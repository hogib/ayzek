"""Runs ayzek over a grid of settings on several datasets and tabulates the results.

Each --set gives one ayzek option and its values, separated by spaces; the grid
is every combination. For an option without a value, use the values on and off.
Magnitude estimation is disabled, since it does not affect detection or
association.

A dataset is NAME:CACHE:FILES. CACHE is a probability cache from an earlier run
with --scores (the runs then use --scores-in, which gives the same result as
running the detector); leave it empty for STA/LTA. FILES is a comma-separated
list of miniSEED files or globs.

    # learned detector, from cached probabilities
    python3 tools/sweep_trigger.py --catalog CAT \\
        --dataset 'marmara:data/runs/cache_marmara:data/marmara/*.mseed' \\
        --set 'threshold=0.9 0.85 0.8' --set 'trigger-windows=1 4 8' --set 'instant-threshold=2 0.9' --out sweep.csv

    # STA/LTA
    python3 tools/sweep_trigger.py --catalog CAT --args '--detector stalta' \\
        --dataset 'marmara::data/marmara/*.mseed' \\
        --set 'stalta-on=3 4 6' --set 'lta=10 30' --set 'stalta-3c=off on' --out stalta.csv

With --shifts S1,S2,..., every setting is also scored against copies of the
catalogue with all times shifted by +S and -S seconds. The mean detection rate
against the shifted catalogues is the chance level c, and the chance-corrected
detection rate is (o - c) / (1 - c) for the observed rate o.

Per setting the table gives: alarms, alarms without a catalogue event, detected
catalogue events, chance level and corrected count (with --shifts), median
alarm delay after origin, station triggers, detections per magnitude band, and
the alarm time of the largest detected event.

Standard library only.
"""
import argparse
import concurrent.futures
import csv
import datetime as dt
import glob
import itertools
import re
import statistics
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

PATTERNS = {
    "alarms": r"Alarms raised\s+(\d+)",
    "correct": r"correct\s+\(match an AFAD event\)\s+(\d+)",
    "false": r"false\s+\(no AFAD event\)\s+(\d+)",
    "catalogue": r"AFAD events within \d+ km\s+(\d+)",
    "detected": r"^\s+detected\s+(\d+)",
    "missed": r"^\s+missed\s+(\d+)",
    "alarm_delay_s": r"alarm after origin time\s+median ([\d.]+) s",
}
BANDS = ["M < 2", "M 2-3", "M 3-4", "M 4-5", "M 5+"]


def parse_report(out):
    row = {}
    for key, pat in PATTERNS.items():
        m = re.search(pat, out, re.M)
        row[key] = m.group(1) if m else ""
    for band in BANDS:
        m = re.search(rf"^\s+{re.escape(band)}\s+(\d+)\s+(\d+)\s+(\d+)", out, re.M)
        row[band] = f"{m.group(2)}/{m.group(1)}" if m else ""
    # Largest catalogue event that raised an alarm: its alarm delay and stations.
    best = None
    for m in re.finditer(r"^  alarm\s+\S+ UTC by (.*)\n(?:.*\n){0,3}?  AFAD\s+\w+ ([\d.]+) at .*?alarm ([\d.]+) s after origin",
                         out, re.M):
        if best is None or float(m.group(2)) > best[0]:
            best = (float(m.group(2)), m.group(3), m.group(1).replace("; later", " +").replace(",", ""))
    row["main_M"], row["main_alarm_s"], row["main_stations"] = best if best else ("", "", "")
    table = out[out.find("Station processing"):]
    row["triggers"] = sum(int(m.group(1)) for m in re.finditer(r"^  [A-Z0-9]+\s+\d+\s+\d+\s+(\d+)", table, re.M))
    row["dropped_warnings"] = out.count("dropped")
    return row


def record_epoch(rec):
    year, doy, hour, minute, sec, _, frac = struct.unpack(">HHBBBBH", rec[20:30])
    t = dt.datetime(year, 1, 1, tzinfo=dt.UTC) + dt.timedelta(
        days=doy - 1, hours=hour, minutes=minute, seconds=sec, microseconds=frac * 100)
    return t.timestamp()


def data_span(files):
    """First and last record start times over the files (records are in time order)."""
    t0, t1 = float("inf"), float("-inf")
    for f in files:
        with open(f, "rb") as fh:
            first = fh.read(512)
            fh.seek(-512, 2)
            last = fh.read(512)
        t0, t1 = min(t0, record_epoch(first)), max(t1, record_epoch(last))
    return t0, t1


def write_catalogue(src, dst, span, shift):
    """Copies the AFAD catalogue events within `span` (plus a day on either side),
    with all times shifted by `shift` seconds."""
    lo, hi = span[0] - 86400, span[1] + 86400
    with open(src, encoding="utf-8-sig") as f, open(dst, "w") as g:
        g.write(next(f))
        for line in f:
            date, rest = line.split(",", 1)
            try:
                t = dt.datetime.strptime(date, "%d/%m/%Y %H:%M:%S").replace(tzinfo=dt.UTC)
            except ValueError:
                continue
            if lo <= t.timestamp() <= hi:
                g.write((t + dt.timedelta(seconds=shift)).strftime("%d/%m/%Y %H:%M:%S") + "," + rest)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--binary", default="build-release/app/ayzek")
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--dataset", action="append", required=True, help="NAME:CACHE:FILE[,FILE...]; CACHE may be empty")
    ap.add_argument("--set", action="append", default=[], metavar="OPTION=V1 V2 ...",
                    help="ayzek option (without --) and values; 'on'/'off' for options without a value")
    ap.add_argument("--args", default="", help="options passed to every run")
    ap.add_argument("--shifts", default="", help="catalogue time shifts in seconds for the chance level, e.g. 127,263,391")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    options = []
    for spec in a.set:
        name, values = spec.split("=", 1)
        options.append((name.strip(), values.split()))
    settings = [list(zip([o for o, _ in options], combo)) for combo in itertools.product(*[v for _, v in options])]
    shifts = [s * sign for s in (int(x) for x in a.shifts.split(",") if x) for sign in (1, -1)]

    tmpdir = tempfile.TemporaryDirectory(prefix="sweep_")
    tmp = Path(tmpdir.name)
    datasets = []
    for spec in a.dataset:
        name, cache, files = spec.split(":", 2)
        paths = sorted(p for f in files.split(",") for p in glob.glob(f))
        if not paths:
            sys.exit(f"{name}: no files")
        span = data_span(paths)
        catalogues = {}
        for shift in [0, *shifts]:
            catalogues[shift] = tmp / f"{name}_{shift}.csv"
            write_catalogue(a.catalog, catalogues[shift], span, shift)
        datasets.append((name, cache, paths, catalogues))

    def run(dataset, setting, shift):
        name, cache, files, catalogues = dataset
        cmd = [a.binary, "--speed", "0", "--no-color", "--no-magnitude", "--catalog", str(catalogues[shift]), *a.args.split()]
        if cache:
            cmd += ["--scores-in", cache]
        for option, value in setting:
            if value == "on":
                cmd.append(f"--{option}")
            elif value != "off":
                cmd += [f"--{option}", value]
        out = subprocess.run([*cmd, *files], capture_output=True, text=True).stdout
        return name, tuple(setting), shift, parse_report(out)

    jobs = [(d, s, shift) for d in datasets for s in settings for shift in [0, *shifts]]
    results = {}
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as pool:
        for name, setting, shift, row in pool.map(lambda j: run(*j), jobs):
            results[(name, setting, shift)] = row

    tmpdir.cleanup()

    rows = []
    for name, _, _, _ in datasets:
        for setting in settings:
            row = {"dataset": name, **{o: v for o, v in setting}, **results[(name, tuple(setting), 0)]}
            det, cat = int(row["detected"] or 0), int(row["catalogue"] or 0)
            if shifts and cat:
                chance = statistics.mean(int(results[(name, tuple(setting), s)]["detected"] or 0) /
                                         max(1, int(results[(name, tuple(setting), s)]["catalogue"] or 0)) for s in shifts)
                corrected = (det / cat - chance) / (1 - chance) if chance < 1 else 0.0
                row["chance"], row["corrected"] = f"{chance:.3f}", f"{corrected * cat:.1f}"
            else:
                row["chance"], row["corrected"] = "", ""
            rows.append(row)

    fields = ["dataset", *[o for o, _ in options], *PATTERNS, "chance", "corrected", *BANDS,
              "main_M", "main_alarm_s", "main_stations", "triggers", "dropped_warnings"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    label_width = max([len(" ".join(f"{o}={v}" for o, v in s)) for s in settings] + [7])
    for name, _, _, _ in datasets:
        print(f"\n{name}")
        print(f"  {'setting':<{label_width}} {'alarms':>6} {'unmatched':>9} {'detected':>9} {'chance':>6} {'corrected':>9} "
              f"{'delay s':>7} {'trig':>5}   " + "  ".join(f"{b:>6}" for b in BANDS) + "   largest event alarm")
        for r in rows:
            if r["dataset"] != name:
                continue
            label = " ".join(f"{o}={r[o]}" for o, _ in options)
            chance = f"{float(r['chance']):.0%}" if r["chance"] else ""
            print(f"  {label:<{label_width}} {r['alarms']:>6} {r['false']:>9} {r['detected'] + '/' + r['catalogue']:>9} "
                  f"{chance:>6} {r['corrected']:>9} {r['alarm_delay_s']:>7} {r['triggers']:>5}   "
                  + "  ".join(f"{r[b]:>6}" for b in BANDS)
                  + (f"   M{r['main_M']} +{r['main_alarm_s']} s ({r['main_stations']})" if r["main_M"] else ""))
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
