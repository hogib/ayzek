"""Runs ayzek over a grid of trigger settings and tabulates the results.

Each dataset needs a probability cache from one earlier run with --scores DIR;
the sweep runs use --scores-in DIR, which gives the same result as running the
detector, in a fraction of the time. Magnitude estimation is disabled, since it
does not affect detection or association.

    python3 tools/sweep_trigger.py --catalog CAT \\
        --dataset sindirgi:data/runs/cache_sindirgi:data/demo/DEMI.mseed,data/demo/MANT.mseed,data/demo/BAND.mseed \\
        --dataset marmara:data/runs/cache_marmara:data/marmara/*.mseed \\
        --thresholds 0.9,0.85,0.8,0.7 --windows 1,4,8 --instant 2,0.9 --out sweep.csv

An instant threshold above 1 disables the single-window rule. The printed table
gives, per setting: alarms, alarms without a catalogue event, detected
catalogue events, median alarm delay, station triggers, detections per
magnitude band, and the alarm time of the largest detected event.

Standard library only.
"""
import argparse
import concurrent.futures
import csv
import glob
import itertools
import re
import subprocess
import sys

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


def run(binary, catalog, name, cache, files, threshold, windows, instant):
    cmd = [binary, "--speed", "0", "--no-color", "--no-magnitude", "--catalog", catalog, "--scores-in", cache,
           "--threshold", str(threshold), "--trigger-windows", str(windows), "--instant-threshold", str(instant), *files]
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    row = {"dataset": name, "threshold": threshold, "windows": windows, "instant": instant}
    for key, pat in PATTERNS.items():
        m = re.search(pat, out, re.M)
        row[key] = m.group(1) if m else ""
    for band in BANDS:
        m = re.search(rf"^\s+{re.escape(band)}\s+(\d+)\s+(\d+)\s+(\d+)", out, re.M)
        row[band] = f"{m.group(2)}/{m.group(1)}" if m else ""
    # Largest catalogue event that raised an alarm: its alarm delay and stations.
    best = None
    for m in re.finditer(r"^  alarm\s+\S+ UTC by (.*)\n(?:.*\n){0,3}?  AFAD\s+\w+ ([\d.]+) at .*?alarm ([\d.]+) s after origin", out, re.M):
        if best is None or float(m.group(2)) > best[0]:
            best = (float(m.group(2)), m.group(3), m.group(1).replace("; later", " +").replace(",", ""))
    row["main_M"], row["main_alarm_s"], row["main_stations"] = best if best else ("", "", "")
    table = out[out.find("Station processing"):]
    row["triggers"] = sum(int(m.group(1)) for m in re.finditer(r"^  [A-Z0-9]+\s+\d+\s+\d+\s+(\d+)", table, re.M))
    row["dropped_warnings"] = out.count("dropped")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default="build-release/app/ayzek")
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--dataset", action="append", required=True, help="NAME:CACHE_DIR:FILE[,FILE...] (globs allowed)")
    ap.add_argument("--thresholds", default="0.9,0.85,0.8,0.7")
    ap.add_argument("--windows", default="1,2,3")
    ap.add_argument("--instant", default="2", help="instant thresholds; values > 1 disable the single-window rule")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    datasets = []
    for spec in a.dataset:
        name, cache, files = spec.split(":", 2)
        paths = sorted(p for f in files.split(",") for p in glob.glob(f))
        if not paths:
            sys.exit(f"{name}: no files")
        datasets.append((name, cache, paths))
    grid = list(itertools.product(datasets, [float(x) for x in a.thresholds.split(",")],
                                  [int(x) for x in a.windows.split(",")], [float(x) for x in a.instant.split(",")]))

    rows = []
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as pool:
        futures = [pool.submit(run, a.binary, a.catalog, name, cache, files, thr, n, inst)
                   for (name, cache, files), thr, n, inst in grid]
        for f in concurrent.futures.as_completed(futures):
            rows.append(f.result())
    rows.sort(key=lambda r: (r["dataset"], -r["threshold"], r["windows"], r["instant"]))

    fields = ["dataset", "threshold", "windows", "instant", *PATTERNS, *BANDS,
              "main_M", "main_alarm_s", "main_stations", "triggers", "dropped_warnings"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    for name, _, _ in datasets:
        print(f"\n{name}")
        print(f"  {'thr':>5} {'N':>2} {'inst':>4} {'alarms':>7} {'false':>6} {'detected':>9} {'missed':>7} {'delay s':>8} {'trig':>5}   "
              + "  ".join(f"{b:>6}" for b in BANDS) + "   main event alarm")
        for r in rows:
            if r["dataset"] != name:
                continue
            print(f"  {r['threshold']:>5} {r['windows']:>2} {('-' if r['instant'] > 1 else r['instant']):>4} {r['alarms']:>7} {r['false']:>6} "
                  f"{r['detected'] + '/' + r['catalogue']:>9} {r['missed']:>7} {r['alarm_delay_s']:>8} {r['triggers']:>5}   "
                  + "  ".join(f"{r[b]:>6}" for b in BANDS) + f"   M{r['main_M']} +{r['main_alarm_s']} s ({r['main_stations']})")
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
