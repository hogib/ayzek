"""Per-station STA/LTA trigger ratios: the per-station baseline of the held-out
protocol (section 4).

    python3 tools/stalta_station_ratios.py --binary data/runs/tune_dev/ayzek \\
        --args '--detector stalta --stalta-band 1,10 --sta 0.5 --lta 10 --stalta-off 1.5' \\
        --on 8 --tuning ../ayzek_heldout/data/tune_* --out data/runs/stalta_stations

Each tuning set (a set directory of the held-out repository: STATION.mseed,
stations.csv, pd_gain.csv, catalog.csv, dataset.json, named tune_GROUP_N; its
first 30 min warm the filters and are not counted) is replayed with the
network-wide setting (--args: the equal-unmatched STA/LTA setting of the
development data) at every ratio of --grid. A station's trigger rate is its
station detections per hour of quiet data over all the sets of its group. The
group's target is the median of its stations' rates at the network-wide ratio
(--on), and each station is given the lowest ratio of the grid at which its
rate does not exceed the target, or the highest ratio if none does: noisy
stations rise above the network-wide ratio and quiet ones fall below it, the
group's overall rate staying near the network-wide setting's.

Quiet data: a station's time within [D/8 - 10 s, D/6 + 120 s] after the origin
of a catalogued event it could have recorded (epicentral distance D in km; the
visibility rule of tools/compare_detectors.py) is left out with its triggers,
since the tuning windows are random 6 h windows, so that the rate counts noise
rather than earthquakes. Hours are the span, not the data present; the tuning
sets' stations have at least 90% coverage.

Writes OUT/targets.csv (each group's stations and target), OUT/rates_GROUP.csv
(each station's quiet hours, its triggers and rate at every ratio, and its
ratio) and OUT/stalta_on_GROUP.csv (station,ratio: ayzek's
--stalta-on-stations).

Standard library only.
"""
import argparse
import concurrent.futures
import csv
import datetime as dt
import json
import math
import re
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

VISIBILITY = ((50.0, 0.0), (150.0, 2.0), (400.0, 3.0), (1500.0, 4.5))  # compare_detectors.py
LEAD_S = 1800.0


def epoch(text):
    """An ISO time (UTC if naive) or the catalogue's DD/MM/YYYY HH:MM:SS[.f], in s."""
    text = text.strip()
    if "/" in text:
        day, clock = text.split()
        d, m, y = day.split("/")
        text = f"{y}-{m}-{d}T{clock}"
    t = dt.datetime.fromisoformat(text)
    return (t if t.tzinfo else t.replace(tzinfo=dt.UTC)).timestamp()


def km(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    h = (math.sin((p2 - p1) / 2) ** 2
         + math.cos(p1) * math.cos(p2) * math.sin(math.radians(lon2 - lon1) / 2) ** 2)
    return 2 * 6371.0 * math.asin(math.sqrt(h))


def catalogue(path):
    return [(epoch(r["Date"]), float(r["Latitude"]), float(r["Longitude"]), float(r["Magnitude"]))
            for r in csv.DictReader(open(path, encoding="utf-8-sig"))]


def arrivals(events, lat, lon, t0, t1):
    """Merged spans of [t0, t1) within the arrivals of the events visible at (lat, lon)."""
    spans = []
    for t, la, lo, mag in events:
        d = km(lat, lon, la, lo)
        if any(d <= dd and mag >= mm for dd, mm in VISIBILITY):
            a, b = max(t0, t + d / 8.0 - 10.0), min(t1, t + d / 6.0 + 120.0)
            if a < b:
                spans.append((a, b))
    merged = []
    for a, b in sorted(spans):
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return merged


def replay(binary, args, files):
    """The stations (code -> lat, lon) and each one's detection times of one ayzek run."""
    with tempfile.TemporaryDirectory(prefix="stalta_stations_") as tmp:
        rec = Path(tmp) / "run.rec"
        cmd = [binary, "--speed", "0", "--no-color", "--record", str(rec), *args, *files]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode or not rec.exists():
            raise SystemExit(f"ayzek failed ({p.returncode}): {' '.join(cmd[:16])} ...\n"
                             f"{p.stderr[-500:]}")
        stations, times = {}, defaultdict(list)
        for line in rec.read_text().splitlines():
            f = line.split("\t")
            if f[0] == "S":
                stations[f[1]] = (float(f[2]), float(f[3]))
            elif f[0] == "D":
                times[f[1]].append(float(f[3]))  # declared_at
    return stations, times


def quiet_counts(stations, times, events, t0, t1):
    """code -> (quiet hours, triggers in them) over [t0, t1)."""
    out = {}
    for code, (lat, lon) in stations.items():
        spans = arrivals(events, lat, lon, t0, t1)
        hours = (t1 - t0 - sum(b - a for a, b in spans)) / 3600.0
        n = sum(1 for t in times.get(code, ())
                if t0 <= t < t1 and not any(a <= t < b for a, b in spans))
        out[code] = (hours, n)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--binary", default="build-release/app/ayzek")
    ap.add_argument("--args", required=True, help="the network-wide setting, without --stalta-on")
    ap.add_argument("--on", type=float, required=True, help="its trigger ratio")
    ap.add_argument("--grid", default="3,4,5,6,8,10,12,15", help="the ratios a station may get")
    ap.add_argument("--tuning", nargs="+", required=True, help="tuning set directories")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    args = a.args.split()
    if "--stalta-on" in args or "--stalta-on-stations" in args:
        sys.exit("--args: the setting without its ratio, which is --on")
    grid = sorted(float(x) for x in a.grid.split(","))
    if a.on not in grid:
        sys.exit(f"--on {a.on:g} must be one of the grid's ratios")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    sets = []
    for d in map(Path, sorted(a.tuning)):
        m = re.fullmatch(r"tune_([a-z]+)_\d+", d.name)
        files = sorted(str(p) for p in d.glob("*.mseed"))
        if not m or not files or not (d / "dataset.json").exists():
            sys.exit(f"{d}: not a fetched tuning set (tune_GROUP_N with miniSEED and dataset.json)")
        span = json.loads((d / "dataset.json").read_text())["span"]
        sets.append((m.group(1), d, files, epoch(span[0]) + LEAD_S, epoch(span[1])))

    def tuning_window(d, files, t0, t1, ratio):
        run_args = [*args, "--stalta-on", f"{ratio:g}", "--stations", str(d / "stations.csv")]
        if (d / "pd_gain.csv").exists():
            run_args += ["--pd-gains", str(d / "pd_gain.csv")]
        stations, times = replay(a.binary, run_args, files)
        return quiet_counts(stations, times, catalogue(d / "catalog.csv"), t0, t1)

    with concurrent.futures.ThreadPoolExecutor(a.jobs) as pool:
        tf = {(g, d.name, r): pool.submit(tuning_window, d, files, t0, t1, r)
              for g, d, files, t0, t1 in sets for r in grid}
        tuning = {k: f.result() for k, f in tf.items()}

    targets = []
    for group in sorted({g for g, *_ in sets}):
        acc = defaultdict(lambda: defaultdict(lambda: [0.0, 0]))  # code -> ratio -> [h, n]
        for (g, _, r), counts in tuning.items():
            if g == group:
                for code, (h, n) in counts.items():
                    acc[code][r][0] += h
                    acc[code][r][1] += n
        rate = {code: {r: acc[code][r][1] / acc[code][r][0] if acc[code][r][0] > 0 else math.inf
                       for r in grid} for code in acc}
        target = statistics.median(v[a.on] for v in rate.values() if math.isfinite(v[a.on]))
        rows = []
        for code in sorted(acc):
            meets = [r for r in grid if rate[code][r] <= target]
            rows.append((code, acc[code][grid[0]][0], min(meets) if meets else grid[-1],
                         bool(meets)))
        with open(out / f"rates_{group}.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["station", "quiet_hours", *[f"triggers_{r:g}" for r in grid],
                        *[f"per_hour_{r:g}" for r in grid], "ratio", "meets_target"])
            for code, hours, ratio, ok in rows:
                w.writerow([code, f"{hours:.2f}", *[acc[code][r][1] for r in grid],
                            *[f"{rate[code][r]:.4f}" for r in grid], f"{ratio:g}", int(ok)])
        with open(out / f"stalta_on_{group}.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["station", "ratio"])
            for code, _, ratio, _ in rows:
                w.writerow([code, f"{ratio:g}"])
        hist = defaultdict(int)
        for *_, ratio, _ in rows:
            hist[ratio] += 1
        n_sets = len({d for g, d, *_ in sets if g == group})
        targets.append((group, len(rows), n_sets, target))
        print(f"{group}: {len(rows)} stations over {n_sets} sets, target {target:.4f} triggers "
              f"per quiet hour at ratio {a.on:g}; ratios "
              + ", ".join(f"{r:g}: {hist[r]}" for r in grid)
              + f"; {sum(1 for *_, ok in rows if not ok)} above the target at every ratio")
    with open(out / "targets.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["group", "stations", "sets", f"median_per_hour_at_{a.on:g}"])
        for group, n, n_sets, target in targets:
            w.writerow([group, n, n_sets, f"{target:.4f}"])


if __name__ == "__main__":
    main()
