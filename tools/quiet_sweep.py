#!/usr/bin/env python3
"""Alarms on the quiet windows for every setting of a trigger grid.

    # STA/LTA
    python3 tools/quiet_sweep.py --args '--detector stalta --stalta-off 1.5' \\
        --set 'stalta-band=1,10 2,10 2,20' --set 'sta=0.5 1' --set 'lta=10 30 60' \\
        --set 'stalta-on=3 4 5 6 8 10' --set 'stalta-3c=off on' --out quiet_stalta.csv
    # learned detector, from window probabilities written by `ayzek --scores ROOT/NAME`
    python3 tools/quiet_sweep.py --scores-root data/runs/quiet100_scores \\
        --set 'threshold=0.7 0.8 0.9' --set 'trigger-windows=1 8' --out quiet_learned.csv

Each window of quiet_index.csv is replayed on its own from the end of its
lead-in, as data/runs/quiet100/run.sh does, against the catalogue excerpt
tests/catalogs/quiet100.csv; every alarm counts. One row per setting: the
alarms summed over the windows, the scored hours, and alarms per day. The grid
syntax is tools/sweep_trigger.py's ('on'/'off' for options without a value).
Standard library only.
"""
import argparse
import concurrent.futures
import csv
import glob
import itertools
import re
import subprocess
from pathlib import Path

ALARMS = re.compile(r"Alarms raised\s+(\d+)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--binary", default="build-release/app/ayzek")
    ap.add_argument("--data", default="data/quiet100")
    ap.add_argument("--catalog", default="tests/catalogs/quiet100.csv")
    ap.add_argument("--scores-root", help="ROOT/NAME probability caches for the learned detector")
    ap.add_argument("--set", action="append", default=[], metavar="OPTION=V1 V2 ...")
    ap.add_argument("--args", default="", help="options passed to every run")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    options = [(n.strip(), v.split()) for n, v in (spec.split("=", 1) for spec in a.set)]
    settings = [list(zip([o for o, _ in options], combo)) for combo in itertools.product(*[v for _, v in options])]
    windows = list(csv.DictReader(open(f"{a.data}/quiet_index.csv")))
    hours = sum(float(w["scored_hours"]) for w in windows)

    def run(setting, w):
        cmd = [a.binary, "--speed", "0", "--no-color", "--no-magnitude", "--catalog", a.catalog,
               "--from", w["score_from_utc"].split("+")[0], *a.args.split()]
        if a.scores_root:
            cmd += ["--scores-in", f"{a.scores_root}/{w['name']}"]
        for option, value in setting:
            if value == "on":
                cmd.append(f"--{option}")
            elif value != "off":
                cmd += [f"--{option}", value]
        out = subprocess.run([*cmd, *sorted(glob.glob(f"{a.data}/{w['name']}/*.mseed"))],
                             capture_output=True, text=True).stdout
        m = ALARMS.search(out)
        if not m:
            raise SystemExit(f"no report for {w['name']} at {setting}")
        return tuple(setting), int(m.group(1))

    jobs = [(s, w) for s in settings for w in windows]
    alarms = {}
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as pool:
        for setting, n in pool.map(lambda j: run(*j), jobs):
            alarms[setting] = alarms.get(setting, 0) + n

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([o for o, _ in options] + ["alarms", "hours", "per_day"])
        for s in settings:
            n = alarms[tuple(s)]
            w.writerow([v for _, v in s] + [n, f"{hours:.2f}", f"{24 * n / hours:.2f}"])
    print(f"{len(settings)} settings x {len(windows)} windows -> {a.out}")


if __name__ == "__main__":
    main()
