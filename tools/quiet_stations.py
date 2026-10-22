#!/usr/bin/env python3
"""Detections per station on the quiet windows, for two detectors (supplement Table S3).

    python3 tools/quiet_stations.py --runs learned=data/runs/quiet100 \\
        stalta=data/runs/quiet100_stalta --data data/quiet100 > table_s3.csv

For each station: the scored hours of the windows it is part of (each window
from the end of its lead-in, `scored_hours` of quiet_index.csv), and per
detector the station detections after the lead-in with bursts within 40 s
counted once (as tools/misfire_review.py does), their rate per day, and the
alarms the station declared (named on the ALARM line of the report).
Standard library only.
"""
import argparse
import collections
import csv
import glob
import re
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from misfire_review import cluster  # noqa: E402

ALARM = re.compile(r"^ALARM\s+#\d+\s+\S+\s+earthquake detected by ([A-Z0-9, ]+?)\s*(?:\(|$)", re.M)


def count(run, index):
    det, alarms = collections.Counter(), collections.Counter()
    for w in index:
        t0 = datetime.fromisoformat(w["score_from_utc"]).timestamp()
        times = collections.defaultdict(list)
        for line in open(f"{run}/{w['name']}.rec"):
            f = line.split("\t")
            if f[0] == "D" and float(f[3]) >= t0:
                times[f[1]].append(float(f[3]))
        for s, ts in times.items():
            det[s] += len(cluster(sorted(ts), 40.0))
        for named in ALARM.findall(Path(f"{run}/{w['name']}.txt").read_text()):
            for s in {x.strip() for x in named.split(",")}:
                alarms[s] += 1
    return det, alarms


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", nargs="+", required=True, metavar="LABEL=DIR")
    ap.add_argument("--data", default="data/quiet100", help="NAME/STATION.mseed and quiet_index.csv")
    a = ap.parse_args()

    index = list(csv.DictReader(open(f"{a.data}/quiet_index.csv")))
    hours = collections.Counter()
    for w in index:
        for p in glob.glob(f"{a.data}/{w['name']}/*.mseed"):
            hours[Path(p).stem] += float(w["scored_hours"])
    runs = [x.split("=", 1) for x in a.runs]
    counts = {label: count(d, index) for label, d in runs}

    out = csv.writer(sys.stdout)
    out.writerow(["station", "hours"] + [f"{label}_{k}" for label, _ in runs
                                         for k in ("detections", "per_day", "alarms")])
    for s in sorted(hours):
        row = [s, f"{hours[s]:.1f}"]
        for label, _ in runs:
            det, alarms = counts[label]
            row += [det[s], f"{24 * det[s] / hours[s]:.1f}", alarms[s]]
        out.writerow(row)
    out.writerow(["total", ""] + [x for label, _ in runs
                                  for x in (sum(counts[label][0].values()), "", "")])


if __name__ == "__main__":
    main()
