"""Quiet-window alarm rates with stations left out, for every detector alike.

    python3 tools/quiet_without_stations.py --index data/quiet100/quiet_index.csv \\
        --runs learned=data/runs/quiet100 stalta=data/runs/quiet100_stalta \\
        --exclude-top 2 --by stalta

Each quiet window's recording (ayzek --record, one NAME.rec per window in each
run directory) is replayed through the network stage without the left-out
stations (tools/network_subsets --quiet), which is the same as running ayzek
without them. The stations left out are given (--exclude CMH,KIRK) or are the
--exclude-top N stations with the most station detections in the --by run over
all windows (held-out protocol, section 7: the two noisiest stations by
STA/LTA's triggers). A window without a left-out station is replayed as it is.
Prints, per detector, the alarms and the rate per day with all stations and
without the left-out ones.

Standard library only.
"""
import argparse
import collections
import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def stations_of(rec):
    return [l.split("\t")[1] for l in rec.read_text().splitlines() if l.startswith("S\t")]


def alarms(binary, rec, keep):
    """Alarms of the network stage on `rec` with exactly the stations `keep`."""
    k = len(keep)
    p = subprocess.run([str(binary), str(rec), "--quiet", "--min-k", str(k), "--max-k", str(k)],
                       capture_output=True, text=True)
    if p.returncode:
        raise SystemExit(f"{rec}: network_subsets failed: {p.stderr[-300:]}")
    for r in csv.DictReader(p.stdout.splitlines()):
        if set(r["stations"].split("+")) == set(keep):
            return int(r["unmatched"])
    raise SystemExit(f"{rec}: no subset {sorted(keep)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--index", required=True, help="quiet_index.csv (name, scored_hours)")
    ap.add_argument("--runs", nargs="+", required=True, metavar="LABEL=DIR")
    ap.add_argument("--exclude", default="", help="comma-separated station codes")
    ap.add_argument("--exclude-top", type=int, default=0)
    ap.add_argument("--by", help="the LABEL whose station detections rank the stations")
    ap.add_argument("--binary", default=str(ROOT / "build-release/tools/network_subsets"))
    a = ap.parse_args()

    index = list(csv.DictReader(open(a.index)))
    runs = dict(x.split("=", 1) for x in a.runs)
    hours = sum(float(w["scored_hours"]) for w in index)
    out = set(s for s in a.exclude.split(",") if s)
    if a.exclude_top:
        if a.by not in runs:
            raise SystemExit("--exclude-top needs --by LABEL of one of --runs")
        count = collections.Counter()
        for w in index:
            for l in (Path(runs[a.by]) / f"{w['name']}.rec").read_text().splitlines():
                if l.startswith("D\t"):
                    count[l.split("\t")[1]] += 1
        top = [s for s, _ in count.most_common(a.exclude_top)]
        print(f"station detections of {a.by}: " +
              ", ".join(f"{s} {n}" for s, n in count.most_common()))
        out |= set(top)
    print(f"left out: {', '.join(sorted(out))}; {len(index)} windows, {hours:.1f} h")
    for label, d in runs.items():
        n_all = n_out = 0
        for w in index:
            rec = Path(d) / f"{w['name']}.rec"
            st = stations_of(rec)
            n_all += alarms(a.binary, rec, st)
            keep = [s for s in st if s not in out]
            n_out += alarms(a.binary, rec, keep) if len(keep) >= 2 else 0
        print(f"  {label:10s} all stations {n_all:5d} alarms, {24 * n_all / hours:7.2f} per day; "
              f"without them {n_out:5d}, {24 * n_out / hours:7.2f} per day")


if __name__ == "__main__":
    sys.exit(main())
