#!/usr/bin/env python3
"""Chance level per detector, replay and shift, and event-avoiding shifts.

The chance correction (paper Equation 4) measures how often an alarm sequence
matches a catalogue shifted in time. Inside an aftershock sequence a shifted
event can land next to a real one, and a true detection of that real event then
"matches" it (referee round 5, M5). This tool reads the per-event outcome files
of tools/event_outcomes.py and

  report   prints the chance level c of every detector on every replay at every
           shift, with the share of shifted events that land within --near
           seconds of a real catalogue event;
  filter   writes a copy of the outcome files from which every shifted event
           that lands within --near seconds of a real one is removed, for
           tools/significance.py. The observed rows are unchanged, so only c
           moves. The removal depends on times alone, so it is the same for
           every detector.

Real events are the catalogue events of the unshifted run (those within 250 km
in the replay); a shifted event at origin o and shift s lands at o + s.

    python3 tools/chance_sensitivity.py report --outcomes data/runs/r5/outcomes --near 60
    python3 tools/chance_sensitivity.py filter --outcomes data/runs/r5/outcomes \\
        --near 60 --out data/runs/r5/outcomes_avoid60

Standard library only.
"""
import argparse
import collections
import csv
import glob
import os
from pathlib import Path


def seconds(origin):
    h, m, s = origin.split(":")
    return int(h) * 3600 + int(m) * 60 + float(s)


def load(d):
    files = sorted(glob.glob(f"{d}/*.csv"))
    rows = {f: list(csv.DictReader(open(f))) for f in files}
    real = collections.defaultdict(set)
    for rs in rows.values():
        for r in rs:
            if int(r["shift"]) == 0:
                real[r["dataset"]].add(seconds(r["origin"]))
    return rows, {k: sorted(v) for k, v in real.items()}


def near(t, times, w):
    return any(abs(t - x) < w for x in times)


def report(a):
    rows, real = load(a.outcomes)
    print(f"chance level c per shift; [n near]: shifted events landing within {a.near:g} s "
          f"of a real event")
    for f, rs in rows.items():
        if not rs:
            continue
        by = collections.defaultdict(list)
        for r in rs:
            by[int(r["shift"])].append(r)
        ds, det = rs[0]["dataset"], rs[0]["detector"]
        o = by[0]
        cells, cs, cs_far = [], [], []
        for s in sorted(k for k in by if k):
            g = by[s]
            close = [r for r in g if near(seconds(r["origin"]) + s, real[ds], a.near)]
            far = [r for r in g if r not in close]
            c = sum(int(r["detected"]) for r in g) / len(g)
            cf = sum(int(r["detected"]) for r in far) / len(far) if far else float("nan")
            cs.append(c)
            cs_far.append(cf)
            cells.append(f"{s:+d} {c:.3f} [{len(close)}/{len(g)}]")
        ob = sum(int(r["detected"]) for r in o) / len(o)
        c, cf = sum(cs) / len(cs), sum(cs_far) / len(cs_far)
        corr = (ob - c) / (1 - c) * len(o)
        corr_f = (ob - cf) / (1 - cf) * len(o)
        print(f"{ds:11s} {det:15s} o {ob:.3f} c {c:.3f} (far {cf:.3f}) corrected {corr:5.1f} "
              f"(far {corr_f:5.1f}) | " + "  ".join(cells))


def filter_(a):
    rows, real = load(a.outcomes)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    dropped = collections.Counter()
    for f, rs in rows.items():
        keep = [r for r in rs if int(r["shift"]) == 0
                or not near(seconds(r["origin"]) + int(r["shift"]), real[r["dataset"]], a.near)]
        if rs:
            dropped[rs[0]["dataset"], rs[0]["detector"]] = len(rs) - len(keep)
        with open(out / os.path.basename(f), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=rs[0].keys() if rs else ["dataset"])
            w.writeheader()
            w.writerows(keep)
    for (ds, det), n in sorted(dropped.items()):
        if det == "learned":
            print(f"{ds:11s} {n} shifted events within {a.near:g} s of a real one removed")
    print(f"-> {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("report", "filter"):
        q = sub.add_parser(name)
        q.add_argument("--outcomes", required=True)
        q.add_argument("--near", type=float, default=60.0,
                       help="seconds between a shifted event and a real one below which it is removed")
        if name == "filter":
            q.add_argument("--out", required=True)
    a = ap.parse_args()
    (report if a.cmd == "report" else filter_)(a)


if __name__ == "__main__":
    main()
