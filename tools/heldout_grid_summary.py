"""Summary of an exploratory grid of tools/heldout_grid.py: EXPLORATORY, outside the
held-out protocol.

    python3 tools/heldout_grid_summary.py --grid tradeoff [--out DIR]      (from ayzek_code/)
    python3 tools/heldout_grid_summary.py --grid instant  [--out DIR]

DIR is the grid's output directory (default: data/runs/heldout_grid for tradeoff,
data/runs/heldout_instant for instant; the former drivers logs/grid_summary.py and
logs/instant_summary.py, which differed in that directory alone).

Inputs: DIR/grid_settings.json (the settings and their order), DIR/<net>/{runs,quiet,
outcomes} of the grid, and of the release run data/runs/heldout/<net>/quiet/
quiet_index.csv and outcomes/learned_*.csv.

Per network and grid setting, from the grid's reports: detected and scorable events and
unmatched alarms over the sequence sets, alarms and scored hours over the quiet windows.
Then tools/significance.py, unchanged, compares the frozen learned detector (release
outcome tables, linked beside the grid's) with every grid setting, --against all of
them in the order of grid_settings.json: chance-corrected counts and the paired
bootstrap interval of the difference.

Writes DIR/summary.csv (through summary.csv.tmp) and DIR/<net>/analysis/
significance.{txt,json}. One line per finished network. Exits nonzero when
significance.py fails or nothing is produced. Standard library only.
"""
import argparse
import csv
import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

DEFAULT_OUT = {"tradeoff": "data/runs/heldout_grid", "instant": "data/runs/heldout_instant"}   # tools/heldout_grid.py
R = Path("data/runs/heldout")
NETS = ("greece_now", "nz_now", "ko_now")
SIG = re.compile(r"^learned vs (\S+): (\d+) events\n  detected\s+(\d+) vs (\d+).*\n  chance-corrected\s+([\d.]+) vs ([\d.]+), "
                 r"difference ([+-][\d.]+) \(95% bootstrap ([+-][\d.]+) to ([+-][\d.]+)", re.M)


def num(pat, t):
    return int(m.group(1)) if (m := re.search(pat, t, re.M)) else 0


def parse_significance(text):
    """{setting: (setting, events, detected learned, detected setting, corrected
    learned, corrected setting, difference, interval low, interval high)}, the
    figures as significance.py prints them."""
    return {m.group(1): m.groups() for m in SIG.finditer(text)}


def report_counts(seq, qt):
    """The counts of one setting from its reports' texts (sequence sets, quiet windows)."""
    return {"alarms": sum(num(r"^  Alarms raised\s+(\d+)", t) for t in seq),
            "unmatched_alarms": sum(num(r"^    false\s+\(no \w+ event\)\s+(\d+)", t) for t in seq),
            "quiet_alarms": sum(num(r"^  Alarms raised\s+(\d+)", t) for t in qt)}


def one(net, G, settings):
    hours = sum(float(r["scored_hours"]) for r in csv.DictReader(open(R / net / "quiet/quiet_index.csv")))
    link = G / net / "analysis" / "outcomes"
    link.mkdir(parents=True, exist_ok=True)
    for f in list((G / net / "outcomes").glob("*.csv")) + list((R / net / "outcomes").glob("learned_*.csv")):
        if not (link / f.name).exists():
            (link / f.name).symlink_to(f.resolve())
    seqs = sorted(f.stem[len("learned_"):] for f in (R / net / "outcomes").glob("learned_*.csv"))
    out = G / net / "analysis" / "significance"
    p = subprocess.run([sys.executable, "tools/significance.py", "--outcomes", str(link),
                        "--groups", f"{net}={','.join(seqs)}", "--against", ",".join(settings),
                        "--json", str(out) + ".json"], capture_output=True, text=True)
    if p.returncode:
        sys.exit(p.stderr[-2000:])
    Path(str(out) + ".txt").write_text(p.stdout)
    sig = parse_significance(p.stdout)
    rows = []
    for s in settings:
        seq = [(f.read_text()) for f in sorted((G / net / "runs" / s).glob("*.txt"))]
        qt = [(f.read_text()) for f in sorted((G / net / "quiet" / s).glob("*.txt"))]
        g = sig[s]
        c = report_counts(seq, qt)
        rows.append({"network": net, "setting": s, "sequence_sets": len(seq), "events": g[1], "detected": g[3],
                     "corrected": g[5], "alarms": c["alarms"],
                     "unmatched_alarms": c["unmatched_alarms"],
                     "quiet_windows": len(qt), "quiet_alarms": c["quiet_alarms"],
                     "quiet_hours": f"{hours:.1f}", "frozen_learned_corrected": g[4],
                     "difference": g[6], "ci_low": g[7], "ci_high": g[8]})
    print(net, "done", flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--grid", required=True, choices=list(DEFAULT_OUT))
    ap.add_argument("--out", help="the grid's output directory (default: data/runs/heldout_grid for tradeoff, "
                                  "data/runs/heldout_instant for instant)")
    a = ap.parse_args()
    G = Path(a.out or DEFAULT_OUT[a.grid])
    settings = [k for k in json.loads((G / "grid_settings.json").read_text()) if not k.startswith("_")]
    with ThreadPoolExecutor(3) as pool:
        rows = [r for rs in pool.map(lambda n: one(n, G, settings), NETS) for r in rs]
    if not rows:
        sys.exit("nothing produced")
    tmp = G / "summary.csv.tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    tmp.rename(G / "summary.csv")
    print("wrote", G / "summary.csv", len(rows), "rows")


if __name__ == "__main__":
    main()
