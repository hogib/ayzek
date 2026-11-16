#!/usr/bin/env python3
"""Alarm rates on the stratified quiet windows (EXPLORATORY, outside the held-out protocol).

    python3 tools/heldout_quietstrat.py --draw data/runs/heldout_quietstrat/greece_strat \\
        --run data/runs/heldout_quietstrat/greece --settings data/runs/heldout_quietstrat/settings_all20.json \\
        --groups "Crete and the southern Aegean=4,13,17;mainland=2,3,8,9,10,11,14,20"

The stratified windows were drawn after the tests, the same number on every selected
sequence's station set (ayzek_heldout/tools/quiet_draw.py), and replayed with every
frozen setting by tools/run_heldout.py (stage runs). This tool only counts: per window
the alarms of each detector (the report's "Alarms raised" line, as tools/significance.py
reads it) and the scored hours (the window less its first 30 min, as the protocol), and
from them, with tools/significance.py's own functions,

  - per detector: alarms, rate per day, and for every detector but the first of
    --reference (default learned) the ratio of the reference's rate to its rate with the
    exact conditional 95 % interval and the interval from resampling windows;
  - the same counts per group of station sets (--groups NAME=SEQ,SEQ;...; a window
    belongs to the sequence whose station set it was drawn on) and per station set.

Writes RUN/analysis/quietstrat_windows.csv (one row per window), quietstrat_summary.csv,
quietstrat_groups.csv and quietstrat.json (with argv, the commit and the input hashes),
each through a temporary name. Exits nonzero when a report is missing or incomplete.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import significance  # noqa: E402

QUIET_LEAD_S = 1800.0      # tools/run_heldout.py


def windows_of(draw):
    sel = json.loads((Path(draw) / "selection.json").read_text())
    tag = f"{sel['group']}_{sel['name']}"
    out = []
    for j, w in enumerate(sel["quiet"]["windows"], 1):
        t0, t1 = dt.datetime.fromisoformat(w["start"]), dt.datetime.fromisoformat(w["end"])
        out.append({"name": f"{tag}_quiet{j:02d}", "sequence": int(w["sequence"]),
                    "start": w["start"], "hours": (t1 - t0).total_seconds() / 3600 - QUIET_LEAD_S / 3600})
    return out


def parse_groups(text):
    groups = {}
    for part in filter(None, (p.strip() for p in text.split(";"))):
        name, seqs = part.split("=")
        groups[name.strip()] = {int(s) for s in seqs.split(",")}
    return groups


def rates(windows, counts, labels, reference, rng):
    """One row per label: alarms, hours, per day, and against the reference the ratio
    reference / label with its exact and window-bootstrap intervals."""
    hours = {w["name"]: w["hours"] for w in windows}
    T = sum(hours.values())
    rows = []
    for lab in labels:
        n = sum(counts[lab][w] for w in hours)
        row = {"detector": lab, "windows": len(hours), "hours": round(T, 2), "alarms": n,
               "per_day": round(24 * n / T, 2) if T else ""}
        if lab != reference and hours:
            q = significance.quiet(hours, {reference: {w: counts[reference][w] for w in hours},
                                           lab: {w: counts[lab][w] for w in hours}}, rng)
            row.update({"ratio_reference_over_this": q["ratio"], "exact_low": q["ratio_exact_ci"][0],
                        "exact_high": q["ratio_exact_ci"][1], "p_equal_rates": q["p_equal_rates"],
                        "boot_low": q["ratio_boot_ci"][0], "boot_high": q["ratio_boot_ci"][1]})
        rows.append(row)
    return rows


def write(path, rows):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    tmp.rename(path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--draw", required=True, help="the draw's directory (selection.json)")
    ap.add_argument("--run", required=True, help="the run_heldout.py output directory (quiet/LABEL/WINDOW.txt)")
    ap.add_argument("--settings", required=True, help="the settings file of the run (its labels, in order)")
    ap.add_argument("--reference", default="learned")
    ap.add_argument("--groups", default="", help="NAME=SEQ,SEQ;NAME=SEQ,...")
    ap.add_argument("--seed", type=int, default=20261022, help="as tools/significance.py")
    a = ap.parse_args(argv)
    run = Path(a.run)
    labels = [k for k in json.loads(Path(a.settings).read_text()) if not k.startswith("_")]
    if a.reference not in labels:
        sys.exit(f"{a.reference}: not a label of {a.settings}")
    windows = windows_of(a.draw)
    if not windows:
        sys.exit(f"{a.draw}: no windows")
    index = run / "analysis" / "quietstrat_index.csv"
    index.parent.mkdir(parents=True, exist_ok=True)
    write(index, [{"name": w["name"], "scored_hours": f"{w['hours']:.4f}"} for w in windows])
    _, counts = significance.quiet_counts(index, {lab: run / "quiet" / lab for lab in labels})
    per_window = [{"window": w["name"], "sequence": w["sequence"], "start": w["start"], "scored_hours": w["hours"],
                   **{lab: counts[lab][w["name"]] for lab in labels}} for w in windows]
    write(run / "analysis" / "quietstrat_windows.csv", per_window)
    summary = rates(windows, counts, labels, a.reference, random.Random(a.seed))
    write(run / "analysis" / "quietstrat_summary.csv", summary)
    groups = parse_groups(a.groups)
    named = set().union(*groups.values()) if groups else set()
    for s in sorted({w["sequence"] for w in windows}):
        groups[f"set {s:03d}"] = {s}
    if named and {w["sequence"] for w in windows} - named:
        sys.exit(f"--groups leaves out the sets {sorted({w['sequence'] for w in windows} - named)}")
    grouped = []
    for name, seqs in groups.items():
        ws = [w for w in windows if w["sequence"] in seqs]
        for r in rates(ws, counts, labels, a.reference, random.Random(a.seed)):
            grouped.append({"group": name, **r})
    write(run / "analysis" / "quietstrat_groups.csv", grouped)
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    rec = {"note": "exploratory, outside the held-out protocol", "argv": sys.argv,
           "commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                                    cwd=Path(__file__).resolve().parent).stdout.strip(),
           "inputs_sha256": {"selection.json": sha(Path(a.draw) / "selection.json"), "settings": sha(a.settings)},
           "windows": len(windows), "labels": labels, "summary": summary}
    tmp = run / "analysis" / "quietstrat.json.tmp"
    tmp.write_text(json.dumps(rec, indent=1))
    tmp.rename(run / "analysis" / "quietstrat.json")
    for r in summary:
        extra = (f"  ratio {r['ratio_reference_over_this']:.2f} ({r['exact_low']:.2f} to {r['exact_high']:.2f})"
                 if "ratio_reference_over_this" in r else "")
        print(f"{r['detector']:28s} {r['alarms']:5d} alarms  {r['per_day']:7.2f} per day{extra}")
    print(f"-> {run / 'analysis'}")


if __name__ == "__main__":
    main()
