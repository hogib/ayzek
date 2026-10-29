"""Select trigger settings from tuning sweeps by the paper's two rules.

    python3 tools/select_setting.py --sweeps 'data/runs/sb/eqtransformer_*_causal' \\
        --datasets demo,marmara --unmatched 26 --corrected 48.7

Settings are the trigger options (--threshold, --trigger-windows,
--instant-threshold), FilterPicker's (--fp-windows, --fp-thresholds) or
STA/LTA's (--stalta-band, --sta, --lta, --stalta-on, --stalta-3c).
Each directory holds one tools/sweep_trigger.py table per tuning replay
(sweep_DATASET.csv), for one detector or one set of weights. Settings are summed
over the replays, as §4.2 of the paper does, and the rules applied over every
directory and setting together:

- equal unmatched alarms: the most chance-corrected detections among the
  settings with at most --unmatched unmatched alarms; ties to fewer unmatched
  alarms, then the higher threshold;
- equal detections: the fewest unmatched alarms among the settings with at
  least --corrected chance-corrected detections; ties to more detections;
- equal quiet alarms (with --quiet-alarms and --quiet-table): the most
  chance-corrected detections among the settings that raise at most
  --quiet-alarms alarms on the quiet windows (tools/quiet_sweep.py tables, one
  per source directory); ties to fewer quiet alarms, then fewer unmatched
  alarms.

Standard library only.
"""
import argparse
import csv
import glob
from pathlib import Path

KEYS = ("threshold", "trigger-windows", "instant-threshold", "fp-windows", "fp-thresholds",
        "stalta-band", "sta", "lta", "stalta-on", "stalta-3c")


def threshold(setting):
    """The tie-break of the first rule: the trigger threshold, FilterPicker's
    (threshold1, threshold2), or STA/LTA's trigger ratio."""
    s = dict(setting)
    if "fp-thresholds" in s:
        return tuple(float(x) for x in s["fp-thresholds"].split(","))
    if "stalta-on" in s:
        return (float(s["stalta-on"]),)
    return (float(s["threshold"]),)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sweeps", required=True, help="glob of directories")
    ap.add_argument("--datasets", default="demo,marmara")
    ap.add_argument("--unmatched", type=float, required=True)
    ap.add_argument("--corrected", type=float, required=True)
    ap.add_argument("--quiet-alarms", type=float,
                    help="the learned detector's alarms on the quiet windows (third rule)")
    ap.add_argument("--quiet-table", action="append", default=[], metavar="SOURCE=CSV",
                    help="tools/quiet_sweep.py table of the source directory SOURCE")
    a = ap.parse_args()

    rows = []
    for d in sorted(glob.glob(a.sweeps)):
        tables = {}
        for ds in a.datasets.split(","):
            f = Path(d) / f"sweep_{ds}.csv"
            if not f.exists():
                raise SystemExit(f"missing {f}")
            tables[ds] = {tuple((k, r[k]) for k in KEYS if k in r): r for r in csv.DictReader(open(f))}
        for setting in set.intersection(*(set(t) for t in tables.values())):
            rs = [tables[ds][setting] for ds in tables]
            rows.append({
                "source": Path(d).name, "setting": setting, "quiet": None,
                "detected": sum(int(r["detected"]) for r in rs),
                "catalogue": sum(int(r["catalogue"]) for r in rs),
                "corrected": sum(float(r["corrected"]) for r in rs),
                "unmatched": sum(int(r["false"]) for r in rs),
            })
    if not rows:
        raise SystemExit("no settings")
    for spec in a.quiet_table:
        src, path = spec.split("=", 1)
        q = {tuple((k, r[k]) for k in KEYS if k in r): int(r["alarms"]) for r in csv.DictReader(open(path))}
        for r in rows:
            if r["source"] == src:
                r["quiet"] = q.get(r["setting"])

    def show(label, r):
        if r is None:
            print(f"{label}: none satisfies the rule")
            return
        s = " ".join(f"{k}={v}" for k, v in r["setting"])
        print(f"{label}: {r['source']}  {s}  detected {r['detected']}/{r['catalogue']}, "
              f"corrected {r['corrected']:.1f}, unmatched {r['unmatched']}")

    fa = [r for r in rows if r["unmatched"] <= a.unmatched]
    best_fa = max(fa, key=lambda r: (r["corrected"], -r["unmatched"], threshold(r["setting"]))) if fa else None
    dt = [r for r in rows if r["corrected"] >= a.corrected]
    best_dt = min(dt, key=lambda r: (r["unmatched"], -r["corrected"])) if dt else None
    show("equal unmatched", best_fa)
    show("equal detections", best_dt)
    if a.quiet_alarms is not None:
        qa = [r for r in rows if r["quiet"] is not None and r["quiet"] <= a.quiet_alarms]
        best_q = max(qa, key=lambda r: (r["corrected"], -r["quiet"], -r["unmatched"])) if qa else None
        show("equal quiet alarms", best_q)
        if best_q:
            print(f"  (quiet windows: {best_q['quiet']} alarms, at most {a.quiet_alarms:g})")
    print(f"\n{len(rows)} settings over {len({r['source'] for r in rows})} sources; best per source at <= {a.unmatched:g} unmatched:")
    for src in sorted({r["source"] for r in rows}):
        c = [r for r in fa if r["source"] == src]
        if c:
            show(f"  {src}", max(c, key=lambda r: (r["corrected"], -r["unmatched"])))
        else:
            print(f"  {src}: none at <= {a.unmatched:g} unmatched")


if __name__ == "__main__":
    main()
