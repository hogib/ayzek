#!/usr/bin/env python3
"""Recall and event AUC on the 728-day MANT record with the S wave kept out.

The paper's Table 12 scores an event by the highest window score in its guard,
10 s before the predicted P (minus the window length) to 60 s after it. Windows
that late contain the S wave and coda of events out to ~500 km, so splitting the
events by distance does not isolate detection on P energy (referee round 5,
M4). Here an event is scored only from windows that END at least MARGIN seconds
before its predicted S (IASP91, from the per-event table), and only events whose
S-P interval leaves every detector such a window are counted: S-P >= 6 s (the
longest window, and the scan step of the 6 s detectors) + MARGIN.

Thresholds are unchanged: they come from the background, which does not depend
on how events are scored. The tool reproduces the report's thresholds and full-
guard recalls first and stops if they disagree.

The event AUC is recomputed the same way: each event's P-only maximum against
the maximum over one random stretch of record of the same length, away from
every catalogued event (seeded).

    uv run --with numpy --with scipy python tools/s_wave_control.py \\
        --runs ~/Projects/sismokaos/cnn_earthquake/runs_mant \\
        --margins 0 1 2 --out data/runs/s_wave_control.csv
"""
import argparse
import csv
import glob
import re
import sys
from pathlib import Path

# name in runs_mant, window seconds, label (as in plot_recall_budget.py)
ARMS = [
    ("pw", 6.0, "pipeline 6 s detector"),
    ("ps", 6.0, "pipeline STA/LTA"),
    ("6s", 6.0, "6 s detector, station scaling"),
    ("pnat", 3.4, "P-wave, natural negatives"),
    ("ponly", 3.4, "P-wave, amplitude-matched"),
    ("stalta", 6.0, "STA/LTA, classical"),
]
GUARD_PRE, GUARD_POST, SNR_MIN = 10.0, 60.0, 3.0
BUDGETS = (100.0, 10.0, 1.0)
SP_FLOOR = 6.0  # the longest window; every arm then has a window ending in (P, S - margin]


def auc(pos, neg):
    """Mann-Whitney AUC, ties counted half."""
    import numpy as np
    from scipy.stats import rankdata

    r = rankdata(np.r_[pos, neg])
    n1, n0 = len(pos), len(neg)
    return (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def main():
    import numpy as np

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", required=True, help="cnn_earthquake runs_mant")
    ap.add_argument("--margins", type=float, nargs="+", default=[1.0],
                    help="seconds a window must end before the predicted S")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    runs = Path(a.runs).expanduser()

    rows = []
    for name, win, label in ARMS:
        files = sorted(glob.glob(str(runs / "scores_mant" / name / "*.npz")))
        if not files:
            sys.exit(f"no scores in {runs}/scores_mant/{name}")
        t = np.concatenate([np.load(f)["t"] for f in files])  # window starts
        p = np.concatenate([np.load(f)["p"] for f in files])
        order = np.argsort(t)
        t, p = t[order], p[order]
        step = float(np.median(np.diff(t[:100000])))
        days = len(t) * step / 86400.0

        ev = list(csv.DictReader(open(runs / f"final_{name}_events.csv")))
        pe_all = np.array([float(r["p_epoch"]) for r in ev])
        explained = np.zeros(len(t), bool)
        for i, j in zip(np.searchsorted(t, pe_all - GUARD_PRE - win),
                        np.searchsorted(t, pe_all + GUARD_POST, side="right")):
            explained[i:j] = True
        bg = p[~explained]

        ref = {float(r["target_per_day"]): (float(r["threshold"]), float(r["recall_snr3"]))
               for r in csv.DictReader(open(runs / f"final_{name}_thresholds.csv"))}
        thr = {}
        for b in BUDGETS:
            thr[b] = float(np.quantile(bg, 1.0 - b * days / len(bg)))
            if abs(thr[b] - ref[b][0]) > 1e-6:
                sys.exit(f"{name} {b}/day: threshold {thr[b]:.6f}, the report says {ref[b][0]:.6f}")

        # Events with a waveform (SNR >= 3) and a score, as in Table 12.
        sel = [r for r in ev if r["best_prob"] and r["snr"] and float(r["snr"]) >= SNR_MIN]
        best = np.array([float(r["best_prob"]) for r in sel])
        pe = np.array([float(r["p_epoch"]) for r in sel])
        se = np.array([float(r["s_epoch"]) for r in sel])
        dist = np.array([float(r["dist"]) for r in sel])
        for b in BUDGETS:
            if abs((best > thr[b]).mean() - ref[b][1]) > 1e-9:
                sys.exit(f"{name} {b}/day: recall {(best > thr[b]).mean():.6f}, "
                         f"the report says {ref[b][1]:.6f}")

        # The random full-length guards of the report, for the full-guard AUC on a subset.
        m = re.search(r"^\s+3\s+[\d,]+\s+([\d.]+)", (runs / f"report_{name}.log").read_text(), re.M)
        auc_report = float(m.group(1)) if m else float("nan")

        pe_sorted = np.sort(pe_all)

        def random_max(lengths, rng):
            """Max score over a random stretch of each length, clear of every event."""
            out = np.empty(len(lengths))
            for k, L in enumerate(lengths):
                while True:
                    c = t[rng.integers(len(t))]
                    lo, hi = c - GUARD_POST - 60, c + L + GUARD_PRE + win + 60
                    if np.searchsorted(pe_sorted, hi) - np.searchsorted(pe_sorted, lo) == 0:
                        i, j = np.searchsorted(t, c), np.searchsorted(t, c + L, side="right")
                        if j > i:
                            out[k] = p[i:j].max()
                            break
            return out

        for margin in a.margins:
            keep = (se - pe) >= SP_FLOOR + margin
            lo = np.searchsorted(t, pe[keep] - GUARD_PRE - win)
            hi = np.searchsorted(t, se[keep] - margin - win, side="right")  # start + win <= S - margin
            ponly = np.array([p[i:j].max() if j > i else np.nan for i, j in zip(lo, hi)])
            ok = ~np.isnan(ponly)  # a gap can still leave an event without a window
            full = best[keep][ok]
            ponly = ponly[ok]
            far = dist[keep][ok] >= 110.0

            rng = np.random.default_rng(a.seed)
            lengths = (se[keep][ok] - margin - win) - (pe[keep][ok] - GUARD_PRE - win)
            rnd_p = random_max(lengths, rng)
            rnd_full = random_max(np.full(ok.sum(), GUARD_PRE + win + GUARD_POST), rng)

            row = {"detector": name, "label": label, "margin_s": margin,
                   "events": int(ok.sum()), "events_ge_110km": int(far.sum()),
                   "auc_report": f"{auc_report:.4f}",
                   "auc_full_subset": f"{auc(full, rnd_full):.4f}",
                   "auc_p_only": f"{auc(ponly, rnd_p):.4f}"}
            for b in BUDGETS:
                row[f"recall_full_{b:g}"] = f"{(full > thr[b]).mean():.4f}"
                row[f"recall_p_only_{b:g}"] = f"{(ponly > thr[b]).mean():.4f}"
                row[f"recall_p_only_ge110_{b:g}"] = f"{(ponly[far] > thr[b]).mean():.4f}"
            rows.append(row)
            print(f"{name:7s} margin {margin:g} s: {ok.sum():,} events (S-P >= {SP_FLOOR + margin:g} s); "
                  f"AUC full {row['auc_full_subset']} P-only {row['auc_p_only']}; "
                  f"1/day full {row['recall_full_1']} P-only {row['recall_p_only_1']}")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
