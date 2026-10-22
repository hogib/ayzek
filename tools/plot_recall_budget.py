#!/usr/bin/env python3
"""The paper's Figure 5: recall against the alarm budget on the 728-day MANT record.

Two steps, so the figure can be redrawn without the archive:

    # 1. curves, from the continuous scan (cnn_earthquake runs_mant: window
    #    scores, the per-event table and the report of each detector)
    uv run --with numpy python tools/plot_recall_budget.py curves \\
        --runs ~/Projects/sismokaos/cnn_earthquake/runs_mant \\
        --out ../ayzek_paper/figures/fig5_recall_budget.csv

    # 2. the figure, from that table alone
    uv run --with matplotlib python tools/plot_recall_budget.py plot \\
        --csv ../ayzek_paper/figures/fig5_recall_budget.csv \\
        --out ../ayzek_paper/figures/fig5_recall_budget

The curves follow `sk falsealarm report` exactly, so that they pass through the
numbers of Table 10: a window is background when it overlaps no catalogued
event's guard (10 s before the predicted P minus the window length, to 60 s
after); the threshold for a budget of B a day is the background quantile that
B windows a day exceed; recall is the share of events with SNR >= 3 and a
waveform whose highest score in the guard is above it. `curves` checks itself
against the report's thresholds table at 1, 10 and 100 a day and stops if it
disagrees. The event AUC in the legend is read from the report.
"""
import argparse
import csv
import glob
import re
import sys
from pathlib import Path

# name in runs_mant, window seconds, label, group
ARMS = [
    ("pw", 6.0, "6 s detector", "deployed"),
    ("ps", 6.0, "STA/LTA, causal 2–20 Hz", "deployed"),
    ("6s", 6.0, "6 s detector, station scaling", "offline"),
    ("pnat", 3.4, "P-wave, natural negatives", "offline"),
    ("ponly", 3.4, "P-wave, amplitude-matched", "offline"),
    ("stalta", 6.0, "STA/LTA, classical", "offline"),
]
GUARD_PRE, GUARD_POST, SNR_MIN = 10.0, 60.0, 3.0
# 20 budgets a decade from 0.1 to 100 a day; 0.1, 1, 10 and 100 among them.
BUDGETS = [round(10 ** (k / 20), 6) for k in range(-20, 41)]


def curves(a):
    import numpy as np

    runs = Path(a.runs).expanduser()
    rows = []
    for name, win, label, group in ARMS:
        files = sorted(glob.glob(str(runs / "scores_mant" / name / "*.npz")))
        if not files:
            sys.exit(f"no scores in {runs}/scores_mant/{name}")
        t = np.concatenate([np.load(f)["t"] for f in files])
        p = np.concatenate([np.load(f)["p"] for f in files])
        order = np.argsort(t)
        t, p = t[order], p[order]
        step = float(np.median(np.diff(t[:100000])))
        days = len(t) * step / 86400.0

        pe, best = [], []
        with open(runs / f"final_{name}_events.csv") as fh:
            for r in csv.DictReader(fh):
                pe.append(float(r["p_epoch"]))
                if r["best_prob"] and r["snr"] and float(r["snr"]) >= SNR_MIN:
                    best.append(float(r["best_prob"]))
        pe, best = np.asarray(pe), np.asarray(best)
        explained = np.zeros(len(t), bool)
        lo = np.searchsorted(t, pe - GUARD_PRE - win)
        hi = np.searchsorted(t, pe + GUARD_POST, side="right")
        for i, j in zip(lo, hi):
            explained[i:j] = True
        bg = p[~explained]

        m = re.search(r"^\s+3\s+[\d,]+\s+([\d.]+)", (runs / f"report_{name}.log").read_text(), re.M)
        auc = float(m.group(1)) if m else float("nan")

        ref = {float(r["target_per_day"]): float(r["recall_snr3"])
               for r in csv.DictReader(open(runs / f"final_{name}_thresholds.csv"))}
        for b in BUDGETS:
            want = b * days
            if want >= len(bg):
                continue
            thr = float(np.quantile(bg, 1.0 - want / len(bg)))
            rec = float((best > thr).mean())
            if b in ref and abs(rec - ref[b]) > 1e-9:
                sys.exit(f"{name} at {b}/day: {rec:.6f}, the report says {ref[b]:.6f}")
            rows.append({"detector": name, "label": label, "group": group, "event_auc": f"{auc:.4f}",
                         "budget_per_day": b, "threshold": f"{thr:.6g}", "recall_snr3": f"{rec:.6f}",
                         "events": len(best)})
        print(f"{name:7s} {len(best):,} events, {days:.1f} days, AUC {auc:.4f}: "
              + ", ".join(f"{b:g}/day {ref[b]:.3f}" for b in (1.0, 10.0, 100.0) if b in ref)
              + "  (matches the report)")
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"-> {a.out}")


def plot(a):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e5e7eb"
    # The detectors keep the hues of the paper's other figures: learned #9a3412,
    # STA/LTA #1d4ed8 (validated pair); the offline counterparts take the same
    # hue dashed, the P-wave detectors two greys.
    STYLE = {
        "pw": dict(color="#9a3412", lw=2.2, ls="-"),
        "ps": dict(color="#1d4ed8", lw=2.2, ls="-"),
        "6s": dict(color="#9a3412", lw=1.3, ls=(0, (5, 2))),
        "stalta": dict(color="#1d4ed8", lw=1.3, ls=(0, (5, 2))),
        "pnat": dict(color="#6b7280", lw=1.3, ls="-"),
        "ponly": dict(color="#9ca3af", lw=1.3, ls=(0, (1.5, 1.5))),
    }
    plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "pdf.fonttype": 42})
    rows = list(csv.DictReader(open(a.csv)))
    fig, ax = plt.subplots(figsize=(6.0, 3.9))
    fig.subplots_adjust(left=0.1, right=0.98, top=0.97, bottom=0.13)
    handles = {}
    for name, _, _, _ in ARMS:
        r = [x for x in rows if x["detector"] == name]
        if not r:
            continue
        x = [float(v["budget_per_day"]) for v in r]
        y = [float(v["recall_snr3"]) for v in r]
        (h,) = ax.plot(x, y, solid_capstyle="round", **STYLE[name])
        marks = [(xi, yi) for xi, yi in zip(x, y) if xi in (1.0, 10.0, 100.0)]
        ax.plot(*zip(*marks), "o", ms=4.5 if r[0]["group"] == "deployed" else 3.5,
                color=STYLE[name]["color"], mec="white", mew=0.8, zorder=4, clip_on=False)
        handles[name] = (h, f"{r[0]['label']}  (AUC {float(r[0]['event_auc']):.3f})", r[0]["group"])
    ax.axvline(1.0, color=MUTED, lw=0.6, ls=":", zorder=0)
    ax.set_xscale("log")
    ax.set_xlim(0.1, 100)
    ax.set_xticks([0.1, 1, 10, 100])
    ax.set_xticklabels(["0.1", "1", "10", "100"])
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("alarm budget (false alarms per day at MANT)", color=MUTED)
    ax.set_ylabel("recall, events with SNR ≥ 3", color=MUTED)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("right", "top"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
        ax.spines[side].set_linewidth(0.6)
    ax.tick_params(colors=MUTED, length=3, width=0.6)
    blank = plt.Line2D([], [], lw=0)
    hs, ls = [], []
    for group, title in (("deployed", "Run by the pipeline"), ("offline", "Offline")):
        hs.append(blank)
        ls.append(r"$\bf{" + title.replace(" ", r"\ ") + "}$")
        for h, lab, g in handles.values():
            if g == group:
                hs.append(h)
                ls.append(lab)
    leg = ax.legend(hs, ls, loc="upper left", frameon=False, fontsize=7, handlelength=2.6)
    for text in leg.get_texts():
        text.set_color(INK)
    out = Path(a.out)
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"-> {out}.pdf, {out}.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("curves", help="recall at each budget, from the scan")
    c.add_argument("--runs", required=True, help="cnn_earthquake runs_mant directory")
    c.add_argument("--out", required=True, help="CSV")
    q = sub.add_parser("plot", help="the figure, from the curves CSV")
    q.add_argument("--csv", required=True)
    q.add_argument("--out", required=True, help="path without extension")
    a = ap.parse_args()
    {"curves": curves, "plot": plot}[a.cmd](a)


if __name__ == "__main__":
    main()
