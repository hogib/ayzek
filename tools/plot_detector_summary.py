#!/usr/bin/env python3
"""The paper's summary figure: every detector's held-out detections against its
alarm rate on the quiet windows, at the settings selected on the tuning data.

Two steps, so the figure can be redrawn without the runs:

    # 1. the table, from the paired statistics and the quiet-window replays
    python3 tools/plot_detector_summary.py data \\
        --significance data/runs/r5/base.json --out ../ayzek_paper/figures/fig_detectors.csv

    # 2. the figure, from that table alone
    uv run --with matplotlib python tools/plot_detector_summary.py plot \\
        --csv ../ayzek_paper/figures/fig_detectors.csv --out ../ayzek_paper/figures/fig3_detectors

Detections are chance-corrected (tools/significance.py, groups `sindirgi61` and
`simav26`); quiet-window alarms are counted as tools/significance.py counts them,
over the 102.8 h of the 23 AFAD quiet windows. The equal-detections STA/LTA
setting was replayed on the quiet windows only within its grid
(data/runs/r4/quiet_stalta.csv, R13), which gives its count.

Design (dataviz skill): one panel per held-out sequence on a shared log x axis;
colour for the three detector families that matter to the argument (learned
#9a3412 and STA/LTA #1d4ed8, the paper's validated pair, and FilterPicker
#0f7a5a; all-pairs CVD dE >= 8.5 on white), the published models in neutral ink;
marker shape names the published model, a hollow marker the look-ahead
(offline) modes; every point carries a direct label, so no identity rests on
colour alone.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

# label in the outcome files, quiet-window run directory (None: from the grid),
# family, mode, rule, short label on the figure
DETECTORS = [
    ("learned", "data/runs/quiet100", "learned", "causal", "", "learned"),
    ("stalta_eqfa", "data/runs/quiet100_stalta", "stalta", "causal", "equal unmatched", "STA/LTA, eq. unm."),
    ("stalta_eqdet", None, "stalta", "causal", "equal detections", "STA/LTA, eq. det."),
    ("fp_fa", "data/runs/fp/eval/fp_fa/quiet", "filterpicker", "causal", "equal unmatched", "FilterPicker, eq. unm."),
    ("fp_det", "data/runs/fp/eval/fp_det/quiet", "filterpicker", "causal", "equal detections", "FilterPicker, eq. det."),
    ("pn_causal", "data/runs/sb/eval/pn_causal/quiet", "phasenet", "causal", "both", "PhaseNet"),
    ("eqt_causal_fa", "data/runs/sb/eval/eqt_causal_fa/quiet", "eqtransformer", "causal", "equal unmatched", "EQT, eq. unm."),
    ("eqt_causal_det", "data/runs/sb/eval/eqt_causal_det/quiet", "eqtransformer", "causal", "equal detections", "EQT, eq. det."),
    ("pn_look_fa", "data/runs/sb/eval/pn_look_fa/quiet", "phasenet", "look-ahead", "equal unmatched", "PhaseNet LA, eq. unm."),
    ("pn_look_det", "data/runs/sb/eval/pn_look_det/quiet", "phasenet", "look-ahead", "equal detections", "PhaseNet LA, eq. det."),
    ("eqt_look", "data/runs/sb/eval/eqt_look/quiet", "eqtransformer", "look-ahead", "both", "EQT LA"),
]
STALTA_EQDET = ('"2,10"', "0.5", "10", "6", "off")  # its row in quiet_stalta.csv


def data(a):
    from significance import quiet_counts

    summ = json.load(open(a.significance))["summary"]
    hours, counts = quiet_counts(a.quiet_index, {d[0]: d[1] for d in DETECTORS if d[1]})
    total_h = sum(hours.values())
    grid = {}
    for line in open(a.quiet_grid).read().splitlines()[1:]:
        f = line.replace('"', "").split(",")
        grid[('"' + f[0] + "," + f[1] + '"', f[2], f[3], f[4], f[5])] = int(f[6])
    rows = []
    for label, _, family, mode, rule, short in DETECTORS:
        n = sum(counts[label].values()) if label in counts else grid[STALTA_EQDET]
        row = {"detector": label, "family": family, "mode": mode, "rule": rule, "short": short,
               "quiet_alarms": n, "quiet_per_day": f"{n / total_h * 24:.2f}"}
        for g, col in (("sindirgi61", "heldout_corrected"), ("simav26", "simav_corrected")):
            key = f"{g}: learned corrected" if label == "learned" else f"{g}: {label} corrected"
            row[col] = summ[key]
        rows.append(row)
        print(f"{label:15s} held out {row['heldout_corrected']:5} Simav {row['simav_corrected']:5} "
              f"quiet {n} ({row['quiet_per_day']}/day)")
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
    COLOR = {"learned": "#9a3412", "stalta": "#1d4ed8", "filterpicker": "#0f7a5a",
             "phasenet": "#3f3f3c", "eqtransformer": "#3f3f3c"}
    MARKER = {"learned": "*", "stalta": "s", "filterpicker": "D", "phasenet": "o", "eqtransformer": "^"}
    SIZE = {"learned": 150, "stalta": 46, "filterpicker": 40, "phasenet": 50, "eqtransformer": 56}
    plt.rcParams.update({"font.size": 7.5, "font.family": "DejaVu Sans", "pdf.fonttype": 42})
    rows = list(csv.DictReader(open(a.csv)))
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.5), sharex=True)
    fig.subplots_adjust(left=0.075, right=0.99, top=0.9, bottom=0.14, wspace=0.18)
    # label placement per panel: (dx, dy) in points and alignment, set against collisions
    R, L, C = "left", "right", "center"
    OFF = {
        "heldout_corrected": {"learned": (8, 4, R), "stalta_eqfa": (8, -5, R), "stalta_eqdet": (8, 4, R),
                              "fp_fa": (8, -5, R), "fp_det": (8, -5, R), "pn_causal": (6, -9, R),
                              "eqt_causal_fa": (8, 0, R), "eqt_causal_det": (8, -4, R),
                              "pn_look_fa": (8, 3, R), "pn_look_det": (0, 9, C), "eqt_look": (-8, 0, L)},
        "simav_corrected": {"learned": (8, -4, R), "stalta_eqfa": (8, -5, R), "stalta_eqdet": (8, 4, R),
                            "fp_fa": (8, -5, R), "fp_det": (-8, -7, L), "pn_causal": (-8, 2, L),
                            "eqt_causal_fa": (8, 0, R), "eqt_causal_det": (8, -5, R),
                            "pn_look_fa": (8, -7, R), "pn_look_det": (0, 9, C), "eqt_look": (-8, 0, L)},
    }
    for ax, col, title in ((axes[0], "heldout_corrected", "(a)  Sındırgı Mw 6.1, 66 events"),
                           (axes[1], "simav_corrected", "(b)  Simav Mw 4.8, 56 events")):
        for r in rows:
            x, y, fam = float(r["quiet_per_day"]), float(r[col]), r["family"]
            hollow = r["mode"] == "look-ahead"
            ax.scatter([x], [y], s=SIZE[fam], marker=MARKER[fam], zorder=3,
                       facecolor="white" if hollow else COLOR[fam], edgecolor=COLOR[fam],
                       linewidth=1.3 if hollow else 0.8)
            dx, dy, ha = OFF[col][r["detector"]]
            ax.annotate(r["short"], (x, y), xytext=(dx, dy), textcoords="offset points",
                        fontsize=6.5, color=INK if fam != "learned" else COLOR[fam], ha=ha, va="center")
        ax.set_xscale("log")
        ax.set_xlim(0.8, 2000)
        ax.set_ylim(5, 48)
        ax.set_xticks([1, 3, 10, 30, 100, 300, 1000])
        ax.set_xticklabels(["1", "3", "10", "30", "100", "300", "1000"])
        ax.minorticks_off()
        ax.grid(color=GRID, lw=0.6, zorder=0)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(MUTED)
        ax.tick_params(colors=MUTED)
        ax.set_title(title, loc="left", fontweight="bold", fontsize=8)
        ax.set_xlabel("alarms per day on the 102.8 h of quiet data", color=MUTED)
    axes[0].set_ylabel("chance-corrected detections", color=MUTED)
    fig.text(0.075, 0.965, "filled: causal (real time); hollow: look-ahead (offline only)",
             color=MUTED, fontsize=6.5)
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=300)
    print(f"-> {a.out}.pdf, {a.out}.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("data")
    d.add_argument("--significance", required=True, help="tools/significance.py --json output with "
                   "groups sindirgi61 and simav26 and every detector")
    d.add_argument("--quiet-index", default="data/quiet100/quiet_index.csv")
    d.add_argument("--quiet-grid", default="data/runs/r4/quiet_stalta.csv")
    d.add_argument("--out", required=True)
    p = sub.add_parser("plot")
    p.add_argument("--csv", required=True)
    p.add_argument("--out", required=True, help="path without extension")
    a = ap.parse_args()
    (data if a.cmd == "data" else plot)(a)


if __name__ == "__main__":
    main()
