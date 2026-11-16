#!/usr/bin/env python3
"""Trade-off curves on the three test networks: detections against alarms as one
threshold of each detector varies (exploratory, outside the held-out protocol).

    # 1. one row per network and grid setting, from the grid's summary table
    python3 tools/plot_tradeoff.py data --summary data/runs/heldout_grid/summary.csv \\
        --out ../ayzek_paper/figures/fig_tradeoff.csv
    # 2. the figure, from that table alone
    uv run --with matplotlib python tools/plot_tradeoff.py plot \\
        --csv ../ayzek_paper/figures/fig_tradeoff.csv --out ../ayzek_paper/figures/fig_tradeoff

One column per network. Top row: chance-corrected detections against unmatched alarms
over the network's sequences. Bottom row: the same detections against alarms per day on
its quiet windows. One line per detector, through its grid settings in threshold order;
the frozen setting of the tests is the large ringed marker. Both alarm axes are
logarithmic above 1 and linear below, so that settings without an alarm can be drawn.

Colours as in the paper's other figures: learned #9a3412, STA/LTA #1d4ed8, FilterPicker
#0f7a5a, the published models in neutral ink. The two published models share the ink and
differ by dash and marker, and every detector has its own marker, so no pair is told
apart by colour alone.
"""
import argparse
import csv
from pathlib import Path

# family: (label, the frozen value of the varied threshold, what is varied)
FAMILY = {
    "learned": ("learned detector", "0.8", "probability threshold of the eight-window rule"),
    "stalta": ("STA/LTA", "8", "trigger ratio"),
    "fp": ("FilterPicker", "20", "first threshold"),
    "phasenet": ("PhaseNet, causal", "0.8", "probability threshold"),
    "eqt": ("EQTransformer, causal", "0.7", "probability threshold"),
}
NETS = {"greece_now": "Greece", "nz_now": "New Zealand", "ko_now": "KO"}


def data(a):
    rows = []
    for r in csv.DictReader(open(a.summary)):
        _, fam, val = r["setting"].split("_", 2)
        if fam not in FAMILY:
            raise SystemExit(f"{r['setting']}: unknown detector family")
        val = val.replace("p", ".")
        hours = float(r["quiet_hours"])
        rows.append({"network": r["network"], "detector": fam, "value": val,
                     "frozen": int(val == FAMILY[fam][1]), "events": int(r["events"]),
                     "corrected": float(r["corrected"]), "unmatched_alarms": int(r["unmatched_alarms"]),
                     "quiet_alarms": int(r["quiet_alarms"]), "quiet_hours": hours,
                     "quiet_per_day": round(int(r["quiet_alarms"]) / hours * 24.0, 2)})
    if not rows:
        raise SystemExit(f"{a.summary}: no rows")
    for net in NETS:
        for fam in FAMILY:
            n = sum(1 for r in rows if r["network"] == net and r["detector"] == fam and r["frozen"])
            if n != 1:
                raise SystemExit(f"{net} {fam}: {n} frozen settings, expected 1")
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} rows -> {a.out}")


def plot(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e5e7eb"
    STYLE = {
        "learned": dict(color="#9a3412", marker="o", ls="-", lw=1.8),
        "stalta": dict(color="#1d4ed8", marker="s", ls="-", lw=1.8),
        "fp": dict(color="#0f7a5a", marker="^", ls="-", lw=1.4),
        "phasenet": dict(color="#3f3f3c", marker="D", ls="-", lw=1.2),
        "eqt": dict(color="#3f3f3c", marker="v", ls=(0, (4, 2)), lw=1.2),
    }
    plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "pdf.fonttype": 42})
    rows = list(csv.DictReader(open(a.csv)))
    for r in rows:
        for k in ("corrected", "quiet_per_day"):
            r[k] = float(r[k])
        for k in ("unmatched_alarms", "frozen", "events"):
            r[k] = int(r[k])

    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.0), sharey="col")
    fig.subplots_adjust(left=0.075, right=0.965, top=0.86, bottom=0.09, wspace=0.2, hspace=0.34)
    for j, (net, name) in enumerate(NETS.items()):
        rs = [r for r in rows if r["network"] == net]
        events = rs[0]["events"]
        for i, (key, xlabel) in enumerate((("unmatched_alarms", "unmatched alarms"),
                                           ("quiet_per_day", "quiet-window alarms per day"))):
            ax = axes[i][j]
            for fam in ("eqt", "phasenet", "fp", "stalta", "learned"):      # learned drawn last, on top
                fs = [r for r in rs if r["detector"] == fam]                 # in grid (threshold) order
                st = STYLE[fam]
                ax.plot([r[key] for r in fs], [r["corrected"] for r in fs], color=st["color"], ls=st["ls"],
                        lw=st["lw"], marker=st["marker"], ms=2.6, mec=st["color"], mfc=st["color"], zorder=3)
                for r in fs:
                    if r["frozen"]:
                        ax.plot(r[key], r["corrected"], marker=st["marker"], ms=7.5, color=st["color"],
                                mec="white", mew=1.1, zorder=5, clip_on=False)
            ax.set_xscale("symlog", linthresh=1, linscale=0.5)
            right = max(float(r[key]) for r in rs) * 1.6
            ticks = [t for t in (0, 1, 10, 100, 1000, 10000) if t <= right]
            ax.set_xlim(-0.15, right)
            ax.set_xticks(ticks)
            ax.set_xticklabels([f"{t:,}".replace(",", " ") for t in ticks])
            ax.set_ylim(0, events)
            ax.set_xlabel(xlabel, color=MUTED)
            if j == 0:
                ax.set_ylabel("chance-corrected detections", color=MUTED)
            if i == 0:
                ax.set_title(f"{name} ({events} events)", loc="left", fontweight="bold", color=INK)
            ax.grid(color=GRID, lw=0.6)
            ax.set_axisbelow(True)
            for side in ("right", "top"):
                ax.spines[side].set_visible(False)
            for side in ("left", "bottom"):
                ax.spines[side].set_color(MUTED)
                ax.spines[side].set_linewidth(0.6)
            ax.tick_params(colors=MUTED, length=3, width=0.6, labelleft=True)
    order = ("learned", "stalta", "fp", "phasenet", "eqt")
    h = [plt.Line2D([], [], color=STYLE[f]["color"], ls=STYLE[f]["ls"], lw=STYLE[f]["lw"],
                    marker=STYLE[f]["marker"], ms=4) for f in order]
    h.append(plt.Line2D([], [], ls="", marker="o", ms=7.5, color=MUTED, mec="white", mew=1.1))
    leg = fig.legend(h, [FAMILY[f][0] for f in order] + ["setting of the tests (frozen)"], loc="upper center",
                     ncols=3, frameon=False, fontsize=7.5, handlelength=2.6, columnspacing=1.8,
                     bbox_to_anchor=(0.5, 1.0))
    for t in leg.get_texts():
        t.set_color(INK)
    out = Path(a.out)
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"-> {out}.pdf, {out}.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("data", help="the figure's table from the grid summary")
    d.add_argument("--summary", required=True)
    d.add_argument("--out", required=True)
    d.set_defaults(fn=data)
    p = sub.add_parser("plot", help="the figure from its table")
    p.add_argument("--csv", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(fn=plot)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
