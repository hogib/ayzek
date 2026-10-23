#!/usr/bin/env python3
"""Every trigger setting of both detectors: held-out operating curve and quiet-window alarms.

    # 1. one row per setting, from the grids (data/runs/r4: tools/sweep_trigger.py
    #    on the three replays, tools/quiet_sweep.py on the quiet windows)
    python3 tools/plot_operating_curves.py data --runs data/runs/r4 \\
        --out ../ayzek_paper/figures/fig7_operating_curves.csv
    # 2. the figure, from that table alone
    uv run --with matplotlib python tools/plot_operating_curves.py plot \\
        --csv ../ayzek_paper/figures/fig7_operating_curves.csv \\
        --out ../ayzek_paper/figures/fig7_operating_curves

(a) Held-out replay: chance-corrected detections against unmatched alarms for
all 60 settings of the learned detector and 216 of STA/LTA (off level 1.5),
with each detector's upper envelope. Post hoc: no setting was chosen on it.
(b) Quiet windows: alarms per day against chance-corrected detections on the
tuning replays, with each detector's lower envelope (the fewest alarms among
settings with at least that many detections). Settings with no alarm are drawn
at the bottom edge. The selected settings are marked in both panels.
"""
import argparse
import csv
from pathlib import Path

LEARNED = {"adopted": ("0.8", "8", "0.9")}
STALTA = {"equal unmatched": ("2,20", "1", "30", "8", "off"),
          "equal detections": ("2,10", "0.5", "10", "6", "off"),
          "rule at 26 unmatched": ("1,10", "0.5", "10", "8", "off")}
LKEYS = ("threshold", "trigger-windows", "instant-threshold")
SKEYS = ("stalta-band", "sta", "lta", "stalta-on", "stalta-3c")


def data(a):
    runs = Path(a.runs)
    rows = []
    for det, keys, marks in (("learned", LKEYS, LEARNED), ("stalta", SKEYS, STALTA)):
        ev = {}
        for n in ("demo", "marmara", "sindirgi61"):
            for r in csv.DictReader(open(runs / f"{det}_{n}.csv")):
                ev.setdefault(tuple(r[k] for k in keys), {})[n] = (float(r["corrected"]), int(r["false"]))
        quiet = {tuple(r[k] for k in keys): float(r["per_day"])
                 for r in csv.DictReader(open(runs / f"quiet_{det}.csv"))}
        names = {v: k for k, v in marks.items()}
        for k, v in ev.items():
            rows.append({"detector": det, "setting": " ".join(k), "selected": names.get(k, ""),
                         "heldout_corrected": v["sindirgi61"][0], "heldout_unmatched": v["sindirgi61"][1],
                         "tuning_corrected": round(v["demo"][0] + v["marmara"][0], 1),
                         "tuning_unmatched": v["demo"][1] + v["marmara"][1], "quiet_per_day": quiet[k]})
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} settings -> {a.out}")


def plot(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e5e7eb"
    COLOR = {"learned": "#9a3412", "stalta": "#1d4ed8"}   # the paper's validated pair
    LABEL = {"learned": "learned detector", "stalta": "STA/LTA"}
    MARK = {"adopted": "o", "equal unmatched": "s", "equal detections": "^", "rule at 26 unmatched": "D"}
    plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "pdf.fonttype": 42})
    rows = list(csv.DictReader(open(a.csv)))
    for r in rows:
        for k in ("heldout_corrected", "tuning_corrected", "quiet_per_day"):
            r[k] = float(r[k])
        for k in ("heldout_unmatched", "tuning_unmatched"):
            r[k] = int(r[k])
    FLOOR = 0.3   # settings without an alarm on the quiet windows

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.8))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.80, bottom=0.13, wspace=0.27)
    for det in ("stalta", "learned"):
        rs = [r for r in rows if r["detector"] == det]
        c = COLOR[det]
        ax1.scatter([r["heldout_unmatched"] for r in rs], [r["heldout_corrected"] for r in rs],
                    s=7, color=c, alpha=0.35, linewidths=0)
        env, best = [], -1
        for r in sorted(rs, key=lambda r: (r["heldout_unmatched"], -r["heldout_corrected"])):
            if r["heldout_corrected"] > best:
                best = r["heldout_corrected"]
                env.append((r["heldout_unmatched"], best))
        ax1.step([x for x, _ in env], [y for _, y in env], where="post", color=c, lw=1.8, label=LABEL[det])
        ax2.scatter([r["tuning_corrected"] for r in rs], [max(r["quiet_per_day"], FLOOR) for r in rs],
                    s=7, color=c, alpha=0.35, linewidths=0)
        xs = sorted({r["tuning_corrected"] for r in rs})
        low = [(x, max(min(r["quiet_per_day"] for r in rs if r["tuning_corrected"] >= x), FLOOR)) for x in xs]
        ax2.step([x for x, _ in low], [y for _, y in low], where="pre", color=c, lw=1.8, label=LABEL[det])
        for r in rs:
            if r["selected"]:
                kw = dict(marker=MARK[r["selected"]], ms=6.5, color=c, mec="white", mew=0.9, zorder=5,
                          clip_on=False)
                if r["selected"] == "rule at 26 unmatched":
                    kw.update(mfc="white", mec=c, mew=1.3)
                ax1.plot(r["heldout_unmatched"], r["heldout_corrected"], **kw)
                ax2.plot(r["tuning_corrected"], max(r["quiet_per_day"], FLOOR), **kw)
    ax1.set_xscale("symlog", linthresh=10)
    ax1.set_xlim(0, 100)
    ax1.set_xticks([0, 2, 5, 10, 20, 50, 100])
    ax1.set_xticklabels(["0", "2", "5", "10", "20", "50", "100"])
    ax1.set_ylim(0, 55)
    ax1.set_xlabel("unmatched alarms, held-out replay", color=MUTED)
    ax1.set_ylabel("chance-corrected detections of 66", color=MUTED)
    ax1.set_title("(a)  Held-out replay (post hoc)", loc="left", fontweight="bold", color=INK)
    ax2.set_yscale("log")
    ax2.set_ylim(FLOOR * 0.9, 300)
    ax2.set_yticks([FLOOR, 1, 3, 10, 30, 100])
    ax2.set_yticklabels(["0", "1", "3", "10", "30", "100"])
    ax2.set_xlim(0, 60)
    ax2.set_xlabel("chance-corrected detections of 66, tuning replays", color=MUTED)
    ax2.set_ylabel("alarms per day, 102.8 h of quiet data", color=MUTED)
    ax2.set_title("(b)  Quiet windows", loc="left", fontweight="bold", color=INK)
    for ax in (ax1, ax2):
        ax.grid(color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("right", "top"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(MUTED)
            ax.spines[side].set_linewidth(0.6)
        ax.tick_params(colors=MUTED, length=3, width=0.6)
    h = [plt.Line2D([], [], color=COLOR["learned"], lw=1.8), plt.Line2D([], [], color=COLOR["stalta"], lw=1.8)]
    h += [plt.Line2D([], [], ls="", marker=m, color=MUTED, mec="white" if k != "rule at 26 unmatched" else MUTED,
                     mfc=MUTED if k != "rule at 26 unmatched" else "white", ms=6) for k, m in MARK.items()]
    fig.legend(h, ["learned detector", "STA/LTA", "adopted rule", "STA/LTA, equal unmatched",
                   "STA/LTA, equal detections", "STA/LTA, the rule at 26 unmatched"],
               loc="upper center", ncols=3, frameon=False, fontsize=7, handlelength=1.8,
               columnspacing=1.6, bbox_to_anchor=(0.5, 1.0))
    out = Path(a.out)
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"-> {out}.pdf, {out}.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("data")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    p = sub.add_parser("plot")
    p.add_argument("--csv", required=True)
    p.add_argument("--out", required=True)
    a = ap.parse_args()
    {"data": data, "plot": plot}[a.cmd](a)


if __name__ == "__main__":
    main()
