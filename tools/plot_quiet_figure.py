#!/usr/bin/env python3
"""The paper's figure of the quiet-window alarms: an earthquake, a misfire, and why.

    uv run --with matplotlib --with obspy --with numpy python tools/plot_quiet_figure.py \\
        --out ../ayzek_paper/figures/fig6_quiet_alarms

(a) an uncatalogued earthquake behind an alarm of the learned detector, (b) a
KIRK-pulse alarm: KIRK's instrumental pulse coinciding with background at KAND, (c) where
KIRK's detections fall in a 10-minute cycle, (d) how long the signal behind each
alarm lasts, for both detectors (tools/transient_duration.py).

Design: traces near-black, components stacked and labelled, one amplitude scale
per station; picks dark grey, P solid and S dashed, so that colour is left to
the two detectors, which keep the hues of the paper's other figures (learned
#9a3412, STA/LTA #1d4ed8; validated pair, CVD dE 30). Vector PDF and 300 dpi PNG.
"""
import argparse
import csv
import glob
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime, read

sys.path.insert(0, str(Path(__file__).parent))
from misfire_review import cluster, read_record  # noqa: E402
from plot_candidates import picks_for  # noqa: E402

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e5e7eb"
PICK = "#3a3a3a"
LEARNED, STALTA = "#9a3412", "#1d4ed8"
VP, VS = 6.0, 3.5
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8.5, "font.family": "DejaVu Sans",
                     "pdf.fonttype": 42, "svg.fonttype": "none"})


def traces(ax, path, t_ref, pre, post, marks, title):
    """Three stacked components band-passed 1-20 Hz, one scale; marks = [(t, style, label)]."""
    t0, t1 = UTCDateTime(t_ref - pre), UTCDateTime(t_ref + post)
    st = read(str(path), starttime=t0 - 5, endtime=t1 + 5)
    st.merge(method=1, fill_value=None)
    st = st.split()
    st.detrend("demean")
    st.taper(0.05)
    st.filter("bandpass", freqmin=1, freqmax=20, corners=4, zerophase=True)
    st.trim(t0, t1)
    peak = max(float(np.abs(tr.data).max()) for tr in st) or 1.0
    for i, comp in enumerate("ZNE"):
        sel = [tr for tr in st if tr.stats.channel.endswith(comp)]
        if not sel:
            continue
        tr = max(sel, key=lambda t: len(t.data))
        x = tr.times() + (tr.stats.starttime.timestamp - t_ref)
        ax.plot(x, tr.data / peak - 2.3 * i, color=INK, lw=0.45)
        ax.text(-pre + 0.3, -2.3 * i + 0.75, f"HH{comp}", color=MUTED, fontsize=7, va="bottom")
    for t, style, label in marks:
        ax.axvline(t, color=PICK, ls=style, lw=1.0, zorder=3)
        ax.text(t + 0.25, 1.45, label, color=PICK, fontsize=7.5, ha="left", va="center", fontweight="bold")
    ax.set_xlim(-pre, post)
    ax.set_ylim(-2.3 * 2 - 1.3, 1.8)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(axis="x", colors=MUTED, labelsize=7, length=3, width=0.6)
    ax.set_title(f"{title}   peak {peak:,.0f} counts", loc="left", color=INK, pad=3, fontsize=8)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", default="data/runs")
    ap.add_argument("--data", default="data/quiet100")
    ap.add_argument("--out", required=True, help="path without extension")
    a = ap.parse_args()
    R, D = Path(a.runs), Path(a.data)

    fig = plt.figure(figsize=(7.0, 8.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.95], hspace=0.75, wspace=0.22,
                          left=0.06, right=0.98, top=0.925, bottom=0.07)

    def heading(ax, text):
        """A row heading just above the row's axes titles."""
        y = ax.get_position().y1 + 0.032
        fig.text(0.06, y, text, fontsize=9, fontweight="bold", color=INK, va="bottom")

    # (a) The earthquake: BAND-CMH, 2024-05-23 04:43 UTC, EMSC M2.4.
    rec = R / "quiet100/quiet_20240523_0352.rec"
    t = 1716439378.66
    pk = picks_for(rec, {"BAND", "CMH"}, t)
    for j, s in enumerate(("CMH", "BAND")):
        pt, st_, _, _ = pk[s]
        d = (st_ - pt) / (1 / VS - 1 / VP)
        ax = fig.add_subplot(gs[0, j])
        traces(ax, D / "quiet_20240523_0352" / f"{s}.mseed", pt, 4, 26,
               [(0.0, "-", "P"), (st_ - pt, "--", "S")], f"{s}: S$-$P {st_ - pt:.1f} s, {d:.0f} km")
        ax.set_xlabel("seconds from the station's P pick", color=MUTED, fontsize=7.5)
        if j == 0:
            heading(ax, "(a)  Alarm 2: uncataloged earthquake, EMSC M2.4, origin 2024-05-23 04:42:52 UTC")

    # (b) The misfire: KAND-KIRK, 2024-06-09 13:50 UTC, KIRK's pulse.
    rec = R / "quiet100/quiet_20240609_1210.rec"
    _, det, _, _ = read_record(str(rec))
    t = 1717941053.50
    # One time base for both stations, so the coincidence is visible.
    ref = min(det["KIRK"], key=lambda x: abs(x - t))
    for j, s in enumerate(("KIRK", "KAND")):
        ax = fig.add_subplot(gs[1, j])
        traces(ax, D / "quiet_20240609_1210" / f"{s}.mseed", ref, 4, 26, [], s)
        ax.set_xlabel("seconds from KIRK's reported P", color=MUTED, fontsize=7.5)
        if j == 0:
            heading(ax, "(b)  Alarm 12, KIRK pulse, 2024-06-09 13:51 UTC: pulse at KIRK, noise at KAND")

    # (c) KIRK's detections in the 10-minute cycle.
    phases = []
    for path in sorted(glob.glob(str(R / "quiet100/*.rec"))):
        _, det, _, _ = read_record(path)
        phases += [x % 600 for x in cluster(det.get("KIRK", []), 40.0)]
    ax = fig.add_subplot(gs[2, 0])
    counts, edges = np.histogram(phases, bins=np.arange(0, 610, 10))
    ax.bar(edges[:-1], counts, width=8, align="edge", color=LEARNED, linewidth=0)
    top = int(np.argmax(counts))
    ax.annotate(f"{counts[top]} of {len(phases)}", (edges[top] + 4, counts[top]),
                xytext=(edges[top] + 90, counts[top] * 0.85), fontsize=7.5, color=INK,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    ax.set_xlim(0, 600)
    ax.set_xticks([0, 120, 240, 360, 480, 600])
    ax.set_xlabel("seconds after each 10-minute mark (UTC)", color=MUTED, fontsize=7.5)
    ax.set_ylabel("KIRK detections", color=MUTED, fontsize=7.5)
    ax.set_title("(c)  KIRK's detections in the 10-minute cycle", loc="left", fontweight="bold", pad=14)

    # (d) Duration of the signal behind each alarm, both detectors.
    ax = fig.add_subplot(gs[2, 1])
    n = {}
    for run, color, label in (("quiet100", LEARNED, "learned"), ("quiet100_stalta", STALTA, "STA/LTA")):
        v = np.sort([max(float(r["seconds_above"]), 0.02)
                     for r in csv.DictReader(open(R / run / "durations.csv"))])
        n[run] = len(v)
        y = np.arange(1, len(v) + 1) / len(v)
        ax.step(v, y, where="post", color=color, lw=2)
    for x0 in (0.5, 3.0):
        ax.axvline(x0, color=MUTED, lw=0.6, ls=":")
    for xm, name in ((0.11, "spike"), (1.2, "transient"), (8, "extended")):
        ax.text(xm, 1.035, name, color=MUTED, fontsize=7, ha="center", va="bottom")
    ax.set_xscale("log")
    ax.set_xlim(0.02, 25)
    ax.set_xticks([0.03, 0.1, 0.3, 1, 3, 10])
    ax.set_xticklabels(["0.03", "0.1", "0.3", "1", "3", "10"])
    ax.minorticks_off()
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("seconds the signal stays above 4× noise", color=MUTED, fontsize=7.5)
    ax.set_ylabel("fraction of station detections", color=MUTED, fontsize=7.5)
    ax.set_title("(d)  Signal behind each alarm", loc="left", fontweight="bold", pad=14)
    ax.legend([plt.Line2D([], [], color=LEARNED, lw=2), plt.Line2D([], [], color=STALTA, lw=2)],
              [f"learned detector ({n['quiet100']})", f"STA/LTA ({n['quiet100_stalta']})"],
              frameon=False, fontsize=7,
              loc="upper left")
    for axx in (fig.axes[-2], fig.axes[-1]):
        axx.grid(axis="y", color=GRID, lw=0.6)
        axx.set_axisbelow(True)
        for side in ("right", "top"):
            axx.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            axx.spines[side].set_color(MUTED)
            axx.spines[side].set_linewidth(0.6)
        axx.tick_params(colors=MUTED, labelsize=7, length=3, width=0.6)

    out = Path(a.out)
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"-> {out}.pdf, {out}.png")


if __name__ == "__main__":
    main()
