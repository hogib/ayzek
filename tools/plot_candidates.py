#!/usr/bin/env python3
"""Publication figures for the quiet-window alarms that pass the S-P test.

    uv run python tools/plot_candidates.py \
        --csv data/runs/quiet100/misfires.csv --data data/quiet100 \
        --records data/runs/quiet100 --out data/runs/quiet100/figures

One figure per candidate (two stations side by side, three components each) and
one composite over all of them, written at 300 dpi PNG and as vector PDF.

Design choices, so the figure survives a journal:

- Components are **stacked and labelled**, so identity is positional, not
  colour-coded: the traces stay near-black and the figure reads in greyscale.
- The P and S picks are the only coloured marks, and they also differ in line
  style (P solid, S dashed) so colour is never the sole encoding. The two hues
  are the validated categorical slots 1 and 2 (blue #2a78d6, orange #eb6834),
  which keep dE >= 24 under protanopia and deuteranopia.
- Each panel is annotated with what makes the candidate an earthquake rather
  than a coincidence: the S-P interval, the distance it implies, and the pick
  confidences. The caption carries the cross-station test.
- Axes are recessive, there is no grid behind the traces, and the per-panel
  peak amplitude is stated in counts instead of a tick-laden amplitude axis.
  All three components share one scale, so their relative sizes are real.
"""

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
# TrueType, not Type 3: many journals reject Type-3 fonts outright, and this
# keeps the text selectable and editable in the vector file. Nothing in the
# figure is rasterised, so the PDF is vector throughout at any zoom.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["pdf.compression"] = 6
matplotlib.rcParams["savefig.transparent"] = False
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from obspy import read, UTCDateTime

P_COLOR, S_COLOR = "#2a78d6", "#eb6834"
INK, MUTED = "#0b0b0b", "#6b6b6b"
VP, VS = 6.0, 3.5
COMPS = ("Z", "N", "E")


def hav(a, b, c, d):
    r = 6371.0
    p = math.radians
    x = (math.sin(p(c - a) / 2) ** 2
         + math.cos(p(a)) * math.cos(p(c)) * math.sin(p(d - b) / 2) ** 2)
    return 2 * r * math.asin(math.sqrt(x))


def picks_for(rec_path, codes, t, span=90.0):
    """Best (highest P+S confidence) pick per station within `span` of t."""
    best = {}
    with open(rec_path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] != "P" or len(p) < 7:
                continue
            code, pt, st_, pp, sp = p[1], float(p[3]), float(p[4]), float(p[5]), float(p[6])
            if code not in codes or abs(pt - t) > span or st_ <= pt:
                continue
            if code not in best or pp + sp > best[code][2] + best[code][3]:
                best[code] = (pt, st_, pp, sp)
    return best


def panel(ax, path, station, pick, freqmin, freqmax, pre, post):
    """Three stacked components for one station, with its picks marked."""
    pt, st_, pp, sp = pick
    t0 = UTCDateTime(pt - pre)
    t1 = UTCDateTime(pt + post)
    stream = read(str(path), starttime=t0 - 5, endtime=t1 + 5)
    stream.merge(method=1, fill_value=None)
    stream = stream.split()
    stream.detrend("demean")
    stream.taper(0.05)
    stream.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
    stream.trim(t0, t1)

    # One shared scale over the three components, so relative amplitude is real.
    peak = max((float(np.abs(tr.data).max()) for tr in stream if len(tr.data)), default=1.0)
    peak = peak or 1.0

    drawn = 0
    for i, comp in enumerate(COMPS):
        sel = [tr for tr in stream if tr.stats.channel.endswith(comp)]
        offset = -i * 2.4
        if not sel:
            ax.text(pre * -1 + 0.3, offset, f"HH{comp}: no data", color=MUTED,
                    fontsize=7, va="center")
            continue
        tr = max(sel, key=lambda t: len(t.data))
        x = np.arange(len(tr.data)) * tr.stats.delta - pre
        ax.plot(x, tr.data / peak + offset, color=INK, lw=0.45,
                solid_joinstyle="bevel", rasterized=False, antialiased=True)
        ax.text(-pre + 0.2, offset + 0.95, f"HH{comp}", color=MUTED, fontsize=7,
                va="top", ha="left")
        drawn += 1

    for when, color, style, label in ((0.0, P_COLOR, "-", "P"),
                                      (st_ - pt, S_COLOR, "--", "S")):
        ax.axvline(when, color=color, ls=style, lw=1.1, zorder=3)
        ax.text(when, 1.55, label, color=color, fontsize=8, fontweight="bold",
                ha="center", va="bottom")

    sp_t = st_ - pt
    dist = sp_t / (1 / VS - 1 / VP)
    ax.set_title(f"{station}   S−P = {sp_t:.1f} s  →  {dist:.0f} km"
                 f"        P {pp:.2f} / S {sp:.2f}        peak {peak:,.0f} counts",
                 fontsize=8.5, color=INK, pad=6, loc="left")
    ax.set_xlim(-pre, post)
    ax.set_ylim(-2.4 * 2 - 1.4, 2.0)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(axis="x", colors=MUTED, labelsize=7.5, length=3, width=0.6)
    return drawn


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--data", required=True, help="the quiet dataset root")
    ap.add_argument("--records", required=True, help="directory of .rec files")
    ap.add_argument("--stations-csv", default="models/stations.csv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--freqmin", type=float, default=1.0)
    ap.add_argument("--freqmax", type=float, default=20.0)
    ap.add_argument("--pre", type=float, default=5.0)
    ap.add_argument("--dpi", type=int, default=600,
                   help="raster dpi for the PNG; the PDF is vector regardless")
    ap.add_argument("--verdict", default="earthquake",
                   help="which rows to plot: earthquake (default), misfire, unclassified")
    ap.add_argument("--limit", type=int, default=0, help="plot at most N of them")
    ap.add_argument("--label", default="Candidate", help="word used in titles")
    a = ap.parse_args()

    st = {}
    with open(a.stations_csv, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            st[r.get("code") or r.get("station")] = (
                float(r.get("lat") or r["latitude"]), float(r.get("lon") or r["longitude"]))

    cands = [r for r in csv.DictReader(open(a.csv)) if r["verdict"] == a.verdict]
    if not cands:
        raise SystemExit(f"no rows with verdict {a.verdict!r} in {a.csv}")
    if a.limit:
        # Spread the sample over the run rather than taking the first N.
        step = max(1, len(cands) // a.limit)
        cands = cands[::step][:a.limit]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    print(f"{len(cands)} candidate(s)")

    made = []
    for n, r in enumerate(cands, 1):
        rec = Path(a.records) / r["record"]
        window = r["record"].replace(".rec", "")
        t = float(r["time_utc"])
        A, B = r["station_a"], r["station_b"]
        pk = picks_for(rec, {A, B}, t)
        if len(pk) < 2:
            print(f"  {n}: picks missing, skipped")
            continue
        sp_max = max(pk[A][1] - pk[A][0], pk[B][1] - pk[B][0])
        post = max(20.0, sp_max + 12.0)

        fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9), dpi=a.dpi,
                                 constrained_layout=True)
        for ax, code in zip(axes, (A, B)):
            panel(ax, Path(a.data) / window / f"{code}.mseed", code, pk[code],
                  a.freqmin, a.freqmax, a.pre, post)
            ax.set_xlabel("seconds from P", fontsize=8, color=MUTED)

        da = (pk[A][1] - pk[A][0]) / (1 / VS - 1 / VP)
        db = (pk[B][1] - pk[B][0]) / (1 / VS - 1 / VP)
        sep = hav(*st[A], *st[B])
        pred = (da - db) / VP
        obs = float(r["dt_s"])
        err = min(abs(pred - obs), abs(-pred - obs))
        when = UTCDateTime(pk[A][0]).strftime("%Y-%m-%d %H:%M:%S")
        fig.suptitle(
            f"{a.label} {n}  ·  {when} UTC  ·  {A}–{B}, {sep:.0f} km apart  ·  "
            f"band-pass {a.freqmin:g}–{a.freqmax:g} Hz\n"
            f"distances from S−P: {da:.0f} and {db:.0f} km  →  predicted Δt "
            f"{abs(pred):.1f} s, observed {abs(obs):.1f} s, error {err:.1f} s",
            fontsize=9.5, color=INK, y=1.10)
        fig.legend(handles=[Line2D([], [], color=P_COLOR, ls="-", lw=1.1, label="P pick"),
                            Line2D([], [], color=S_COLOR, ls="--", lw=1.1, label="S pick")],
                   loc="lower right", frameon=False, fontsize=8, ncols=2,
                   bbox_to_anchor=(1.0, -0.06))
        stem = out / f"{a.label.lower().replace(' ', '_')}_{n}_{window[6:]}_{A}_{B}"
        for ext in ("png", "pdf"):
            fig.savefig(f"{stem}.{ext}", dpi=a.dpi, bbox_inches="tight",
                        facecolor="white")
        plt.close(fig)
        made.append((n, r, pk, da, db, sep, pred, obs, err, when, window, A, B, post))
        print(f"  {n}: {stem.name}.png/.pdf")

    # Composite: one row per candidate, for a single paper figure.
    if made:
        fig, axes = plt.subplots(len(made), 2, dpi=a.dpi,
                                 figsize=(10.5, 3.5 * len(made)),
                                 constrained_layout=True, squeeze=False)
        for row, (n, r, pk, da, db, sep, pred, obs, err, when, window, A, B, post) in enumerate(made):
            for ax, code in zip(axes[row], (A, B)):
                panel(ax, Path(a.data) / window / f"{code}.mseed", code, pk[code],
                      a.freqmin, a.freqmax, a.pre, post)
                if row == len(made) - 1:
                    ax.set_xlabel("seconds from P", fontsize=8, color=MUTED)
            axes[row][0].text(-0.055, 0.5, f"({chr(96 + n)})", transform=axes[row][0].transAxes,
                              fontsize=11, fontweight="bold", va="center", ha="right", color=INK)
            axes[row][1].text(1.0, 1.28, f"{when} UTC · Δt predicted {abs(pred):.1f} s, "
                              f"observed {abs(obs):.1f} s (error {err:.1f} s)",
                              transform=axes[row][1].transAxes, fontsize=8,
                              color=MUTED, ha="right", va="bottom")
        # Below the axes, not over panel (a)'s title.
        fig.legend(handles=[Line2D([], [], color=P_COLOR, ls="-", lw=1.1, label="P pick"),
                            Line2D([], [], color=S_COLOR, ls="--", lw=1.1, label="S pick")],
                   loc="lower center", frameon=False, fontsize=8, ncols=2,
                   bbox_to_anchor=(0.5, -0.012))
        for ext in ("png", "pdf"):
            fig.savefig(out / f"{a.label.lower().replace(' ', '_')}s_all.{ext}", dpi=a.dpi,
                        bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"  composite: {a.label.lower()}s_all.png/.pdf")
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
