#!/usr/bin/env python3
"""Excess coincidence over the time-slide null against station separation.

    uv run python tools/plot_coincidence.py \
        --arm model:data/runs/quiet100 --arm STA/LTA:data/runs/quiet100_stalta \
        --out data/runs/quiet100/figures/coincidence_excess

This is the form of the paper's Figure 7, built from the 102.8 h quiet set
rather than from the 2,977 station-day network run, which is not on this
machine. Two differences matter and are printed on the figure:

- **One operating point per detector, not three alarm budgets.** The quiet
  replays were run at the shipped trigger settings, so there is no budget axis.
  The two detectors take its place, which is arguably more informative here:
  the model and STA/LTA sit at opposite ends of the same axis.
- **Far fewer station-days**, so per-pair counts are small and most pairs are
  unresolved. They are drawn hollow, as in the paper.

For each station pair the observed coincidences are counted within
`d/vp + slack`, and the null is the mean over `--slides` circular shifts of
each station's own detection series inside its window — preserving each
station's rate and clustering while destroying real coincidences. Excess is
observed ÷ null; a pair at 1.0 is indistinguishable from independence.
"""

import argparse
import itertools
import math
import statistics
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from misfire_review import read_record, cluster, hav  # noqa: E402

INK, MUTED, GRID = "#0b0b0b", "#6b6b6b", "#d8d8d6"
# The paper's own series colours, sampled from the v3 figures: the learned
# detector is rust and STA/LTA is blue throughout, so this figure must not
# invert them. Marker shape carries the same distinction, so identity survives
# greyscale printing and does not rest on colour alone.
COLORS = ["#993412", "#1d4ed7", "#1baf7a", "#eda100"]
MARKERS = ["o", "^", "s", "D"]


def pair_counts(paths, vp, slack, coda, slide, slides):
    """-> {(a, b): (observed, null_mean, null_sd, separation_km)}"""
    obs, nulls, seps = {}, {}, {}
    for path in paths:
        stations, det, span, _ = read_record(path)
        if not span:
            continue
        det = {s: cluster(ts, coda) for s, ts in det.items()}
        codes = sorted(det)
        for a, b in itertools.combinations(sorted(stations), 2):
            seps.setdefault((a, b), hav(*stations[a], *stations[b]))

        def count(shift=None):
            got = {}
            series = {}
            for s, ts in det.items():
                off = (shift or {}).get(s, 0.0)
                if off:
                    lo, hi = span
                    n = hi - lo
                    ts = sorted(lo + ((t - lo + off) % n) for t in ts)
                series[s] = ts
            for a, b in itertools.combinations(codes, 2):
                if a not in stations or b not in stations:
                    continue
                w = hav(*stations[a], *stations[b]) / vp + slack
                ta, tb, j, n = series[a], series[b], 0, 0
                for t in ta:
                    while j < len(tb) and tb[j] < t - w:
                        j += 1
                    k = j
                    while k < len(tb) and tb[k] <= t + w:
                        n += 1
                        k += 1
                got[(a, b)] = got.get((a, b), 0) + n
            return got

        for k, v in count().items():
            obs[k] = obs.get(k, 0) + v
        if len(codes) > 1:
            for i in range(1, slides + 1):
                sh = {c: (j + 1) * i * slide for j, c in enumerate(codes)}
                sh[codes[0]] = 0.0
                for k, v in count(sh).items():
                    nulls.setdefault(k, []).append(v)
    out = {}
    for k in obs:
        vals = nulls.get(k, [])
        # Sum across windows happens above for obs; the null list holds one
        # entry per (window, slide), so the per-slide network total is the mean
        # times the number of windows that contributed.
        mu = statistics.fmean(vals) * (len(vals) / max(1, slides)) if vals else 0.0
        sd = (statistics.pstdev(vals) * math.sqrt(len(vals) / max(1, slides))
              if len(vals) > 1 else 0.0)
        out[k] = (obs[k], mu, sd, seps.get(k, float("nan")))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", action="append", required=True,
                    help="LABEL:directory of .rec files; repeatable")
    ap.add_argument("--out", required=True, help="output path without extension")
    ap.add_argument("--vp", type=float, default=6.0)
    ap.add_argument("--slack", type=float, default=3.0)
    ap.add_argument("--coda", type=float, default=40.0)
    ap.add_argument("--slide", type=float, default=600.0)
    ap.add_argument("--slides", type=int, default=50)
    ap.add_argument("--dpi", type=int, default=600)
    a = ap.parse_args()

    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=a.dpi, constrained_layout=True)
    summary = []
    for i, spec in enumerate(a.arm):
        label, d = spec.split(":", 1)
        paths = sorted(Path(d).glob("*.rec"))
        if not paths:
            sys.exit(f"no .rec files in {d}")
        counts = pair_counts(paths, a.vp, a.slack, a.coda, a.slide, a.slides)
        color = COLORS[i % len(COLORS)]
        marker = MARKERS[i % len(MARKERS)]
        xs_r, ys_r, xs_u, ys_u = [], [], [], []
        for (sa, sb), (o, mu, sd, sep) in sorted(counts.items(), key=lambda kv: kv[1][3]):
            if o == 0 and mu == 0:
                continue
            exc = o / mu if mu > 0 else float("nan")
            if not np.isfinite(exc) or exc <= 0:
                continue
            resolved = sd > 0 and abs(o - mu) > 3 * sd
            (xs_r if resolved else xs_u).append(sep)
            (ys_r if resolved else ys_u).append(exc)
        ax.scatter(xs_u, ys_u, s=38, marker=marker, facecolors="none",
                   edgecolors=color, linewidths=1.0, zorder=3)
        ax.scatter(xs_r, ys_r, s=38, marker=marker, color=color, zorder=4)
        tot_o = sum(v[0] for v in counts.values())
        tot_n = sum(v[1] for v in counts.values())
        summary.append(f"{label}: {tot_o} observed, {tot_n:.1f} null "
                       f"({tot_o / tot_n:.2f}x)" if tot_n else f"{label}: {tot_o} observed")
        print(f"{label}: {len(counts)} pairs, {tot_o} observed, null {tot_n:.1f}, "
              f"{len(xs_r)} resolved")

    ax.axhline(1.0, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.text(ax.get_xlim()[1], 1.0, " independence", color=MUTED, fontsize=8,
            va="center", ha="left")
    ax.set_yscale("log")
    ax.set_xlabel("station separation (km)", fontsize=9.5, color=INK)
    ax.set_ylabel("coincidences ÷ time-slide null", fontsize=9.5, color=INK)
    ax.tick_params(colors=MUTED, labelsize=8.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
        ax.spines[s].set_linewidth(0.6)
    ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)

    handles = []
    for i, spec in enumerate(a.arm):
        handles.append(Line2D([], [], marker=MARKERS[i % len(MARKERS)], ls="none",
                              color=COLORS[i % len(COLORS)], markersize=6,
                              label=spec.split(":", 1)[0]))
    handles.append(Line2D([], [], marker="o", ls="none", markerfacecolor="none",
                          markeredgecolor=MUTED, markersize=6,
                          label="unresolved (within 3σ of the null)"))
    ax.legend(handles=handles, frameon=False, fontsize=8.5, loc="upper right")
    # No title inside the figure -- it belongs in the caption. The numbers it
    # used to carry are printed to stdout so they can be pasted there.
    print("caption numbers: " + "; ".join(summary))
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=a.dpi, bbox_inches="tight", facecolor="white")
    print(f"-> {a.out}.pdf/.png")


if __name__ == "__main__":
    main()
