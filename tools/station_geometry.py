"""Scores station subsets by network geometry, without waveforms.

For a grid of possible epicentres over the region around the stations, weighted
by catalogue seismicity (event counts, or the cells with any event), and for
every subset of the given stations, it computes:

  coverage      fraction of the grid weight with at least two subset
                stations within `reach` km (the minimum for a declaration)
  alert_s       P travel time to the second-nearest station within reach, plus
                `latency` (time from P arrival to declaration)
  blind_km      radius the S wave has travelled at the alert time
  location_km   95% horizontal error of a P+S location, from the linearised
                travel-time equations with pick errors sigma_p and sigma_s
  gap_deg       azimuthal gap of the stations within reach
  worst_*       coverage and alert time with the most damaging single station
                removed

Metrics other than coverage are averaged over covered grid nodes (location
error: weighted median). Uniform velocities, fixed depth.

With --at LAT,LON[,DEPTH] --compare SUBSETS.csv, it also predicts the alert
time at one epicentre for every subset in a network_subsets output and
compares it with the observed alert time of the main event.

    uv run --with numpy --with pandas --with matplotlib python tools/station_geometry.py \\
        --stations ARNA,BAND,CATL,DEMI,ELBA,KIRK,MANT,SEMS --catalog CATALOG.csv --max-k 4 --plot best.png
"""
import argparse
import itertools
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

VP, VS = 6.0, 3.5
CHI2_95_2D = 5.991


def load_stations(path, codes):
    table = pd.read_csv(path).drop_duplicates("code").set_index("code")
    missing = [c for c in codes if c not in table.index]
    if missing:
        sys.exit(f"no coordinates for {missing}")
    return table.loc[codes, ["lat", "lon"]].to_numpy(float)


def haversine(lat1, lon1, lat2, lon2):
    r = np.pi / 180
    a = np.sin((lat2 - lat1) * r / 2) ** 2 + np.cos(lat1 * r) * np.cos(lat2 * r) * np.sin((lon2 - lon1) * r / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def grid_and_weights(stations, catalog, margin, step, since, min_mag, weighting):
    lat0, lat1 = stations[:, 0].min() - margin, stations[:, 0].max() + margin
    lon0, lon1 = stations[:, 1].min() - margin, stations[:, 1].max() + margin
    lats = np.arange(lat0, lat1 + 1e-9, step)
    lons = np.arange(lon0, lon1 + 1e-9, step)
    glat, glon = [a.ravel() for a in np.meshgrid(lats, lons, indexing="ij")]
    if weighting == "uniform":
        return glat, glon, np.ones_like(glat)
    cat = pd.read_csv(catalog, encoding="utf-8-sig")
    year = pd.to_numeric(cat.Date.str.slice(6, 10), errors="coerce")
    cat = cat[(year >= since) & (cat.Magnitude >= min_mag)
              & cat.Latitude.between(lat0, lat1) & cat.Longitude.between(lon0, lon1)]
    i = np.clip(np.round((cat.Latitude.to_numpy() - lat0) / step).astype(int), 0, len(lats) - 1)
    j = np.clip(np.round((cat.Longitude.to_numpy() - lon0) / step).astype(int), 0, len(lons) - 1)
    counts = np.zeros((len(lats), len(lons)))
    np.add.at(counts, (i, j), 1.0)
    if weighting == "footprint":
        counts = (counts > 0).astype(float)
    return glat, glon, counts.ravel()


class Geometry:
    """Per-node distances, azimuths and travel-time derivatives for all stations."""

    def __init__(self, glat, glon, stations, depth):
        self.depth = depth
        self.epi = haversine(glat[:, None], glon[:, None], stations[None, :, 0], stations[None, :, 1])  # (nodes, n)
        self.hyp = np.hypot(self.epi, depth)
        dx = (stations[None, :, 1] - glon[:, None]) * 111.32 * np.cos(np.radians(glat[:, None]))
        dy = (stations[None, :, 0] - glat[:, None]) * 110.57
        self.azimuth = np.degrees(np.arctan2(dx, dy)) % 360
        # d(travel time)/d(epicentre x, y) = -(station - epicentre) / (v * hypocentral distance)
        self.gx, self.gy = -dx / self.hyp, -dy / self.hyp


def score(geo, w, idx, reach, latency, sigma_p, sigma_s):
    """Metrics for the station subset `idx` (list of column indices)."""
    epi = geo.epi[:, idx]
    in_reach = epi <= reach
    n_in = in_reach.sum(axis=1)
    covered = n_in >= 2
    wsum = w.sum()
    coverage = w[covered].sum() / wsum if wsum > 0 else 0.0
    out = {"coverage": coverage}
    if not covered.any():
        return out | {"alert_s": np.nan, "blind_km": np.nan, "location_km": np.nan, "gap_deg": np.nan}

    tp = np.where(in_reach, geo.hyp[:, idx] / VP, np.inf)
    alert = np.sort(tp, axis=1)[:, 1] + latency
    blind = np.sqrt(np.maximum((VS * alert) ** 2 - geo.depth ** 2, 0))
    wc = w[covered]
    out["alert_s"] = np.average(alert[covered], weights=wc) if wc.sum() > 0 else np.nan
    out["blind_km"] = np.average(blind[covered], weights=wc) if wc.sum() > 0 else np.nan

    # Location covariance from P and S rows: unknowns (x, y, origin time).
    mask = in_reach[covered].astype(float)
    rows = []
    for v, sigma in ((VP, sigma_p), (VS, sigma_s)):
        rows.append(np.stack([geo.gx[covered][:, idx] / v, geo.gy[covered][:, idx] / v, np.ones_like(mask)], -1)
                    * (mask / sigma)[..., None])
    J = np.concatenate(rows, axis=1)                                   # (nodes, 2k, 3)
    N = np.einsum("nij,nik->njk", J, J)
    cov = np.linalg.inv(N + np.eye(3) * 1e-12)
    horizontal = cov[:, :2, :2]
    lam = np.linalg.eigvalsh(horizontal)[:, -1]
    err = np.sqrt(CHI2_95_2D * np.maximum(lam, 0))
    order = np.argsort(err)
    cw = np.cumsum(wc[order])
    out["location_km"] = err[order][np.searchsorted(cw, cw[-1] / 2)] if cw[-1] > 0 else np.nan

    az = np.where(in_reach[covered], geo.azimuth[covered][:, idx], np.nan)
    az = np.sort(az, axis=1)
    gaps = np.diff(np.concatenate([az, az[:, :1] + 360], axis=1), axis=1)
    out["gap_deg"] = np.average(np.nanmax(np.where(np.isnan(gaps), -1, gaps), axis=1), weights=wc) if wc.sum() > 0 else np.nan
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stations", required=True, help="comma-separated station codes")
    ap.add_argument("--station-table", default="models/stations.csv")
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--since", type=int, default=2010)
    ap.add_argument("--min-mag", type=float, default=2.5, help="catalogue events used for weights")
    ap.add_argument("--weighting", choices=["seismicity", "footprint", "uniform"], default="seismicity",
                    help="seismicity: catalogue event count per cell; footprint: 1 for every cell with "
                         "an event, 0 otherwise (removes the dominance of aftershock sequences); uniform")
    ap.add_argument("--margin", type=float, default=1.0, help="degrees around the stations")
    ap.add_argument("--step", type=float, default=0.1, help="grid step in degrees")
    ap.add_argument("--depth", type=float, default=10.0)
    ap.add_argument("--reach", type=float, default=150.0, help="km within which a station contributes")
    ap.add_argument("--latency", type=float, default=2.0, help="s from P arrival to declaration")
    ap.add_argument("--sigma-p", type=float, default=0.3)
    ap.add_argument("--sigma-s", type=float, default=0.5)
    ap.add_argument("--min-k", type=int, default=2)
    ap.add_argument("--max-k", type=int, default=None)
    ap.add_argument("--out", default=None, help="CSV of all subsets")
    ap.add_argument("--plot", default=None, help="PNG map of the best subset for each k")
    ap.add_argument("--at", default=None, help="LAT,LON[,DEPTH] for --compare")
    ap.add_argument("--compare", default=None, help="network_subsets CSV")
    a = ap.parse_args()

    codes = a.stations.split(",")
    st = load_stations(a.station_table, codes)
    max_k = a.max_k or len(codes)

    if a.compare:
        lat, lon, *rest = [float(x) for x in a.at.split(",")]
        depth = rest[0] if rest else a.depth
        obs = pd.read_csv(a.compare)
        obs = obs[obs.main_declared == 1]
        hyp = np.hypot(haversine(lat, lon, st[:, 0], st[:, 1]), depth)
        pos = {c: i for i, c in enumerate(codes)}
        pred = [np.sort(hyp[[pos[s] for s in names.split("+")]])[1] / VP for names in obs.stations]
        obs = obs.assign(predicted_p2_s=pred, latency_s=obs.main_alert_s - np.array(pred))
        rho = obs.predicted_p2_s.rank().corr(obs.main_alert_s.rank())
        print(f"{len(obs)} subsets that declared the main event")
        print(f"observed alert = P travel time to the 2nd station + latency: "
              f"latency median {obs.latency_s.median():.2f} s, IQR {obs.latency_s.quantile(.25):.2f}..{obs.latency_s.quantile(.75):.2f} s")
        print(f"rank correlation of predicted P2 time and observed alert: {rho:.3f}")
        worst = obs.reindex(obs.latency_s.abs().sort_values(ascending=False).index).head(5)
        print("largest residuals:\n" + worst[["stations", "predicted_p2_s", "main_alert_s", "latency_s"]].to_string(index=False))
        return

    glat, glon, w = grid_and_weights(st, a.catalog, a.margin, a.step, a.since, a.min_mag, a.weighting)
    geo = Geometry(glat, glon, st, a.depth)
    print(f"{len(codes)} stations, {len(glat)} grid nodes, weighting {a.weighting}, total weight {w.sum():.0f} "
          f"(catalogue M>={a.min_mag} since {a.since})")

    rows = []
    for k in range(a.min_k, max_k + 1):
        for combo in itertools.combinations(range(len(codes)), k):
            idx = list(combo)
            m = score(geo, w, idx, a.reach, a.latency, a.sigma_p, a.sigma_s)
            worst_cov, worst_alert = m["coverage"], m["alert_s"]
            if k > 2:
                for drop in idx:
                    rest = [i for i in idx if i != drop]
                    r = score(geo, w, rest, a.reach, a.latency, a.sigma_p, a.sigma_s)
                    worst_cov = min(worst_cov, r["coverage"])
                    worst_alert = np.nanmax([worst_alert, r["alert_s"]])
            else:
                worst_cov, worst_alert = 0.0, np.nan
            rows.append({"k": k, "stations": "+".join(codes[i] for i in idx), **m,
                         "worst_coverage": worst_cov, "worst_alert_s": worst_alert})
    df = pd.DataFrame(rows)
    if a.out:
        df.to_csv(a.out, index=False)

    # Ranking: coverage first (rounded to 1%), then alert time, then location error.
    df["rank_key"] = list(zip(-df.coverage.round(2), df.alert_s.round(1), df.location_km))
    pd.set_option("display.width", 200)
    for k, g in df.groupby("k"):
        best = g.sort_values("rank_key").head(3)
        print(f"\nk = {k} ({len(g)} subsets), best 3:")
        print(best[["stations", "coverage", "alert_s", "blind_km", "location_km", "gap_deg", "worst_coverage",
                    "worst_alert_s"]].round(2).to_string(index=False))

    if a.plot:
        plot(a, codes, st, glat, glon, w, geo, df)


def plot(a, codes, st, glat, glon, w, geo, df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ks = sorted(df.k.unique())
    show = [k for k in (2, 3, 4) if k in ks] + ([ks[-1]] if ks[-1] not in (2, 3, 4) else [])
    fig, axes = plt.subplots(1, len(show), figsize=(5 * len(show), 5), squeeze=False)
    for ax, k in zip(axes[0], show):
        best = df[df.k == k].sort_values("rank_key").iloc[0]
        idx = [codes.index(c) for c in best.stations.split("+")]
        in_reach = geo.epi[:, idx] <= a.reach
        tp = np.where(in_reach, geo.hyp[:, idx] / VP, np.inf)
        alert = np.sort(tp, axis=1)[:, 1] + a.latency
        blind = np.where(in_reach.sum(1) >= 2, np.sqrt(np.maximum((VS * alert) ** 2 - a.depth ** 2, 0)), np.nan)
        sc = ax.scatter(glon, glat, c=blind, s=6, marker="s", cmap="viridis_r", vmin=0, vmax=120)
        quake = w > 0
        ax.scatter(glon[quake], glat[quake], s=np.clip(w[quake], 1, 30) * 0.4, c="k", alpha=0.25, lw=0)
        others = [i for i in range(len(codes)) if i not in idx]
        ax.scatter(st[others, 1], st[others, 0], marker="^", s=60, facecolors="none", edgecolors="red")
        ax.scatter(st[idx, 1], st[idx, 0], marker="^", s=80, c="red", edgecolors="k")
        for i, c in enumerate(codes):
            ax.annotate(c, (st[i, 1], st[i, 0]), xytext=(3, 3), textcoords="offset points", fontsize=7)
        ax.set_title(f"k={k}: {best.stations}\ncoverage {best.coverage:.0%}, mean blind zone {best.blind_km:.0f} km",
                     fontsize=9)
        ax.set_aspect(1 / math.cos(math.radians(float(np.mean(glat)))))
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")
    fig.colorbar(sc, ax=axes[0].tolist(), label="blind-zone radius at alert (km); blank: < 2 stations in reach")
    fig.savefig(a.plot, dpi=130, bbox_inches="tight")
    print(f"\nmap: {a.plot}")


if __name__ == "__main__":
    main()
