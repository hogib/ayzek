"""Magnitude from peak P displacement: fit and evaluation.

For each Pd window tau (pd.py) the relation

    log10 Pd_ij = alpha + beta M_j + gamma log10 R_ij + s_i + e_ij,   e ~ N(0, sigma^2)

is fitted on the train split, where i is the station, j the event,
R_ij = sqrt(D_ij^2 + DEPTH_KM^2) with D the epicentral distance, and s_i a
station term (sum zero over the stations that have one), which absorbs site
amplification and errors in the station's response metadata. The fixed
depth matches deployment, where D comes from the geometry head and the depth
is unknown.

**Censoring.** Displacement emphasises the long periods at which the
microseism is strongest, and the pre-P noise level `pd_noise` is comparable
to the Pd of a small regional event: at tau = 3 s only 2-5% of the M <= 3
values on fdsn_wide_x exceed MIN_SNR times it, against 75-90% above M 4. A
least-squares fit restricted to the values above that threshold therefore
sees only the brightest small events; it overestimates log10 Pd at small M,
which flattens the fitted slope, and the flattened slope, extrapolated,
overestimates the largest events (by about one magnitude unit above M 5 on
fdsn_wide_x). The relation is therefore fitted by censored maximum
likelihood (a Tobit model): a value below the threshold enters as the
upper bound log10(MIN_SNR pd_noise) on log10 Pd, with likelihood
Phi((c - mu) / sigma), and a value above it with the normal density.

**Event magnitude.** An event's magnitude is the M that maximises the same
likelihood over all its stations, those whose Pd exceeds the threshold and
those where it does not. The latter bound the estimate from above, so that a
small event recorded above the noise at a single station is not assigned
that station's magnitude. An event with no station above the threshold has
no estimate. The likelihood is maximised over a grid of M (M_GRID), which
ayzek evaluates the same way.

Two further estimators are fitted for comparison, both on the values above
the threshold only and both with an event's estimate the median over its
stations: the same relation by weighted least squares, inverted for M
(`ols`), and the regression of M on log10 Pd and log10 R (`direct`). The
latter minimises the error in M and so shrinks the slope by the ratio of the
magnitude variance to the total variance (regression dilution), which pulls
the largest events towards the mean magnitude.

**Channels that do not record.** A vertical channel that has failed, or
records far less than it should, reports every value as censored, and the
spurious upper bounds pull the estimates down: on the 2025-04-23 Marmara
replay, one such channel (MRMT, 25 times below its usual noise level) moved
the Mw 6.2 mainshock's estimate from M6.0 to M5.5. A value whose pre-P noise
lies more than DEAD_NOISE_FACTOR below its station's median (`station_noise`,
over all splits, since it uses no magnitude and the splits are
station-disjoint) is therefore left out, in the fit and in the estimates.

Each event carries a total weight of one in every fit, shared among its
stations. Station terms are estimated for stations with at least
MIN_STATION_ROWS values; others get zero. A station whose term exceeds
BAD_STATION_LOG10 in magnitude is reported as a probable response error and
left out of the estimates.

Evaluation is by event, on the validation and test splits, by magnitude band
and window, with the number of events that have an estimate, the bias (mean
signed error) and the mean absolute error.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .pd import PD_WINDOWS_S

DEPTH_KM = 10.0
MIN_SNR = 3.0
MIN_FIT_M = 2.0
MIN_STATION_ROWS = 20
BAD_STATION_LOG10 = 1.0
DEAD_NOISE_FACTOR = 10.0   # pre-P noise this far below the station's median: not recording
SUM_ZERO_PENALTY = 1e3
M_GRID = np.arange(1.0, 8.5001, 0.01)
B_VALUE_MC = 2.5          # completeness magnitude for the b-value estimate
DIST_KNOTS_KM = (70.0, 140.0)
M_KNOTS = (4.0,)
BANDS = ((0.0, 3.0), (3.0, 4.0), (4.0, 5.0), (5.0, 9.9))


def log_r(distance_km) -> np.ndarray:
    return np.log10(np.hypot(np.asarray(distance_km, float), DEPTH_KM))


def dist_basis(lr: np.ndarray, knots=DIST_KNOTS_KM) -> np.ndarray:
    """(n, 1 + len(knots)): log10 R and a hinge max(0, log10 R/knot) per knot,
    so that the attenuation is piecewise linear in log10 R."""
    lr = np.asarray(lr, float)
    return np.column_stack([lr] + [np.maximum(0.0, lr - np.log10(k)) for k in knots])


def mag_basis(m, knots=M_KNOTS) -> np.ndarray:
    """(..., 1 + len(knots)): M and a hinge max(0, M - knot) per knot, so that
    log10 Pd is piecewise linear in M."""
    m = np.asarray(m, float)
    return np.stack([m] + [np.maximum(0.0, m - k) for k in knots], axis=-1)


def mag_term(f: dict, m) -> np.ndarray:
    """The magnitude part of a censored fit's relation at M."""
    return mag_basis(m, f["m_knots"]) @ np.asarray(f["beta"])


def station_noise(df: pd.DataFrame) -> dict:
    """Each station's median pre-P noise level (m), its typical background."""
    d = df[df.pd_noise > 0]
    return {s: float(v) for s, v in d.groupby("station").pd_noise.median().items()}


def recording(df: pd.DataFrame, noise: dict) -> pd.Series:
    """False for a value whose pre-P noise lies more than DEAD_NOISE_FACTOR
    below its station's median: a channel that is not recording ground motion
    reports every value as censored, and those spurious upper bounds would
    pull an estimate down. Stations without a median are kept."""
    typical = df.station.map(noise)
    return (df.pd_noise > 0) & ~(df.pd_noise * DEAD_NOISE_FACTOR < typical)


def usable(df: pd.DataFrame, tau: float) -> pd.Series:
    """Pd above MIN_SNR times the pre-P noise level."""
    v = df[f"pd_{tau:g}s"]
    return (v > 0) & (v >= MIN_SNR * df["pd_noise"]) & df["magnitude"].notna()


def _rows(df: pd.DataFrame, tau: float):
    """(log10 Pd or its upper bound, censored flag, log10 R, M, event weight)."""
    above = usable(df, tau).to_numpy()
    y = np.where(above, np.log10(df[f"pd_{tau:g}s"].clip(lower=1e-30)),
                 np.log10((MIN_SNR * df["pd_noise"]).clip(lower=1e-30)))
    w = 1.0 / df.groupby("event_id").event_id.transform("size").to_numpy()
    return y, ~above, log_r(df.distance_km), df.magnitude.to_numpy(float), w


def _stations(df: pd.DataFrame):
    counts = df.groupby("station").size()
    stations = sorted(counts[counts >= MIN_STATION_ROWS].index)
    col = {s: k for k, s in enumerate(stations)}
    idx = df.station.map(col).fillna(-1).astype(int).to_numpy()
    return stations, idx


def _loglik(mu, y, cens, sigma):
    """Per-row log-likelihood and its derivatives in mu and log sigma."""
    from scipy.special import log_ndtr
    z = (y - mu) / sigma
    ll = np.where(cens, log_ndtr(z), -0.5 * z * z - np.log(sigma) - 0.5 * np.log(2 * np.pi))
    lam = np.exp(-0.5 * z * z - 0.5 * np.log(2 * np.pi) - log_ndtr(z))   # phi / Phi
    d_mu = np.where(cens, -lam / sigma, z / sigma)
    d_ls = np.where(cens, -lam * z, z * z - 1.0)
    return ll, d_mu, d_ls


def fit_censored(train: pd.DataFrame, tau: float, knots=DIST_KNOTS_KM,
                 m_knots=M_KNOTS, noise: dict | None = None) -> dict:
    """The relation by censored maximum likelihood (module docstring).
    `noise`: the stations' median noise levels (station_noise), by default
    from `train`."""
    from scipy.optimize import minimize
    noise = station_noise(train) if noise is None else noise
    d = train[(train.magnitude >= MIN_FIT_M) & train.magnitude.notna()
              & recording(train, noise)]
    y, cens, lr, m, w = _rows(d, tau)
    B, Mb = dist_basis(lr, knots), mag_basis(m, m_knots)
    nb, nm = B.shape[1], Mb.shape[1]
    stations, sidx = _stations(d)
    k = len(stations)
    has = sidx >= 0
    start = fit_ols(train, tau, knots=knots)
    x0 = np.r_[start["alpha"], start["beta"], np.zeros(nm - 1), start["gamma"],
               [start["station"].get(s, 0.0) for s in stations], np.log(0.3)]
    ib, ig, is_ = slice(1, 1 + nm), slice(1 + nm, 1 + nm + nb), slice(1 + nm + nb, 1 + nm + nb + k)

    def f(x):
        a, b, g, s, ls = x[0], x[ib], x[ig], x[is_], x[-1]
        mu = a + Mb @ b + B @ g + np.where(has, s[np.clip(sidx, 0, None)], 0.0)
        ll, d_mu, d_ls = _loglik(mu, y, cens, np.exp(ls))
        # Sum-zero station terms, as a quadratic penalty.
        val = -(w * ll).sum() + 0.5 * SUM_ZERO_PENALTY * s.sum() ** 2
        gw = -(w * d_mu)
        grad = np.empty_like(x)
        grad[0], grad[ib], grad[ig] = gw.sum(), Mb.T @ gw, B.T @ gw
        grad[is_] = (np.bincount(sidx[has], weights=gw[has], minlength=k)
                     + SUM_ZERO_PENALTY * s.sum())
        grad[-1] = -(w * d_ls).sum()
        return val, grad

    r = minimize(f, x0, jac=True, method="L-BFGS-B", options={"maxiter": 5000})
    x = r.x
    return {"tau_s": tau, "alpha": float(x[0]), "beta": list(map(float, x[ib])),
            "m_knots": list(m_knots), "gamma": list(map(float, x[ig])), "knots_km": list(knots),
            "sigma": float(np.exp(x[-1])),
            "station": dict(zip(stations, map(float, x[is_]))),
            "station_noise": noise,
            "n_values": int(len(d)), "n_censored": int(cens.sum()),
            "n_events": int(d.event_id.nunique()), "converged": bool(r.success)}


def fit_ols(train: pd.DataFrame, tau: float, inverse: bool = False,
            knots=DIST_KNOTS_KM) -> dict:
    """Weighted least squares on the values above the threshold: the relation
    (or, with `inverse`, M = a log10 Pd + b . dist_basis(log10 R) + c + s_i)."""
    d = train[usable(train, tau) & (train.magnitude >= MIN_FIT_M)]
    stations, sidx = _stations(d)
    lp, m = np.log10(d[f"pd_{tau:g}s"].to_numpy()), d.magnitude.to_numpy()
    B = dist_basis(log_r(d.distance_km), knots)
    nb = B.shape[1]
    X = np.zeros((len(d), 2 + nb + len(stations)))
    X[:, 0] = 1.0
    X[:, 1] = lp if inverse else m
    X[:, 2:2 + nb] = B
    X[np.flatnonzero(sidx >= 0), 2 + nb + sidx[sidx >= 0]] = 1.0
    y = m if inverse else lp
    w = 1.0 / d.groupby("event_id").event_id.transform("size").to_numpy()
    # Sum-zero station terms, as one heavily weighted extra equation.
    X = np.vstack([X, np.r_[np.zeros(2 + nb), np.ones(len(stations))]])
    y, w = np.r_[y, 0.0], np.r_[w, 1e3 * w.sum()]
    sw = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    out = {"tau_s": tau, "n_values": int(len(d)), "n_events": int(d.event_id.nunique()),
           "knots_km": list(knots), "station": dict(zip(stations, map(float, coef[2 + nb:])))}
    g = list(map(float, coef[2:2 + nb]))
    if inverse:
        out.update(a=float(coef[1]), b=g, c=float(coef[0]))
    else:
        out.update(alpha=float(coef[0]), beta=float(coef[1]), gamma=g)
    return out


def _flagged(df: pd.DataFrame, f: dict) -> np.ndarray:
    return df.station.map(lambda c: abs(f["station"].get(c, 0.0)) > BAD_STATION_LOG10).to_numpy()


def estimate_median(df: pd.DataFrame, f: dict, inverse: bool = False) -> pd.DataFrame:
    """Event estimates: the median of the stations' inverted values above the
    threshold (estimators `ols` and `direct`)."""
    tau = f["tau_s"]
    v = df[f"pd_{tau:g}s"]
    lp = np.log10(v.where(v > 0))
    B = dist_basis(log_r(df.distance_km), f["knots_km"])
    s = df.station.map(f["station"]).fillna(0.0)
    m = (f["a"] * lp + B @ np.asarray(f["b"]) + f["c"] + s if inverse
         else (lp - f["alpha"] - B @ np.asarray(f["gamma"]) - s) / f["beta"])
    m = m.where(usable(df, tau) & ~_flagged(df, f))
    e = df.assign(est=m).dropna(subset=["est"])
    return e.groupby("event_id").agg(magnitude=("magnitude", "first"),
                                     est=("est", "median"), n=("est", "size"))


def b_value(magnitudes, mc: float = B_VALUE_MC, dm: float = 0.1) -> float:
    """Gutenberg-Richter b-value by maximum likelihood (Aki 1965), with the
    half-bin correction for magnitudes reported to `dm`."""
    m = np.asarray(magnitudes, float)
    m = m[m >= mc]
    return float(np.log10(np.e) / (m.mean() - (mc - dm / 2)))


def estimate_censored(df: pd.DataFrame, f: dict, b: float | None = None) -> pd.DataFrame:
    """Event estimates from the censored likelihood over all of the event's
    stations, on M_GRID. With a b-value (`b`, default the fit's `b_value`),
    the posterior under the Gutenberg-Richter prior p(M) ~ 10^(-b M): its mean
    `est` and standard deviation `sd`. With b = 0, the maximum-likelihood M.
    None for an event no station records above the threshold."""
    b = f.get("b_value", 0.0) if b is None else b
    tau = f["tau_s"]
    d = df[recording(df, f.get("station_noise", {})) & ~_flagged(df, f)]
    y, cens, lr, _, _ = _rows(d, tau)
    s = d.station.map(f["station"]).fillna(0.0).to_numpy()
    base = f["alpha"] + dist_basis(lr, f["knots_km"]) @ np.asarray(f["gamma"]) + s
    out = []
    for eid, rows in d.groupby("event_id").indices.items():
        if cens[rows].all():
            continue
        mu = base[rows, None] + mag_term(f, M_GRID)[None, :]
        ll, _, _ = _loglik(mu, y[rows, None], cens[rows, None], f["sigma"])
        lp = ll.sum(0)
        if b == 0:
            est, sd = float(M_GRID[lp.argmax()]), np.nan
        else:
            lp = lp - b * np.log(10.0) * M_GRID
            post = np.exp(lp - lp.max())
            post /= post.sum()
            est = float((post * M_GRID).sum())
            sd = float(np.sqrt((post * (M_GRID - est) ** 2).sum()))
        out.append((eid, float(d.magnitude.iloc[rows[0]]), est, sd, int((~cens[rows]).sum())))
    return pd.DataFrame(out, columns=["event_id", "magnitude", "est", "sd", "n"]).set_index("event_id")


def calibration(df: pd.DataFrame, f: dict) -> pd.DataFrame:
    """How well the fitted model describes the data it was fitted on, by
    magnitude and distance band: the observed share of censored values against
    the share the model predicts, and the mean residual of the values above
    the threshold against the mean the model predicts for them (the mean of a
    normal truncated at the threshold). Under a correct model the pairs agree;
    an excess residual marks a band whose Pd the relation underpredicts."""
    from scipy.stats import norm
    tau = f["tau_s"]
    d = df[(df.magnitude >= MIN_FIT_M) & recording(df, f.get("station_noise", {}))]
    y, cens, lr, m, _ = _rows(d, tau)
    s = d.station.map(f["station"]).fillna(0.0).to_numpy()
    mu = f["alpha"] + mag_term(f, m) + dist_basis(lr, f["knots_km"]) @ np.asarray(f["gamma"]) + s
    z = (np.log10(MIN_SNR * d.pd_noise.to_numpy()) - mu) / f["sigma"]
    t = d.assign(cens=cens, pred_cens=norm.cdf(z),
                 res=np.where(cens, np.nan, y - mu),
                 pred_res=np.where(cens, np.nan,
                                   f["sigma"] * norm.pdf(z) / np.clip(norm.sf(z), 1e-12, None)))
    out = []
    for col, edges in (("magnitude", [2, 2.5, 3, 3.5, 4, 4.5, 5, 7]),
                       ("distance_km", [0, 30, 60, 100, 150, 250, 600])):
        for band, g in t.groupby(pd.cut(t[col], edges), observed=True):
            out.append({"tau_s": tau, "by": col, "band": str(band), "n": len(g),
                        "censored": g.cens.mean(), "pred_censored": g.pred_cens.mean(),
                        "residual": g.res.mean(), "pred_residual": g.pred_res.mean()})
    return pd.DataFrame(out)


def band_table(ev: pd.DataFrame) -> list[dict]:
    rows = []
    for lo, hi in BANDS + ((0.0, 9.9),):
        g = ev[(ev.magnitude >= lo) & (ev.magnitude < hi)]
        err = g.est - g.magnitude
        name = "all" if (lo, hi) == (0.0, 9.9) else (f">={lo:g}" if hi > 9 else f"{lo:g}-{hi:g}")
        rows.append({"band": name, "n": int(len(g)),
                     "bias": float(err.mean()) if len(g) else np.nan,
                     "mae": float(err.abs().mean()) if len(g) else np.nan})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(prog="onset fit-pd", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pd", required=True, help="table from onset measure-pd")
    ap.add_argument("--out", required=True, help="directory for fit.json and eval.csv")
    a = ap.parse_args(argv)

    df = pd.read_csv(a.pd)
    train = df[df.split == "train"]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    fits, rows, cal = {}, [], []
    b = b_value(train.drop_duplicates("event_id").magnitude)
    # Every station's median noise level, from all splits: it uses no
    # magnitude, and the splits are station-disjoint, so the validation and
    # test stations (and those ayzek runs on) need their own.
    noise = station_noise(df)
    print(f"Gutenberg-Richter b-value of the train events above M{B_VALUE_MC:g}: {b:.2f}")
    for tau in PD_WINDOWS_S:
        cen = fit_censored(train, tau, noise=noise)
        ols, inv = fit_ols(train, tau), fit_ols(train, tau, True)
        cen["b_value"] = b
        fits[f"{tau:g}"] = cen
        cal.append(calibration(train, cen))
        bad = sorted(s for s, v in cen["station"].items() if abs(v) > BAD_STATION_LOG10)
        print(f"\ntau {tau:g} s, censored ML: log10 Pd = {cen['alpha']:.3f} + magnitude "
              f"{' '.join(f'{b:+.3f}' for b in cen['beta'])} + distance {' '.join(f'{g:+.3f}' for g in cen['gamma'])} + s_i, "
              f"sigma {cen['sigma']:.2f}   "
              f"({cen['n_events']} events, {cen['n_values']} values, {cen['n_censored']} censored"
              f"{'' if cen['converged'] else ', NOT CONVERGED'})"
              f"   [OLS above threshold: beta {ols['beta']:.3f}]"
              + (f"\n  flagged stations: {', '.join(bad)}" if bad else ""))
        for split in ("val", "test"):
            part = df[df.split == split]
            for name, ev in (("censored", estimate_censored(part, cen)),
                             ("censored_ml", estimate_censored(part, cen, b=0.0)),
                             ("ols", estimate_median(part, ols)),
                             ("direct", estimate_median(part, inv, inverse=True))):
                for r in band_table(ev):
                    rows.append({"tau_s": tau, "split": split, "estimator": name, **r})
        t = pd.DataFrame([r for r in rows if r["tau_s"] == tau and r["split"] == "test"])
        for band, g in t.groupby("band", sort=False):
            print(f"  test {band:>5}: " + "   ".join(
                f"{x.estimator} n {x.n} bias {x.bias:+.2f} MAE {x.mae:.2f}" for x in g.itertuples()))
    (out / "fit.json").write_text(json.dumps(
        {"relation": "log10 Pd = alpha + beta . [M, max(0, M - k) for k in m_knots] + gamma . "
                     "[lr, max(0, lr - log10 k) for k in knots_km] + s_i, "
                     "lr = log10 sqrt(D^2 + depth^2)",
         "estimator": "censored maximum likelihood; event M the posterior mean under a "
                      "Gutenberg-Richter prior with b_value",
         "depth_km": DEPTH_KM, "min_snr": MIN_SNR, "bad_station_log10": BAD_STATION_LOG10,
         "m_grid": [float(M_GRID[0]), float(M_GRID[-1]), 0.01], "windows": fits}, indent=1))
    pd.DataFrame(rows).to_csv(out / "eval.csv", index=False)
    cal = pd.concat(cal)
    cal.to_csv(out / "calibration.csv", index=False)
    print("\ncalibration on the train split, tau 3 s (observed against predicted):")
    print(cal[cal.tau_s == 3.0].drop(columns="tau_s").round(3).to_string(index=False))
    print(f"\n-> {out}/fit.json, eval.csv, calibration.csv")
