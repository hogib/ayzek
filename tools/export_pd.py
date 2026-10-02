"""Exports onset's Pd magnitude relation and the stations' gains for ayzek.

    uv run --project tools python tools/export_pd.py \\
        --fit ../../onset/runs/pd_v4/fit.json \\
        --inventory ~/Projects/sismokaos/data_downloader/raw/data/station_inventory

Writes (docs/impl/17-pd-magnitude.md):

  models/pd_relation.csv    one row per Pd window: the relation's coefficients,
                            its knots, residual sd, b-value, threshold and
                            fixed depth, and the high-pass second-order
                            section of the displacement chain
  models/pd_gain.csv        each station's vertical channel and its
                            sensitivity (counts per m/s) for every epoch of
                            its StationXML response
  models/pd_terms.csv       the station terms per window
  models/pd_noise.csv       each station's median pre-P noise level: a value
                            whose noise lies more than dead_noise_factor below
                            it is taken not to record, more than
                            coda_noise_factor above it to lie in a coda
  data/fixtures/pd_chain.csv, pd_estimate.csv
                            reference outputs of onset's pd.py and pd_fit.py
                            for tests/test_pd.cpp

The relation is fitted in onset (`onset fit-pd`, onset's docs/MAGNITUDE.md);
its code is copied into tools/reference/onset/ (pd.py, pd_fit.py) so that the
fixtures come from the code the relation was fitted with.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reference.onset import pd as pdm          # noqa: E402
from reference.onset import pd_fit             # noqa: E402

BANDS = ("HH", "BH", "EH")


def gains(inventory: Path, stations: set[str]) -> list[dict]:
    """Every epoch of each station's preferred vertical channel."""
    from obspy import read_inventory
    rows = []
    for f in sorted(inventory.glob("*.xml")):
        for net in read_inventory(str(f)):
            for sta in net:
                if sta.code not in stations:
                    continue
                for band in BANDS:
                    chans = [c for c in sta.channels if c.code == band + "Z"
                             and c.response is not None
                             and c.response.instrument_sensitivity is not None]
                    if not chans:
                        continue
                    best_loc = min(c.location_code for c in chans)
                    for c in chans:
                        if c.location_code != best_loc:
                            continue
                        s = c.response.instrument_sensitivity
                        if (s.input_units or "").upper().replace(" ", "") != "M/S":
                            continue
                        rows.append({"station": sta.code, "channel": c.code,
                                     "start": c.start_date.timestamp if c.start_date else 0.0,
                                     "end": c.end_date.timestamp if c.end_date else "",
                                     "sensitivity": float(s.value)})
                    break
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", type=Path, required=True, help="onset fit-pd's fit.json")
    ap.add_argument("--inventory", type=Path, required=True, help="StationXML directory")
    ap.add_argument("--stations", type=Path, default=Path("models/stations.csv"))
    a = ap.parse_args()

    fit = json.loads(a.fit.read_text())
    windows = [fit["windows"][k] for k in sorted(fit["windows"], key=float)]
    for w in windows:
        if len(w["beta"]) != 2 or len(w["gamma"]) != 3:
            sys.exit("ayzek reads one magnitude knot and two distance knots")
    sos = np.asarray(pdm.highpass_sos(pdm.FS))
    if sos.shape != (1, 6):
        sys.exit("ayzek reads a single second-order section")
    b0, b1, b2, a0, a1, a2 = sos[0]

    out = Path("models")
    with open(out / "pd_relation.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["tau_s", "alpha", "beta_m", "beta_above", "m_knot", "gamma_r", "gamma_k1",
                    "gamma_k2", "k1_km", "k2_km", "sigma", "b_value", "min_snr", "depth_km",
                    "flag_log10", "dead_noise_factor", "coda_noise_factor", "hp_b0", "hp_b1",
                    "hp_b2", "hp_a1", "hp_a2"])
        for x in windows:
            w.writerow([x["tau_s"], x["alpha"], *x["beta"], x["m_knots"][0], *x["gamma"],
                        *x["knots_km"], x["sigma"], x["b_value"], fit["min_snr"],
                        fit["depth_km"], fit["bad_station_log10"], pd_fit.DEAD_NOISE_FACTOR,
                        pd_fit.CODA_NOISE_FACTOR,
                        b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0])
    with open(out / "pd_terms.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["station", "tau_s", "term"])
        for x in windows:
            for s, t in sorted(x["station"].items()):
                w.writerow([s, x["tau_s"], t])
    with open(out / "pd_noise.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["station", "noise_m"])
        for s, v in sorted(windows[0]["station_noise"].items()):
            w.writerow([s, v])
    codes = {r["code"] for r in csv.DictReader(open(a.stations, encoding="utf-8-sig"))}
    g = gains(a.inventory, codes)
    with open(out / "pd_gain.csv", "w", newline="") as f:
        w = csv.DictWriter(f, ["station", "channel", "start", "end", "sensitivity"],
                           lineterminator="\n")
        w.writeheader()
        w.writerows(g)
    flagged = sorted({s for x in windows for s, t in x["station"].items()
                      if abs(t) > fit["bad_station_log10"]})
    print(f"  models/pd_relation.csv  {len(windows)} windows")
    print(f"  models/pd_terms.csv     {sum(len(x['station']) for x in windows)} terms; "
          f"flagged: {', '.join(flagged)}")
    print(f"  models/pd_noise.csv     {len(windows[0]['station_noise'])} stations")
    print(f"  models/pd_gain.csv      {len(g)} epochs, {len({r['station'] for r in g})} "
          f"of {len(codes)} stations")

    # --- fixtures -----------------------------------------------------------
    fx = Path("data/fixtures")
    fx.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(7)
    n, p = 9000, 6000
    t = np.arange(n) / pdm.FS
    counts = 1e5 + 300 * rng.normal(size=n) + 2e3 * np.sin(2 * np.pi * 0.2 * t)
    counts[p:] += 4e4 * np.exp(-(t[p:] - t[p]) / 4.0) * np.sin(2 * np.pi * 3.0 * (t[p:] - t[p]))
    counts = np.round(counts)
    sens = 1.5e9
    d = pdm.displacement(counts, sens, "M/S")
    peaks = pdm.peak_displacement(d, p)
    with open(fx / "pd_chain.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["counts", "displacement"])
        w.writerows(zip(counts.astype(np.int64), d))
    with open(fx / "pd_chain_peaks.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["sensitivity", "p_index"] + list(peaks))
        w.writerow([sens, p] + list(peaks.values()))
    # The event estimator on hand-made cases, at the 5 s window.
    rel = fit["windows"]["5"]
    stations = sorted(rel["station"])[:4]
    cases = [
        [(stations[0], 30.0, 2e-5, 1e-7), (stations[1], 80.0, 3e-6, 1e-7),
         (stations[2], 120.0, 1e-6, 2e-7)],
        [(stations[0], 40.0, 5e-7, 1e-7), (stations[1], 90.0, 2e-7, 1.5e-7),
         (stations[3], 150.0, 1e-7, 2e-7)],
        [(stations[2], 25.0, 4e-4, 1e-7), (stations[3], 60.0, 9e-5, 1e-7)],
        # the last station's noise is far below its median: left out
        [(stations[0], 30.0, 2e-5, 1e-7), (stations[1], 80.0, 3e-6, 1e-7),
         (stations[2], 60.0, 2e-10, 1e-10)],
        # the last station's noise is far above its median (a coda): left out
        [(stations[0], 30.0, 2e-5, 1e-7), (stations[1], 80.0, 3e-6, 1e-7),
         (stations[2], 60.0, 5e-3, 1e-4)],
    ]
    import pandas as pd
    with open(fx / "pd_estimate.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["case", "station", "distance_km", "pd", "pd_noise", "est", "sd"])
        for k, rows in enumerate(cases):
            df = pd.DataFrame([{"event_id": k, "station": s, "distance_km": r, "magnitude": 0.0,
                                "pd_5s": v, "pd_noise": nz} for s, r, v, nz in rows])
            e = pd_fit.estimate_censored(df, rel)
            est, sd = (float(e.est.iloc[0]), float(e.sd.iloc[0])) if len(e) else (np.nan, np.nan)
            for s, r, v, nz in rows:
                w.writerow([k, s, r, v, nz, est, sd])
    print(f"  {fx}/pd_chain.csv, pd_chain_peaks.csv, pd_estimate.csv")


if __name__ == "__main__":
    main()
