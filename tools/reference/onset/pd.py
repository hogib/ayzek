"""Peak P-wave displacement (Pd) for magnitude estimation.

The detector's input carries no absolute amplitude: each component is divided
by the station's noise RMS and band-passed at 1–45 Hz (conditioning.py,
dsp.py). Both are deliberate for detection, and both remove what magnitude
depends on: the first cancels the instrument gain together with the site's
noise level, and the second removes most of the energy below 1 Hz, where
the spectra of large earthquakes differ from those of small ones. Magnitude
is therefore estimated from a separate measurement on the raw vertical
component, the peak displacement of the P wave (Wu & Kanamori 2005):

    counts / sensitivity -> ground velocity (m/s)
    -> causal 2-pole Butterworth high-pass at HP_HZ
    -> integration -> displacement (m)
    -> the same high-pass again, against integration drift
    Pd(tau) = max |d(t)| for P <= t < P + tau

Every step is causal, so the value at a time t can also be computed at t in
real time; ayzek computes the same quantity from the transformer's dated P.
Pd is measured for several windows `tau` (PD_WINDOWS_S): a magnitude
estimate then becomes available one second after P and is refined as the
window grows. The noise level `pd_noise`, the same measurement over the
NOISE_S before P, qualifies each value.

The magnitude follows from a scaling relation between log10 Pd, magnitude
and distance, piecewise linear in both and with station terms, fitted to the
stored event traces by censored maximum likelihood (pd_fit.py,
docs/MAGNITUDE.md). `onset measure-pd` writes the table it is fitted to: one
row per event trace of a store, with the catalogue magnitude, the distance,
and Pd for every window.
"""
from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

HP_HZ = 0.075
HP_ORDER = 2
PD_WINDOWS_S = (1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0)
NOISE_S = 10.0          # the noise window ends 1 s before P
SETTLE_S = 30.0         # filter run-in before the noise window
FS = 100.0
FDSN_ROOT = os.environ.get("ONSET_FDSN_ROOT", "/home/oguzb/Projects/sismokaos/data_downloader")


@lru_cache(maxsize=None)
def highpass_sos(fs: float = FS, hz: float = HP_HZ, order: int = HP_ORDER):
    from scipy import signal
    return signal.butter(order, hz, btype="highpass", fs=fs, output="sos")


def _highpass(x: np.ndarray, fs: float) -> np.ndarray:
    from scipy import signal
    sos = highpass_sos(fs)
    zi = signal.sosfilt_zi(sos) * x[0]
    return signal.sosfilt(sos, x, zi=zi)[0]


def displacement(counts: np.ndarray, sensitivity: float, units: str, fs: float = FS):
    """Raw counts of one gap-free component -> causal displacement in metres.

    `units` are the response's input units: velocity (M/S) is integrated once,
    acceleration (M/S**2) twice, each integration preceded by the high-pass.
    """
    x = (counts - counts[0]) / sensitivity
    n = {"M/S": 1, "M/S**2": 2}.get(units.upper().replace(" ", ""))
    if n is None:
        raise ValueError(f"unsupported response input units {units!r}")
    for _ in range(n):
        x = np.cumsum(_highpass(x, fs)) / fs
    return _highpass(x, fs)


def peak_displacement(d: np.ndarray, p: int, fs: float = FS) -> dict:
    """Pd for every window after sample `p`, and the noise level before it."""
    out = {f"pd_{t:g}s": float(np.abs(d[p: p + int(t * fs)]).max()) for t in PD_WINDOWS_S}
    a, b = p - int((NOISE_S + 1.0) * fs), p - int(fs)
    out["pd_noise"] = float(np.abs(d[a:b]).max())
    return out


# -- the measurement over a store --------------------------------------------

_W: dict = {}


def _init(inventory_dir: str):
    from obspy import Inventory, read_inventory
    inv = Inventory()
    for f in sorted(Path(inventory_dir).glob("*.xml")):
        inv += read_inventory(str(f))
    _W["inv"] = inv


def _response(net, sta, loc, cha, t):
    """(sensitivity counts per input unit, input units) at time t, or None."""
    sel = _W["inv"].select(network=net, station=sta, location=loc, channel=cha, time=t)
    for n in sel:
        for s in n:
            for c in s:
                r = c.response
                if r is not None and r.instrument_sensitivity is not None:
                    return (float(r.instrument_sensitivity.value),
                            r.instrument_sensitivity.input_units or "")
    return None


def _event(job):
    """All of one event file's traces: [row dicts] and {reason: count}."""
    from obspy import UTCDateTime, read
    from onset.build_fdsn import instruments, to_grid
    path, rows = job
    drops: dict = {}
    out = []

    def drop(why):
        drops[why] = drops.get(why, 0) + 1

    try:
        st = read(str(path))
    except Exception:
        drop("unreadable")
        return out, drops
    insts = instruments(st)
    for r in rows:
        comps = insts.get((r["network"], r["station"]))
        if comps is None:
            drop("no_trace")
            continue
        z = comps[0]
        t0 = UTCDateTime(r["start_time"])
        n = int(r["n_samples"])
        counts, missing = to_grid(z, t0, n)
        p = int(round(r["p_sample"]))
        a = p - int((SETTLE_S + NOISE_S + 1.0) * FS)
        b = p + int(max(PD_WINDOWS_S) * FS)
        if a < 0 or b > n:
            drop("window_short")
            continue
        if missing[a:b].any():
            drop("gap")
            continue
        tr = z[0].stats
        resp = _response(tr.network, tr.station, tr.location, tr.channel,
                         t0 + p / FS)
        if resp is None:
            drop("no_response")
            continue
        try:
            d = displacement(counts[a:b], resp[0], resp[1])
        except ValueError:
            drop("units")
            continue
        out.append({**{k: r[k] for k in ("key", "split", "network", "station", "event_id",
                                         "magnitude", "depth_km", "distance_km",
                                         "p_source")},
                    "channel": tr.channel, "sensitivity": resp[0], "units": resp[1],
                    **peak_displacement(d, p - a)})
    return out, drops


def main(argv=None):
    ap = argparse.ArgumentParser(prog="onset measure-pd", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="dataset store whose event traces to measure")
    ap.add_argument("--events-dir", help="raw event pulls; default: the store's build.json")
    ap.add_argument("--inventory", default=f"{FDSN_ROOT}/raw/data/station_inventory",
                    help="StationXML files with full responses")
    ap.add_argument("--catalog", help="for the magnitude type; default: the store's build.json")
    ap.add_argument("--out", help="default: <data>/pd.csv")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args(argv)

    root = Path(a.data)
    args = json.loads((root / "build.json").read_text())["args"]
    events_dir = Path(a.events_dir or args["events_dir"])
    idx = pd.read_csv(root / "index.csv", low_memory=False)
    ev = idx[(idx.kind == "event") & (idx.source == "fdsn") & idx.p_sample.notna()]
    jobs = [(events_dir / f"event_{int(eid)}_raw.mseed", g.to_dict("records"))
            for eid, g in ev.groupby("event_id")]
    rows, drops = [], {}
    with ProcessPoolExecutor(a.workers, initializer=_init, initargs=(a.inventory,)) as ex:
        for got, why in ex.map(_event, jobs, chunksize=16):
            rows += got
            for k, v in why.items():
                drops[k] = drops.get(k, 0) + v
    out = pd.DataFrame(rows)
    cat = pd.read_csv(a.catalog or args["catalog"], encoding="utf-8-sig")
    out["magnitude_type"] = out.event_id.map(dict(zip(cat.EventID, cat.Type)))
    path = Path(a.out) if a.out else root / "pd.csv"
    out.to_csv(path, index=False)
    print(f"{len(out)} of {len(ev)} event traces measured -> {path}")
    print("not measured: " + ", ".join(f"{k} {v}" for k, v in sorted(drops.items())))
