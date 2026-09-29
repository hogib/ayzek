"""Exports the onset transformer detector and reference outputs for its C++ test.

    uv run --project tools python tools/export_transformer.py \\
        --run ../../onset/runs/fdsn_v1 [--checkpoint best.pt] [--dt-reset 1,5]

Writes (AYZW, docs/impl/02-weights.md):

  models/transformer.ayzw          weights, geometry, the causal band-pass as
                                   second-order sections, the validation
                                   operating threshold and the default dt
                                   restart (--dt-reset); `geo_head.*` too for a
                                   run trained with --geometry 1, which makes
                                   ayzek locate from it (--locate geometry)
  data/fixtures/transformer.ayzw   PyTorch/scipy outputs on DEMI around the
                                   2025-11-10 M4.9, for tests/test_transformer.cpp

The model code is a copy of onset's in `tools/reference/onset/` (provenance in
its `__init__.py`), so the export needs only the run directory, not an onset
checkout. The copy must match the code the run was trained with: re-copy it
when onset's model, conditioning or filter changes.
Run from the ayzek root; needs data/demo/DEMI.mseed.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import torch
from obspy import UTCDateTime, read
from scipy import signal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ayzw import write_ayzw  # noqa: E402

FS = 100.0
RELEASE_RATIO = 0.5          # onset.metrics: re-arm below half the threshold


def load_onset():
    from reference.onset import conditioning, dsp  # noqa: F401
    from reference.onset.config import load_run_config
    from reference.onset.model import OnsetDetector
    return conditioning, dsp, load_run_config, OnsetDetector


def geometry(out):
    """(T, 2) float32: log distance and its sd, as ayzek's `Geometry`."""
    return torch.stack([out["log_dist"][0].float(),
                        torch.exp(0.5 * out["log_dist_var"][0].float())], dim=1).numpy()


def demi(t0, seconds):
    """(n, 3) float64 Z N E counts and (n, 3) missing, from data/demo/DEMI.mseed."""
    st = read("data/demo/DEMI.mseed")
    n = int(seconds * FS)
    x, m = np.zeros((n, 3)), np.ones((n, 3), bool)
    for c, comp in enumerate("ZNE"):
        s = st.select(component=comp).copy()
        s.merge(method=0, fill_value=None)
        tr = s[0]
        tr.trim(t0, t0 + (n - 1) / FS, pad=True, nearest_sample=True, fill_value=None)
        d = np.ma.asarray(tr.data, dtype=np.float64)
        x[:, c] = d.filled(0.0)[:n]
        m[:, c] = np.ma.getmaskarray(d)[:n]
    return x, m


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", type=Path, required=True, help="onset run directory")
    ap.add_argument("--checkpoint", default="best.pt")
    ap.add_argument("--threshold", type=float, default=None,
                    help="Operating threshold; default: the run's val_best.json.")
    ap.add_argument("--dt-reset", default=None, metavar="BELOW,FROM",
                    help="ayzek's default dt restart for this model; default: the "
                         "run's training rule (dt_reset_below, dt_reset_from)")
    a = ap.parse_args()
    conditioning, dsp, load_run_config, OnsetDetector = load_onset()
    torch.set_grad_enabled(False)

    mcfg, dcfg, tcfg, _ = load_run_config(a.run)
    ckpt = a.run / a.checkpoint
    model = OnsetDetector(mcfg).eval()
    model.load_state_dict(torch.load(ckpt, map_location="cpu", weights_only=True))
    best = json.loads((a.run / "val_best.json").read_text())
    thr = a.threshold if a.threshold is not None else best["summary"]["threshold"]

    # --- weights ------------------------------------------------------------
    tensors = {k: v.numpy().astype(np.float32) for k, v in model.state_dict().items()
               if not k.startswith("context.stem.")}          # the shared stem, once
    tensors["config.ints"] = np.array([mcfg.d_model, mcfg.n_heads, mcfg.n_layers,
                                       mcfg.window_tokens, mcfg.context_tokens,
                                       mcfg.context_layers], np.int64)
    tensors["config.stem"] = np.array(mcfg.stem, np.int64)
    tensors["config.floats"] = np.array([mcfg.max_dt_s, dcfg.ctx_seconds,
                                         dcfg.fallback_scale_s, thr, RELEASE_RATIO,
                                         mcfg.sample_rate], np.float64)
    if mcfg.geometry:
        # The geometry head's output clamp (onset model.geometry): the
        # log-variance of log distance.
        tensors["config.geo"] = np.array([2 * np.log(getattr(mcfg, "geo_min_sd", np.exp(-4.0))),
                                          6.0], np.float64)
    # ayzek's default dt restart: [below, from, tokens] (trigger.hpp).
    below, frm = (map(float, a.dt_reset.split(",")) if a.dt_reset else
                  (getattr(tcfg, "dt_reset_below", 2.0), getattr(tcfg, "dt_reset_from", 5.0)))
    tensors["config.trigger"] = np.array([below, frm, getattr(tcfg, "dt_reset_tokens", 2)],
                                         np.float64)
    sos = dsp.bandpass_sos(FS)
    tensors["filter.sos"] = np.asarray(sos, np.float64)
    tensors["filter.zi"] = signal.sosfilt_zi(sos).astype(np.float64)
    write_ayzw(Path("models/transformer.ayzw"), tensors, {
        "model": "onset-transformer", "source": str(ckpt),
        "sha256_16": hashlib.sha256(ckpt.read_bytes()).hexdigest()[:16],
        "threshold": thr, "dt_reset": [below, frm], "val_epoch": best["epoch"],
        "val_recall_1s": best["summary"]["recall@1.0s"],
        "val_false_per_hour": best["summary"]["false_per_hour"],
        "geometry": bool(mcfg.geometry),
        "filter": "butter(4, [1, 45], bandpass, fs=100) as sos, causal, restart after gaps"})

    # --- fixtures -----------------------------------------------------------
    fx = {}
    # 1. The filter, on the vertical with a 2 s gap and a 0.5 s run after it
    #    (too short to filter, so missing), then data again.
    raw, miss = demi(UTCDateTime("2025-11-10T18:20:30"), 40.0)
    z, zm = raw[:, 0].copy(), miss[:, 0].copy()
    zm[1500:1700] = True
    zm[1750:1800] = True
    y, ym = dsp.causal_filter(z, zm, FS)
    fx.update({"filter.raw": z, "filter.missing": zm.astype(np.uint8),
               "filter.out": y.astype(np.float64), "filter.out_missing": ym.astype(np.uint8)})

    # 2. The network on conditioned input: 30 s over the P arrival, with the
    #    minute from 18:10 as context, and again with no context.
    ctx_raw, ctx_m = demi(UTCDateTime("2025-11-10T18:10:00"), dcfg.ctx_seconds)
    cw, cm = dsp.filter_components(ctx_raw, ctx_m, FS)
    scale = conditioning.channel_scale(cw, cm)
    xr, xm = demi(UTCDateTime("2025-11-10T18:20:40"), 30.0)
    xw, xmask = dsp.filter_components(xr, xm, FS)
    x = conditioning.condition(xw, xmask, scale)
    c = conditioning.condition(cw, cm, scale)
    xt, ct = torch.from_numpy(x)[None], torch.from_numpy(c)[None]
    with_ctx = model(xt, ct, torch.tensor([True]))
    without = model(xt)
    fx.update({"net.x": x, "net.ctx": c,
               "net.p_ctx": torch.sigmoid(with_ctx["logit"])[0].numpy(),
               "net.dt_ctx": with_ctx["dt"][0].numpy(),
               "net.p_null": torch.sigmoid(without["logit"])[0].numpy(),
               "net.dt_null": without["dt"][0].numpy()})
    if mcfg.geometry:
        fx.update({"net.geo_ctx": geometry(with_ctx), "net.geo_null": geometry(without)})

    # 3. The whole front end: 50 s of raw counts with a 1.5 s horizontal gap,
    #    through filter, first-second scale and the null context. Shorter than
    #    the 60 s before the first context refresh, so one forward pass is the
    #    reference.
    sr, sm = demi(UTCDateTime("2025-11-10T18:20:20"), 50.0)
    sm[2000:2150, 1] = True
    sw, smask = dsp.filter_components(sr, sm, FS)
    n0 = int(dcfg.fallback_scale_s * FS)
    s_scale = conditioning.channel_scale(sw[:n0], smask[:n0])
    sx = conditioning.condition(sw, smask, s_scale)
    out = model(torch.from_numpy(sx)[None])
    fx.update({"stream.raw": sr, "stream.missing": sm.astype(np.uint8),
               "stream.p": torch.sigmoid(out["logit"])[0].numpy(),
               "stream.dt": out["dt"][0].numpy()})
    if mcfg.geometry:
        fx["stream.geo"] = geometry(out)
    write_ayzw(Path("data/fixtures/transformer.ayzw"), fx,
               {"source": "data/demo/DEMI.mseed", "checkpoint": str(ckpt)})

    p = fx["stream.p"]
    j = int(np.argmax(p >= thr)) if (p >= thr).any() else None
    print(f"  threshold {thr:.6f}; stream fixture: max p {p.max():.4f}"
          + (f", first crossing at {(j * 10 + 9) / FS:.2f} s into the 18:20:20 window, "
             f"P dated {(j * 10 + 9) / FS - fx['stream.dt'][j]:.2f} s" if j is not None else ""))
    if mcfg.geometry and j is not None:
        # DEMI and the 2025-11-10 Sındırgı M4.9 (AFAD 39.2250N 28.1715E), for a
        # sanity check of the head on data it may never have seen.
        g = fx["stream.geo"]
        sta = next(l.split(",") for l in Path("models/stations.csv").read_text().splitlines()
                   if l.startswith("DEMI,"))
        slat, slon, elat, elon = np.radians([float(sta[1]), float(sta[2]), 39.2250, 28.1715])
        true_km = 2 * 6371.0 * np.arcsin(np.sqrt(
            np.sin((elat - slat) / 2) ** 2
            + np.cos(slat) * np.cos(elat) * np.sin((elon - slon) / 2) ** 2))
        print(f"  DEMI to the AFAD epicentre: {true_km:.0f} km")
        for k in (j, min(j + 50, len(g) - 1), min(j + 100, len(g) - 1)):
            print(f"  geometry {(k - j) / 10:4.1f} s after the trigger: "
                  f"{np.exp(g[k, 0]):.0f} km (x/{np.exp(g[k, 1]):.2f})")


if __name__ == "__main__":
    main()
