"""Window probabilities from a published SeisBench model, as an ayzek score cache.

    uv run --with seisbench==0.12.6 python tools/seisbench_scores.py \\
        --model eqtransformer --weights instance --mode causal \\
        --out data/runs/sb/eqt_instance_causal/marmara data/marmara/*.mseed

ayzek scores a 6 s window every 0.5 s on an absolute grid; with --scores-in DIR
it takes each window's probability from DIR/STATION.csv instead of running its
detector, and everything downstream (trigger rule, STA/LTA anchoring, picker,
association, location, catalogue scoring) is unchanged. This writes such a
cache from EQTransformer or PhaseNet, so that a published model goes through
the identical evaluation as the pipeline's own detector.

The score of the window [s, s + 6) is the model's largest output over those
6 s: the detection trace for EQTransformer, the P probability for PhaseNet.
The model sees one input per window:

- causal: the input ends at s + 6, the end of the window, so the score uses
  only data a live pipeline would have when it scores that window (60 s for
  EQTransformer, 30.01 s for PhaseNet);
- lookahead: the input is centred on the window, the way these models are
  normally applied offline, and the score is an upper bound that a real-time
  system cannot reach.

With --ahead S (causal mode only), the input ends S seconds after the window
instead, so that the window sits S seconds from the input's edge: the test of
how much of the causal mode's weakness is the input edge (referee round 5, M6a;
ayzek_paper revisions2.md R18). A live pipeline would score such a window S
seconds late, which the replay does not show: S has to be added to its alarm
times. --ahead 0 is the causal mode.

A score is written for every grid point the record covers, because ayzek falls
back to its own detector for any window missing from the cache. Each model's
own normalisation is applied per input (`annotate_batch_pre`); a model whose
weights carry a band-pass (EQTransformer "original") gets it causally in the
causal mode and zero-phase in the lookahead mode. Gaps are filled with zeros
after the mean is removed. PhaseNet's annotation blinds the outer 2.5 s of
each input; here the raw output is used, and reliability is left to the
trigger rule, as for the pipeline's own detector.
"""
import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
import seisbench.models as sbm
import torch
from obspy import read

MODELS = {"eqtransformer": sbm.EQTransformer, "phasenet": sbm.PhaseNet}
# The weight versions the candidates were selected with (held-out protocol,
# section 6): SeisBench's from_pretrained takes the latest version by default,
# which a later upload would change under the frozen method. SHA-256 of each
# .pt file in the release notes.
VERSIONS = {("eqtransformer", "original"): "3", ("eqtransformer", "stead"): "2",
            ("eqtransformer", "instance"): "2", ("eqtransformer", "ethz"): "2",
            ("phasenet", "original"): "2", ("phasenet", "stead"): "2",
            ("phasenet", "instance"): "2", ("phasenet", "ethz"): "2"}
STEP, WINDOW = 0.5, 6.0  # ayzek's grid and window, seconds


def station_array(path, fs):
    """(3, n) float32 ZNE at fs from one station file, gaps as NaN, and the
    epoch of its first sample."""
    st = read(str(path))
    st.merge(method=0, fill_value=None)
    comps = {}
    for tr in st:
        c = tr.stats.channel[-1]
        c = {"1": "N", "2": "E"}.get(c, c)
        if tr.stats.sampling_rate != fs:
            sys.exit(f"{path}: {tr.stats.sampling_rate} Hz, the model wants {fs}")
        comps[c] = tr
    if set(comps) < {"Z", "N", "E"}:
        sys.exit(f"{path}: components {sorted(comps)}, need Z, N and E")
    t0 = max(comps[c].stats.starttime.timestamp for c in "ZNE")
    t1 = min(comps[c].stats.endtime.timestamp for c in "ZNE")
    n = int(round((t1 - t0) * fs)) + 1
    out = np.full((3, n), np.nan, dtype=np.float32)
    for i, c in enumerate("ZNE"):
        tr = comps[c]
        off = int(round((t0 - tr.stats.starttime.timestamp) * fs))
        d = np.ma.filled(np.ma.masked_invalid(tr.data.astype(np.float64)), np.nan)
        out[i] = d[off:off + n]
    return out, t0


def bandpass(x, fs, args, zerophase):
    """The model's band-pass, per contiguous run of samples."""
    from scipy.signal import butter, sosfilt, sosfiltfilt
    sos = butter(args.get("corners", 2), [args["freqmin"], args["freqmax"]], "bandpass", fs=fs, output="sos")
    y = np.full_like(x, np.nan)
    for c in range(x.shape[0]):
        ok = np.isfinite(x[c])
        edges = np.flatnonzero(np.diff(np.r_[0, ok.astype(np.int8), 0]))
        for a, b in zip(edges[::2], edges[1::2]):
            seg = x[c, a:b].astype(np.float64)
            seg -= seg.mean()
            y[c, a:b] = (sosfiltfilt if zerophase else sosfilt)(sos, seg) if b - a > 30 else 0.0
    return y


def scores(model, x, t0, fs, mode, label, device, batch, ahead=0.0):
    n_in = model.in_samples
    w = int(round(WINDOW * fs))
    first = math.ceil(t0 / STEP) * STEP
    starts = np.arange(first, t0 + (x.shape[1] - w) / fs + 1e-9, STEP)
    # Sample offset of each window start, and of its model input.
    s_idx = np.round((starts - t0) * fs).astype(np.int64)
    if mode == "causal":
        in_start = s_idx + w + int(round(ahead * fs)) - n_in
    else:
        in_start = s_idx + w // 2 - n_in // 2
    pad = max(0, -in_start.min()), max(0, (in_start + n_in).max() - x.shape[1])
    xp = np.pad(np.nan_to_num(x - np.nanmean(x, axis=1, keepdims=True)), ((0, 0), pad))
    in_start = in_start + pad[0]
    lo = s_idx + pad[0] - in_start  # window start within the input
    out = np.empty(len(starts), dtype=np.float32)
    with torch.no_grad():
        for b in range(0, len(starts), batch):
            idx = in_start[b:b + batch]
            xb = np.stack([xp[:, i:i + n_in] for i in idx])
            tb = torch.from_numpy(xb).to(device)
            tb = model.annotate_batch_pre(tb, {})
            y = model(tb)
            if isinstance(y, (tuple, list)):
                y = y[label]  # EQTransformer: (detection, P, S)
            else:
                y = y[:, label]  # PhaseNet: (batch, labels, samples)
            y = y.float().cpu().numpy()
            for k in range(len(idx)):
                a = lo[b + k]
                out[b + k] = y[k, a:a + w].max()
    return starts, out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--model", choices=MODELS, required=True)
    ap.add_argument("--weights", required=True)
    ap.add_argument("--weights-version", help="SeisBench weight version (default: VERSIONS)")
    ap.add_argument("--mode", choices=("causal", "lookahead"), required=True)
    ap.add_argument("--out", required=True, help="cache directory, one STATION.csv per file")
    ap.add_argument("--ahead", type=float, default=0.0,
                    help="causal mode: seconds of input after the window's end")
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    if a.ahead and a.mode != "causal":
        ap.error("--ahead needs --mode causal")

    device = "cuda" if torch.cuda.is_available() else "cpu"
    version = a.weights_version or VERSIONS.get((a.model, a.weights))
    if version is None:
        ap.error(f"no pinned version for {a.model} {a.weights}; give --weights-version")
    model = MODELS[a.model].from_pretrained(a.weights, version_str=version).to(device).eval()
    fs = model.sampling_rate
    if a.model == "eqtransformer":
        label = 0  # detection
    else:
        label = model.labels.index("P")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for f in sorted(a.files):
        station = Path(f).stem
        x, t0 = station_array(f, fs)
        if getattr(model, "filter_args", None):
            x = bandpass(x, fs, model.filter_kwargs, zerophase=a.mode == "lookahead")
        starts, p = scores(model, x, t0, fs, a.mode, label, device, a.batch, a.ahead)
        with open(out / f"{station}.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["station", "window_start", "probability"])
            for s, v in zip(starts, p):
                w.writerow([station, f"{s:.2f}", f"{v:.6f}"])
        print(f"{station}: {len(starts)} windows, max {p.max():.3f}, median {np.median(p):.3f}")


if __name__ == "__main__":
    main()
