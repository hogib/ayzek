"""Exports detector and picker weights, filter coefficients, the station table
and reference outputs for the C++ tests.

Output files (AYZW format, docs/impl/02-weights.md):

  models/detector_s{42,43,44}.ayzw   the 3-seed per-window 6 s detector
  models/spicker.ayzw                the sphase `wave` S/P picker
  models/bandpass.ayzw               butter(4, [1, 45] Hz) as b, a, and lfilter_zi
  models/stations.csv                code, lat, lon, elevation of AFAD stations
  data/fixtures/dsp.ayzw             raw -> cleaned -> standardised windows (scipy)
  data/fixtures/detector.ayzw        inputs, per-seed logits, intermediates (torch)
  data/fixtures/picker.ayzw          inputs, logits, intermediates (torch)

Run from the ayzek root. Requires data/demo/DEMI.mseed (tools/make_demo_data.py):

    uv run --project tools python tools/export_models.py \
        --detector-dir DIR --picker CKPT [--stations-csv CSV]

--detector-dir holds the three seed checkpoints of one detector run
(`best_..._seed42.pth` and so on, as `cascade-impl detect` writes them);
--picker is an sphase `wave` checkpoint. The model classes are in
tools/reference/. models/stations.csv is tracked, so --stations-csv is only
needed to rebuild it from an AFAD station list (Code, Latitude, Longitude,
Height).
"""
import argparse
import hashlib
import re
import sys
from pathlib import Path

import numpy as np
import torch
from obspy import UTCDateTime, read
from scipy import signal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ayzw import write_ayzw  # noqa: E402
from reference.conditioning import clean_block, taper_vector  # noqa: E402
from reference.detector import WaveformDetector  # noqa: E402
from reference.picker import PhasePicker  # noqa: E402

DEMO = Path("data/demo")
FS, FMIN, FMAX = 100.0, 1.0, 45.0

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def seed_checkpoints(ckpt_dir):
    """(seed, path) for every seed checkpoint of the one run in `ckpt_dir`.

    The seeds are averaged, so they must come from the same run: the filename
    minus its `_pid<n>_seed<n>.pth` tail has to agree.
    """
    tail = re.compile(r"_pid\d+_seed(\d+)\.pth$")
    found = [(int(m.group(1)), p) for p in sorted(Path(ckpt_dir).glob("*.pth"))
             if (m := tail.search(p.name))]
    if not found:
        sys.exit(f"no *_pid<n>_seed<n>.pth checkpoints in {ckpt_dir}")
    runs = {tail.sub("", p.name) for _, p in found}
    if len(runs) > 1:
        sys.exit(f"{ckpt_dir} holds checkpoints of several runs: {sorted(runs)}")
    return sorted(found)


def state_arrays(sd):
    return {k: v.detach().cpu().numpy() for k, v in sd.items() if v.dtype.is_floating_point}


def component_trace(st, comp, t):
    """The contiguous trace of one component covering epoch `t`."""
    for tr in st:
        if tr.stats.channel.endswith(comp) and tr.stats.starttime <= t < tr.stats.endtime:
            return tr
    raise LookupError(f"no {comp} trace covers {t}")


def cut(st, t, n):
    """(n, 3) float64 Z/N/E at epoch t, plus the absolute start position in samples."""
    out = np.empty((n, 3))
    pos = None
    for k, comp in enumerate("ZNE"):
        tr = component_trace(st, comp, t)
        i = int(round((t - tr.stats.starttime) * FS))
        start_pos = int(round(tr.stats.starttime.ns / 1e7)) + i
        if pos is not None and start_pos != pos:
            raise ValueError("components are not aligned")
        pos = start_pos
        out[:, k] = tr.data[i:i + n]
    return out, pos


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--detector-dir", type=Path, required=True)
    ap.add_argument("--picker", type=Path, required=True)
    ap.add_argument("--stations-csv", type=Path)
    args = ap.parse_args()
    torch.set_grad_enabled(False)
    models = Path("models")

    # --- detector ---------------------------------------------------------
    print("detector")
    dets = []
    for seed, c in seed_checkpoints(args.detector_dir):
        m = WaveformDetector(3, hidden=48, fusion_dim=96, branch1d="cnn-lstm")
        m.load_state_dict(torch.load(c, weights_only=True, map_location="cpu"), strict=True)
        m.eval()
        dets.append((seed, m))
        write_ayzw(models / f"detector_s{seed}.ayzw", state_arrays(m.state_dict()),
                   {"model": "detector-cnn-lstm", "window": 600, "channels": "ZNE",
                    "standardize": "perwindow", "transform": "asinh", "seed": seed,
                    "source": str(c), "sha256_16": sha(c)})

    # --- picker -----------------------------------------------------------
    print("picker")
    ck = torch.load(args.picker, weights_only=False, map_location="cpu")
    picker = PhasePicker(n_chunks=250)
    picker.load_state_dict(ck["state"], strict=True)
    picker.eval()
    write_ayzw(models / "spicker.ayzw", state_arrays(picker.state_dict()),
               {"model": "sphase-wave", "window": 6000, "n_chunks": 250, "channels": "ZNE",
                "standardize": "per-trace z-score", "source": str(args.picker),
                "sha256_16": sha(args.picker)})

    # --- filter and stations ----------------------------------------------
    b, a = signal.butter(4, [FMIN, FMAX], btype="bandpass", fs=FS)
    write_ayzw(models / "bandpass.ayzw", {"b": b, "a": a, "zi": signal.lfilter_zi(b, a)},
               {"design": "scipy.signal.butter(4, [1, 45], 'bandpass', fs=100)"})
    if args.stations_csv:
        import csv
        rows = list(csv.DictReader(open(args.stations_csv, encoding="utf-8-sig")))
        with open(models / "stations.csv", "w") as f:
            f.write("code,lat,lon,elev_m\n")
            for r in rows:
                f.write(f"{r['Code']},{r['Latitude']},{r['Longitude']},{r['Height']}\n")
        print(f"  {models / 'stations.csv'}  ({len(rows)} stations)")

    # --- fixtures from real data ------------------------------------------
    print("fixtures")
    st = read(str(DEMO / "DEMI.mseed"))
    st.merge(method=1, fill_value=None)
    st = st.split()
    # Windows: noise, the P onset of the 2025-11-10 18:20:51 M4.9 at several offsets, S coda.
    starts = ["2025-11-10T18:10:00", "2025-11-10T18:20:54", "2025-11-10T18:20:56",
              "2025-11-10T18:20:58", "2025-11-10T18:21:05", "2025-11-10T18:28:00"]
    raw, pos = [], []
    for s in starts:
        w, p = cut(st, UTCDateTime(s), 600)
        raw.append(w)
        pos.append(p)
    raw = np.stack(raw)                                          # (k, 600, 3)
    taper = taper_vector(600)
    cleaned = np.stack([clean_block(raw[:, :, c], FS, FMIN, FMAX, taper) for c in range(3)], axis=-1)
    mu = cleaned.mean(axis=1, keepdims=True)
    sd = np.maximum(cleaned.std(axis=1, keepdims=True), 1e-12)
    standardized = ((cleaned - mu) / sd).astype(np.float32)
    write_ayzw(Path("data/fixtures/dsp.ayzw"),
               {"raw": raw, "cleaned": cleaned, "standardized": standardized,
                "position": np.array(pos, dtype=np.int64)},
               {"source": "data/demo/DEMI.mseed", "starts": starts})

    x = torch.asinh(torch.from_numpy(standardized))
    out = {"input": x.numpy()}
    probs = []
    for seed, m in dets:
        logit = m(x).squeeze(1)
        out[f"logit_s{seed}"] = logit.numpy()
        probs.append(torch.sigmoid(logit))
    out["prob"] = torch.stack(probs).mean(0).numpy()
    # Intermediate activations of the first seed, for layer-by-layer comparison.
    m = dets[0][1]
    h = m.b1.conv(x.transpose(1, 2)).transpose(1, 2)
    out["conv"] = h.numpy()
    l, _ = m.b1.lstm(h)
    out["lstm"] = l.numpy()
    a_, _ = m.b1.attn(l, l, l)
    out["attn"] = a_.numpy()
    out["pooled"] = m.b1.norm(l + a_).mean(1).numpy()
    write_ayzw(Path("data/fixtures/detector.ayzw"), out, {"intermediates_seed": dets[0][0]})
    print("  ensemble probabilities:", np.round(out["prob"], 4).tolist())

    praw = np.stack([cut(st, UTCDateTime(s), 6000)[0] for s in
                     ["2025-11-10T18:20:55", "2025-11-10T18:08:00"]])     # (k, 6000, 3)
    pc = np.stack([clean_block(praw[:, :, c], FS, FMIN, FMAX, taper_vector(6000))
                   for c in range(3)], axis=-1)
    pc = pc - pc.mean(axis=1, keepdims=True)
    psd = pc.std(axis=1, keepdims=True)
    pin = torch.from_numpy((pc / np.where(psd > 0, psd, 1.0)).astype(np.float32)).transpose(1, 2)
    logits = picker(pin)
    w = picker.wave
    stem = w.stem(pin)
    pooled = torch.nn.functional.adaptive_avg_pool1d(stem, 250).transpose(1, 2)
    pl, _ = w.lstm(pooled)
    write_ayzw(Path("data/fixtures/picker.ayzw"),
               {"input": pin.numpy(), "logits": logits.numpy(), "stem": stem.numpy(),
                "pooled": pooled.numpy(), "lstm": pl.numpy(), "wave": w(pin).numpy()},
               {"starts": ["2025-11-10T18:20:55", "2025-11-10T18:08:00"]})
    prob = torch.softmax(logits, -1)
    for i in range(len(pin)):
        p_idx, s_idx = int(prob[i, :, 1].argmax()), int(prob[i, :, 2].argmax())
        print(f"  picker window {i}: P at {(p_idx + .5) * 24 / 100:.2f} s (p={float(prob[i, p_idx, 1]):.2f}), "
              f"S at {(s_idx + .5) * 24 / 100:.2f} s (p={float(prob[i, s_idx, 2]):.2f})")


if __name__ == "__main__":
    main()
