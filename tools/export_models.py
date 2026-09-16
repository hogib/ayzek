"""Export weights and reference outputs for ayzek.

Writes, in the AYZW format (docs/impl/02-weights.md):

  models/detector_s{42,43,44}.ayzw   the 3-seed per-window 6 s detector
  models/spicker.ayzw                the sphase `wave` S/P picker
  models/bandpass.ayzw               butter(4, [1, 45] Hz) as b, a
  models/stations.csv                code, lat, lon for the demo stations
  data/fixtures/dsp.ayzw             raw -> cleaned -> standardised windows (scipy)
  data/fixtures/detector.ayzw        inputs, per-seed logits, intermediates (torch)
  data/fixtures/picker.ayzw          inputs, logits, intermediates (torch)

Run from the ayzek root with the archive pipeline's environment, which has
torch, scipy and obspy:

    uv run --project ~/Projects/sismokaos/archive_pipeline python tools/export_models.py
"""
import hashlib
import json
import struct
import sys
from pathlib import Path

import numpy as np
import torch
from obspy import UTCDateTime, read
from scipy import signal

SISMO = Path.home() / "Projects/sismokaos"
SPHASE = Path.home() / "Projects/Codings/sphase"
sys.path.insert(0, str(SISMO / "archive_pipeline/src"))
sys.path.insert(0, str(SPHASE / "src"))

from archive_pipeline.archive import clean_block, taper_vector  # noqa: E402
from archive_pipeline.products.detector import WaveformDetector, find_checkpoints  # noqa: E402
from sphase.model import PhasePicker  # noqa: E402

DETECTOR_DIR = SISMO / "cnn_earthquake/trained_model_perwindow_6s"
PICKER_CKPT = SPHASE / "runs/wave_n250.pt"
STATIONS_CSV = SISMO / "data_downloader/catalogs/istasyon_katalog.csv"
DEMO = Path("data/demo")
FS, FMIN, FMAX = 100.0, 1.0, 45.0

DTYPES = {np.dtype(np.float32): 0, np.dtype(np.float64): 1, np.dtype(np.int64): 2, np.dtype(np.uint8): 3}


def write_ayzw(path, tensors, meta):
    """magic, version, count, then (name, dtype, dims, bytes) per tensor."""
    path.parent.mkdir(parents=True, exist_ok=True)
    items = dict(tensors)
    items["__meta__"] = np.frombuffer(json.dumps(meta, indent=1).encode(), dtype=np.uint8)
    with open(path, "wb") as f:
        f.write(b"AYZW" + struct.pack("<II", 1, len(items)))
        for name, arr in items.items():
            arr = np.ascontiguousarray(arr)
            if arr.dtype not in DTYPES:
                arr = arr.astype(np.float32)
            nb = name.encode()
            f.write(struct.pack("<H", len(nb)) + nb)
            f.write(struct.pack("<BB", DTYPES[arr.dtype], arr.ndim))
            f.write(struct.pack(f"<{arr.ndim}I", *arr.shape))
            f.write(struct.pack("<Q", arr.nbytes))
            f.write(arr.astype(arr.dtype.newbyteorder("<")).tobytes())
    print(f"  {path}  ({len(items)} tensors, {path.stat().st_size / 1e3:.0f} kB)")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


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
    torch.set_grad_enabled(False)
    models = Path("models")

    # --- detector ---------------------------------------------------------
    print("detector")
    ckpts = find_checkpoints(DETECTOR_DIR, "1d", "linear", "cnn-lstm")
    dets = []
    for c in ckpts:
        seed = int(str(c).rsplit("seed", 1)[1].split(".")[0])
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
    ck = torch.load(PICKER_CKPT, weights_only=False, map_location="cpu")
    picker = PhasePicker(arm="wave", n_chunks=250)
    picker.load_state_dict(ck["state"], strict=True)
    picker.eval()
    write_ayzw(models / "spicker.ayzw", state_arrays(picker.state_dict()),
               {"model": "sphase-wave", "window": 6000, "n_chunks": 250, "channels": "ZNE",
                "standardize": "per-trace z-score", "source": str(PICKER_CKPT),
                "sha256_16": sha(PICKER_CKPT)})

    # --- filter and stations ----------------------------------------------
    b, a = signal.butter(4, [FMIN, FMAX], btype="bandpass", fs=FS)
    write_ayzw(models / "bandpass.ayzw", {"b": b, "a": a, "zi": signal.lfilter_zi(b, a)},
               {"design": "scipy.signal.butter(4, [1, 45], 'bandpass', fs=100)"})
    import csv
    rows = list(csv.DictReader(open(STATIONS_CSV, encoding="utf-8-sig")))
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
    # Noise, the P onset of the 2025-11-10 18:20:51 M4.9 at a few offsets, and the S coda.
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
    # Intermediates of the first seed, for localising a mismatch.
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
