"""Exports the magnitude regressor weights and reference outputs for the C++ tests.

The reference inputs are computed with the encoder functions that generated
the training data (tools/reference/conditioning.py). Run from the ayzek root:

    uv run --project tools python tools/export_magnitude.py \
        --partition P0 --partition P1 --partition P2 --dataset DATASET

Each --partition is a regressor checkpoint (or a directory holding one), in
partition order; --dataset is the regressor's training set
(`dataset_magreg_fdsn_10s`), whose test split supplies three corpus windows for
the fixture. Both come from cascade_impl (`REPRODUCING.md`, recipe A).

Writes:
  models/magnitude_p{0,1,2}.ayzw   the three station- and event-disjoint partition models
  data/fixtures/magnitude.ayzw     noise baselines, spectrograms, normalised
                                   inputs and model outputs on real windows
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from obspy import UTCDateTime, read

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ayzw import write_ayzw  # noqa: E402
from reference.conditioning import SpectrogramEncoder, clean_and_filter_1d, standardize  # noqa: E402
from reference.magnitude import DualChannelNet  # noqa: E402

FS, FMIN, FMAX, N_FFT, HOP, TOP_DB, WINDOW = 100.0, 1.0, 45.0, 128, 32, 80.0, 10.0


def build(path):
    m = DualChannelNet(3, 3, aux_dim=0, hidden=64, fusion_dim=128, dropout=0.3, channels="1d+2d",
                       fusion="linear", n_classes=1, squeeze_output=True, branch1d="lstm")
    m.load_state_dict(torch.load(path, weights_only=True, map_location="cpu"), strict=True)
    return m.eval()


def cut(st, t, n):
    out = np.empty((n, 3))
    for k, comp in enumerate("ZNE"):
        tr = next(tr for tr in st if tr.stats.channel.endswith(comp) and tr.stats.starttime <= t < tr.stats.endtime)
        i = int(round((t - tr.stats.starttime) * FS))
        out[:, k] = tr.data[i:i + n]
    return out


def seq_transform(path):
    """How the waveform channel was fed in training, from the run's model.json.

    The regressor's own spec records it under `protocol`; a checkpoint trained
    before the flag existed has none, which means "none". The C++ estimator
    reads the same value out of the weights file and applies it, so a model
    cannot be scored with the wrong input convention.
    """
    spec = (path if path.is_dir() else path.parent) / "model.json"
    if not spec.exists():
        return "none"
    return json.loads(spec.read_text()).get("protocol", {}).get("seq_transform", "none")


def checkpoint(path):
    """A checkpoint file, or the single one inside a directory."""
    if path.is_file():
        return path
    found = sorted(path.glob("*.pth"))
    if len(found) != 1:
        sys.exit(f"{path}: expected one *.pth, found {len(found)}")
    return found[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--partition", type=Path, action="append", required=True,
                    help="checkpoint or directory, once per partition, in order")
    ap.add_argument("--dataset", type=Path, required=True)
    args = ap.parse_args()
    DATASET = args.dataset
    torch.set_grad_enabled(False)
    enc = SpectrogramEncoder(n_fft=N_FFT, hop_length=HOP, top_db=TOP_DB, nominal_fs=FS,
                             window_seconds=WINDOW, normalize="station")
    _, spec_tf, db_tf = enc._transforms()

    transforms = {seq_transform(d) for d in args.partition}
    if len(transforms) != 1:
        sys.exit(f"partitions disagree on seq_transform: {sorted(transforms)}")
    transform = transforms.pop()
    print(f"seq_transform: {transform}")

    models = {}
    for p, d in enumerate(args.partition):
        ckpt = checkpoint(d)
        models[p] = build(ckpt)
        sd = {k: v.numpy() for k, v in models[p].state_dict().items() if v.dtype.is_floating_point}
        write_ayzw(Path(f"models/magnitude_p{p}.ayzw"), sd,
                   {"model": "magnitude-1d+2d-lstm", "window": 1000, "pre_p_seconds": 2.0,
                    "seq_transform": transform,
                    "spectrogram": {"n_fft": N_FFT, "hop": HOP, "top_db": TOP_DB},
                    "normalize": "station noise sigma (seq), station median noise dB profile (img)",
                    "partition": p, "source": str(ckpt), "reported_mae": "0.4203 +- 0.0165 over p0-p2"})

    # --- noise baseline, computed as in the encoder's baseline functions --------------
    st = read("data/demo/DEMI.mseed")
    st.merge(method=1, fill_value=None)
    st = st.split()
    noise_starts = [f"2025-11-10T18:1{m}:{s:02d}" for m in (0, 1) for s in (0, 30)] + \
                   ["2025-11-10T18:12:00", "2025-11-10T18:12:30"]
    noise_raw = np.stack([cut(st, UTCDateTime(t), 1000) for t in noise_starts])     # (6, 1000, 3)
    s_acc = np.zeros(3)
    ss_acc = np.zeros(3)
    n_acc = np.zeros(3)
    frames = {c: [] for c in range(3)}
    for w in noise_raw:
        for c in range(3):
            cl = clean_and_filter_1d(w[:, c].copy(), FS, FMIN, FMAX)
            s_acc[c] += cl.sum()
            ss_acc[c] += (cl ** 2).sum()
            n_acc[c] += cl.size
            t = torch.from_numpy(np.ascontiguousarray(cl)).float().unsqueeze(0)
            frames[c].append(db_tf(spec_tf(t))[0])                                   # (freq, time)
    mu = s_acc / n_acc
    sigma = np.sqrt(np.maximum(ss_acc / n_acc - mu ** 2, 0))
    profiles = np.stack([torch.cat(frames[c], dim=1).median(dim=1).values.numpy() for c in range(3)])

    # --- model inputs for an event window and a noise window (SpectrogramDualEncoder steps) ---
    starts = ["2025-11-10T18:20:59.80", "2025-11-10T18:16:30.00"]     # P at DEMI - 2 s; quiet
    raw = np.stack([cut(st, UTCDateTime(t), 1000) for t in starts])
    cleaned = np.stack([np.stack([clean_and_filter_1d(w[:, c].copy(), FS, FMIN, FMAX) for c in range(3)], -1)
                        for w in raw])
    spec_db = torch.stack([db_tf(spec_tf(torch.from_numpy(np.ascontiguousarray(w.T)).float())) for w in cleaned])
    enc.noise_profiles = {("TU.DEMI", c): profiles[i] for i, c in enumerate("ZNE")}
    img = torch.stack([enc.normalize_spec(s, "TU.DEMI", ("Z", "N", "E")) for s in spec_db])
    enc.noise_profiles = {}
    img_fallback = torch.stack([enc.normalize_spec(s, "TU.DEMI", ("Z", "N", "E")) for s in spec_db])
    seq = np.stack([np.stack([standardize(w[:, c], mu=mu[c], sigma=sigma[c]) for c in range(3)], -1)
                    for w in cleaned]).astype(np.float32)

    # --- training-corpus tensors with their catalogue magnitudes -------------------------
    man = pd.read_csv(DATASET / "manifest.csv")
    rows = man[man.split == "test"].sort_values("magnitude").iloc[[0, len(man[man.split == "test"]) // 2, -1]]
    corpus = [torch.load(DATASET / r.split / r.filename, weights_only=True) for r in rows.itertuples()]
    c_seq = torch.stack([d["seq"].float() for d in corpus])
    if transform == "asinh":
        # Both the fixture's `seq` (compared against the C++ estimator's own
        # input) and the tensors fed to the model below are post-transform, as
        # the trainer's dataset applies it on load.
        seq = np.arcsinh(seq).astype(np.float32)
        c_seq = torch.asinh(c_seq)
    c_img = torch.stack([d["img"].float() for d in corpus])

    out = {"noise_raw": noise_raw, "noise_mu": mu, "noise_sigma": sigma, "noise_profile": profiles.astype(np.float32),
           "raw": raw, "cleaned": cleaned, "spec_db": spec_db.numpy(), "img": img.numpy(),
           "img_fallback": img_fallback.numpy(), "seq": seq,
           "corpus_seq": c_seq.numpy(), "corpus_img": c_img.numpy(),
           "corpus_label": rows.magnitude.to_numpy(np.float32)}
    x_seq = torch.cat([torch.from_numpy(seq), c_seq])
    x_img = torch.cat([img, c_img])
    for p, m in models.items():
        out[f"mag_p{p}"] = m(x_seq, x_img, None).numpy()
    m = models[0]
    lstm, _ = m.b1.lstm(x_seq[:1])
    out["p0_lstm"] = lstm.numpy()
    out["p0_b1"] = m.b1(x_seq).numpy()
    out["p0_b2"] = m.b2(x_img).numpy()
    write_ayzw(Path("data/fixtures/magnitude.ayzw"), out,
               {"windows": starts, "noise_windows": noise_starts, "corpus_files": rows.filename.tolist()})

    ens = np.mean([out[f"mag_p{p}"] for p in models], axis=0)
    print("DEMI M4.9 window (catalogue Mw 4.9): %.2f   quiet window: %.2f" % (ens[0], ens[1]))
    for label, est in zip(out["corpus_label"], ens[2:]):
        print(f"corpus test window, label {label:.1f}: {est:.2f}")
    print("DEMI noise sigma (Z, N, E):", np.round(sigma, 1))


if __name__ == "__main__":
    main()
