"""Causal band-pass filtering. Needs scipy (the `build` extra).

The filter is causal, and that is the point. A zero-phase filter (`filtfilt`)
runs backwards as well as forwards and smears onset energy *ahead* of the
arrival. A model trained on that learns to "see" P before it has arrived, and
deployment, which can only filter forwards, takes that ability away. Training,
replay and the C++ runtime all filter with this: 4th-order Butterworth, 1–45 Hz,
second-order sections, state carried forward.

Each contiguous run of samples is filtered on its own, starting from
steady-state initial conditions, and missing samples are 0. A gap resets the
filter, just as a real-time filter has to restart after a gap.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy import signal

from .config import SAMPLE_RATE

BAND_HZ = (1.0, 45.0)
ORDER = 4
MIN_RUN_S = 1.0          # shorter runs are marked missing: too short to filter


@lru_cache(maxsize=None)
def bandpass_sos(fs: float = SAMPLE_RATE, band=BAND_HZ, order: int = ORDER):
    return signal.butter(order, band, btype="bandpass", fs=fs, output="sos")


def runs(valid: np.ndarray) -> list[tuple[int, int]]:
    """Half-open [start, end) spans where `valid` is True."""
    v = np.concatenate([[False], valid, [False]]).astype(np.int8)
    d = np.diff(v)
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def causal_filter(x: np.ndarray, missing: np.ndarray, fs: float = SAMPLE_RATE):
    """Filters one component.

    Args:
        x: (T,) samples; values under `missing` are ignored.
        missing: (T,) bool.

    Returns:
        (filtered (T,) float32, missing (T,) bool). The returned mask also
        covers runs too short to filter.
    """
    sos = bandpass_sos(fs)
    out = np.zeros(len(x), np.float32)
    missing = missing.copy()
    for a, b in runs(~missing):
        if b - a < MIN_RUN_S * fs:
            missing[a:b] = True
            continue
        seg = x[a:b].astype(np.float64)
        seg = seg - seg[: int(fs)].mean()        # DC from the run's first second: causal
        zi = signal.sosfilt_zi(sos) * seg[0]
        out[a:b] = signal.sosfilt(sos, seg, zi=zi)[0]
    return out, missing


def filter_components(x: np.ndarray, missing: np.ndarray, fs: float = SAMPLE_RATE):
    """(T, 3) samples and (T, 3) per-component missing -> filtered (T, 3) and a
    combined (T,) mask that is True wherever any component is missing."""
    y = np.zeros(x.shape, np.float32)
    m = np.zeros(x.shape, bool)
    for c in range(x.shape[1]):
        y[:, c], m[:, c] = causal_filter(x[:, c], missing[:, c], fs)
    combined = m.any(axis=1)
    y[combined] = 0.0
    return y, combined
