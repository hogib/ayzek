"""From filtered counts to the model's input. Numpy only, no scipy.

The model never sees counts. Each component is divided by that station's noise
RMS and compressed with asinh:

    x_in = asinh(x / sigma_station)

- Dividing by the station's own noise cancels instrument gain and leaves
  "how far above this station's background", which is the quantity a detector
  should trigger on.
- asinh is linear near zero and logarithmic in the tail, so an M6 at 10 km
  stays inside fp16/bf16 range while the noise keeps its resolution.

`sigma_station` comes from the station context when there is one, and from the
first `fallback_scale_s` of the stream when there is not (a station that has
just come up). Training uses both, in the ratio set by `DataConfig.ctx_drop`.

A missing sample is 0 on every component and 1 on the fourth, the gap channel.
"""
from __future__ import annotations

import numpy as np


def channel_scale(x: np.ndarray, missing: np.ndarray | None = None) -> np.ndarray:
    """Per-component RMS over the samples that are present; (T, 3) -> (3,)."""
    if missing is not None and missing.any():
        x = x[~missing]
    if len(x) == 0:
        return np.ones(x.shape[1], dtype=np.float32)
    rms = np.sqrt(np.mean(np.square(x, dtype=np.float64), axis=0))
    return np.where(rms > 0, rms, 1.0).astype(np.float32)


def condition(x: np.ndarray, missing: np.ndarray, scale: np.ndarray) -> np.ndarray:
    """(T, 3) filtered counts + (T,) missing mask -> (T, 4) model input."""
    out = np.empty((len(x), 4), dtype=np.float32)
    out[:, :3] = np.arcsinh(x / scale)
    out[missing, :3] = 0.0
    out[:, 3] = missing
    return out


def token_missing(missing: np.ndarray, stride: int) -> np.ndarray:
    """True for every token that has at least one missing sample in its span."""
    n = len(missing) // stride
    return missing[: n * stride].reshape(n, stride).any(axis=1)
