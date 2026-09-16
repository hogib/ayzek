# 03 · Window conditioning

A network with fixed weights expects inputs on exactly the scale and in exactly
the passband it was trained on. So `src/dsp.cpp` is a transcription of the
training pipeline (`archive_pipeline.archive.clean_block` and the `perwindow`
branch of `products/scan.py`), not a new design:

| step | C++ | Python it matches |
|---|---|---|
| remove linear trend | `detrend_linear` (closed-form least squares) | `signal.detrend(type="linear")` |
| remove mean | `detrend_constant` | `signal.detrend(type="constant")` |
| taper 5% at each end | Hann halves, `k = n // 20` | `taper_vector` |
| 1–45 Hz zero-phase bandpass | `FiltFilt`: DF-II transposed, forward then backward, odd extension of 27 samples, `lfilter_zi` initial state | `signal.filtfilt(b, a)` |
| standardise | per component, by the window's own mean and std | `(c - mean) / max(std, 1e-12)` |
| amplitude transform | `asinh`, inside `DetectorEnsemble` | `torch.asinh` in `score_block` |

Everything runs in double, as scipy does, and is cast to float only for the
network.

## Why per window, not continuous

`DESIGN.md` planned a continuous causal filter feeding a second ring. The
detector was trained on windows filtered *individually* with a zero-phase
filter, and a causal filter shifts phase and changes the onset shape the
network learned. Matching training beats the tidier architecture.

The cost is small. A window is 600 samples, and filtering it takes tens of
microseconds against about a millisecond for one model.

## Agreement with scipy

On six real DEMI windows (noise, several P onsets, S coda):

- cleaned signal: at most **1.1e-8 relative**. An 8th-order IIR in
  transfer-function form amplifies rounding, so this is expected, not a bug.
- standardised float window: at most **1.9e-9 absolute**, float32 rounding.

## The picker's conditioning

The sphase picker was trained on STEAD traces, z-scored per trace and per
component, with no filter of its own. ayzek cleans each 60 s window with the
same detrend, taper and bandpass as the detector, then z-scores it. That is a
choice, not a transcription: STEAD's waveforms are distributed pre-processed,
and filtering live data brings it closer to that. It is also a place a
domain-shift problem could hide, and the picker's transfer from STEAD to
Turkish broadbands is untested.
