# 01 · SIMD kernels

`src/simd.hpp` is the only file that knows which instruction set is in use.
Every layer in `nn.cpp` reduces to a handful of kernels, so vectorising those
vectorises the whole of inference.

| kernel | used by | NEON |
|---|---|---|
| `dot(float*)` | conv (per output position and channel), linear, LSTM gates, attention scores | `vmlaq_f32`, 4 accumulators, 16 elements a step |
| `gemv` | linear layers, LSTM input and recurrent products, attention projections | rows of `dot` |
| `axpy` | attention's weighted sum over values | `vmlaq_f32` |
| `add`, `affine`, `mul`, `relu` | residuals, folded batch norm, ReLU | `vaddq` / `vmulq` / `vmaxq` |
| `sum`, `dot(double*)`, `affine(double*)`, `mul_d` | window mean, variance, detrend, taper | `float64x2_t` |

## Backends

- **`neon`** is compiled whenever `__aarch64__` or `__ARM_NEON` is defined,
  which covers every 64-bit Raspberry Pi. It is written explicitly with
  intrinsics, so the Pi's speed does not depend on the compiler's
  auto-vectoriser.
- **`scalar`** is plain loops for x86 development builds. `dot` keeps four
  independent accumulators. A single `s += a[i] * b[i]` cannot be vectorised
  without `-ffast-math`, because that would reorder a float sum; four separate
  sums have no order to preserve, so GCC is free to vectorise them. That one
  change took a detector window from **15 ms to 5.2 ms** (three models).

`simd::kBackend` names the active backend at runtime. `test_models` prints it.

## Why not everything

- **GELU needs `erf`**, which has no NEON form. It runs once per activation,
  against one dot product per weight row, so it is a small share of the time.
- **LSTM gate nonlinearities** (`exp`, `tanh`) are scalar for the same reason.
- **The Butterworth recurrence is sequential by construction**: each output
  sample depends on the previous one. It is also about 1% of a window's work.
  Vectorising it across the three components is possible but not worth doing
  yet.

## Agreement

The two backends sum in different orders, so they agree to float rounding
rather than bit for bit. Both are held to PyTorch by `test_models`: the NEON
build runs under `qemu-aarch64` and matches the reference logits to 2.4e-7,
the scalar build to 6e-7.
