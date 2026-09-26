# The environment every number here was produced in

Recorded 2026-09-25. A reproduction does not need these exact versions, but a
disagreement in the last digit of a benchmark usually traces to one of them, so
they are written down rather than remembered.

Regenerate this table with `tools/env_report.sh`.

## Build toolchain

| component | version | where it is used |
|---|---|---|
| **zig** | **0.16.0** | cross-compiles all three non-native targets |
| ↳ clang inside zig | **21.1.8** | what actually compiles the cross builds |
| ↳ musl | **1.2.5** | static Linux and Raspberry Pi binaries |
| ↳ mingw-w64 | **13.0** | Windows target |
| **gcc / g++** | **16.2.1** (20260810) | native `build-release`, `build-tsan`, `build-asan` |
| clang (system) | 22.1.8 | **not** used for builds |
| **glibc** | **2.44** | native builds link against this |
| meson | 1.12.1 | |
| ninja | 1.13.2 | |
| qemu-aarch64 | **11.1.1** | `exe_wrapper`, runs the aarch64 tests on x86 |

Two distinctions are easy to get wrong and both matter in print:

- the cross builds report **clang 21.1.8**, which is the LLVM bundled inside
  zig, not the system clang 22.1.8;
- the static binaries are **musl 1.2.5**; only the native build uses glibc 2.44.

Exact distribution strings (Arch, rolling): `zig 0.16.0-1`,
`gcc 16.2.1+r23+gd564253eb6c8-1`, `clang 22.1.8-1`,
`glibc 2.44+r24+g16be1518495f-1`, `meson 1.12.1-1`, `ninja 1.13.2-3`,
`qemu-user 11.1.1-4`.

## Targets

| build dir | triple | libc | CPU | tests run by |
|---|---|---|---|---|
| `build-release` | x86_64-linux-gnu | glibc 2.44 | native | directly |
| `build-pi` | aarch64-linux-musl | musl 1.2.5 | `cortex_a72` (Pi 4; `cortex_a76` for Pi 5) | qemu-aarch64 |
| `build-linux` | x86_64-linux-musl | musl 1.2.5 | generic, static | directly |
| `build-win` | x86_64-windows-gnu | mingw-w64 13.0 | generic | wine, or a Windows host |

## Python, and the two environments that are not the same

The detector scans and the false-alarm runs use `cnn_earthquake`; the magnitude
corpus and training use `cascade_impl`. They have drifted, which is why the
column is split rather than averaged.

| | `cnn_earthquake` | `cascade_impl` |
|---|---|---|
| python | 3.12.13 | 3.12.13 |
| **torch** | **2.13.0** | **2.14.0** |
| numpy | 2.5.1 | 2.5.3 |
| scipy | 1.18.0 | 1.18.1 |
| obspy | 1.5.0 | 1.5.1 |
| matplotlib | 3.11.1 | 3.11.2 |
| pandas | 3.0.3 | 3.0.5 |
| scikit-learn | 1.9.0 | 1.9.0 |
| h5py | 3.16.0 | — |

The system interpreter is 3.14.7 and `uv` is 0.12.18; neither runs any of the
science. The C++ pipeline has no Python dependency at all — the tools around it
do.

## Machine

| | |
|---|---|
| CPU | AMD Ryzen 5 5600, 6 cores / 12 threads |
| RAM | 31 GB |
| GPU | **NVIDIA RTX 3060 Ti, 8 GB**, driver 615.71.09 |
| CUDA / cuDNN (through torch) | 13.0 / 9.24.00 |
| OS | Arch Linux, kernel 7.2.6-arch2-1, x86_64 |

**The 8 GB card is a real constraint on reproduction, not a detail.** Three
concurrent scan processes fit; four exhaust it (`ayzek_paper/revisions2.md`
R2.2), and a scan cannot share the card with a training run. Both failures
observed during this work were that, not a bug.

Timings in `docs/impl/10-benchmarks.md` are single-threaded per station on this
CPU unless the text says otherwise; the replay throughput figures assume the
stations run in parallel across these 12 threads.
