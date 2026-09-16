# 06 · Raspberry Pi

The target is a 64-bit Raspberry Pi 4 (Cortex-A72) or 5 (Cortex-A76). The build
is a single statically linked aarch64 binary, so deploying it means copying
one file and the `models/` directory. Nothing needs installing on the Pi.

## Toolchain

No aarch64 cross-compiler is packaged here, so the build uses **zig as a C++
compiler**. It bundles clang and libc++ with musl for every target, and uv
fetches and caches it:

```
cross/zig-c++   uv run --with ziglang python -m ziglang c++ -target aarch64-linux-musl
cross/zig-cc    the same for C
cross/zig-ar    the archiver
cross/aarch64-pi.ini   the Meson cross file
```

The wrappers drop two flags Meson adds that zig's linker driver rejects or
warns about: `-Wl,--fatal-warnings` and `-Wl,-O1`. Compile warnings still fail
the build under `werror`.

```bash
meson setup build-pi --cross-file cross/aarch64-pi.ini --buildtype=release
ninja -C build-pi
meson test -C build-pi --timeout-multiplier 10     # runs under qemu-aarch64
```

For a Pi 5, change `cpu = 'cortex_a72'` to `'cortex_a76'` in the cross file.
zig spells CPU names with underscores.

## Portability fixes the Pi build forced

- **`std::chrono::parse`** is not implemented in libc++. Both time formats are
  parsed by hand in `parse_utc` (`src/pipeline/common.hpp`).
- **libc++ headers** trip clang's `-Wnullability-completeness` under `-Werror`.
  It is disabled for the cross build only.
- **`assert` in the stage-1 tests** compiled away in release builds, so those
  tests checked nothing there. clang's unused-variable warning exposed this.
  They use an aborting `CHECK` now, on every build.

## Verification

Every test runs on the aarch64 binary under `qemu-aarch64` with the NEON
backend, including the comparison against PyTorch and the end-to-end M4.9 run.
`test_models` prints `neon backend`, and its detector logit matches PyTorch to
2.4e-7.

Timing under qemu is emulation, not the Pi. Real Pi numbers come from running
`tools/demo.sh` on the device. The summary prints per-window and per-pick
milliseconds, and the throughput in station-seconds per second.

## Deploying

```bash
scp build-pi/app/ayzek pi@raspberrypi:~/ayzek/
scp -r models data/demo tools/demo.sh pi@raspberrypi:~/ayzek/
ssh pi@raspberrypi 'cd ayzek && ./ayzek --speed 10 --from 2025-11-10T18:20:00 data/demo/DEMI.mseed data/demo/MANT.mseed data/demo/BAND.mseed'
```

Budget: three models per window every 0.5 s is 6 inferences a second per
station. At x86's 5.2 ms per window, a station costs about 1% of one core. A Pi
core is several times slower, which still leaves room for tens of stations.
