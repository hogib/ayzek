# 11 · Windows build

A 64-bit Windows build, cross-compiled from Linux with zig, as the Raspberry Pi
build is (`06-raspberry-pi.md`). It produces `ayzek.exe` and every tool and
test, each self-contained: the only DLLs they import are `KERNEL32` and the
Universal CRT, which ships with Windows 10 and 11.

```bash
meson setup build-win --cross-file cross/x86_64-windows.ini --buildtype=release
ninja -C build-win
```

| | |
|---|---|
| `cross/x86_64-windows.ini` | the cross file |
| `cross/zig-c++-windows`, `zig-cc-windows`, `zig-ar-windows` | zig with `-target x86_64-windows-gnu`, fetched by uv, dropping the two linker flags Meson adds that zig rejects |

`build-win/app/ayzek.exe` is 1.6 MB and runs from any directory that has
`models/` beside the data:

```
ayzek.exe --speed 0 --catalog tests\catalogs\demo.csv data\demo\DEMI.mseed ...
```

## What the port needed

Two calls, both in terminal and machine reporting:

- `isatty(STDOUT_FILENO)`, which decides whether to colour the output, is
  `_isatty(_fileno(stdout))` on Windows. `app/ayzek.cpp` wraps it in
  `stdout_is_tty()`. The first attempt defined `STDOUT_FILENO` as a macro,
  which MinGW's headers already define, and `-Werror` caught it.
- `gethostname` in `ayzek_bench` needs winsock on Windows, which is not worth
  linking for a label; it reads `%COMPUTERNAME%` instead.

Nothing else was platform-specific. In particular:

- **Every binary read already opens with `std::ios::binary`** — the weights
  loader, the miniSEED ingest, `mseed_dump` and `replay_check`. A stream left
  in text mode would have corrupted Steim frames on Windows only, which is the
  kind of fault that reaches a user before a test.
- **Paths are built with `/`**, which Windows accepts.
- `ayzek_bench` reads the CPU name from `/proc`, absent on Windows; it reports
  an empty name rather than failing.
- The pipeline's threading is `std::thread` and `std::atomic`, and the SIMD
  kernels are aarch64-only, so Windows takes the same scalar path as x86 Linux.

## Not yet verified by running

The binaries are **built and linked, not executed**: no Windows machine and no
wine on this one. The cross file names `wine` as its `exe_wrapper`, so once
wine is installed the suite runs the same way the Pi build runs under qemu:

```bash
meson test -C build-win            # needs wine
```

Until then, the Windows build proves only that the code compiles and links for
the target. The numerical agreement tests (`test_models`, `test_magnitude`,
`test_dsp`, `test_stalta`) are the ones worth running there first: they compare
against the PyTorch and scipy fixtures and would catch a floating-point or
library difference.
