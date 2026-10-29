# 10 · Benchmarks and the scorecard

Two tools, one for speed and one for quality, so that a change to the code or
the models is judged on both.

| | what | command |
|---|---|---|
| `ayzek_bench` | time per stage and per layer, real-time capacity | `build-release/bench/ayzek_bench` |
| `tools/scorecard.py` | detection, false alarms, magnitude and location on fixed datasets | `python3 tools/scorecard.py` |
| `tools/bench_compare.py` | two `ayzek_bench` results side by side | `python3 tools/bench_compare.py A.json B.json` |

## Speed: `ayzek_bench`

```bash
meson test -C build-release --benchmark        # all groups, JSON in build-release/bench/
build-release/bench/ayzek_bench --group stages,layers,capacity --json run.json
build-release/bench/ayzek_bench --layer lstm:T=125,in=96,h=64 --layer attn:T=125,d=128
```

| group | what it times |
|---|---|
| `stages` | the pipeline's work with the shipped models: conditioning, one detector and the ensemble, 50 STA/LTA samples, the per-window total (every 0.5 s per station), the picker, the Pd chain and estimate, the transformer (until v0.1.0-beta.6 also the retired magnitude regressor and its noise baseline) |
| `layers` | every network layer at the shape a shipped model runs it at, with multiply-accumulates per call and the rate achieved |
| `shapes` | layers at shapes given with `--layer KIND:K=V,...` (kinds `conv1d`, `conv2d`, `lstm`, `attn`, `linear`), to price an architecture before it exists in C++ |
| `capacity` | detector windows per second on 1 to N threads, the stations that fit (two windows per station per second), and the magnitude estimate's latency while the other cores run detectors |

Inputs are synthetic and seeded. Inference time does not depend on the sample
values, so no waveform data is needed and the benchmark runs on a bare board.
Each benchmark warms up, then times calls for `--min-time` seconds (1 s by
default, at least 10 calls) and reports p50, p90 and p99. `--quick` does a few
calls of everything; it runs as the ordinary test `bench_quick`, so the
benchmark keeps building on every build type.

**It refuses to run on an unoptimised build** unless given `--allow-debug`.
Every result records the CPU, core count, compiler, SIMD backend, whether the
build was optimised, and the commit (`git describe --dirty`). The guard exists
because `build-release` was for a while a `-O0` build left over from before the
repository moved, and ran the models about 17 times slower without anything
saying so.

### Results on the development machine

AMD Ryzen 5 5600 (6 cores, 12 threads), GCC 16.2, release build, scalar
backend (the NEON kernels are aarch64 only):

| stage | p50 |
|---|---:|
| condition a 6 s window | 0.04 ms |
| detector, one model | 1.7 ms |
| detector ensemble (3 models) | 5.2 ms |
| STA/LTA, 50 samples (one 0.5 s step) | 0.0005 ms |
| per-window total (conditioning, ensemble, STA/LTA) | 5.2 ms |
| picker, 60 s window | 12.0 ms |
| magnitude estimate (3 models) | 294 ms |
| magnitude, one model | 96 ms |
| transformer, one 0.1 s token | 0.094 ms |
| transformer stream, 0.5 s of one station (5 tokens, filtering) | 0.49 ms |
| transformer context refresh (once a quiet minute) | 8.6 ms |
| Pd displacement chain, 50 samples | 0.0007 ms |
| Pd magnitude estimate, 10 stations, every window | 0.26 ms |

Measured 21 October 2026 at 7748168 on an idle machine (load average below
1.6 at the start of each group).

- **The detector costs about 1% of a core per station** (5.2 ms every 0.5 s).
  One thread sustains 191 windows a second, 95 stations; four threads 379
  stations; all twelve 653, where SMT and the shared caches take their share.
- **The transformer is about ten times cheaper**: 0.49 ms per 0.5 s of a
  station, because it steps one 0.1 s token against cached keys and values
  instead of rescoring a 6 s window. One thread sustains 882 stations, all
  twelve 5,956. Pd magnitude costs nothing measurable beside it.
- **STA/LTA costs about 10,000 times less than the ensemble** (0.0005 ms per
  0.5 s step).
- **The magnitude estimate is the expensive stage**, and 91% of one model is a
  single layer: self-attention over the 1000 raw samples of the 1D branch
  (85.6 of 94 ms). With the other cores busy running detectors it takes
  480 ms (p99 602 ms), and it runs on the station's own thread.
- **A lighter 1D branch, priced with `--layer`:** the detector's design
  (three strided convolutions to 125 steps, BiLSTM, attention over 125 steps)
  costs 4.0 ms against 92 ms, which would make a magnitude model about 7 ms.

### Raspberry Pi 4

Raspberry Pi 4 Model B (Cortex-A72, 4 cores at 1.8 GHz), NEON backend, the
static `build-pi` binary, 21 October 2026. Each group started below 55 °C with
the clock and throttle flags logged every second; the clock held 1.8 GHz
throughout. A Pi 4 without a heatsink reaches 80–85 °C under the capacity
group and its accuracy runs, and throttles there; timings taken hot read slow.

| stage | p50 | x86 |
|---|---:|---:|
| condition a 6 s window | 0.14 ms | 0.04 ms |
| detector, one model | 11.8 ms | 1.7 ms |
| detector ensemble (3 models) | 36.3 ms | 5.2 ms |
| STA/LTA, 50 samples | 0.0019 ms | 0.0005 ms |
| picker, 60 s window | 75.7 ms | 12.0 ms |
| magnitude estimate (3 models) | 3.2 s | 294 ms |
| transformer stream, 0.5 s of one station | 3.9 ms | 0.49 ms |
| Pd magnitude estimate, 10 stations | 2.0 ms | 0.26 ms |

- **One core runs 14 stations** with the 6 s detector, four cores 52.
  The transformer runs 114 on one core but only 197 on four.
- **The results are the x86 results**: on demo_ko and quiet_ko the scorecards
  of the 6 s detector and of the transformer (with either magnitude) agree with
  x86 field for field.

`tools/bench_pi.sh HOST` runs the benchmark and the scorecard there. The JSON
of either machine can be compared with the other:

```bash
python3 tools/bench_compare.py bench/results/x86.json bench/results/pi.json --threshold 0
```

`bench_compare.py` flags every benchmark more than `--threshold` percent
(default 10) slower in the second file and exits 1 if any is. Timing on a busy
machine is not a reference: record reference results with nothing else
running.

## Quality: the scorecard

```bash
python3 tools/scorecard.py --json score.json                   # shipped models
python3 tools/scorecard.py --models NEW_MODELS --json new.json # candidate models
python3 tools/scorecard.py --args '--detector stalta' --json stalta.json
python3 tools/scorecard.py --only demo,marmara                 # a subset
python3 tools/scorecard.py --compare score.json new.json
```

| dataset | what | role |
|---|---|---|
| `demo` | Sındırgı 2025-11-10, 30 min, DEMI MANT BAND, the Mw 4.9 | tuning |
| `marmara` | Marmara Sea 2025-04-23, 3 h, 8 stations, the Mw 6.2 | tuning |
| `sindirgi61` | Sındırgı 2025-08-10, 3 h, 6 stations, the Mw 6.1 | held out |
| `noise_a1`, `a2`, `b`, `c` | four windows with no AFAD event within 250 km of any station, 33.3 h in total (`data/runs/noise`) | false alarms |

"Tuning" and "held out" follow `05-pipeline.md`: the trigger settings were
chosen on the first two, and the third was not looked at until they were fixed.
The scorecard keeps them apart so that a change tuned on it can still be
checked on data it has not seen.

Per dataset it reports catalogue events detected, alarms, alarms without a
catalogue event, the median alarm delay after origin and the mean magnitude
error of the matched alarms; on the quiet windows, false alarms per day. For
every matched event of AFAD M ≥ 4 it lists the alarm delay, the magnitude at
the alarm and final, and the epicentre error. The summary pools these into the
headline numbers: detections and unmatched alarms (tuning and held out
separately), false alarms per day, magnitude error over all matched events and
over M ≥ 4, and the median epicentre error. `--bench` adds the p50 of the
per-window, picker and magnitude stages from `ayzek_bench`.

- **Catalogues are fixed.** Each dataset is scored against its excerpt in
  `tests/catalogs/`, the AFAD events within a day of the data, cut from the
  2026-08-30 export. A score therefore does not move when AFAD revises its
  catalogue. On the demo the excerpt gives the same report as the full export.
- **The data is not tracked.** `data/` is built with `tools/make_demo_data.py`
  (commands in the README and `05-pipeline.md`); a dataset whose files are
  missing is skipped and listed as such.
- **Magnitude error is rounded per event.** ayzek prints each event's error to
  0.1, and the scorecard averages those, so its mean can differ from ayzek's
  own summary line in the second decimal.

The scorecard for the shipped models is in `bench/results/scorecard_baseline.json`.
