# ayzek — implementation

A real-time earthquake early-warning pipeline in C++23 with hand-written
inference. It replays AFAD miniSEED through a trained detector, a P/S picker,
a magnitude regressor, network association and location, and it runs on a
Raspberry Pi.

`DESIGN.md` is the reasoning behind the architecture. The documents here
describe what was built.

| | |
|---|---|
| [01 · SIMD kernels](impl/01-simd.md) | NEON on aarch64, scalar elsewhere; where vectorising pays and where it does not |
| [02 · Weights and fixtures](impl/02-weights.md) | exporting PyTorch checkpoints, the AYZW format, reference outputs |
| [03 · Window conditioning](impl/03-dsp.md) | detrend, taper, scipy-exact `filtfilt`, per-window standardisation |
| [04 · Models](impl/04-models.md) | the layers, the 3-seed detector, the picker, and their agreement with PyTorch |
| [05 · Pipeline](impl/05-pipeline.md) | threads, rings, trigger, association, location, the M4.9 demo and its results |
| [06 · Raspberry Pi](impl/06-raspberry-pi.md) | cross-compiling with zig, running under qemu, deploying |
| [07 · Magnitude](impl/07-magnitude.md) | the regressor, station noise baselines, spectrogram, results |
| [08 · Station selection](impl/08-station-selection.md) | recording and subset replay, geometry model, 2025-04-23 Mw 6.2 Marmara results |

## In one screen

```
AFAD miniSEED ─> decode (Steim2) ─> reorder ─> ring Z/N/E ─┐   one pair of threads per station
                                                           ├─> 6 s windows every 0.5 s
                                                           │   detrend · taper · 1–45 Hz · z-score
                                                           │   3-seed CNN-BiLSTM-attention (NEON)
                                                           │   trigger ─> 60 s P/S picker (NEON)
                                                           │   trigger ─> 10 s magnitude, station noise baseline (NEON)
                                                           └─> network: associate ─> declare ─> locate ─> median M
```

On the 2025-11-10 Sındırgı M4.9, from DEMI, MANT and BAND:

- **Declared** 18.0 s after origin, carrying **M4.3**. That is the second
  station's P travel time, not compute.
- **Magnitude** refined to **M4.6** once the picks are in, 63 s after origin.
- **Located** 4.9 km from AFAD's epicentre, origin within 1.1 s, rms 0.17 s.
- **Aftershocks:** four of them (ML 1.4, 1.9, 3.6, 2.7) located 4–6 km from
  AFAD. Magnitude error averages 0.40 units over all five, and 0.17 over the
  three inside the model's training range (M ≥ 2).
- **False alarms:** one declared event in 30 minutes is not in the catalogue,
  at threshold 0.9.
- **Speed:** 5.2 ms per detector window (three models) and 300 ms per
  magnitude estimate on x86. Seven stations with everything on run at 67× real
  time per station.
- **Agreement with training:** all three networks match PyTorch to about 1e-6
  on their outputs, the conditioning matches scipy to float rounding, the
  spectrogram matches torchaudio to 2e-3 dB, and the NEON build passes the same
  tests under qemu.

## First run

```bash
# 1. data: 30 minutes around the M4.9 from seven stations (stdlib Python)
python3 tools/make_demo_data.py --out data/demo \
    --start 2025-11-10T18:05:00 --end 2025-11-10T18:35:00 \
    ~/Projects/sismokaos/tdvms/afad_raw/{DEMI,MANT,BAND,KAND,KIRK}/*_2025-10-29.zip \
    ~/Projects/sismokaos/tdvms/afad_raw/{ELBA,SEMS}/*_2025-10-21.zip

# 2. weights and test fixtures from the PyTorch checkpoints
uv run --project ~/Projects/sismokaos/archive_pipeline python tools/export_models.py
uv run --project ~/Projects/sismokaos/data_downloader python tools/export_magnitude.py

# 3. build and test
meson setup build-release --buildtype=release
meson test -C build-release          # includes the end-to-end M4.9 check

# 4. demo
tools/demo.sh                        # 3 stations, 10x real time
tools/demo.sh all                    # 7 stations, full speed
```

## `ayzek` options

```
ayzek [options] STATION.mseed...
  --models DIR          weights directory (default: models)
  --speed X             replay speed; 1 = real time, 0 = as fast as possible
  --from TIME           start detecting at this UTC time
  --threshold P         detector trigger probability (default 0.9)
  --step N              samples between detector windows (default 50 = 0.5 s)
  --min-stations N      detections needed to declare an event (default 2)
  --no-pick             detector only
  --no-magnitude        skip the magnitude regressor
  --catalog CSV         AFAD catalogue export to score events against
  --scores DIR          write every window's probability to DIR/STATION.csv
  --record FILE         write all station outputs, for tools/network_subsets
```

## Layout

```
src/            ring, mseed, reorder (stages 1–2); simd, nn, models, dsp, weights, spectrogram, magnitude
src/pipeline/   ingest, processor, network
app/            the ayzek binary
tests/          unit tests, PyTorch/scipy agreement, end-to-end demo check
tools/          export_models.py, export_magnitude.py, make_demo_data.py, demo.sh, mseed_dump, replay_check,
                network_subsets, station_geometry.py, scan_encodings.py
cross/          zig toolchain wrappers and the Raspberry Pi cross file
```

## Status against DESIGN.md's stages

| stage | status |
|---|---|
| 1 ring buffer | done |
| 2 file replay | done (miniSEED, Steim2, reordering) |
| 3 filtering | done, but per window rather than continuous. Matching training wins (03). |
| 4 windowing | done, gap-aware |
| 5 kernels | done, NEON + scalar |
| 6 detector | done, matches PyTorch |
| 7 live SeedLink | not started; `ReplaySource` is the seam |
| 8 cascade | detector → picker → magnitude → location, done |
| 9 location | done: grid search, P+S, uniform half-space, drops ill-fitting stations |
