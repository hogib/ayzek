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
- **Located** 5.3 km from AFAD's epicentre, origin within 1.0 s, rms 0.16 s.
- **Aftershocks:** 6 of the other 8 catalogue events declared. Four (ML 1.4,
  1.9, 3.6, 2.7) are located 3.7–6.1 km from AFAD; two located from two
  distant stations are 100 and 199 km off. Magnitude error averages 0.40
  units over the seven declared events.
- **Unmatched alarms:** two of nine alarms in 30 minutes have no catalogue
  event.
- **Trigger rule:** 8 consecutive windows at p ≥ 0.8, or one at p ≥ 0.9. The
  earlier rule (one window at p ≥ 0.9) declared 5 of 9 events here, and 24
  instead of 37 of 66 on a held-out 3 h replay of the Sındırgı Mw 6.1
  (`impl/05-pipeline.md`, *Threshold*).
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
  --threshold P         trigger when N consecutive windows reach P (default 0.8) ...
  --trigger-windows N   ... (default 8, which adds 3.5 s)
  --instant-threshold P or when one window reaches P (default 0.9; > 1 disables)
  --release P           probability below which a trigger resets (default 0.3)
  --step N              samples between detector windows (default 50 = 0.5 s)
  --min-stations N      detections needed to declare an event (default 2)
  --no-pick             detector only
  --no-magnitude        skip the magnitude regressor
  --catalog CSV         AFAD catalogue export to score events against
  --scores DIR          write every window's probability to DIR/STATION.csv
  --scores-in DIR       use probabilities from an earlier --scores run instead of the detector
  --record FILE         write all station outputs, for tools/network_subsets
  --site NAME,LAT,LON   also report the S-wave warning time at this place (repeatable)
  --verbose             print every station detection, pick and magnitude estimate
  --no-color
```

## Output

While running, one line per alarm and one per magnitude or location update:

```
ALARM    #3   18:21:09.00  earthquake detected by DEMI, MANT  (AFAD MW 4.9 at 18:20:51.00, 18.0 s ago)
MAG      #3   18:21:09.00  magnitude M4.3 from 1 station (DEMI 4.3)
LOCATE   #3   18:22:05.00  located at 39.221N 28.088E, origin 18:20:52.02, rms 0.16 s from 3 stations
```

At the end, a report with one block per alarm, a summary, and station
processing statistics:

```
Event #3: detected M4.6 earthquake at 39.221N 28.088E, origin 18:20:52.02 UTC
  alarm      18:21:09.00 UTC by DEMI, MANT; later BAND
  magnitude  M4.3 at the alarm, M4.6 final from 3 stations
  location   rms 0.16 s from 3 stations
  AFAD       MW 4.9 at 18:20:51.00 UTC: alarm 18.0 s after origin, magnitude -0.3, epicentre 5.3 km off, origin +1.0 s
  S wave     had travelled 62 km from the epicentre when the alarm sounded (AFAD hypocentre)
  warning               distance   S wave arrives            warning
             DEMI          54 km   18:21:08.76 UTC (picked)      -0.2 s  after S
             MANT          91 km   18:21:17.88 UTC (picked)      +8.9 s  before S
             BAND         124 km   18:21:28.16 UTC (picked)     +19.2 s  before S
             Istanbul     209 km   18:21:50.80 UTC (predicted)   +41.8 s  before S
```

The warning time is the S-wave arrival minus the alarm time: positive means
the alarm came before the S wave. At stations with an S pick for the event the
picked time is used, otherwise the arrival is predicted from the AFAD
hypocentre (or ayzek's location without a catalogue) with Vs = 3.5 km/s.

The summary counts alarms that match an AFAD event (correct) and alarms without
one (false; some may be real earthquakes too small for the catalogue), AFAD
events detected and missed, detection by magnitude band, and lists each missed
event and each false alarm.

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
