# ayzek — implementation

A real-time earthquake early-warning pipeline in C++23 with hand-written
inference. It replays AFAD miniSEED through a trained detector, a P/S picker,
network association and location, and a magnitude from peak P displacement,
and it runs on a Raspberry Pi.

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
| [07 · Magnitude (retired)](impl/07-magnitude.md) | the learned regressor of v0.1.0-beta.6 and earlier, station noise baselines, spectrogram, results |
| [08 · Station selection](impl/08-station-selection.md) | recording and subset replay, geometry model, 2025-04-23 Mw 6.2 Marmara results |
| [09 · STA/LTA benchmark](impl/09-sta-lta.md) | a recursive STA/LTA in place of the detector, tuned and compared on the same replays |
| [10 · Benchmarks and the scorecard](impl/10-benchmarks.md) | time per stage and layer, real-time capacity, and one scorecard of detection, false alarms and magnitude for every model change |
| [11 · Windows build](impl/11-windows.md) | cross-compiling a self-contained ayzek.exe with zig |
| [12 · Environment](impl/12-environment.md) | the toolchain and machines every number here was produced on |
| [13 · Claims](impl/13-claims.md) | every published number, and the command that reproduces it |
| [14 · Streaming transformer](impl/14-transformer.md) | the onset transformer as a detector: causal band-pass, station context, token trigger with dt restarts, the picker placed on its P |
| [15 · Location from the geometry head](impl/15-geometry-location.md) | per-station distance in place of S-P picks, the locator, its settings and geo_v1's numbers |
| [16 · Alarm assessment](impl/16-alarm-assessment.md) | `--assess`: S-P consistency and independent evidence to tell uncatalogued earthquakes from misfires |
| [17 · Magnitude from peak P displacement](impl/17-pd-magnitude.md) | the Pd chain on raw counts, the censored Pd–M–R estimator calibrated in onset, its agreement with onset and its results on the replays |
| [18 · FilterPicker](impl/18-filterpicker.md) | Lomax's FilterPicker5 ported from his C code, checked against it, in place of the detector |

## In one screen

```
AFAD miniSEED ─> decode (Steim2) ─> reorder ─> ring Z/N/E ─┐   one pair of threads per station
                                                           ├─> 6 s windows every 0.5 s
                                                           │   detrend · taper · 1–45 Hz · z-score
                                                           │   3-seed CNN-BiLSTM-attention (NEON)
                                                           │   trigger ─> 60 s P/S picker (NEON)
                                                           │   trigger ─> peak P displacement at the onset (Pd)
                                                           └─> network: associate ─> declare ─> locate ─> Pd M
```

On the 2025-11-10 Sındırgı M4.9, from DEMI, MANT and BAND:

- **Declared** 18.0 s after origin. That is the second station's P travel
  time, not compute.
- **Magnitude:** none here. The Pd magnitude needs the stations' gains, which
  `models/pd_gain.csv` does not yet hold for the AFAD stations. On the six KO
  stations of `demo_ko` the M4.9 gets M4.7 3.5 s after the alarm and M4.8
  final (`impl/17-pd-magnitude.md`).
- **Located** 5.3 km from AFAD's epicentre, origin within 1.0 s, rms 0.16 s.
- **Aftershocks:** 6 of the other 8 catalogue events declared. Four (ML 1.4,
  1.9, 3.6, 2.7) are located 3.7–6.1 km from AFAD; two located from two
  distant stations are 100 and 199 km off.
- **Unmatched alarms:** two of nine alarms in 30 minutes have no catalogue
  event.
- **Trigger rule:** 8 consecutive windows at p ≥ 0.8, or one at p ≥ 0.9. The
  earlier rule (one window at p ≥ 0.9) declared 5 of 9 events here, and 24
  instead of 37 of 66 on a held-out 3 h replay of the Sındırgı Mw 6.1
  (`impl/05-pipeline.md`, *Threshold*).
- **Speed:** 5.2 ms per detector window (three models) on x86.
- **Agreement with training:** all three networks match PyTorch to about 1e-6
  on their outputs, the conditioning matches scipy to float rounding, the Pd
  chain matches onset's reference to 1e-9, and the NEON build passes the same
  tests under qemu.

## First run

```bash
# 1. data: 30 minutes around the M4.9 from seven stations (stdlib Python)
python3 tools/make_demo_data.py --out data/demo \
    --start 2025-11-10T18:05:00 --end 2025-11-10T18:35:00 \
    ~/Projects/sismokaos/tdvms/afad_raw/{DEMI,MANT,BAND,KAND,KIRK}/*_2025-10-29.zip \
    ~/Projects/sismokaos/tdvms/afad_raw/{ELBA,SEMS}/*_2025-10-21.zip

# 2. weights and test fixtures from the PyTorch checkpoints
uv run --project tools python tools/export_models.py \
    --detector-dir CKPT_DIR --picker wave_n250.pt

# 3. build and test
meson setup build-release --buildtype=release
meson test -C build-release          # includes the end-to-end M4.9 checks

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
  --detector KIND       model (default) or stalta, the reference STA/LTA trigger
  --sta S, --lta S      STA/LTA averaging lengths in seconds (default 1, 30)
  --stalta-on R         STA/LTA ratio that triggers (default 8)
  --stalta-off R        STA/LTA ratio below which a trigger resets (default 1.5)
  --stalta-band LO,HI   STA/LTA pass band in Hz (default 2,20)
  --stalta-3c           STA/LTA on the energy of all three components (default: vertical)
  --no-anchor           do not date model triggers by the STA/LTA onset
  --anchor-on R         STA/LTA ratio taken as the onset when anchoring (default 3)
  --anchor-picker       place the picker window from the anchored P too
  --require-onset       a detector trigger also needs an STA/LTA onset
  --step N              samples between detector windows (default 50 = 0.5 s)
  --min-stations N      detections needed to declare an event (default 2;
                        transformer: 3)
  --slack S             tolerance on the inter-station P travel time, seconds (default 3)
  --no-pick             detector only
  --no-magnitude        no magnitude estimate (the Pd magnitude)
  --catalog CSV         AFAD catalogue export to score events against
  --scores DIR          write every window's probability to DIR/STATION.csv
  --scores-in DIR       use probabilities from an earlier --scores run instead of the detector
  --record FILE         write all station outputs, for tools/network_subsets
  --site NAME,LAT,LON   also report the S-wave warning time at this place (repeatable)
  --verbose             print every station detection, pick and Pd value
  --no-color
```

## Output

While running, one line per alarm and one per magnitude or location update:

```
ALARM    #3   18:21:09.00  earthquake detected by KULA, SIMA  (AFAD MW 4.9 at 18:20:51.00, 18.0 s ago)
PLOC     #3   18:21:10.50  P times locate it at 39.240N 28.119E +-23 km, origin 18:20:52.18, rms 0.01 s from 3 stations
PDMAG    #3   18:21:12.50  Pd magnitude M4.7 ± 0.2 (M<4 0%, 4-5.5 100%, M>=5.5 0%) from 3 stations (3 above noise), windows up to 5 s, distances from the P-time location
LOCATE   #3   18:22:04.50  located at 39.206N 28.137E, origin 18:20:52.66, rms 0.21 s from 5 stations
```

(demo_ko; `README.md` lists every option.)

At the end, a report with one block per alarm, a summary, and station
processing statistics:

```
Event #3: detected M4.8 earthquake at 39.206N 28.137E, origin 18:20:52.66 UTC
  alarm      18:21:09.00 UTC by KULA, SIMA; later DKL, KCTX, TVSB, YAYO
  Pd magn.   M4.7 3.5 s after the alarm, M4.8 ± 0.1 (M<4 0%, 4-5.5 100%, M>=5.5 0%) final from 6 stations (6 above noise)
  location   rms 0.21 s from 5 stations, first 51.5 s after the alarm (P and S picks)
  P loc.     39.247N 28.121E +-8 km from 6 stations, first 1.5 s after the alarm, epicentre 2.2 km off
  AFAD       MW 4.9 at 18:20:51.00 UTC: alarm 18.0 s after origin, Pd magnitude -0.1, epicentre 3.7 km off, origin +1.7 s
  S wave     had travelled 62 km from the epicentre when the alarm sounded (AFAD hypocentre)
  warning               distance   S wave arrives            warning
             SIMA          74 km   18:21:14.42 UTC (picked)      +5.4 s  before S
             KULA          92 km   18:21:18.14 UTC (picked)      +9.1 s  before S
             YAYO         101 km   18:21:22.02 UTC (picked)     +13.0 s  before S
             DKL          109 km   18:21:23.54 UTC (picked)     +14.5 s  before S
             KCTX         115 km   18:21:23.97 UTC (predicted)   +15.0 s  before S
             TVSB         116 km   18:21:25.98 UTC (picked)     +17.0 s  before S
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
src/            ring, mseed, reorder (stages 1–2); simd, nn, models, dsp, weights, stalta, filterpicker, transformer, pd_magnitude
src/pipeline/   ingest, processor, network
app/            the ayzek binary
tests/          unit tests, PyTorch/scipy agreement, end-to-end demo check
tools/          export_models.py, export_pd.py, export_stalta.py, make_demo_data.py, demo.sh, mseed_dump,
                replay_check, network_subsets, station_geometry.py, scan_encodings.py, sweep_trigger.py
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
| 8 cascade | detector → picker → magnitude → location, done; the magnitude at the picked P runs only for declared events |
| 9 location | done: grid search, P+S, uniform half-space, drops ill-fitting stations |
