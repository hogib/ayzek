# ayzek

An earthquake early-warning (EEW) pipeline in C++23. It reads continuous
three-component waveforms from several stations and, for each earthquake:

1. detects it at each station (a CNN–BiLSTM–attention detector, 3-seed
   ensemble; a recursive STA/LTA is available for comparison)
2. declares an event when stations agree within the P travel time between them
3. estimates the magnitude from 10 s of waveform per station
4. picks P and S and locates the event by grid search

All inference is hand-written (no ML runtime), with NEON kernels on aarch64.
The target platform is a Raspberry Pi 4 or 5. The input is currently a replay
of AFAD miniSEED files at real time or faster; a live SeedLink client is not
implemented yet.

Named after my pet bird.

## Requirements

- Meson and Ninja, and a C++23 compiler and standard library with
  `std::expected`, `std::print` and `std::mdspan`. Tested with Meson 1.12,
  GCC 16.2, and Clang with libc++ through zig 0.16.
- For exporting the models: [uv](https://docs.astral.sh/uv/) and the
  `archive_pipeline` and `data_downloader` environments that hold the PyTorch
  checkpoints.
- For the Raspberry Pi build: uv (it fetches zig) and, to run the tests on
  x86, `qemu-aarch64`.

## Build and test

```bash
meson setup build-release --buildtype=release
ninja -C build-release
meson test -C build-release
```

The tests that compare against PyTorch, scipy, torchaudio and ObsPy need the fixtures
in `data/fixtures/` and are skipped if these are missing. The end-to-end test
(`demo_m49`) replays the 2025-11-10 Sındırgı M4.9. It requires the event to be
detected, located within 10 km and given a magnitude between 3.9 and 5.9.

Raspberry Pi (static aarch64 binary, `docs/impl/06-raspberry-pi.md`):

```bash
meson setup build-pi --cross-file cross/aarch64-pi.ini --buildtype=release
ninja -C build-pi
meson test -C build-pi --timeout-multiplier 10     # under qemu-aarch64
```

## Data and models

`data/` and `models/` are not tracked. To create them:

```bash
# waveforms: cut a time range from the AFAD archive chunks, one file per station
python3 tools/make_demo_data.py --out data/demo \
    --start 2025-11-10T18:05:00 --end 2025-11-10T18:35:00 \
    ~/Projects/sismokaos/tdvms/afad_raw/{DEMI,MANT,BAND,KAND,KIRK}/*_2025-10-29.zip \
    ~/Projects/sismokaos/tdvms/afad_raw/{ELBA,SEMS}/*_2025-10-21.zip

# weights (AYZW format), filter coefficients, station table, test fixtures
uv run --project ~/Projects/sismokaos/archive_pipeline python tools/export_models.py
uv run --project ~/Projects/sismokaos/data_downloader python tools/export_magnitude.py
uv run --project ~/Projects/sismokaos/archive_pipeline python tools/export_stalta.py
```

Archive chunks are named by their start date, so a given day is usually in a
chunk that starts some days earlier.

## Usage

```
ayzek [options] STATION.mseed...
```

Each file holds one station's HH? channels. The station code is read from the
records and its coordinates from `models/stations.csv`.

```bash
CAT=~/Projects/sismokaos/data_downloader/catalogs/catalog_afad_full_2026-08-30.csv

# three stations at 10x real time, triggers from 18:20 on
build-release/app/ayzek --speed 10 --from 2025-11-10T18:20:00 --catalog $CAT \
    data/demo/DEMI.mseed data/demo/MANT.mseed data/demo/BAND.mseed

# all stations as fast as possible, with S-wave warning times for two cities
build-release/app/ayzek --speed 0 --catalog $CAT \
    --site Istanbul,41.0082,28.9784 --site Bursa,40.1828,29.0667 data/demo/*.mseed

# the same through the demo script
tools/demo.sh          # 3 stations, 10x
tools/demo.sh all      # 7 stations, full speed
```

### Options

| option | default | meaning |
|---|---|---|
| `--models DIR` | `models` | weights directory |
| `--speed X` | 1 | replay speed; 1 = real time, 0 = as fast as possible |
| `--from TIME` | start of data | no triggers before this UTC time; earlier data only updates the station noise baselines |
| `--threshold P` | 0.8 | trigger when `--trigger-windows` consecutive windows reach P |
| `--trigger-windows N` | 8 | see `--threshold`; 8 windows add 3.5 s to a trigger |
| `--instant-threshold P` | 0.9 | also trigger on a single window at or above P; a value > 1 disables this |
| `--release P` | 0.3 | a trigger re-arms after two windows below P |
| `--detector KIND` | `model` | `stalta` replaces the detector with a recursive STA/LTA; picking, association, location and magnitude are unchanged |
| `--sta S`, `--lta S` | 1, 30 | STA/LTA averaging lengths in seconds |
| `--stalta-on R`, `--stalta-off R` | 8, 1.5 | STA/LTA ratio that triggers, and below which a trigger re-arms |
| `--stalta-band LO,HI` | 2,20 | STA/LTA pass band in Hz (4th-order Butterworth high- and low-pass) |
| `--stalta-3c` | | STA/LTA on the energy of all three components instead of the vertical |
| `--no-anchor` | | do not take the P time of a model trigger from the STA/LTA onset |
| `--anchor-on R` | 3 | STA/LTA ratio taken as the onset when anchoring |
| `--anchor-picker`, `--anchor-magnitude` | | also place the picker and early magnitude windows from the anchored P; both make their models worse (`docs/impl/09-sta-lta.md`) |
| `--require-onset` | | a detector trigger also needs an STA/LTA onset; measured, and worse than not requiring one (`docs/impl/09-sta-lta.md`) |
| `--mag-lead S` | 2 | magnitude window start, seconds before P |
| `--step N` | 50 | samples between detector windows (0.5 s at 100 Hz) |
| `--min-stations N` | 2 | station detections needed to declare an event |
| `--no-pick` | | no P/S picker, and so no location |
| `--no-magnitude` | | no magnitude regressor. Otherwise each trigger gets an early estimate, and triggers of declared events a second one at the picked P |
| `--catalog CSV` | | AFAD catalogue export; each alarm is compared with it |
| `--site NAME,LAT,LON` | | also report the S-wave warning time at this place (repeatable) |
| `--scores DIR` | | write every window's detector probability to `DIR/STATION.csv` |
| `--scores-in DIR` | | use the probabilities of an earlier `--scores` run instead of running the detector |
| `--record FILE` | | write all station outputs, for `tools/network_subsets` |
| `--verbose` | | print every station detection, pick and magnitude estimate |
| `--no-color` | | |

### Output

While running, one line per alarm and per magnitude or location update:

```
ALARM    #3   18:21:09.00  earthquake detected by DEMI, MANT  (AFAD MW 4.9 at 18:20:51.00, 18.0 s ago)
MAG      #3   18:21:09.00  magnitude M4.3 from 1 station (DEMI 4.3)
LOCATE   #3   18:22:05.00  located at 39.221N 28.088E, origin 18:20:52.02, rms 0.16 s from 3 stations
```

At the end:

- **One block per alarm:** alarm time and stations, magnitude at the alarm and
  final magnitude, location, comparison with the AFAD event, and the S-wave
  warning time at each station and `--site` (S arrival minus alarm time;
  positive means the alarm came first).
- **Summary:** alarms that match an AFAD event and alarms that do not, AFAD
  events detected and missed (in total and by magnitude band), alarm delay,
  and lists of missed events and unmatched alarms.
- **Station processing:** windows, triggers, picks and magnitude estimates per
  station, and the time per window, pick and estimate.

An alarm matches an AFAD event if its location agrees (epicentre within 30 km,
origin within 5 s), or if the event's predicted P arrival at a station that
raised the alarm is within 6 s of that station's detection. An unmatched alarm
is not necessarily false: the catalogue is incomplete for small events,
especially during aftershock sequences.

## Evaluation tools

| tool | purpose |
|---|---|
| `tools/sweep_trigger.py` | runs a grid of any ayzek options on several datasets and tabulates detections, unmatched alarms, alarm delay and, with `--shifts`, detections corrected for chance matches; uses cached probabilities for the model |
| `tools/network_subsets` | runs the network stage for every subset of stations from one `--record` file; identical to running ayzek on each subset |
| `tools/station_geometry.py` | coverage, alert time, blind zone and location error of station subsets over a grid of epicentres, without waveforms |
| `tools/mseed_dump`, `tools/validate_mseed.py` | record-level inspection, and a decoding check against ObsPy |
| `tools/replay_check`, `tools/gaps_vs_events.py` | gap statistics of the archive, and their relation to earthquake times |
| `tools/scan_encodings.py` | counts the miniSEED data encodings in archive chunks |

Trigger sweeps on one dataset:

```bash
# model: store the detector probabilities once, then sweep the trigger rule
build-release/app/ayzek --speed 0 --no-magnitude --catalog $CAT --scores data/runs/cache data/demo/*.mseed
python3 tools/sweep_trigger.py --catalog $CAT --dataset 'demo:data/runs/cache:data/demo/*.mseed' \
    --set 'threshold=0.9 0.85 0.8' --set 'trigger-windows=1 4 8' --set 'instant-threshold=2 0.9' \
    --shifts 127,263,391 --out sweep.csv

# STA/LTA (no cache)
python3 tools/sweep_trigger.py --catalog $CAT --args '--detector stalta' --dataset 'demo::data/demo/*.mseed' \
    --set 'stalta-on=4 6 8' --set 'lta=10 30' --set 'stalta-3c=off on' --shifts 127,263,391 --out stalta.csv
```

## Results

Replays of AFAD recordings with the default settings, scored against the AFAD
catalogue (events within 250 km):

| replay | stations | catalogue events detected | alarms without a catalogue event | main event |
|---|---:|---|---|---|
| Sındırgı, 2025-11-10, 30 min | 3 | 7 of 9 | 2 of 9 | Mw 4.9: alarm 18.0 s after origin with M4.3, final M4.6, located 5.3 km off |
| Marmara Sea, 2025-04-23, 3 h | 8 | 44 of 57 | 24 of 68 | Mw 6.2: alarm 11.0 s after origin, M4.6 at +14 s, final M4.8, located 3.1 km off |
| Sındırgı, 2025-08-10, 3 h | 6 | 39 of 66 | 7 of 46 | Mw 6.1: alarm 17.0 s after origin, M4.9 at +17.5 s, final M4.7, located 32.5 km off |

- **Trigger settings.** Chosen on the first two replays. The third was not used
  for the choice. There, the earlier rule (one window at p ≥ 0.9) detected 24
  of 66 events with 2 alarms without a catalogue event
  (`docs/impl/05-pipeline.md`).
- **Alarm time.** Set by the P travel time to the second station, not by
  computation.
- **Magnitude.** Saturates for large events (final M4.7 for the Mw 6.1, M4.8
  for the Mw 6.2). The regressor was trained on 10 s windows with few events
  above M5.5 (`docs/impl/07-magnitude.md`).
- **Unmatched alarms.** Some are probably uncatalogued aftershocks. Of the 24
  in the Marmara replay, 9 were located and 5 of those within 15 km of the
  Mw 6.2.
- **Against STA/LTA.** STA/LTA settings were chosen on the first two replays
  to match the model's unmatched alarms. On the held-out replay that setting
  detects 25 events against the model's 39, with 9 unmatched alarms against 7.
  To detect as many events as the model, STA/LTA needs about three times the
  unmatched alarms. It alarms 1–1.3 s earlier for the Mw 6.1 and 6.2, and
  costs about 0.001 ms per update against 5 ms (`docs/impl/09-sta-lta.md`).
- **Speed.** On a 12-thread x86 laptop: 5–6 ms per detector window (three
  models, every 0.5 s per station), 13 ms per picker run, and 300–400 ms per
  magnitude estimate. The aarch64 NEON build passes the same agreement tests
  under qemu. Timing on a Raspberry Pi has not been measured yet.

## Documentation

- `docs/DESIGN.md`: the architecture and the reasoning behind it
- `docs/IMPLEMENTATION.md`: an overview of what was built, with an index of
  `docs/impl/01`–`09` (SIMD kernels, weights format, conditioning, models,
  pipeline, Raspberry Pi, magnitude, station selection, STA/LTA benchmark)

## Layout

```
src/            ring buffer, miniSEED decoding, reordering, SIMD kernels, layers, models, DSP, magnitude
src/pipeline/   ingest, per-station processor, network association and location, recording
app/            the ayzek binary
tests/          unit tests, agreement with PyTorch/scipy/torchaudio, end-to-end check
tools/          model export, data preparation, evaluation
cross/          zig toolchain wrappers and the Raspberry Pi cross file
docs/           design and implementation notes
```
