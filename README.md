# ayzek

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23672495.svg)](https://doi.org/10.5281/zenodo.23672495)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/license-GPL--3.0--or--later-blue.svg)](LICENSE)

ayzek is a real-time earthquake detection pipeline for early warning, written
in C++23 with no machine-learning runtime. It reads continuous three-component
waveforms from a set of stations, detects earthquakes at each station with a
learned detector, declares an event when stations agree within the P travel
time between them, and reports a first magnitude from the peak displacement of
the P wave within seconds of the alarm. All inference is hand-written and
verified against the PyTorch models it was exported from. The same source
builds for x86-64 Linux, the Raspberry Pi 4 and 5 (static aarch64 binary with
NEON kernels) and Windows.

The repository also holds the evaluation method: a replay protocol that scores
detectors at their operating point with chance correction, paired tests and
catalog-quiet data, and the held-out test protocol under which the method was
frozen as release v1.0.0 and then tested on three networks it had not been
tuned on. Every number in the accompanying paper is produced by a script in
`tools/` and checked against `bench/claims.json`.

The name is that of the author's pet bird.

## Contents

- [What it does](#what-it-does)
- [Status](#status)
- [Building](#building)
- [Running](#running)
- [Evaluation](#evaluation)
- [Results on the development data](#results-on-the-development-data)
- [Platforms](#platforms)
- [Data sources and licences](#data-sources-and-licences)
- [Documentation and layout](#documentation-and-layout)
- [Citing](#citing)
- [License](#license)

## What it does

For each station, ayzek ingests miniSEED records (from files replayed at any
speed, or live from a SeedLink server), conditions the three components, and
runs one of four detectors over a 0.5 s grid:

| detector | what it is |
|---|---|
| `6s` (default) | a convolutional–recurrent network with self-attention, 142 059 parameters, three seeds averaged, scoring 6 s windows; its P time is anchored to a recursive STA/LTA onset within the window |
| `transformer` | a streaming onset transformer with one output per 0.1 s and, with its geometry head, a per-station distance (`docs/impl/14-transformer.md`) |
| `stalta` | a recursive STA/LTA trigger, numerically identical to ObsPy's |
| `filterpicker` | a port of Lomax's FilterPicker5, giving identical picks |

Station detections are associated when their P times differ by no more than
the inter-station P travel time plus a slack, and an event is declared when two
stations agree. The event is then located, first from the detections' P times
alone and later from P and S picks (or from the transformer's distances), and
its magnitude is estimated from Pd with the probability of each of three size
classes (`docs/impl/17-pd-magnitude.md`). When a catalog is given, every alarm
is scored against it, with chance matches corrected by shifting the catalog in
time.

Replays are deterministic: repeated and concurrent runs of the same data give
byte-identical reports, and a recorded run can be replayed through the network
stage for every subset of stations.

## Status

Release **v1.0.0** (30 October 2026) is the method freeze. The held-out
protocol (`docs/protocol-heldout.md`) fixes, by that release, every default of
the pipeline, the models by their SHA-256, the settings of every baseline
detector it is compared with, and the rules that select the test data. The
freeze is timestamped by the Zenodo deposit of the release, and no waveform
from a test pool was fetched before it. The test pools are the Greek networks
served by the National Observatory of Athens, New Zealand's GeoNet, and the
KOERI network the detector was trained on, held out in time. The selections
made after the freeze are in `docs/heldout-selections/`; the test runs and
their results are reported in the paper, and the scripts that produce them are
`tools/run_heldout.py` and the tools it calls.

The tool that selects and fetches the test data lives in a separate repository
and is attached to the release as `ayzek-heldout-440a2d9.tar.gz`.

## Building

Requirements: Meson, Ninja and a C++23 compiler and standard library with
`std::expected`, `std::print` and `std::mdspan` (tested with Meson 1.12,
GCC 16.2, and Clang with libc++ through zig 0.16). The model weights are in the
repository, so neither Python nor an ML runtime is needed to build or run. The
cross builds use [uv](https://docs.astral.sh/uv/) to fetch zig.

```bash
meson setup build-release --buildtype=release
ninja -C build-release
meson test -C build-release
```

The agreement tests against PyTorch, SciPy, torchaudio and ObsPy read fixtures
from `data/fixtures/` and are skipped when those are absent. The end-to-end
test replays the 2025-11-10 Sındırgı M4.9 and requires the event to be
detected, located within 10 km and given a magnitude between 3.9 and 5.9.

The Python analysis tools have their own tests, which need neither a build nor
data:

```bash
uv run --project tools --with pytest pytest tools/tests -q
```

## Running

```
ayzek [options] STATION.mseed...
```

Each file holds one station's HH channels at 100 Hz; the station code is read
from the records and its coordinates from the station table. The quickest data
to try is the public KOERI set attached to the releases (4 MB):

```bash
mkdir -p data/demo_ko
curl -L https://github.com/hogib/ayzek/releases/download/v0.1.0-beta.6/ayzek-demo_ko.tar.gz \
    | tar xz -C data/demo_ko
build-release/app/ayzek --speed 10 --from 2025-11-10T18:20:00 \
    --catalog data/demo_ko/demo_ko.csv data/demo_ko/*.mseed
```

While running, ayzek prints one line per alarm and per magnitude or location
update:

```
ALARM    #3   18:21:09.00  earthquake detected by KULA, SIMA  (AFAD MW 4.9 at 18:20:51.00, 18.0 s ago)
PLOC     #3   18:21:10.50  P times locate it at 39.240N 28.119E +-23 km, origin 18:20:52.18, rms 0.01 s from 3 stations
PDMAG    #3   18:21:12.50  Pd magnitude M4.7 ± 0.2 (M<4 0%, 4-5.5 100%, M>=5.5 0%) from 3 stations (3 above noise), windows up to 5 s, distances from the P-time location
LOCATE   #3   18:22:04.50  located at 39.206N 28.137E, origin 18:20:52.66, rms 0.21 s from 5 stations
```

The final report gives one block per alarm (stations, magnitude at the alarm
and at the end, location, the matching catalog event, S-wave warning time at
each station and at every `--site`), a summary of matched and unmatched alarms
and of detected and missed catalog events by magnitude band, and per-station
processing counts and timings. An alarm matches a catalog event when its
epicentre lies within 30 km and its origin within 5 s, or when the event's
predicted P arrival at a station that raised the alarm lies within 6 s of that
station's detection. An unmatched alarm is not necessarily false: catalogs are
incomplete for small events, especially in aftershock sequences.

The options most often needed:

| option | meaning |
|---|---|
| `--speed X` | replay speed: 1 real time, 10 ten times faster, 0 as fast as possible |
| `--from TIME` | no triggers before this UTC time; earlier data still warms the filters |
| `--catalog CSV` | score every alarm against this catalog (AFAD export format) |
| `--detector KIND` | `6s`, `transformer`, `stalta` or `filterpicker` |
| `--stations FILE`, `--pd-gains FILE` | the station table and Pd gains of another network |
| `--site NAME,LAT,LON` | report the S-wave warning time at this place |
| `--scores DIR`, `--scores-in DIR` | cache the detector's window probabilities, and replay from the cache |
| `--record FILE` | record every station output for `tools/network_subsets` |
| `--seedlink HOST:PORT --streams NET.STA,...` | live input instead of files (`docs/impl/19-seedlink.md`) |

The full table, including every detector's trigger settings, the association
and location parameters and the live-input options, is in
[`docs/OPTIONS.md`](docs/OPTIONS.md).

## Evaluation

The replay protocol and the tools that implement it are described in the
paper and in `docs/impl/05-pipeline.md`, `09-sta-lta.md`, `10-benchmarks.md`
and `13-claims.md`. In brief: trigger settings are chosen on tuning data by
rules stated in advance; a declared event is credited to the one catalog event
it fits best; the chance level is measured by shifting the catalog by ±127,
±263 and ±391 s and removed from the detection rate; detectors are compared on
per-event outcomes by an exact McNemar test and a paired bootstrap; and false
alarms are counted on windows with no cataloged event within 250 km of any
station.

| tool | purpose |
|---|---|
| `tools/sweep_trigger.py` | a grid of ayzek options over several datasets, with detections, unmatched alarms, delays and chance-corrected counts; uses cached probabilities |
| `tools/select_setting.py` | applies the selection rules (equal unmatched alarms, equal detections, equal quiet alarms) to a sweep |
| `tools/event_outcomes.py`, `tools/significance.py` | per-event outcome tables and the paired tests, equivalence test and depth-restricted repeats |
| `tools/quiet_sweep.py`, `tools/quiet_without_stations.py` | alarm rates on quiet windows, with and without the noisiest stations |
| `tools/stalta_station_ratios.py` | the per-station STA/LTA baseline |
| `tools/seisbench_scores.py` | EQTransformer and PhaseNet (SeisBench 0.12.6) scored on ayzek's window grid, for `--scores-in` |
| `tools/network_subsets` | the network stage for every subset of stations of one recorded run |
| `tools/run_heldout.py` | a held-out test end to end, from a selection to every endpoint of the protocol, with a manifest |
| `tools/reproduce_all.sh`, `tools/verify_claims.py` | every result of the paper, checked against `bench/claims.json` |
| `tools/scorecard.py`, `bench/ayzek_bench` | detection, magnitude and location on the fixed datasets; time per stage and per layer |

Smaller utilities (miniSEED inspection and validation, gap scans, a SeedLink
test server, KOERI data fetching, map and figure scripts) are listed in
`docs/IMPLEMENTATION.md`.

## Results on the development data

Every setting of the pipeline was chosen on replays of the AFAD archive and the
public KOERI sets; these are development results, not tests. With the default
settings, scored against the AFAD catalog within 250 km of the stations:

| replay | stations | catalog events detected | unmatched alarms | main event |
|---|---:|---|---|---|
| Sındırgı, 2025-11-10, 30 min | 3 | 7 of 9 | 2 | Mw 4.9: alarm 18.0 s after origin |
| Marmara Sea, 2025-04-23, 3 h | 8 | 44 of 57 | 24 | Mw 6.2: alarm 11.0 s after origin, S wave at the nearest station 2.1 s later |
| Sındırgı, 2025-08-10, 3 h | 6 | 39 of 66 | 7 | Mw 6.1: alarm 17.0 s after origin |

- Against STA/LTA tuned to the same number of unmatched alarms on the tuning
  replays, the learned detector detected 36.4 chance-corrected events of the
  66 of the Sındırgı Mw 6.1 sequence against 21.7 (difference 14.7, 95 %
  interval 5.9–23.7). To detect as many events as the learned detector,
  STA/LTA needed about three times its unmatched alarms. STA/LTA alarms
  1.0–1.3 s earlier on the two largest events and costs about 0.0005 ms per
  update against 5.2 ms for the detector ensemble.
- Over 102.8 h of recordings with no cataloged event within 250 km of any
  station, the learned detector raised 16 alarms (3.74 a day) and the STA/LTA
  setting above 323 (75.4 a day). Nine of the 16 came from one station's
  periodic instrumental pulse.
- Causal PhaseNet, run through the same pipeline with pre-trained weights,
  matched the learned detector on both development sequences; causal
  EQTransformer did not.
- On the public KOERI sets: the Sındırgı Mw 4.9 is declared 18.0 s after
  origin and receives M4.7 at +3.5 s and M4.8 at the end; the Marmara Mw 6.2
  is declared 10.5 s after origin and receives M6.0; on 5.8 h of quiet data the
  learned detector raised 1 alarm against 25 for STA/LTA. The Marmara set
  overlaps the detector's training data and is a demonstration only.
- Speed: on a Ryzen 5 5600, 5.2 ms per detector window (three models, one
  window every 0.5 s per station), and one thread sustains 95 stations. On a
  Raspberry Pi 4 a window costs 36.3 ms and four cores sustain 52 stations,
  with the same alarms at the same times.

The complete development record, with the trigger grids, the selection of
every baseline and the quiet-window classification, is in the paper's
supplement and in `docs/impl/`.

## Platforms

- **Linux x86-64**: the development build above. The tests of the held-out
  protocol run the static build from `cross/x86_64-linux.ini` (zig and musl),
  whose reports equal the development build's on every development replay.
- **Raspberry Pi 4 and 5**: `meson setup build-pi --cross-file cross/aarch64-pi.ini
  --buildtype=release`; the test suite runs under `qemu-aarch64` with
  `--timeout-multiplier 10` (`docs/impl/06-raspberry-pi.md`).
- **Windows 10 and 11**: `ayzek.exe` is cross-built on Linux with
  `cross/x86_64-windows.ini` and needs no runtime DLLs. The build is linked and
  attached to the releases, but the numerical tests have not yet been run on a
  Windows machine; live SeedLink input is not available in this build
  (`docs/impl/11-windows.md`, `docs/OPTIONS.md`).

## Data sources and licences

`models/` is tracked (exported weights in the AYZW format, filter coefficients,
the station table and the Pd calibration, about 5 MB). `data/` is not tracked.

The public datasets attached to the releases are KO network waveforms of the
Kandilli Observatory and Earthquake Research Institute (KOERI), Boğaziçi
University, open under citation:

> Kandilli Observatory And Earthquake Research Institute (KOERI), Boğaziçi
> University (1971). *Kandilli Observatory And Earthquake Research Institute
> (KOERI)*. International Federation of Digital Seismograph Networks.
> https://doi.org/10.7914/SN/KO

The AFAD waveforms of the development replays are available from AFAD's data
portal (https://tdvms.afad.gov.tr) on registration and cannot be
redistributed; the catalog excerpts in `tests/catalogs/` are AFAD's
(https://deprem.afad.gov.tr/). The test waveforms of the held-out protocol
(the Greek networks HL, HT, HA, HP and HC through the NOA EIDA node, NZ through
GeoNet, KO through KOERI) are not redistributed either; each set's request
record and file hashes are published with the protocol, so that the same data
can be fetched again and checked. `tools/fetch_ko_data.py` fetches KO
waveforms into the replay layout; `tools/make_demo_data.py` cuts replays from
AFAD archive chunks.

## Documentation and layout

- `docs/DESIGN.md`: the architecture and the reasoning behind it.
- `docs/IMPLEMENTATION.md`: what was built, with an index of the notes in
  `docs/impl/` (SIMD kernels, weights format, conditioning, models, pipeline,
  Raspberry Pi, magnitude, station selection, STA/LTA, benchmarks and the
  scorecard, Windows build, environment, claims, the streaming transformer,
  location from its geometry head, alarm assessment, Pd magnitude,
  FilterPicker, live SeedLink input).
- `docs/protocol-heldout.md`: the frozen held-out test protocol;
  `docs/release-notes-freeze.md`: what v1.0.0 fixes, with every hash;
  `docs/heldout-selections/`: the test sets selected after the freeze.
- `docs/OPTIONS.md`: every command-line option.

```
src/            ring buffer, miniSEED decoding, reordering, SeedLink client, SIMD kernels,
                layers, models, DSP, magnitude
src/pipeline/   ingest (replay and live), per-station processing, association, location,
                alarm assessment, recording
app/            the ayzek binary
tests/          unit tests, agreement with PyTorch/SciPy/torchaudio/ObsPy, end-to-end check
tools/          model export, data preparation, evaluation, the held-out runner
bench/          ayzek_bench, claims.json and reference results
cross/          zig toolchain wrappers and the Linux, Raspberry Pi and Windows cross files
docs/           design, implementation notes, protocol and release notes
```

## Citing

Every release is archived on Zenodo. The DOI below resolves to the latest
release; each release's Zenodo page gives its own DOI, which is the one to cite
for a specific version (v1.0.0: 10.5281/zenodo.23806743). `CITATION.cff` holds
the same metadata.

> Bolat, H. O., Çimen, S., Bolat, B., & Yıldırım, T. ayzek: a self-contained
> real-time earthquake detection and early-warning pipeline. Zenodo.
> https://doi.org/10.5281/zenodo.23672495

```bibtex
@software{ayzek,
  author    = {Bolat, H. O{\u{g}}uz and {\c{C}}imen, Sibel and Bolat, B{\"u}lent and Y{\i}ld{\i}r{\i}m, T{\"u}lay},
  title     = {ayzek: a self-contained real-time earthquake detection and early-warning pipeline},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23672495},
  url       = {https://github.com/hogib/ayzek}
}
```

The paper describing ayzek and the held-out tests is under preparation and
will be cited here once published.

## License

GNU General Public License, version 3 or (at your option) any later version;
the full text is in `LICENSE`. The FilterPicker port derives from Anthony
Lomax's C code (GPL-2.0-or-later). Releases up to and including v0.1.0-beta.3
were published under the Apache License 2.0 and remain available under it.
