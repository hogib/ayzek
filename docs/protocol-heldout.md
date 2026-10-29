# Held-out tests: protocol (DRAFT)

Status: **draft, not frozen.** Values marked *[...]* are still open. The
protocol is frozen when this file, the release it names and the hashes below
are committed publicly and timestamped (Zenodo). No waveform from a test pool
is fetched before then.

## 1. What is frozen

- **ayzek**: the release tagged at the freeze (commit and tag recorded here),
  with every default as it is in that commit. The 6 s detector keeps its rule
  (0.8 for 8 windows, or 0.9 once). The baselines' settings are those selected
  on the development data (section 6).
- **Models**: the SHA-256 of every file in `models/` (detector, picker, Pd
  relation, gains, terms, noise levels, transformer) is listed in the release
  notes.
- **The transformer** (its own paper): its weights, threshold, dt restart,
  station count and P-time uncertainty as in that release.
- **The binary**: the tests run the release's static x86-64 binary
  (`build-linux`, zig/musl). On every development replay, for the learned
  detector, STA/LTA, FilterPicker and the transformer, its reports equal the
  development build's apart from compute times and the order of
  station-thread lines (checked 2026-10-29; repeated at the tag).
- **Analysis code**, at the frozen commits: `tools/run_heldout.py` (the test
  end to end, from a selection to every endpoint, with its manifest) and the
  tools it runs, `tools/sweep_trigger.py`, `tools/event_outcomes.py`,
  `tools/significance.py`, `tools/quiet_sweep.py`,
  `tools/quiet_without_stations.py` with `tools/network_subsets`,
  `tools/seisbench_scores.py` with its pinned weight versions,
  `tools/select_setting.py` and `tools/stalta_station_ratios.py` (the
  selection of the settings); `tools/heldout_settings.json` (every frozen
  setting) and `tools/detector_station_split.csv` (the 6 s detector's station
  split, endpoint 6); and in the held-out repository `heldout select`,
  `probe`, `pull` and `catalog`.

## 2. Development and test data

- **Development** (used for every choice; no test result is reported from it):
  the AFAD (TU) archive replays (Sındırgı 30 min, Marmara 3 h, Sındırgı
  Mw 6.1 3 h, Simav Mw 4.8 3 h, the 102.8 h of quiet windows) and the KOERI sets
  demo_ko, marmara_ko and quiet_ko. marmara_ko overlaps the detector's
  training data and is never used to compare detectors.
- **Test networks**, each a separate result, never pooled across networks:
  - **Greece**: all five Greek networks served by the NOA EIDA node, HL (NOA),
    HT (AUTH), HA (University of Athens), HP (University of Patras) and HC
    (Crete), as one network group;
  - **New Zealand**: NZ (GeoNet);
  - **KO** (KOERI): held out in time only; both models trained on its stations.
    Its retrospective pool is two months, which the pre-pool rate of
    qualifying sequences (37 from 2024-01-01 to 2026-08-29, 1.16 a month, by
    EMSC's catalogue and the stations usable in July-August 2026) puts at about
    two sequences, one of them "now": the KO retrospective result is expected
    to be far under its target and is reported as such, the Greek and New
    Zealand results carrying the headline.
  - **Prospective**: New Zealand and KO again, on data recorded after the
    freeze (section 10), each a result of its own.
- **Stations**: three HH channels at 100 Hz, open over the pool period, with at
  least 90% data coverage over it (from the availability services, before any
  waveform is fetched), whose components' orientations can be rotated to north
  and east (the determinant of their three unit vectors, from azimuth and dip,
  at least 0.5 in absolute value; NZ.MKTZ from 2025-09-30 fails it).
  Coordinates from each set's StationXML (`--stations`).
- **Pd gains, station terms and noise levels**:
  - Greece and New Zealand: the gains of each set's StationXML (`pd_gain.csv`,
    `--pd-gains`). A vertical channel epoch gets a gain only if its velocity
    response divided by the stated sensitivity lies within 0.8-1.25 at 0.1,
    0.2, 0.5, 1, 2 and 5 Hz (Pd divides by that one number); otherwise the
    station is not measured for Pd, and the held-out tool records the reason
    (16 of the 166 Greek vertical epochs of the tuning sets fail it, none of
    KO's or New Zealand's). With `--pd-gains` the Pd station terms are zero and
    no station noise levels apply (no dead-channel or coda exclusion, section
    11).
  - KO: runs without `--pd-gains` (`tools/run_heldout.py`). It keeps the
    gains, station terms and noise levels of its calibration in `models/`
    (`pd_gain.csv`, `pd_terms.csv`, `pd_noise.csv`, keyed by station code),
    with the exclusions that calibration flagged; a KO station without a gain
    epoch there covering the replay is not measured for Pd.

## 3. Pools and the reserve

- **Pool periods**: KO, from 2026-08-30 (after the last training
  window of the transformer and the Pd calibration, 2026-08-29) to the freeze
  date; Greece and New Zealand, from 2023-01-01 to the freeze date.
- **Catalogues**: EMSC's fdsnws/event for KO (KOERI runs no event service, and
  AFAD's event API, the first choice, answered every request of 2026-10-29
  with an empty list), NOA's fdsnws/event for the Greek networks and GeoNet's
  for New Zealand; each catalogue's preferred origin and magnitude as given.
  Before each catalogue request the fetch tool asks the same service for a
  fixed pre-pool window known to contain events (its canary) and stops if
  fewer than about 40% of them come back; `dataset.json` records the check.
- **Usable stations**: those of section 2 over the whole pool period (three
  HH channels at 100 Hz, at least 90% coverage over the period).
- **Sequences**: in time order, every catalogue event of M >= 4.5 within the
  pool period with at least five usable stations within 150 km of its
  epicentre starts a sequence, unless it lies within 3 h after the first event
  of an earlier sequence with which its stations share one (it then belongs to
  that sequence). A sequence is a replay of 3 h starting 10 min before its
  first event, on the 8 nearest usable stations within 150 km of it (one per
  station code). Catalogue events are scored within 250 km of the centroid of
  the replay's stations, over the replayed span, as the development replays
  were.
- **Reserve split**: the qualifying sequences of each network, in time order,
  are assigned to "now" or "reserve" by a coin with seed 42: Python's
  `random.Random(42)`, one draw per sequence, below 0.5 is "now". Reserve
  sequences are not fetched until a later protocol, committed in the same way,
  names them.
- **Stopping rule**: "now" sequences are taken in time order until the
  scorable catalogue events (as scored above, counted from the catalogue and
  the station list before any waveform is fetched) reach **240 on KO and 150
  on Greece and on New Zealand**; the sequence that crosses the target is
  included whole. If a network's "now" sequences do not reach it, all of them
  are used and the result is reported as under the target.
- **Quiet data**: 24 windows of 6 h per network, drawn with Python's
  `random.Random(43)`: each draw is a start time on the minute, uniform over
  the pool period less 6 h and a 30 min margin at either end, then one of the
  selected "now" sequences' station sets, uniform. A draw is accepted if no
  catalogue event lies within 250 km of any of its stations during the window
  or 30 min on either side, it does not overlap an accepted window that shares
  a station, and the availability service gives every station at least 90%
  of the window; draws continue until 24 are accepted. If 10,000 draws give
  fewer (as expected for New Zealand), the draws start again from seed 43
  with catalogue events below M2 farther than 100 km from every station
  disregarded. The first 30 min of each window warm the filters and are not
  scored.
- **Tool**: `heldout select` of the held-out repository applies these rules
  and writes the selection (`selection.json`: every qualifying sequence with
  its coin, stations and scorable count, the stopping point, every quiet pass)
  and the set files of the selected sequences and quiet windows; reserve
  sequences get none. It refuses a pool period without `--allow-test-period`,
  which is used only after the freeze. Its dry runs on 2022 (Greece, New
  Zealand) and 2026-07-01 to 2026-08-29 (KO), before the pools, are kept with
  the protocol.

## 4. Per-station STA/LTA (secondary baseline)

Each test station's own STA/LTA trigger ratio, from its own data before the
pool: four 6 h windows per network group, drawn with seed 20221231
(`tools/make_tuning_sets.py` of the held-out repository; Greece and New Zealand
in 2022, KO from 2026-07-01 to 2026-08-29), with every station of the group
that has three 100 Hz HH channels and 90% coverage. The setting is the
network-wide equal-unmatched one (section 6) with the trigger ratio replaced by
the lowest ratio of the grid (3, 4, 5, 6, 8, 10, 12, 15) at which the station's
trigger rate on its quiet data over the four windows does not exceed its
network group's target: the median station rate of the network-wide setting
(ratio 8) on the same windows, 0.845 triggers per quiet hour for Greece, 0.548
for KO and 0.318 for New Zealand. A station above the target at every ratio
gets 15; a test station without tuning data keeps the network-wide ratio.

Changed on 2026-10-29, before the freeze: the target was first the median
station rate of the network-wide setting on the AFAD quiet windows (0.307 per
quiet hour), with the grid ending at 10. Those windows were chosen for having
no catalogued event within 250 km, the tuning windows are random, and against
that target 117 of 150 Greek and 93 of 144 KO stations ended at the grid's top
ratio, so that the baseline measured how quiet the AFAD windows were chosen to
be rather than per-station tuning (`data/runs/stalta_stations_afad_target`).

Quiet data: on each tuning window, the time outside [D/8 - 10 s, D/6 + 120 s]
after the origin of every catalogued event
visible at the station (epicentral distance D in km; visible within 50 km at
any magnitude, 150 km from M2, 400 km from M3 and 1500 km from M4.5, as in
`tools/compare_detectors.py`), after the first 30 min of the window. A rate is
the station's
detections in its quiet time per quiet hour. `tools/stalta_station_ratios.py`
computes the target and the ratios (ayzek's `--stalta-on-stations`). These
pre-pool data are fetched before the freeze, for this tuning only.

## 5. Scoring

As on the development data: chance correction with catalogue shifts of +-127,
+-263 and +-391 s; an alarm credited to the catalogue event it fits best and an
event to the first alarm credited with it; paired comparisons on per-event
outcomes, with the exact McNemar test (on raw outcomes) and a paired bootstrap
of the difference in chance-corrected detections, 10 000 resamples within each
sequence.

## 6. Detectors and their settings

Selected on the development event replays (the four AFAD sequences and demo_ko)
at the learned detector's counts on those replays (`data/runs/tune_dev/run.sh`,
`tools/select_setting.py`), by three rules:

- **equal unmatched alarms**: the most chance-corrected detections among the
  settings with at most the learned detector's unmatched alarms;
- **equal detections**: the fewest unmatched alarms among the settings with at
  least its chance-corrected detections;
- **equal quiet alarms**: the most chance-corrected detections among the
  settings that raise at most its alarms on the 23 AFAD quiet windows (102.8 h;
  `data/runs/tune_dev/run_quiet.sh`). Unmatched alarms in an aftershock sequence
  are partly uncatalogued earthquakes, more so for a sensitive detector, so the
  first rule places a baseline at a higher true false-alarm rate than the learned
  detector; this rule matches the rate on data without catalogued seismicity.
  For the published models it is applied within the weights selected by the
  first two rules.


Selected 2026-10-29 on the five development event replays (197 catalogued
events; the learned detector: 128.2 chance-corrected detections, 52 unmatched
alarms, 16 alarms on the 23 AFAD quiet windows, 3.74 a day). "Corrected" and
"unmatched" are the setting's totals over the five replays; "quiet" its alarms
on the quiet windows (102.8 h). The published models' rules 1 and 2 ran over
the original, STEAD, INSTANCE and ETHZ weights of EQTransformer and the STEAD,
INSTANCE and ETHZ weights of PhaseNet; rule 3 within the weights they chose.
`tools/heldout_settings.json` holds every setting as ayzek options.

| detector | rule | setting | corrected | unmatched | quiet |
|---|---|---|---|---|---|
| learned (6 s) | fixed | 0.8 x 8 windows, or 0.9 once | 128.2 | 52 | 16 |
| STA/LTA | equal unmatched | 1-10 Hz, STA 0.5 s, LTA 10 s, ratio 8, vertical | 110.3 | 47 | 63 |
| STA/LTA | equal detections | 1-10 Hz, STA 1 s, LTA 10 s, ratio 4, vertical | 131.6 | 113 | 224 |
| STA/LTA | equal quiet alarms | 1-10 Hz, STA 0.5 s, LTA 10 s, ratio 10, three components | 64.8 | 27 | 12 |
| STA/LTA, per station | section 4 | the equal-unmatched setting with each station's own ratio (3 to 15) | | | |
| FilterPicker | equal unmatched | windows 6/10/0.4 s, thresholds 20 and 10 | 102.5 | 44 | 291 |
| FilterPicker | equal detections | windows 6/10/0.4 s, thresholds 15 and 8 | 130.9 | 101 | 954 |
| FilterPicker | equal quiet alarms | windows 6/10/0.4 s, thresholds 6 and 15 | 52.0 | 9 | 9 |
| PhaseNet, causal | equal unmatched | INSTANCE, 0.8 x 1 | 119.7 | 52 | 49 |
| PhaseNet, causal | equal detections | INSTANCE, 0.5 x 4 | 128.3 | 70 | 29 |
| PhaseNet, causal | equal quiet alarms | INSTANCE, 0.7 x 2 | 126.1 | 62 | 15 |
| EQTransformer, causal | equal unmatched | INSTANCE, 0.7 x 8 | 62.1 | 49 | 7 |
| EQTransformer, causal | equal detections | INSTANCE, 0.3 x 8 | 128.5 | 183 | 26 |
| EQTransformer, causal | equal quiet alarms | INSTANCE, 0.6 x 2 | 94.4 | 130 | 14 |
| PhaseNet, look-ahead (offline) | equal unmatched | INSTANCE, 0.7 x 4 | 123.2 | 52 | 12 |
| PhaseNet, look-ahead (offline) | equal detections | INSTANCE, 0.6 x 8 | 135.9 | 61 | 14 |
| PhaseNet, look-ahead (offline) | equal quiet alarms | INSTANCE, 0.7 x 1 | 138.1 | 74 | 16 |
| EQTransformer, look-ahead (offline) | equal unmatched | INSTANCE, 0.6 x 4 | 113.0 | 50 | 11 |
| EQTransformer, look-ahead (offline) | equal detections | original, 0.7 x 8 | 140.0 | 88 | 13 |
| EQTransformer, look-ahead (offline) | equal quiet alarms | original, 0.6 x 2 | 153.2 | 116 | 16 |
| transformer (its own paper) | fixed | the model's defaults; geometry locator sd scale 2.1, sigma_p 2.75 s | 121.9 | 70 | 6 |

The published models' "threshold x windows" is the trigger rule on the model's
window probability, the instant branch off. Quiet alarms of the learned
detector and the baselines are network alarms on the 23 windows.

## 7. Endpoints

**Primary, per network**: the difference in chance-corrected detections between
the learned detector and STA/LTA at equal unmatched alarms, over the network's
"now" sequences, with its 95% paired bootstrap interval. The same difference
for STA/LTA at equal quiet alarms is the first secondary endpoint.

**Secondary, per network**:
1. every other detector against the learned detector, as in the primary;
2. PhaseNet (causal) against the learned detector: equivalence by two
   one-sided tests at 5% with a margin of +-10% of the scorable events
   (`tools/significance.py --equivalence`);
3. alarms per day on the quiet windows, every detector, with exact Poisson
   intervals; and the same with the network's two stations of most STA/LTA
   (equal unmatched alarms) station detections on its quiet windows left out,
   for every detector alike (`tools/quiet_without_stations.py --exclude-top 2
   --by stalta`), so that the comparison does not rest on two noisy stations;
4. the neighbourhood of the learned detector's rule: the adjacent settings of
   its grid (threshold 0.75 and 0.85; 6 and 10 windows), reported beside it;
5. detection before S: the share of the detections a match uses (the stations
   that raised the alarm, declared by its time) that were declared before the
   predicted S at their station, and the primary comparison restricted to
   events detected before S, those whose every such detection was. S is
   predicted from the catalogue hypocentre in the uniform half-space of the P
   prediction (3.5 and 6 km/s). The match itself keeps its 6 s tolerance on
   the predicted P, which the half-space's error at regional distances needs.
   ayzek prints `pre-S  N of K alarm detections ...` after each matched event
   and the count of events detected before S in its summary;
   `tools/event_outcomes.py` writes the per-event outcome (`detected_pre_s`)
   and `tools/significance.py --outcome detected_pre_s` the comparison;
6. KO only: the primary comparison split by the 6 s detector's station split
   (training against test-split stations);
7. magnitude (Pd): error and bias by band, class accuracy at the first estimate
   and the final one, and the share of alarmed events with an estimate, for
   events no deeper than 40 km (deeper ones reported apart);
8. latency: alarm delay after origin, and warning time before S.

## 8. Deviations

A bug found after the test data are fetched is fixed, and the results are
reported both ways, labelled as a deviation. Nothing is re-tuned on test data.

## 9. Records

The fetched sets (miniSEED, StationXML, catalogue excerpts, `dataset.json`) are
published where the networks' licences allow, each under its own licence (CC BY
4.0 for the NOA node's services, CC BY 3.0 NZ for GeoNet's) with the networks'
citations. Every run's command, commit and output is kept. Nothing is published
before the data centre has confirmed in writing that its licence covers the
redistribution: for the Greek networks this confirmation from NOA
(eida@noa.gr), requested on 2026-10-29, is still outstanding, and so is a check of EMSC's terms for the
KO catalogue excerpts.

## 10. Prospective tests (New Zealand and KO)

Tests on data that did not exist when the method was frozen, one per network,
each reported as its own confirmatory result and never pooled with the
retrospective test of the same network or with any other network.

- **Data**: NZ (GeoNet) and KO (KOERI), recorded from the freeze (the UTC time
  of this protocol's public commit) on. Stations by the rule of section 2,
  their availability checked over the span before any waveform is fetched.
- **Sequences**: by the rule of section 3 (M >= 4.5 in GeoNet's catalogue, at
  least five usable stations within 150 km, not within 3 h after an earlier
  qualifying event; 3 h from 10 min before), every qualifying sequence in time
  order, with no reserve split.
- **Span and stopping rule**: from the freeze until the scorable catalogue
  events of the network's sequences reach its target (150 for New Zealand, 240
  for KO), the sequence that crosses the target included whole, or 12 months,
  whichever comes first; if 12 months pass under the target, every sequence is
  used and the result is reported as under it.
- **Fetching**: once a month after the freeze, the qualifying sequences whose
  first event is at least 4 weeks old, so that the catalogue's review is in,
  are identified from the network's catalogue (GeoNet's; EMSC's for KO) and
  fetched with the frozen tool, each with its catalogue excerpt at that time
  and the canary check. A sequence is fetched once; later catalogue revisions
  are not taken up.
- **Quiet data**: after the span closes, 24 windows of 6 h drawn over it by
  the rule of section 3 with seed 44 (`heldout select nz --prospective`,
  `heldout select ko --prospective`), on the station sets of the network's
  sequences, and fetched then.
- **Detectors, settings, scoring and endpoints**: those of sections 5 to 7 for
  the network, with the frozen release and settings; nothing is retuned. The
  primary endpoint is that of section 7 on the prospective sequences.
- **Not live**: the archive delivers the post-freeze data intact, so the test
  replays it like the retrospective ones. A live SeedLink run, if one is made,
  is reported apart and is not part of this test.
