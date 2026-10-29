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
- **Analysis code**: `tools/sweep_trigger.py`, `tools/event_outcomes.py`,
  `tools/significance.py` and the fetch tool of the held-out repository, at the
  frozen commits.

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
  - **Prospective**: New Zealand again, on data recorded after the freeze
    (section 10), a result of its own.
- **Stations**: three HH channels at 100 Hz, open over the pool period, with at
  least 90% data coverage over it (from the availability services, before any
  waveform is fetched). Coordinates and gains from each network's StationXML
  (`--stations`, `--pd-gains`); Pd station terms are zero for Greece and New
  Zealand.

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
- **Sequences**: every catalogue event of M >= 4.5 in the network's catalogue
  within the pool period, with at least five usable
  stations within 150 km, that is not within 3 h after an earlier qualifying
  event (which belongs to that event's sequence). Each sequence is a replay of
  3 h starting 10 min before its first event. Catalogue events
  are scored within 250 km of the station centroid, at stations where they are
  visible (SNR >= 3, as on the development data).
- **Reserve split**: the qualifying sequences of each network, in time order,
  are assigned to "now" or "reserve" by a coin with seed *[seed]*. Reserve
  sequences are not fetched until a later protocol, committed in the same way,
  names them.
- **Stopping rule**: "now" sequences are replayed in time order until the
  scorable catalogue events reach **240 on KO and 150 on Greece and on New
  Zealand**; the sequence that crosses the target is included whole. If a
  network's "now" sequences do not reach it, all of them are used and the
  result is reported as under the target.
- **Quiet data**: 24 windows of 6 h per network, start times
  drawn uniformly over the pool period with seed *[seed]*, on the station
  clusters of the "now" sequences, each accepted only if no catalogue event
  lies within 250 km of any of its stations during the window or 30 min on
  either side. Where fewer than 24 windows qualify (New Zealand), catalogue
  events below M2 farther than 100 km are disregarded. The first 30 min of each
  window warm the filters and are not scored.

## 4. Per-station STA/LTA (secondary baseline)

Each test station's own STA/LTA trigger ratio, from its own data before the
pool: four 6 h windows per network group, drawn with seed 20221231
(`tools/make_tuning_sets.py` of the held-out repository; Greece and New Zealand
in 2022, KO from 2026-07-01 to 2026-08-29), with every station of the group
that has three 100 Hz HH channels and 90% coverage. The setting is the
network-wide equal-unmatched one (section 6) with the trigger ratio replaced by
the lowest ratio of the grid (3, 4, 5, 6, 8, 10) at which the station's trigger
rate on its quiet data over the four windows does not exceed the median
per-station rate of the network-wide setting on the AFAD quiet windows. A
station above the target at every ratio gets 10; a test station without tuning
data keeps the network-wide ratio.

Quiet data, on the tuning windows and the AFAD quiet windows alike: the time
outside [D/8 - 10 s, D/6 + 120 s] after the origin of every catalogued event
visible at the station (epicentral distance D in km; visible within 50 km at
any magnitude, 150 km from M2, 400 km from M3 and 1500 km from M4.5, as in
`tools/compare_detectors.py`), after the first 30 min of each tuning window
(the quiet windows' own lead-in is not scored). A rate is the station's
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


| detector | setting |
|---|---|
| learned (6 s) | 0.8 x 8 windows, or 0.9 once |
| STA/LTA, equal unmatched alarms | *[from the retuning]* |
| STA/LTA, equal detections | *[from the retuning]* |
| STA/LTA, equal quiet alarms | *[from the retuning]* |
| FilterPicker, the three rules | *[from the retuning]* |
| PhaseNet, EQTransformer, causal; and look-ahead (offline) | *[weights and settings from the retuning]* |
| per-station STA/LTA | section 4 |

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
citations. Every run's command, commit and output is kept.

## 10. Prospective test (New Zealand)

A test on data that did not exist when the method was frozen, reported as its
own confirmatory result and never pooled with the retrospective New Zealand
test or any other network.

- **Data**: NZ (GeoNet), recorded from the freeze (the UTC time of this
  protocol's public commit) on. Stations by the rule of section 2, their
  availability checked over the span before any waveform is fetched.
- **Sequences**: by the rule of section 3 (M >= 4.5 in GeoNet's catalogue, at
  least five usable stations within 150 km, not within 3 h after an earlier
  qualifying event; 3 h from 10 min before), every qualifying sequence in time
  order, with no reserve split.
- **Span and stopping rule**: from the freeze until the scorable catalogue
  events of its sequences reach 150, the sequence that crosses the target
  included whole, or 12 months, whichever comes first; if 12 months pass under
  the target, every sequence is used and the result is reported as under it.
- **Fetching**: once a month after the freeze, the qualifying sequences whose
  first event is at least 4 weeks old, so that GeoNet's review is in, are
  identified from its catalogue and fetched with the frozen tool, each with its
  catalogue excerpt at that time and the canary check. A sequence is fetched
  once; later catalogue revisions are not taken up.
- **Quiet data**: after the span closes, 24 windows of 6 h drawn uniformly over
  it with seed *[seed]* on the station clusters of its sequences, accepted by
  the rule of section 3 (with its New Zealand fallback), and fetched then.
- **Detectors, settings, scoring and endpoints**: those of sections 5 to 7 for
  New Zealand, with the frozen release and settings; nothing is retuned. The
  primary endpoint is that of section 7 on the prospective sequences.
- **Not live**: the archive delivers the post-freeze data intact, so the test
  replays it like the retrospective ones. A live SeedLink run, if one is made,
  is reported apart and is not part of this test.
