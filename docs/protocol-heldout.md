# Held-out tests: protocol (DRAFT)

Status: **draft, not frozen.** Values marked *[to confirm]* are proposals. The
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
  - **Greece**: HL (NOA) and HT (AUTH) *[networks to confirm with the survey]*;
  - **New Zealand**: NZ (GeoNet);
  - **KO** (KOERI): held out in time only; both models trained on its stations.
- **Stations**: three HH channels at 100 Hz, open over the pool period, with at
  least 90% data coverage over it (from the availability services, before any
  waveform is fetched). Coordinates and gains from each network's StationXML
  (`--stations`, `--pd-gains`); Pd station terms are zero for Greece and New
  Zealand.

## 3. Pools and the reserve

- **Pool periods** *[to confirm]*: KO, from 2026-08-30 (after the last training
  window of the transformer and the Pd calibration, 2026-08-29) to the freeze
  date; Greece and New Zealand, from 2023-01-01 to the freeze date.
- **Sequences**: every catalogue event of M >= 4.5 *[to confirm]* in the
  network's own catalogue within the pool period, with at least five usable
  stations within 150 km, that is not within 3 h after an earlier qualifying
  event (which belongs to that event's sequence). Each sequence is a replay of
  3 h starting 10 min before its first event *[to confirm]*. Catalogue events
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
- **Quiet data**: 24 windows of 6 h per network *[to confirm]*, start times
  drawn uniformly over the pool period with seed *[seed]*, on the station
  clusters of the "now" sequences, each accepted only if no catalogue event
  lies within 250 km of any of its stations during the window or 30 min on
  either side. Where fewer than 24 windows qualify (New Zealand), catalogue
  events below M2 farther than 100 km are disregarded. The first 30 min of each
  window warm the filters and are not scored.

## 4. Per-station STA/LTA (secondary baseline)

Each test station's own STA/LTA trigger ratio, from its own data before the
pool *[period to confirm, e.g. 2017]*: the network-wide equal-unmatched setting
(section 6) with the trigger ratio replaced by the lowest ratio of the grid
(3, 4, 5, 6, 8, 10) at which the station's trigger rate on its own quiet data
does not exceed the median per-station rate of the network-wide setting on the
AFAD quiet windows (`--stalta-on-stations`). These pre-pool data are fetched
before the freeze, for this tuning only.

## 5. Scoring

As on the development data: chance correction with catalogue shifts of +-127,
+-263 and +-391 s; an alarm credited to the catalogue event it fits best and an
event to the first alarm credited with it; paired comparisons on per-event
outcomes, with the exact McNemar test (on raw outcomes) and a paired bootstrap
of the difference in chance-corrected detections, 10 000 resamples within each
sequence.

## 6. Detectors and their settings

Selected on the development event replays (the four AFAD sequences and demo_ko)
by the two rules, at the learned detector's counts on those replays
(`data/runs/tune_dev/run.sh`, `tools/select_setting.py`):

| detector | setting |
|---|---|
| learned (6 s) | 0.8 x 8 windows, or 0.9 once |
| STA/LTA, equal unmatched alarms | *[from the retuning]* |
| STA/LTA, equal detections | *[from the retuning]* |
| FilterPicker, both rules | *[from the retuning]* |
| PhaseNet, EQTransformer, causal; and look-ahead (offline) | *[weights and settings from the retuning]* |
| per-station STA/LTA | section 4 |

## 7. Endpoints

**Primary, per network**: the difference in chance-corrected detections between
the learned detector and STA/LTA at equal unmatched alarms, over the network's
"now" sequences, with its 95% paired bootstrap interval.

**Secondary, per network**:
1. every other detector against the learned detector, as in the primary;
2. PhaseNet (causal) against the learned detector: equivalence by two
   one-sided tests at 5% with a margin of +-10% of the scorable events
   (`tools/significance.py --equivalence`);
3. alarms per day on the quiet windows, every detector, with exact Poisson
   intervals;
4. the neighbourhood of the learned detector's rule: the adjacent settings of
   its grid (threshold 0.75 and 0.85; 6 and 10 windows), reported beside it;
5. detection before S: the share of matched detections whose triggering window
   ended before the predicted S, and the primary comparison restricted to them
   *[matcher tolerance to check before the freeze]*;
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
published where the networks' licences allow. Every run's command, commit and
output is kept.
