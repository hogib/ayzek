# 05 · The real-time pipeline

```
               per station, 2 threads                                   once
 ┌──────────────────────────────────────────────────────────────┐
 │ ingest     miniSEED record -> Steim2 decode -> Reorderer(0)  │
 │               -> SpscRing<int32> x 3 (Z, N, E), gaps marked  │
 │ processor  6 s window every 0.5 s -> condition -> 3-seed     │──Bus──> network
 │               detector -> trigger -> 60 s picker window      │         associate · declare
 └──────────────────────────────────────────────────────────────┘         locate · score
```

| file | stage |
|---|---|
| `src/pipeline/common.hpp` | stream clock, messages, bus, logger |
| `src/pipeline/station.hpp` | three component rings on one absolute sample clock |
| `src/pipeline/ingest.cpp` | replay source, decode, reorder, ring writes |
| `src/pipeline/processor.cpp` | windows, detector, trigger, picker |
| `src/pipeline/network.cpp` | association, declaration, location, catalogue scoring |
| `app/ayzek.cpp` | CLI, threads, stream-time ordering, summary |

## Ingest

- **Source.** `ReplaySource` reads a miniSEED file and orders records by the
  time of their last sample, the earliest a datalogger could have sent them.
  Each record is released when the stream clock reaches that time. A SeedLink
  client would replace this class and nothing downstream changes.
- **Decoding.** `src/mseed.hpp` decodes Steim2 and uncompressed 16- and 32-bit
  integer records, the latter in either byte order. During the 2025-04-23
  Mw 6.2 Marmara earthquake the dataloggers at CATL (57 km) and ARNA (68 km)
  switched from Steim2 to 32-bit integer records; ELBA (38 km) did not. The
  first version accepted only Steim2, so CATL and ARNA lost every record of
  the strongest shaking and produced no detection. Records that cannot be
  decoded are now reported with the reason. CATL (little-endian integer
  records) and ARNA (big-endian) decode bit-exactly against ObsPy.
  `tools/scan_encodings.py` over the whole archive (13 stations, 286 chunks,
  398 million records, 2024–2026) found 304 non-Steim2 records, all at ARNA
  and CATL on 2025-04-23 between 09:49 and 12:13 UTC.
- **Reordering.** Each component has its own `Reorderer` at
  `max_lateness = 0` (`DESIGN.md`): a hole is a gap at once and a late record
  is stale.
- **Gaps are written into the ring** as `INT32_MIN` samples, so ring index and
  absolute position stay in lockstep and a window over a gap is recognisable by
  content.
- **Replay backpressure.** A full ring makes ingest wait for the processor
  rather than drop, which is right for replay and would be wrong live. The wait
  is capped at a minute of wall time. Beyond that it drops, counts and warns.
  The cap exists because a component silent for longer than the ring holds
  would otherwise stall its siblings forever.

## Processor

- **Windows.** A window is 600 samples, one every `--step` samples (default
  50, i.e. 0.5 s), on an absolute grid so every run sees the same windows.
- **Conditioning and scoring.** Each window is conditioned per component
  (`03-dsp.md`) and scored by the 3-seed ensemble (`04-models.md`).
- **Trigger.**
  - It fires when 8 consecutive windows reach p ≥ 0.8 (`--threshold`,
    `--trigger-windows`), or earlier if a window of that run reaches p ≥ 0.9
    (`--instant-threshold`). Windows within 15 s of the previous trigger do
    not start a run. The reasons for this rule are under *Threshold* below.
  - The detection carries a P time: the STA/LTA onset in that window if there
    is one, otherwise 3.5 s after the start of the first window of the run
    (`09-sta-lta.md`). Association and catalogue matching use it.
  - It re-arms after two consecutive windows below 0.3.
  - A trigger sends a `Detection` and queues a picker window starting 2 s
    before the first window of the run, and an early magnitude window.
- **Picker.** Every trigger gets a pick. A queued job runs when the scored
  windows reach the end of its data, which for the picker is about a minute
  after P. The picks serve location, not the alert. Jobs are run by the
  position of the scored windows rather than by how much data ingest has
  written, because in a fast replay ingest can be far ahead, and the result
  would then depend on thread timing (see *Two mistakes worth recording*).
- **Magnitude.** An early estimate runs 11.5 s after every trigger. A
  confident pick carries the 10 s window at its P, and the main thread
  estimates the magnitude from it if the trigger belongs to a declared event.
  Quiet windows keep the station's noise baseline current (`07-magnitude.md`).
- **STA/LTA.** One runs alongside the detector and gives the P time of a
  trigger, in place of the convention that P lies 3.5 s into the window; it
  never triggers anything. With `--detector stalta` it replaces the detector
  instead, for comparison (`09-sta-lta.md`).
- **Ring floor.** Every 256 windows the processor raises the floor to whatever
  no future window or queued job can need. That floor is what lets ingest
  write. (Raising it only after scoring everything buffered, up to 22 minutes
  of data, blocked ingest past its one-minute limit whenever a processor ran
  slower than about 20x real time, for example under a parallel sweep.)
- **Progress.** Every processor also sends `Progress`: a promise that nothing
  it sends later is declared before a given stream time.

## Network

- **Ordering.** `app/ayzek.cpp` holds messages in a priority queue and releases
  them in stream-time order against the minimum of all stations' promises.
  Without this, at full replay speed stations race each other: one station's
  detection can reach the network after later detections from other stations,
  which changes association. Messages are released only for times strictly
  before every station's promise, since a station may still send a message at
  exactly the promised time. Messages with the same time, which are common
  because detections end on the 0.5 s window grid, are ordered by station code
  and then by the order the station sent them.
- **Association.** Two detections can belong to one event only if their
  window times differ by at most the P travel time between the two stations
  plus 3 s.
- **Declaration.** An event is declared when `--min-stations` (default 2)
  agree.
- **Coda.** Later detections at a station within 40 s of a declared event's
  first are absorbed as its S wave and coda.
- **Pick quality control.** A pick is used only if its P and S probabilities
  are both at least 0.5, S follows P by 0–60 s, and P lies near the triggering
  window. The picker always returns an argmax, even in noise, so this matters.
- **Location.**
  - Epicentre grid search (±3° at 0.05°, then ±0.1° at 0.005°), depth fixed at
    10 km, uniform half-space with Vp 6.0 and Vs 3.5 km/s, using P and S from
    every accepted pick.
  - The origin time at each node is the mean residual.
  - While rms exceeds 2 s, the worst-fitting station is dropped and the event
    relocated. Stations beyond ~150 km are the usual casualties, because a
    half-space ignores the faster upper-mantle P path.
- **Report.** At the end, `Network::report` prints one block per alarm
  (magnitude, alarm time, location, AFAD comparison, S-wave warning time per
  station and `--site`) and a summary of correct and false alarms and detected
  and missed catalogue events; see `IMPLEMENTATION.md`, "Output".
- **Scoring.** `--catalog` takes an AFAD CSV export. A declared event matches a
  catalogue event if its location agrees (≤ 30 km, origin ≤ 5 s). Otherwise it
  matches if a detecting window fired where that event's P should have arrived
  (within 6 s). The summary lists every catalogue event within 250 km with its
  alert delay and location error.

## Demo: the 2025-11-10 Sındırgı M4.9

```bash
python3 tools/make_demo_data.py --out data/demo \
    --start 2025-11-10T18:05:00 --end 2025-11-10T18:35:00 \
    ~/Projects/sismokaos/tdvms/afad_raw/{DEMI,MANT,BAND,KAND,KIRK}/*_2025-10-29.zip \
    ~/Projects/sismokaos/tdvms/afad_raw/{ELBA,SEMS}/*_2025-10-21.zip
uv run --project ~/Projects/sismokaos/archive_pipeline python tools/export_models.py

CAT=~/Projects/sismokaos/data_downloader/catalogs/catalog_afad_full_2026-08-30.csv
build-release/app/ayzek --speed 10 --from 2025-11-10T18:20:00 --catalog $CAT \
    data/demo/DEMI.mseed data/demo/MANT.mseed data/demo/BAND.mseed
```

Or `tools/demo.sh`. The chunks are named by their *start* date, which is why
10 November is in the `2025-10-29` files.

### What it found

30 minutes of the Sındırgı aftershock sequence. DEMI, MANT and BAND are 54, 91
and 124 km from the M4.9. Default trigger settings, full-speed replay:

| AFAD catalogue | ayzek |
|---|---|
| 18:07:17 ML 1.4 | alert +17.0 s, M2.3, located 6.1 km off, origin +1.9 s |
| 18:10:29 ML 1.0 | missed |
| 18:16:19 ML 1.9 | alert +17.5 s, M2.3, located 3.7 km off, origin +0.9 s |
| **18:20:51 Mw 4.9** | **alert +18.0 s with M4.3, final M4.6, located 5.3 km off, origin +1.0 s** |
| 18:23:21 ML 3.1 | alert +27.0 s (BAND, MANT), M3.6; two-station location 199 km off, origin −21 s |
| 18:25:51 ML 3.6 | alert +18.0 s, M3.6, located 4.7 km off, origin +1.2 s |
| 18:29:11 ML 2.4 | missed |
| 18:30:02 ML 2.7 | alert +18.0 s, M2.9, located 4.4 km off, origin +1.5 s |
| 18:32:40 ML 2.4 | alert +22.5 s, M2.9; two-station location 100 km off, origin +3.4 s |

That is 7 of 9 declared. Two alarms have no catalogue event: one at
18:29:15.5, 4.5 s after the ML 2.4's origin and so before its P wave reached
DEMI, and one at 18:34:24.5. An ML 2.6 at 18:34:51 is excluded from scoring,
because its waves arrive after the replay ends. With the earlier single-window
rule at 0.9, the ML 3.1 and the 18:32:40 ML 2.4 were missed and one alarm had
no catalogue event.

Five of the seven are 3.7–6.1 km from AFAD's epicentres (rms 0.16–0.26 s),
with origins 1–2 s late, which is about what a uniform velocity model costs.
The ML 3.1 and the 18:32:40 ML 2.4 were located from two stations each (BAND
and MANT, DEMI and MANT) and are 199 and 100 km off. With P and S at only two
stations the misfit can have two minima, one on each side of the line through
the stations; this is a possible cause, not checked. The ML 1.4 was also
located from two stations and is 6.1 km off.

**Alert delay is set by the network, not the software.** The alert waits for
the second-nearest station's P wave (MANT, 91 km, ~15 s) plus the ~2 s until the
P onset is well inside a window. Compute adds milliseconds. A denser network,
or `--min-stations 1` for the nearest station alone, shortens it.

### Threshold

**The saturation plateau.** The detector was trained with label smoothing
(targets 0.1 and 0.9), and its outputs saturate between 0.90 and 0.907. The
maximum at each of eight Marmara stations over three hours lies in that range.
The original rule, a single window at p ≥ 0.9, therefore sits on the plateau:
whether a clearly recorded event triggers can depend on the third decimal. KIRK
reached 0.8992 for the 2025-04-23 Mw 6.2 and did not trigger. The same rule
detected 39 of 57 catalogue events in the Marmara replay.

**Rule.** A threshold below the plateau passes many short excursions of noise.
Requiring the probability to stay above it for several consecutive windows
removes most of these, but it delays every trigger by `(N − 1) × 0.5 s`. A
single window at the plateau still triggers at once, so events that saturate
the detector are not delayed.

**Evaluation.** The detector's probabilities do not depend on the trigger
settings. One run with `--scores DIR` stores them, and runs with `--scores-in
DIR` reuse them, which gives an identical report (checked by diffing) in
under a second for 30 minutes of three stations. `tools/sweep_trigger.py` runs a
grid of settings on the cached probabilities.

Two corrections were needed before the counts could be compared:

- **Chance matches.** An alarm matches a catalogue event by location, or when
  the event's predicted P arrival at an alarm station is within 6 s of the
  detection. With many alarms, some events match by coincidence. The chance
  level was measured by shifting all catalogue times by ±127, ±263 and ±391 s
  and repeating the scoring. With observed detection rate *o* and chance rate
  *c*, the corrected rate is (*o* − *c*) / (1 − *c*).
- **Stations that joined later.** Originally, a station that joined an event
  after the alarm could also satisfy the P criterion. In one case an alarm
  raised 4.5 s after an ML 2.4's origin time (before its P wave could reach
  either alarm station) was credited to it through BAND's later detection.
  Only the detections that raised the alarm are now used. This did not change
  the result of the original rule on any dataset.

Settings were chosen on two datasets and then checked on a third, which was
not used in the choice:

| dataset | stations | duration | catalogue events |
|---|---:|---:|---:|
| Sındırgı 2025-11-10 18:05 (tuning) | 3 | 30 min | 9 |
| Marmara 2025-04-23 08:30 (tuning) | 8 | 3 h | 57 |
| Sındırgı Mw 6.1 2025-08-10 16:00 (held out) | 6 (CMH, DEMI, KAND, KIRK, SEMS, VIZE) | 3 h | 66 |

Detected / catalogue events, chance-corrected detection rate, and alarms
without a catalogue event:

| threshold | windows | instant | Sındırgı 30 min | Marmara | held out |
|---:|---:|---:|---|---|---|
| 0.9 | 1 | – | 5/9, 56%, 1 | 39/57, 64%, 10 | 24/66, 33%, 2 |
| 0.85 | 4 | 0.9 | 6/9, 64%, 6 | 46/57, 75%, 52 | 51/66, 73%, 20 |
| 0.85 | 6 | 0.9 | 6/9, 64%, 5 | 47/57, 79%, 29 | 42/66, 58%, 13 |
| 0.85 | 8 | 0.9 | 6/9, 67%, 2 | 44/57, 74%, 16 | 32/66, 45%, 2 |
| 0.8 | 4 | 0.9 | 7/9, 73%, 9 | 47/57, 74%, 80 | 52/66, 73%, 44 |
| **0.8** | **8** | **0.9** | **7/9, 78%, 2** | **44/57, 73%, 24** | **39/66, 55%, 7** |
| 0.8 | 12 | 0.9 | 5/9, 56%, 1 | 42/57, 70%, 10 | 25/66, 35%, 2 |
| 0.7 | 12 | 0.9 | 6/9, 67%, 1 | 42/57, 69%, 15 | 25/66, 35%, 2 |

(Full grid: `data/runs/v6/model_sweep.csv`, not tracked.)

- **Choice.** 0.8 with 8 windows, plus the single-window rule at 0.9, was
  chosen on the two tuning datasets before the held-out run. (The correction
  for stations that joined later was made after it; the choice was not
  revisited, and the table uses the corrected scoring.) On the held-out
  data it detects 39 instead of 24 of 66 events, with 7 instead of 2 alarms
  without a catalogue event. Over all three, the chance-corrected count rises
  from 64 to 85 of 132 events.
- **Where the gain is.** On the held-out data, M2–3 detections go from 2 to 12
  of 27, and M3–4 from 16 to 21 of 32. Both settings detect 5 of 6 M4–5 and the
  Mw 6.1.
- **Alarm time.** The main events alarm at the same time as before (Mw 6.2
  +11.0 s, Mw 6.1 +17.0 s, M4.9 +18.0 s), because they reach the plateau. The
  median alarm delay does not change on any dataset.
- **Unmatched alarms.** Marmara goes from 10 to 24. 2 of the 24 fall in the 79
  minutes before the Mw 6.2, when the catalogue should be complete; the others
  follow the Mw 6.2 and 5 of them are located within 15 km of it
  (`08-station-selection.md`).
- **Sensitivity.** On the tuning data, settings around the chosen one give
  similar counts. The held-out data is more sensitive to the number of windows:
  at 0.8, 8 windows detect 39 events, 10 windows 29 and 12 windows 25. A
  likely reason, not checked, is that only two of its stations (DEMI 60 km,
  CMH 85 km) are within 100 km of the sequence, and small events at these
  distances stay above the threshold for fewer windows. 0.85 with 6 windows
  detects 42 there, with 13 unmatched alarms.
- **Limits.** 132 catalogue events in three aftershock-rich windows, two of
  them around the same sequence. A quiet period, to measure the false-alarm
  rate without uncatalogued aftershocks, has not been evaluated.

### Association

Two detections can belong to one event if their P times differ by no more than
the P travel time between the stations plus `--slack` (3 s). The slack was sized
for the old P estimate, which was good to about 2.5 s; anchored P times are good
to 0.1 s (`09-sta-lta.md`), so the obvious question is whether a tighter slack
removes alarms that pair by chance. It does not.

Chance-corrected detections and alarms without a catalogue event, from the
same replays (`data/runs/v7/slack_sweep.csv`, not tracked):

| slack | Sındırgı 30 min | Marmara | held out |
|---:|---|---|---|
| 3 s | 7.0, 2 unmatched | 41.7, 24 | 36.4, 7 |
| 1.5 s | 7.0, 2 | 41.7, 25 | 35.4, 7 |
| 0.5 s | 7.0, 2 | 43.0, 25 | 35.4, 7 |

Nothing moves by more than one event or one alarm. The reason is that the
slack is the small term: the tolerance is dominated by the inter-station travel
time, which is 7 s for ARNA–ELBA and about 50 s for the widest pairs in these
networks. Anchoring the P times cannot help with that.

The lever that does work is requiring a third station (`--min-stations 3`),
and it costs warning time:

| | Marmara | held out |
|---|---|---|
| detected, 2 stations | 44 of 57 | 39 of 66 |
| detected, 3 stations | 38 | 25 |
| unmatched, 2 → 3 stations | 24 → 12 | 7 → 1 |
| median alarm delay | 12.0 → 13.2 s | 17.5 → **42.0 s** |
| largest event | Mw 6.2 +11.0 → +11.5 s | Mw 6.1 +17.0 → **+37.0 s** |

In the dense Marmara network a third station halves the unmatched alarms for
half a second of delay. In the held-out network it is ruinous: the third
station is 230–290 km away, so the Mw 6.1 alarm arrives 20 s later, which is
most of the warning. The default stays at two stations; three is worth setting
only where the network is dense.

A real reduction in unmatched alarms would need the declaration to depend on
whether the detections fit one origin, rather than on pairwise time
differences. That is not implemented: the location runs after the alarm, from
the picks, about a minute later.

### Scaling

All seven stations (up to 296 km) at full speed on a 12-thread x86 laptop:

- **Detector and picker only:** 25,924 windows in 24.6 s wall. That is 528
  station-seconds per second, 75× real time per station, at 5.2–6.6 ms per
  window (three models) and 13–16 ms per 60 s pick.
- **With magnitude:** 27.8 s wall, 67× real time per station.

The M4.9 was detected at all seven.

## Two mistakes worth recording

**Dropped data in a parallel sweep.** An early threshold sweep ran four replays
in parallel. Oversubscribed, the processors stalled past the backpressure cap
of the time (2 s), and ingest silently dropped samples. The results changed
from run to run. One of them showed the ML 3.6 with an origin 10 s off, which
briefly read as an AFAD timing error. It was data loss. Since then drops are
counted and printed, and the cap is a minute.

**Results that depended on thread timing.** Two replays were diffed after that
fix and agreed, which was taken as evidence of determinism. Much later, three
identical Marmara replays gave three slightly different reports: the same alarm
and detection counts, but different stations listed as raising some alarms,
different picks, and so different locations. There were two causes:

- Messages with equal times were released in the order they reached the main
  thread, and messages were released at, not only before, the minimum promise.
- A station kept a single pending picker window. Whether a trigger within 60 s
  of the previous one got a pick depended on whether that pick had already run,
  and picks ran as soon as ingest had written their data, which in a fast
  replay depends on how far ingest is ahead of the detector.

Both are fixed (*Network* and *Processor* above). Checked by five cached
Marmara replays, three of them run concurrently, which now give identical
reports, and by a full-inference replay of the Sındırgı demo that equals a
replay from cached probabilities. With the second fix every trigger gets a
pick; the results in this document were recomputed with it.

## Not done yet

- **Early magnitude off the detector thread.** The at-pick estimate runs in
  the main thread, but the early estimate still pauses the station's windows
  for 300 ms. It should be its own worker (`07-magnitude.md`).
- **Live SeedLink.** `ReplaySource` is the seam.
- **Velocity model.** A 1D model with Moho and Pn would keep distant stations
  in the location instead of dropping them.
- **Late records.** Replay delivers records on time. The archive's own late
  records (`DESIGN.md`) are not reproduced.
