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
  integer records, the latter in either byte order. AFAD dataloggers switch
  from Steim2 to 32-bit integer records during strong shaking. The first
  version accepted only Steim2, so during the 2025-04-23 Mw 6.2 Marmara
  earthquake CATL (57 km) and ARNA (68 km) lost every record written during
  the strongest shaking and produced no detection. Records that cannot be
  decoded are now reported with the reason. CATL (little-endian integer
  records) and ARNA (big-endian) decode bit-exactly against ObsPy.
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
  - It fires when p ≥ `--threshold`, if no trigger fired in the last 15 s.
  - It re-arms after two consecutive windows below 0.3.
  - A trigger sends a `Detection` and schedules a picker window starting 2 s
    before the triggering window.
- **Picker.** The picker runs as soon as its 60 s of data is in. The P and S
  picks therefore arrive about a minute after P. They serve location, not the
  alert.
- **Magnitude.** An early estimate runs 11.5 s after the trigger, and another
  at the picked P. Quiet windows keep the station's noise baseline current
  (`07-magnitude.md`).
- **Ring floor.** Each loop raises the floor to whatever no future window or
  pending pick can need. That floor is what lets ingest write.
- **Progress.** Every processor also sends `Progress`: a promise that nothing
  it sends later is declared before a given stream time.

## Network

- **Ordering.** `app/ayzek.cpp` holds messages in a priority queue and releases
  them in stream-time order against the minimum of all stations' promises.
  Without this, at full replay speed stations race each other: one station's
  detection can reach the network after later detections from other stations,
  which changes association. With it, a full-speed replay and a real-time
  replay produce identical events. Two runs were diffed to check.
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
and 124 km from the M4.9. Threshold 0.9, full-speed replay:

| AFAD catalogue | ayzek |
|---|---|
| 18:07:17 ML 1.4 | alert +17.0 s, M2.3, located 6.1 km off, origin +1.9 s |
| 18:10:29 ML 1.0 | missed |
| 18:16:19 ML 1.9 | alert +17.5 s, M2.5, located 5.9 km off, origin +1.8 s |
| **18:20:51 Mw 4.9** | **alert +18.0 s with M4.3, final M4.6, located 4.9 km off, origin +1.1 s** |
| 18:23:21 ML 3.1 | missed |
| 18:25:51 ML 3.6 | alert +18.0 s, M3.6, located 4.7 km off, origin +1.2 s |
| 18:29:11 ML 2.4 | missed |
| 18:30:02 ML 2.7 | alert +18.0 s, M2.9, located 4.0 km off, origin +1.7 s |
| 18:32:40 ML 2.4 | missed |

That is 5 of 9 declared, all 5 located within 6.1 km. One declared event is
not in the catalogue. An ML 2.6 at 18:34:51 is excluded from scoring, because
its waves arrive after the replay ends.

Every declared event was *located*, not only detected: P and S picks from two
or three stations, rms 0.13–0.25 s. The locations sit 4–6 km from AFAD's, with
origins 1–2 s late, which is about what a uniform velocity model costs.

**Alert delay is set by the network, not the software.** The alert waits for
the second-nearest station's P wave (MANT, 91 km, ~15 s) plus the ~2 s until the
P onset is well inside a window. Compute adds milliseconds. A denser network,
or `--min-stations 1` for the nearest station alone, shortens it.

### Threshold

The detector was trained with label smoothing (targets 0.1 and 0.9), so its
probabilities saturate near 0.91 and useful thresholds sit near the top of the
range. On this 30 minutes, run sequentially so no data was dropped:

| threshold | catalogue events declared | located | declared, not in catalogue |
|---:|---:|---:|---:|
| 0.5 | 6 of 9 | 4 | 23 |
| 0.7 | 6 | 6 | 18 |
| 0.8 | 6 | 5 | 13 |
| 0.85 | 7 | 6 | 10 |
| **0.9** | 5 | 5 | 1 |

"Not in the catalogue" is not the same as false. This is an active aftershock
sequence, and AFAD's catalogue is incomplete for small events. It is not
checked here either way.

The missed events are all ML 1.0–3.1 inside the sequence's busiest minutes,
where detections at one station fall into another event's coda window or a
refractory period.

### Scaling

All seven stations (up to 296 km) at full speed on a 12-thread x86 laptop:

- **Detector and picker only:** 25,924 windows in 24.6 s wall. That is 528
  station-seconds per second, 75× real time per station, at 5.2–6.6 ms per
  window (three models) and 13–16 ms per 60 s pick.
- **With magnitude:** 27.8 s wall, 67× real time per station.

The M4.9 was detected at all seven.

## A mistake worth recording

An early threshold sweep ran four replays in parallel. Oversubscribed, the
processors stalled past the backpressure cap of the time (2 s), and ingest
silently dropped samples. The results changed from run to run. One of them
showed the ML 3.6 with an origin 10 s off, which briefly read as an AFAD timing
error. It was data loss. Since then drops are counted and printed, the cap is a
minute, and determinism is checked by diffing two runs.

## Not done yet

- **Magnitude off the detector thread.** A 300 ms estimate pauses that
  station's windows. It should be its own worker (`07-magnitude.md`).
- **Live SeedLink.** `ReplaySource` is the seam.
- **Velocity model.** A 1D model with Moho and Pn would keep distant stations
  in the location instead of dropping them.
- **Late records.** Replay delivers records on time. The archive's own late
  records (`DESIGN.md`) are not reproduced.
