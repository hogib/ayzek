# 08 · Station selection

Which stations should a deployment use when more are available than can be
streamed or processed, and how much does each one contribute? Two methods are
used, and they answer different parts of the question.

| method | tool | inputs | answers |
|---|---|---|---|
| replay | `ayzek --record`, `tools/network_subsets` | recorded waveforms, AFAD catalogue | what each subset of real stations actually detects, how fast, and how accurately, for the events in the recording |
| geometry | `tools/station_geometry.py` | station coordinates, AFAD catalogue | expected alert time, coverage and location error for any source position in the region |

The replay is exact for the recorded events but limited to them. The geometry
model covers the whole region but relies on assumptions (uniform velocities, a
fixed detection reach and latency). The replay is used to calibrate and check
the geometry model.

## Runtime selection

At runtime ayzek does not choose stations: it processes every station it is
given, and the network stage uses every station that detects an event. The
detector costs about 5 ms per window per station on x86, so compute does not
limit the number of stations. Per event, the network stage already:

- associates only detections compatible with the inter-station P travel time
- removes the worst-fitting station from the location while rms > 2 s

## Replay evaluation

Station processing (detection, picking, magnitude estimation) does not depend
on the other stations; only the network stage combines them. `ayzek --record`
writes every station's outputs in processing order
(`src/pipeline/recording.hpp`). `tools/network_subsets` reads the recording
once and runs the network stage for every subset of stations, which gives the
same result as running ayzek on each subset, without reprocessing waveforms.

```bash
build-release/app/ayzek --speed 0 --catalog $CAT --record run.tsv data/marmara/*.mseed
build-release/tools/network_subsets run.tsv --catalog $CAT > subsets.csv
```

Check: for the full station set of the Sındırgı demo, the subset tool
reproduces ayzek's own summary (5 of 9 events declared, 5 located, 1 not in
the catalogue, median location error 4.9 km).

Per subset it reports the number of catalogue events declared and located,
declared events without a catalogue match, mean alert delay, median location
error, mean absolute magnitude error, and for the main (largest) event the
alert delay, the S-wave blind-zone radius at the alert, the first and final
magnitude, and the location error. A summary lists the best subset per size
and the marginal gain of each station: for every evaluated subset without the
station, the change when the station is added, averaged over those subsets.
(A comparison of subsets with and without a station at fixed size is not
used, because including one station excludes another.)

## Geometry model

`tools/station_geometry.py` evaluates every subset on a 0.1° grid of
epicentres around the stations (1° margin), at 10 km depth:

| metric | definition |
|---|---|
| coverage | fraction of the grid weight with ≥ 2 subset stations within `--reach` km |
| alert time | P travel time to the second-nearest station in reach + `--latency` |
| blind zone | distance travelled by the S wave at the alert time |
| location error | 95% horizontal error ellipse (semi-major axis) of a P+S location, from the linearised travel-time equations with pick errors `--sigma-p`, `--sigma-s` |
| azimuthal gap | largest gap between the azimuths of the stations in reach |
| worst case | coverage and alert time with the most damaging single station removed |

Grid weights (`--weighting`):

- `seismicity`: number of catalogue events (M ≥ `--min-mag`, since `--since`)
  per cell
- `footprint`: 1 for each cell with at least one such event
- `uniform`

Subsets are ranked by coverage, then alert time, then location error.

**Weighting changes the answer.** With event counts (M ≥ 2.5 since 2010), the
weight is dominated by aftershock sequences in the south-west (Sındırgı, Simav,
the Aegean), and the best subsets are southern stations (DEMI, MANT, GCAM).
With the footprint of M ≥ 4 events, BAND and CMH enter every best subset and
coverage extends towards the Marmara region. Neither weighting includes
population or exposure; for early warning, weighting by exposure would be the
more relevant choice.

## Results: 2025-04-23 Mw 6.2 Marmara Sea

Data: 2025-04-23 08:30–11:30 UTC from eight stations with continuous data at
that time (KURT's chunk has no data in the window). Distances from the Mw 6.2
epicentre: ELBA 38 km, BAND 55, CATL 57, ARNA 68, SEMS 128, KIRK 140, DEMI 205,
MANT 263. The window contains the Mw 6.2 (09:49:10), an Mw 5.9 (09:51:19),
Mw 4.8, 4.5 and 4.9 aftershocks around 10:01–10:03, and 52 smaller catalogue
events. Replay of the three hours for eight stations took 3 minutes.

### Decoder defect found by this replay

In the first run, CATL and ARNA produced no detection of the Mw 6.2. From the
Mw 6.2 onwards their dataloggers wrote some uncompressed 32-bit integer records
(encoding 3; CATL little-endian, ARNA big-endian) instead of Steim2, and the
decoder rejected every such record: 104 at CATL and 155 at ARNA in the replay
window, at 09:49–09:51 and 10:02–10:03. ELBA, closer to the epicentre, wrote
only Steim2. A scan of the whole archive found no integer records at any other
station or time. The decoder now supports integer records
(`05-pipeline.md`); both files decode bit-exactly against ObsPy. The results
below are from the corrected run.

### All eight stations

| | |
|---|---|
| alarms | 49: 38 matching an AFAD event, 11 without one |
| AFAD events within 250 km | 57: 38 detected, 19 missed |
| located | 30 of the 38 |
| **Mw 6.2 alert** | **11.0 s after origin** (ARNA, ELBA), S-wave blind zone 38 km |
| Mw 6.2 location | 11.8 km from AFAD's epicentre, origin +0.4 s |
| Mw 6.2 magnitude | M4.5 at +15 s, at most M4.9 at +32 s, final M4.6 |
| Mw 5.9 | alert +12.5 s, M4.3 |

S-wave warning time for the Mw 6.2 (S arrival minus alarm time; picked S at
stations with a pick, otherwise predicted from the AFAD hypocentre with
Vs = 3.5 km/s; `--site` for the three cities):

| place | distance | warning |
|---|---:|---:|
| Silivri | 26 km | −3.5 s (S wave before the alarm) |
| ELBA | 38 km | +2.1 s |
| BAND | 55 km | +6.6 s |
| CATL | 57 km | +7.0 s |
| Istanbul (Fatih) | 67 km | +8.1 s |
| ARNA | 68 km | +8.0 s |
| Bursa | 103 km | +18.4 s |

Over all 38 correctly detected events, the alarm preceded the S wave at 257 of
304 station arrivals (median warning 25 s). This count is dominated by the
distant stations.

The magnitude estimates are 1.3–1.7 units low for the Mw 6.2 and Mw 5.9. The
regressor's training data has 160 windows of M ≥ 5.5 among 13,150, and a 10 s
window saturates (`07-magnitude.md`).

### Subsets

Best subset per size, by the alert time for the Mw 6.2:

| k | stations | Mw 6.2 alert | events declared | not in catalogue | subsets of this size that declare the Mw 6.2 |
|---:|---|---:|---:|---:|---:|
| 2 | ARNA+ELBA | 11.0 s | 16 | 0 | 21 of 28 |
| 3 | ARNA+ELBA+SEMS | 11.0 s | 27 | 0 | 56 of 56 |
| 4 | ARNA+CATL+ELBA+SEMS | 11.0 s | 35 | 8 | 70 of 70 |
| 8 | all | 11.0 s | 38 | 11 | 1 of 1 |

Marginal gain of adding a station to a subset without it (mean over the 120
subsets of 2–7 stations that exclude it):

| station | + events declared | + declared, not in catalogue | Mw 6.2 alert change |
|---|---:|---:|---:|
| ARNA | 2.94 | 0.04 | −7.9 s |
| ELBA | 6.57 | 0.64 | −7.9 s |
| CATL | 9.07 | 4.31 | −7.6 s |
| BAND | 6.08 | 1.88 | −7.4 s |
| SEMS | 10.49 | 4.34 | −2.0 s |
| DEMI | 0.90 | 0.54 | −0.8 s |
| KIRK | 1.65 | 1.08 | 0 |
| MANT | 0.10 | 0.04 | 0 |

- The four stations within 70 km each shorten the Mw 6.2 alert by 7–8 s.
- CATL and SEMS add the most declared events, and also the most declared
  events that are not in the catalogue. This sequence is active, so some of
  these may be real uncatalogued events; this was not checked.
- ARNA and ELBA add declarations with almost no unmatched events.

### Geometry model against the replay

For every subset that declared the Mw 6.2, the observed alert time was
compared with the P travel time from the AFAD epicentre to the subset's
second-nearest station:

- observed alert − predicted P time: median **2.05 s**, IQR 0.6–2.8 s, which
  matches the model's `--latency` default of 2 s
- rank correlation between predicted and observed alert times: 0.68

The largest residuals are subsets containing KIRK. KIRK's highest detector
probability for the Mw 6.2 was 0.8992, below the 0.9 threshold, so it did not
trigger and those subsets waited for a more distant station.

**Threshold at the saturation level.** Because the detector was trained with
label smoothing, its outputs saturate between 0.90 and 0.907 (the maximum at
each of the eight stations over the three hours). A 0.9 threshold lies on this
plateau, so whether a clearly recorded event triggers can depend on the third
decimal. A threshold below the plateau, combined with a requirement for several
consecutive windows above it, would avoid this; it has not been evaluated yet.

### Geometry model over the region

Eight Marmara stations, footprint of M ≥ 4 events since 2010, reach 150 km,
latency 2.05 s (`data/runs/geometry_marmara8.png`):

| k | best subset | coverage | mean alert | mean blind zone | median location error (95%) |
|---:|---|---:|---:|---:|---:|
| 2 | DEMI+MANT | 0.31 | 19.6 s | 68 km | 9.4 km |
| 3 | BAND+DEMI+MANT | 0.36 | 19.7 s | 68 km | 9.3 km |
| 4 | BAND+DEMI+ELBA+MANT | 0.53 | 18.8 s | 65 km | 9.7 km |
| 5 | BAND+DEMI+ELBA+MANT+SEMS | 0.54 | 18.7 s | 65 km | 9.0 km |

Over the whole region the southern stations DEMI and MANT rank highly, because
the region around these eight stations includes the seismically active area
south of them. The replay ranks the northern stations highest, because its
events are in the Marmara Sea. The two results answer different questions:
the geometry model gives regional averages, the replay the performance for
one source area.
