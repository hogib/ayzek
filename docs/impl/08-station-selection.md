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
below are from the corrected run, with the default trigger rule (8 windows at
p ≥ 0.8 or one at p ≥ 0.9, `05-pipeline.md`). The earlier rule, one window at
p ≥ 0.9, is given for comparison where it differs.

### All eight stations

| | |
|---|---|
| | default rule | earlier rule |
|---|---|---|
| alarms | 68: 45 matching an AFAD event, 23 without one | 49: 38 and 11 |
| AFAD events within 250 km | 57: 45 detected, 12 missed | 38 detected, 19 missed |
| located | 34 of the 45 | 30 of the 38 |
| **Mw 6.2 alert** | **11.0 s after origin** (ARNA, ELBA), S-wave blind zone 38 km | the same |
| Mw 6.2 location | 3.1 km from AFAD's epicentre, origin +1.3 s, 6 stations | 11.8 km, origin +0.4 s |
| Mw 6.2 magnitude | M4.6 at +14 s, at most M4.9 at +31 s, final M4.8 | M4.5 at +15 s, final M4.6 |
| Mw 5.9 | alert +12.5 s, M4.4, located 21.1 km off | alert +12.5 s, M4.3, 59.6 km off |

Of the 23 alarms without an AFAD event, 5 were located within 15 km of the
Mw 6.2 epicentre (four within 7 km) and are probably uncatalogued aftershocks. Two others follow a missed catalogue event by 7 and
8.5 s (10:42:34 ML 2.8, 10:47:05 ML 3.2) and may be the same events. Only 2
of the 23 fall before the Mw 6.2. None of this was checked on the waveforms.

S-wave warning time for the Mw 6.2 (S arrival minus alarm time; picked S at
stations with a pick, otherwise predicted from the AFAD hypocentre with
Vs = 3.5 km/s; `--site` for the three cities):

| place | distance | warning |
|---|---:|---:|
| Silivri | 26 km | −3.5 s (S wave before the alarm) |
| ELBA | 38 km | +2.1 s |
| BAND | 55 km | +6.5 s |
| CATL | 57 km | +7.0 s |
| Istanbul (Fatih) | 67 km | +8.1 s |
| ARNA | 68 km | +7.9 s |
| Bursa | 103 km | +18.4 s |

Over all 45 detected events, the alarm preceded the S wave at 319 of 360
station arrivals (median warning 25 s). This count is dominated by the
distant stations.

The final magnitude estimates are 1.4 and 1.5 units low for the Mw 6.2 and
Mw 5.9. The
regressor's training data has 160 windows of M ≥ 5.5 among 13,150, and a 10 s
window saturates (`07-magnitude.md`).

### Subsets

Best subset per size, by the alert time for the Mw 6.2:

| k | stations | Mw 6.2 alert | events declared | not in catalogue | subsets of this size that declare the Mw 6.2 |
|---:|---|---:|---:|---:|---:|
| 2 | ARNA+ELBA | 11.0 s | 18 | 1 | 28 of 28 |
| 3 | ARNA+ELBA+SEMS | 11.0 s | 36 | 5 | 56 of 56 |
| 4 | ARNA+ELBA+KIRK+SEMS | 11.0 s | 41 | 9 | 70 of 70 |
| 5 | ARNA+BAND+ELBA+KIRK+SEMS | 11.0 s | 44 | 11 | 56 of 56 |
| 8 | all | 11.0 s | 45 | 23 | 1 of 1 |

With the earlier rule, 21 of the 28 pairs declared the Mw 6.2, and the best
subsets of 2, 3 and 4 stations declared 16, 27 and 35 events.

Marginal gain of adding a station to a subset without it (mean over the 120
subsets of 2–7 stations that exclude it):

| station | + events declared | + declared, not in catalogue | Mw 6.2 alert change |
|---|---:|---:|---:|
| ARNA | 3.40 | 0.48 | −7.4 s |
| ELBA | 7.36 | 2.11 | −7.4 s |
| CATL | 7.93 | 8.53 | −7.1 s |
| BAND | 7.20 | 2.54 | −6.9 s |
| SEMS | 12.13 | 8.45 | −1.5 s |
| KIRK | 4.49 | 3.17 | −0.9 s |
| DEMI | 1.55 | 1.20 | −0.4 s |
| MANT | 0.13 | 0.57 | 0 |

- The four stations within 70 km each shorten the Mw 6.2 alert by about 7 s.
- SEMS and CATL add the most declared events, and also the most declared
  events that are not in the catalogue. This sequence is active, so some of
  these may be real uncatalogued events; this was not checked.
- ARNA adds declarations with few unmatched events.
- KIRK contributes with the default rule. With the earlier rule its highest
  probability for the Mw 6.2 was 0.8992, below 0.9, and it added nothing to
  the Mw 6.2 alert.

### Geometry model against the replay

For every subset that declared the Mw 6.2, the observed alert time was
compared with the P travel time from the AFAD epicentre to the subset's
second-nearest station:

- observed alert − predicted P time: median **2.05 s**, IQR 0.6–2.8 s, which
  matches the model's `--latency` default of 2 s
- rank correlation between predicted and observed alert times: 0.71 (0.68
  with the earlier rule)

The largest residuals (5.7 s) are subsets in which KIRK, at 140 km, is the
second station. KIRK does not reach the 0.9 plateau for the Mw 6.2 and triggers
through the 8-window rule, which adds 3.5 s. With the earlier rule KIRK did
not trigger at all, and those subsets waited for a more distant station. The
effect of the saturation plateau on the trigger is described in
`05-pipeline.md`, *Threshold*.

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
