# 15 · Location from the transformer's geometry head

With `--detector transformer` and a model trained with onset's geometry head
(`onset train --geometry 1`), ayzek locates events from the transformer's own
per-station estimates of how far away the event is, instead of from the S-P
picker.
It is the default whenever the loaded `transformer.ayzw` contains `geo_head.*`;
`--locate picks` restores the picker, `--locate geometry` insists on the head.

## Why

The S-P locator needs a 60 s picker window after each trigger and an S pick
at two or more stations, so an event is located a minute after its alarm at
best. The geometry head reads distance from the waveform the detector
already has, at every 0.1 s token after P:

| output | meaning |
|---|---|
| `log_dist`, `log_dist_sd` | epicentral distance, Gaussian in log km |

It is trained with exactly that likelihood, so its uncertainties are the
weights of the location. Before S it infers distance from the P wave alone;
once S is inside the ~32 s lookback it can read the S-P time itself, and its
stated uncertainty should shrink.

## What runs

1. **Station** (`processor.cpp`). A trigger is followed by a `StationGeometry`
   message at the trigger token, then every second
   (`ProcessorConfig::geometry_every`) until 20 s after P
   (`geometry_seconds`), or until p falls below the release level. Each
   carries the detector's P time (trigger token end − dt) and the token's
   geometry.
2. **Network** (`network.cpp`, `locate.cpp`). Each event keeps the latest
   estimate per station and relocates on every new one, including the ones
   that arrived before the event was declared. The location minimises

   ```
   J(x) = Σ (log max(D_i(x), 1) − log_dist_i)² / (2 sd_i²)
        + Σ (r_i(x) − origin(x))² / (2 σ_p²)            σ_p = 0.5 s
   ```

   on a grid (±3° at 0.05° around the station centroid, then ±0.1° at
   0.005°), depth fixed at 10 km, `r_i` the P time minus the uniform
   half-space travel time. Each station is a ring around it; two rings cross
   in two mirror points, so it takes three stations (`kMinGeometryStations`).
3. **Too few stations.** Until a third station reports, the event is not
   located. Instead a `RANGE` line gives each station's distance and its
   one-sd range, and prints again when a station joins or a distance moves
   by more than 10%:

   ```
   RANGE    #1   18:07:38.30  83 km from KULA (75-92 km), 93 km from SIMA (76-113 km); 2 of 3 stations to locate
   ```

   An event that never gets a third station ends the final report with
   `not located (2 of 3 stations)` and a `distance` line in the same form.
4. **Quality control.** As for picks: while the P-time rms exceeds 2 s and
   more than three stations remain, the station contributing most to J is
   left out. A solution is printed when it moves by 2 km or more or its station
   count changes, with the radius of its 68% region (`+-N km`, the coarse
   grid points within ΔJ = 1.15 of the minimum).

The final report's location line says which locator ran and how long after
the alarm the first location came. `--record` writes the estimates as `G`
records, and `tools/network_subsets` locates from them when present.

## Agreement with onset

`src/pipeline/locate.cpp` transcribes `onset/src/onset/locate.py`.
`tests/test_locate.cpp` runs the scenarios of onset's `tests/test_locate.py`
and checks the same answers (for example 0.35 km and 15.40 km for an
honestly-uncertain and an overconfident bad station; 68% radii of 2.78 and
36.19 km with loose P times), plus the network stage locating a declared
event, leaving out a station with a bad P time, and reporting ranges until a
third station arrives.

`tests/test_transformer.cpp` compares the head's outputs with PyTorch at every
token when the export has it (`tools/export_transformer.py` writes
`net.geo_*` and `stream.geo`). On a randomly initialised model the worst
difference was 1e-5 (log distance and relative sd).

**Older models.** The head once also gave a back-azimuth (von Mises), which
let one station locate. It was dropped: the network always has three
stations, and rings with P times fix the epicentre. A model exported with the
five-output head still loads: ayzek reads its first two outputs, the same
distance, and ignores the direction, as `test_transformer` does with its
four-column fixtures. Recordings with `baz` and `kappa` in their `G` records
still load too.

## Ship a geometry model

```bash
cd ~/Projects/Codings/onset
uv run onset train --data datasets/fdsn_v1 --out runs/geo_v1 --geometry 1
uv run onset evaluate runs/geo_v1 --data datasets/fdsn_v1 --split test   # geometry table
cd ~/Projects/Codings/ayzek/ayzek_code
uv run --project tools python tools/export_transformer.py --run ../../onset/runs/geo_v1
meson test -C build-release test_transformer test_locate
build-release/app/ayzek --speed 0 --detector transformer --verbose \
    --catalog tests/catalogs/demo_ko.csv data/demo_ko/*.mseed
```

The export prints the head's estimate on DEMI 0, 5 and 10 s after the first
trigger, next to its distance from the AFAD epicentre of the Sındırgı M4.9, as a first sanity
check. `--verbose` prints every station estimate as a `geo` line.

## Settings and what geo_v1 says about them

| flag | default | effect |
|---|---|---|
| `--geo-sd-scale X` | 1 | multiplies every station's distance sd |
| `--geo-max-z Z` | off | while a station's distance is more than Z (scaled) sds from the solution and more than three stations remain, the worst station is left out; a solution still failing is not reported |
| `--geo-max-err-km KM` | off | a solution whose 68% radius exceeds KM is not reported |

`tools/network_subsets` takes the same flags and re-runs the network stage
from a `--record` file in seconds, and now counts located catalogue events
within 30 km and beyond 50 km. On the geo_v1 recordings of demo_ko (6
catalogue events declared) and marmara_ko (19), measured when the locator
still used the back-azimuth:

| settings | demo_ko: median, ≤30 / >50 km | marmara_ko: median, ≤30 / >50 km |
|---|---|---|
| defaults | 6.6 km, 4 / 2 | 11.7 km, 14 / 3 |
| sd × 2 | 3.4 km, 4 / 2 | 19.5 km, 14 / 3 |
| sd × 2, z 3 | 15.6 km, 4 / 2 | 8.7 km, 14 / 3 |
| sd × 2, z 3, 30 km | 6.8 km, 3 / 1 (5 located) | 7.6 km, 13 / 2 (17 located) |

None of them removes the far-off solutions: those are distant events on
which every station's head says about 45 km (the training set ends at 55
km), so the stations agree with each other and there is no disagreement to
gate on. The medians move in opposite directions on the two sets, which is
noise at 6 and 19 events. So the defaults stay off, and the fix is the head
itself: a model trained on wider distances.

The head's distance sd is overconfident on these sets (29–42% of true
distances within one sd instead of 68%), which is what `--geo-sd-scale` is
for once a model's calibration has been measured.

## Open questions

- σ_p = 0.5 s is a guess at the dt-dated P error; `compare_detectors.py`
  measures that error and should set it.
- The KO training set is regional: stations beyond the training distances
  get estimates the head has not learned, and their stated uncertainty is the
  only guard.
