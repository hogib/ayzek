# 15 · Location from the transformer's geometry head

With `--detector transformer` and a model trained with onset's geometry head
(`onset train --geometry 1`), ayzek locates events from the transformer's own
per-station estimates of where the event is, instead of from the S-P picker.
It is the default whenever the loaded `transformer.ayzw` contains `geo_head.*`;
`--locate picks` restores the picker, `--locate geometry` insists on the head.

## Why

The S-P locator needs a 60 s picker window after each trigger and an S pick
at two or more stations, so an event is located a minute after its alarm at
best. The geometry head reads distance and direction from the waveform the
detector already has, at every 0.1 s token after P:

| output | meaning |
|---|---|
| `log_dist`, `log_dist_sd` | epicentral distance, Gaussian in log km |
| `baz`, `kappa` | back-azimuth (station → event), von Mises |

It is trained with exactly those likelihoods, so its uncertainties are the
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
        + Σ kappa_i (1 − cos(AZ_i(x) − baz_i))
        + Σ (r_i(x) − origin(x))² / (2 σ_p²)            σ_p = 0.5 s
   ```

   on a grid (±3° at 0.05° around the mean of the single-station epicentres,
   then ±0.1° at 0.005°), depth fixed at 10 km, `r_i` the P time minus the
   uniform half-space travel time. One station with a back-azimuth already
   fixes an epicentre; without back-azimuths, two stations are needed.
3. **Quality control.** As for picks: while the P-time rms exceeds 2 s and
   more than two stations remain, the station contributing most to J is left
   out. A solution is printed when it moves by 2 km or more or its station
   count changes, with the radius of its 68% region (`+-N km`, the coarse
   grid points within ΔJ = 1.15 of the minimum).

The final report's location line says which locator ran and how long after
the alarm the first location came. `--record` writes the estimates as `G`
records, and `tools/network_subsets` locates from them when present.

## Agreement with onset

`src/pipeline/locate.cpp` transcribes `onset/src/onset/locate.py`.
`tests/test_locate.cpp` runs the scenarios of onset's `tests/test_locate.py`
and checks the same answers (for example 0.37 km and 15.78 km for an
honestly-uncertain and an overconfident bad station; 68% radii of 4.31 and
37.56 km), plus the network stage locating a declared event and leaving out
a station with a bad P time.

`tests/test_transformer.cpp` compares the head's outputs with PyTorch at every
token when the export has it (`tools/export_transformer.py` writes
`net.geo_*` and `stream.geo`). On a randomly initialised model the worst
difference was 1e-5 (log distance, relative sd and kappa, angle in radians).

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
trigger, next to the AFAD epicentre of the Sındırgı M4.9, as a first sanity
check. `--verbose` prints every station estimate as a `geo` line.

## Settings and what geo_v1 says about them

| flag | default | effect |
|---|---|---|
| `--geo-sd-scale X` | 1 | multiplies every station's distance sd |
| `--geo-max-z Z` | off | while a station's distance is more than Z (scaled) sds from the solution and more than two stations remain, the worst station is left out; a solution still failing is not reported |
| `--geo-max-err-km KM` | off | a solution whose 68% radius exceeds KM is not reported |

`tools/network_subsets` takes the same flags and re-runs the network stage
from a `--record` file in seconds, and now counts located catalogue events
within 30 km and beyond 50 km. On the geo_v1 recordings of demo_ko (6
catalogue events declared) and marmara_ko (19):

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
