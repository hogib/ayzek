# 14 · The streaming transformer detector

`--detector transformer` replaces the 3-seed 6 s window detector with the
onset project's streaming transformer. Association, magnitude and catalogue
scoring are the same code as for the 6 s detector. Three things differ: the
trigger (below), the picker, which is placed on the transformer's own P, and
location, which a model with the geometry head does itself instead of the
picker (`15-geometry-location.md`). A comparison with the 6 s detector
therefore measures the detection stage and those two choices together.

## What runs

`src/transformer.hpp`, per station, one sample at a time:

1. **Causal band-pass.** 4th-order Butterworth, 1–45 Hz, as second-order
   sections.
   - Each run of samples is filtered on its own, from steady-state initial
     conditions, after subtracting the mean of the run's first second. That
     first second is released one second late, all at once.
   - A run shorter than a second, or any gap, makes those samples missing.
2. **Scale.** Each component is divided by its noise RMS, then passed through
   asinh.
   - At the start the RMS comes from the stream's first second.
   - After every minute whose outputs all stayed below 0.3 and whose data is at
     least 95% complete, that minute becomes the station context and its RMS
     the new scale.
3. **Model.** Every 10 samples (0.1 s), one step:
   - the causal conv stem runs over the last 30 samples;
   - four layers of attention over the last 80 tokens (keys and values cached),
     with ALiBi distance bias;
   - cross-attention to the 8 context tokens;
   - output: p (an event is under way) and dt (seconds since it began).
4. **Trigger** (`src/pipeline/trigger.hpp`). Either of:
   - **rising edge**: the first token with p ≥ threshold after p last fell
     below the release level;
   - **dt restart** (`--dt-reset BELOW,FROM`, default 2,5; `--no-dt-reset`
     turns it off): p still ≥ threshold, and dt back at ≤ 2 s for two tokens
     in a row after it had reached ≥ 5 s since the last trigger. In an
     aftershock sequence p stays high between events, so the next onset has
     no rising edge; a restarted dt says it is a new event, not the coda of
     the old one.

   The 15 s minimum between triggers applies to both. The other detectors'
   trigger rules are unchanged. A detection from a dt restart is flagged
   (`Detection::restart`), and the network stage does not absorb it as coda
   of the event the station last detected: it can join or declare a new
   event, as an aftershock inside the 40 s coda window should.
   - P is dated at the token's end minus dt.
   - The Detection carries a window start 3.5 s before that P, the network
     stage's convention. The picker and early magnitude windows are therefore
     placed from the transformer's own P, with no STA/LTA anchor.
   - The picker searches for P only within 3 s of that P, and for S only
     after it and before this station's next accepted trigger
     (`Picker::Search`). In an aftershock sequence the 60 s window often holds
     a later event, and an unconstrained search picks its P or S instead:
     on marmara_ko, 39 of 201 picks had P outside the trigger window and 26
     an S not after P. `--pick-anywhere` restores the whole-window search.
   - The decision time is the newest sample read when the token completed.
   - Missing samples are fed as missing, not skipped: the model has a gap
     channel and runs through gaps. The 6 s detector instead drops every
     window that touches a gap.

**Threshold.** `--threshold` and `--release` default to the model's
validation operating point, stored in `models/transformer.ayzw`: the lowest
threshold with at most one false trigger per hour on validation noise, and half
of it for release. The 6 s defaults (0.8, 0.3, 8 windows) do not apply.

**Noise baseline** (magnitude). Each token counts as the 6 s window ending with
it, so `maybe_add_noise` gates on the same probability rule.

## Export and agreement

```bash
uv run --project tools python tools/export_transformer.py --run ../../onset/runs/fdsn_v1
meson test -C build-release test_transformer
```

The export imports the model from the onset checkout rather than a copy.
`tests/test_transformer.cpp` compares against PyTorch and scipy on DEMI around
the 2025-11-10 M4.9:

| stage | worst error |
|---|---:|
| band-pass across a gap and a too-short run | 3e-8 of the signal (the float32 reference) |
| network, with station context, 300 tokens | p 4.8e-7, dt 3.8e-6 s |
| network, null context | p 3.7e-7, dt 3.3e-6 s |
| raw counts → tokens (filter, first scale, gap), 500 tokens | p 6.6e-7, dt 1.9e-6 s |

### What the dt restart recovers

Measured with `tools/sweep_transformer_trigger.py` on the v2 model's token
streams (`--scores`, which for the transformer now includes dt), threshold
0.989:

| set | rising edges only | + dt restart 2,5 | unmatched/h | restarts |
|---|---:|---:|---:|---:|
| demo_ko | 18 / 36 | 21 / 36 | 8.71 → 8.71 | 3 |
| marmara_ko | 112 / 296 | 114 / 296 | 6.07 → 6.21 | 4 |
| quiet_ko | – | – | 0.19 → 0.19 | 0 |

Why the rest are missed, marmara_ko, with the restart on:

| | arrivals |
|---|---:|
| detected | 114 |
| p never reaches the threshold within 20 s | 66 |
| p reaches it, but the trigger is not re-armed and dt does not restart | 116 |

On those 116, dt stays in the coda: its median minimum after P is 9.4 s,
against 10 s before, so no setting of the rule sees them. The v2 model reads an
aftershock inside a coda as more coda. Its training crops hold one event
each, so the dt head has never seen a second onset whose dt should restart.
A training example with a second event spliced into the coda, labelled to
restart dt at its P, is what would teach it that. The trigger rule then turns
those restarts into detections without further change.

The saved `compare_detectors.py` runs cannot be re-scored for this: their
records hold only each trigger's Detection. Replay the transformer with
`--scores` once, then sweep from the CSVs:

```bash
build-release/app/ayzek --speed 0 --detector transformer --no-dt-reset --no-pick \
    --no-magnitude --scores cmp/scores_marmara_ko data/marmara_ko/*.mseed
uv run --project tools python tools/sweep_transformer_trigger.py --threshold 0.989013 \
    cmp/scores_marmara_ko:tests/catalogs/marmara_ko.csv
```

## Comparing the two detectors

```bash
uv run --project tools python tools/compare_detectors.py --json compare.json
uv run --project tools python tools/compare_detectors.py --only demo_ko,marmara_ko,quiet_ko
```

Each scorecard dataset is run with both detectors, all else equal. The script
reports two levels:

- **Station level** measures each detector on its own stream:
  - recall within 0.5–8 s of P, and median latency;
  - P-dating error, which is what places the picker window;
  - unmatched triggers per station-hour.

  The reference P is TauP, refined by an AIC pick on the station's own vertical
  where one is accepted.
- **Network level** is the scorecard's numbers: catalogue events detected,
  unmatched alarms, alarm delay, false alarms per day and magnitude error.

`--workdir DIR --reuse` rescores saved runs without replaying them.

A third table compares location: the locator each run used, how many
catalogue-matched events it located, the median time from alarm to first
location, the epicentre error, and the error on only the events both runs
located (15-geometry-location.md).

## Results

See the end of this file (filled from `compare_detectors.py`).
