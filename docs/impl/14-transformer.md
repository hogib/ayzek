# 14 · The streaming transformer detector

`--detector transformer` replaces the 3-seed 6 s window detector with the
onset project's streaming transformer. Everything after the trigger is the
same code: the picker, association, location, magnitude and catalogue scoring.
A model with the geometry head also locates events itself, replacing the
picker (`15-geometry-location.md`).
So, as with STA/LTA (`09-sta-lta.md`), a comparison measures the detection
stage alone.

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
4. **Trigger.** The first token with p ≥ threshold after p last fell below the
   release level.
   - P is dated at the token's end minus dt.
   - The Detection carries a window start 3.5 s before that P, the network
     stage's convention. The picker and early magnitude windows are therefore
     placed from the transformer's own P, with no STA/LTA anchor.
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

## Results

See the end of this file (filled from `compare_detectors.py`).
