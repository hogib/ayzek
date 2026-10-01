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
   - **dt restart** (`--dt-reset BELOW,FROM`; `--no-dt-reset` turns it
     off): p still ≥ threshold, and dt back at ≤ BELOW s for two tokens in a
     row after it had reached ≥ FROM s since the last trigger. The default
     comes with the model (`config.trigger`, written by
     `export_transformer.py --dt-reset`; 2,5 for a model without it); see
     [The restart setting per model](#the-restart-setting-per-model). In an
     aftershock sequence p stays high between events, so the next onset has
     no rising edge; a restarted dt says it is a new event, not the coda of
     the old one.

   A minimum of 5 s between triggers (by their P dates) applies to both;
   `--retrigger S` changes it. It was 15 s, as for the other detectors, but a
   restart already waits for dt to reach 5 s, so the longer gap only dropped
   aftershocks 5–15 s behind the last event: on marmara_ko, 5 s took the
   arrivals caught within 1 s from 43% to 49% and the median epicentre error
   from 10.1 to 7.8 km. The other detectors keep 15 s and their trigger
   rules are unchanged.

   **One earthquake, one alarm.** With 5 s a station can trigger twice in
   one earthquake (early, on a stray or a coda, and again on its P), and the
   network split such an earthquake into two alarms: 21 unmatched alarms on
   marmara_ko instead of 14, of which 6 were a second alarm for a catalogued
   event. So with the transformer the network stage
   (`NetworkConfig::merge_duplicates`, off with `--no-merge`):
   - lets a detection that fits several events join the declared one with
     the most stations, not the oldest, so a stray trigger cannot hold the
     rest of an earthquake's stations in an event of its own;
   - before declaring an event, tries it against the events declared in the
     last 30 s: if its P dates and the declared event's fit one source within
     1° of the stations (one detection per station, each choice tried, at
     most one left out, four stations or more, rms ≤ 1 s), it joins that
     event instead of raising a second alarm (a `MERGE` line). A declared
     event whose own four or more P dates fit one source loses at most one
     detection to this, so two real earthquakes seconds apart stay two.

   | marmara_ko | gap 15 s | gap 5 s | gap 5 s, merge |
   |---|---|---|---|
   | catalogued events alarmed | 34/37 | 34/37 | 35/37 |
   | unmatched alarms | 14 | 21 | 17 |
   | station arrivals alarmed before S | 240/272 | 245/272 | 256/280 |

   Against 15 s, the gap with the merge adds 4 unmatched alarms and removes
   1. The 4, by their waveforms
   (`data/runs/marmara_gap5_alarms`): two are earthquakes the catalogue
   does not list (10:11:17, 10:16:12), one a second real event 12 s after an
   uncatalogued one that took the catalogue match (10:26:17), and one a
   duplicate of the ML 2.2 at 10:24:10 that no fit can merge (its first
   alarm is mostly stray triggers). 10:16:37 is no longer raised.
   `tools/network_subsets --merge --log` re-runs the network stage from a
   recording with the merge and prints its log. A detection from a dt restart is flagged
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
5. **Declaring** (network stage). With the transformer an event needs three
   stations, not two (`--min-stations`; the other detectors keep 2). See
   [Stations needed to declare](#stations-needed-to-declare).

**Threshold.** `--threshold` and `--release` default to the model's
validation operating point, stored in `models/transformer.ayzw`: the lowest
threshold with at most one false trigger per hour on validation noise, and half
of it for release. The 6 s defaults (0.8, 0.3, 8 windows) do not apply.

**Noise baseline** (magnitude). Each token counts as the 6 s window ending with
it, so `maybe_add_noise` gates on the same probability rule.

## Export and agreement

```bash
uv run --project tools python tools/export_transformer.py --run ../../onset/runs/fdsn_wide --dt-reset 1,5
ninja -C build-release && meson test -C build-release test_transformer
```

The export needs only the onset run directory (`config.json`, the checkpoint,
`val_best.json`). The model code is a copy of onset's in
`tools/reference/onset/`, from onset commit c1f02fc with only the package
import changed; re-exporting `runs/fdsn_wide` with it gives the deployed
`models/transformer.ayzw` and the test fixtures bit for bit. Re-copy it when
onset's model, conditioning or filter changes. `--dt-reset` sets the restart
level stored in the model file; without it the export takes the rule the run
was trained with.

Rebuild after pulling. A binary built before a change to the model file's
contents ignores what it does not know: one built before `config.trigger` ran
the current model with the old 2 s restart.
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

### Stations needed to declare

In a dense aftershock sequence the transformer triggers in codas, both on
rising edges and on dt restarts. Two stations ringing at the same time then
satisfy a two-station rule. Measured by replaying `--record` recordings with
`tools/network_subsets` against each run's AFAD catalogue:

| stations | dt restart | marmara_ko found / false | Sındırgı M6.1, first 2 h: found / false |
|---:|---|---:|---:|
| 2 | on | 31 / 30 | 38 / 38 |
| **3** | **on** | **29 / 9** | **28 / 7** |
| 4 | on | 25 / 2 | 13 / 3 |
| 2 | off | 19 / 17 | 31 / 23 |
| 3 | off | 17 / 5 | 16 / 4 |

Found is catalogue events detected (37 on marmara_ko, 66 on Sındırgı).
False is alarms with no catalogue event of their own, including second
alarms for an event already declared.

Three stations is the default:

- It removes most of the false alarms and costs few events: 2 of 37 on
  marmara_ko, both small. The median alarm delay rises from 10.4 s to 11.3 s.
- Four stations suits the dense Marmara network but not Sındırgı, where only
  two stations are within 100 km.
- The dt restart stays on. It adds 12 events on each set for a few false
  alarms.

What the Sındırgı false alarms were (six stations, four of them 230–300 km
away; an M ≥ 2.5 every two minutes on average):

- 10 were second alarms for a catalogued event: a stray trigger took the
  station's slot in the event, and its real P then paired with a far
  station.
- About 20 were pairs of far stations triggering in the coda of an event
  within the previous 200 s.
- 4 triggered on the S of an event that was otherwise missed.
- 4 had no catalogued event within 200 s.

Per-station rules did not help. A longer retrigger gate, or a coda window
scaled by the station's magnitude, removed about two false alarms for every
event lost.

### The restart setting per model

How far dt falls at a new onset differs between models, so the restart
setting travels with the exported weights. The model before the stronger
second-onset training (onset `781a36c`) brought dt down only to 2–3 s at
most missed aftershocks; the one after it restarts fully, and at 2 s it
fires on coda too. Network replays at 3 stations (found / false):

| model, restart at dt ≤ | Sındırgı | marmara_ko | demo_ko | quiet_ko, 6 h |
|---|---:|---:|---:|---:|
| before, 2 s | 28/66 +7 | 29/37 +9 | 8/9 +4 | 1 alarm |
| before, 2.5 s | 33/66 +11 | 31/37 +14 | 8/9 +4 | 1 alarm |
| after, 1 s | 33/66 +10 | 32/37 +12 | 8/9 +5 | 1 alarm |
| after, 1.5 s | 44/66 +24 | 33/37 +13 | 8/9 +8 | 1 alarm |
| after, 2 s | 47/66 +35 | 33/37 +13 | 8/9 +8 | 1 alarm |

The current model ships with 1,5: about the old false-alarm level with a
few more events. 1.5 finds many more aftershocks for many more alarms.
Replays of the Sındırgı alarms show only a minority of the extra ones fit a
source in the aftershock zone as well as catalogued alarms do.

## Comparing the two detectors

```bash
uv run --project tools python tools/compare_detectors.py --json compare.json
uv run --project tools python tools/compare_detectors.py --only demo_ko,marmara_ko,quiet_ko
```

Each scorecard dataset is run with both detectors, all else equal except each
detector's defaults: the transformer declares from 3 stations, the 6 s
detector from 2 (`--args '--min-stations 2'` equalises). The script reports
two levels:

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

The deployed model: onset `runs/fdsn_wide`, epoch 15, threshold 0.9630,
restart level 1 s.

**Held-out stations** (`onset evaluate --split test`, 13,544 event traces,
283 h of noise, at 1.0 false trigger per station-hour):

| | within 0.5 s | within 1 s | within 2 s | ever |
|---|---:|---:|---:|---:|
| first onsets | 43.2 % | 55.5 % | 61.6 % | 89.2 % |
| second onsets in a coda (3,387) | 37.6 % | 46.5 % | 51.3 % | 53.2 % |

- Median latency 0.53 s; median error of the P dated from dt 0.40 s. At a
  second onset dt falls to a median 0.79 s, below the 1 s restart level.
- By label source: 78.8 % within 1 s on the 9,442 traces whose P is an AIC
  pick, 1.8 % on the 4,102 whose P is only a TauP prediction. Most of the
  gap is the label: a TauP P is often off by more than the 1 s window.
- By magnitude, within 1 s: M < 2 32 %, M 2–3 53 %, M 3–4 77 %, M 4–5 90 %,
  M ≥ 5 96 %.

**Continuous replays** (whole ayzek pipeline, each detector's defaults; the
transformer declares from 3 stations, the 6 s detector from 2). Found is
catalogue events detected; false is alarms with no catalogue event of their
own, second alarms for one event included:

| set | 6 s: found / false, median alarm | transformer: found / false, median alarm |
|---|---|---|
| marmara_ko (37 events) | 28 / 23, 11.5 s | 32 / 12, 10.7 s |
| demo_ko (9 events) | 8 / 2, 17.8 s | 8 / 5, 20.3 s |
| quiet_ko (6 h, none) | 1 alarm | 1 alarm |
| Sındırgı M6.1, first 2 h (66 events)¹ | 39 / 7, 17.5 s | 33 / 10, 42.6 s |

¹ Replayed from recordings of the station outputs, without picks and
magnitudes; the others are full runs.

- Marmara M<sub>w</sub> 6.2: alarm 9.9 s after origin, epicentre 2.1 km off
  (6 s: 10.5 s, 10.9 km).
- The transformer is the better detector where many stations are near the
  source (Marmara), the 6 s detector where only two are (Sındırgı: four of
  six stations 230–300 km away). There the 6 s detector can declare from the
  two near stations with 7 false alarms, the transformer with 42, because it
  also triggers in their codas; needing a third, distant station makes its
  alarms late.
- Combining the two (the 6 s detector confirming the transformer's
  triggers, or either counting as a station) did not beat the better single
  detector on both sequences: confirmation cut false alarms by about as many
  events as it lost.
