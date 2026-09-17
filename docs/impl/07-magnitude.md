# 07 · Magnitude

Each triggered station sizes the event from a 10 s window, and the network
reports the median. The regressor is cnn_earthquake's waveform-only magnitude
model, three partitions averaged. It is transcribed like the detector and held
to PyTorch the same way.

## The model

`trained_model_magreg_fdsn10s_nonaux_p{0,1,2}`: `DualChannelNet` with
`channels="1d+2d"`, 245k parameters each.

```
waveform (1000, 3)  -> BiLSTM 3 -> 2x64 -> self-attention (4 heads) + residual -> LayerNorm -> mean -> Linear 128  x w1 ┐
spectrogram (3, 65, 32) -> conv 3->32 · conv 32->64 s2 · conv 64->128 s2 (BN + GELU each) -> mean -> Linear 128 x w2 ┘
   sum -> LayerNorm -> Linear 128 -> GELU -> Linear 1 = magnitude
```

**Why these checkpoints.**

- Training data: 13,150 FDSN windows from 161 KOERI stations. Each window
  starts 2 s before a TauP-predicted P and runs 10 s.
- Partitions: the three are station-disjoint *and* event-disjoint. Reported
  MAE is **0.42 ± 0.02** on stations the model never saw, against 0.61 for a
  ridge fit on log SNR (`docs/experiment_magnitude_fdsn10s_2026-09-06.md`).
- Inputs: waveform only, no distance. A fresh detection has no hypocentre, so
  a model that needs one cannot run in a cascade.
- Rejected alternatives: the `alarm10s` and `cont10s` models were trained on
  MANT alone and make no claim about other stations.

Averaging the three partitions is an ensemble of models that each saw a
different subset of stations.

## Station noise normalisation

This is the part the detector did not need. The regressor's inputs are scaled
by the *station's* noise:

| input | training (seismic_cli) | ayzek |
|---|---|---|
| waveform | `(cleaned - mu) / sigma`, from the station's noise files | the same, from `NoiseBaseline` |
| spectrogram | dB minus the station's median noise dB per frequency bin | the same, from `NoiseBaseline` |
| no baseline | per-window z-score, and spectrogram `(s - mean) / (std + 1e-6)` | the same fallback |

**`NoiseBaseline` is built online** from the station's own data:

- **Selection.** Every 30 s, the 10 s ending at the newest scored window
  becomes a noise window, if no detector window overlapping it scored 0.3 or
  more. It is cleaned exactly as a training noise trace is.
- **σ** pools the sums over the last 30 noise windows (5 minutes).
- **Profile.** Each window's spectrogram is computed one component at a time,
  because the training generator applied `top_db` per trace. The profile is
  the median over all stored frames, using torch's lower-median convention.
- **Readiness.** It is trusted after 6 windows (60 s), matching the generator's
  `min_seconds`.

**The fallback costs accuracy.** Without a baseline, per-window scaling throws
away absolute amplitude. On the M4.9 with only a minute of warm-up, the
estimate was M4.0; with a full baseline it was M4.6. So with `--from`, stations
still score the data before that time silently, to warm their baselines.

## Two estimates per station

| | window | runs for | computed in | why |
|---|---|---|---|---|
| early | trigger window + 1.5 s (P assumed 3.5 s in, minus 2 s) | every station trigger | station thread, trigger + 11.5 s | goes out with the alert |
| at pick | picked P - 2 s | triggers that belong to a declared event | main thread (network stage), when the 60 s pick completes | the window the model was trained on |

The at-pick estimate replaces the early one for that station. The network
magnitude is the median across stations, reprinted whenever it changes.

**Why the two are gated differently.** A station trigger need not become an
alarm: it can remain a single-station detection, or be an S-wave or coda
re-trigger that the network assigns to an earlier event. The early estimate has
to run before the network has decided, because the event is declared only when
a second station detects it, and waiting would delay the first magnitude by
about 10 s. The at-pick estimate comes about a minute after the trigger, when
the decision is known, so it runs only for triggers that belong to a declared
event.

**Mechanism.** When a pick is confident and its P lies near the trigger, the
station thread copies the 10 s window starting 2 s before P and the station's
current noise baseline into the `Pick` message (`MagnitudeWindow`,
`src/pipeline/processor.hpp`). The main thread releases messages in stream-time
order; for a pick with a window it asks the network whether the trigger belongs
to a declared event and, if so, runs the regressor and passes the estimate on
with the pick's time. The estimate is the same as one made in the station
thread, and the decision does not depend on thread timing. The regressor
runtime moves from the station thread to the main thread.

With `--record`, every pick with a window is estimated and written, and the
network ignores at-pick estimates for events not declared at that time. A
subset replay (`tools/network_subsets`) can then declare events that the full
network did not and still has their at-pick estimates.

**Effect.** Estimates computed with and without the gate:

| replay | triggers (early estimates) | picks with a window | at-pick estimates made | skipped | all estimates |
|---|---:|---:|---:|---:|---|
| Sındırgı 30 min | 33 | 23 | 20 | 3 | 56 → 53 |
| Marmara 3 h | 389 | 256 | 202 | 54 | 645 → 591 |
| Sındırgı Mw 6.1 3 h | 188 | 125 | 96 | 29 | 313 → 284 |

The saving is 5–9% of all estimates, much less than expected when this was
proposed. Most picks with a confident P near the trigger belong to declared
events, and the early estimate, which runs for every trigger, accounts for
two thirds of the work. Gating the early estimate as well would save more but
would delay the first magnitude until the second station detects the event.

The event magnitudes do not change: the Marmara replay with `--record`, which
estimates every pick, gives a report identical to the one without.

## The spectrogram

`src/spectrogram.cpp` matches torchaudio's
`Spectrogram(n_fft=128, hop_length=32, power=2)` and
`AmplitudeToDB(stype="power", top_db=80)`:

- periodic Hann window, `center=True` with reflect padding (edge sample not
  repeated), unnormalised one-sided power
- `10 log10(max(p, 1e-10))`, then a floor at the maximum minus 80 dB. The
  maximum is taken over everything passed together: all three components for
  an event window, one component for a noise trace.

The STFT is a small radix-2 FFT in double. torchaudio runs in float32, so the
two agree to about 2e-3 dB.

## Agreement

`tests/test_magnitude.cpp` against `tools/export_magnitude.py`. That script
runs in the data_downloader environment and calls seismic_cli's own
`clean_and_filter_1d`, `standardize`, `SpectrogramEncoder.spec_db` and
`normalize_spec`, the functions that built the training corpus.

| stage | max abs error |
|---|---:|
| noise σ | 8e-15 relative |
| noise dB profile | 1.3e-3 dB |
| spectrogram | 1.8e-3 dB |
| normalised waveform | 1.5e-5 |
| BiLSTM over 1000 steps | 1.9e-5 |
| 1D / 2D branch embeddings | 4.3e-6 / 4.8e-7 |
| magnitude | 9.5e-7 |

On the real DEMI window of the M4.9, the C++ ensemble gives 4.429, the same as
PyTorch to three decimals.

## Results on the demo

DEMI, MANT and BAND, the full 30 minutes, so baselines are warm:

| AFAD | at the alert (+17–18 s) | final (+63 s) | error |
|---|---:|---:|---:|
| **Mw 4.9** | **M4.3** (DEMI) | **M4.6** (3 stations) | **−0.3** |
| ML 3.6 | M3.5 | M3.6 | 0.0 |
| ML 2.7 | M2.9 | M2.9 | +0.2 |
| ML 1.9 | M2.5 | M2.5 | +0.6 |
| ML 1.4 | M2.4 | M2.3 | +0.9 |

Final errors average 0.40 over all five events, in line with the 0.42 reported
on unseen stations. The errors are not uniform, and the training corpus
explains why:

| corpus magnitude band | 2.0 | 2.5 | 3.0 | 3.5 | 4.0 | 4.5 | 5.0 | 5.5 | ≥ 6.0 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| windows | 2,229 | 2,180 | 2,142 | 2,111 | 2,061 | 1,794 | 473 | 116 | 44 |

- **Below M2.0 the model has seen nothing.** The ML 1.4 and 1.9 are outside
  its training range, and it answers M2.3–2.5, near the bottom of what it
  knows. Inside the range (ML 2.7, 3.6, Mw 4.9) the mean error is **0.17**.
- **Large events are underestimated.** Only 160 windows are M5.5 or above, and
  a 10 s window saturates. The model scores a corpus window labelled M6.2 at
  4.3. Expect an M6 to be sized low.

With all seven stations, stations 200–300 km away pull the M4.9's median down
to M4.4, so distance matters even for a model given none.

**Caveat worth saying out loud:** the model was trained on KOERI stations and
runs here on AFAD stations. Station normalisation is designed to cancel
instrument gain, and these five events suggest the transfer holds. Five events
are not a validation.

## Cost

About 100 ms per model and 300 ms per three-model estimate on x86 (scalar
build). Most of it is attention over 1000 steps: four heads of 1000 × 1000
scores. It runs on the station's processing thread, so an estimate pauses that
station's detector for that long. Replay backpressure absorbs this. A live
deployment on a Pi should move magnitude to its own worker, which is the next
change.
