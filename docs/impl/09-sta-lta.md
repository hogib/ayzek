# 09 · STA/LTA benchmark

A recursive STA/LTA trigger, the standard classical method, replaces the
learned detector with `--detector stalta`. Everything after the trigger is
unchanged: the picker, association, location, magnitude and the catalogue
scoring are the same code. The comparison therefore measures the detection
stage alone.

## Method

`src/stalta.hpp`, one sample at a time per station:

1. subtract the first sample after a (re)start, so the high-pass does not ring
   on the DC offset of the raw counts
2. 4th-order Butterworth high-pass, then 4th-order Butterworth low-pass
   (`--stalta-band LO,HI`), each as two biquads (bilinear transform with
   prewarping, as `scipy.signal.butter(..., output="sos")`)
3. energy: the square of the filtered vertical component, or the sum over all
   three components (`--stalta-3c`)
4. recursive averages over `--sta` and `--lta` seconds (Withers et al. 1998,
   as `obspy.signal.trigger.recursive_sta_lta`); the ratio is 0 for the first
   LTA length after a start

A station triggers at the first sample whose ratio reaches `--stalta-on` and
re-arms when it falls below `--stalta-off`. Any gap restarts the filter and the
averages.

**Integration.** The processor feeds the STA/LTA the samples of each window
that it has not yet seen, so within continuous data each 0.5 s step adds 50
samples. A trigger at sample *t* is declared at *t* (no window to wait for) and
dated 3.5 s earlier, the detector windows' convention for where P lies. With
that, the picker window, the early magnitude window (which starts 2 s before
*t*) and the catalogue matching are placed as for the model. The 15 s minimum
between triggers is the same. The progress promise to the network stage is the
next sample to be fed, not the end of the next window.

**Agreement.** `tests/test_stalta.cpp` compares with scipy and ObsPy on 4
minutes of DEMI around the Sındırgı M4.9 (`tools/export_stalta.py`):

| | max difference |
|---|---:|
| biquad denominators vs `butter(4, ..., output="sos")` | 2e-16 |
| filtered trace (peak 9e5 counts) | 1e-8 counts |
| STA/LTA ratio, vertical, vs `recursive_sta_lta` | 5e-14 relative |
| STA/LTA ratio, three components | 8e-14 relative |

## Evaluation

The procedure is the one used to choose the model's trigger settings
(`05-pipeline.md`, *Threshold*): tune on the Sındırgı 30 min and Marmara 3 h
replays, then evaluate once on the held-out Sındırgı Mw 6.1 3 h replay.
Detections are also given corrected for chance matches, measured with
catalogues shifted by ±127, ±263 and ±391 s.

**Grid.** 648 settings on the tuning data (`tools/sweep_trigger.py`):

| parameter | values |
|---|---|
| band | 1–10, 2–10, 2–20 Hz |
| STA | 0.5, 1 s |
| LTA | 10, 30, 60 s |
| on | 3, 4, 5, 6, 8, 10 |
| off | 1, 1.5, 2 |
| energy | vertical, three components |

The off level changes the counts by at most 2 detections and 7 unmatched alarms
for any setting on either dataset, and by at most one alarm for the two chosen
below; it was fixed at 1.5. The chance correction needs 7 runs per setting, so
it was computed for 20 candidates chosen from the uncorrected counts, which are
never lower than the corrected ones. For both rules below, every setting that
could have been selected after correction was among the candidates.

**Selection rules, fixed before the held-out run.** The model's setting gives
50.0 chance-corrected detections and 25 unmatched alarms over the two tuning
datasets.

- **Equal unmatched alarms:** the STA/LTA setting with the most
  chance-corrected detections among those with at most 25 unmatched alarms.
  Result: 2–20 Hz, STA 1 s, LTA 30 s, on 8, vertical. These are the defaults of
  `--detector stalta`.
- **Equal detections:** the setting with the fewest unmatched alarms among
  those with at least 50.0 chance-corrected detections (ties: more
  detections). Result: 2–10 Hz, STA 0.5 s, LTA 10 s, on 6, vertical.

## Results

Tuning data (66 catalogue events):

| detector | detected | chance-corrected | unmatched alarms | Marmara median alarm delay | Mw 6.2 alarm |
|---|---:|---:|---:|---:|---:|
| model: 0.8 × 8 windows, or 0.9 | 52 | 50.0 | 25 | 12.0 s | +11.0 s |
| STA/LTA, equal unmatched alarms | 38 | 33.6 | 25 | 10.6 s | +10.0 s |
| STA/LTA, equal detections | 59 | 56.3 | 63 | 10.1 s | +9.9 s |

Held out (Sındırgı 2025-08-10, 6 stations, 66 catalogue events):

| detector | detected | chance-corrected | unmatched alarms | M2–3 | M3–4 | M4–5 | median alarm delay | Mw 6.1 alarm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| model | 37 | 34.4 | 8 | 10/27 | 21/32 | 5/6 | 17.5 s | +17.0 s |
| STA/LTA, equal unmatched alarms | 25 | 21.7 | 9 | 4/27 | 14/32 | 6/6 | 17.8 s | +15.7 s |
| STA/LTA, equal detections | 42 | 37.9 | 26 | 15/27 | 20/32 | 6/6 | 17.0 s | +15.7 s |

- **At the same number of unmatched alarms** the model detects more: 50.0
  against 33.6 chance-corrected on the tuning data, and 34.4 against 21.7 on
  the held-out data (8 and 9 unmatched alarms). The difference is in M2–4.
- **At the same number of detections** STA/LTA raises about three times as many
  alarms without a catalogue event: 63 against 25 on the tuning data, 26 against
  8 on the held-out data (with 37.9 against 34.4 corrected detections).
- **STA/LTA alarms earlier for the large events**: by 1.0–1.1 s for the Mw 6.2
  and 1.3 s for the Mw 6.1. It triggers on the onset sample, while the model
  needs P inside a 6 s window. For the median event the difference is small on
  the held-out data (17.0–17.8 s against 17.5 s) and 1.4–1.9 s on Marmara.
- **Cost.** The STA/LTA update takes about 0.001 ms per station per 0.5 s
  step, against 5.2 ms for the three detector models on x86. A 3 h, 8-station
  Marmara replay without magnitude takes 2.3 s of wall time with STA/LTA
  (the rest of the pipeline: decoding, picking, association) and 3 min with
  the model.

**Limits.**

- Unmatched alarms are not necessarily false. Both detectors are scored against
  the same incomplete catalogue.
- STA/LTA had more settings to choose from (648) than the model (60), which
  favours it on the tuning data. The held-out comparison is not affected
  by this.
- Only the basic recursive STA/LTA was tested. Variants used in operational
  systems (several frequency bands, adaptive thresholds, Z-detect,
  FilterPicker) were not, and could narrow the gap.
- 132 catalogue events in total, most of them in two aftershock sequences.
