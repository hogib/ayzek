# 17 · Magnitude from peak P-wave displacement

With the Pd model in `models/` (exported by `tools/export_pd.py`), every
declared event receives its magnitude estimate from it, whichever detector
runs; until v0.1.0-beta.6 it ran alongside the learned regressor of
`07-magnitude.md`, since retired. It is computed from
the peak displacement of the P wave (Pd) at the stations that triggered and
from their distances, and comes with the probability of each size class
(section 2.3). The estimator, its calibration and its
evaluation are described in onset's `docs/MAGNITUDE.md`. This document
describes the implementation in ayzek, its agreement with onset and its
results on the continuous replays. `--no-pd-magnitude` (or `--no-magnitude`)
disables it. Sections 4 and 5 compare it with the regressor as it ran in
v0.1.0-beta.6 and earlier.

## 1. Motivation

The regressor of `07-magnitude.md` estimates magnitude from a 10 s window of
waveform. It was trained on few large events, and it underestimates them
(onset `docs/MAGNITUDE.md`, section 1). The Pd estimator replaces learning
from examples with a scaling relation between log10 Pd, magnitude and
distance. That relation is calibrated on 96,908 values from 11,806 events,
and its extrapolation to large magnitudes rests on fitted slopes rather than
on examples of large events. It requires the absolute ground motion. The
transformer's input does not carry it, because it is normalised by the
station's noise level and band-passed at 1–45 Hz. The estimator therefore
runs on a separate processing chain on the raw counts.

## 2. Implementation

### 2.1 Station processing (`Processor`, `pd_magnitude.hpp`)

Every raw sample of the vertical component also enters a causal
displacement chain (`Displacement`): a two-pole Butterworth high-pass at
0.075 Hz, integration, and the same high-pass again. The high-pass
coefficients are exported from scipy, not recomputed. Each contiguous run
of samples starts from rest at its first sample, and a gap ends the run, as
in onset's `pd.displacement`. The processor keeps |displacement| for the last
64 s.

When a station triggers, a measurement is scheduled at the detection's onset
(`Detection::p_time`): for the transformer the P it dates (token end minus
dt); for the 6 s window detector the pick of a FilterPicker running alongside
it that lies nearest the trigger's P within 3 s, or that P (the anchoring
STA/LTA onset, or 3.5 s into the first window of the run) without one
(section 2.4); for STA/LTA its trigger sample and for FilterPicker its pick. As each window τ of the model
(1, 2, 3, 4, 5, 7 and 10 s) ends, the station sends a `PdEstimate` with

    Pd(τ) = max |d(t)| over [P, P + τ),
    noise = max |d(t)| over [P − 11 s, P − 1 s),

both divided by the sensitivity of the station's vertical channel in force
at P (`models/pd_gain.csv`; StationXML epochs, so that replays of older data
use the gain of their time). A measurement is sent only if the chain has run
without a gap from 41 s before P, the run-in and noise window of onset's
measurement. A trigger dated so far back that its noise window has left the
64 s buffer is not measured.

The transformer feeds the chain one sample at a time and runs it through
gaps in the horizontals. The window detectors (6 s detector, STA/LTA,
FilterPicker) feed it the vertical samples of each scored window not yet
seen, and a window that crosses a gap in any component is not scored, so
there the chain's run ends at a gap in any component and starts again at
the first window after it. Their values are sent at the end of the window
in which τ ends, at most 0.5 s after the transformer would send them.

### 2.2 Network stage (`Network::report_pd`)

Each event keeps, for each station, the value of the longest window
received. The estimate is recomputed whenever a value arrives, when the event
is declared, merged or relocated. The station's epicentral distance is
taken from the event's location (the S-P picks, or the transformer's
geometry) once the location's 68% radius is within 20 km
(`NetworkConfig::pd_max_err_km`); before that, for the window detectors,
from the P-time location (section 2.4) once its radius is within 20 km, and
for the transformer from the station's own geometry estimate. A poorly
constrained early location can lie hundreds of kilometres from the
epicentre, and Pd depends on distance.

The estimate (`pd_magnitude`) is the posterior mean, under a
Gutenberg–Richter prior with the exported b-value, of the censored likelihood
of the event's values under each value's own window relation, evaluated on a
grid of M from 1.0 to 8.5 in steps of 0.01. Its standard deviation is
reported as the uncertainty. A value below `min_snr` (3) times its noise
level enters as an upper bound. Values are excluded:

- from stations whose station term exceeds ±1 in log10 Pd (CTKS, DAT, TOKT
  and YEDI), taken to have erroneous response metadata;
- whose noise level lies more than a factor of 10 below the station's median
  (`models/pd_noise.csv`). Such a channel does not record ground motion, and
  its values would all enter as censored upper bounds that the event did not
  satisfy (section 4);
- whose noise level lies more than a factor of 10 above the station's median.
  Such a value was measured inside the coda of another event, and the peak
  after P includes that coda (section 4).

An event with no value above its threshold receives no estimate.

A `PDMAG` line is printed when the estimate moves by 0.1 or more or when a
station is added:

```
PDMAG    #3   18:22:00.50  Pd magnitude M4.8 ± 0.1 (M<4 0%, 4-5.5 100%, M>=5.5 0%) from 6 stations (6 above noise), windows up to 10 s, distances from the location
```

The final report gives, for each event, the first and the final estimate
(`Pd magn.`) and the error against the catalogue. The summary gives the mean
absolute error and the bias over the catalogued events that have an estimate.
`--record` writes the values as `Q` records, and `tools/network_subsets
--pd models` recomputes the estimates from a recording.

### 2.3 Size classes

With each estimate the network stage reports the posterior probability of
three size classes, small (M < 4), medium (4 ≤ M < 5.5) and large
(M ≥ 5.5), the bounds of `kPdClassBounds` (`pd_magnitude.hpp`). They are
sums of the same posterior over the same grid as the mean, a grid point at
a bound counting in the class above it. For events whose rupture outlasts
the 10 s window the estimate is a lower bound (section 5), and so is the
class.

### 2.4 Early distances for the window detectors

The S-P picks come from a 60 s picker window, so a location from them, and
with it the first Pd estimate, came about 52 s after the alarm. The window
detectors therefore also locate each declared event from its detections'
onsets alone (`locate_p_times`, `locate.hpp`): the geometry locator's grid
search with the distance term left out, from three stations on, the station
contributing most to the misfit removed while the rms exceeds 2 s and more
than three remain. This P-time location (`PLOC` lines, `P loc.` in the
report) serves the Pd distances only, until a pick location replaces it
there. Catalogue matching and the warning times keep the pick location, so
it changes no detection result.

The onsets come from a FilterPicker at its default settings running alongside
the 6 s detector on the vertical (`ProcessorConfig::onset_fp`). It triggers
nothing and does not move the detection's window start, which the
association and the catalogue matching use. A trigger takes the pick nearest
its own P within 3 s, and keeps its own P without one. Measured on the six
development replays (demo_ko, marmara_ko, demo, marmara, sindirgi61,
simav26) against TauP P refined by an AIC pick on the station's vertical,
for the 540 detected arrivals whose AIC pick was accepted:

| onset | arrivals | median error | within 0.5 s | within 1 s | robust sd |
|---|---:|---:|---:|---:|---:|
| the detection's own P, all | 540 | 0.69 s | 44% | 59% | 0.90 s |
| the detection's own P, where a pick exists | 418 | 0.51 s | 50% | 65% | 0.56 s |
| FilterPicker pick nearest it within 3 s | 418 | 0.06 s | 86% | 93% | 0.09 s |

A pick exists for 77% of the detected arrivals. Windows of 2 and 1.5 s find
fewer (375 and 305) and are no more accurate; the earliest pick in the window
instead of the nearest gives the same result.

The P-time uncertainty in the locator's misfit (`NetworkConfig::p_sigma`,
`--p-sigma`) sets its 68% radius, and so whether Pd uses the location. It is
1.35 s, the value at which the radius contains the catalogue epicentre for
69% of the 118 catalogue-matched events the P-time location locates on the
same replays (0.5 s gave 31%, 1.75 s 85%). It is much larger than the onsets'
own error because the uniform half-space at a fixed 10 km depth dominates.
Against the pick location (median epicentre error; first location after the
alarm, median):

| replay | stations | pick location | P-time location | P-time located |
|---|---:|---|---|---:|
| demo_ko | 6 | 3.7 km, 52 s | 3.3 km, 3.8 s | 8 of 8 |
| marmara_ko | 8 | 6.1 km, 52 s | 7.0 km, 1.0 s | 23 of 28 |
| marmara | 8 | 10.6 km, 52 s | 16.3 km, 2.0 s | 37 of 44 |
| simav26 | 6 | 6.1 km, 51 s | 16.1 km, 10.5 s | 23 of 37 |
| sindirgi61 | 6 | 32.2 km, 52 s | 46.2 km, 24 s | 23 of 39 |
| demo | 3 | 5.3 km, 52 s | 109 km, 6.8 s | 4 of 7 |

Within a dense network it matches the picks about 50 s earlier. Where most
stations lie far from the source (sindirgi61, four of six at 230-300 km) or
there are only three (demo), it is poor, its radius says so, and Pd waits
for the pick location. On the KO replays, the only ones with gains, the four
events that receive a Pd estimate get their first one a median 3.0 s after
the alarm (from 52 s with the pick location alone), with a mean absolute
error of 0.30 units at the first and at the final estimate.

## 3. Agreement with onset

`tools/export_pd.py` writes, from onset's `runs/pd_v6/fit.json` and the
StationXML responses:

| file | content |
|---|---|
| `models/pd_relation.csv` | per window: α, the two magnitude slopes and knot, the three distance slopes and two knots, σ, the b-value, the threshold, the fixed depth, the station-term and noise limits, and the high-pass section |
| `models/pd_terms.csv` | the station terms per window (875) |
| `models/pd_noise.csv` | the stations' median noise levels (148) |
| `models/pd_gain.csv` | the vertical channels' sensitivity epochs (154 epochs, 121 stations) |

It also writes reference outputs computed with a copy of onset's `pd.py` and
`pd_fit.py` (`tools/reference/onset/`). `tests/test_pd.cpp` holds the C++
implementation to them:

- the displacement of a 90 s synthetic record, sample by sample: the largest
  difference is 2 × 10⁻¹⁸ m, on a signal of 2 × 10⁻⁶ m;
- Pd for every window and the noise level, to a relative 10⁻⁹;
- the event estimator on five cases (three stations above the noise; one of
  three; two; two with a third that does not record; and two with a third
  inside a coda), to 10⁻⁶ in mean and standard deviation;
- the network stage, whose estimate for a declared event equals the
  estimator's at the stations' distances and ignores a value for another
  trigger.

The values that ayzek's processor measures on the Marmara replay were
recomputed independently with onset's code at the same P times, from the
same miniSEED files; Pd and the noise level agree to every printed digit.

## 4. Results

The two continuous replays used throughout (marmara_ko, the 2025-04-23
Marmara sequence; demo_ko, the 2025-11-10 Sındırgı M4.9), with the transformer
(`runs/fdsn_sp_1`), its 5 s retrigger gap and duplicate merge, and the Pd
model of onset's `runs/pd_v5` (outputs in `data/runs/pd_v2/`). The regressor
of `07-magnitude.md` runs in the same replays. Errors are against the AFAD
catalogue; "above" counts the stations whose final value exceeds three times
its noise level.

| replay | event (AFAD) | regressor | Pd, first | Pd, final | stations above |
|---|---|---:|---:|---:|---:|
| marmara_ko | ML 3.9, 09:13:05 | −0.5 | M4.0 | M3.8 ± 0.1 (−0.1) | 6 of 6 |
| marmara_ko | Mw 6.2, 09:49:10 | −0.8 | M6.4 | M6.0 ± 0.1 (−0.2) | 6 of 6 |
| marmara_ko | Mw 4.9, 10:02:32 | −0.4 | M4.9 | M5.1 ± 0.1 (+0.2) | 2 of 5 |
| marmara_ko | ML 3.4, 10:07:01 | +0.2 | M4.5 | M4.5 ± 0.2 (+1.1) | 1 of 5 |
| marmara_ko | ML 3.4, 10:16:23 | 0.0 | M3.7 | M3.8 ± 0.2 (+0.4) | 2 of 6 |
| marmara_ko | ML 2.2, 10:43:53 | +0.2 | M2.9 | M2.9 ± 0.3 (+0.7) | 1 of 4 |
| demo_ko | Mw 4.9, 18:20:51 | −0.2 | M4.5 | M4.8 ± 0.1 (−0.1) | 6 of 6 |

Seven of the 42 catalogued events that were alarmed (35 on marmara_ko and 7
on demo_ko) receive a Pd estimate, as do two uncatalogued events. The
remainder have no value above the threshold. On marmara_ko this is almost
always because the event occurred in the coda of the Mw 6.2 mainshock or of
a preceding aftershock. Of the 355 values at τ = 5 s on that replay, 17
exceed three times their noise level. After the mainshock, the noise level at
the triggering stations rises from about 10⁻⁷ m to as much as 4 × 10⁻⁴ m.

The results divide according to the number of stations above the threshold.
For the three events recorded above the noise at all six stations, the Pd
estimate is within 0.2 units of the catalogue (−0.1, −0.2, −0.1), against
−0.5, −0.8 and −0.2 for the regressor. In particular, the Mw 6.2 mainshock is
estimated at M 6.0 ± 0.1. The estimate at the alarm (M 6.4) was made from
two stations at the distances of their geometry estimates; the location was
still poorly constrained at that time. The events recorded above the noise
at one or two stations, all of them aftershocks within the coda of the
mainshock, are overestimated by 0.2 to 1.1 units. Two effects contribute.
First, the peak after P includes the still-decaying coda of the preceding
event, which the noise window before P captures only in part. Second, the
stations that did not trigger contribute no censored bound (section 5).
Over all six catalogued Marmara events, the mean absolute error is 0.44 with
a bias of +0.34, against 0.40 for the regressor over all 35 of its events.

**The channel that does not record.** In the first version of this replay,
the estimate of the Mw 6.2 mainshock was M 5.5. The vertical channel of MRMT
in these files records a noise level of about 5 × 10⁻⁹ m, 25 times below its
median of 1.3 × 10⁻⁷ m; its raw counts have a standard deviation of 6. Every
value from it is therefore censored. For the mainshock, it asserted that Pd
was below 1.6 × 10⁻⁸ m at 64 km, and that one bound lowered the posterior
from M 6.04 to M 5.50. The rule of section 2.2, which excludes values more
than a factor of 10 below the station's median noise, removes it. It was
added to onset's fit and evaluation as well, so the two implementations
remain identical. The values that ayzek measured on this replay were
recomputed with onset's code at the same P times and agree exactly. The
discrepancy was therefore a property of the data, not of the implementation.

**Values inside a coda.** The overestimated events of the table above were
aftershocks within the coda of the mainshock. Relative to each station's
median noise level, the values from which they were estimated had been
measured at 13–340 times that level (the Mw 4.9 at 10:02 at 162–344 times,
the ML 3.4 at 10:07 at 37–108 times, the ML 3.4 at 10:16 at 13–32 times),
whereas the events within 0.2 units of the catalogue were measured at about
the median. The peak after P of such a value includes the coda of the
preceding event. Requiring two stations above the threshold would not
address this, since offline, where most events have one test-split station
above the noise, single-station estimates are not biased (bias −0.12 to
+0.08 at τ = 5 s). Values more than a factor of 10 above their station's
median are therefore excluded (section 2.2; onset `runs/pd_v6`, in which the
same rule changes no offline test result by more than 0.01). The replay
outputs with the rule are in `data/runs/pd_v3/`:

| replay | event (AFAD) | Pd, final, without the rule | Pd, final, with the rule |
|---|---|---:|---:|
| marmara_ko | ML 3.9, 09:13:05 | M3.8 (−0.1) | M3.8 (−0.1) |
| marmara_ko | Mw 6.2, 09:49:10 | M6.0 (−0.2) | M6.0 (−0.2) |
| marmara_ko | Mw 4.9, 10:02:32 | M5.1 (+0.2) | no estimate |
| marmara_ko | ML 3.4, 10:07:01 | M4.5 (+1.1) | no estimate |
| marmara_ko | ML 3.4, 10:16:23 | M3.8 (+0.4) | no estimate |
| marmara_ko | ML 2.2, 10:43:53 | M2.9 (+0.7) | M2.9 (+0.7) |
| demo_ko | Mw 4.9, 18:20:51 | M4.8 (−0.1) | M4.8 (−0.1) |

Over the catalogued Marmara events, the mean absolute error falls from 0.44
(six events, bias +0.34) to 0.32 (three events, bias +0.12). The three
aftershocks within the coda now receive no estimate instead of an incorrect
one, at the cost of the Mw 4.9, whose estimate had happened to be close.
The remaining error of +0.7 is that of a small event (ML 2.2) recorded above
the noise at one station at an ordinary noise level.

**Conclusion.** For events recorded above the noise at most of their
stations, the Pd estimate is substantially better than the regressor's,
most clearly for the largest event. Within the coda of a larger event the
estimator now declines rather than overestimates. Contributing censored
bounds from stations within range that did not trigger is the evident next
step for small events.

## 5. Limitations

- Only stations that triggered contribute. In onset's evaluation every
  station of the event pull contributes, including those at which the event
  remained below the noise. The censored bounds of the stations that did not
  trigger are absent here, which biases the estimates of small events
  upward.
- Inside the coda of a larger event the noise level is the coda itself: the
  values of later events are censored or, beyond ten times the station's
  median noise, excluded. Most aftershocks of a large event therefore receive
  no Pd estimate (section 4).
- For events whose rupture outlasts the 10 s window the estimate is a lower
  bound (onset `docs/MAGNITUDE.md`, section 6.1).
- The relation and the station terms come from KOERI's network. A station
  without a gain epoch in its StationXML is not measured, and one without a
  station term is assumed to have a term of zero.
