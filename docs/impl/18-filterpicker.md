# 18 · FilterPicker

FilterPicker (Lomax, Satriano & Vassallo 2012) is a picker in operational use:
Earthworm's `pick_FP`, PRESTo and Early-est run it. `--detector filterpicker`
puts it where the learned detector would be, as `--detector stalta` does for
the STA/LTA (`09-sta-lta.md`): picking, association, location, magnitude and
the catalogue scoring are the same code, so a comparison measures the
detection stage alone. Background, parameter values in the literature and the
other implementations: `docs/research/filterpicker.md`.

## Method

`src/filterpicker.hpp`, one vertical sample at a time per station. It is a
port of Lomax's `FilterPicker5.c` v5.0.3 (GPL-2.0-or-later, from
`picker_20130212.tgz`; provenance in `tests/reference/filterpicker/README.md`),
kept statement for statement so that the two round the same way.

1. the first difference of the signal goes through N + 1 bands, N the power
   of two at or above filter window / dt (10 bands for 3 s at 100 Hz, 11 for
   6 s); band k is two one-pole high-passes and a one-pole low-pass with time
   constant 2^k dt / 2π, so the bands' periods run from dt to 2^N dt
2. per band, E = y²; its mean and variance decay over the long-term window,
   and the band's characteristic function (CF) is (E − mean) / sd, clipped at
   5 × threshold1 (E clipped with it, so a large event does not hold the
   statistics up for long); the CF of a sample is the maximum over the bands
3. a pick is triggered where the CF rises past threshold1, and declared once
   the clipped CF summed from there reaches threshold2 × tUpEvent samples, at
   most tUpEvent after the trigger
4. the pick time is the last sample before the trigger at which the
   triggering band's CF rose past its uncertainty threshold (a slow average of
   its own CF, kept between 0.5 and threshold1 / 2); the uncertainty interval
   is symmetric about it and at least a twentieth of the band's period wide;
   the polarity comes from the band's derivative summed from that point
5. no pick in the first long-term window after a start, and after a pick none
   until the CF of every band is back at 2 or below

The one difference from the C function is the start. `Pick()` begins its first
difference from the mean of the first long-term window of samples of the first
packet it is given; a stream does not have them yet, so the port begins from
the first sample. The C function does exactly that when its first packet holds
one sample, which is how the test feeds it.

**Integration.** As for the STA/LTA: the processor feeds the samples of each
window not yet seen, 50 per 0.5 s step in continuous data, and any gap
restarts the picker. Every pick triggers, subject to `--retrigger` (15 s). The
detection is dated at the pick time and declared at the sample that accepted
the pick, which is tUpEvent (0.2 s) or more after the trigger; it carries a
window start 3.5 s before P, as the STA/LTA's do, so the picker and magnitude
windows are placed as for the other detectors. The value written per window
(`--scores`) is the largest clipped CF among its new samples, and a window is
taken as noise for the magnitude baseline when that is below 5
(`fp_noise_below`). STA/LTA anchoring is off.

## Parameters

| option | default | |
|---|---|---|
| `--fp-windows F,L,U` | 3,5,0.2 | filter window, long-term window, tUpEvent, seconds |
| `--fp-thresholds T1,T2` | 10,10 | trigger threshold, mean CF over tUpEvent to accept |

The defaults are Lomax et al. (2012) Table 1 at 100 Hz: 300, 500 and 20 sample
intervals. Earthworm's `pick_FP` floors dt at 0.02 s when it sets them, which
at 100 Hz doubles the windows to 6, 10 and 0.4 s. Neither has been tuned on
this network; the comparison selects a setting on the tuning replays, as for
the STA/LTA and the published models.

## Agreement with the C function

`tests/test_filterpicker.cpp` builds Lomax's files unchanged
(`tests/reference/filterpicker/`, the test only; nothing of it is in `ayzek`)
and runs both on the same data: 20 minutes of synthetic noise with seven P/S
pairs from 1.5 to 20 Hz, a spike and a step, and the 4 minutes of DEMI that
`test_stalta` uses when that fixture exists. Three settings: the defaults in
10 s packets, Earthworm's windows in packets of 37 samples, and thresholds
5/5. Measured 26 October 2026 (x86-64, GCC, `-O3`):

| data | setting | picks | differing | CF, max rel. difference |
|---|---|---|---|---|
| synthetic | defaults | 13 | 0 | 0 |
| synthetic | Earthworm windows | 12 | 0 | 0 |
| synthetic | thresholds 5/5 | 20 | 0 | 0 |
| DEMI | defaults | 1 | 0 | 0 |
| DEMI | Earthworm windows | 2 | 0 | 0 |
| DEMI | thresholds 5/5 | 14 | 0 | 0 |

Pick time, interval, polarity and band agree exactly; strength, polarity
weight and the CF of every sample agree to the last bit. The test allows
1e-5 on the CF (the C function returns it as floats) and 1e-9 on strength and
weight, for compilers that contract to FMA.

## First run

The demo half hour (7 KO stations, `tests/catalogs/demo.csv`), defaults,
magnitude on:

| detector | alarms | match AFAD | no AFAD event | AFAD events detected |
|---|---|---|---|---|
| `stalta` (8, 1.5; 1/30 s; 2-20 Hz) | 8 | 8 | 0 | 8 of 9 |
| `filterpicker` (3/5/0.2 s, 10/10) | 11 | 8 | 3 | 8 of 9 |

0.006 ms per 0.5 s step per station (scalar backend). One replay, untuned: a
smoke test, not a comparison.
