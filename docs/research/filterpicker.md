# FilterPicker: reference note for an ayzek comparison

Research note, 2026-10-26. It answers a referee's suggestion to compare ayzek with a
standard operational picker. No code was changed. Bracketed numbers such as [1] point
to the reference list at the end. **Unverified** marks a claim I could not check
against a primary source.

Sources read in full for this note:
- the Lomax et al. (2012) and Vassallo et al. (2012) author PDFs [1, 2];
- the C distribution `picker_20130212.tgz` from Lomax's archive [3], which contains
  FilterPicker5.c v5.0.3;
- the current Earthworm `pick_FP` tree on GitLab [5];
- the Early-est mirror's picker [7];
- the `mbagagli/filterpicker` Python repository [10];
- PhasePApy's FBpicker [12].

I did not compile-and-run any of this code. Every run time below is an estimate.

**Update, 2026-10-26.** ayzek has since moved to GPL-3.0-or-later, which removes
the licence obstacle discussed below, and FilterPicker5 is ported into the pipeline
(`docs/impl/18-filterpicker.md`), checked against Lomax's C code.

---

## 1. Algorithm

FilterPicker (FP) is a single-channel phase detector and picker. It runs one sample
at a time with recursive state, so it suits continuous real-time streams [1].
According to its authors it is "loosely based" on the Baer–Kradolfer (1987) and Allen
(1978) pickers [1].

### 1.1 Filter bank

- Input is the raw broadband stream y(i) with sample interval Δt. "Little or no
  preprocessing" is needed: no filtering and no mean removal [1].
- The first previous sample is initialised to the mean of the first T_long of data
  (code: `init_filterPicker5_Memory`) [3].
- The bank has N_band bands, n = 0 … N_band−1, with corner period T_n = 2ⁿ Δt and time
  constant w_n = T_n / 2π [1].
- Each band applies two one-pole high-pass filters, with C_HP = w_n / (w_n + Δt),
  followed by one one-pole low-pass filter, with C_LP = Δt / (w_n + Δt) [1, 3].
- The paper's eq. (5), the low-pass, has an erratum. The corrected form is
  Y_n(i) = Y_n(i−1) + C_LP [Y_n^HP2(i) − Y_n(i−1)] [4].
- The paper describes the result as "similar to the raw signal twice differentiated".
  It emphasises high frequencies, so impulsive onsets keep their sharpness [1].
- In code, the HP recursion is y[i] = α(y[i−1] + x[i] − x[i−1]), applied to the raw
  sample differences [3]. The PDF renders eqs. 2–4 poorly. The description here
  follows the C code.
- Number of bands:
  - Paper: N_band = ceil(log2(T_filter / Δt)) [1, Table 1].
  - Code: it doubles a counter until it reaches T_filter / Δt, then adds one band [3].
  - At 100 Hz with T_filter = 3 s, this gives 10 bands with T_n = 0.01 … 5.12 s.
  - Cost grows linearly with N_band, so "N_band should not be set higher than
    necessary" [1].

### 1.2 Characteristic function

- Per band, the envelope is E_n = Y_n² (code: `xRec[k] = dy*dy`) [3].
- The band CF is F_n(i) = (E_n(i) − ⟨E_n⟩(i−1)) / ⟨σ_En⟩(i−1).
  - The mean and standard deviation are recursive averages with decay
    C_long = 1 − Δt / T_long [1, eq. 1; 3].
  - Table 1 of [1] prints the note as "C_long = 1 − T_long/ΔT". The text and the code
    both use 1 − Δt/T_long, so the table note looks like a typo.
- Each band CF is clipped at 5·S1. The clipped value is also written back into the
  statistics, "to avoid long recovery time after strong events" [3].
- The summary CF is F_C(i) = max_n F_n(i). The paper prefers the maximum to the sum
  because a sum "would fold noise into the summary CF" and miss narrow-band onsets [1].
- Because of the first difference and the adaptive statistics, the CF does not depend
  on absolute amplitude once T_long has passed. Offsets and long-period noise do not
  matter [1].

### 1.3 Trigger and pick declaration

Parameter names are given as paper / code.

1. **Trigger.** F_C ≥ S1 (`threshold1`). The trigger band k is the highest-frequency
   (lowest-n) band with F_n ≥ S1 [1, 3].
2. **Pick.** A pick is declared once the integral of the clipped F_C, from the trigger
   onward, reaches S2·T_up (`threshold2` × `tUpEvent`) [1, 3].
   - In the code, the integral runs over a ring buffer of n_up = round(T_up/Δt) + 1
     samples. The critical value is n_up·S2 [3].
   - A strong onset can therefore be declared before T_up has elapsed. The latest a
     pick can arrive is T_up after the trigger. "For a 100 samples/sec data stream the
     default value of T_up is 0.2 sec" [1].
3. **Spike guard.** Each sample's contribution to the integral is limited to 5·S1 [1, 3].
4. **Re-arm.** After a pick, no new pick is allowed until every band CF has dropped
   below 2.0 (`maxAllowNewPickThreshold`). This is the "avoids excessive picking
   during large events" feature [1, 3].
5. **Warm-up.** No picks are made during the first T_long after the start of data, or
   after a restart [3, 6].

### 1.4 Pick time, uncertainty, polarity and strength

- **Pick time.** Each band keeps a running "uncertainty threshold": the long-term mean
  of its CF, clamped to [0.5, S1/2]. Each time the band CF rises past it, that sample
  is stored [1, 3]. At declaration:
  - t_pick is the stored time for band k;
  - the uncertainty σ_pick = t_trig − t_pick;
  - the pick is reported as [t_pick − σ, t_trig], so t_pick ± σ_pick [1, 3].
- **Minimum uncertainty.** If the interval is narrower than T_k/20, it is widened
  symmetrically. The paper states the same rule as T_k/40 for t_trig − t_pick [1, 3].
- **Polarity.**
  - Over [t_pick, t_trig], FP sums the first differences of the band-k filtered
    signal, and also their absolute values.
  - If |Σ| / Σ|·| > 0.667, the polarity is the sign of Σ. Otherwise it is unknown [1, 3].
  - v5.0.3 (2012-10-19) adds a polarity weight, (ratio − 0.667) / (1 − 0.667) [3].
    The Earthworm copy does not have this field [5].
- **Strength and period.** The CF value at the trigger is reported as pick "amplitude"
  (strength). The trigger band's period is also reported [1, 3].
- **Phase type.** FP does not label phases. Every pick is an onset; P/S assignment is
  left to an associator.

### 1.5 Parameters

Defaults from [1, Table 1] and [2, Table 3]. The Earthworm column is what `pick_FP`
uses when a station-file value is negative [6, 5].

| Paper | Code / Earthworm name | Controls | Default [1] | Earthworm autoset [5] |
|---|---|---|---|---|
| T_filter | `filterWindow` | Longest band period, and so the number of bands (cost) and the longest onset period FP can pick | 300 Δt | 300 × max(Δt, 0.02 s) |
| T_long | `longTermWindow` | Decay constant of all running statistics (how fast the CF adapts to noise changes), and the warm-up delay | 500 Δt | 500 × max(Δt, 0.02 s) |
| S1 | `threshold1` | Trigger level on the summary CF. The CF is clipped at 5·S1 and the uncertainty threshold is capped at S1/2 | 10 | from the station file |
| S2 | `threshold2` | Mean CF over T_up needed to accept the pick | 10 | from the station file |
| T_up | `tUpEvent` | Validation window. It sets the maximum pick latency after the trigger | 20 Δt | 20 × max(Δt, 0.02 s) |

Further constants are hard-coded in the C code, not exposed as parameters [3]:
- clip at 5·S1;
- re-arm at CF < 2.0;
- uncertainty-threshold clamp to [0.5, S1/2];
- polarity ratio 0.667.

Earthworm floors Δt at 0.02 s for high-rate data [5]. Lomax's own test harness has the
same line commented out [3]. So at 100 Hz the "defaults" differ by implementation:
- Lomax table: 3.0 s / 5.0 s / 0.2 s;
- Earthworm autoset: 6.0 s / 10.0 s / 0.4 s (11 bands).

### 1.6 Computational cost

- The paper says the algorithm uses "almost exclusively logical and arithmetic
  operations (i.e., uses the square root but no exponentials, logarithms, or
  transform algorithms) and so is computationally highly efficient" [1].
- Per sample and band, the code does:
  - three one-pole filters;
  - a square and a division;
  - two recursive-average updates and one `sqrt`;
  - ring-buffer bookkeeping [3].
- Memory per channel is O(N_band × n_up) doubles [3].
- **Estimate (not measured):** at 100 Hz with 10 bands, a channel needs 1,000
  band-updates per second, which is a negligible fraction of one core. The 102.8 h
  quiet set is about 3.7 × 10⁷ samples per channel, about 3.7 × 10⁸ band-updates,
  which should take seconds per channel in C.
- **Only published timing:** Mousavi et al. (2020) ran a Python FilterPicker over the
  STEAD test set (>113k one-minute three-component traces) in 3 h 25 min on a 2.7 GHz
  i7. Kurtosis took 62 h 12 min and AIC 31 h 18 min [13]. That run used the PhasePApy
  implementation, not Lomax's code (see §4).

### 1.7 Versions and variants

- **TestPicker4_1 → FilterPicker5.** Early-est replaced TestPicker4_1 with
  FilterPicker5 in v0.1.1 (2010-01-14) [8]. Lomax's Java class is
  `net.alomax.timedom.FilterPicker5` [4].
- **FilterPicker5 C history** (from the file header [3]):

  | Version | Date | Change |
  |---|---|---|
  | 5.0.0 | 2008-07-14 | Ported from Java |
  | 5.0.1 | 2011-11-04 | Release |
  | 5.0.2 | 2011-11-14 | Release |
  | 5.0.3 | 2012-10-19 | Polarity weight added |

  Lomax's archive has four tarballs, 2011-10-24 to 2013-02-12 [3].
- **Early-est mirror.** It adds a 2022-02-02 fix (a `ppick_list_size` argument for the
  pick-list allocation) and makes `message_str` static [7].
- **"FilterPicker 6".** Not found anywhere (searched the web and GitHub code).
  **Unverified / probably does not exist.**
- **Vassallo et al. (2012).** This is the companion paper on parameter optimisation by
  genetic algorithm, not a new version [2].
- **Reimplementations** with different numerics:
  - Kamer (MATLAB);
  - Bagagli (Python, derived from the MATLAB version);
  - PhasePApy's FBpicker (a different algorithm in the same spirit).

  See §3.

## 2. Parameter settings in practice

| Setting | Sampling | T_filter | T_long | S1 | S2 | T_up | Context | Source |
|---|---|---|---|---|---|---|---|---|
| Paper default | any | 300 Δt (3.0 s at 100 Hz) | 500 Δt (5.0 s) | 10 | 10 | 20 Δt (0.2 s) | "generally work well" | [1] Table 1 |
| Paper Fig. 1–2 | 20 Hz | 15.0 s | 25.0 s | 10 | 10 | 1.0 s | synthetic; teleseismic M7.6 Izmit at G.RER | [1] |
| Paper Fig. 4 | 125 Hz | 2.4 s | 4.0 s | 10 | 10 | 0.16 s | small regional event, ISNet | [1] |
| Vassallo optimised | 125 Hz | 0.865 s (108 Δt) | 12.0 s (1500 Δt) | 9.36 | 9.21 | 0.388 s (49 Δt) | ISNet (southern Italy) local/regional events and noise; genetic algorithm on 105 traces | [1] Table 1; [2] Table 3 |
| Earthworm `pick_FP` autoset | 100 Hz (Δt floored to 0.02) | 6.0 s | 10.0 s | — | — | 0.4 s | used when the station file gives −1 | [5] `call_FilterPicker.c` |
| Earthworm sample station file | ISNet channels | −1 (auto) | −1 (auto) | 8.6 | 17.2 | −1 (auto) | example `pick_FP.sta` | [5, 6] |
| Earthworm / Lomax harness comment | — | 4.0 s | 10.0 s | 8.0 | 8.0 | 0.2 s | "defaults" comment in code | [3, 5] |
| Early-est, raw broadband | teleseismic | 100 s | 50 s | 10 | 10 | 2.0 s | global tsunami-warning monitor | [7] source defaults |
| Early-est, HF stream | teleseismic | 10 s | 25 s | 10 | 8 | 2.0 s | same | [7] |
| Python port defaults | any | 3.0 s | 5.0 s | 10 | 10 | 0.2 s | `filterpicker` 1.1.0 | [10] |
| PRESTo in Liaoning, China | ? | ? | ? | 10 | 8 | 0.5 s | regional EEW | [17] **unverified** (search snippet only; full text not accessible) |
| PRESToPlus, schools study (ISNet / SOSEWIN) | 100 Hz? | 1 s? | 55 s? | 15 / 11 | ? | ? | EEW feasibility, southern Italy | [18] **unverified** (search snippet only) |

Notes:
- **Turkey.** I found no published FP parameter values for AFAD or KOERI stations.
  - Clinton et al. (2016) report that PRESTo, which uses FP, was installed at KOERI
    within the REAKT project, alongside the Istanbul EEW system. They say it was not
    integrated with the existing system [19]. **Unverified**: seen only in a
    search-result excerpt; full text not accessible.
  - KOERI hosts a copy of the Early-est user guide under
    `koeri.boun.edu.tr/earlyest/` [8]. This suggests, but does not prove, that KOERI
    runs an Early-est (FP-based) instance.
- **Search range.** Vassallo et al. searched each FP parameter between 0 and twice its
  default. Their fitness function depended "significantly" on all five parameters [2].
  Their optimised S1/S2 values are close to 10. T_filter and T_up differ by up to a
  factor of 3 from the defaults [1].
- **Long windows.** A 100 Hz local network would plausibly use T_long = 10–12 s and
  T_up = 0.2–0.4 s. In the sources above, longer T_filter values (≥10 s) appear only
  in teleseismic configurations (Early-est).

## 3. Implementations

| Implementation | Language | License | Where | Maintained? | Real-time (sample/packet) | Status |
|---|---|---|---|---|---|---|
| Lomax FilterPicker5, C | C (≈1,000 lines incl. headers and PickData) | GPL-2.0-or-later (file headers + `COPYING`) | http://alomax.free.fr/FilterPicker/ → `archive/picker_20130212.tgz` [3, 4] | Last release 2013; no public VCS | **Yes.** `Pick()` takes packets and a per-channel memory struct (`useMemory=1`). The test harness feeds 10 s packets [3] | Confirmed (downloaded, read) |
| Lomax Java (`net.alomax.timedom.FilterPicker5`, SeisGram2K) | Java | GNU GPL (version not stated on page) | http://alomax.net/projects/java/ [20] | Unknown | Yes in SeisGram2K (assumed from [1]) | Page confirmed; source not read |
| Earthworm `pick_FP` | C | Module files GPL-2.0-or-later. The repository root `LICENSE` is MIT, but the file headers govern | https://gitlab.com/seismic-software/earthworm/-/tree/master/src/seismic_processing/pick_FP [5] | **Yes.** Commits 2024-11 (overlap handling) and 2026-10 | **Yes.** Reads TRACEBUF2 from a ring, picks each channel independently, releases picks immediately. Interpolates gaps up to `MaxGap` samples, otherwise restarts. Standalone tracebuf-file mode [5, 6] | Confirmed. Its picker copy predates v5.0.3 (no polarity weight) |
| Early-est | C | GPL-3.0-or-later (picker sub-directory GPL-2.0) | Third-party mirror https://github.com/comoglu/early-est (v1.2.9xDEV, 2023) [7] | Mirror pushed 2026-07 | Yes (SeedLink monitor) | Mirror confirmed. Official upstream download location **unverified** |
| PRESTo | C++ | "Free (GPL source code)" per website | http://www.prestoews.org/ (download behind a registration form) [16] | v1.0 (2016-03-11), nothing newer listed | Yes (SeedLink EEW platform) | Website confirmed. Source not inspected, so the FP version inside is **unverified** |
| SeisComP | — | — | `scautopick` lists AIC, BK and GFZ (re)pickers [21]. A code search of the SeisComP GitHub org found nothing | — | — | **No FilterPicker found** |
| ObsPy | — | — | `obspy.signal.trigger` has STA/LTA variants, Z-detector, `ar_pick`, `pk_baer` [22] | — | — | **Not in ObsPy** |
| SeisBench | — | — | The model list has no FilterPicker [23] | — | — | **Not in SeisBench** |
| `filterpicker` (M. Bagagli) | Python + NumPy | MIT | https://github.com/mbagagli/filterpicker, `pip install filterpicker` (1.1.0, 2022-05-02), Zenodo doi:10.5281/zenodo.6511985 [10] | Last push 2022-05 | **No.** "At the moment the picker works on the entire trace, not suitable for online picking" [10] | Confirmed (read source). Not a port of the C code; see below |
| Kamer, MATLAB | MATLAB | Not read (File Exchange "View License") | https://www.mathworks.com/matlabcentral/fileexchange/69211 [11] | v1.0.0, 2018-10-25 | Whole trace (assumed; the Python port derives from it) | Page confirmed; license **unverified** |
| PhasePApy FBpicker | Python | **No license file** (GitHub reports none, so all rights reserved by default) | https://github.com/austinholland/PhasePApy [12] | Last push 2018 | No (whole trace, NumPy) | Confirmed. A **different algorithm** (see below) |
| SEISAN `AUTOPHASE` | ? | ? | SEISAN manual §21.2 [24] | ? | Also used in RTQUAKE | **Unverified**. The manual says it is "based on FilterPicker with modifications". Page now redirects to a 404 |
| `calum-chamberlain/rt2detection` | C (copy of FilterPicker5) | No license file | https://github.com/calum-chamberlain/rt2detection | Last push 2017 | Packet mode | Exists; not relevant |

### Where the reimplementations differ from the C reference

**Bagagli `filterpicker`** [10] is a pure-Python per-sample loop:
- It loops over bands `n = 1 … N−1`, so it **skips band 0** (T = Δt). For T_filter = 3 s
  at 100 Hz it uses 8 bands (0.02–2.56 s), against 10 in C.
- It sets CF values below 1 to 0.
- It tests S2 with a sliding average (`np.convolve`) instead of the
  integral-from-trigger.
- It reports no polarity.
- It stores three float64 arrays of N_band × N samples, about 216 bytes per sample. That
  is about 78 MB per channel-hour, or about 8 GB per channel for the whole 102.8 h quiet
  set, so it would have to be chunked, with a T_long overlap at each boundary.

Expect picks close to the C reference, but not identical.

**PhasePApy FBpicker** [12, 25] is a different picker inspired by FP:
- 4-pole Butterworth octave band-passes, starting at `freqmin`;
- moving-window RMS or standard-deviation statistics;
- a dynamic threshold (moving average + n_sigma);
- noise-ratio and polarity heuristics.

Its parameters (`t_long`, `t_ma`, `nsigma`, `t_up`, …) do not map onto S1/S2.

**Reference code quirks**, relevant to vendoring [3, 5, 7]:
- `ew_bridge.h` *defines* `char message_str[1024]` in a header included by three
  translation units. Earthworm changed it to `extern` and Early-est to `static`. With
  GCC ≥ 10 (`-fno-common` by default) the original is expected to fail at link time
  (expected, not tested).
- `fmax`/`fmin` are `#define`d as macros, which clash with `<cmath>` in C++.
- `resultType` is a mutable global.
- v5.0.3 calls `info(message_str)` twice per accepted pick, printing a stale buffer to
  stdout. Earthworm and Early-est comment these calls out.
- Pick indices are relative to the current packet and can be negative (the onset lies
  in an earlier packet). The caller converts them using the packet start time [3, 5].
- Earthworm's `call_FilterPicker.c` copies `nsamp + 1` samples (`i <= length`), an
  off-by-one read [5].
- Per-channel state is entirely in the memory struct. Apart from the globals above, the
  code is reentrant.

## 4. Benchmarks

**Vassallo et al. (2012), optimised FP vs optimised and suggested `pick_ew` (Allen
picker)** [1, 2]:
- *Data:* ISNet, southern Italy, 125 Hz. 5,048 manually picked vertical traces of
  local/regional events, mostly M_L < 2.5, plus 948 high-noise traces.
- *Correct picks* (within 2 s of the manual pick): FP_OPT 90% on both broadband and
  short-period traces, against 74% / 87% for `pick_ew` optimised and 59% / 53% with
  suggested parameters [2, Table 4].
- *As summarised in [1]:*
  - picks matching manual picks: FP 4559 vs 4271;
  - manual picks missed: 489 vs 777;
  - false picks: 537 vs 644;
  - by noise class: FP missed 1% vs 6% in the lowest-noise class, 22% vs 29% in the
    highest.
- *Picks per trace:*
  - noise-only traces with more than 4 picks: 56 (6%) for FP_OPT, against 299 (32%)
    and 287 (30%) for `pick_ew` with suggested and optimised parameters [2];
  - event traces with more than 4 picks: 358 for FP_OPT, 330 for `pick_ew`
    optimised (7% each) and 530 (10%) for `pick_ew` suggested [2].
- *Timing residuals* are similar for both pickers: mean near 0, σ ≈ 0.1–0.5 s depending
  on class [2, Fig. 6].
- *Full chain* (picks → Earthworm association → events) [1]:

  | Dataset | FP | `pick_ew` |
  |---|---|---|
  | 301 events inside the network | 282 detected | 265 detected |
  | 104 events outside the network | 76 detected | 71 detected |
  | 49 storm-noise "false event" periods | **10** declarations | **31** declarations |

  All events with M > 2.5 were detected by both. Below M 2.5, FP detected 5–10% more.

  This is the closest published analogue to ayzek's quiet-data false-alarm measure.

**Mousavi et al. (2020), EQTransformer paper** [13]:
- *Data:* STEAD test set, over 113k one-minute three-component waveforms. A pick is
  correct within 0.5 s.
- *Implementation:* the "FilterPicker" baseline was "based on the python implementation
  in PhasePApy", i.e. FBpicker, not Lomax's code. No parameters are reported.
- *Results* (Pr = precision, Re = recall):

  | Phase | Model | Pr | Re | F1 | μ ± σ (s) |
  |---|---|---|---|---|---|
  | P | FilterPicker | 0.95 | 0.82 | 0.88 | −0.01 ± 0.08 |
  | P | Kurtosis | 0.94 | 0.79 | 0.86 | |
  | P | AIC | 0.92 | 0.83 | 0.87 | |
  | P | PhaseNet | 0.96 | 0.96 | 0.96 | |
  | P | EQTransformer | 0.99 | 0.99 | 0.99 | |
  | S | FilterPicker | 0.61 | 0.41 | 0.49 | |
  | S | PhaseNet | | | 0.94 | |
  | S | EQTransformer | | | 0.98 | |

- *Run time:* see §1.6.
- These are per-trace picking scores on pre-cut event windows. They are not detection
  or false-alarm rates on continuous data.

**Wang et al. (2023), Frontiers in Earth Science** [14]:
- *Gansu (China) 100 Hz dataset*, 283,467 records:
  - FilterPicker: P precision 0.84 / recall 0.69; S 0.46 / 0.22;
  - PhaseNet: P 0.91 / 0.88; S 0.87 / 0.84;
  - STA/LTA event detection: 0.73 / 0.65.
- *Their STEAD numbers* for FilterPicker are identical to Mousavi's, so they were
  probably copied.
- *Missing details:* no FP implementation, parameters or match tolerance are given.

**Other benchmarks:**
- The SeisBench benchmark (Münchmeyer et al., 2022) evaluates deep-learning pickers.
  I could not access the full text to check its classical baselines. FilterPicker is
  not in the SeisBench model zoo [23]. **Unverified** whether it is used there at all.
- I found no published head-to-head of Lomax's FilterPicker against PhaseNet or
  EQTransformer on continuous data with network association. A comparison inside ayzek
  would be new.

## 5. Effort in our case

### What ayzek already provides

- `--detector stalta` swaps the learned detector for a sample-by-sample trigger
  (`src/stalta.hpp`). Everything downstream is unchanged:
  - picker window;
  - association with `min_stations = 2`, travel-time slack and coda window;
  - location, magnitude, catalogue scoring (`docs/impl/09-sta-lta.md`).
- `--scores-in DIR` replaces the detector with a per-window score cache. This is how
  PhaseNet and EQTransformer were brought in (`tools/seisbench_scores.py`).

FilterPicker can enter through either door.

### Route A: in the C++ pipeline

**A1. Vendor Lomax's C.**
- *License.* The code is GPL-2.0-or-later and ayzek is Apache-2.0 (`LICENSE`). The FSF
  regards Apache-2.0 as compatible with GPLv3 but not GPLv2, and "GPLv3 software cannot
  be included in Apache projects" [15]. A binary linking the vendored picker would
  therefore have to be distributed as GPLv3 as a whole.
  - That can be fine for a public research repo: keep the GPL files with their notices
    in e.g. `third_party/filterpicker/` with a provenance README, behind a meson option
    that is off by default.
  - It conflicts with keeping ayzek a standalone Apache-2.0 product, and it needs an
    explicit licensing decision by the authors. Not legal advice.
- *Code work:*
  - fix `message_str`, the `fmax`/`fmin` macros and the stray `info()` calls;
  - wrap `Pick()` in a C++ class fed per sample or per 0.5 s chunk;
  - restart on gaps, as `pick_FP` does;
  - convert packet-relative indices to absolute times;
  - map a declared pick to an ayzek detection.
- *Estimate:* about 1–2 days of code and tests, plus the licensing decision.

**A2. Independent reimplementation under Apache-2.0.**
- *Basis:* the paper plus erratum [1, 4]. Every rule in §1 is stated in the paper,
  including:
  - S1/2 clamp;
  - 5·S1 spike limit;
  - re-arm below 2;
  - T_k/40 minimum uncertainty;
  - 66% polarity ratio.
- *Size:* about 200–300 lines in the style of `src/stalta.hpp`.
- *Validation:* compare against picks produced offline by the GPL reference. Lomax's
  tarball includes test SAC files and reference picks from both the Java and the C
  versions (`test_data/seisgram2k.pick`) [3]. Running a GPL program to produce test
  data does not put its license on the data.
- *Caveat:* write it from the paper, not by translating the C file.
- *Estimate:* about 3–5 days including tests. This keeps ayzek standalone and gives a
  real-time FP inside the pipeline, latency included (pick ≤ T_up after the trigger).

### Route B: offline on the replays

**B1 (preferred offline): Lomax's C behind a thin wrapper.**
- *Steps:*
  - compile it as a small shared library or CLI kept outside ayzek's build;
  - read miniSEED with ObsPy (already a tool dependency);
  - feed the vertical channel in packets, restarting at gaps;
  - write per-station picks (time, uncertainty, polarity, band period, CF strength).
- *Getting into ayzek:* either
  - write an ayzek score cache for `--scores-in`, with a window score of 1 when a
    declared pick falls where the window convention expects P, else 0, so association
    and scoring stay identical to the other detectors; or
  - write a separate associator applying the same 2-station rule.
- *License:* keep the wrapper and the C files in a separately licensed directory or
  repo, or fetch the tarball at run time with a checksum. Then ayzek's own code and
  binaries never link GPL code.
- *Estimate:* about 1–2 days, plus about 1 day for the S1/S2 sweep. Compute is
  negligible (§1.6).

**B2: Bagagli's Python package.**
- *Pros:* `pip install`, MIT license.
- *Cons:*
  - whole-trace only;
  - pure-Python loop, which I estimate at minutes to tens of minutes per channel for the
    quiet set (not measured);
  - about 8 GB of arrays per quiet-set channel unless chunked;
  - not numerically the reference algorithm (§3).
- *Use:* a cross-check only.

### Pitfalls

1. **A picker is not a detector.**
   - FP declares onsets per channel; ayzek alarms on two or more stations.
   - The fair comparison is "FP picks → the same association rule".
   - FP's re-arm rule (CF < 2 in all bands) and ayzek's 15 s minimum between triggers
     interact. Decide which applies.
   - FP will also pick S, later phases and coda bursts. ayzek's association has to
     absorb them as it does STA/LTA re-triggers.
2. **The time convention differs.** The STA/LTA integration dates a trigger 3.5 s
   before the trigger sample, to match the detector windows (`09-sta-lta.md`). FP
   already returns the onset estimate. The mapping to window scores or detections must
   be chosen so that the downstream picker window and catalogue matching sit where they
   do for the other detectors.
3. **Broadband sensitivity on quiet data.**
   - FP picks teleseismic P (paper Fig. 2) and regional events. On "quiet" windows
     these become coincident multi-station picks.
   - Whether they count as false alarms must be decided exactly as for the other
     detectors (catalogue radius 250 km).
   - At 100 Hz, band 0 has T = 0.01 s. Cultural noise and glitches can trigger the
     highest bands. FP has no band-limit parameter other than T_filter; a pre-filter
     would depart from "standard" FP.
4. **Parameters.**
   - Report the defaults (S1 = S2 = 10, T_filter 3 s, T_long 5 s, T_up 0.2 s at
     100 Hz).
   - State which implementation's defaults are used, since Earthworm's autoset gives
     6 / 10 / 0.4 s at 100 Hz.
   - Sweep S1 (and S2) for an operating curve, as for STA/LTA. Tuning on the
     evaluation replays risks overfitting; keep tuning and test sequences separate.
5. **Warm-up and gaps.** No picks during the first T_long after start or after each
   restart. Gap handling (interpolate short gaps vs restart) changes pick counts near
   gaps.
6. **Benchmark wording.** The widely quoted EQTransformer-paper FilterPicker numbers
   come from PhasePApy's FBpicker, a different algorithm, on pre-cut STEAD traces. Do
   not present them as Lomax-FP performance on continuous data.
7. **Components.** FP is single-channel. `pick_FP` and Vassallo et al. pick on vertical
   channels. Use Z only to match standard practice; ayzek's STA/LTA defaults to
   vertical only too.

## Recommendation for ayzek

- **Do it offline first (Route B1)**, with Lomax's own C FilterPicker5 v5.0.3 [3], or
  the Earthworm copy [5], as the reference implementation.
  - Wrap it outside ayzek's Apache-2.0 tree.
  - Run it on the Z channel of the replays and the 102.8 h quiet set with paper
    defaults, plus an S1/S2 sweep.
  - Feed the declared picks through ayzek's existing association and scoring via a
    `--scores-in` cache, so only the detection stage changes, as was done for
    STA/LTA, PhaseNet and EQTransformer.
  - About 2–3 days in total, negligible compute.
  - This answers the referee with the actual operational algorithm (as run in Earthworm
    `pick_FP`, PRESTo and Early-est) without putting GPL code into ayzek.
- **Real-time FP inside the C++ pipeline** (e.g. for latency figures) is optional, if
  the paper needs it. Prefer a clean Apache-2.0 reimplementation from the paper (A2,
  about 3–5 days), validated against the B1 outputs, over vendoring the GPL source (A1),
  which would make distributed ayzek binaries GPLv3.
- **Do not use** the Python port or PhasePApy as "FilterPicker" in the paper, except as
  a disclosed cross-check.

---

## References

1. Lomax, A., C. Satriano, and M. Vassallo (2012). Automatic picker developments and
   optimization: FilterPicker—a robust, broadband picker for real-time seismic
   monitoring and earthquake early warning. *Seismological Research Letters* 83(3),
   531–540. https://doi.org/10.1785/gssrl.83.3.531. Author PDF:
   http://alomax.free.fr/articles/2012__Lomax_et__Automatic_picker_developments_and_optimization_FilterPicker_a_robust_broadband_picker_for_real-time_seismic_monitoring_and_earthquake_early-warning__SRL.pdf
2. Vassallo, M., C. Satriano, and A. Lomax (2012). Automatic picker developments and
   optimization: A strategy for improving the performances of automatic phase pickers.
   *Seismological Research Letters* 83(3), 541–554. https://doi.org/10.1785/gssrl.83.3.541.
   Author PDF:
   http://alomax.free.fr/articles/2012__Vassallo_et__Automatic_picker_developments_and_optimization_a_strategy_for_improving_the_performances_of_automatic_phase_pickers__SRL.pdf
3. Lomax, A. FilterPicker C source, `picker_20130212.tgz` (FilterPicker5.c v5.0.3;
   `FilterPicker5_Memory.c`, `PickData.c`, `DOC_FP5.txt`, `COPYING` = GPL v2, test data).
   http://alomax.free.fr/FilterPicker/archive/
4. Lomax, A. FilterPicker home page (downloads, Java class name, implementations list)
   and erratum to eq. (5). http://alomax.free.fr/FilterPicker/ ;
   http://alomax.free.fr/FilterPicker/erratum.html
5. Earthworm, `pick_FP` source (README, `call_FilterPicker.c`, `pick_FP.sta`,
   `picker/`), GitLab `seismic-software/earthworm`, master branch, accessed 2026-10-26.
   https://gitlab.com/seismic-software/earthworm/-/tree/master/src/seismic_processing/pick_FP ;
   repository license: https://gitlab.com/seismic-software/earthworm/-/blob/master/LICENSE
6. Earthworm documentation: pick_FP overview and commands.
   https://folkworm.ceri.memphis.edu/ew-doc/ovr/pick_FP_ovr.html ;
   https://folkworm.ceri.memphis.edu/ew-doc/cmd/pick_FP_cmd.html
7. Early-est source mirror (v1.2.9xDEV, GPL-3.0; `picker/`,
   `timedomain_processing/timedomain_processing.c`). https://github.com/comoglu/early-est
8. Lomax, A. Users Guide for Early-est (2016-06-27), copy hosted by KOERI.
   http://www.koeri.boun.edu.tr/earlyest/warnings/early-est_users_guide.pdf
9. Lancieri, M., A. Fuenzalida, S. Ruiz, and R. Madariaga (2011). Magnitude scaling of
   early-warning parameters for the Mw 7.8 Tocopilla, Chile, earthquake and its
   aftershocks. *BSSA* 101(2), 447–463. https://doi.org/10.1785/0120100045
   (cited in [1] as an FP application on continuous records).
10. Bagagli, M. (2022). filterpicker: general purpose, broad-band, phase detector and
    picker, v1.1.0. https://doi.org/10.5281/zenodo.6511985 ;
    https://github.com/mbagagli/filterpicker ; https://pypi.org/project/filterpicker/ ;
    https://filterpicker.readthedocs.io/
11. Kamer, Y. (2018). FilterPicker: a robust, broadband phase detector and picker
    (MATLAB), v1.0.0. MATLAB Central File Exchange 69211.
    https://www.mathworks.com/matlabcentral/fileexchange/69211-filterpicker-a-robust-broadband-phase-detector-and-picker
12. PhasePApy source (FBpicker: `fbpicker.py`, `cf_fb.py`).
    https://github.com/austinholland/PhasePApy
13. Mousavi, S. M., W. L. Ellsworth, W. Zhu, L. Y. Chuang, and G. C. Beroza (2020).
    Earthquake transformer—an attentive deep-learning model for simultaneous earthquake
    detection and phase picking. *Nature Communications* 11, 3952.
    https://doi.org/10.1038/s41467-020-17591-w
14. Wang, S., F. Liu, X. Yin, K. Chen, and R. Cai (2023). Employing
    convolution-enhanced attention mechanisms for earthquake detection and phase picking
    models. *Frontiers in Earth Science* 11, 1283857.
    https://doi.org/10.3389/feart.2023.1283857
15. Apache Software Foundation. Apache License v2.0 and GPL compatibility.
    https://www.apache.org/licenses/GPL-compatibility.html
16. PRESTo project website (v1.0, license "Free (GPL source code)").
    http://www.prestoews.org/ — and Satriano, C., L. Elia, C. Martino, M. Lancieri,
    A. Zollo, and G. Iannaccone (2011). PRESTo, the earthquake early warning system for
    Southern Italy: Concepts, capabilities and future perspectives. *Soil Dynamics and
    Earthquake Engineering* 31(2), 137–153. https://doi.org/10.1016/j.soildyn.2010.06.008
17. Xu, et al. (2020). Earthquake early warning system in Liaoning, China based on
    PRESTo. *Earthquake Science* 33, 281–292. https://doi.org/10.29382/eqs-2020-0281-01
    (parameter values **unverified**: full text not accessible).
18. Picozzi, M., A. Emolo, C. Martino, A. Zollo, et al. (2015). Earthquake early warning
    system for schools: A feasibility study in southern Italy. *Seismological Research
    Letters* 86(2A), 398–412. https://doi.org/10.1785/0220140194 (parameter values
    **unverified**: full text not accessible).
19. Clinton, J., A. Zollo, A. Marmureanu, C. Zulfikar, et al. (2016). State-of-the art
    and future of earthquake early warning in the European region. *Bulletin of
    Earthquake Engineering* 14, 2441–2458. https://doi.org/10.1007/s10518-016-9922-7
    (KOERI/PRESTo statement **unverified**: full text not accessible).
20. Lomax, A. net.alomax and related Java source code (GNU GPL).
    http://alomax.net/projects/java/
21. SeisComP documentation, `scautopick`. https://www.seiscomp.de/doc/apps/scautopick.html
22. ObsPy documentation, `obspy.signal` (trigger module).
    https://docs.obspy.org/packages/obspy.signal.html
23. SeisBench model implementations.
    https://github.com/seisbench/seisbench/tree/main/seisbench/models
24. SEISAN manual §21.2, AUTOPHASE (now redirects; content seen only in search results,
    **unverified**). https://seis.geus.net/software/seisan/node153.html
25. Chen, C., and A. A. Holland (2016). PhasePApy: A robust pure Python package for
    automatic identification of seismic phases. *Seismological Research Letters* 87(6),
    1384–1396. https://doi.org/10.1785/0220160019
