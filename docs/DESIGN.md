# ayzek — a real-time EEW pipeline

Detect, pick and size earthquakes on live continuous streams from three
stations, in C++23, with hand-written inference.

## What this is for

Two goals, honestly ranked. Learning C++ on a problem where the language's
actual strengths matter, and producing something that runs. Every stage is
checkable against Python that already works, so debugging is never "is it the
code or the science".

---

## The architecture question

The proposal was: one ring buffer per station reading all channels at once,
then another buffer reading from the first and doing each model's preparation,
and so on. Three of those four ideas are right.

### Right: three stations

Three is the minimum that locates an event. Differential arrival times cancel
the unknown origin time, so N stations give N-1 equations, and three give
latitude and longitude with depth fixed. A fourth adds a residual but not a new
answer.

Two consequences, both already worked out in the sibling project's
`locate-from-station-triplets` note and worth carrying here: with exactly three
stations the fit has no residual, so it cannot report its own quality and must
be validated externally; and two hyperbolae can intersect twice, which arrival
times alone cannot disambiguate.

### Right: do not split the components

E, N and Z are sampled on the same clock at the same instants and every stage
consumes all three together. Splitting them into separate buffers would mean
three indices to keep synchronised and no benefit. One sample is

```cpp
struct Sample { float z, n, e; };   // 12 bytes, one index
```

and the ring is `Ring<Sample, N>`. Reading a window touches contiguous memory
and pulls all three components into cache together, which is exactly the access
pattern.

### Right: one ring per station

Stations have independent clocks, independent packet arrival, independent gaps
and independent event arrival times. They cannot share a buffer, and trying to
would make the gap logic intractable.

### Needs refinement: a ring between *every* stage

A ring buffer exists to decouple a producer and a consumer that run at
different rates or on different threads. Where both sides run in the same
thread at the same rate, a ring is just a slower function call: the window
extractor and the model do not need one between them, they need a view.

So the pipeline is **not** ring -> ring -> ring -> ring. It is two rings and
then a chain of transforms over views.

**But two rings is correct**, and for a reason worth stating because it is the
part that is easy to get wrong.

### Why the second ring genuinely earns its place

The obvious design filters each window as it is extracted. That is wrong twice.
Overlapping windows refilter the same samples, which is waste; and far more
importantly, an IIR filter applied to a short isolated window has **edge
transients**. The filter's state has to spin up, so the first fraction of a
second of every window is corrupted by the filter itself -- right where a P
arrival might be.

Filtering continuously, as samples arrive, with the filter state carried across
packet boundaries, removes both problems. The filtered samples then need
somewhere to live, and that somewhere is the second ring.

```
SeedLink ─→ decode ─→ [RING 1: raw @ 100 Hz]
                            │  continuous IIR, state carried across packets
                            ↓
                      [RING 2: filtered @ target rate]
                            │  mdspan views, zero copy
                            ↓
                      detector ─→ trigger ─→ picker + magnitude
```

Everything after RING 2 is a view. `std::mdspan` over the ring's storage gives
a window as a `(n_samples, 3)` matrix with no allocation and no copy.

---

## Threading

Data rate is trivial -- three stations at 100 Hz three-component is 900 floats
a second. **Latency is the constraint, not throughput.** So the threading model
should minimise hops, not maximise parallelism.

```
thread per station:  socket read -> decode -> RING 1 -> filter -> RING 2
single processor:    poll all three RING 2s, extract windows, run models
```

Each ring stays **single-producer, single-consumer**, which is the one lock-free
structure worth hand-writing: two atomics, acquire/release ordering, no mutex,
no allocation. Filtering on the reader thread rather than in a fourth stage
keeps the hop count down and keeps the filter state owned by exactly one thread.

---

## The part that will actually bite: gaps

A ring indexed by sample number silently diverges from wall-clock time the
moment a packet is missing. Every timing bug in this kind of system traces back
to that.

Rules, adopted up front:

- Every ring carries `(sample_index, timestamp)` anchors, not a start time plus
  a count.
- A detected gap is **filled explicitly** with a marked value and the anchor is
  re-established, so index arithmetic stays exact.
- **A window that spans a gap is never emitted.** The sibling cascade paper
  already holds this line ("windows never span a gap in the record"); breaking
  it here would silently feed the model fabricated samples.
- Overlapping or out-of-order packets are resolved by timestamp, not arrival
  order.

---

## Measured: a few horizontal records arrive very late

A full 21-day DEMI chunk, decoded and replayed.

**The file is sorted by channel, not by time** -- all HHZ, then HHE, then HHN --
so replay feeds each channel's records in their own file order.

**The vertical is never late. The horizontals rarely are, but by a lot:**

| channel | records | late | steps back | smallest | median | largest |
|---|---:|---:|---:|---:|---:|---:|
| HHZ | 488,059 | 0 | 0 | -- | -- | -- |
| HHE | 498,768 | 1,528 (0.31%) | 757 | 7.7 s | 39.2 s | 615.9 s |
| HHN | 497,247 | 1,585 (0.32%) | 811 | 7.5 s | 43.5 s | 609.2 s |

"Steps back" are places where file order jumps backwards in time, sized in the
last three columns. Of the late records, 3 of 3,113 trail the newest data by
under 5 s, about 9% by
5-10 s, 31% by 10-30 s, 21% by 30-120 s, and 40% by more than two minutes. They
come in ~700 short runs (at most 18 records) spread over ~225 of the 504 hours.
That shape is not network jitter, which would crowd near zero. It fits packets
lost and recovered minutes later, then stored in the order they arrived -- an
inference from the distribution, not something verified.

**Genuine gaps are rare once sorted:** 25 on HHZ, 41 on HHE, 58 on HHN, every
one at least a second long, with no sub-sample jitter and no overlaps.

ObsPy assembles 26, 2,052 and 2,176 traces from these channels because it
builds traces in file order, and each late record breaks the current trace.

### Correction

An earlier version of this section gave the horizontals' **largest** backward
jump as 7.7 s and 7.5 s. Those are the **smallest**: the analysis took the
maximum of negative time differences, which is the one closest to zero. The
largest are about ten minutes.

The argument built on the wrong number -- that aligning three components means
waiting up to 7.7 s for the latest horizontal, against a warning budget of a
few seconds -- led to a decision to run detection on Z alone. That decision is
reversed below.

### What the replay shows

`tools/replay_check` feeds each channel through the reorder buffer across a
range of `max_lateness` settings and checks every emitted sample and gap against
the time-sorted records.

- **HHZ at `max_lateness = 0` is exact:** all 181,427,887 samples match, all 25
  declared gaps sit on genuine holes, and nothing is delayed.
- **No usable window rescues late horizontals.** A 10 s window recovers about a
  third of them (HHE stale records 1,528 -> 1,070); a lossless one would need
  about ten minutes. Data that late is useless for early warning whatever the
  design, so no design waits for it.

### Decided: keep the three-component detector, for now

Since no design waits for late horizontals, the latency cost attributed to
aligning three components does not exist. What late packets actually cost is
**real-time gaps on the horizontals**: 485 on HHE and 513 on HHN over the 21 days,
against 25 on HHZ, losing about 0.3% of the horizontal data. They are not short: a
late run blanks the horizontals for a median of 10 s (90th percentile 14.5 s,
longest 85 s), and only 13% last 5 s or less. A three-component detector loses
every window those gaps touch -- with the detector's 6 s window, about 16 s of
windows per gap. A Z-only detector would recover that ~0.3%,
at the cost of retraining -- worth doing only if the loss turns out to matter, so
it is deferred until measured.

Consequences:

1. **One ring per component, all indexed by the same absolute position.** Kept
   despite the reversal, on independent grounds: gaps are per component -- Z runs
   on while E has a hole -- and a shared `Sample {z, n, e}` would force every hole
   into all three. `SpscRing<float, N>` needs no change.
2. **Every channel passes through the reorder buffer with `max_lateness = 0`.**
   Never wait; any hole is an immediate gap; any late record is counted as stale.
   The class still earns its place with no reordering to do, because it declares
   gaps, trims overlaps and catches duplicates.
3. **A three-component window is emitted only when all three components are
   contiguous across it.** Windows touching a horizontal gap are skipped and
   **counted**. That count, set against real detections, is what decides whether
   a Z-only detector is worth building.

   Bridging gaps with interpolated samples instead was considered. Bridging a
   gap costs no latency -- its far edge is known the moment it is declared -- but
   a typical gap is longer than the detector's whole window. A bridged window
   is really a Z-only window with fabricated horizontals, and the detector was
   trained with no gap augmentation. Whether it tolerates that is measurable on
   its existing test set without retraining, and is worth measuring before
   either bridging or a Z-only model.
4. **Late records are not lost to offline use.** Training data and the coda work
   keep reading the time-sorted archive; only the real-time path drops them.

Whether live SeedLink behaves like this archive is unverified until stage 7.
The design holds either way: if live horizontals arrive on time, these gaps
simply do not occur.

### Measured: positions are exact

Every record on every channel starts exactly on a 10 ms boundary -- one distinct
phase, 0.000 ms, across all 1,484,074 records of the DEMI chunk. All three
components share a sample clock locked to the epoch, so absolute position is
integer arithmetic rather than an estimate:

```cpp
std::uint64_t position = start_ns / 10'000'000;   // samples since 1970 at 100 Hz
```

Indexing three rings by one position is therefore exact, a gap is a jump in
position, and position times the sample period *is* the timestamp -- which
satisfies the anchor rule above without storing anchors at all. Ingest verifies
the alignment per record and counts violations, since other stations or live
data may not be locked this way.

## Inference, hand-written

No LibTorch, no ONNX. The models are small enough that writing the forward pass
is tractable and is most of the point:

| model | size | what it needs |
|---|---|---|
| detector | 142k params | conv1d stack, BiLSTM, attention |
| S picker | 95k params | conv1d stack, BiLSTM, attention, dense head |
| magnitude | small | conv1d + dense |

All three reduce to the same handful of kernels: 1D convolution, batch norm,
an LSTM cell, scaled dot-product attention, and a couple of activations.
Weights export from PyTorch as a flat binary plus a shape manifest.

Every kernel is validated against the Python it came from, to single-precision
agreement, on the same inputs. That check is not optional -- a transposed
weight matrix produces plausible-looking wrong answers rather than a crash.

---

## Stages

Each stage is a complete, useful, testable thing. Do not start the next until
the current one has a passing test.

**1. Ring buffer.** Fixed capacity, power of two, SPSC, no allocation in the hot
path, `std::span` and `std::mdspan` views. Tests: wraparound correctness,
overwrite detection, a producer and consumer thread hammering it under
ThreadSanitizer. *No dependencies. This is where most of the C++ is learned.*

**2. File replay.** Read a miniSEED file and feed the ring at real-time rate.
Realistic load, no network. Byte-compare the ring's output against the same
file read in Python. *Makes stage 1 provably correct on real data.*

**3. Continuous filtering.** Butterworth as second-order sections, state carried
across packet boundaries, into RING 2. Validate against `scipy.signal.sosfilt`
to single precision. *Where cache layout starts to matter.*

**4. Windowing.** Sliding `mdspan` views, gap-aware, never spanning a gap.

**5. Kernels.** conv1d, batchnorm, LSTM cell, attention, dense. Each one tested
against PyTorch output on the same weights and inputs, independently, before
anything is composed.

**6. The detector**, composed from stage 5, matching Python end to end.

**7. Live SeedLink.** `libslink` via a Meson wrap. *Deliberately late -- network
code is the least interesting part and the easiest to defer.*

**8. The cascade.** Detector triggers, picker and magnitude run, alert emitted
with latency measured at every hop.

**9. Location** from three stations, once picks exist. See the caveats above.

---

## Measuring it

The number that matters is **latency from sample arrival to alert**, broken
down per hop, reported as a distribution rather than a mean. `std::chrono`
timestamps at each boundary, written to a ring of their own.

A pipeline that is correct but takes four seconds is not an early warning
system.

---

## Build notes

Meson, C++23, no modules. Modules were tested on this toolchain: GCC 16.2 and
Clang 22.1 both compile them correctly by hand, but Meson 1.12 only wires
dyndep scanning for MSVC, and forcing it on via `-fmodules` reaches ninja with
a malformed dyndep. Revisit when that is fixed upstream; nothing here depends
on it.

What C++23 does give this project, and it is a lot: `std::mdspan` for zero-copy
windows over the ring, `std::span` for buffer views, `std::expected` for error
handling with no exceptions on the hot path, and `std::print` so iostream never
appears.
