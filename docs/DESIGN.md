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
