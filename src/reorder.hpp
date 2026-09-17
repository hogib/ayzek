#pragma once

// Reorders records of one channel into a contiguous stream of samples.
//
// Input: records in arrival order, each with an absolute start position in
// samples (start_ns / 10'000'000 at 100 Hz). Output: contiguous runs in
// strictly increasing position, and explicit gaps.
//
// Rules:
//   - a record starting at the frontier (next expected position) is emitted
//     immediately without copying
//   - a record starting beyond the frontier is held until the preceding hole
//     is filled
//   - a hole is declared a gap once data more than `max_lateness` samples
//     beyond it has been received; the frontier then moves past it
//   - a record arriving after its hole was declared a gap is counted as stale
//     and discarded, because the samples after it have already been emitted
//   - samples overlapping already-emitted data are trimmed
//
// About 0.3% of horizontal records in the TDVMS archive arrive 7.7 s to ~10 min
// late (docs/DESIGN.md). The real-time pipeline uses max_lateness = 0; larger
// values are used for offline replay analysis.

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>

namespace ayzek {

struct ReorderStats {
  std::uint64_t records = 0;
  std::uint64_t in_order = 0; // arrived exactly at the frontier
  std::uint64_t held = 0;     // arrived early and waited for a hole to fill
  std::uint64_t stale = 0;    // wholly behind the frontier, dropped
  std::uint64_t trimmed_samples =
      0; // leading samples behind the frontier, dropped
  std::uint64_t gaps = 0;
  std::uint64_t gap_samples = 0;
  std::uint64_t forced_gaps = 0; // declared early because the pool was full
  std::uint64_t rejected = 0;    // longer than a slot can hold
  // Largest amount, in samples, by which a record's end trailed the newest
  // data already received. Same quantity as the max_lateness setting.
  std::uint64_t max_lateness = 0;
};

// `MaxHeld`: records that can be held at once. `MaxRecordSamples`: samples per
// record. Storage is a fixed array; offer() does not allocate.
template <std::size_t MaxHeld, std::size_t MaxRecordSamples> class Reorderer {
  struct Slot {
    std::uint64_t pos = 0;
    std::uint32_t n = 0;
    bool used = false;
    std::array<std::int32_t, MaxRecordSamples> data{};
  };

  std::array<Slot, MaxHeld> pool_{};
  std::size_t count_ = 0;
  std::uint64_t max_lateness_;

  // The frontier is set once max_lateness worth of data has been seen, so that
  // a late-arriving first record is not dropped.
  bool started_ = false;
  std::uint64_t next_ = 0;     // position of the next sample to emit
  std::uint64_t seen_end_ = 0; // highest record end ever offered

  ReorderStats stats_{};

public:
  // `max_lateness_samples`: how much later data may arrive before a hole is
  // declared a gap. Zero declares every hole a gap immediately.
  explicit Reorderer(std::uint64_t max_lateness_samples) noexcept
      : max_lateness_(max_lateness_samples) {}

  [[nodiscard]] const ReorderStats &stats() const noexcept { return stats_; }
  [[nodiscard]] std::size_t held() const noexcept { return count_; }
  [[nodiscard]] std::uint64_t frontier() const noexcept { return next_; }

  // Offers one record. `sink(pos, samples)` receives contiguous runs in
  // increasing position; `gap(from, to)` reports [from, to) as missing. The
  // callbacks are template parameters, so calls are resolved at compile time.
  template <typename Sink, typename Gap>
    requires std::invocable<Sink &, std::uint64_t,
                            std::span<const std::int32_t>> &&
             std::invocable<Gap &, std::uint64_t, std::uint64_t>
  void offer(std::uint64_t pos, std::span<const std::int32_t> samples,
             Sink &&sink, Gap &&gap) {
    const std::uint64_t n = samples.size();
    if (n == 0)
      return;
    ++stats_.records;

    const std::uint64_t end = pos + n;
    if (seen_end_ > end)
      stats_.max_lateness = std::max(stats_.max_lateness, seen_end_ - end);
    seen_end_ = std::max(seen_end_, end);

    if (started_) {
      if (end <= next_) { // entirely behind: too late
        ++stats_.stale;
        return;
      }
      if (pos <= next_) { // at, or straddling, the frontier
        const std::uint64_t skip = next_ - pos;
        stats_.trimmed_samples += skip;
        if (skip == 0)
          ++stats_.in_order;
        sink(next_, samples.subspan(static_cast<std::size_t>(skip)));
        next_ = end;
        drain(sink);
        resolve(sink, gap, false);
        return;
      }
    }

    // Ahead of the frontier, or frontier not yet established: hold it.
    if (n > MaxRecordSamples) {
      // Cannot happen while MaxRecordSamples matches the decoder's limit;
      // counted if it does.
      ++stats_.rejected;
      return;
    }
    if (count_ == MaxHeld)
      resolve(sink, gap, true);
    Slot &s = free_slot();
    s.pos = pos;
    s.n = static_cast<std::uint32_t>(n);
    s.used = true;
    std::copy(samples.begin(), samples.end(), s.data.begin());
    ++count_;
    ++stats_.held;

    resolve(sink, gap, false);
  }

  // End of stream: emits all held records, declaring remaining holes as gaps.
  template <typename Sink, typename Gap> void flush(Sink &&sink, Gap &&gap) {
    while (count_ > 0) {
      if (!started_)
        start_at_earliest();
      Slot *m = earliest();
      if (m->pos > next_) {
        declare_gap(gap, m->pos);
      }
      drain(sink);
    }
  }

private:
  Slot &free_slot() noexcept {
    for (Slot &s : pool_)
      if (!s.used)
        return s;
    return pool_[0]; // guarded by the caller's count_ check
  }

  // Held record with the lowest position. Linear scan; MaxHeld is small.
  Slot *earliest() noexcept {
    Slot *best = nullptr;
    for (Slot &s : pool_)
      if (s.used && (best == nullptr || s.pos < best->pos))
        best = &s;
    return best;
  }

  void release(Slot &s) noexcept {
    s.used = false;
    --count_;
  }

  void start_at_earliest() noexcept {
    next_ = earliest()->pos;
    started_ = true;
  }

  template <typename Gap> void declare_gap(Gap &gap, std::uint64_t to) {
    gap(next_, to);
    ++stats_.gaps;
    stats_.gap_samples += to - next_;
    next_ = to;
  }

  // Emits, in order, every held record that now starts at or before the
  // frontier.
  template <typename Sink> void drain(Sink &sink) {
    for (;;) {
      Slot *m = earliest();
      if (m == nullptr || m->pos > next_)
        return;
      const std::uint64_t end = m->pos + m->n;
      if (end <= next_) {
        ++stats_.stale; // duplicate of emitted data
      } else {
        const std::uint64_t skip = next_ - m->pos;
        stats_.trimmed_samples += skip;
        sink(next_,
             std::span<const std::int32_t>(m->data.data() + skip, m->n - skip));
        next_ = end;
      }
      release(*m);
    }
  }

  // Declares gaps for holes that have waited longer than max_lateness, or for
  // the oldest hole when the pool is full.
  template <typename Sink, typename Gap>
  void resolve(Sink &sink, Gap &gap, bool pool_full) {
    while (count_ > 0) {
      if (!started_) {
        const std::uint64_t first = earliest()->pos;
        if (!pool_full && seen_end_ - first <= max_lateness_)
          return;
        start_at_earliest();
        drain(sink);
        continue;
      }
      Slot *m = earliest();
      // Waiting time is measured from the end of the hole (the start of the
      // next held record), which is the quantity stats().max_lateness
      // reports. Measuring from the start of the hole would add the hole's
      // own length to its lateness.
      const bool waited_out = seen_end_ - m->pos > max_lateness_;
      if (!waited_out && !pool_full)
        return;
      if (m->pos > next_) {
        if (pool_full && !waited_out)
          ++stats_.forced_gaps;
        declare_gap(gap, m->pos);
      }
      drain(sink);
      pool_full = false; // one freed slot suffices
    }
  }
};

} // namespace ayzek
