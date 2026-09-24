#pragma once

// Writes one component's samples into its ring, in position order, without
// ever blocking inside a write and without ever dropping a sample.
//
// Why it exists. A ring write can only proceed while the processor keeps
// raising the floor, and the processor only advances while *all three*
// components of the station have data (it scores the minimum of their ends).
// Ingest runs one thread per station, so a component that waits for room
// inside its write stops the other two from being fed, and the station
// deadlocks until something gives. A gap longer than the ring -- 2^17 samples,
// 21.8 minutes at 100 Hz -- always does this, because filling it needs more
// room than the ring can ever offer at once.
//
// So writes that do not fit are queued and retried by `drain()`, which the
// ingest loop calls as it goes; the component that is holding the processor
// back is by construction the one whose ring has room, so the backlog always
// drains. A gap is queued as a count rather than as samples: an outage of any
// length costs nothing to hold and is never waited on.
//
// Nothing is dropped. The absolute position of a sample is its ring index plus
// the component's base (station.hpp), so dropping one would move every later
// sample of that component earlier in time -- a silent timing error rather
// than a loss. Backpressure is applied to the source instead.

#include "station.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <span>
#include <vector>

namespace ayzek::pipeline {

// Upper bound on the samples ingest may hold back per station before it stops
// reading from the source. Gaps are excluded, since they are counts.
inline constexpr std::size_t kMaxQueuedSamples = kRingCapacity;

class ComponentWriter {
public:
  explicit ComponentWriter(ComponentStream &cs) noexcept : cs_(&cs) {}

  // Appends samples, writing as many as fit and queueing the rest.
  void write(std::span<const std::int32_t> s) {
    if (q_.empty()) { // steady state: straight into the ring, no allocation
      const std::size_t k = std::min(room(), s.size());
      cs_->ring.push_bulk(s.first(k));
      s = s.subspan(k);
      if (s.empty())
        return;
    }
    q_.push_back(Chunk{0, std::vector<std::int32_t>(s.begin(), s.end()), 0});
    queued_ += q_.back().samples.size();
  }

  // Appends `n` missing samples as kGap. Held as a count, so length is free.
  void write_gap(std::uint64_t n) {
    if (n == 0)
      return;
    if (!q_.empty() && q_.back().samples.empty())
      q_.back().gap += n; // merge with the gap at the tail
    else
      q_.push_back(Chunk{n, {}, 0});
    drain();
  }

  // Writes as much of the backlog as the ring accepts. True if it is empty
  // afterwards, i.e. the component is fully written up to its ingested end.
  bool drain() {
    while (!q_.empty()) {
      const std::size_t r = room();
      if (r == 0)
        return false;
      Chunk &c = q_.front();
      if (c.gap > 0) {
        const auto k = static_cast<std::size_t>(
            std::min<std::uint64_t>({c.gap, r, gap_fill().size()}));
        cs_->ring.push_bulk(
            std::span<const std::int32_t>(gap_fill().data(), k));
        c.gap -= k;
      } else {
        const std::size_t k = std::min(r, c.samples.size() - c.done);
        cs_->ring.push_bulk(
            std::span<const std::int32_t>(c.samples.data() + c.done, k));
        c.done += k;
        queued_ -= k;
      }
      if (c.gap == 0 && c.done == c.samples.size())
        q_.pop_front();
    }
    return true;
  }

  // Samples waiting to be written. Queued gaps are not counted: they cost no
  // memory, and waiting for one to fit is what deadlocked ingest.
  [[nodiscard]] std::size_t queued() const noexcept { return queued_; }

  [[nodiscard]] bool empty() const noexcept { return q_.empty(); }

  // Samples the ring will accept right now.
  [[nodiscard]] std::size_t room() const noexcept {
    return kRingCapacity -
           static_cast<std::size_t>(cs_->ring.written() - cs_->ring.floor());
  }

private:
  // One pending write: either `gap` missing samples or `samples` from
  // `done` on. The two are never mixed in one chunk.
  struct Chunk {
    std::uint64_t gap = 0;
    std::vector<std::int32_t> samples;
    std::size_t done = 0;
  };

  static const std::array<std::int32_t, 1024> &gap_fill() {
    static const std::array<std::int32_t, 1024> f = [] {
      std::array<std::int32_t, 1024> a{};
      a.fill(kGap);
      return a;
    }();
    return f;
  }

  ComponentStream *cs_;
  std::deque<Chunk> q_;
  std::size_t queued_ = 0;
};

} // namespace ayzek::pipeline
