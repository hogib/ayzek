// Tests for ComponentWriter: a write must never wait for room, and the
// absolute position of every sample must survive a gap longer than the ring.
//
// The case that matters is the one that used to deadlock ingest: a gap of
// more than kRingCapacity samples (21.8 minutes at 100 Hz) cannot be written
// in one go, and the processor -- which only advances when all three
// components have data -- cannot make room while ingest waits inside the
// write. Here the consumer is driven by hand, so the property is checked
// without threads: the writer returns at once, the backlog drains as the
// floor rises, and position + base still names the sample it should.

#include "pipeline/writer.hpp"

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <print>
#include <vector>

#define CHECK(cond)                                                            \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::println(stderr, "CHECK failed: {}  ({}:{})", #cond, __FILE__,       \
                   __LINE__);                                                  \
      std::abort();                                                            \
    }                                                                          \
  } while (0)

using namespace ayzek::pipeline;

namespace {

// Each sample carries its own absolute position, so a shift of the stream
// shows up as a wrong value rather than as missing data.
std::vector<std::int32_t> ramp(std::uint64_t from, std::size_t n) {
  std::vector<std::int32_t> v(n);
  for (std::size_t i = 0; i < n; ++i)
    v[i] = static_cast<std::int32_t>(from + i);
  return v;
}

// Stands in for Processor::run: reads everything the ring holds, checks it,
// and raises the floor to release it. Returns the number of samples read.
std::size_t consume(ComponentStream &cs, std::uint64_t &next,
                    std::size_t &gaps_seen) {
  const std::uint64_t end = cs.ring.written();
  std::size_t read = 0;
  std::vector<std::int32_t> buf;
  while (next < end) {
    const auto n = static_cast<std::size_t>(std::min<std::uint64_t>(4096, end - next));
    buf.assign(n, 0);
    auto got = cs.ring.copy_out(next, std::span<std::int32_t>(buf.data(), n));
    CHECK(got.has_value());
    for (std::size_t i = 0; i < n; ++i) {
      if (buf[i] == kGap)
        ++gaps_seen;
      else
        CHECK(buf[i] == static_cast<std::int32_t>(next + i));
    }
    next += n;
    read += n;
  }
  cs.ring.set_floor(next);
  return read;
}

// A gap longer than the ring: queued as a count, written as the floor rises,
// and the samples after it land at their own positions.
void gap_longer_than_the_ring() {
  auto cs = std::make_unique<ComponentStream>();
  ComponentWriter w(*cs);
  cs->base.store(0);

  const auto head = ramp(0, 1000);
  w.write(head);
  CHECK(cs->ring.written() == 1000);

  // Three rings and a bit of missing data. The call must return without
  // waiting for room it can never have.
  const std::uint64_t gap_len = 3 * kRingCapacity + 12'345;
  w.write_gap(gap_len);
  CHECK(!w.empty());          // most of it is still owed
  CHECK(w.queued() == 0);     // but a queued gap costs no memory

  const auto tail = ramp(1000 + gap_len, 1000);
  w.write(tail);

  // Drain against a consumer that only reads what is there. Bounded: each
  // round writes at least one ring's worth unless the backlog is done.
  std::uint64_t next = 0;
  std::size_t gaps_seen = 0, rounds = 0;
  while (!w.drain()) {
    CHECK(++rounds < 1000);
    consume(*cs, next, gaps_seen);
  }
  consume(*cs, next, gaps_seen);

  CHECK(cs->ring.dropped() == 0);
  CHECK(gaps_seen == gap_len);
  CHECK(cs->end() == 1000 + gap_len + 1000);
  CHECK(next == cs->end());
}

// Real data beyond the ring's capacity is held, not dropped: the position of
// everything after it depends on that.
void data_beyond_the_ring_is_held() {
  auto cs = std::make_unique<ComponentStream>();
  ComponentWriter w(*cs);
  cs->base.store(0);

  const std::size_t n = kRingCapacity + 50'000;
  const auto all = ramp(0, n);
  for (std::size_t i = 0; i < n; i += 500)
    w.write(std::span<const std::int32_t>(all.data() + i,
                                          std::min<std::size_t>(500, n - i)));
  CHECK(w.queued() == n - kRingCapacity);

  std::uint64_t next = 0;
  std::size_t gaps_seen = 0, rounds = 0;
  while (!w.drain()) {
    CHECK(++rounds < 1000);
    consume(*cs, next, gaps_seen);
  }
  consume(*cs, next, gaps_seen);

  CHECK(cs->ring.dropped() == 0);
  CHECK(gaps_seen == 0);
  CHECK(w.queued() == 0);
  CHECK(cs->end() == n);
  CHECK(next == n);
}

// Gaps arriving back to back merge into one pending count, so an outage made
// of many missing records is as cheap as one long one.
void consecutive_gaps_merge() {
  auto cs = std::make_unique<ComponentStream>();
  ComponentWriter w(*cs);
  cs->base.store(0);

  w.write(ramp(0, kRingCapacity)); // fills the ring, nothing queued
  CHECK(w.empty());
  for (int i = 0; i < 1000; ++i)
    w.write_gap(1000);
  CHECK(w.queued() == 0);

  std::uint64_t next = 0;
  std::size_t gaps_seen = 0, rounds = 0;
  while (!w.drain()) {
    CHECK(++rounds < 1000);
    consume(*cs, next, gaps_seen);
  }
  consume(*cs, next, gaps_seen);
  CHECK(gaps_seen == 1000 * 1000);
  CHECK(cs->end() == kRingCapacity + 1000 * 1000);
  CHECK(cs->ring.dropped() == 0);
}

} // namespace

int main() {
  gap_longer_than_the_ring();
  std::println("  gap longer than the ring  ok");
  data_beyond_the_ring_is_held();
  std::println("  data held, never dropped  ok");
  consecutive_gaps_merge();
  std::println("  consecutive gaps merge    ok");
  std::println("all writer tests passed");
}
