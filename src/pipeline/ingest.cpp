#include "ingest.hpp"

#include "mseed.hpp"
#include "reorder.hpp"
#include "writer.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <thread>

namespace ayzek::pipeline {

namespace {

std::size_t component_of(std::string_view cha) {
  switch (cha.empty() ? '?' : cha.back()) {
  case 'Z':
    return Z;
  case 'N':
  case '1':
    return N;
  case 'E':
  case '2':
    return E;
  default:
    return 3;
  }
}

} // namespace

ReplaySource::ReplaySource(const std::string &path) {
  std::ifstream f(path, std::ios::binary);
  if (!f)
    throw std::runtime_error("cannot open " + path);
  std::array<std::byte, 512> raw{};
  std::map<std::string, std::size_t> skipped; // reason -> records
  while (f.read(reinterpret_cast<char *>(raw.data()), raw.size())) {
    auto h = mseed::parse_header(raw);
    if (!h) {
      ++skipped[std::string(mseed::to_string(h.error()))];
      continue;
    }
    const std::size_t comp = component_of(h->cha());
    if (comp > E || h->cha().substr(0, 2) != "HH" || h->sample_rate != 100.0) {
      ++skipped["not a 100 Hz HH? channel"];
      continue;
    }
    if (station_.empty())
      station_ = std::string(h->sta());
    const double start = static_cast<double>(h->start_ns) / 1e9;
    const double end = start + h->num_samples / h->sample_rate;
    records_.push_back({end, comp, raw});
  }
  if (records_.empty())
    throw std::runtime_error(path + ": no usable 100 Hz HH? records");
  std::stable_sort(
      records_.begin(), records_.end(),
      [](const Record &a, const Record &b) { return a.end_time < b.end_time; });
  first_ = records_.front().end_time;
  last_ = records_.back().end_time;
  for (const auto &[reason, count] : skipped)
    Log::get().line("replay", "33", "{}: skipped {} records ({})", station_,
                    count, reason);
}

void run_ingest(const ReplaySource &src, Station &st, const StreamClock &clock,
                const std::atomic<bool> &stop) {
  // One reorderer per component with max_lateness 0: holes become gaps at once
  // and late records are discarded as stale (docs/DESIGN.md).
  std::array<std::unique_ptr<Reorderer<8, mseed::kMaxDiffs>>, 3> reorder;
  for (auto &r : reorder)
    r = std::make_unique<Reorderer<8, mseed::kMaxDiffs>>(0);
  std::array<std::int32_t, mseed::kMaxDiffs> buf{};

  // One writer per component. A write never waits: what does not fit is
  // queued and retried by `drain`, so a component whose ring is full cannot
  // stop the other two from being fed (writer.hpp).
  std::array<ComponentWriter, 3> writer{ComponentWriter(st.comp[Z]),
                                        ComponentWriter(st.comp[N]),
                                        ComponentWriter(st.comp[E])};
  auto drain_all = [&] {
    bool empty = true;
    for (auto &w : writer)
      empty = w.drain() && empty;
    return empty;
  };
  auto queued = [&] {
    std::size_t n = 0;
    for (const auto &w : writer)
      n += w.queued();
    return n;
  };
  auto written = [&] {
    std::uint64_t n = 0;
    for (const auto &c : st.comp)
      n += c.ring.written();
    return n;
  };
  // Waits while `held` is true, writing whatever the rings accept meanwhile.
  // Gives up, returning false, after a minute in which not one sample reached
  // a ring: a station whose third component never arrives has a processor
  // that can never advance, and the run must end rather than hang. A long
  // wait is not by itself a stall -- across a gap longer than the ring the
  // other components are being written the whole time.
  auto wait_draining = [&](auto held, std::string_view what) {
    std::uint64_t mark = written();
    for (int idle_ms = 0;;) {
      drain_all();
      if (!held())
        return true;
      if (stop.load(std::memory_order_relaxed))
        return false;
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
      if (const auto now = written(); now != mark) {
        mark = now;
        idle_ms = 0;
        continue;
      }
      if (++idle_ms >= 60'000) {
        Log::get().line("warn", "33",
                        "{}: {} for 60 s with no progress, {} samples held",
                        st.code, what, queued());
        return false;
      }
    }
  };
  bool hold_source = true; // cleared for the rest of the run if it gives up

  for (const auto &rec : src.records()) {
    if (stop.load(std::memory_order_relaxed))
      break;
    while (rec.end_time > clock.now() &&
           !stop.load(std::memory_order_relaxed)) {
      drain_all();
      std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    // Backpressure on the source rather than inside a write: the replay is
    // held back while the processor is more than a ring's worth behind. Live
    // data would be dropped here by the socket instead, which is the one
    // place where dropping does not shift later samples in time. If the
    // processor turns out not to be advancing at all, the rest of the replay
    // is buffered instead, which costs memory but keeps every position right.
    if (hold_source && queued() > kMaxQueuedSamples)
      hold_source = wait_draining([&] { return queued() > kMaxQueuedSamples; },
                                  "processor behind, holding the source");

    auto h = mseed::parse_header(rec.bytes);
    auto n = mseed::decode(rec.bytes, *h, buf);
    if (!n)
      continue;
    const auto pos = static_cast<std::uint64_t>(h->start_ns) / 10'000'000;
    auto &cs = st.comp[rec.comp];
    auto &w = writer[rec.comp];

    auto sink = [&](std::uint64_t p, std::span<const std::int32_t> s) {
      if (cs.base.load(std::memory_order_relaxed) == kUnset)
        cs.base.store(p, std::memory_order_release);
      w.write(s);
    };
    auto gap = [&](std::uint64_t a, std::uint64_t b) {
      cs.gaps.fetch_add(1, std::memory_order_relaxed);
      cs.gap_samples.fetch_add(b - a, std::memory_order_relaxed);
      w.write_gap(b - a);
    };
    reorder[rec.comp]->offer(pos, std::span<const std::int32_t>(buf.data(), *n),
                             sink, gap);
    cs.stale.store(reorder[rec.comp]->stats().stale, std::memory_order_relaxed);
    st.records.fetch_add(1, std::memory_order_relaxed);
    drain_all();
  }
  // The stream is not finished until everything held has been written, or the
  // processor would stop at a truncated one. If it cannot be written the tail
  // is given up, which loses the end of the replay but leaves the position of
  // every sample before it right.
  if (!wait_draining([&] { return !drain_all(); }, "flushing the last samples") &&
      queued() > 0)
    Log::get().line("warn", "33", "{}: {} samples at the end not written",
                    st.code, queued());
  st.ingest_done.store(true, std::memory_order_release);
}

} // namespace ayzek::pipeline
