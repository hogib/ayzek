#pragma once

// Ingest stage: miniSEED records -> decoding -> reordering -> component rings.
//
// ReplaySource reads records from a file, and run_ingest feeds them, one
// thread per station. Each record is released when the stream clock reaches
// the time of its last sample, the earliest time a datalogger could transmit
// it. The live input (live.hpp) takes the place of both: one thread receives
// every station's records from a SeedLink server and releases each on
// arrival. Both go through StationFeed, so decoding, reordering and the ring
// writes are the same code.

#include "common.hpp"
#include "mseed.hpp"
#include "reorder.hpp"
#include "station.hpp"
#include "writer.hpp"

#include <array>
#include <atomic>
#include <cstddef>
#include <expected>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace ayzek::pipeline {

class ReplaySource {
public:
  // Reads and indexes every 512-byte record. Throws on an unreadable file.
  explicit ReplaySource(const std::string &path);

  [[nodiscard]] const std::string &station() const noexcept { return station_; }
  [[nodiscard]] double first_time() const noexcept { return first_; }
  [[nodiscard]] double last_time() const noexcept { return last_; }

  struct Record {
    double end_time;
    std::size_t comp;
    std::array<std::byte, 512> bytes;
  };
  // Records sorted by end time, i.e. the delivery order of a stream without
  // late records. Late arrivals present in the archive are not reproduced.
  [[nodiscard]] const std::vector<Record> &records() const noexcept {
    return records_;
  }

private:
  std::string station_;
  double first_ = 0, last_ = 0;
  std::vector<Record> records_;
};

// Positions [t0, t0 + n / 100 Hz) of the three components from the file,
// (n, 3) interleaved counts, bypassing the pipeline (for --assess, after a
// replay). False if any sample in the range is missing.
bool read_window(const ReplaySource &src, double t0, std::size_t n,
                 std::vector<double> &out);

// One station's records into its rings: decoding, one reorderer per component
// and one writer per component. Used by one thread at a time.
class StationFeed {
public:
  explicit StationFeed(Station &st);

  // Decodes record `rec`, whose header is `h`, as component `comp` (Z, N or
  // E), reorders it and writes it. A record that does not decode is counted
  // in Station::undecodable, and nothing else happens.
  std::expected<void, mseed::Error> offer(std::span<const std::byte> rec,
                                          const mseed::Header &h,
                                          std::size_t comp);

  // Writes whatever the rings accept of the samples held back. True if none
  // is left.
  bool drain();
  // Samples held back for want of room in a ring (gaps excluded): all three
  // components', or component `comp`'s.
  [[nodiscard]] std::size_t queued() const;
  [[nodiscard]] std::size_t queued(std::size_t comp) const;

  // Waits while `held()` is true, writing whatever the rings accept meanwhile.
  // Gives up, returning false, when `stop` is set or after a minute in which
  // not one sample reached a ring: a station whose third component never
  // arrives has a processor that can never advance, and the run must end
  // rather than hang. A long wait is not by itself a stall -- across a gap
  // longer than the ring the other components are being written the whole
  // time.
  template <typename Held>
  bool wait_draining(Held held, std::string_view what,
                     const std::atomic<bool> &stop);

private:
  [[nodiscard]] std::uint64_t written() const;
  bool wait_step(const std::atomic<bool> &stop, std::uint64_t &mark,
                 int &idle_ms, std::string_view what);

  Station &st_;
  // max_lateness 0: holes become gaps at once and late records are discarded
  // as stale (docs/DESIGN.md).
  std::array<std::unique_ptr<Reorderer<8, mseed::kMaxDiffs>>, 3> reorder_;
  // A write never waits: what does not fit is queued and retried by `drain`,
  // so a component whose ring is full cannot stop the other two from being
  // fed (writer.hpp).
  std::array<ComponentWriter, 3> writer_;
  std::array<std::int32_t, mseed::kMaxDiffs> buf_{};
};

template <typename Held>
bool StationFeed::wait_draining(Held held, std::string_view what,
                                const std::atomic<bool> &stop) {
  std::uint64_t mark = written();
  for (int idle_ms = 0;;) {
    drain();
    if (!held())
      return true;
    if (!wait_step(stop, mark, idle_ms, what))
      return false;
  }
}

// Feeds all records of `src` into the station's rings, then sets
// st.ingest_done. When a ring is full it waits for the processor (replay only;
// see ingest.cpp).
void run_ingest(const ReplaySource &src, Station &st, const StreamClock &clock,
                const std::atomic<bool> &stop);

} // namespace ayzek::pipeline
