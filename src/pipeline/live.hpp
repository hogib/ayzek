#pragma once

// Live input: one thread receives every station's records from a SeedLink
// server (seedlink.hpp) and releases each on arrival into its station's rings,
// through the StationFeed the replay uses (ingest.hpp), so decoding,
// reordering and the ring writes are the same code.
//
// The pipeline takes HHZ, HHN and HHE at exactly 100 Hz. Every other record --
// another channel or rate, HH1/HH2 (which would need a rotation that live
// input does not do), a second location code, an encoding the decoder does
// not support -- is counted by reason and reported, the first of each reason
// as it happens. A record that does not decode is counted in
// Station::undecodable.
//
// The reader never waits for a processor. While a component holds more than
// a ring's worth of samples back (kMaxQueuedSamples, writer.hpp), its new
// records are dropped and counted. A dropped record becomes a gap when the
// next one arrives, so no later sample moves in time. A component that was
// out holds nothing back, so its records are taken again at once and its
// outage becomes a gap.
//
// With an archive directory, every record received is appended, as received,
// to DIR/<STATION>.mseed, which the replay reads like any other file.

#include "ingest.hpp"
#include "seedlink.hpp"
#include "station.hpp"

#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace ayzek::pipeline {

struct LiveConfig {
  seedlink::ClientConfig client;
  std::string archive_dir; // "": no archive
  double duration_s = 0;   // 0: until stopped
};

struct LiveStationStats {
  std::uint64_t received = 0; // records the client delivered
  std::uint64_t accepted = 0; // decoded into the rings
  std::uint64_t dropped = 0;  // not offered: too much held back
  std::map<std::string, std::uint64_t> rejected; // reason -> records
  // End times of the HHZ, HHN and HHE records at 100 Hz, dropped ones
  // included, as ReplaySource's first_time and last_time.
  double first_end = std::numeric_limits<double>::infinity();
  double last_end = -std::numeric_limits<double>::infinity();
  // Latency: wall-clock arrival minus the end time of the record's last
  // sample, in 0.1 s bins up to two minutes.
  static constexpr std::size_t kLatencyBins = 1200;
  std::uint64_t latency_n = 0;
  double latency_max = -std::numeric_limits<double>::infinity();
  std::array<std::uint32_t, kLatencyBins + 1> latency_hist{}; // last: beyond
  void add_latency(double seconds);
  // Median latency to the bin; +inf beyond the bins, NaN without records.
  [[nodiscard]] double latency_median() const;
};

class LiveInput {
public:
  // stations[i] receives the records of cfg.client.streams[i]. Opens the
  // archive files; throws if one cannot be opened.
  LiveInput(LiveConfig cfg, std::vector<Station *> stations);

  // Receives until `stop` is set or the duration has passed, writes what is
  // held back and sets every station's ingest_done. `hard_stop` ends the
  // writing too. An error that ends the input early (the server refused a
  // station) is kept in error().
  void run(const std::atomic<bool> &stop, const std::atomic<bool> &hard_stop);

  // One record of station i as the client delivers it, arriving at `arrival`
  // (epoch seconds). Public for the tests.
  void on_record(std::size_t i, seedlink::Record rec, double arrival);
  // Writes what the rings accept of the samples held back.
  void drain();

  [[nodiscard]] const std::string &error() const noexcept { return error_; }
  // The earliest first_end and the latest last_end over the stations; NaN
  // before any record.
  [[nodiscard]] double first_time() const;
  [[nodiscard]] double last_time() const;
  [[nodiscard]] const LiveStationStats &stats(std::size_t i) const {
    return slots_[i]->stats;
  }
  [[nodiscard]] const seedlink::ClientStats &client_stats() const noexcept {
    return client_.stats();
  }
  // The input's section of the end-of-run report.
  void report() const;

private:
  struct Slot {
    explicit Slot(Station &s) : st(&s), feed(s) {}
    Station *st;
    StationFeed feed;
    std::ofstream archive;
    bool archive_failed = false;
    std::optional<std::string> loc; // any-location streams: the one taken
    LiveStationStats stats;
  };
  void reject(Slot &s, const std::string &reason);

  LiveConfig cfg_;
  std::vector<std::unique_ptr<Slot>> slots_;
  seedlink::Client client_;
  std::string error_;
};

} // namespace ayzek::pipeline
