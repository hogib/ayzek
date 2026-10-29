#include "live.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <format>
#include <stdexcept>

namespace ayzek::pipeline {

namespace {

double wall_now() {
  return std::chrono::duration<double>(
             std::chrono::system_clock::now().time_since_epoch())
      .count();
}

} // namespace

void LiveStationStats::add_latency(double seconds) {
  ++latency_n;
  latency_max = std::max(latency_max, seconds);
  const double bin = std::floor(std::max(seconds, 0.0) * 10.0);
  ++latency_hist[bin < kLatencyBins ? static_cast<std::size_t>(bin)
                                    : kLatencyBins];
}

double LiveStationStats::latency_median() const {
  if (latency_n == 0)
    return std::numeric_limits<double>::quiet_NaN();
  std::uint64_t below = 0;
  for (std::size_t b = 0; b < kLatencyBins; ++b) {
    below += latency_hist[b];
    if (2 * below >= latency_n)
      return (static_cast<double>(b) + 0.5) / 10.0;
  }
  return std::numeric_limits<double>::infinity();
}

LiveInput::LiveInput(LiveConfig cfg, std::vector<Station *> stations)
    : cfg_(std::move(cfg)),
      client_(
          cfg_.client,
          [this](std::size_t i, seedlink::Record rec) {
            on_record(i, rec, wall_now());
          },
          [this] { drain(); }) {
  if (stations.size() != cfg_.client.streams.size())
    throw std::invalid_argument("LiveInput: one station per stream");
  if (!cfg_.archive_dir.empty())
    std::filesystem::create_directories(cfg_.archive_dir);
  for (Station *st : stations) {
    auto s = std::make_unique<Slot>(*st);
    if (!cfg_.archive_dir.empty()) {
      const auto path =
          std::filesystem::path(cfg_.archive_dir) / (st->code + ".mseed");
      s->archive.open(path, std::ios::binary | std::ios::app);
      if (!s->archive)
        throw std::runtime_error("cannot write " + path.string());
    }
    slots_.push_back(std::move(s));
  }
}

void LiveInput::reject(Slot &s, const std::string &reason) {
  if (s.stats.rejected[reason]++ == 0)
    Log::get().line("warn", "33", "{}: records rejected: {}", s.st->code,
                    reason);
}

void LiveInput::on_record(std::size_t i, seedlink::Record rec,
                          double arrival) {
  Slot &s = *slots_[i];
  ++s.stats.received;
  if (s.archive.is_open() && !s.archive_failed) {
    s.archive.write(reinterpret_cast<const char *>(rec.data()),
                    static_cast<std::streamsize>(rec.size()));
    s.archive.flush();
    if (!s.archive) {
      s.archive_failed = true;
      Log::get().line("warn", "33",
                      "{}: cannot write to the archive; archiving stopped",
                      s.st->code);
    }
  }

  const auto h = mseed::parse_header(rec);
  if (!h) {
    reject(s, std::string(mseed::to_string(h.error())));
    return;
  }
  const auto &spec = cfg_.client.streams[i];
  const std::string loc(h->loc());
  if (spec.loc && loc != *spec.loc) {
    reject(s, std::format("location '{}', not the '{}' requested", loc,
                          *spec.loc));
    return;
  }
  const std::string_view cha = h->cha();
  std::size_t comp = 3;
  if (cha == "HHZ")
    comp = Z;
  else if (cha == "HHN")
    comp = N;
  else if (cha == "HHE")
    comp = E;
  if (comp == 3) {
    reject(s, cha == "HH1" || cha == "HH2"
                  ? std::format("{}: live input takes HHZ, HHN and HHE only "
                                "and does not rotate HH1/HH2",
                                cha)
                  : std::format("{}: live input takes HHZ, HHN and HHE only",
                                cha));
    return;
  }
  if (h->sample_rate != kFs) {
    reject(s, std::format("{} at {:g} Hz: live input takes 100 Hz only", cha,
                          h->sample_rate));
    return;
  }
  // Any location: the first HHZ, HHN or HHE record at 100 Hz fixes it, so a
  // channel the pipeline does not take cannot.
  if (!spec.loc && !s.loc)
    s.loc = loc;
  if (!spec.loc && loc != *s.loc) {
    reject(s, std::format("location '{}': HH? is taken from location '{}', "
                          "the first one received (NET.STA.LOC chooses)",
                          loc, *s.loc));
    return;
  }
  const double end = static_cast<double>(h->start_ns) / 1e9 +
                     h->num_samples / h->sample_rate;
  s.stats.first_end = std::min(s.stats.first_end, end);
  s.stats.last_end = std::max(s.stats.last_end, end);
  s.stats.add_latency(arrival - end);

  // The processor is more than a ring's worth behind this component: drop
  // rather than wait (live.hpp). Only this component's backlog counts: one
  // that has been out holds nothing back, and its records must reach the
  // reorderer, which makes the outage a gap so the processor can advance.
  if (s.feed.queued(comp) > kMaxQueuedSamples) {
    if (s.stats.dropped++ % 100 == 0)
      Log::get().line("warn", "33",
                      "{}: processor more than a ring behind; {} record(s) "
                      "dropped so far",
                      s.st->code, s.stats.dropped);
    return;
  }
  if (const auto ok = s.feed.offer(rec, *h, comp); !ok) {
    if (s.st->undecodable.load(std::memory_order_relaxed) == 1)
      Log::get().line("warn", "33", "{}: a {} record did not decode ({})",
                      s.st->code, cha, mseed::to_string(ok.error()));
    return;
  }
  ++s.stats.accepted;
}

void LiveInput::drain() {
  for (auto &s : slots_)
    s->feed.drain();
}

void LiveInput::run(const std::atomic<bool> &stop,
                    const std::atomic<bool> &hard_stop) {
  const auto t0 = std::chrono::steady_clock::now();
  const auto should_stop = [&] {
    return stop.load(std::memory_order_relaxed) ||
           hard_stop.load(std::memory_order_relaxed) ||
           (cfg_.duration_s > 0 &&
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
                    .count() >= cfg_.duration_s);
  };
  try {
    client_.run(should_stop);
  } catch (const std::exception &e) {
    error_ = e.what();
  }
  // As at the end of a replay (run_ingest): the stream is not finished until
  // everything held back has been written.
  for (auto &s : slots_)
    if (!s->feed.wait_draining([&] { return !s->feed.drain(); },
                               "flushing the last samples", hard_stop) &&
        s->feed.queued() > 0)
      Log::get().line("warn", "33", "{}: {} samples at the end not written",
                      s->st->code, s->feed.queued());
  for (auto &s : slots_)
    s->st->ingest_done.store(true, std::memory_order_release);
}

double LiveInput::first_time() const {
  double t = std::numeric_limits<double>::infinity();
  for (const auto &s : slots_)
    t = std::min(t, s->stats.first_end);
  return std::isinf(t) ? std::numeric_limits<double>::quiet_NaN() : t;
}

double LiveInput::last_time() const {
  double t = -std::numeric_limits<double>::infinity();
  for (const auto &s : slots_)
    t = std::max(t, s->stats.last_end);
  return std::isinf(t) ? std::numeric_limits<double>::quiet_NaN() : t;
}

void LiveInput::report() const {
  auto &log = Log::get();
  const auto &c = client_.stats();
  log.plain(false, "");
  log.plain(true, "SeedLink input from {}:{}", cfg_.client.server.host,
            cfg_.client.server.port);
  log.plain(false,
            "  {} connections, {} sessions; {} packets, {} duplicates "
            "discarded, {} of no requested station",
            c.connections, c.sessions, c.packets, c.duplicates, c.unknown);
  log.plain(false,
            "  {} keepalives sent, {} INFO packets and {} ERROR replies "
            "received{}",
            c.keepalives, c.info_packets, c.server_errors,
            c.seq_restarts ? std::format(", {} packets after a sequence restart",
                                         c.seq_restarts)
                           : std::string());
  log.plain(false, "  {:<5} {:>8} {:>8} {:>8} {:>8} {:>8}   {:>20}", "",
            "received", "accepted", "rejected", "dropped", "bad data",
            "latency median (max)");
  for (const auto &s : slots_) {
    const auto &st = s->stats;
    std::uint64_t rejected = 0;
    for (const auto &[_, n] : st.rejected)
      rejected += n;
    const double med = st.latency_median();
    const std::string latency =
        st.latency_n == 0 ? std::string("-")
        : std::isinf(med)
            ? std::format("> {:g} s",
                          LiveStationStats::kLatencyBins / 10.0)
            : std::format("{:.1f} s ({:.1f} s)", med, st.latency_max);
    log.plain(false, "  {:<5} {:>8} {:>8} {:>8} {:>8} {:>8}   {:>20}",
              s->st->code, st.received, st.accepted, rejected, st.dropped,
              s->st->undecodable.load(), latency);
  }
  for (const auto &s : slots_)
    for (const auto &[reason, n] : s->stats.rejected)
      log.plain(false, "  {} rejected {} records: {}", s->st->code, n, reason);
  if (!cfg_.archive_dir.empty())
    log.plain(false, "  every record received appended to {}/STATION.mseed",
              std::filesystem::absolute(cfg_.archive_dir).string());
}

} // namespace ayzek::pipeline
