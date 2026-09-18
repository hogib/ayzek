// ayzek: earthquake detection, phase picking, magnitude estimation and location
// on miniSEED data replayed in real time, one file per station.
//
//   ayzek [options] STATION.mseed...
//
// See docs/IMPLEMENTATION.md.

#include "pipeline/common.hpp"
#include "pipeline/ingest.hpp"
#include "pipeline/network.hpp"
#include "pipeline/processor.hpp"
#include "pipeline/recording.hpp"
#include "pipeline/station.hpp"
#include "simd.hpp"

#include <algorithm>
#include <chrono>
#include <csignal>
#include <cstdio>
#include <format>
#include <fstream>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <ranges>
#include <sstream>
#include <thread>
#include <tuple>
#include <unistd.h>

using namespace ayzek;
using namespace ayzek::pipeline;

namespace {

std::atomic<bool> g_stop{false};

const char *kUsage = R"(usage: ayzek [options] STATION.mseed...

  --models DIR          weights directory (default: models)
  --speed X             replay speed; 1 = real time, 0 = as fast as possible (default 1)
  --from TIME           start detecting at this UTC time, e.g. 2025-11-10T18:19:00
  --threshold P         trigger when N consecutive windows reach P (default 0.8) ...
  --trigger-windows N   ... (default 8, which adds 3.5 s)
  --instant-threshold P or when one window reaches P (default 0.9; > 1 disables)
  --release P           probability below which a trigger resets (default 0.3)
  --detector KIND       model (default) or stalta, the reference STA/LTA trigger
  --no-anchor           do not date model triggers by the STA/LTA onset
  --anchor-on R         STA/LTA ratio taken as the onset when anchoring (default 3)
  --require-onset       a detector trigger also needs an STA/LTA onset
  --anchor-picker       place the picker window from the anchored P too
  --anchor-magnitude    place the early magnitude window from the anchored P too
  --mag-lead S          magnitude window start, seconds before P (default 2)
  --sta S, --lta S      STA/LTA averaging lengths in seconds (default 1, 30)
  --stalta-on R         STA/LTA ratio that triggers (default 8)
  --stalta-off R        STA/LTA ratio below which a trigger resets (default 1.5)
  --stalta-band LO,HI   STA/LTA pass band in Hz (default 2,20)
  --stalta-3c           STA/LTA on the energy of all three components (default: vertical)
  --step N              samples between detector windows (default 50 = 0.5 s)
  --min-stations N      detections needed to declare an event (default 2)
  --no-pick             detector only
  --no-magnitude        skip the magnitude regressor
  --catalog CSV         AFAD catalogue export to score events against
  --scores DIR          write every window's probability to DIR/STATION.csv
  --scores-in DIR       use probabilities from an earlier --scores run instead of the detector
  --record FILE         write all station outputs for tools/network_subsets
  --site NAME,LAT,LON   report the S-wave warning time at this place (repeatable)
  --verbose             print every station detection, pick and magnitude estimate
  --no-color
)";

double parse_time(const std::string &s) {
  const auto t = parse_utc(s, "ymdhms");
  if (!t)
    throw std::runtime_error("bad time: " + s + " (want YYYY-MM-DDTHH:MM:SS)");
  return *t;
}

std::map<std::string, StationInfo> load_stations(const std::string &path) {
  std::map<std::string, StationInfo> out;
  std::ifstream f(path);
  if (!f)
    throw std::runtime_error("cannot open " + path);
  std::string line;
  std::getline(f, line);
  while (std::getline(f, line)) {
    std::stringstream ss(line);
    std::string code, lat, lon;
    std::getline(ss, code, ',');
    std::getline(ss, lat, ',');
    std::getline(ss, lon, ',');
    // Some station codes occur more than once in the AFAD table; the first is
    // used.
    if (!code.empty() && !out.contains(code))
      out[code] = {std::stod(lat), std::stod(lon)};
  }
  return out;
}

struct Percentiles {
  double mean, p50, p99, max;
};
Percentiles percentiles(std::vector<float> v) {
  if (v.empty())
    return {0, 0, 0, 0};
  std::ranges::sort(v);
  const auto at = [&](double q) {
    return static_cast<double>(
        v[static_cast<std::size_t>(q * static_cast<double>(v.size() - 1))]);
  };
  return {std::accumulate(v.begin(), v.end(), 0.0) /
              static_cast<double>(v.size()),
          at(0.5), at(0.99), static_cast<double>(v.back())};
}

} // namespace

int main(int argc, char **argv) try {
  std::string models = "models", scores_dir, scores_in_dir, catalog_path,
              record_path;
  double speed = 1.0;
  ProcessorConfig pcfg;
  NetworkConfig ncfg;
  std::vector<std::string> files;

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc)
        throw std::runtime_error(a + " needs a value");
      return argv[++i];
    };
    if (a == "--models")
      models = next();
    else if (a == "--speed")
      speed = std::stod(next());
    else if (a == "--from")
      pcfg.from = parse_time(next());
    else if (a == "--threshold")
      pcfg.threshold = std::stof(next());
    else if (a == "--trigger-windows")
      pcfg.trigger_windows = std::stoul(next());
    else if (a == "--instant-threshold")
      pcfg.instant_threshold = std::stof(next());
    else if (a == "--release")
      pcfg.release = std::stof(next());
    else if (a == "--detector") {
      const auto kind = next();
      if (kind == "model")
        pcfg.detector = DetectorKind::Model;
      else if (kind == "stalta")
        pcfg.detector = DetectorKind::StaLta;
      else
        throw std::runtime_error("--detector: model or stalta");
    } else if (a == "--sta")
      pcfg.stalta.sta_seconds = std::stod(next());
    else if (a == "--lta")
      pcfg.stalta.lta_seconds = std::stod(next());
    else if (a == "--stalta-on")
      pcfg.stalta_on = std::stod(next());
    else if (a == "--stalta-off")
      pcfg.stalta_off = std::stod(next());
    else if (a == "--stalta-band") {
      const auto band = next();
      const auto comma = band.find(',');
      if (comma == std::string::npos)
        throw std::runtime_error("--stalta-band: LO,HI");
      pcfg.stalta.f_lo = std::stod(band.substr(0, comma));
      pcfg.stalta.f_hi = std::stod(band.substr(comma + 1));
    } else if (a == "--stalta-3c")
      pcfg.stalta.three_component = true;
    else if (a == "--no-anchor")
      pcfg.anchor = false;
    else if (a == "--anchor-on")
      pcfg.anchor_on = std::stod(next());
    else if (a == "--anchor-picker")
      pcfg.anchor_picker = true;
    else if (a == "--anchor-magnitude")
      pcfg.anchor_magnitude = true;
    else if (a == "--require-onset")
      pcfg.require_onset = true;
    else if (a == "--mag-lead")
      pcfg.magnitude_lead = std::stod(next());
    else if (a == "--step")
      pcfg.step = std::stoul(next());
    else if (a == "--min-stations")
      ncfg.min_stations = std::stoul(next());
    else if (a == "--no-pick")
      pcfg.pick = false;
    else if (a == "--no-magnitude")
      pcfg.magnitude = false;
    else if (a == "--scores")
      scores_dir = next();
    else if (a == "--scores-in")
      scores_in_dir = next();
    else if (a == "--no-color")
      Log::get().color = false;
    else if (a == "--verbose")
      ncfg.verbose = true;
    else if (a == "--site") {
      std::stringstream ss(next());
      Site site;
      std::string lat, lon;
      std::getline(ss, site.name, ',');
      std::getline(ss, lat, ',');
      std::getline(ss, lon, ',');
      site.lat = std::stod(lat);
      site.lon = std::stod(lon);
      ncfg.sites.push_back(site);
    } else if (a == "--catalog")
      catalog_path = next();
    else if (a == "--record")
      record_path = next();
    else if (a == "-h" || a == "--help") {
      std::print("{}", kUsage);
      return 0;
    } else if (a.starts_with("--")) {
      throw std::runtime_error("unknown option " + a);
    } else {
      files.push_back(a);
    }
  }
  if (files.empty()) {
    std::print(stderr, "{}", kUsage);
    return 2;
  }
  if (pcfg.detector == DetectorKind::StaLta && !scores_in_dir.empty())
    throw std::runtime_error("--scores-in applies to the model detector only");
  if (!isatty(STDOUT_FILENO))
    Log::get().color = false;
  std::signal(SIGINT, [](int) { g_stop = true; });

  // --- load
  // ---------------------------------------------------------------------
  std::vector<Weights> detector;
  for (int seed : {42, 43, 44})
    detector.push_back(
        Weights::load(std::format("{}/detector_s{}.ayzw", models, seed)));
  const Weights picker = Weights::load(models + "/spicker.ayzw");
  std::vector<Weights> magnitude;
  if (pcfg.magnitude) {
    for (int p = 0; p < 3; ++p) {
      const auto path = std::format("{}/magnitude_p{}.ayzw", models, p);
      if (std::ifstream(path).good())
        magnitude.push_back(Weights::load(path));
    }
    if (magnitude.empty()) {
      Log::get().line("warn", "33",
                      "no magnitude models in {} (tools/export_magnitude.py); "
                      "magnitude off",
                      models);
      pcfg.magnitude = false;
    }
  }
  const auto bp = dsp::Bandpass::load(Weights::load(models + "/bandpass.ayzw"));
  auto coords = load_stations(models + "/stations.csv");

  std::vector<std::unique_ptr<ReplaySource>> sources;
  std::vector<std::unique_ptr<Station>> stations;
  std::map<std::string, StationInfo> used;
  double t_first = 1e18, t_last = 0;
  for (const auto &f : files) {
    auto src = std::make_unique<ReplaySource>(f);
    auto st = std::make_unique<Station>();
    st->code = src->station();
    if (auto it = coords.find(st->code); it != coords.end()) {
      st->lat = it->second.lat;
      st->lon = it->second.lon;
      used[st->code] = it->second;
    } else {
      Log::get().line("warn", "33",
                      "{}: no coordinates, excluded from location", st->code);
    }
    t_first = std::min(t_first, src->first_time());
    t_last = std::max(t_last, src->last_time());
    sources.push_back(std::move(src));
    stations.push_back(std::move(st));
  }

  const double start =
      pcfg.from > 0 ? pcfg.from - 70.0
                    : t_first; // 70 s covers a picker window before `from`
  const std::string detector_name =
      pcfg.detector == DetectorKind::Model
          ? std::string(pcfg.anchor ? "3-seed detector + STA/LTA anchor"
                                    : "3-seed detector")
          : std::format(
                "STA/LTA {:g}/{:g} s on {:g} off {:g}, {:g}-{:g} Hz, {}",
                pcfg.stalta.sta_seconds, pcfg.stalta.lta_seconds,
                pcfg.stalta_on, pcfg.stalta_off, pcfg.stalta.f_lo,
                pcfg.stalta.f_hi, pcfg.stalta.three_component ? "Z+N+E" : "Z");
  Log::get().line(
      "ayzek", "1",
      "{} stations, {} backend, {}{}{}, replay {} to {} UTC at {}",
      stations.size(), simd::kBackend, detector_name,
      pcfg.pick ? " + picker" : "",
      pcfg.magnitude ? std::format(" + {}-model magnitude", magnitude.size())
                     : "",
      ymd_hms(t_first), hms(t_last),
      speed > 0 ? std::format("{:g}x", speed) : std::string("full speed"));
  if (!catalog_path.empty()) {
    // Catalogue events whose waves can arrive within the replayed data.
    ncfg.catalog = load_afad_catalog(
        catalog_path, std::max(t_first, pcfg.from) - 30, t_last - 30);
    Log::get().line("ayzek", "1", "catalogue: {} events in the replay span",
                    ncfg.catalog.size());
  }

  // --- run
  // ----------------------------------------------------------------------
  Bus bus;
  StreamClock clock(start, speed);
  std::vector<std::unique_ptr<Processor>> procs;
  for (auto &st : stations) {
    const std::string path =
        scores_dir.empty() ? "" : scores_dir + "/" + st->code + ".csv";
    const std::string in_path =
        scores_in_dir.empty() ? "" : scores_in_dir + "/" + st->code + ".csv";
    procs.push_back(std::make_unique<Processor>(
        *st, detector, picker, magnitude, bp, pcfg, bus, path, in_path));
  }
  std::vector<std::jthread> threads;
  for (std::size_t i = 0; i < stations.size(); ++i) {
    threads.emplace_back(
        [&, i] { run_ingest(*sources[i], *stations[i], clock, g_stop); });
    threads.emplace_back([&, i] { procs[i]->run(g_stop); });
  }

  // Messages are held in a priority queue and processed in declared_at order,
  // for times strictly before the minimum Progress time over the stations still
  // running: a station may still send a message declared exactly at the time it
  // promised. Messages with equal times are ordered by station code, then by
  // the order the station sent them, so the order does not depend on thread
  // timing.
  Network net(used, ncfg);
  std::unique_ptr<Recorder> recorder;
  if (!record_path.empty())
    recorder = std::make_unique<Recorder>(record_path, used, t_first, t_last);
  std::map<std::string, double> watermark;
  for (auto &st : stations)
    watermark[st->code] = 0;
  struct Held {
    double t;
    std::string station;
    std::uint64_t seq;
    Message m;
  };
  auto later = [](const Held &a, const Held &b) {
    return std::tie(a.t, a.station, a.seq) > std::tie(b.t, b.station, b.seq);
  };
  std::priority_queue<Held, std::vector<Held>, decltype(later)> held(later);
  std::map<std::string, std::uint64_t> sent;

  // Magnitude at the picked P, for picks of declared events only. The pick
  // carries the window and the station's noise baseline, so the estimate
  // equals one made in the station thread. It runs when the pick is released,
  // in stream-time order, so the decision does not depend on thread timing.
  // With --record every pick is estimated, for tools/network_subsets.
  std::unique_ptr<MagnitudeEstimator> pick_magnitude;
  if (pcfg.magnitude)
    pick_magnitude = std::make_unique<MagnitudeEstimator>(magnitude, bp);
  struct PickMagnitudeStats {
    std::size_t estimated = 0, skipped = 0;
    std::vector<float> ms;
  };
  std::map<std::string, PickMagnitudeStats> pick_stats;

  auto release = [&] {
    const double safe = watermark.empty()
                            ? 1e18
                            : std::ranges::min(watermark | std::views::values);
    while (!held.empty() && held.top().t < safe) {
      const Message m = held.top().m;
      held.pop();
      if (auto *d = std::get_if<Detection>(&m))
        net.on(*d);
      else if (auto *p = std::get_if<Pick>(&m))
        net.on(*p);
      else if (auto *g = std::get_if<MagnitudeEstimate>(&m))
        net.on(*g);
      if (recorder)
        recorder->write(m);

      const auto *p = std::get_if<Pick>(&m);
      if (!p || !p->magnitude_window || !pick_magnitude)
        continue;
      auto &ps = pick_stats[p->station];
      if (!recorder && !net.declared(p->station, p->trigger_window)) {
        ++ps.skipped;
        continue;
      }
      const auto &w = *p->magnitude_window;
      const auto t0 = std::chrono::steady_clock::now();
      const float mag = pick_magnitude->estimate(w.raw, w.noise);
      const double ms = std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - t0)
                            .count();
      ++ps.estimated;
      ps.ms.push_back(static_cast<float>(ms));
      const MagnitudeEstimate g{p->station,
                                p->trigger_window,
                                true,
                                w.start,
                                mag,
                                w.noise.ready ? w.noise.windows : 0,
                                p->declared_at,
                                ms};
      net.on(g);
      if (recorder)
        recorder->write(g);
    }
  };
  while (!watermark.empty()) {
    auto m = bus.receive(std::chrono::milliseconds(100));
    if (!m)
      continue;
    if (auto *pr = std::get_if<Progress>(&*m)) {
      if (auto it = watermark.find(pr->station); it != watermark.end())
        it->second = std::max(it->second, pr->until);
    } else if (auto *sd = std::get_if<StationDone>(&*m)) {
      watermark.erase(sd->station);
    } else {
      const std::string station = station_of(*m);
      held.push(Held{declared_at(*m), station, sent[station]++, std::move(*m)});
    }
    release();
  }
  release();
  threads.clear();
  const double wall = clock.wall_seconds();

  // --- report
  // -----------------------------------------------------------------
  const double stream_seconds = t_last - std::max(start, t_first);
  net.report(std::max(start, t_first), t_last);

  auto &log = Log::get();
  log.plain(false, "");
  log.plain(true, "Station processing ({} backend)", simd::kBackend);
  log.plain(false,
            "  {:<5} {:>8} {:>9} {:>8} {:>6} {:>17}   {:>16} {:>9} {:>14}", "",
            "windows", "with gap", "triggers", "picks", "magnitudes early/P",
            "ms/window (p99)", "ms/pick", "ms/magnitude");
  std::uint64_t windows = 0, triggers = 0, anchored = 0;
  std::size_t pick_estimated = 0, pick_skipped = 0;
  for (std::size_t i = 0; i < stations.size(); ++i) {
    const auto &s = procs[i]->stats();
    const auto &ps = pick_stats[stations[i]->code];
    pick_estimated += ps.estimated;
    pick_skipped += ps.skipped;
    const auto w = percentiles(s.window_ms);
    const auto pk = percentiles(s.pick_ms);
    std::vector<float> magnitude_ms = s.magnitude_ms;
    magnitude_ms.insert(magnitude_ms.end(), ps.ms.begin(), ps.ms.end());
    const auto mg = percentiles(magnitude_ms);
    std::uint64_t stale = 0, dropped = 0;
    for (auto &c : stations[i]->comp) {
      stale += c.stale.load();
      dropped += c.ring.dropped();
    }
    if (dropped)
      Log::get().line("warn", "33", "{}: {} samples dropped by the ring",
                      stations[i]->code, dropped);
    windows += s.windows;
    triggers += s.detections;
    anchored += s.anchored;
    if (stale)
      log.line("warn", "33", "{}: {} late records discarded", stations[i]->code,
               stale);
    log.plain(false,
              "  {:<5} {:>8} {:>9} {:>8} {:>6} {:>17}   {:>8.3f} ({:>5.3f}) "
              "{:>9.0f} {:>14.0f}",
              stations[i]->code, s.windows, s.gap_windows, s.detections,
              s.picks, std::format("{}/{}", s.magnitudes, ps.estimated), w.mean,
              w.p99, pk.mean, mg.mean);
  }
  if (pcfg.detector == DetectorKind::Model && pcfg.anchor)
    log.plain(false,
              "  P time from an STA/LTA onset for {} of {} triggers ({} from "
              "the window start)",
              anchored, triggers, triggers - anchored);
  if (pcfg.magnitude)
    log.plain(false,
              "  magnitude at the picked P: {} estimated, {} not (trigger not "
              "part of a declared event)",
              pick_estimated, pick_skipped);
  log.plain(false, "  {} windows in {:.1f} s: {:.0f}x real time per station",
            windows, wall, stream_seconds / wall);
  return 0;
} catch (const std::exception &e) {
  std::println(stderr, "ayzek: {}", e.what());
  return 1;
}
