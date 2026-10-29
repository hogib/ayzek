// ayzek: earthquake detection, phase picking, magnitude estimation and location
// on miniSEED data replayed in real time, one file per station, or received
// live from a SeedLink server.
//
//   ayzek [options] STATION.mseed...
//   ayzek [options] --seedlink HOST:PORT --streams NET.STA[,NET.STA...]
//
// See docs/IMPLEMENTATION.md.

#include "pipeline/common.hpp"
#include "pipeline/ingest.hpp"
#include "pipeline/live.hpp"
#include "pipeline/network.hpp"
#include "pipeline/processor.hpp"
#include "pipeline/recording.hpp"
#include "pipeline/station.hpp"
#include "seedlink.hpp"
#include "simd.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <ranges>
#include <set>
#include <sstream>
#include <thread>
#include <tuple>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

namespace {
// True when stdout is a terminal, so colour is worth emitting.
bool stdout_is_tty() {
#ifdef _WIN32
  return _isatty(_fileno(stdout)) != 0;
#else
  return isatty(STDOUT_FILENO) != 0;
#endif
}
} // namespace

using namespace ayzek;
using namespace ayzek::pipeline;

namespace {

std::atomic<bool> g_stop{false};
// Live input: the first SIGINT or SIGTERM stops the input, and the run ends
// once what was received has been processed; a second stops at once.
std::atomic<bool> g_input_stop{false};

const char *kUsage = R"(usage: ayzek [options] STATION.mseed...
       ayzek [options] --seedlink HOST:PORT --streams NET.STA[,NET.STA...]

  --models DIR          weights directory (default: models)
  --stations FILE       station coordinates (code,lat,lon[,...]) in place of
                        MODELS/stations.csv, for a replay of another network
  --pd-gains FILE       the replay's own Pd gains (pd_gain.csv's columns) in
                        place of MODELS/pd_gain.csv; KO's station terms and
                        noise levels are then not used
  --speed X             replay speed; 1 = real time, 0 = as fast as possible (default 1)
  --from TIME           start detecting at this UTC time, e.g. 2025-11-10T18:19:00
  --threshold P         trigger when N consecutive windows reach P (default 0.8) ...
                        (transformer: when one token reaches P; default from the model)
  --trigger-windows N   ... (default 8, which adds 3.5 s)
  --instant-threshold P or when one window reaches P (default 0.9; > 1 disables)
  --release P           probability below which a trigger resets (default 0.3;
                        transformer: the model's release ratio times the threshold)
  --detector KIND       6s (default; also "model"): the 3-seed 6 s window detector
                        transformer: the streaming onset transformer, every 0.1 s
                        stalta: the reference STA/LTA trigger
                        filterpicker: Lomax's FilterPicker5 on the vertical
  --transformer FILE    transformer weights (default: MODELS/transformer.ayzw)
  --dt-reset BELOW,FROM transformer: also trigger while p stays high when dt falls
                        back to BELOW s or less after reaching FROM s since the
                        last trigger: a new onset in the coda of the last one
                        (default: the model's, config.trigger; 2,5 without)
  --no-dt-reset         transformer: rising edges only
  --retrigger S         minimum time between a station's triggers (default 15;
                        transformer: 5, so an aftershock 5-15 s into the coda
                        of the last event can still trigger)
  --no-pd-magnitude     no magnitude from peak P displacement (default: on
                        when MODELS/pd_relation.csv exists)
  --no-merge            transformer: do not merge an event into a declared one
                        whose P times it fits (one earthquake, one alarm)
  --pick-anywhere       transformer: let the picker search the whole 60 s window
                        (default: P within 3 s of the transformer's P, S before
                        the station's next trigger)
  --no-anchor           do not date model triggers by the STA/LTA onset
  --no-onset-fp         do not date the P-time location and Pd of model
                        triggers by a FilterPicker pick (default: the pick
                        within 3 s of the trigger's P, if any)
  --no-p-locate         no location from the detections' P times alone
  --p-sigma S           P-time uncertainty of that location, s (default 1.35)
  --anchor-on R         STA/LTA ratio taken as the onset when anchoring (default 3)
  --require-onset       a detector trigger also needs an STA/LTA onset
  --anchor-picker       place the picker window from the anchored P too
  --sta S, --lta S      STA/LTA averaging lengths in seconds (default 1, 30)
  --stalta-on R         STA/LTA ratio that triggers (default 8)
  --stalta-on-stations FILE
                        per-station STA/LTA trigger ratios (CSV station,ratio);
                        stations not listed keep --stalta-on
  --stalta-off R        STA/LTA ratio below which a trigger resets (default 1.5)
  --stalta-band LO,HI   STA/LTA pass band in Hz (default 2,20)
  --stalta-3c           STA/LTA on the energy of all three components (default: vertical)
  --fp-windows F,L,U    FilterPicker filter window, long-term window and tUpEvent
                        in seconds (default 3,5,0.2: Lomax et al. 2012 at 100 Hz)
  --fp-thresholds T1,T2 FilterPicker trigger and acceptance thresholds (default 10,10)
  --step N              samples between detector windows (default 50 = 0.5 s)
  --min-stations N      detections needed to declare an event (default 2;
                        transformer: 3)
  --slack S             tolerance on the inter-station P travel time, seconds (default 3)
  --locate KIND         geometry: from the transformer's per-station distance and
                        back-azimuth, from the trigger on; the picker is not run
                        picks: grid search on P and S picks (60 s picker window)
                        (default: geometry when the transformer has the geometry
                        head, otherwise picks)
  --geo-sd-scale X      geometry: multiply every station's distance sd by X
                        (default 2.1)
  --geo-sigma-p S       geometry: P-time uncertainty of its misfit, s (default 2.75)
  --geo-max-z Z         geometry: leave out, or else refuse, a station whose
                        distance is more than Z sds from the solution (default off)
  --geo-max-err-km KM   geometry: refuse a location whose 68% radius exceeds KM
                        (default off)
  --assess              diagnostic: judge each alarm earthquake, possible, misfire
                        or unclassified from its S-P picks (the picker runs even
                        with the geometry locator); an ASSESS line as picks
                        arrive and a section in the report
  --assess-csv FILE     with --assess: one row per alarm, readable by
                        tools/plot_candidates.py
  --no-pick             detector only
  --no-magnitude        no magnitude estimate (as --no-pd-magnitude); with
                        --from, reading also starts there where no STA/LTA
                        or FilterPicker needs earlier data
  --catalog CSV         AFAD catalogue export to score events against
  --scores DIR          write every window's probability to DIR/STATION.csv
                        (DIR is created if missing)
  --scores-in DIR       use probabilities from an earlier --scores run instead of the detector
  --record FILE         write all station outputs for tools/network_subsets
  --site NAME,LAT,LON   report the S-wave warning time at this place (repeatable)
  --verbose             print every station detection, pick and magnitude estimate
  --no-color

live input, in place of the files (docs/impl/19-seedlink.md):
  --seedlink HOST:PORT  receive from this SeedLink server (port 18000 if omitted)
  --streams LIST        NET.STA[.LOC],...: each station's HHZ, HHN and HHE at
                        100 Hz; coordinates from the station table
  --archive DIR         append every record received to DIR/STATION.mseed
  --duration S          stop after S seconds (default: at SIGINT or SIGTERM)
  --live-wait S         about how long a station's message waits for the
                        other stations' data before it is processed
                        (default 10)
  --seedlink-timeout S  reconnect after S seconds with nothing received
                        (default 120)
  --seedlink-keepalive S
                        send INFO ID after S seconds with nothing received
                        (default 60; 0: never)
)";

// Seconds since the Unix epoch by the system clock.
double wall_epoch() {
  return std::chrono::duration<double>(
             std::chrono::system_clock::now().time_since_epoch())
      .count();
}

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
  std::string stations_path, pd_gains_path, stalta_on_stations;
  std::string models = "models", scores_dir, scores_in_dir, catalog_path,
              record_path, transformer_path, assess_csv;
  bool threshold_set = false, release_set = false, min_stations_set = false,
       dt_reset_set = false, retrigger_set = false, no_merge = false,
       no_pd = false;
  std::string locate_kind; // "", "geometry" or "picks"
  double speed = 1.0;
  ProcessorConfig pcfg;
  NetworkConfig ncfg;
  std::optional<PdModel> pd_model; // outlives the processors and the network
  std::vector<std::string> files;
  // Live input (docs/impl/19-seedlink.md).
  std::string seedlink_server, streams;
  LiveConfig lcfg;
  double live_wait = 10.0;
  bool speed_set = false, live_option = false;

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc)
        throw std::runtime_error(a + " needs a value");
      return argv[++i];
    };
    if (a == "--stations")
      stations_path = next();
    else if (a == "--pd-gains")
      pd_gains_path = next();
    else if (a == "--models")
      models = next();
    else if (a == "--speed")
      speed = std::stod(next()), speed_set = true;
    else if (a == "--from")
      pcfg.from = parse_time(next());
    else if (a == "--threshold") {
      pcfg.threshold = std::stof(next());
      threshold_set = true;
    } else if (a == "--transformer")
      transformer_path = next();
    else if (a == "--dt-reset") {
      const auto v = next();
      const auto comma = v.find(',');
      if (comma == std::string::npos)
        throw std::runtime_error("--dt-reset: BELOW,FROM");
      pcfg.dt_reset_below = std::stof(v.substr(0, comma));
      pcfg.dt_reset_from = std::stof(v.substr(comma + 1));
      if (!(pcfg.dt_reset_below < pcfg.dt_reset_from))
        throw std::runtime_error("--dt-reset: BELOW must be less than FROM");
      pcfg.dt_reset = true, dt_reset_set = true;
    } else if (a == "--no-dt-reset")
      pcfg.dt_reset = false, dt_reset_set = true;
    else if (a == "--retrigger")
      pcfg.retrigger_seconds = std::stod(next()), retrigger_set = true;
    else if (a == "--no-merge")
      no_merge = true;
    else if (a == "--no-pd-magnitude")
      no_pd = true;
    else if (a == "--pick-anywhere")
      pcfg.pick_on_onset = false;
    else if (a == "--trigger-windows")
      pcfg.trigger_windows = std::stoul(next());
    else if (a == "--instant-threshold")
      pcfg.instant_threshold = std::stof(next());
    else if (a == "--release") {
      pcfg.release = std::stof(next());
      release_set = true;
    }
    else if (a == "--detector") {
      const auto kind = next();
      if (kind == "6s" || kind == "model")
        pcfg.detector = DetectorKind::Model;
      else if (kind == "transformer")
        pcfg.detector = DetectorKind::Transformer;
      else if (kind == "stalta")
        pcfg.detector = DetectorKind::StaLta;
      else if (kind == "filterpicker")
        pcfg.detector = DetectorKind::FilterPicker;
      else
        throw std::runtime_error(
            "--detector: 6s, transformer, stalta or filterpicker");
    } else if (a == "--sta")
      pcfg.stalta.sta_seconds = std::stod(next());
    else if (a == "--lta")
      pcfg.stalta.lta_seconds = std::stod(next());
    else if (a == "--stalta-on-stations")
      stalta_on_stations = next();
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
    else if (a == "--fp-windows") {
      const auto v = next();
      const auto c1 = v.find(','), c2 = v.find(',', c1 + 1);
      if (c1 == std::string::npos || c2 == std::string::npos)
        throw std::runtime_error("--fp-windows: FILTER,LONG,TUP");
      pcfg.fp.filter_window = std::stod(v.substr(0, c1));
      pcfg.fp.long_term_window = std::stod(v.substr(c1 + 1, c2 - c1 - 1));
      pcfg.fp.t_up = std::stod(v.substr(c2 + 1));
    } else if (a == "--fp-thresholds") {
      const auto v = next();
      const auto comma = v.find(',');
      if (comma == std::string::npos)
        throw std::runtime_error("--fp-thresholds: T1,T2");
      pcfg.fp.threshold1 = std::stod(v.substr(0, comma));
      pcfg.fp.threshold2 = std::stod(v.substr(comma + 1));
    }
    else if (a == "--no-onset-fp")
      pcfg.onset_fp = false;
    else if (a == "--no-p-locate")
      ncfg.p_locate = false;
    else if (a == "--p-sigma")
      ncfg.p_sigma = std::stod(next());
    else if (a == "--no-anchor")
      pcfg.anchor = false;
    else if (a == "--anchor-on")
      pcfg.anchor_on = std::stod(next());
    else if (a == "--anchor-picker")
      pcfg.anchor_picker = true;
    else if (a == "--require-onset")
      pcfg.require_onset = true;
    else if (a == "--step")
      pcfg.step = std::stoul(next());
    else if (a == "--slack")
      ncfg.slack_seconds = std::stod(next());
    else if (a == "--min-stations")
      ncfg.min_stations = std::stoul(next()), min_stations_set = true;
    else if (a == "--locate") {
      locate_kind = next();
      if (locate_kind != "geometry" && locate_kind != "picks")
        throw std::runtime_error("--locate: geometry or picks");
    } else if (a == "--geo-sd-scale")
      ncfg.geo_sd_scale = std::stod(next());
    else if (a == "--geo-sigma-p")
      ncfg.sigma_p = std::stod(next());
    else if (a == "--geo-max-z")
      ncfg.geo_max_dist_z = std::stod(next());
    else if (a == "--geo-max-err-km")
      ncfg.geo_max_err_km = std::stod(next());
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
    else if (a == "--assess")
      ncfg.assess = true;
    else if (a == "--assess-csv")
      assess_csv = next();
    else if (a == "--seedlink")
      seedlink_server = next();
    else if (a == "--streams")
      streams = next(), live_option = true;
    else if (a == "--archive")
      lcfg.archive_dir = next(), live_option = true;
    else if (a == "--duration")
      lcfg.duration_s = std::stod(next()), live_option = true;
    else if (a == "--live-wait")
      live_wait = std::stod(next()), live_option = true;
    else if (a == "--seedlink-timeout")
      lcfg.client.timeout_s = std::stod(next()), live_option = true;
    else if (a == "--seedlink-keepalive")
      lcfg.client.keepalive_s = std::stod(next()), live_option = true;
    else if (a == "-h" || a == "--help") {
      std::print("{}", kUsage);
      return 0;
    } else if (a.starts_with("--")) {
      throw std::runtime_error("unknown option " + a);
    } else {
      files.push_back(a);
    }
  }
  const bool live = !seedlink_server.empty();
  if (live) {
    if (!files.empty())
      throw std::runtime_error(
          "--seedlink takes the place of the miniSEED files; give one or the "
          "other");
    if (streams.empty())
      throw std::runtime_error("--seedlink needs --streams NET.STA[,...]");
    if (speed_set)
      throw std::runtime_error(
          "--speed applies to replays; live input arrives at its own pace");
    if (ncfg.assess)
      throw std::runtime_error(
          "--assess reads the replay files after the run; replay an "
          "--archive instead");
    if (!seedlink::Client::available())
      throw std::runtime_error(
          "live SeedLink input is not available in this build (it needs "
          "POSIX sockets, which the Windows build does not have)");
    // NaN fails every comparison, so each check is written to fail on it.
    if (!(lcfg.duration_s >= 0) || std::isinf(lcfg.duration_s))
      throw std::runtime_error("--duration must be 0 (until stopped) or more");
    if (!(live_wait >= 0) || std::isinf(live_wait))
      throw std::runtime_error("--live-wait must be 0 or more");
    if (!(lcfg.client.timeout_s > 0) || std::isinf(lcfg.client.timeout_s))
      throw std::runtime_error("--seedlink-timeout must be more than 0");
    if (!(lcfg.client.keepalive_s >= 0) || std::isinf(lcfg.client.keepalive_s))
      throw std::runtime_error("--seedlink-keepalive must be 0 (none) or more");
    if (lcfg.client.keepalive_s > 0 &&
        lcfg.client.keepalive_s >= lcfg.client.timeout_s)
      throw std::runtime_error(
          "--seedlink-keepalive must be shorter than --seedlink-timeout");
    lcfg.client.server = seedlink::parse_endpoint(seedlink_server);
    lcfg.client.streams = seedlink::parse_streams(streams);
  } else if (live_option) {
    throw std::runtime_error("--streams, --archive, --duration, --live-wait "
                             "and --seedlink-* need --seedlink");
  }
  if (files.empty() && !live) {
    std::print(stderr, "{}", kUsage);
    return 2;
  }
  if (pcfg.detector != DetectorKind::Model && !scores_in_dir.empty())
    throw std::runtime_error("--scores-in applies to the 6 s detector only");
  // FilterPicker dates P itself, from its uncertainty point.
  if (pcfg.detector == DetectorKind::FilterPicker) {
    pcfg.anchor = false;
    pcfg.require_onset = false;
  }
  if (!stdout_is_tty())
    Log::get().color = false;
  if (live) {
    const auto on_signal = [](int) {
      if (g_input_stop.exchange(true))
        g_stop = true;
    };
    std::signal(SIGINT, on_signal);
    std::signal(SIGTERM, on_signal);
  } else {
    std::signal(SIGINT, [](int) { g_stop = true; });
  }

  // --- load
  // ---------------------------------------------------------------------
  std::vector<Weights> detector;
  for (int seed : {42, 43, 44})
    detector.push_back(
        Weights::load(std::format("{}/detector_s{}.ayzw", models, seed)));
  const Weights picker = Weights::load(models + "/spicker.ayzw");
  // The transformer dates P itself (dt), so the STA/LTA anchor is not used,
  // and its threshold is the validation operating point it was exported with.
  std::unique_ptr<Weights> transformer;
  if (pcfg.detector == DetectorKind::Transformer) {
    transformer = std::make_unique<Weights>(Weights::load(
        transformer_path.empty() ? models + "/transformer.ayzw" : transformer_path));
    const auto f = transformer->at("config.floats").f64();
    if (!threshold_set)
      pcfg.threshold = static_cast<float>(f[3]);
    if (!release_set)
      pcfg.release = static_cast<float>(f[4] * pcfg.threshold);
    // The dt restart the model was exported with: how far dt must fall
    // depends on how deep a given model's restarts go.
    if (!dt_reset_set && transformer->has("config.trigger")) {
      const auto g = transformer->at("config.trigger", {3}).f64();
      pcfg.dt_reset_below = static_cast<float>(g[0]);
      pcfg.dt_reset_from = static_cast<float>(g[1]);
      pcfg.dt_reset_tokens = static_cast<std::size_t>(g[2]);
    }
    pcfg.anchor = false;
    pcfg.require_onset = false;
    // A dt restart dates its own P, so a short gap between triggers is safe
    // and keeps aftershocks close behind the last event: on onset's wide
    // validation set, second onsets within 1 s go from 38% to 59% at 15 -> 5
    // s, for a few more triggers in codas.
    if (!retrigger_set)
      pcfg.retrigger_seconds = 5.0;
    // With triggers 5 s apart a station can trigger twice in one earthquake,
    // which split events into duplicate alarms (14-transformer.md).
    ncfg.merge_duplicates = !no_merge;
    // Two stations are too easy to satisfy for a detector that also fires in
    // codas (dt restarts): on the Sindirgi and Marmara sequences two far
    // stations ringing together made most of the false alarms, and a third
    // station removes them at the cost of a few small events
    // (docs/impl/14-transformer.md).
    if (!min_stations_set)
      ncfg.min_stations = 3;
  }
  // Magnitude from peak P displacement at every detector's onset
  // (docs/impl/17-pd-magnitude.md).
  if (pcfg.magnitude && !no_pd && (pd_model = PdModel::load(models))) {
    if (!pd_gains_path.empty())
      pd_model->use_gains(pd_gains_path);
    pcfg.pd = &*pd_model;
    ncfg.pd = &*pd_model;
  }
  // Location from the transformer's geometry head replaces the S-P picks: it
  // is available from the trigger on instead of after the 60 s picker window.
  const bool has_geometry = transformer && transformer->has("geo_head.weight");
  if (locate_kind == "geometry" && !has_geometry)
    throw std::runtime_error(
        "--locate geometry needs --detector transformer with a model exported "
        "with the geometry head");
  if (locate_kind == "geometry" || (locate_kind.empty() && has_geometry)) {
    ncfg.locator = Locator::Geometry;
    pcfg.geometry = true;
    pcfg.pick = false;
  }
  // The assessment needs S-P: the picker runs for it whatever locates.
  if (!assess_csv.empty() && !ncfg.assess)
    throw std::runtime_error("--assess-csv needs --assess");
  if (ncfg.assess)
    pcfg.pick = true;
  const auto bp = dsp::Bandpass::load(Weights::load(models + "/bandpass.ayzw"));
  const std::string station_table =
      stations_path.empty() ? models + "/stations.csv" : stations_path;
  auto coords = load_stations(station_table);

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
  // Live input: one station per stream, and every one in the station table,
  // since a station that cannot be located from would be a silent hole in a
  // live network.
  for (const auto &spec : lcfg.client.streams) {
    auto st = std::make_unique<Station>();
    st->code = spec.sta;
    const auto it = coords.find(st->code);
    if (it == coords.end())
      throw std::runtime_error(std::format(
          "{}: not in the station table {}; live input needs the "
          "coordinates of every station (--stations FILE)",
          spec.id(), station_table));
    st->lat = it->second.lat;
    st->lon = it->second.lon;
    used[st->code] = it->second;
    stations.push_back(std::move(st));
  }
  std::unique_ptr<LiveInput> live_in;
  if (live) {
    std::vector<Station *> fed;
    for (auto &st : stations)
      fed.push_back(st.get());
    live_in = std::make_unique<LiveInput>(lcfg, std::move(fed));
  }

  // 70 s covers a picker window before `from`. Live input knows t_first only
  // at the end, where `start` is set again.
  double start = pcfg.from > 0 ? pcfg.from - 70.0 : t_first;
  const std::string detector_name =
      pcfg.detector == DetectorKind::Transformer
          ? std::format("onset transformer (p >= {:g}, release {:g}, {}P from dt)",
                        pcfg.threshold, pcfg.release,
                        pcfg.dt_reset
                            ? std::format("dt restart {:g}/{:g} s, ",
                                          pcfg.dt_reset_below, pcfg.dt_reset_from)
                            : std::string())
      : pcfg.detector == DetectorKind::Model
          ? std::string(pcfg.anchor ? "3-seed 6 s detector + STA/LTA anchor"
                                    : "3-seed 6 s detector")
      : pcfg.detector == DetectorKind::FilterPicker
          ? std::format(
                "FilterPicker5 {:g}/{:g}/{:g} s, thresholds {:g}/{:g}, Z",
                pcfg.fp.filter_window, pcfg.fp.long_term_window, pcfg.fp.t_up,
                pcfg.fp.threshold1, pcfg.fp.threshold2)
          : std::format(
                "STA/LTA {:g}/{:g} s on {:g} off {:g}, {:g}-{:g} Hz, {}",
                pcfg.stalta.sta_seconds, pcfg.stalta.lta_seconds,
                pcfg.stalta_on, pcfg.stalta_off, pcfg.stalta.f_lo,
                pcfg.stalta.f_hi, pcfg.stalta.three_component ? "Z+N+E" : "Z");
  const std::string input =
      live ? std::format("live from {}:{}", lcfg.client.server.host,
                         lcfg.client.server.port)
           : std::format("replay {} to {} UTC at {}", ymd_hms(t_first),
                         hms(t_last),
                         speed > 0 ? std::format("{:g}x", speed)
                                   : std::string("full speed"));
  Log::get().line("ayzek", "1", "{} stations, {} backend, {}{}{}, {}",
                  stations.size(), simd::kBackend, detector_name,
                  pcfg.geometry ? " + geometry locator"
                  : pcfg.pick   ? " + picker"
                                : "",
                  pcfg.pd ? " + Pd magnitude" : "", input);
  if (!catalog_path.empty() && live) {
    // The span is known at the end, when the report keeps the events a
    // replay of the same data would (set_catalog below).
    ncfg.catalog =
        load_afad_catalog(catalog_path, -std::numeric_limits<double>::infinity(),
                          std::numeric_limits<double>::infinity());
    Log::get().line("ayzek", "1",
                    "catalogue: {} events; the report keeps those in the "
                    "span received",
                    ncfg.catalog.size());
  } else if (!catalog_path.empty()) {
    // Catalogue events whose waves can arrive within the replayed data.
    ncfg.catalog = load_afad_catalog(
        catalog_path, std::max(t_first, pcfg.from) - 30, t_last - 30);
    Log::get().line("ayzek", "1", "catalogue: {} events in the replay span",
                    ncfg.catalog.size());
    // For the assessment's aftershock-zone context: the 30 days before too.
    if (ncfg.assess)
      ncfg.context_catalog = load_afad_catalog(
          catalog_path, std::max(t_first, pcfg.from) - 30 * 86400.0, t_last);
  }

  // --- run
  // ----------------------------------------------------------------------
  if (!scores_dir.empty())
    std::filesystem::create_directories(scores_dir);
  Bus bus;
  // Live input runs on the wall clock, and releases records on arrival.
  StreamClock clock(live ? wall_epoch() : start, live ? 1.0 : speed);
  std::vector<std::unique_ptr<Processor>> procs;
  // --stalta-on-stations: the per-station STA/LTA baseline, each station's own
  // trigger ratio.
  std::map<std::string, double> station_on;
  if (!stalta_on_stations.empty()) {
    std::ifstream f(stalta_on_stations);
    if (!f)
      throw std::runtime_error("cannot open " + stalta_on_stations);
    std::string line;
    std::getline(f, line); // header
    while (std::getline(f, line)) {
      const auto comma = line.find(',');
      if (comma != std::string::npos)
        station_on[line.substr(0, comma)] = std::stod(line.substr(comma + 1));
    }
  }
  for (auto &st : stations) {
    const std::string path =
        scores_dir.empty() ? "" : scores_dir + "/" + st->code + ".csv";
    const std::string in_path =
        scores_in_dir.empty() ? "" : scores_in_dir + "/" + st->code + ".csv";
    auto cfg = pcfg;
    if (auto it = station_on.find(st->code); it != station_on.end())
      cfg.stalta_on = it->second;
    procs.push_back(std::make_unique<Processor>(
        *st, detector, picker, bp, cfg, bus, path, in_path, transformer.get()));
  }
  // Made before the threads start, so that an output that cannot be written
  // ends the run before it begins.
  Network net(used, ncfg);
  std::unique_ptr<Recorder> recorder;
  if (!record_path.empty())
    recorder = live ? std::make_unique<Recorder>(record_path, used,
                                                 ncfg.min_stations, pcfg.from)
                    : std::make_unique<Recorder>(record_path, used, t_first,
                                                 t_last, ncfg.min_stations,
                                                 pcfg.from);
  std::vector<std::jthread> threads;
  // An exception from here on unwinds through `threads`, whose joins would
  // wait for threads that end only on a signal, --duration or the end of the
  // data (live input: never). Destroyed before them, this stops them first.
  struct StopOnUnwind {
    int at = std::uncaught_exceptions();
    ~StopOnUnwind() {
      if (std::uncaught_exceptions() > at) {
        g_input_stop = true;
        g_stop = true;
      }
    }
  } stop_on_unwind;
  if (live)
    threads.emplace_back([&] { live_in->run(g_input_stop, g_stop); });
  for (std::size_t i = 0; i < stations.size(); ++i) {
    if (!live)
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
  std::map<std::string, double> watermark;
  for (auto &st : stations)
    watermark[st->code] = 0;
  struct Held {
    double t;
    std::string station;
    std::uint64_t seq;
    Message m;
    std::chrono::steady_clock::time_point received{}; // live input only
  };
  auto later = [](const Held &a, const Held &b) {
    return std::tie(a.t, a.station, a.seq) > std::tie(b.t, b.station, b.seq);
  };
  std::priority_queue<Held, std::vector<Held>, decltype(later)> held(later);
  std::map<std::string, std::uint64_t> sent;

  auto process = [&](const Message &m) {
    if (auto *d = std::get_if<Detection>(&m))
      net.on(*d);
    else if (auto *p = std::get_if<Pick>(&m))
      net.on(*p);
    else if (auto *o = std::get_if<StationGeometry>(&m))
      net.on(*o);
    else if (auto *q = std::get_if<PdEstimate>(&m))
      net.on(*q);
    else if (auto *t = std::get_if<PdTrace>(&m))
      net.on(*t);
    if (recorder)
      recorder->write(m);
  };
  auto release = [&] {
    const double safe = watermark.empty()
                            ? 1e18
                            : std::ranges::min(watermark | std::views::values);
    while (!held.empty() && held.top().t < safe) {
      const Message m = held.top().m;
      held.pop();
      process(m);
    }
  };

  // Live input: a message waits for every station's promise as in a replay,
  // but for --live-wait seconds at most. The stations that held it back that
  // long are marked lagging and not waited for until they have caught up with
  // what has been processed (a station that stops sending would otherwise
  // hold every alarm back). What a lagging station sends meanwhile for
  // earlier times is processed out of order, and counted: with none out of
  // order, the network stage saw the messages in the replay's order, so a
  // replay of the same records gives the same alarms.
  const auto wait_limit = std::chrono::duration<double>(live_wait);
  std::map<std::string, bool> lagging;
  std::tuple<double, std::string, std::uint64_t> last_key{
      -std::numeric_limits<double>::infinity(), "", 0};
  std::uint64_t out_of_order = 0, waited_out = 0;
  // Alarms as they are issued: id, stream time, wall-clock time.
  std::set<int> issued_ids;
  std::vector<std::tuple<int, double, double>> issued;
  auto release_live = [&] {
    const auto now = std::chrono::steady_clock::now();
    while (!held.empty()) {
      double safe = 1e18;
      for (const auto &[code, w] : watermark)
        if (!lagging[code])
          safe = std::min(safe, w);
      const Held &top = held.top();
      if (top.t >= safe) {
        if (now - top.received < wait_limit)
          break;
        for (const auto &[code, w] : watermark)
          if (!lagging[code] && w <= top.t) {
            lagging[code] = true;
            ++waited_out;
            Log::get().line("warn", "33",
                            "{}: held a message back {:g} s ({}); not waited "
                            "for until it catches up",
                            code, live_wait,
                            w > 0 ? "processed to " + hms(w)
                                  : std::string("nothing processed yet"));
          }
        continue;
      }
      auto key = std::tie(top.t, top.station, top.seq);
      if (key < last_key)
        ++out_of_order;
      else
        last_key = key;
      const Message m = top.m;
      held.pop();
      process(m);
    }
    for (auto &[code, lag] : lagging)
      if (const auto w = watermark.find(code);
          lag && w != watermark.end() && w->second > std::get<0>(last_key)) {
        lag = false;
        Log::get().line("live", "36", "{}: caught up, waited for again", code);
      }
    for (const auto &e : net.events())
      if (e.declared && issued_ids.insert(e.id).second) {
        const double wall = wall_epoch();
        issued.emplace_back(e.id, e.declared_at, wall);
        Log::get().line("ISSUED", "1;31",
                        "#{:<3} {}  at {} UTC by the wall clock, {:.1f} s after "
                        "the data",
                        e.id, hms(e.declared_at), hms(wall),
                        wall - e.declared_at);
      }
  };

  while (!watermark.empty()) {
    auto m = bus.receive(std::chrono::milliseconds(100));
    if (!m) {
      if (live) // the wait limit runs out without messages too
        release_live();
      continue;
    }
    if (auto *pr = std::get_if<Progress>(&*m)) {
      if (auto it = watermark.find(pr->station); it != watermark.end())
        it->second = std::max(it->second, pr->until);
    } else if (auto *sd = std::get_if<StationDone>(&*m)) {
      watermark.erase(sd->station);
    } else {
      const std::string station = station_of(*m);
      held.push(Held{declared_at(*m), station, sent[station]++, std::move(*m),
                     std::chrono::steady_clock::now()});
    }
    if (live)
      release_live();
    else
      release();
  }
  if (live)
    release_live();
  else
    release();
  threads.clear();

  if (live) {
    if (!live_in->error().empty())
      throw std::runtime_error(live_in->error());
    // The span received, as a replay of the same records would read it; the
    // wall-clock span of the run if nothing arrived.
    t_first = live_in->first_time();
    t_last = live_in->last_time();
    if (std::isnan(t_first)) {
      t_last = wall_epoch();
      t_first = t_last - clock.wall_seconds();
    }
    start = pcfg.from > 0 ? pcfg.from - 70.0 : t_first;
    if (!catalog_path.empty())
      net.set_catalog(load_afad_catalog(
          catalog_path, std::max(t_first, pcfg.from) - 30, t_last - 30));
    if (recorder)
      recorder->write_span(t_first, t_last);
  }

  // --assess: picks at the stations that did not trigger, from the replay
  // files, in the window each alarm's source predicts. Made here, after the
  // replay and on this thread, so they do not depend on thread timing; a live
  // client would have to keep that data instead (docs/impl/16-alarm-assessment.md).
  if (ncfg.assess) {
    Picker extra_picker(picker);
    dsp::Conditioner cond(bp, Picker::kWindow);
    std::map<std::string, const ReplaySource *> source_of;
    for (const auto &src : sources)
      source_of[src->station()] = src.get();
    std::vector<double> raw;
    std::vector<float> planar(Picker::kWindow * 3);
    std::size_t made = 0, missing = 0;
    constexpr double kLead = 5.5; // window start before the predicted P, as for triggers
    for (const auto &r : net.pick_requests(std::max(start, t_first), t_last)) {
      auto it = source_of.find(r.station);
      const double t0 = r.p_pred - kLead;
      if (it == source_of.end() || !read_window(*it->second, t0, Picker::kWindow, raw)) {
        ++missing;
        continue;
      }
      const auto c0 = std::chrono::steady_clock::now();
      cond.condition_planar(raw, 3, planar);
      const auto pk = extra_picker.pick(
          planar, {kLead - 4.0, kLead + 4.0,
                   std::min(Picker::kWindow / kFs, r.s_cap - t0), 0.5});
      const double ms = std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - c0)
                            .count();
      Pick p{r.station, r.p_pred - 3.5,  t0 + pk.p_seconds, t0 + pk.s_seconds,
             pk.p_prob, pk.s_prob,       t0 + Picker::kWindow / kFs,
             ms,        r.event_id,      r.event_alarm};
      net.on(p);
      if (recorder)
        recorder->write(p);
      ++made;
    }
    Log::get().line("assess", "1",
                    "{} picks at stations that did not trigger{}", made,
                    missing ? std::format(", {} windows skipped for missing data",
                                          missing)
                            : std::string());
  }
  const double wall = clock.wall_seconds();

  // --- report
  // -----------------------------------------------------------------
  const double stream_seconds = t_last - std::max(start, t_first);
  net.report(std::max(start, t_first), t_last);
  if (!assess_csv.empty()) {
    const auto slash = record_path.find_last_of("/\\");
    if (!net.write_assessments(
            assess_csv, std::max(start, t_first), t_last,
            slash == std::string::npos ? record_path : record_path.substr(slash + 1)))
      throw std::runtime_error("cannot write " + assess_csv);
  }
  if (recorder) {
    recorder.reset(); // flush and close
    Log::get().plain(true, "Recording written to {}",
                     std::filesystem::absolute(record_path).string());
  }

  auto &log = Log::get();
  log.plain(false, "");
  log.plain(true, "Station processing ({} backend)", simd::kBackend);
  log.plain(false,
            "  {:<5} {:>8} {:>9} {:>8} {:>6} {:>9}   {:>16} {:>9}", "",
            "windows", "with gap", "triggers", "picks", "Pd values",
            "ms/window (p99)", "ms/pick");
  std::uint64_t windows = 0, triggers = 0, anchored = 0, onset_picks = 0;
  for (std::size_t i = 0; i < stations.size(); ++i) {
    const auto &s = procs[i]->stats();
    const auto w = percentiles(s.window_ms);
    const auto pk = percentiles(s.pick_ms);
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
    onset_picks += s.onset_picks;
    if (stale)
      log.line("warn", "33", "{}: {} late records discarded", stations[i]->code,
               stale);
    if (const auto bad = stations[i]->undecodable.load())
      log.line("warn", "33", "{}: {} records did not decode", stations[i]->code,
               bad);
    log.plain(false,
              "  {:<5} {:>8} {:>9} {:>8} {:>6} {:>9}   {:>8.3f} ({:>5.3f}) "
              "{:>9.0f}",
              stations[i]->code, s.windows, s.gap_windows, s.detections,
              s.picks, s.pd_values, w.mean, w.p99, pk.mean);
  }
  if (pcfg.detector == DetectorKind::Transformer) {
    std::uint64_t refreshes = 0, geometry = 0, dt_resets = 0;
    for (const auto &pr : procs) {
      refreshes += pr->stats().context_refreshes;
      geometry += pr->stats().geometry;
      dt_resets += pr->stats().dt_resets;
    }
    log.plain(false,
              "  transformer: \"windows\" are 0.1 s tokens; {} station context "
              "refreshes; {} of {} triggers from a dt restart{}",
              refreshes, dt_resets, triggers,
              (pcfg.geometry
                   ? std::format("; {} geometry estimates sent", geometry)
                   : std::string()));
  }
  if (pcfg.pd) {
    std::uint64_t pd_values = 0;
    for (const auto &pr : procs)
      pd_values += pr->stats().pd_values;
    log.plain(false, "  Pd: {} values sent", pd_values);
  }
  if (pcfg.detector == DetectorKind::Model && pcfg.anchor)
    log.plain(false,
              "  P time from an STA/LTA onset for {} of {} triggers ({} from "
              "the window start)",
              anchored, triggers, triggers - anchored);
  if (pcfg.detector == DetectorKind::Model && pcfg.onset_fp)
    log.plain(false,
              "  onset from a FilterPicker pick for {} of {} triggers (P-time "
              "location and Pd)",
              onset_picks, triggers);
  log.plain(false, "  {} windows in {:.1f} s: {:.0f}x real time per station",
            windows, wall, stream_seconds / wall);
  if (live) {
    live_in->report();
    log.plain(false,
              "  network stage: {} messages processed out of order; {} times "
              "a station was not waited for (--live-wait {:g} s)",
              out_of_order, waited_out, live_wait);
    for (const auto &[id, at, wall_at] : issued)
      log.plain(false,
                "  alarm #{:<3} {} UTC issued at {} UTC by the wall clock, "
                "{:.1f} s after the data",
                id, hms(at), hms(wall_at), wall_at - at);
  }
  return 0;
} catch (const std::exception &e) {
  std::println(stderr, "ayzek: {}", e.what());
  return 1;
}
