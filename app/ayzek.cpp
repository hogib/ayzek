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
#include "pipeline/station.hpp"
#include "simd.hpp"

#include <algorithm>
#include <chrono>
#include <csignal>
#include <cstdio>
#include <format>
#include <fstream>
#include <memory>
#include <numeric>
#include <queue>
#include <ranges>
#include <sstream>
#include <thread>
#include <unistd.h>

using namespace ayzek;
using namespace ayzek::pipeline;

namespace {

std::atomic<bool> g_stop{false};

const char* kUsage = R"(usage: ayzek [options] STATION.mseed...

  --models DIR          weights directory (default: models)
  --speed X             replay speed; 1 = real time, 0 = as fast as possible (default 1)
  --from TIME           start detecting at this UTC time, e.g. 2025-11-10T18:19:00
  --threshold P         detector trigger probability (default 0.9; outputs saturate near 0.91)
  --step N              samples between detector windows (default 50 = 0.5 s)
  --min-stations N      detections needed to declare an event (default 2)
  --no-pick             detector only
  --no-magnitude        skip the magnitude regressor
  --catalog CSV         AFAD catalogue export to score events against
  --scores DIR          write every window's probability to DIR/STATION.csv
  --no-color
)";

double parse_time(const std::string& s) {
    const auto t = parse_utc(s, "ymdhms");
    if (!t) throw std::runtime_error("bad time: " + s + " (want YYYY-MM-DDTHH:MM:SS)");
    return *t;
}

std::map<std::string, StationInfo> load_stations(const std::string& path) {
    std::map<std::string, StationInfo> out;
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::string line;
    std::getline(f, line);
    while (std::getline(f, line)) {
        std::stringstream ss(line);
        std::string code, lat, lon;
        std::getline(ss, code, ',');
        std::getline(ss, lat, ',');
        std::getline(ss, lon, ',');
        // Some station codes occur more than once in the AFAD table; the first is used.
        if (!code.empty() && !out.contains(code)) out[code] = {std::stod(lat), std::stod(lon)};
    }
    return out;
}

struct Percentiles {
    double mean, p50, p99, max;
};
Percentiles percentiles(std::vector<float> v) {
    if (v.empty()) return {0, 0, 0, 0};
    std::ranges::sort(v);
    const auto at = [&](double q) { return static_cast<double>(v[static_cast<std::size_t>(q * static_cast<double>(v.size() - 1))]); };
    return {std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size()), at(0.5), at(0.99), static_cast<double>(v.back())};
}

}  // namespace

int main(int argc, char** argv) try {
    std::string models = "models", scores_dir, catalog_path;
    double speed = 1.0;
    ProcessorConfig pcfg;
    NetworkConfig ncfg;
    std::vector<std::string> files;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(a + " needs a value");
            return argv[++i];
        };
        if (a == "--models") models = next();
        else if (a == "--speed") speed = std::stod(next());
        else if (a == "--from") pcfg.from = parse_time(next());
        else if (a == "--threshold") pcfg.threshold = std::stof(next());
        else if (a == "--step") pcfg.step = std::stoul(next());
        else if (a == "--min-stations") ncfg.min_stations = std::stoul(next());
        else if (a == "--no-pick") pcfg.pick = false;
        else if (a == "--no-magnitude") pcfg.magnitude = false;
        else if (a == "--scores") scores_dir = next();
        else if (a == "--no-color") Log::get().color = false;
        else if (a == "--catalog") catalog_path = next();
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
    if (!isatty(STDOUT_FILENO)) Log::get().color = false;
    std::signal(SIGINT, [](int) { g_stop = true; });

    // --- load ---------------------------------------------------------------------
    std::vector<Weights> detector;
    for (int seed : {42, 43, 44}) detector.push_back(Weights::load(std::format("{}/detector_s{}.ayzw", models, seed)));
    const Weights picker = Weights::load(models + "/spicker.ayzw");
    std::vector<Weights> magnitude;
    if (pcfg.magnitude) {
        for (int p = 0; p < 3; ++p) {
            const auto path = std::format("{}/magnitude_p{}.ayzw", models, p);
            if (std::ifstream(path).good()) magnitude.push_back(Weights::load(path));
        }
        if (magnitude.empty()) {
            Log::get().line("warn", "33", "no magnitude models in {} (tools/export_magnitude.py); magnitude off", models);
            pcfg.magnitude = false;
        }
    }
    const auto bp = dsp::Bandpass::load(Weights::load(models + "/bandpass.ayzw"));
    auto coords = load_stations(models + "/stations.csv");

    std::vector<std::unique_ptr<ReplaySource>> sources;
    std::vector<std::unique_ptr<Station>> stations;
    std::map<std::string, StationInfo> used;
    double t_first = 1e18, t_last = 0;
    for (const auto& f : files) {
        auto src = std::make_unique<ReplaySource>(f);
        auto st = std::make_unique<Station>();
        st->code = src->station();
        if (auto it = coords.find(st->code); it != coords.end()) {
            st->lat = it->second.lat;
            st->lon = it->second.lon;
            used[st->code] = it->second;
        } else {
            Log::get().line("warn", "33", "{}: no coordinates, excluded from location", st->code);
        }
        t_first = std::min(t_first, src->first_time());
        t_last = std::max(t_last, src->last_time());
        sources.push_back(std::move(src));
        stations.push_back(std::move(st));
    }

    const double start = pcfg.from > 0 ? pcfg.from - 70.0 : t_first;   // 70 s covers a picker window before `from`
    Log::get().line("ayzek", "1", "{} stations, {} backend, 3-seed detector{}{}, replay {} to {} UTC at {}",
                    stations.size(), simd::kBackend, pcfg.pick ? " + picker" : "",
                    pcfg.magnitude ? std::format(" + {}-model magnitude", magnitude.size()) : "", ymd_hms(t_first), hms(t_last),
                    speed > 0 ? std::format("{:g}x", speed) : std::string("full speed"));
    if (!catalog_path.empty()) {
        // Catalogue events whose waves can arrive within the replayed data.
        ncfg.catalog = load_afad_catalog(catalog_path, std::max(t_first, pcfg.from) - 30, t_last - 30);
        Log::get().line("ayzek", "1", "catalogue: {} events in the replay span", ncfg.catalog.size());
    }

    // --- run ----------------------------------------------------------------------
    Bus bus;
    StreamClock clock(start, speed);
    std::vector<std::unique_ptr<Processor>> procs;
    for (auto& st : stations) {
        const std::string path = scores_dir.empty() ? "" : scores_dir + "/" + st->code + ".csv";
        procs.push_back(std::make_unique<Processor>(*st, detector, picker, magnitude, bp, pcfg, bus, path));
    }
    std::vector<std::jthread> threads;
    for (std::size_t i = 0; i < stations.size(); ++i) {
        threads.emplace_back([&, i] { run_ingest(*sources[i], *stations[i], clock, g_stop); });
        threads.emplace_back([&, i] { procs[i]->run(g_stop); });
    }

    // Messages are held in a priority queue and processed in declared_at order,
    // up to the minimum Progress time over the stations still running.
    Network net(used, ncfg);
    std::map<std::string, double> watermark;
    for (auto& st : stations) watermark[st->code] = 0;
    auto later = [](const Message& a, const Message& b) { return declared_at(a) > declared_at(b); };
    std::priority_queue<Message, std::vector<Message>, decltype(later)> held(later);
    auto release = [&] {
        const double safe = watermark.empty() ? 1e18 : std::ranges::min(watermark | std::views::values);
        while (!held.empty() && declared_at(held.top()) <= safe) {
            const Message m = held.top();
            held.pop();
            if (auto* d = std::get_if<Detection>(&m)) net.on(*d);
            else if (auto* p = std::get_if<Pick>(&m)) net.on(*p);
            else if (auto* g = std::get_if<MagnitudeEstimate>(&m)) net.on(*g);
        }
    };
    while (!watermark.empty()) {
        auto m = bus.receive(std::chrono::milliseconds(100));
        if (!m) continue;
        if (auto* pr = std::get_if<Progress>(&*m)) {
            if (auto it = watermark.find(pr->station); it != watermark.end()) it->second = std::max(it->second, pr->until);
        } else if (auto* sd = std::get_if<StationDone>(&*m)) {
            watermark.erase(sd->station);
        } else {
            held.push(std::move(*m));
        }
        release();
    }
    release();
    threads.clear();
    const double wall = clock.wall_seconds();

    // --- summary ------------------------------------------------------------------
    Log::get().line("summary", "1", "{:<5} {:>8} {:>6} {:>7} {:>6} {:>5} {:>5} {:>6}   window ms: {:>5} {:>5} {:>5}   pick ms  mag ms",
                    "sta", "windows", "gaps", "stale", "detect", "picks", "mags", "noise", "mean", "p99", "max");
    std::uint64_t windows = 0;
    for (std::size_t i = 0; i < stations.size(); ++i) {
        const auto& s = procs[i]->stats();
        const auto w = percentiles(s.window_ms);
        const auto pk = percentiles(s.pick_ms);
        const auto mg = percentiles(s.magnitude_ms);
        std::uint64_t stale = 0, dropped = 0;
        for (auto& c : stations[i]->comp) {
            stale += c.stale.load();
            dropped += c.ring.dropped();
        }
        if (dropped) Log::get().line("warn", "33", "{}: {} samples dropped by the ring", stations[i]->code, dropped);
        windows += s.windows;
        Log::get().line("summary", "1", "{:<5} {:>8} {:>6} {:>7} {:>6} {:>5} {:>5} {:>6}   {:>15.2f} {:>5.2f} {:>5.2f}   {:>7.0f} {:>7.0f}",
                        stations[i]->code, s.windows, s.gap_windows, stale, s.detections, s.picks, s.magnitudes,
                        s.noise_windows, w.mean, w.p99, w.max, pk.mean, mg.mean);
    }
    const double stream_seconds = t_last - std::max(start, t_first);
    net.summary();
    Log::get().line("summary", "1", "{} windows in {:.1f} s wall: {:.0f} station-seconds per second ({:.0f}x real time per station)",
                    windows, wall, stream_seconds * static_cast<double>(stations.size()) / wall, stream_seconds / wall);
    return 0;
} catch (const std::exception& e) {
    std::println(stderr, "ayzek: {}", e.what());
    return 1;
}
