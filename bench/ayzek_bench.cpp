// ayzek_bench: how long each part of the pipeline takes on this machine, and how
// many stations it can therefore run in real time. docs/impl/10-benchmarks.md.
//
//   ayzek_bench [--models DIR] [--group stages,layers,shapes,capacity]
//               [--layer SPEC]... [--min-time S] [--threads N] [--quick]
//               [--json FILE] [--allow-debug]
//
// Groups:
//   stages    the pipeline's per-window, per-pick and per-estimate work, with
//             the shipped models
//   layers    each network layer at the shape a shipped model runs it at
//   shapes    layers at the shapes given with --layer, to price a candidate
//             architecture before it exists in C++, e.g.
//               --layer lstm:T=125,in=96,h=64 --layer attn:T=125,d=128,heads=4
//   capacity  detector throughput on 1..N threads, the number of stations that
//             fits, and the magnitude estimate's latency while they run
//
// Inputs are synthetic and seeded: inference time here does not depend on the
// sample values, so no waveform data is needed.

#include "bench_version.hpp"
#include "dsp.hpp"
#include "magnitude.hpp"
#include "models.hpp"
#include "nn.hpp"
#include "simd.hpp"
#include "stalta.hpp"
#include "weights.hpp"

#include <algorithm>
#include <atomic>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <memory>
#include <numeric>
#include <print>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#ifndef _WIN32
#include <unistd.h>
#endif

using namespace ayzek;
using Clock = std::chrono::steady_clock;

namespace {

// --- options --------------------------------------------------------------

struct Options {
  std::string models = "models";
  std::set<std::string> groups{"stages", "layers", "capacity"};
  std::vector<std::string> layers;
  double min_time = 1.0; // seconds of timed calls per benchmark
  std::size_t threads = 0;
  bool quick = false, allow_debug = false;
  std::string json;
};

[[noreturn]] void usage(const char *why) {
  std::println(stderr, "ayzek_bench: {}", why);
  std::println(stderr,
               "usage: ayzek_bench [--models DIR] [--group G[,G...]] "
               "[--layer SPEC]... [--min-time S] [--threads N] [--quick] "
               "[--json FILE] [--allow-debug]");
  std::exit(2);
}

Options parse(int argc, char **argv) {
  Options o;
  bool groups_given = false;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc)
        usage(("missing value for " + a).c_str());
      return argv[++i];
    };
    if (a == "--models")
      o.models = next();
    else if (a == "--group") {
      if (!groups_given)
        o.groups.clear();
      groups_given = true;
      std::stringstream ss(next());
      for (std::string g; std::getline(ss, g, ',');)
        o.groups.insert(g);
    } else if (a == "--layer") {
      o.layers.push_back(next());
      if (!groups_given)
        o.groups = {"shapes"};
    } else if (a == "--min-time")
      o.min_time = std::stod(next());
    else if (a == "--threads")
      o.threads = std::stoul(next());
    else if (a == "--json")
      o.json = next();
    else if (a == "--quick")
      o.quick = true;
    else if (a == "--allow-debug")
      o.allow_debug = true;
    else if (a == "-h" || a == "--help")
      usage("benchmarks the pipeline; see the comment at the top of "
            "bench/ayzek_bench.cpp");
    else
      usage(("unknown option " + a).c_str());
  }
  for (const auto &g : o.groups)
    if (g != "stages" && g != "layers" && g != "shapes" && g != "capacity")
      usage(("unknown group " + g).c_str());
  if (o.groups.contains("shapes") && o.layers.empty())
    usage("--group shapes needs at least one --layer");
  if (o.quick)
    o.min_time = 0.05;
  return o;
}

// --- machine --------------------------------------------------------------

constexpr bool kOptimized =
#ifdef __OPTIMIZE__
    true;
#else
    false;
#endif

std::string first_line_with(const std::string &path, const std::string &key) {
  std::ifstream f(path);
  for (std::string line; std::getline(f, line);)
    if (line.starts_with(key)) {
      const auto c = line.find(':');
      auto v = c == std::string::npos ? line : line.substr(c + 1);
      v.erase(0, v.find_first_not_of(" \t"));
      return v;
    }
  return {};
}

std::string cpu_name() {
  // x86: "model name". A Raspberry Pi reports no model name per core; its
  // board is in the device tree.
  if (auto m = first_line_with("/proc/cpuinfo", "model name"); !m.empty())
    return m;
  std::ifstream dt("/proc/device-tree/model");
  std::string board;
  std::getline(dt, board, '\0');
  if (!board.empty())
    return board;
  return first_line_with("/proc/cpuinfo", "Hardware");
}

std::string hostname() {
#ifdef _WIN32
  // No gethostname without winsock; the machine name is in the environment.
  const char *name = std::getenv("COMPUTERNAME");
  return name ? name : "windows";
#else
  char buf[256] = {};
  gethostname(buf, sizeof buf - 1);
  return buf;
#endif
}

// --- timing ---------------------------------------------------------------

struct Result {
  std::string group, name, shape;
  std::size_t iters = 0;
  double min_ms = 0, p50_ms = 0, p90_ms = 0, p99_ms = 0, mean_ms = 0;
  double macs = 0; // multiply-accumulates per call, 0 if not counted
};

std::vector<Result> g_results;
volatile float g_sink = 0; // keeps results observable to the optimiser

double pct(std::vector<double> v, double q) {
  std::ranges::sort(v);
  const auto i = static_cast<std::size_t>(q * static_cast<double>(v.size() - 1) + 0.5);
  return v[std::min(i, v.size() - 1)];
}

// Times `f` per call: warm-up, then calls until `min_time` seconds and at least
// ten samples. `reps` calls make one sample, for calls too short to time singly.
template <class F>
Result time_it(const Options &o, std::string group, std::string name,
               std::string shape, F &&f, double macs = 0, std::size_t reps = 1) {
  for (int i = 0; i < 3; ++i)
    f();
  std::vector<double> ms;
  const auto t_end = Clock::now() + std::chrono::duration<double>(o.min_time);
  while (ms.size() < 10 || (Clock::now() < t_end && ms.size() < 1'000'000)) {
    const auto t0 = Clock::now();
    for (std::size_t r = 0; r < reps; ++r)
      f();
    ms.push_back(std::chrono::duration<double, std::milli>(Clock::now() - t0).count() /
                 static_cast<double>(reps));
  }
  Result r{std::move(group), std::move(name), std::move(shape), ms.size() * reps,
           *std::ranges::min_element(ms), pct(ms, 0.5), pct(ms, 0.9), pct(ms, 0.99),
           std::accumulate(ms.begin(), ms.end(), 0.0) / static_cast<double>(ms.size()),
           macs};
  std::println("  {:<34} {:>10.4f} {:>10.4f} {:>10.4f} {:>9}{}", r.name, r.p50_ms,
               r.p90_ms, r.p99_ms, r.iters,
               macs > 0 ? std::format("  {:7.2f} GMAC/s", macs / (r.p50_ms * 1e6)) : "");
  g_results.push_back(r);
  return r;
}

void header(const char *title) {
  std::println("\n{}", title);
  std::println("  {:<34} {:>10} {:>10} {:>10} {:>9}", "", "p50 ms", "p90 ms", "p99 ms",
               "calls");
}

// --- synthetic inputs -----------------------------------------------------

std::mt19937 g_rng(42);

std::vector<double> raw_counts(std::size_t n, std::size_t channels) {
  // Broadband noise of a few thousand counts plus an offset, like a quiet
  // HH channel.
  std::normal_distribution<double> d(0.0, 2000.0);
  std::vector<double> x(n * channels);
  for (auto &v : x)
    v = 1.5e4 + d(g_rng);
  return x;
}

std::vector<float> randn(std::size_t n, float scale = 1.0f) {
  std::normal_distribution<float> d(0.0f, scale);
  std::vector<float> x(n);
  for (auto &v : x)
    v = d(g_rng);
  return x;
}

// --- synthetic layers, for any shape --------------------------------------

nn::Conv1d conv1d(std::size_t cin, std::size_t cout, std::size_t k, std::size_t s) {
  nn::Conv1d c;
  c.cin = cin, c.cout = cout, c.k = k, c.stride = s, c.pad = k / 2;
  c.w = randn(cout * cin * k, 0.1f);
  c.col.resize(cin * k);
  return c;
}

nn::Conv2d conv2d(std::size_t cin, std::size_t cout, std::size_t k, std::size_t s) {
  nn::Conv2d c;
  c.cin = cin, c.cout = cout, c.k = k, c.stride = s, c.pad = k / 2;
  c.w = randn(cout * cin * k * k, 0.1f);
  c.col.resize(cin * k * k);
  return c;
}

nn::LstmDirection lstm_dir(std::size_t in, std::size_t h) {
  nn::LstmDirection d;
  d.in = in, d.hidden = h;
  d.w_ih = randn(4 * h * in, 0.1f);
  d.w_hh = randn(4 * h * h, 0.1f);
  d.bias = randn(4 * h, 0.1f);
  d.gates.resize(4 * h);
  d.h.resize(h);
  d.c.resize(h);
  return d;
}

nn::BiLstm bilstm(std::size_t in, std::size_t h) {
  nn::BiLstm b;
  b.fwd = lstm_dir(in, h);
  b.bwd = lstm_dir(in, h);
  return b;
}

nn::SelfAttention attention(std::size_t d, std::size_t heads) {
  nn::SelfAttention a;
  a.d = d, a.heads = heads;
  a.in_w = randn(3 * d * d, 0.1f);
  a.in_b = randn(3 * d, 0.1f);
  a.out_w = randn(d * d, 0.1f);
  a.out_b = randn(d, 0.1f);
  return a;
}

nn::Linear linear(std::size_t in, std::size_t out) {
  nn::Linear l;
  l.in = in, l.out = out;
  l.w = randn(out * in, 0.1f);
  l.b = randn(out, 0.1f);
  return l;
}

// Multiply-accumulates per forward call.
double macs_conv1d(std::size_t L, std::size_t cin, std::size_t cout, std::size_t k,
                   std::size_t s) {
  return static_cast<double>((L + 2 * (k / 2) - k) / s + 1) * cout * cin * k;
}
double macs_conv2d(std::size_t H, std::size_t W, std::size_t cin, std::size_t cout,
                   std::size_t k, std::size_t s) {
  const auto o = [&](std::size_t n) { return (n + 2 * (k / 2) - k) / s + 1; };
  return static_cast<double>(o(H) * o(W)) * cout * cin * k * k;
}
double macs_lstm(std::size_t T, std::size_t in, std::size_t h) {
  return 2.0 * T * 4 * h * (in + h);
}
double macs_attn(std::size_t T, std::size_t d) {
  return 4.0 * T * d * d + 2.0 * T * T * d;
}

// One layer benchmark from a kind and its sizes.
void bench_layer(const Options &o, const std::string &group, const std::string &name,
                 const std::string &kind, std::map<std::string, std::size_t> p) {
  auto need = [&](const char *k) {
    if (!p.contains(k))
      usage(std::format("layer {} needs {}=", kind, k).c_str());
    return p[k];
  };
  std::string shape;
  for (const auto &[k, v] : p)
    shape += std::format("{}{}={}", shape.empty() ? "" : ",", k, v);
  if (kind == "conv1d") {
    const auto L = need("L"), cin = need("cin"), cout = need("cout"), k = need("k");
    const auto s = p.contains("s") ? p["s"] : 1;
    auto c = conv1d(cin, cout, k, s);
    auto x = randn(cin * L);
    std::vector<float> y(cout * c.out_len(L));
    time_it(o, group, name, "conv1d:" + shape, [&] {
      c.forward(x.data(), L, y.data());
      g_sink = y[0];
    }, macs_conv1d(L, cin, cout, k, s));
  } else if (kind == "conv2d") {
    const auto H = need("H"), W = need("W"), cin = need("cin"), cout = need("cout"),
               k = need("k");
    const auto s = p.contains("s") ? p["s"] : 1;
    auto c = conv2d(cin, cout, k, s);
    auto x = randn(cin * H * W);
    std::vector<float> y(cout * c.out_len(H) * c.out_len(W));
    time_it(o, group, name, "conv2d:" + shape, [&] {
      c.forward(x.data(), H, W, y.data());
      g_sink = y[0];
    }, macs_conv2d(H, W, cin, cout, k, s));
  } else if (kind == "lstm") {
    const auto T = need("T"), in = need("in"), h = need("h");
    auto l = bilstm(in, h);
    auto x = randn(T * in);
    std::vector<float> y(T * 2 * h);
    time_it(o, group, name, "lstm:" + shape, [&] {
      l.forward(x.data(), T, y.data());
      g_sink = y[0];
    }, macs_lstm(T, in, h));
  } else if (kind == "attn") {
    const auto T = need("T"), d = need("d");
    const auto heads = p.contains("heads") ? p["heads"] : 4;
    auto a = attention(d, heads);
    auto x = randn(T * d);
    std::vector<float> y(T * d);
    time_it(o, group, name, "attn:" + shape, [&] {
      a.forward(x.data(), T, y.data());
      g_sink = y[0];
    }, macs_attn(T, d));
  } else if (kind == "linear") {
    const auto T = need("T"), in = need("in"), out = need("out");
    auto l = linear(in, out);
    auto x = randn(T * in);
    std::vector<float> y(T * out);
    time_it(o, group, name, "linear:" + shape, [&] {
      l.forward(x.data(), T, y.data());
      g_sink = y[0];
    }, static_cast<double>(T) * in * out);
  } else {
    usage(("unknown layer kind " + kind + " (conv1d, conv2d, lstm, attn, linear)").c_str());
  }
}

// "lstm:T=125,in=96,h=64" -> kind and sizes.
void bench_spec(const Options &o, const std::string &spec) {
  const auto colon = spec.find(':');
  if (colon == std::string::npos)
    usage(("layer spec needs KIND:K=V,...: " + spec).c_str());
  std::map<std::string, std::size_t> p;
  std::stringstream ss(spec.substr(colon + 1));
  for (std::string kv; std::getline(ss, kv, ',');) {
    const auto eq = kv.find('=');
    std::size_t v = 0;
    if (eq == std::string::npos ||
        std::from_chars(kv.data() + eq + 1, kv.data() + kv.size(), v).ec != std::errc{})
      usage(("bad size in layer spec: " + kv).c_str());
    p[kv.substr(0, eq)] = v;
  }
  bench_layer(o, "shapes", spec, spec.substr(0, colon), p);
}

// --- the shipped models ---------------------------------------------------

struct Models {
  std::vector<Weights> detectors, magnitude;
  std::unique_ptr<Weights> picker, bandpass;
};

Models load_models(const std::string &dir) {
  Models m;
  for (int seed : {42, 43, 44})
    m.detectors.push_back(Weights::load(std::format("{}/detector_s{}.ayzw", dir, seed)));
  for (int p = 0; p < 3; ++p)
    m.magnitude.push_back(Weights::load(std::format("{}/magnitude_p{}.ayzw", dir, p)));
  m.picker = std::make_unique<Weights>(Weights::load(dir + "/spicker.ayzw"));
  m.bandpass = std::make_unique<Weights>(Weights::load(dir + "/bandpass.ayzw"));
  return m;
}

// A noise baseline with enough windows to be used, as a station has after its
// first minutes.
StationNoise ready_noise(const dsp::Bandpass &bp) {
  NoiseBaseline nb(bp);
  for (int i = 0; i < 8; ++i)
    nb.add(raw_counts(kMagWindow, 3));
  return nb.noise();
}

// Work per detector window, as the station processor does it every 0.5 s:
// condition 6 s, run the ensemble, and feed the STA/LTA the 50 new samples.
struct WindowWork {
  dsp::Conditioner cond;
  DetectorEnsemble ens;
  dsp::StaLta sl;
  std::vector<double> raw;
  std::vector<float> x;
  std::size_t step = 0;

  WindowWork(const dsp::Bandpass &bp, const std::vector<Weights> &seeds)
      : cond(bp, Detector::kWindow), ens(seeds), sl(dsp::StaLtaConfig{}),
        raw(raw_counts(Detector::kWindow, 3)), x(Detector::kWindow * 3) {}

  float operator()() {
    cond.condition(raw, 3, x);
    const float p = ens.probability(x);
    for (std::size_t i = 0; i < 50; ++i) {
      const std::size_t t = (step + i) % Detector::kWindow;
      const std::array<double, 3> s{raw[t * 3], raw[t * 3 + 1], raw[t * 3 + 2]};
      g_sink = static_cast<float>(sl.step(s));
    }
    step += 50;
    return p;
  }
};

void bench_stages(const Options &o, const Models &m) {
  header("stages (shipped models)");
  const auto bp = dsp::Bandpass::load(*m.bandpass);

  // Detector window.
  {
    dsp::Conditioner cond(bp, Detector::kWindow);
    auto raw = raw_counts(Detector::kWindow, 3);
    std::vector<float> x(Detector::kWindow * 3);
    time_it(o, "stages", "condition 6 s window", "600x3", [&] {
      cond.condition(raw, 3, x);
      g_sink = x[0];
    });
    Detector one(m.detectors[0]);
    std::vector<float> in(x.size());
    std::ranges::transform(x, in.begin(), [](float v) { return std::asinh(v); });
    time_it(o, "stages", "detector, one model", "600x3", [&] { g_sink = one.logit(in); });
    DetectorEnsemble ens(m.detectors);
    time_it(o, "stages", "detector ensemble (3 models)", "600x3",
            [&] { g_sink = ens.probability(x); });
  }
  {
    dsp::StaLta sl(dsp::StaLtaConfig{});
    auto raw = raw_counts(Detector::kWindow, 3);
    std::size_t t = 0;
    time_it(o, "stages", "STA/LTA, 50 samples", "50x3", [&] {
      for (int i = 0; i < 50; ++i, t = (t + 1) % Detector::kWindow) {
        const std::array<double, 3> s{raw[t * 3], raw[t * 3 + 1], raw[t * 3 + 2]};
        g_sink = static_cast<float>(sl.step(s));
      }
    }, 0, 100);
  }
  {
    WindowWork w(bp, m.detectors);
    time_it(o, "stages", "per-window total", "every 0.5 s", [&] { g_sink = w(); });
  }

  // Picker: 60 s, planar.
  {
    dsp::Conditioner cond(bp, Picker::kWindow);
    auto raw = raw_counts(Picker::kWindow, 3);
    std::vector<float> x(Picker::kWindow * 3);
    time_it(o, "stages", "condition 60 s window", "6000x3", [&] {
      cond.condition_planar(raw, 3, x);
      g_sink = x[0];
    });
    Picker pk(*m.picker);
    time_it(o, "stages", "picker", "3x6000",
            [&] { g_sink = static_cast<float>(pk.pick(x).p_seconds); });
  }

  // Magnitude: the estimate includes conditioning and the spectrogram.
  {
    const auto noise = ready_noise(bp);
    auto raw = raw_counts(kMagWindow, 3);
    MagnitudeEstimator est(m.magnitude, bp);
    time_it(o, "stages", "magnitude estimate (3 models)", "1000x3",
            [&] { g_sink = est.estimate(raw, noise); });
    MagnitudeRegressor one(m.magnitude[0]);
    const auto seq = est.seq, img = est.img;
    time_it(o, "stages", "magnitude, one model", "1000x3 + 3x65x32",
            [&] { g_sink = one.predict(seq, img); });
    NoiseBaseline nb(bp);
    time_it(o, "stages", "noise baseline update", "1000x3", [&] {
      nb.add(raw);
      g_sink = static_cast<float>(nb.noise().windows);
    });
  }
}

void bench_layers(const Options &o) {
  header("layers (shapes of the shipped models)");
  using P = std::map<std::string, std::size_t>;
  const std::vector<std::tuple<std::string, std::string, P>> L = {
      {"detector conv 1", "conv1d", P{{"L", 600}, {"cin", 3}, {"cout", 24}, {"k", 7}, {"s", 2}}},
      {"detector conv 2", "conv1d", P{{"L", 300}, {"cin", 24}, {"cout", 48}, {"k", 5}, {"s", 2}}},
      {"detector conv 3", "conv1d", P{{"L", 150}, {"cin", 48}, {"cout", 96}, {"k", 5}, {"s", 2}}},
      {"detector BiLSTM", "lstm", P{{"T", 75}, {"in", 96}, {"h", 48}}},
      {"detector attention", "attn", P{{"T", 75}, {"d", 96}, {"heads", 4}}},
      {"picker conv 1", "conv1d", P{{"L", 6000}, {"cin", 3}, {"cout", 32}, {"k", 9}, {"s", 2}}},
      {"picker conv 2", "conv1d", P{{"L", 3000}, {"cin", 32}, {"cout", 64}, {"k", 9}, {"s", 2}}},
      {"picker conv 3", "conv1d", P{{"L", 1500}, {"cin", 64}, {"cout", 64}, {"k", 7}, {"s", 2}}},
      {"picker BiLSTM", "lstm", P{{"T", 250}, {"in", 64}, {"h", 32}}},
      {"picker attention", "attn", P{{"T", 250}, {"d", 64}, {"heads", 4}}},
      {"magnitude BiLSTM (1D)", "lstm", P{{"T", 1000}, {"in", 3}, {"h", 64}}},
      {"magnitude attention (1D)", "attn", P{{"T", 1000}, {"d", 128}, {"heads", 4}}},
      {"magnitude conv 1 (2D)", "conv2d", P{{"H", 65}, {"W", 32}, {"cin", 3}, {"cout", 32}, {"k", 3}, {"s", 1}}},
      {"magnitude conv 2 (2D)", "conv2d", P{{"H", 65}, {"W", 32}, {"cin", 32}, {"cout", 64}, {"k", 3}, {"s", 2}}},
      {"magnitude conv 3 (2D)", "conv2d", P{{"H", 33}, {"W", 16}, {"cin", 64}, {"cout", 128}, {"k", 3}, {"s", 2}}},
  };
  for (const auto &[name, kind, p] : L)
    bench_layer(o, "layers", name, kind, p);
}

void bench_shapes(const Options &o) {
  header("shapes (--layer)");
  for (const auto &s : o.layers)
    bench_spec(o, s);
}

// Detector throughput on 1..N threads, each a station with its own models.
// The stations that fit are the windows per second the machine sustains,
// divided by the two windows each station needs per second.
struct Capacity {
  std::size_t threads = 0;
  double windows_per_s = 0, stations = 0;
};
std::vector<Capacity> g_capacity;
double g_mag_idle_p50 = 0, g_mag_loaded_p50 = 0, g_mag_loaded_p99 = 0;

void bench_capacity(const Options &o, const Models &m) {
  const auto bp = dsp::Bandpass::load(*m.bandpass);
  const std::size_t hw = std::max(1u, std::thread::hardware_concurrency());
  const std::size_t n = o.threads ? o.threads : hw;
  const double secs = std::max(o.min_time, o.quick ? 0.1 : 2.0);
  std::println("\ncapacity (detector windows per second, all threads together)");
  std::println("  {:>7} {:>12} {:>14} {:>16}", "threads", "windows/s", "per thread",
               "stations (0.5 s)");
  std::vector<std::size_t> counts;
  for (std::size_t t = 1; t <= n; t = (t < 4 ? t + 1 : t * 2)) {
    counts.push_back(t);
    if (t != n && t * 2 > n && t >= 4)
      counts.push_back(n);
  }
  counts.erase(std::unique(counts.begin(), counts.end()), counts.end());
  for (std::size_t t : counts) {
    std::atomic<bool> go{false}, stop{false};
    std::vector<std::size_t> done(t, 0);
    std::vector<std::thread> pool;
    for (std::size_t i = 0; i < t; ++i)
      pool.emplace_back([&, i] {
        WindowWork w(bp, m.detectors);
        w(); // warm up before the clock starts
        while (!go.load())
          std::this_thread::yield();
        while (!stop.load(std::memory_order_relaxed)) {
          g_sink = w();
          ++done[i];
        }
      });
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    const auto t0 = Clock::now();
    go = true;
    std::this_thread::sleep_for(std::chrono::duration<double>(secs));
    stop = true;
    for (auto &th : pool)
      th.join();
    const double el = std::chrono::duration<double>(Clock::now() - t0).count();
    const double wps = static_cast<double>(std::accumulate(done.begin(), done.end(), 0ul)) / el;
    g_capacity.push_back({t, wps, wps / 2.0});
    std::println("  {:>7} {:>12.0f} {:>14.0f} {:>16.0f}", t, wps, wps / static_cast<double>(t),
                 wps / 2.0);
  }

  // The magnitude estimate runs on a station's own thread and holds up that
  // station while it runs, so its latency with the other cores busy is the
  // figure that decides how late the magnitude arrives.
  const auto noise = ready_noise(bp);
  auto raw = raw_counts(kMagWindow, 3);
  MagnitudeEstimator est(m.magnitude, bp);
  auto sample = [&](std::size_t k) {
    std::vector<double> ms;
    for (std::size_t i = 0; i < k; ++i) {
      const auto t0 = Clock::now();
      g_sink = est.estimate(raw, noise);
      ms.push_back(std::chrono::duration<double, std::milli>(Clock::now() - t0).count());
    }
    return ms;
  };
  const std::size_t k = o.quick ? 3 : 20;
  g_mag_idle_p50 = pct(sample(k), 0.5);
  std::atomic<bool> stop{false};
  std::vector<std::thread> load;
  for (std::size_t i = 0; i + 1 < n; ++i)
    load.emplace_back([&] {
      WindowWork w(bp, m.detectors);
      while (!stop.load(std::memory_order_relaxed))
        g_sink = w();
    });
  const auto loaded = sample(k);
  stop = true;
  for (auto &th : load)
    th.join();
  g_mag_loaded_p50 = pct(loaded, 0.5);
  g_mag_loaded_p99 = pct(loaded, 0.99);
  std::println("\nmagnitude estimate (3 models): {:.1f} ms alone, {:.1f} ms (p99 {:.1f}) "
               "with {} detector threads running",
               g_mag_idle_p50, g_mag_loaded_p50, g_mag_loaded_p99, n - 1);
}

// --- output ---------------------------------------------------------------

std::string esc(const std::string &s) {
  std::string o;
  for (char c : s) {
    if (c == '"' || c == '\\')
      o += '\\';
    o += c;
  }
  return o;
}

void write_json(const Options &o) {
  std::ofstream f(o.json);
  if (!f) {
    std::println(stderr, "ayzek_bench: cannot write {}", o.json);
    std::exit(1);
  }
  const auto now = std::chrono::system_clock::now();
  f << "{\n";
  f << std::format("  \"host\": \"{}\",\n  \"cpu\": \"{}\",\n  \"cores\": {},\n",
                   esc(hostname()), esc(cpu_name()), std::thread::hardware_concurrency());
  f << std::format("  \"backend\": \"{}\",\n  \"optimized\": {},\n  \"compiler\": \"{}\",\n",
                   simd::kBackend, kOptimized, esc(__VERSION__));
  f << std::format("  \"version\": \"{}\",\n  \"models\": \"{}\",\n  \"date\": \"{:%FT%TZ}\",\n",
                   AYZEK_BENCH_VERSION, esc(o.models),
                   std::chrono::floor<std::chrono::seconds>(now));
  f << std::format("  \"min_time_s\": {},\n  \"quick\": {},\n", o.min_time, o.quick);
  f << "  \"results\": [\n";
  for (std::size_t i = 0; i < g_results.size(); ++i) {
    const auto &r = g_results[i];
    f << std::format("    {{\"group\": \"{}\", \"name\": \"{}\", \"shape\": \"{}\", "
                     "\"iters\": {}, \"min_ms\": {:.6g}, \"p50_ms\": {:.6g}, "
                     "\"p90_ms\": {:.6g}, \"p99_ms\": {:.6g}, \"mean_ms\": {:.6g}, "
                     "\"macs\": {:.6g}}}{}\n",
                     r.group, esc(r.name), esc(r.shape), r.iters, r.min_ms, r.p50_ms,
                     r.p90_ms, r.p99_ms, r.mean_ms, r.macs,
                     i + 1 < g_results.size() ? "," : "");
  }
  f << "  ],\n  \"capacity\": [\n";
  for (std::size_t i = 0; i < g_capacity.size(); ++i) {
    const auto &c = g_capacity[i];
    f << std::format("    {{\"threads\": {}, \"windows_per_s\": {:.6g}, \"stations\": {:.6g}}}{}\n",
                     c.threads, c.windows_per_s, c.stations,
                     i + 1 < g_capacity.size() ? "," : "");
  }
  f << std::format("  ],\n  \"magnitude_latency_ms\": {{\"idle_p50\": {:.6g}, "
                   "\"loaded_p50\": {:.6g}, \"loaded_p99\": {:.6g}}}\n}}\n",
                   g_mag_idle_p50, g_mag_loaded_p50, g_mag_loaded_p99);
  std::println("\nresults -> {}", o.json);
}

} // namespace

int main(int argc, char **argv) {
  const auto o = parse(argc, argv);
  std::println("ayzek_bench {}  |  {}  |  {} cores  |  {} backend  |  {}", AYZEK_BENCH_VERSION,
               cpu_name(), std::thread::hardware_concurrency(), simd::kBackend,
               kOptimized ? "optimised build" : "UNOPTIMISED build");
  if (!kOptimized && !o.allow_debug) {
    // A -O0 build runs the models about 17 times slower; its numbers would be
    // read as the pipeline's.
    std::println(stderr, "ayzek_bench: this build is not optimised; configure with "
                         "--buildtype=release, or pass --allow-debug to run anyway");
    return 1;
  }
  try {
    std::unique_ptr<Models> m;
    if (o.groups.contains("stages") || o.groups.contains("capacity"))
      m = std::make_unique<Models>(load_models(o.models));
    if (o.groups.contains("stages"))
      bench_stages(o, *m);
    if (o.groups.contains("layers"))
      bench_layers(o);
    if (o.groups.contains("shapes"))
      bench_shapes(o);
    if (o.groups.contains("capacity"))
      bench_capacity(o, *m);
  } catch (const std::exception &e) {
    std::println(stderr, "ayzek_bench: {}", e.what());
    return 1;
  }
  if (!o.json.empty())
    write_json(o);
  return 0;
}
