// Compares dsp::FilterPicker with Lomax's FilterPicker5 C function
// (tests/reference/filterpicker), fed in packets as in his
// picker_func_test_memory.c: the picks (time, uncertainty interval, polarity,
// strength, band) and the characteristic function of every sample. Data: a
// synthetic 20 minutes of noise, earthquakes, a spike and a step, and the 4
// minutes of DEMI of tools/export_stalta.py when that fixture exists.

#include "filterpicker.hpp"
#include "fixture.hpp"
#include "reference/filterpicker/fp_reference.h"

#include <array>
#include <cmath>
#include <filesystem>
#include <numbers>
#include <random>
#include <string>
#include <tuple>
#include <vector>

using namespace ayzek;

namespace {

// Raw counts, as integers in floats (the C function takes floats): an offset,
// noise whose level drifts, and P/S wave trains of several frequencies and
// sizes, a one-sample spike and a step.
std::vector<float> synthetic() {
  constexpr double fs = 100.0;
  const std::size_t n = static_cast<std::size_t>(1200 * fs);
  std::mt19937 rng(20261026);
  std::normal_distribution<double> noise(0.0, 1.0);
  std::vector<double> x(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double t = static_cast<double>(i) / fs;
    x[i] = 12345.0 + 40.0 * (1.0 + 0.5 * std::sin(t / 90.0)) * noise(rng);
  }
  struct Quake {
    double t, f, amp, sp; // P time, frequency, P amplitude, S - P
  };
  for (const auto &q :
       {Quake{60, 8.0, 400, 4.0}, Quake{150, 3.0, 5000, 12.0},
        Quake{260, 15.0, 150, 2.5}, Quake{400, 1.5, 60000, 30.0},
        Quake{520, 6.0, 120, 6.0}, Quake{700, 20.0, 2000, 1.5},
        Quake{900, 4.0, 900, 9.0}}) {
    for (const auto &[t0, amp, f, decay] :
         {std::tuple{q.t, q.amp, q.f, 3.0},
          std::tuple{q.t + q.sp, 3.0 * q.amp, 0.6 * q.f, 6.0}}) {
      for (std::size_t i = static_cast<std::size_t>(t0 * fs); i < n; ++i) {
        const double s = static_cast<double>(i) / fs - t0;
        if (s > 8 * decay)
          break;
        x[i] += amp * (1.0 - std::exp(-s * 20.0)) * std::exp(-s / decay) *
                std::sin(2 * std::numbers::pi * f * s);
      }
    }
  }
  x[static_cast<std::size_t>(330 * fs)] += 1e5; // spike
  for (std::size_t i = static_cast<std::size_t>(800 * fs); i < n; ++i)
    x[i] += 2000.0; // step
  std::vector<float> out(n);
  for (std::size_t i = 0; i < n; ++i)
    out[i] = static_cast<float>(std::round(x[i]));
  return out;
}

struct Case {
  const char *name;
  dsp::FilterPickerConfig cfg;
  int packet;
};

void compare(const std::string &data, const std::vector<float> &x,
             const Case &c, bool expect_picks) {
  const auto &k = c.cfg;
  std::vector<std::array<double, 7>> ref(4096);
  std::vector<double> ref_cf(x.size());
  const int nref = fp_reference(
      1.0 / k.fs, x.data(), static_cast<int>(x.size()), c.packet,
      k.filter_window, k.long_term_window, k.threshold1, k.threshold2, k.t_up,
      reinterpret_cast<double (*)[7]>(ref.data()), static_cast<int>(ref.size()),
      ref_cf.data());
  CHECK(nref >= 0);

  dsp::FilterPicker fp(k);
  std::vector<dsp::FilterPick> picks;
  std::vector<std::int64_t> declared;
  double cf_err = 0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (auto p = fp.step(x[i])) {
      picks.push_back(*p);
      declared.push_back(static_cast<std::int64_t>(i));
    }
    // The C function returns the CF as floats.
    const double ours = static_cast<float>(fp.cf());
    cf_err = std::max(cf_err, std::abs(ours - ref_cf[i]) /
                                  std::max(1.0, std::abs(ref_cf[i])));
  }

  std::size_t mismatched = 0;
  double amp_err = 0, weight_err = 0;
  CHECK(picks.size() == static_cast<std::size_t>(nref));
  for (std::size_t i = 0; i < picks.size(); ++i) {
    const auto &p = picks[i];
    const auto &r = ref[i];
    if (static_cast<double>(p.time) != r[0] ||
        static_cast<double>(p.begin) != r[1] ||
        static_cast<double>(p.end) != r[2] || p.polarity != r[3] ||
        p.period != r[6])
      ++mismatched;
    amp_err = std::max(amp_err, std::abs(p.amplitude - r[5]) / std::abs(r[5]));
    weight_err = std::max(weight_err, std::abs(p.polarity_weight - r[4]));
    // Declared after its trigger, never before the pick it reports.
    CHECK(p.begin <= p.time && p.time <= p.end);
    CHECK(p.time <= p.trigger && p.trigger < declared[i]);
  }
  std::println("  {} {} ({} bands, packets of {}): {} picks, {} differ; "
               "strength {:.1e}, polarity weight {:.1e}, CF {:.1e}",
               data, c.name, fp.bands(), c.packet, picks.size(), mismatched,
               amp_err, weight_err, cf_err);
  CHECK(mismatched == 0);
  CHECK(amp_err < 1e-9 && weight_err < 1e-9);
  CHECK(cf_err < 1e-5); // float rounding of the reference's CF output
  if (expect_picks)
    CHECK(!picks.empty());

  // reset() starts again from scratch.
  fp.reset();
  std::size_t again = 0;
  for (const float v : x)
    again += fp.step(v).has_value();
  CHECK(again == picks.size());
}

} // namespace

int main(int argc, char **argv) {
  const auto root = fixture::root(argc, argv);
  const Case cases[] = {
      {"Lomax et al. defaults", {}, 1000},
      {"Earthworm windows",
       {.filter_window = 6.0, .long_term_window = 10.0, .t_up = 0.4},
       37},
      {"thresholds 5/5", {.threshold1 = 5.0, .threshold2 = 5.0}, 1000},
  };

  const auto syn = synthetic();
  for (const auto &c : cases)
    compare("synthetic", syn, c, true);

  // Real data, if tools/export_stalta.py has been run.
  const auto path = root + "/data/fixtures/stalta.ayzw";
  if (std::filesystem::exists(path)) {
    const auto fx = Weights::load(path);
    const auto &raw = fx.at("raw");
    const auto v = raw.f64();
    std::vector<float> z(raw.shape[0]);
    for (std::size_t i = 0; i < z.size(); ++i)
      z[i] = static_cast<float>(v[i * 3]);
    for (const auto &c : cases)
      compare("DEMI", z, c, false);
  } else {
    std::println("  DEMI: {} not found, synthetic data only", path);
  }
  std::println("ok");
}
