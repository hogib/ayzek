// Compares dsp::StaLta with scipy's Butterworth filters and ObsPy's recursive
// STA/LTA on 4 minutes of DEMI (tools/export_stalta.py).

#include "fixture.hpp"
#include "stalta.hpp"

#include <algorithm>
#include <array>
#include <vector>

using namespace ayzek;

int main(int argc, char **argv) {
  const auto root = fixture::root(argc, argv);
  auto fx = fixture::load_or_skip(root + "/data/fixtures/stalta.ayzw");

  // Filter design: the pole pairs (section denominators) against scipy's. scipy
  // orders the sections differently and puts the whole gain in the first one,
  // so the numerators are checked through the filtered output below.
  const auto hp = dsp::butter4_highpass(1.0, 100.0),
             lp = dsp::butter4_lowpass(10.0, 100.0);
  auto sos_hp = fx.at("sos_hp").f64(),
       sos_lp = fx.at("sos_lp").f64(); // (2, 6): b0 b1 b2 a0 a1 a2
  double design = 0;
  for (const auto &[ours, ref] :
       {std::pair{hp, sos_hp}, std::pair{lp, sos_lp}}) {
    for (const auto &s : ours) {
      double best = 1e300;
      for (std::size_t k = 0; k < 2; ++k)
        best = std::min(
            best, std::max(std::abs(s.a1 - ref[k * 6 + 4] / ref[k * 6 + 3]),
                           std::abs(s.a2 - ref[k * 6 + 5] / ref[k * 6 + 3])));
      design = std::max(design, best);
    }
  }
  std::println("  section denominators: max abs difference {:.1e}", design);
  CHECK(design < 1e-12);

  const auto &raw = fx.at("raw");
  const std::size_t n = raw.shape[0];
  auto x = raw.f64();
  auto filtered = fx.at("filtered").f64();
  auto ratio_z = fx.at("ratio_z").f64(), ratio_zne = fx.at("ratio_zne").f64();

  for (const bool three : {false, true}) {
    dsp::StaLta sl({.sta_seconds = 1.0,
                    .lta_seconds = 30.0,
                    .f_lo = 1.0,
                    .f_hi = 10.0,
                    .three_component = three,
                    .fs = 100.0});
    const auto ref = three ? ratio_zne : ratio_z;
    double filt_err = 0, filt_scale = 0, ratio_err = 0, max_ratio = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const std::array<double, 3> s = {x[i * 3], x[i * 3 + 1], x[i * 3 + 2]};
      const double r = sl.step(s);
      for (std::size_t c = 0; c < (three ? 3u : 1u); ++c) {
        filt_err =
            std::max(filt_err, std::abs(sl.filtered(c) - filtered[i * 3 + c]));
        filt_scale = std::max(filt_scale, std::abs(filtered[i * 3 + c]));
      }
      ratio_err = std::max(ratio_err, std::abs(r - ref[i]) /
                                          std::max(1.0, std::abs(ref[i])));
      max_ratio = std::max(max_ratio, r);
    }
    std::println("  {}: filtered max abs error {:.1e} (peak {:.0f}), ratio max "
                 "relative error {:.1e}, peak ratio {:.1f}",
                 three ? "Z+N+E" : "Z", filt_err, filt_scale, ratio_err,
                 max_ratio);
    CHECK(filt_err < 1e-9 * filt_scale);
    CHECK(ratio_err < 1e-9);
    CHECK(max_ratio > 5); // the M4.9 P onset
  }
  std::println("STA/LTA matches scipy and ObsPy");
}
