// Window conditioning against scipy, on real DEMI windows.

#include "dsp.hpp"
#include "fixture.hpp"

#include <vector>

using namespace ayzek;

int main(int argc, char** argv) {
    const auto root = fixture::root(argc, argv);
    auto fx = fixture::load_or_skip(root + "/data/fixtures/dsp.ayzw");
    auto bp = dsp::Bandpass::load(fixture::load_or_skip(root + "/models/bandpass.ayzw"));

    const auto& raw = fx.at("raw");                 // (k, 600, 3) f64
    const std::size_t k = raw.shape[0], n = raw.shape[1], C = raw.shape[2];
    auto cleaned_ref = fx.at("cleaned").f64();
    auto std_ref = fx.at("standardized").f32();

    dsp::Conditioner cond(bp, n);
    std::vector<float> out(n * C);
    fixture::Diff worst_clean, worst_std;
    for (std::size_t w = 0; w < k; ++w) {
        cond.condition(raw.f64().subspan(w * n * C, n * C), C, out);
        auto dc = fixture::compare<double, double>(cond.cleaned, cleaned_ref.subspan(w * n * C, n * C), 1.0);
        auto ds = fixture::compare<float, float>(out, std_ref.subspan(w * n * C, n * C), 1.0);
        std::println("  window {}: cleaned max rel err {:.2e}, standardised max abs err {:.2e}", w, dc.rel, ds.abs);
        worst_clean.rel = std::max(worst_clean.rel, dc.rel);
        worst_std.abs = std::max(worst_std.abs, ds.abs);
    }
    // An 8th-order IIR in transfer-function form amplifies rounding, so double
    // agreement with scipy is ~1e-8 relative, not machine epsilon. What the
    // model sees, the float32 standardised window, agrees to float rounding.
    CHECK(worst_clean.rel < 1e-7);
    CHECK(worst_std.abs < 1e-5);
    std::println("dsp matches scipy: cleaned {:.1e} relative, standardised {:.1e} absolute", worst_clean.rel, worst_std.abs);
}
