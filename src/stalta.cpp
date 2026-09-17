#include "stalta.hpp"

#include <cmath>
#include <limits>
#include <numbers>

namespace ayzek::dsp {

namespace {

// Quality factors of the two pole pairs of a 4th-order Butterworth filter,
// Q = 1 / (2 sin((2k - 1) pi / 8)), k = 1, 2.
constexpr std::array<double, 2> kQ = {1.0 / (2.0 * 0.3826834323650898), 1.0 / (2.0 * 0.9238795325112867)};

// Second-order low-pass (lowpass = true) or high-pass section by the bilinear
// transform with the cutoff prewarped (Bristow-Johnson's formulas).
Biquad section(double fc, double fs, double q, bool lowpass) {
    const double w0 = 2.0 * std::numbers::pi * fc / fs;
    const double c = std::cos(w0), alpha = std::sin(w0) / (2.0 * q);
    const double a0 = 1.0 + alpha;
    Biquad s;
    if (lowpass) {
        s.b0 = (1.0 - c) / 2.0 / a0;
        s.b1 = (1.0 - c) / a0;
    } else {
        s.b0 = (1.0 + c) / 2.0 / a0;
        s.b1 = -(1.0 + c) / a0;
    }
    s.b2 = s.b0;
    s.a1 = -2.0 * c / a0;
    s.a2 = (1.0 - alpha) / a0;
    return s;
}

}  // namespace

std::array<Biquad, 2> butter4_lowpass(double fc, double fs) {
    return {section(fc, fs, kQ[0], true), section(fc, fs, kQ[1], true)};
}

std::array<Biquad, 2> butter4_highpass(double fc, double fs) {
    return {section(fc, fs, kQ[0], false), section(fc, fs, kQ[1], false)};
}

StaLta::StaLta(const StaLtaConfig& cfg)
    : cfg_(cfg),
      csta_(1.0 / std::round(cfg.sta_seconds * cfg.fs)),
      clta_(1.0 / std::round(cfg.lta_seconds * cfg.fs)),
      nlta_(static_cast<std::size_t>(std::round(cfg.lta_seconds * cfg.fs))) {
    const auto hp = butter4_highpass(cfg.f_lo, cfg.fs);
    const auto lp = butter4_lowpass(cfg.f_hi, cfg.fs);
    for (auto& s : sections_) s = {hp[0], hp[1], lp[0], lp[1]};
    reset();
}

void StaLta::reset() {
    for (auto& comp : sections_)
        for (auto& s : comp) s.reset();
    offset_ = {};
    y_ = {};
    sta_ = 0;
    lta_ = std::numeric_limits<double>::min();
    n_ = 0;
}

double StaLta::step(std::span<const double, 3> x) {
    const std::size_t components = cfg_.three_component ? 3 : 1;
    double e = 0;
    for (std::size_t c = 0; c < components; ++c) {
        if (n_ == 0) offset_[c] = x[c];
        double v = x[c] - offset_[c];
        for (auto& s : sections_[c]) v = s.step(v);
        y_[c] = v;
        e += v * v;
    }
    // ObsPy starts the recursion at the second sample.
    if (n_ > 0) {
        sta_ = csta_ * e + (1.0 - csta_) * sta_;
        lta_ = clta_ * e + (1.0 - clta_) * lta_;
    }
    return n_++ < nlta_ ? 0.0 : sta_ / lta_;
}

}  // namespace ayzek::dsp
