#pragma once

// Recursive STA/LTA trigger on band-passed data, updated one sample at a time.
// It is the reference method the learned detector is compared with
// (docs/impl/09-sta-lta.md).
//
// Filter: 4th-order Butterworth high-pass at `f_lo` followed by a 4th-order
// Butterworth low-pass at `f_hi`, each as two biquads designed by the bilinear
// transform with prewarping (identical to scipy.signal.butter(4, f, fs=...,
// output="sos")). The first sample after a reset is subtracted from the input,
// so that the high-pass does not ring on the DC offset of the raw counts.
//
// Averages (Withers et al. 1998; obspy.signal.trigger.recursive_sta_lta):
//   sta_i = e_i / n_sta + (1 - 1/n_sta) sta_{i-1}
//   lta_i = e_i / n_lta + (1 - 1/n_lta) lta_{i-1}
// with e_i the squared filtered vertical component, or the sum of the squares of
// all three components. The ratio is 0 for the first n_lta samples after a
// reset, as in ObsPy.

#include <array>
#include <cstddef>
#include <span>

namespace ayzek::dsp {

struct Biquad {
    double b0 = 1, b1 = 0, b2 = 0, a1 = 0, a2 = 0;
    double z1 = 0, z2 = 0;   // transposed direct form II state

    double step(double x) {
        const double y = b0 * x + z1;
        z1 = b1 * x - a1 * y + z2;
        z2 = b2 * x - a2 * y;
        return y;
    }
    void reset() { z1 = z2 = 0; }
};

// 4th-order Butterworth sections for cutoff `fc` at sampling rate `fs`.
std::array<Biquad, 2> butter4_lowpass(double fc, double fs);
std::array<Biquad, 2> butter4_highpass(double fc, double fs);

struct StaLtaConfig {
    double sta_seconds = 1.0;
    double lta_seconds = 30.0;
    double f_lo = 1.0, f_hi = 10.0;   // pass band, Hz
    bool three_component = false;     // false: vertical only
    double fs = 100.0;
};

class StaLta {
public:
    explicit StaLta(const StaLtaConfig& cfg);
    // One sample per component, in the order Z, N, E. Returns the STA/LTA ratio.
    double step(std::span<const double, 3> x);
    void reset();
    // Filtered value of component c from the last step, for tests.
    [[nodiscard]] double filtered(std::size_t c) const { return y_[c]; }
    [[nodiscard]] const StaLtaConfig& config() const noexcept { return cfg_; }

private:
    StaLtaConfig cfg_;
    std::array<std::array<Biquad, 4>, 3> sections_;   // per component: high-pass, then low-pass
    std::array<double, 3> offset_{}, y_{};
    double csta_, clta_, sta_ = 0, lta_ = 0;
    std::size_t nlta_, n_ = 0;
};

}  // namespace ayzek::dsp
