#pragma once

// Window preprocessing, the same sequence of operations as the detector's
// training pipeline (archive_pipeline.archive.clean_block followed by the
// per-window standardisation in products/scan.py):
//
//   linear detrend -> mean removal -> 5% Hann taper at each end
//   -> zero-phase 4th-order Butterworth band-pass 1-45 Hz (scipy filtfilt)
//   -> subtract mean, divide by standard deviation (ddof 0)
//
// Computation is in double, as in scipy; the result is converted to float.

#include "weights.hpp"

#include <span>
#include <vector>

namespace ayzek::dsp {

struct Bandpass {
    std::vector<double> b, a, zi;                 // a[0] == 1, zi = lfilter_zi(b, a)
    static Bandpass load(const Weights& w);
};

void detrend_linear(std::span<double> x) noexcept;
void detrend_constant(std::span<double> x) noexcept;
void hann_taper(std::span<double> x) noexcept;   // scipy.signal.windows.hann(2k) halves, k = n // 20

// scipy.signal.filtfilt(b, a, x) with its defaults: odd extension of
// 3 * max(len(a), len(b)) samples, and lfilter_zi initial conditions.
class FiltFilt {
public:
    explicit FiltFilt(Bandpass bp);
    void apply(std::span<double> x);
private:
    Bandpass bp_;
    std::vector<double> ext_, state_;
    void lfilter(std::span<double> y, double x0);
};

// (x - mean) / max(std, 1e-12), converted to float.
void standardize(std::span<const double> x, std::span<float> out) noexcept;

// The full preprocessing of one window of raw counts, per component.
class Conditioner {
public:
    Conditioner(const Bandpass& bp, std::size_t window);
    // raw: (n, C) interleaved counts; out: (n, C) interleaved standardised.
    void condition(std::span<const double> raw, std::size_t channels, std::span<float> out);
    // Same, with (C, n) planar output, the picker's input layout.
    void condition_planar(std::span<const double> raw, std::size_t channels, std::span<float> out);
    // Filtered signal of the last call before standardisation, (n, C) interleaved.
    std::vector<double> cleaned;
private:
    FiltFilt filt_;
    std::size_t n_;
    std::vector<double> taper_;
};

}  // namespace ayzek::dsp
