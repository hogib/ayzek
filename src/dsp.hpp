#pragma once

// Window conditioning, operation for operation what the detector was trained
// on (archive_pipeline.archive.clean_block and the per-window standardisation
// in products/scan.py):
//
//   detrend linear -> detrend constant -> 5% Hann taper each end
//   -> zero-phase 4th-order Butterworth 1-45 Hz (scipy filtfilt, odd padding)
//   -> subtract mean, divide by std (ddof 0)
//
// All in double, like scipy, and cast to float only at the end. The filter
// recurrence is sequential by nature and costs about 1% of a window's compute;
// the window statistics use the NEON sums in simd.hpp.

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

// Mean/std standardisation, into float.
void standardize(std::span<const double> x, std::span<float> out) noexcept;

// The whole chain for one (n, 3) interleaved window of raw counts.
class Conditioner {
public:
    Conditioner(const Bandpass& bp, std::size_t window);
    // raw: (n, C) interleaved counts; out: (n, C) interleaved standardised.
    void condition(std::span<const double> raw, std::size_t channels, std::span<float> out);
    // Same, but writes (C, n) planar -- the layout the picker's convolutions take.
    void condition_planar(std::span<const double> raw, std::size_t channels, std::span<float> out);
    // Exposed for tests: the cleaned signal of the last channel processed per call.
    std::vector<double> cleaned;
private:
    FiltFilt filt_;
    std::size_t n_;
    std::vector<double> taper_;
};

}  // namespace ayzek::dsp
