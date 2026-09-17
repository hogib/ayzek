#include "dsp.hpp"

#include "simd.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace ayzek::dsp {

Bandpass Bandpass::load(const Weights &w) {
  Bandpass bp;
  auto b = w.at("b").f64(), a = w.at("a").f64(), zi = w.at("zi").f64();
  bp.b.assign(b.begin(), b.end());
  bp.a.assign(a.begin(), a.end());
  bp.zi.assign(zi.begin(), zi.end());
  if (bp.a.size() != bp.b.size() || bp.zi.size() + 1 != bp.a.size() ||
      bp.a[0] != 1.0)
    throw std::runtime_error(
        w.path() + ": filter coefficients are not a normalised b, a pair");
  return bp;
}

void detrend_linear(std::span<double> x) noexcept {
  // scipy fits c0 * t + c1, t = (i + 1) / N, by least squares. The fitted line
  // does not depend on the scaling of t, so the closed-form regression is used.
  const std::size_t n = x.size();
  if (n < 2)
    return;
  const double N = static_cast<double>(n);
  const double t_mean = (N + 1.0) / (2.0 * N);
  const double x_mean = simd::sum(x.data(), n) / N;
  double sxy = 0.0, sxx = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double t = static_cast<double>(i + 1) / N - t_mean;
    sxy += t * (x[i] - x_mean);
    sxx += t * t;
  }
  const double slope = sxy / sxx;
  for (std::size_t i = 0; i < n; ++i)
    x[i] -= x_mean + slope * (static_cast<double>(i + 1) / N - t_mean);
}

void detrend_constant(std::span<double> x) noexcept {
  const double m =
      simd::sum(x.data(), x.size()) / static_cast<double>(x.size());
  simd::affine(x.data(), 1.0, -m, x.size());
}

namespace {
std::vector<double> taper_vector(std::size_t n) {
  std::vector<double> t(n, 1.0);
  const std::size_t k = n / 20; // int(n * 0.05)
  if (k == 0)
    return t;
  const std::size_t M = 2 * k; // symmetric hann(M)
  for (std::size_t i = 0; i < k; ++i) {
    const double w =
        0.5 - 0.5 * std::cos(2.0 * std::numbers::pi * static_cast<double>(i) /
                             static_cast<double>(M - 1));
    const double w_end = 0.5 - 0.5 * std::cos(2.0 * std::numbers::pi *
                                              static_cast<double>(k + i) /
                                              static_cast<double>(M - 1));
    t[i] = w;
    t[n - k + i] = w_end;
  }
  return t;
}
} // namespace

void hann_taper(std::span<double> x) noexcept {
  const auto t = taper_vector(x.size());
  for (std::size_t i = 0; i < x.size(); ++i)
    x[i] *= t[i];
}

FiltFilt::FiltFilt(Bandpass bp) : bp_(std::move(bp)) {
  state_.resize(bp_.zi.size());
}

void FiltFilt::lfilter(std::span<double> y, double x0) {
  // Direct form II transposed, as scipy.signal.lfilter, starting from zi * x0.
  const auto &b = bp_.b;
  const auto &a = bp_.a;
  const std::size_t m = state_.size();
  for (std::size_t k = 0; k < m; ++k)
    state_[k] = bp_.zi[k] * x0;
  for (double &v : y) {
    const double in = v;
    const double out = b[0] * in + state_[0];
    for (std::size_t k = 0; k + 1 < m; ++k)
      state_[k] = b[k + 1] * in + state_[k + 1] - a[k + 1] * out;
    state_[m - 1] = b[m] * in - a[m] * out;
    v = out;
  }
}

void FiltFilt::apply(std::span<double> x) {
  const std::size_t n = x.size();
  const std::size_t edge = 3 * bp_.a.size();
  if (n <= edge)
    throw std::invalid_argument("filtfilt: input shorter than the padding");
  ext_.resize(n + 2 * edge);
  for (std::size_t i = 0; i < edge; ++i) {
    ext_[i] = 2.0 * x[0] - x[edge - i];
    ext_[n + edge + i] = 2.0 * x[n - 1] - x[n - 2 - i];
  }
  std::copy(x.begin(), x.end(),
            ext_.begin() + static_cast<std::ptrdiff_t>(edge));

  lfilter(ext_, ext_[0]);
  std::reverse(ext_.begin(), ext_.end());
  lfilter(ext_, ext_[0]);
  std::reverse(ext_.begin(), ext_.end());
  std::copy_n(ext_.begin() + static_cast<std::ptrdiff_t>(edge), n, x.begin());
}

void standardize(std::span<const double> x, std::span<float> out) noexcept {
  const std::size_t n = x.size();
  const double mean = simd::sum(x.data(), n) / static_cast<double>(n);
  const double ss = simd::dot(x.data(), x.data(), n) / static_cast<double>(n);
  const double sd = std::max(std::sqrt(std::max(ss - mean * mean, 0.0)), 1e-12);
  for (std::size_t i = 0; i < n; ++i)
    out[i] = static_cast<float>((x[i] - mean) / sd);
}

Conditioner::Conditioner(const Bandpass &bp, std::size_t window)
    : filt_(bp), n_(window), taper_(taper_vector(window)) {}

void Conditioner::condition(std::span<const double> raw, std::size_t C,
                            std::span<float> out) {
  thread_local std::vector<float> planar;
  planar.resize(C * n_);
  condition_planar(raw, C, planar);
  for (std::size_t c = 0; c < C; ++c)
    for (std::size_t t = 0; t < n_; ++t)
      out[t * C + c] = planar[c * n_ + t];
}

void Conditioner::condition_planar(std::span<const double> raw, std::size_t C,
                                   std::span<float> out) {
  if (raw.size() != n_ * C || out.size() != n_ * C)
    throw std::invalid_argument("conditioner: wrong window size");
  cleaned.resize(n_ * C);
  thread_local std::vector<double> x;
  x.resize(n_);
  for (std::size_t c = 0; c < C; ++c) {
    for (std::size_t t = 0; t < n_; ++t)
      x[t] = raw[t * C + c];
    detrend_linear(x);
    detrend_constant(x);
    simd::mul_d(taper_.data(), x.data(), n_);
    filt_.apply(x);
    for (std::size_t t = 0; t < n_; ++t)
      cleaned[t * C + c] = x[t];
    standardize(x, out.subspan(c * n_, n_));
  }
}

} // namespace ayzek::dsp
