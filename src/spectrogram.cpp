#include "spectrogram.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace ayzek::dsp {

Spectrogram::Spectrogram(std::size_t n_fft, std::size_t hop,
                         std::size_t samples)
    : n_fft_(n_fft), hop_(hop), samples_(samples) {
  if (!std::has_single_bit(n_fft) || samples <= n_fft / 2)
    throw std::invalid_argument("spectrogram geometry");
  window_.resize(n_fft);
  for (std::size_t n = 0; n < n_fft; ++n) // torch.hann_window(periodic=True)
    window_[n] =
        0.5 - 0.5 * std::cos(2.0 * std::numbers::pi * static_cast<double>(n) /
                             static_cast<double>(n_fft));
  padded_.resize(samples + n_fft);
  re_.resize(n_fft);
  im_.resize(n_fft);
}

namespace {

// In-place iterative radix-2 FFT.
void fft(std::vector<double> &re, std::vector<double> &im) {
  const std::size_t n = re.size();
  for (std::size_t i = 1, j = 0; i < n; ++i) {
    std::size_t bit = n >> 1;
    for (; j & bit; bit >>= 1)
      j ^= bit;
    j ^= bit;
    if (i < j) {
      std::swap(re[i], re[j]);
      std::swap(im[i], im[j]);
    }
  }
  for (std::size_t len = 2; len <= n; len <<= 1) {
    const double ang = -2.0 * std::numbers::pi / static_cast<double>(len);
    const double wr = std::cos(ang), wi = std::sin(ang);
    for (std::size_t i = 0; i < n; i += len) {
      double cr = 1.0, ci = 0.0;
      for (std::size_t j = 0; j < len / 2; ++j) {
        const std::size_t u = i + j, v = i + j + len / 2;
        const double tr = re[v] * cr - im[v] * ci;
        const double ti = re[v] * ci + im[v] * cr;
        re[v] = re[u] - tr;
        im[v] = im[u] - ti;
        re[u] += tr;
        im[u] += ti;
        const double ncr = cr * wr - ci * wi;
        ci = cr * wi + ci * wr;
        cr = ncr;
      }
    }
  }
}

} // namespace

void Spectrogram::power(std::span<const double> x, std::span<double> out) {
  if (x.size() != samples_ || out.size() != bins() * frames())
    throw std::invalid_argument("spectrogram: wrong size");
  const std::size_t half = n_fft_ / 2, F = frames(), B = bins();
  // Reflect padding (torch pad_mode="reflect"): the edge sample is not
  // repeated.
  for (std::size_t i = 0; i < half; ++i) {
    padded_[i] = x[half - i];
    padded_[half + samples_ + i] = x[samples_ - 2 - i];
  }
  std::copy(x.begin(), x.end(),
            padded_.begin() + static_cast<std::ptrdiff_t>(half));
  for (std::size_t f = 0; f < F; ++f) {
    for (std::size_t n = 0; n < n_fft_; ++n) {
      re_[n] = padded_[f * hop_ + n] * window_[n];
      im_[n] = 0.0;
    }
    fft(re_, im_);
    for (std::size_t k = 0; k < B; ++k)
      out[k * F + f] = re_[k] * re_[k] + im_[k] * im_[k];
  }
}

void Spectrogram::db(std::span<double> p, double top_db) {
  double mx = -1e300;
  for (double &v : p) {
    v = 10.0 * std::log10(std::max(v, 1e-10));
    mx = std::max(mx, v);
  }
  const double floor = mx - top_db;
  for (double &v : p)
    v = std::max(v, floor);
}

std::vector<float> median_profile(std::span<const double> frames_by_bin,
                                  std::size_t bins) {
  const std::size_t n = frames_by_bin.size() / bins;
  std::vector<float> out(bins);
  std::vector<double> col(n);
  for (std::size_t k = 0; k < bins; ++k) {
    for (std::size_t f = 0; f < n; ++f)
      col[f] = frames_by_bin[f * bins + k];
    const std::size_t mid = (n - 1) / 2;
    std::nth_element(col.begin(),
                     col.begin() + static_cast<std::ptrdiff_t>(mid), col.end());
    out[k] = static_cast<float>(col[mid]);
  }
  return out;
}

} // namespace ayzek::dsp
