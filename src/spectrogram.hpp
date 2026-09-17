#pragma once

// Log-power spectrogram, equivalent to the torchaudio transforms used to build
// the magnitude regressor's training data:
//
//   T.Spectrogram(n_fft=128, hop_length=32, power=2.0)
//       periodic Hann window of n_fft, center=True with reflect padding,
//       one-sided, unnormalised: |sum_n x[n] w[n] e^{-i 2 pi k n / N}|^2
//   T.AmplitudeToDB(stype="power", top_db=80)
//       10 log10(max(p, 1e-10)), then floored at (max over the whole tensor -
//       80)
//
// The top_db floor is computed over all channels passed to `db` together.
// seismic_cli applies it per component for noise profiles and over all three
// components for event windows; callers do the same.

#include <cstddef>
#include <span>
#include <vector>

namespace ayzek::dsp {

class Spectrogram {
public:
  Spectrogram(std::size_t n_fft, std::size_t hop, std::size_t samples);

  [[nodiscard]] std::size_t bins() const noexcept { return n_fft_ / 2 + 1; }
  [[nodiscard]] std::size_t frames() const noexcept {
    return samples_ / hop_ + 1;
  }

  // `x`: one channel of `samples`. Writes (bins, frames) power.
  void power(std::span<const double> x, std::span<double> out);

  // In place over (channels, bins, frames): dB with one top_db floor for all.
  static void db(std::span<double> p, double top_db = 80.0);

private:
  std::size_t n_fft_, hop_, samples_;
  std::vector<double> window_, padded_, re_, im_;
};

// Per-bin median over frames of a (frames x bins) dB matrix. For an even number
// of frames it returns the lower middle value, as torch.median does.
std::vector<float> median_profile(std::span<const double> frames_by_bin,
                                  std::size_t bins);

} // namespace ayzek::dsp
