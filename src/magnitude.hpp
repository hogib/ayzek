#pragma once

// Magnitude estimation from one station's 10 s window starting 2 s before P.
// Model: cnn_earthquake's waveform-only regressor trained on FDSN (KO) windows;
// the three station- and event-disjoint partition models are averaged
// (reported MAE 0.42 +- 0.02 on held-out stations).
//
// Inputs are normalised by the station's noise level:
//   waveform     (cleaned - noise mean) / noise sigma, per component
//   spectrogram  dB minus the median noise dB per frequency bin
// NoiseBaseline estimates both online from windows the detector scored as
// noise. With less than 60 s of noise, both inputs use per-window normalisation
// instead, which is what the training data generator used for stations without
// a noise baseline.

#include "dsp.hpp"
#include "nn.hpp"
#include "spectrogram.hpp"
#include "weights.hpp"

#include <array>
#include <deque>
#include <span>
#include <vector>

namespace ayzek {

inline constexpr std::size_t kMagWindow = 1000; // 10 s
inline constexpr std::size_t kMagBins = 65;     // n_fft 128
inline constexpr std::size_t kMagFrames = 32;   // hop 32

class MagnitudeRegressor {
public:
  explicit MagnitudeRegressor(const Weights &w);
  // seq: (1000, 3) interleaved; img: (3, 65, 32) planar.
  [[nodiscard]] float predict(std::span<const float> seq,
                              std::span<const float> img);

  std::vector<float> lstm_out; // kept for tests
  std::array<float, 128> branch1d{}, branch2d{};

private:
  nn::BiLstm lstm_;
  nn::SelfAttention attn_;
  nn::LayerNorm norm_, head_norm_;
  nn::Linear p1_, p2_, head1_, head2_;
  nn::Conv2d c1_, c2_, c3_;
  nn::BatchNorm1d n1_, n2_, n3_;
  float w1_ = 1, w2_ = 1;
  std::vector<float> attn_out_, b1_, b2_, b3_;
};

struct StationNoise {
  bool ready = false;
  std::array<double, 3> mu{}, sigma{};
  std::array<std::vector<float>, 3> profile; // dB per frequency bin
  std::size_t windows = 0;
};

class NoiseBaseline {
public:
  explicit NoiseBaseline(const dsp::Bandpass &bp, std::size_t keep_windows = 30,
                         std::size_t min_windows = 6);
  // Adds one noise window: (1000, 3) interleaved raw counts.
  void add(std::span<const double> raw);
  [[nodiscard]] const StationNoise &noise() const noexcept { return noise_; }

private:
  struct Entry {
    std::array<double, 3> sum{}, sumsq{};
    std::array<std::vector<double>, 3>
        frames; // (frames, bins) dB, per component
  };
  dsp::Conditioner cond_;
  dsp::Spectrogram spec_;
  std::deque<Entry> entries_;
  std::size_t keep_, min_;
  StationNoise noise_;
};

class MagnitudeEstimator {
public:
  MagnitudeEstimator(const std::vector<Weights> &partitions,
                     const dsp::Bandpass &bp);
  // raw: (1000, 3) interleaved counts starting 2 s before P.
  [[nodiscard]] float estimate(std::span<const double> raw,
                               const StationNoise &noise);
  [[nodiscard]] std::size_t size() const noexcept { return models_.size(); }

  // The model inputs from the last call, for tests.
  std::vector<float> seq, img, spec_db;
  std::vector<double> cleaned;
  std::vector<float> per_model;

private:
  std::vector<MagnitudeRegressor> models_;
  bool asinh_seq_ = false; // from the weights' `seq_transform`
  dsp::Conditioner cond_;
  dsp::Spectrogram spec_;
  std::vector<float> scratch_;
  std::vector<double> channel_, power_;
};

} // namespace ayzek
