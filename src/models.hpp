#pragma once

// The detector and phase picker, built from the layers in nn.hpp. Each class
// follows the structure of its PyTorch module, and the constructor reads
// parameters by their state_dict names.
//
// Instances hold scratch buffers and must be used by one thread at a time; each
// station processor owns its own instances.

#include "nn.hpp"
#include "weights.hpp"

#include <array>
#include <span>
#include <vector>

namespace ayzek {

// archive_pipeline WaveformDetector, branch1d = "cnn-lstm":
//   conv(3->24 k7 s2) conv(24->48 k5 s2) conv(48->96 k5 s2), each BN + GELU
//   -> BiLSTM(96, 48) -> self-attention(4 heads) + residual -> LayerNorm
//   -> time mean -> Linear(96, 96) * w1 -> LayerNorm -> Linear -> GELU ->
//   Linear(96, 1)
class Detector {
public:
  static constexpr std::size_t kWindow = 600; // 6 s at 100 Hz
  static constexpr std::size_t kChannels = 3; // Z, N, E

  explicit Detector(const Weights &w);

  // `x`: (600, 3) interleaved, already standardised and asinh-transformed.
  [[nodiscard]] float logit(std::span<const float> x);

  // Intermediate activations of the last call, compared layer by layer in
  // tests.
  std::vector<float> conv_out, lstm_out, attn_out;
  std::array<float, 96> pooled{};

private:
  nn::Conv1d c1_, c2_, c3_;
  nn::BatchNorm1d n1_, n2_, n3_;
  nn::BiLstm lstm_;
  nn::SelfAttention attn_;
  nn::LayerNorm norm_, head_norm_;
  nn::Linear p1_, head1_, head2_;
  float w1_ = 1.0f;
  std::vector<float> planar_, b1_, b2_, b3_;
};

// Mean of the sigmoid outputs of several detectors (one per training seed).
class DetectorEnsemble {
public:
  explicit DetectorEnsemble(const std::vector<Weights> &seeds);
  // `standardized`: (600, 3) interleaved. Applies asinh before the models, as
  // the training pipeline did.
  [[nodiscard]] float probability(std::span<const float> standardized);
  [[nodiscard]] std::size_t size() const noexcept { return members_.size(); }

private:
  std::vector<Detector> members_;
  std::vector<float> x_;
};

// sphase PhasePicker(arm="wave"): class scores {noise, P, S} for each of 250
// chunks (0.24 s) of a 60 s window.
class Picker {
public:
  static constexpr std::size_t kWindow = 6000;
  static constexpr std::size_t kChunks = 250;
  static constexpr std::size_t kClasses = 3;

  explicit Picker(const Weights &w);

  // `x`: (3, 6000) planar, per-trace per-component z-scored.
  // Returns (250, 3) logits, interleaved.
  [[nodiscard]] std::span<const float> logits(std::span<const float> x);

  struct Picks {
    double p_seconds, p_prob, s_seconds, s_prob; // from the window start
  };
  [[nodiscard]] Picks pick(std::span<const float> x);

  std::vector<float> stem_out, pooled_out, lstm_out, wave_out;

private:
  nn::Conv1d c1_, c2_, c3_;
  nn::BatchNorm1d n1_, n2_, n3_;
  nn::BiLstm lstm_;
  nn::SelfAttention attn_;
  nn::LayerNorm norm_;
  nn::Linear head1_, head2_;
  std::vector<float> b1_, b2_, b3_, logits_, hidden_;
};

} // namespace ayzek
