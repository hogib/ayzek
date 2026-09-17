#include "models.hpp"

#include "simd.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ayzek {

// --- detector ----------------------------------------------------------------

Detector::Detector(const Weights &w)
    : c1_(w, "b1.conv.0.0", 2, 3, false), c2_(w, "b1.conv.1.0", 2, 2, false),
      c3_(w, "b1.conv.3.0", 2, 2, false), n1_(w, "b1.conv.0.1"),
      n2_(w, "b1.conv.1.1"), n3_(w, "b1.conv.3.1"), lstm_(w, "b1.lstm"),
      attn_(w, "b1.attn", 4), norm_(w, "b1.norm"), head_norm_(w, "head.0"),
      p1_(w, "p1"), head1_(w, "head.2"), head2_(w, "head.5") {
  w1_ = w.at("w1").f32()[0];
  if (c1_.cin != kChannels || head2_.out != 1)
    throw std::runtime_error(w.path() + ": not a detector");
  planar_.resize(kChannels * kWindow);
  b1_.resize(c1_.cout * c1_.out_len(kWindow));
  b2_.resize(c2_.cout * c2_.out_len(c1_.out_len(kWindow)));
  b3_.resize(c3_.cout * c3_.out_len(c2_.out_len(c1_.out_len(kWindow))));
}

float Detector::logit(std::span<const float> x) {
  if (x.size() != kWindow * kChannels)
    throw std::invalid_argument("detector wants 600 x 3");
  nn::interleaved_to_planar(x.data(), kWindow, kChannels, planar_.data());

  std::size_t L = kWindow;
  c1_.forward(planar_.data(), L, b1_.data());
  L = c1_.out_len(L);
  n1_.forward(b1_.data(), L);
  nn::gelu(b1_.data(), c1_.cout * L);
  c2_.forward(b1_.data(), L, b2_.data());
  L = c2_.out_len(L);
  n2_.forward(b2_.data(), L);
  nn::gelu(b2_.data(), c2_.cout * L);
  c3_.forward(b2_.data(), L, b3_.data());
  L = c3_.out_len(L);
  n3_.forward(b3_.data(), L);
  nn::gelu(b3_.data(), c3_.cout * L);

  const std::size_t T = L, D = c3_.cout;
  conv_out.resize(T * D);
  nn::planar_to_interleaved(b3_.data(), D, T, conv_out.data());
  lstm_out.resize(T * lstm_.out_dim());
  lstm_.forward(conv_out.data(), T, lstm_out.data());
  attn_out.resize(lstm_out.size());
  attn_.forward(lstm_out.data(), T, attn_out.data());

  // LayerNorm(lstm + attention), then the mean over time steps.
  const std::size_t E = lstm_.out_dim();
  thread_local std::vector<float> resid;
  resid.assign(lstm_out.begin(), lstm_out.end());
  simd::add(attn_out.data(), resid.data(), resid.size());
  norm_.forward(resid.data(), T);
  pooled.fill(0.0f);
  for (std::size_t t = 0; t < T; ++t)
    simd::add(resid.data() + t * E, pooled.data(), E);
  for (auto &v : pooled)
    v /= static_cast<float>(T);

  std::array<float, 96> f{}, g{};
  p1_.forward(pooled.data(), 1, f.data());
  for (auto &v : f)
    v *= w1_;
  head_norm_.forward(f.data(), 1);
  head1_.forward(f.data(), 1, g.data());
  nn::gelu(g.data(), g.size());
  float out = 0.0f;
  head2_.forward(g.data(), 1, &out);
  return out;
}

DetectorEnsemble::DetectorEnsemble(const std::vector<Weights> &seeds) {
  for (const auto &w : seeds)
    members_.emplace_back(w);
  if (members_.empty())
    throw std::runtime_error("detector ensemble has no members");
  x_.resize(Detector::kWindow * Detector::kChannels);
}

float DetectorEnsemble::probability(std::span<const float> standardized) {
  for (std::size_t i = 0; i < x_.size(); ++i)
    x_[i] = std::asinh(standardized[i]);
  double acc = 0.0;
  for (auto &m : members_)
    acc += 1.0 / (1.0 + std::exp(-static_cast<double>(m.logit(x_))));
  return static_cast<float>(acc / static_cast<double>(members_.size()));
}

// --- picker ------------------------------------------------------------------

Picker::Picker(const Weights &w)
    : c1_(w, "wave.stem.0", 2, 4, true), c2_(w, "wave.stem.3", 2, 4, true),
      c3_(w, "wave.stem.6", 2, 3, true), n1_(w, "wave.stem.1"),
      n2_(w, "wave.stem.4"), n3_(w, "wave.stem.7"), lstm_(w, "wave.lstm"),
      attn_(w, "wave.attn", 4), norm_(w, "wave.norm"), head1_(w, "head.0"),
      head2_(w, "head.2") {
  b1_.resize(c1_.cout * c1_.out_len(kWindow));
  b2_.resize(c2_.cout * c2_.out_len(c1_.out_len(kWindow)));
  b3_.resize(c3_.cout * c3_.out_len(c2_.out_len(c1_.out_len(kWindow))));
  logits_.resize(kChunks * kClasses);
}

std::span<const float> Picker::logits(std::span<const float> x) {
  if (x.size() != 3 * kWindow)
    throw std::invalid_argument("picker wants 3 x 6000");
  std::size_t L = kWindow;
  c1_.forward(x.data(), L, b1_.data());
  L = c1_.out_len(L);
  n1_.forward(b1_.data(), L);
  nn::relu(b1_.data(), c1_.cout * L);
  c2_.forward(b1_.data(), L, b2_.data());
  L = c2_.out_len(L);
  n2_.forward(b2_.data(), L);
  nn::relu(b2_.data(), c2_.cout * L);
  c3_.forward(b2_.data(), L, b3_.data());
  L = c3_.out_len(L);
  n3_.forward(b3_.data(), L);
  nn::relu(b3_.data(), c3_.cout * L);
  stem_out.assign(b3_.begin(), b3_.end());

  // adaptive_avg_pool1d(750 -> 250). The length divides exactly, so this is a
  // mean over non-overlapping groups of 3.
  const std::size_t C = c3_.cout;
  if (L % kChunks != 0)
    throw std::runtime_error("picker stem length does not divide into chunks");
  const std::size_t f = L / kChunks;
  pooled_out.assign(kChunks * C, 0.0f);
  for (std::size_t c = 0; c < C; ++c)
    for (std::size_t t = 0; t < kChunks; ++t) {
      float s = 0.0f;
      for (std::size_t j = 0; j < f; ++j)
        s += b3_[c * L + t * f + j];
      pooled_out[t * C + c] = s / static_cast<float>(f);
    }

  lstm_out.resize(kChunks * lstm_.out_dim());
  lstm_.forward(pooled_out.data(), kChunks, lstm_out.data());
  wave_out.resize(lstm_out.size());
  attn_.forward(lstm_out.data(), kChunks, wave_out.data());
  simd::add(lstm_out.data(), wave_out.data(), wave_out.size());
  norm_.forward(wave_out.data(), kChunks);

  hidden_.resize(kChunks * head1_.out);
  head1_.forward(wave_out.data(), kChunks, hidden_.data());
  nn::relu(hidden_.data(), hidden_.size());
  head2_.forward(hidden_.data(), kChunks, logits_.data());
  return logits_;
}

Picker::Picks Picker::pick(std::span<const float> x) {
  auto lg = logits(x);
  const double chunk_seconds = static_cast<double>(kWindow) / kChunks / 100.0;
  Picks out{};
  double best_p = -1.0, best_s = -1.0;
  for (std::size_t t = 0; t < kChunks; ++t) {
    const float *r = lg.data() + t * kClasses;
    const float mx = std::max({r[0], r[1], r[2]});
    const double z =
        std::exp(r[0] - mx) + std::exp(r[1] - mx) + std::exp(r[2] - mx);
    const double pp = std::exp(r[1] - mx) / z, ps = std::exp(r[2] - mx) / z;
    if (pp > best_p) {
      best_p = pp;
      out.p_seconds = (static_cast<double>(t) + 0.5) * chunk_seconds;
    }
    if (ps > best_s) {
      best_s = ps;
      out.s_seconds = (static_cast<double>(t) + 0.5) * chunk_seconds;
    }
  }
  out.p_prob = best_p;
  out.s_prob = best_s;
  return out;
}

} // namespace ayzek
