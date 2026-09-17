#include "nn.hpp"

#include "simd.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ayzek::nn {

void gelu(float *x, std::size_t n) noexcept {
  // Scalar: there is no NEON erf.
  for (std::size_t i = 0; i < n; ++i)
    x[i] = 0.5f * x[i] * (1.0f + std::erf(x[i] * 0.70710678118654752f));
}

void relu(float *x, std::size_t n) noexcept { simd::relu(x, n); }

void sigmoid(float *x, std::size_t n) noexcept {
  for (std::size_t i = 0; i < n; ++i)
    x[i] = 1.0f / (1.0f + std::exp(-x[i]));
}

void planar_to_interleaved(const float *x, std::size_t C, std::size_t L,
                           float *y) noexcept {
  for (std::size_t c = 0; c < C; ++c)
    for (std::size_t t = 0; t < L; ++t)
      y[t * C + c] = x[c * L + t];
}

void interleaved_to_planar(const float *x, std::size_t L, std::size_t C,
                           float *y) noexcept {
  for (std::size_t c = 0; c < C; ++c)
    for (std::size_t t = 0; t < L; ++t)
      y[c * L + t] = x[t * C + c];
}

// --- Conv1d ------------------------------------------------------------------

Conv1d::Conv1d(const Weights &W, const std::string &prefix, std::size_t stride_,
               std::size_t pad_, bool bias)
    : stride(stride_), pad(pad_) {
  const auto &t = W.at(prefix + ".weight");
  if (t.shape.size() != 3)
    throw std::runtime_error(prefix + ".weight is not (cout, cin, k)");
  cout = t.shape[0];
  cin = t.shape[1];
  k = t.shape[2];
  w.assign(t.f32().begin(),
           t.f32().end()); // (cout, cin, k) row-major == (cout, cin * k)
  if (bias)
    b = W.vec(prefix + ".bias", {cout});
  col.resize(cin * k);
}

void Conv1d::forward(const float *x, std::size_t L, float *y) const {
  const std::size_t Lo = out_len(L), ck = cin * k;
  for (std::size_t o = 0; o < Lo; ++o) {
    // Copy the receptive field into `col`, with zeros outside [0, L).
    const std::ptrdiff_t start = static_cast<std::ptrdiff_t>(o * stride) -
                                 static_cast<std::ptrdiff_t>(pad);
    const bool interior =
        start >= 0 && start + static_cast<std::ptrdiff_t>(k) <=
                          static_cast<std::ptrdiff_t>(L);
    for (std::size_t c = 0; c < cin; ++c) {
      const float *src = x + c * L;
      float *dst = col.data() + c * k;
      if (interior) {
        std::copy_n(src + start, k, dst);
      } else {
        for (std::size_t j = 0; j < k; ++j) {
          const std::ptrdiff_t p = start + static_cast<std::ptrdiff_t>(j);
          dst[j] =
              (p >= 0 && p < static_cast<std::ptrdiff_t>(L)) ? src[p] : 0.0f;
        }
      }
    }
    for (std::size_t oc = 0; oc < cout; ++oc)
      y[oc * Lo + o] = simd::dot(w.data() + oc * ck, col.data(), ck) +
                       (b.empty() ? 0.0f : b[oc]);
  }
}

// --- Conv2d
// --------------------------------------------------------------------

Conv2d::Conv2d(const Weights &W, const std::string &prefix, std::size_t stride_,
               std::size_t pad_, bool bias)
    : stride(stride_), pad(pad_) {
  const auto &t = W.at(prefix + ".weight");
  if (t.shape.size() != 4 || t.shape[2] != t.shape[3])
    throw std::runtime_error(prefix + ".weight is not (cout, cin, k, k)");
  cout = t.shape[0];
  cin = t.shape[1];
  k = t.shape[2];
  w.assign(t.f32().begin(), t.f32().end());
  if (bias)
    b = W.vec(prefix + ".bias", {cout});
  col.resize(cin * k * k);
}

void Conv2d::forward(const float *x, std::size_t H, std::size_t Wd,
                     float *y) const {
  const std::size_t Ho = out_len(H), Wo = out_len(Wd), patch = cin * k * k;
  const auto sH = static_cast<std::ptrdiff_t>(H),
             sW = static_cast<std::ptrdiff_t>(Wd);
  for (std::size_t oy = 0; oy < Ho; ++oy) {
    for (std::size_t ox = 0; ox < Wo; ++ox) {
      const std::ptrdiff_t y0 = static_cast<std::ptrdiff_t>(oy * stride) -
                                static_cast<std::ptrdiff_t>(pad);
      const std::ptrdiff_t x0 = static_cast<std::ptrdiff_t>(ox * stride) -
                                static_cast<std::ptrdiff_t>(pad);
      float *dst = col.data();
      for (std::size_t c = 0; c < cin; ++c) {
        const float *plane = x + c * H * Wd;
        for (std::size_t i = 0; i < k; ++i) {
          const std::ptrdiff_t yy = y0 + static_cast<std::ptrdiff_t>(i);
          for (std::size_t j = 0; j < k; ++j) {
            const std::ptrdiff_t xx = x0 + static_cast<std::ptrdiff_t>(j);
            *dst++ = (yy >= 0 && yy < sH && xx >= 0 && xx < sW)
                         ? plane[yy * sW + xx]
                         : 0.0f;
          }
        }
      }
      for (std::size_t oc = 0; oc < cout; ++oc)
        y[oc * Ho * Wo + oy * Wo + ox] =
            simd::dot(w.data() + oc * patch, col.data(), patch) +
            (b.empty() ? 0.0f : b[oc]);
    }
  }
}

// --- BatchNorm1d
// ----------------------------------------------------------------

BatchNorm1d::BatchNorm1d(const Weights &W, const std::string &prefix,
                         float eps) {
  auto g = W.at(prefix + ".weight").f32();
  auto bb = W.at(prefix + ".bias").f32();
  auto mean = W.at(prefix + ".running_mean").f32();
  auto var = W.at(prefix + ".running_var").f32();
  scale.resize(g.size());
  shift.resize(g.size());
  for (std::size_t c = 0; c < g.size(); ++c) {
    // Computed in double to limit rounding in the folded constants.
    const double s = g[c] / std::sqrt(static_cast<double>(var[c]) + eps);
    scale[c] = static_cast<float>(s);
    shift[c] = static_cast<float>(bb[c] - mean[c] * s);
  }
}

void BatchNorm1d::forward(float *x, std::size_t L) const noexcept {
  for (std::size_t c = 0; c < scale.size(); ++c)
    simd::affine(x + c * L, scale[c], shift[c], L);
}

// --- Linear / LayerNorm
// -------------------------------------------------------------

Linear::Linear(const Weights &W, const std::string &prefix) {
  const auto &t = W.at(prefix + ".weight");
  out = t.shape.at(0);
  in = t.shape.at(1);
  w.assign(t.f32().begin(), t.f32().end());
  b = W.vec(prefix + ".bias", {out});
}

void Linear::forward(const float *x, std::size_t T, float *y) const noexcept {
  for (std::size_t t = 0; t < T; ++t)
    simd::gemv(w.data(), out, in, x + t * in, b.data(), y + t * out);
}

LayerNorm::LayerNorm(const Weights &W, const std::string &prefix) {
  g = W.vec(prefix + ".weight");
  b = W.vec(prefix + ".bias", {g.size()});
  d = g.size();
}

void LayerNorm::forward(float *x, std::size_t T) const noexcept {
  for (std::size_t t = 0; t < T; ++t) {
    float *r = x + t * d;
    double mean = 0.0, var = 0.0;
    for (std::size_t i = 0; i < d; ++i)
      mean += r[i];
    mean /= static_cast<double>(d);
    for (std::size_t i = 0; i < d; ++i)
      var += (r[i] - mean) * (r[i] - mean);
    var /= static_cast<double>(d);
    const float inv = static_cast<float>(1.0 / std::sqrt(var + eps));
    for (std::size_t i = 0; i < d; ++i)
      r[i] = (r[i] - static_cast<float>(mean)) * inv * g[i] + b[i];
  }
}

// --- LSTM
// ---------------------------------------------------------------------------

LstmDirection::LstmDirection(const Weights &W, const std::string &prefix,
                             const std::string &suffix) {
  const auto &ih = W.at(prefix + ".weight_ih_l0" + suffix);
  hidden = ih.shape.at(0) / 4;
  in = ih.shape.at(1);
  w_ih.assign(ih.f32().begin(), ih.f32().end());
  w_hh = W.vec(prefix + ".weight_hh_l0" + suffix, {4 * hidden, hidden});
  auto bi = W.at(prefix + ".bias_ih_l0" + suffix, {4 * hidden}).f32();
  auto bh = W.at(prefix + ".bias_hh_l0" + suffix, {4 * hidden}).f32();
  bias.resize(4 * hidden);
  for (std::size_t i = 0; i < bias.size(); ++i)
    bias[i] = bi[i] + bh[i];
  gates.resize(4 * hidden);
  h.resize(hidden);
  c.resize(hidden);
}

void LstmDirection::forward(const float *x, std::size_t T, bool reverse,
                            float *out) const {
  const std::size_t H = hidden;
  std::fill(h.begin(), h.end(), 0.0f);
  std::fill(c.begin(), c.end(), 0.0f);
  for (std::size_t s = 0; s < T; ++s) {
    const std::size_t t = reverse ? T - 1 - s : s;
    simd::gemv(w_ih.data(), 4 * H, in, x + t * in, bias.data(), gates.data());
    for (std::size_t r = 0; r < 4 * H; ++r)
      gates[r] += simd::dot(w_hh.data() + r * H, h.data(), H);
    for (std::size_t j = 0; j < H; ++j) {
      const float i_g = 1.0f / (1.0f + std::exp(-gates[j]));
      const float f_g = 1.0f / (1.0f + std::exp(-gates[H + j]));
      const float g_g = std::tanh(gates[2 * H + j]);
      const float o_g = 1.0f / (1.0f + std::exp(-gates[3 * H + j]));
      c[j] = f_g * c[j] + i_g * g_g;
      h[j] = o_g * std::tanh(c[j]);
    }
    std::copy_n(h.data(), H, out + t * H);
  }
}

BiLstm::BiLstm(const Weights &W, const std::string &prefix)
    : fwd(W, prefix, ""), bwd(W, prefix, "_reverse") {}

void BiLstm::forward(const float *x, std::size_t T, float *y) const {
  const std::size_t H = fwd.hidden;
  thread_local std::vector<float> a, b;
  a.resize(T * H);
  b.resize(T * H);
  fwd.forward(x, T, false, a.data());
  bwd.forward(x, T, true, b.data());
  for (std::size_t t = 0; t < T; ++t) {
    std::copy_n(a.data() + t * H, H, y + t * 2 * H);
    std::copy_n(b.data() + t * H, H, y + t * 2 * H + H);
  }
}

// --- attention
// ----------------------------------------------------------------------

SelfAttention::SelfAttention(const Weights &W, const std::string &prefix,
                             std::size_t heads_)
    : heads(heads_) {
  const auto &t = W.at(prefix + ".in_proj_weight");
  d = t.shape.at(1);
  if (t.shape.at(0) != 3 * d || d % heads != 0)
    throw std::runtime_error(prefix + ": bad attention shape");
  in_w.assign(t.f32().begin(), t.f32().end());
  in_b = W.vec(prefix + ".in_proj_bias", {3 * d});
  out_w = W.vec(prefix + ".out_proj.weight", {d, d});
  out_b = W.vec(prefix + ".out_proj.bias", {d});
}

void SelfAttention::forward(const float *x, std::size_t T, float *y) const {
  const std::size_t hd = d / heads;
  q.resize(T * d);
  k.resize(T * d);
  v.resize(T * d);
  scores.resize(T);
  mixed.resize(T * d);
  for (std::size_t t = 0; t < T; ++t) {
    simd::gemv(in_w.data(), d, d, x + t * d, in_b.data(), q.data() + t * d);
    simd::gemv(in_w.data() + d * d, d, d, x + t * d, in_b.data() + d,
               k.data() + t * d);
    simd::gemv(in_w.data() + 2 * d * d, d, d, x + t * d, in_b.data() + 2 * d,
               v.data() + t * d);
  }
  const float inv_sqrt = 1.0f / std::sqrt(static_cast<float>(hd));
  std::fill(mixed.begin(), mixed.end(), 0.0f);
  for (std::size_t h = 0; h < heads; ++h) {
    for (std::size_t i = 0; i < T; ++i) {
      const float *qi = q.data() + i * d + h * hd;
      float mx = -1e30f;
      for (std::size_t j = 0; j < T; ++j) {
        scores[j] = simd::dot(qi, k.data() + j * d + h * hd, hd) * inv_sqrt;
        mx = std::max(mx, scores[j]);
      }
      double z = 0.0; // softmax, shifted by the max for stability
      for (std::size_t j = 0; j < T; ++j)
        z += std::exp(static_cast<double>(scores[j] - mx));
      float *out = mixed.data() + i * d + h * hd;
      for (std::size_t j = 0; j < T; ++j) {
        const float p = static_cast<float>(
            std::exp(static_cast<double>(scores[j] - mx)) / z);
        simd::axpy(p, v.data() + j * d + h * hd, out, hd);
      }
    }
  }
  for (std::size_t t = 0; t < T; ++t)
    simd::gemv(out_w.data(), d, d, mixed.data() + t * d, out_b.data(),
               y + t * d);
}

} // namespace ayzek::nn
