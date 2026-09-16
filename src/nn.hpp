#pragma once

// Inference-only layers, each matching its PyTorch counterpart in eval mode.
//
// Two memory layouts, chosen per layer to match how PyTorch indexes it, so no
// weight ever needs reordering:
//   planar       (channels, length)  -- Conv1d, BatchNorm1d
//   interleaved  (steps, features)   -- Linear, LayerNorm, LSTM, attention
//
// Layers own their scratch buffers, sized on first use and reused after, so a
// model built for one input length does no allocation per call.

#include "weights.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace ayzek::nn {

void gelu(float* x, std::size_t n) noexcept;     // exact erf form, PyTorch's default
void relu(float* x, std::size_t n) noexcept;
void sigmoid(float* x, std::size_t n) noexcept;

// (C, L) -> (L, C) and back.
void planar_to_interleaved(const float* x, std::size_t C, std::size_t L, float* y) noexcept;
void interleaved_to_planar(const float* x, std::size_t L, std::size_t C, float* y) noexcept;

struct Conv1d {
    std::size_t cin = 0, cout = 0, k = 0, stride = 1, pad = 0;
    std::vector<float> w;                         // (cout, cin * k)
    std::vector<float> b;                         // (cout) or empty
    mutable std::vector<float> col;               // im2col scratch, (cin * k)

    Conv1d() = default;
    Conv1d(const Weights& W, const std::string& prefix, std::size_t stride, std::size_t pad, bool bias);
    [[nodiscard]] std::size_t out_len(std::size_t L) const noexcept { return (L + 2 * pad - k) / stride + 1; }
    void forward(const float* x, std::size_t L, float* y) const;   // (cin, L) -> (cout, out_len)
};

// (cin, H, W) -> (cout, H', W'), square kernel, zero padding, no dilation.
struct Conv2d {
    std::size_t cin = 0, cout = 0, k = 0, stride = 1, pad = 0;
    std::vector<float> w;                         // (cout, cin * k * k)
    std::vector<float> b;
    mutable std::vector<float> col;

    Conv2d() = default;
    Conv2d(const Weights& W, const std::string& prefix, std::size_t stride, std::size_t pad, bool bias);
    [[nodiscard]] std::size_t out_len(std::size_t L) const noexcept { return (L + 2 * pad - k) / stride + 1; }
    void forward(const float* x, std::size_t H, std::size_t W, float* y) const;
};

// BatchNorm1d and BatchNorm2d are the same per-channel affine in eval mode; for
// 2D pass L = H * W.
struct BatchNorm1d {
    std::vector<float> scale, shift;              // folded: w / sqrt(var + eps), b - mean * scale

    BatchNorm1d() = default;
    BatchNorm1d(const Weights& W, const std::string& prefix, float eps = 1e-5f);
    void forward(float* x, std::size_t L) const noexcept;          // (C, L) in place
};

struct Linear {
    std::size_t in = 0, out = 0;
    std::vector<float> w, b;

    Linear() = default;
    Linear(const Weights& W, const std::string& prefix);
    void forward(const float* x, std::size_t T, float* y) const noexcept;   // (T, in) -> (T, out)
};

struct LayerNorm {
    std::size_t d = 0;
    std::vector<float> g, b;
    float eps = 1e-5f;

    LayerNorm() = default;
    LayerNorm(const Weights& W, const std::string& prefix);
    void forward(float* x, std::size_t T) const noexcept;          // (T, d) in place
};

// One direction of a single-layer PyTorch LSTM. Gate order i, f, g, o.
struct LstmDirection {
    std::size_t in = 0, hidden = 0;
    std::vector<float> w_ih, w_hh, bias;          // (4H, in), (4H, H), b_ih + b_hh
    mutable std::vector<float> gates, h, c;

    LstmDirection() = default;
    LstmDirection(const Weights& W, const std::string& prefix, const std::string& suffix);
    // Writes the hidden state of every step into out, (T, H), in time order.
    void forward(const float* x, std::size_t T, bool reverse, float* out) const;
};

struct BiLstm {
    LstmDirection fwd, bwd;

    BiLstm() = default;
    BiLstm(const Weights& W, const std::string& prefix);
    [[nodiscard]] std::size_t out_dim() const noexcept { return 2 * fwd.hidden; }
    void forward(const float* x, std::size_t T, float* y) const;   // (T, in) -> (T, 2H), [fwd, bwd]
};

// nn.MultiheadAttention(batch_first=True) used as self-attention, eval mode.
struct SelfAttention {
    std::size_t d = 0, heads = 0;
    std::vector<float> in_w, in_b, out_w, out_b;
    mutable std::vector<float> q, k, v, scores, mixed;

    SelfAttention() = default;
    SelfAttention(const Weights& W, const std::string& prefix, std::size_t heads);
    void forward(const float* x, std::size_t T, float* y) const;   // (T, d) -> (T, d)
};

}  // namespace ayzek::nn
