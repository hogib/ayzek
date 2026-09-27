#pragma once

// The streaming onset transformer (the `onset` project, OnsetDetector) and the
// front end that feeds it: a causal band-pass that restarts after gaps,
// scaling by the station's noise, and a station context refreshed from quiet
// data. Selected with `--detector transformer`.
//
// Where the 6 s detector scores a fresh window every 0.5 s, this runs one step
// per 0.1 s token and keeps state between steps: the stem's last samples, and
// the keys and values of the last 80 tokens in every layer. Each output
// depends only on samples at or before the end of its token. The outputs are
//
//   p   probability that an event is under way
//   dt  seconds since its onset, so a trigger at time t puts P at t - dt
//
// and, for a model trained with the geometry head (onset ModelConfig.geometry),
// where the event is as seen from this station: epicentral distance and
// back-azimuth, each with the model's own uncertainty (Geometry below). The
// network stage locates events from these (pipeline/locate.hpp).
//
// Weights, geometry and the operating threshold come from
// models/transformer.ayzw (tools/export_transformer.py). docs/impl/14-transformer.md.

#include "nn.hpp"
#include "weights.hpp"

#include <array>
#include <cstdint>
#include <deque>
#include <optional>
#include <span>
#include <vector>

namespace ayzek {

struct TransformerConfig {
  std::size_t d = 0, heads = 0, layers = 0, window = 0, ctx_tokens = 0,
              ctx_layers = 0;
  bool geometry = false; // the export has the geometry head
  // Clamps on the head's log-variance of log distance and log concentration
  // of back-azimuth; exports without `config.geo` used these.
  std::array<float, 4> geo_clamp{-8.0f, 6.0f, -4.0f, 8.0f};
  std::vector<std::array<std::size_t, 3>> stem; // (kernel, stride, channels)
  std::size_t stride = 0;  // samples per token
  std::size_t history = 0; // samples the stem needs for the newest token
  double max_dt = 0, ctx_seconds = 0, fallback_scale_seconds = 0;
  double threshold = 0, release_ratio = 0; // validation operating point
};

// Multi-head attention with separate query and key/value projections,
// `onset.model.Attention`.
struct TransformerAttention {
  nn::Linear q, kv, out;
  std::size_t d = 0, heads = 0;
  TransformerAttention() = default;
  TransformerAttention(const Weights &w, const std::string &prefix,
                       std::size_t heads);
  // (n, d) -> (n, 2d) keys then values.
  void project_kv(const float *mem, std::size_t n, std::vector<float> &k,
                  std::vector<float> &v) const;
  // Queries (m, d) against keys/values (n, d); `bias` (heads, m, n) or null,
  // -inf entries excluded. Writes (m, d).
  void attend(const float *x, std::size_t m, const float *k, const float *v,
              std::size_t n, const float *bias, float *y) const;
};

// `onset.model.OnsetDetector.geometry` at one token: the distance is Gaussian
// in log km, the back-azimuth von Mises. Only meaningful after P.
struct Geometry {
  float log_dist;    // log km
  float log_dist_sd; // standard deviation of log_dist
  float baz;         // station -> event, radians clockwise from north
  float kappa;       // von Mises concentration of baz
};

class OnsetTransformer {
public:
  struct Out {
    float p, dt;
    std::optional<Geometry> geo; // with the geometry head only
  };

  explicit OnsetTransformer(const Weights &w);
  [[nodiscard]] const TransformerConfig &config() const noexcept {
    return cfg_;
  }

  // Forgets every sample and cached token, keeping the context.
  void reset();
  // `ctx`: (n, 4) conditioned station context, or empty for the learned null
  // context.
  void set_context(std::span<const float> ctx);
  // One token: `stride` new conditioned samples, (stride, 4) interleaved.
  Out step(std::span<const float> samples);

private:
  struct Block {
    nn::LayerNorm ln_sa, ln_ca, ln_ff;
    TransformerAttention sa, ca;
    nn::Linear ff1, ff2;
    std::vector<float> ck, cv;       // cross-attention keys/values of the context
    std::deque<std::vector<float>> k, v; // self-attention cache, newest last
  };
  struct ContextLayer {
    nn::LayerNorm ln_sa, ln_ff;
    TransformerAttention sa;
    nn::Linear ff1, ff2;
  };
  struct StemStage {
    nn::Conv1d conv; // pad 0; the causal padding is added by hand
    nn::LayerNorm norm;
    std::size_t k, s;
  };

  // (T, 4) samples -> (T / stride, d) tokens, zero-padded on the left per
  // stage exactly as `onset.model.CausalConv`.
  void stem(std::span<const float> x, std::size_t T, std::vector<float> &out);
  void feed_forward(const nn::Linear &a, const nn::Linear &b, float *x,
                    std::size_t n);
  void encode_context(std::span<const float> ctx);

  TransformerConfig cfg_;
  std::vector<StemStage> stem_;
  nn::Linear stem_proj_;
  std::vector<Block> blocks_;
  // context encoder
  std::vector<float> queries_, null_context_, context_;
  nn::LayerNorm ctx_ln_mem_, ctx_ln_q_, ctx_ln_ff_, ctx_ln_out_;
  TransformerAttention ctx_cross_;
  nn::Linear ctx_ff1_, ctx_ff2_;
  std::vector<ContextLayer> ctx_layers_;
  nn::LayerNorm ln_out_;
  nn::Linear head_, geo_head_;
  std::vector<float> slopes_;
  std::vector<float> samples_; // last `history` samples, (n, 4)
  std::vector<float> e_, h_, a_, hidden_, tokens_, planar_, next_, qbuf_, bias_;
};

// `onset.dsp.causal_filter` as a stream, for one component: each run of
// samples is filtered on its own with second-order sections, starting from
// steady-state conditions after removing the mean of its first second. Runs
// shorter than a second count as missing. The first second of a run is
// therefore released one second late, all at once; after that, one output
// per input.
class CausalBandpass {
public:
  struct Sample {
    double value;
    bool missing;
  };
  // `sos`: (sections, 6); `zi`: (sections, 2), scipy's sosfilt_zi.
  CausalBandpass(std::vector<double> sos, std::vector<double> zi,
                 std::size_t first_run);
  // Appends zero or more outputs, in position order.
  void push(std::optional<double> x, std::vector<Sample> &out);

private:
  double filter(double x);
  std::vector<double> sos_, zi_, state_;
  std::size_t first_run_;
  std::vector<double> run_; // the start of the current run, until it is long enough
  bool running_ = false;
  double mean_ = 0;
};

// One station's transformer detector: raw counts in, token outputs out.
class OnsetStream {
public:
  struct Token {
    std::uint64_t end; // position of the token's last sample
    float p, dt;
    bool refreshed; // the station context was refreshed after this token
    std::optional<Geometry> geo;
  };
  explicit OnsetStream(const Weights &w);
  [[nodiscard]] const TransformerConfig &config() const noexcept {
    return model_.config();
  }
  // Feeds the next position; a missing component is nullopt. Appends the
  // tokens this completes.
  void push(std::uint64_t pos, const std::array<std::optional<double>, 3> &x,
            std::vector<Token> &out);
  [[nodiscard]] std::size_t refreshes() const noexcept { return refreshes_; }

  static constexpr float kQuiet = 0.3f;        // context only from blocks below this
  static constexpr double kBlockSeconds = 60.0; // how often the context may refresh

private:
  void emit(std::uint64_t pos, std::array<double, 3> x, bool missing,
            std::vector<Token> &out);
  void condition(const std::array<double, 3> &x, bool missing);
  void maybe_refresh();

  OnsetTransformer model_;
  std::array<CausalBandpass, 3> filters_;
  std::array<std::deque<CausalBandpass::Sample>, 3> filtered_;
  std::vector<CausalBandpass::Sample> scratch_;
  std::uint64_t next_ = ~std::uint64_t{0}; // position of the next combined sample
  // Before the first scale: the stream's first samples.
  std::vector<std::array<double, 3>> startup_;
  std::vector<bool> startup_missing_;
  bool scaled_ = false;
  std::array<double, 3> scale_{1, 1, 1};
  // The last ctx_seconds of filtered samples, for refreshing the context.
  std::vector<std::array<double, 3>> ring_;
  std::vector<std::uint8_t> ring_missing_;
  std::size_t ring_head_ = 0, ring_count_ = 0;
  std::deque<float> recent_p_;
  std::vector<float> pending_; // conditioned samples of the token being filled
  std::uint64_t tokens_ = 0;
  std::size_t refreshes_ = 0;
  bool refreshed_ = false;
};

} // namespace ayzek
