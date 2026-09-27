#include "transformer.hpp"

#include "simd.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <stdexcept>

namespace ayzek {

namespace {

constexpr float kNegInf = -std::numeric_limits<float>::infinity();

std::vector<float> alibi_slopes(std::size_t heads) {
  std::vector<float> s(heads);
  for (std::size_t i = 0; i < heads; ++i)
    s[i] = static_cast<float>(std::pow(
        2.0, -8.0 * static_cast<double>(i + 1) / static_cast<double>(heads)));
  return s;
}

float sigmoid(float x) { return 1.0f / (1.0f + std::exp(-x)); }

} // namespace

// --- attention
// ----------------------------------------------------------------------

TransformerAttention::TransformerAttention(const Weights &w,
                                           const std::string &prefix,
                                           std::size_t heads_)
    : q(w, prefix + ".q"), kv(w, prefix + ".kv"), out(w, prefix + ".out"),
      d(q.out), heads(heads_) {
  if (kv.out != 2 * d || d % heads != 0)
    throw std::runtime_error(prefix + ": bad attention shape");
}

void TransformerAttention::project_kv(const float *mem, std::size_t n,
                                      std::vector<float> &k,
                                      std::vector<float> &v) const {
  std::vector<float> both(n * 2 * d);
  kv.forward(mem, n, both.data());
  k.resize(n * d);
  v.resize(n * d);
  for (std::size_t t = 0; t < n; ++t) {
    std::copy_n(both.data() + t * 2 * d, d, k.data() + t * d);
    std::copy_n(both.data() + t * 2 * d + d, d, v.data() + t * d);
  }
}

void TransformerAttention::attend(const float *x, std::size_t m,
                                  const float *k, const float *v, std::size_t n,
                                  const float *bias, float *y) const {
  const std::size_t hd = d / heads;
  thread_local std::vector<float> qs, mixed, scores;
  qs.resize(m * d);
  mixed.assign(m * d, 0.0f);
  scores.resize(n);
  q.forward(x, m, qs.data());
  const float inv_sqrt = 1.0f / std::sqrt(static_cast<float>(hd));
  for (std::size_t h = 0; h < heads; ++h) {
    for (std::size_t i = 0; i < m; ++i) {
      const float *qi = qs.data() + i * d + h * hd;
      const float *bi = bias ? bias + (h * m + i) * n : nullptr;
      float mx = kNegInf;
      for (std::size_t j = 0; j < n; ++j) {
        const float b = bi ? bi[j] : 0.0f;
        scores[j] = b == kNegInf
                        ? kNegInf
                        : simd::dot(qi, k + j * d + h * hd, hd) * inv_sqrt + b;
        mx = std::max(mx, scores[j]);
      }
      double z = 0.0;
      for (std::size_t j = 0; j < n; ++j)
        if (scores[j] != kNegInf)
          z += std::exp(static_cast<double>(scores[j] - mx));
      float *o = mixed.data() + i * d + h * hd;
      for (std::size_t j = 0; j < n; ++j) {
        if (scores[j] == kNegInf)
          continue;
        const float p = static_cast<float>(
            std::exp(static_cast<double>(scores[j] - mx)) / z);
        simd::axpy(p, v + j * d + h * hd, o, hd);
      }
    }
  }
  out.forward(mixed.data(), m, y);
}

// --- the network
// ----------------------------------------------------------------------

OnsetTransformer::OnsetTransformer(const Weights &w) {
  const auto ints = w.at("config.ints").i64();
  const auto floats = w.at("config.floats").f64();
  if (ints.size() != 6 || floats.size() != 6)
    throw std::runtime_error(w.path() + ": not an onset transformer export");
  cfg_.d = static_cast<std::size_t>(ints[0]);
  cfg_.heads = static_cast<std::size_t>(ints[1]);
  cfg_.layers = static_cast<std::size_t>(ints[2]);
  cfg_.window = static_cast<std::size_t>(ints[3]);
  cfg_.ctx_tokens = static_cast<std::size_t>(ints[4]);
  cfg_.ctx_layers = static_cast<std::size_t>(ints[5]);
  cfg_.max_dt = floats[0];
  cfg_.ctx_seconds = floats[1];
  cfg_.fallback_scale_seconds = floats[2];
  cfg_.threshold = floats[3];
  cfg_.release_ratio = floats[4];
  if (floats[5] != 100.0)
    throw std::runtime_error(w.path() + ": trained at a sample rate other than 100 Hz");

  const auto &st = w.at("config.stem");
  const auto sv = st.i64();
  cfg_.stride = 1;
  for (std::size_t i = 0; i < st.shape.at(0); ++i) {
    const auto k = static_cast<std::size_t>(sv[3 * i]),
               s = static_cast<std::size_t>(sv[3 * i + 1]),
               c = static_cast<std::size_t>(sv[3 * i + 2]);
    cfg_.stem.push_back({k, s, c});
    cfg_.stride *= s;
    stem_.push_back({nn::Conv1d(w, std::format("stem.net.{}.conv", 3 * i), s, 0, true),
                     nn::LayerNorm(w, std::format("stem.net.{}.norm", 3 * i + 1)), k,
                     s});
  }
  std::size_t rf = 1;
  for (auto it = cfg_.stem.rbegin(); it != cfg_.stem.rend(); ++it)
    rf = (rf - 1) * (*it)[1] + (*it)[0];
  cfg_.history = (rf + cfg_.stride - 1) / cfg_.stride * cfg_.stride;
  stem_proj_ = nn::Linear(w, "stem.proj");

  for (std::size_t i = 0; i < cfg_.layers; ++i) {
    const auto p = std::format("blocks.{}", i);
    Block b;
    b.ln_sa = nn::LayerNorm(w, p + ".ln_sa");
    b.ln_ca = nn::LayerNorm(w, p + ".ln_ca");
    b.ln_ff = nn::LayerNorm(w, p + ".ln_ff");
    b.sa = TransformerAttention(w, p + ".sa", cfg_.heads);
    b.ca = TransformerAttention(w, p + ".ca", cfg_.heads);
    b.ff1 = nn::Linear(w, p + ".ff.0");
    b.ff2 = nn::Linear(w, p + ".ff.3");
    blocks_.push_back(std::move(b));
  }

  queries_ = w.vec("context.queries", {cfg_.ctx_tokens, cfg_.d});
  null_context_ = w.vec("null_context", {cfg_.ctx_tokens, cfg_.d});
  ctx_ln_mem_ = nn::LayerNorm(w, "context.ln_mem");
  ctx_ln_q_ = nn::LayerNorm(w, "context.ln_q");
  ctx_ln_ff_ = nn::LayerNorm(w, "context.ln_ff");
  ctx_ln_out_ = nn::LayerNorm(w, "context.ln_out");
  ctx_cross_ = TransformerAttention(w, "context.cross", cfg_.heads);
  ctx_ff1_ = nn::Linear(w, "context.ff.0");
  ctx_ff2_ = nn::Linear(w, "context.ff.3");
  for (std::size_t i = 0; i < cfg_.ctx_layers; ++i) {
    const auto p = std::format("context.self_layers.{}", i);
    ctx_layers_.push_back({nn::LayerNorm(w, p + ".ln_sa"), nn::LayerNorm(w, p + ".ln_ff"),
                           TransformerAttention(w, p + ".sa", cfg_.heads),
                           nn::Linear(w, p + ".ff.0"), nn::Linear(w, p + ".ff.3")});
  }
  ln_out_ = nn::LayerNorm(w, "ln_out");
  head_ = nn::Linear(w, "head");
  if (head_.out != 2 || stem_.front().conv.cin != 4)
    throw std::runtime_error(w.path() + ": unexpected input or head width");
  if (w.has("geo_head.weight")) {
    geo_head_ = nn::Linear(w, "geo_head");
    if (geo_head_.out != 5 || geo_head_.in != cfg_.d)
      throw std::runtime_error(w.path() + ": unexpected geometry head shape");
    cfg_.geometry = true;
    if (w.has("config.geo")) {
      const auto c = w.at("config.geo", {4}).f64();
      for (std::size_t i = 0; i < 4; ++i)
        cfg_.geo_clamp[i] = static_cast<float>(c[i]);
    }
  }
  slopes_ = alibi_slopes(cfg_.heads);
  set_context({});
}

void OnsetTransformer::reset() {
  samples_.clear();
  for (auto &b : blocks_) {
    b.k.clear();
    b.v.clear();
  }
}

void OnsetTransformer::stem(std::span<const float> x, std::size_t T,
                            std::vector<float> &out) {
  // planar_ holds the current stage's input, (C, L).
  std::size_t C = 4, L = T;
  planar_.resize(C * L);
  nn::interleaved_to_planar(x.data(), L, C, planar_.data());
  for (auto &st : stem_) {
    const std::size_t pad = st.k - st.s, Lp = L + pad;
    thread_local std::vector<float> padded;
    padded.assign(C * Lp, 0.0f);
    for (std::size_t c = 0; c < C; ++c)
      std::copy_n(planar_.data() + c * L, L, padded.data() + c * Lp + pad);
    const std::size_t Lo = st.conv.out_len(Lp), Co = st.conv.cout;
    next_.resize(Co * Lo);
    st.conv.forward(padded.data(), Lp, next_.data());
    // LayerNorm over channels at each step, then GELU: in (L, C) layout.
    thread_local std::vector<float> inter;
    inter.resize(Co * Lo);
    nn::planar_to_interleaved(next_.data(), Co, Lo, inter.data());
    st.norm.forward(inter.data(), Lo);
    nn::gelu(inter.data(), inter.size());
    planar_.resize(Co * Lo);
    nn::interleaved_to_planar(inter.data(), Lo, Co, planar_.data());
    C = Co;
    L = Lo;
  }
  thread_local std::vector<float> inter;
  inter.resize(C * L);
  nn::planar_to_interleaved(planar_.data(), C, L, inter.data());
  out.resize(L * cfg_.d);
  stem_proj_.forward(inter.data(), L, out.data());
}

void OnsetTransformer::feed_forward(const nn::Linear &a, const nn::Linear &b,
                                    float *x, std::size_t n) {
  // x += b(gelu(a(x))), with x already normalised by the caller into h_.
  hidden_.resize(n * a.out);
  a.forward(x, n, hidden_.data());
  nn::gelu(hidden_.data(), hidden_.size());
  a_.resize(n * b.out);
  b.forward(hidden_.data(), n, a_.data());
}

void OnsetTransformer::encode_context(std::span<const float> ctx) {
  const std::size_t d = cfg_.d, K = cfg_.ctx_tokens, n_s = ctx.size() / 4;
  std::vector<float> mem;
  stem(ctx, n_s, mem);
  const std::size_t n = mem.size() / d;
  ctx_ln_mem_.forward(mem.data(), n);
  std::vector<float> km, vm;
  ctx_cross_.project_kv(mem.data(), n, km, vm);

  // Tokens spanning a gap are masked, unless every one does.
  std::vector<bool> gap(n, false);
  bool all = true;
  for (std::size_t t = 0; t < n; ++t) {
    for (std::size_t j = 0; j < cfg_.stride; ++j)
      gap[t] = gap[t] || ctx[(t * cfg_.stride + j) * 4 + 3] > 0.0f;
    all = all && gap[t];
  }
  bias_.assign(cfg_.heads * K * n, 0.0f);
  if (!all)
    for (std::size_t i = 0; i < cfg_.heads * K; ++i)
      for (std::size_t t = 0; t < n; ++t)
        if (gap[t])
          bias_[i * n + t] = kNegInf;

  std::vector<float> z = queries_, h = queries_, o(K * d);
  ctx_ln_q_.forward(h.data(), K);
  ctx_cross_.attend(h.data(), K, km.data(), vm.data(), n, bias_.data(), o.data());
  simd::add(o.data(), z.data(), z.size());
  h = z;
  ctx_ln_ff_.forward(h.data(), K);
  feed_forward(ctx_ff1_, ctx_ff2_, h.data(), K);
  simd::add(a_.data(), z.data(), z.size());
  for (auto &l : ctx_layers_) {
    h = z;
    l.ln_sa.forward(h.data(), K);
    std::vector<float> kk, vv;
    l.sa.project_kv(h.data(), K, kk, vv);
    l.sa.attend(h.data(), K, kk.data(), vv.data(), K, nullptr, o.data());
    simd::add(o.data(), z.data(), z.size());
    h = z;
    l.ln_ff.forward(h.data(), K);
    feed_forward(l.ff1, l.ff2, h.data(), K);
    simd::add(a_.data(), z.data(), z.size());
  }
  ctx_ln_out_.forward(z.data(), K);
  context_ = std::move(z);
}

void OnsetTransformer::set_context(std::span<const float> ctx) {
  if (ctx.empty())
    context_ = null_context_;
  else
    encode_context(ctx);
  for (auto &b : blocks_)
    b.ca.project_kv(context_.data(), cfg_.ctx_tokens, b.ck, b.cv);
}

OnsetTransformer::Out OnsetTransformer::step(std::span<const float> x) {
  if (x.size() != cfg_.stride * 4)
    throw std::invalid_argument("transformer step wants one token of samples");
  const std::size_t d = cfg_.d;
  samples_.insert(samples_.end(), x.begin(), x.end());
  if (samples_.size() > cfg_.history * 4)
    samples_.erase(samples_.begin(),
                   samples_.end() - static_cast<std::ptrdiff_t>(cfg_.history * 4));
  // Until `history` samples exist the stem sees the stream from its first
  // sample, as the full forward pass pads it; after that only the newest
  // `history` samples can reach the newest token.
  stem(samples_, samples_.size() / 4, tokens_);
  e_.assign(tokens_.end() - static_cast<std::ptrdiff_t>(d), tokens_.end());

  for (auto &b : blocks_) {
    h_ = e_;
    b.ln_sa.forward(h_.data(), 1);
    std::vector<float> k, v;
    b.sa.project_kv(h_.data(), 1, k, v);
    b.k.push_back(std::move(k));
    b.v.push_back(std::move(v));
    if (b.k.size() > cfg_.window) {
      b.k.pop_front();
      b.v.pop_front();
    }
    const std::size_t n = b.k.size();
    thread_local std::vector<float> K, V;
    K.resize(n * d);
    V.resize(n * d);
    for (std::size_t j = 0; j < n; ++j) {
      std::copy_n(b.k[j].data(), d, K.data() + j * d);
      std::copy_n(b.v[j].data(), d, V.data() + j * d);
    }
    bias_.resize(cfg_.heads * n);
    for (std::size_t hh = 0; hh < cfg_.heads; ++hh)
      for (std::size_t j = 0; j < n; ++j)
        bias_[hh * n + j] = -slopes_[hh] * static_cast<float>(n - 1 - j);
    qbuf_.resize(d);
    b.sa.attend(h_.data(), 1, K.data(), V.data(), n, bias_.data(), qbuf_.data());
    simd::add(qbuf_.data(), e_.data(), d);

    h_ = e_;
    b.ln_ca.forward(h_.data(), 1);
    b.ca.attend(h_.data(), 1, b.ck.data(), b.cv.data(), cfg_.ctx_tokens, nullptr,
                qbuf_.data());
    simd::add(qbuf_.data(), e_.data(), d);

    h_ = e_;
    b.ln_ff.forward(h_.data(), 1);
    feed_forward(b.ff1, b.ff2, h_.data(), 1);
    simd::add(a_.data(), e_.data(), d);
  }
  h_ = e_;
  ln_out_.forward(h_.data(), 1);
  float o[2];
  head_.forward(h_.data(), 1, o);
  Out out{sigmoid(o[0]), static_cast<float>(cfg_.max_dt) * sigmoid(o[1]), {}};
  if (cfg_.geometry) {
    // log distance, its log-variance, the back-azimuth as an unnormalised
    // (sin, cos) and the log of its concentration; clamps as in onset.
    float g[5];
    geo_head_.forward(h_.data(), 1, g);
    const auto &c = cfg_.geo_clamp;
    out.geo = Geometry{g[0], std::exp(0.5f * std::clamp(g[1], c[0], c[1])),
                       std::atan2(g[2], g[3]),
                       std::exp(std::clamp(g[4], c[2], c[3]))};
  }
  return out;
}

// --- causal band-pass
// ----------------------------------------------------------------

CausalBandpass::CausalBandpass(std::vector<double> sos, std::vector<double> zi,
                               std::size_t first_run)
    : sos_(std::move(sos)), zi_(std::move(zi)), state_(zi_.size()),
      first_run_(first_run) {
  if (sos_.size() % 6 != 0 || zi_.size() * 3 != sos_.size())
    throw std::runtime_error("band-pass: sos must be (n, 6) and zi (n, 2)");
  run_.reserve(first_run_);
}

double CausalBandpass::filter(double x) {
  // scipy.signal.sosfilt: direct form II transposed, one section after another.
  for (std::size_t s = 0; s < sos_.size() / 6; ++s) {
    const double *c = sos_.data() + 6 * s;
    double *z = state_.data() + 2 * s;
    const double y = c[0] * x + z[0];
    z[0] = c[1] * x - c[4] * y + z[1];
    z[1] = c[2] * x - c[5] * y;
    x = y;
  }
  return x;
}

void CausalBandpass::push(std::optional<double> x, std::vector<Sample> &out) {
  if (!x) {
    // A gap ends the run; a run still too short to filter becomes missing.
    for (std::size_t i = 0; i < run_.size(); ++i)
      out.push_back({0.0, true});
    run_.clear();
    running_ = false;
    out.push_back({0.0, true});
    return;
  }
  if (running_) {
    out.push_back({filter(*x - mean_), false});
    return;
  }
  run_.push_back(*x);
  if (run_.size() < first_run_)
    return;
  double m = 0;
  for (double v : run_)
    m += v;
  mean_ = m / static_cast<double>(run_.size());
  for (std::size_t i = 0; i < state_.size(); ++i)
    state_[i] = zi_[i] * (run_[0] - mean_);
  for (double v : run_)
    out.push_back({filter(v - mean_), false});
  run_.clear();
  running_ = true;
}

// --- the station stream
// ----------------------------------------------------------------

namespace {
CausalBandpass make_filter(const Weights &w) {
  const auto sos = w.at("filter.sos").f64(), zi = w.at("filter.zi").f64();
  return CausalBandpass({sos.begin(), sos.end()}, {zi.begin(), zi.end()}, 100);
}

std::array<double, 3> rms(const std::vector<std::array<double, 3>> &x,
                          const std::vector<std::uint8_t> &missing) {
  std::array<double, 3> s{0, 0, 0};
  std::size_t n = 0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (missing[i])
      continue;
    for (std::size_t c = 0; c < 3; ++c)
      s[c] += x[i][c] * x[i][c];
    ++n;
  }
  for (auto &v : s)
    v = n ? std::sqrt(v / static_cast<double>(n)) : 0.0;
  for (auto &v : s)
    if (!(v > 0))
      v = 1.0;
  return s;
}
} // namespace

OnsetStream::OnsetStream(const Weights &w)
    : model_(w), filters_{make_filter(w), make_filter(w), make_filter(w)} {
  const auto &c = model_.config();
  const std::size_t ctx_len =
      static_cast<std::size_t>(std::lround(c.ctx_seconds * 100.0 / static_cast<double>(c.stride))) *
      c.stride;
  ring_.resize(ctx_len);
  ring_missing_.resize(ctx_len);
  pending_.reserve(c.stride * 4);
}

void OnsetStream::push(std::uint64_t pos,
                       const std::array<std::optional<double>, 3> &x,
                       std::vector<Token> &out) {
  if (next_ == ~std::uint64_t{0})
    next_ = pos;
  for (std::size_t c = 0; c < 3; ++c) {
    scratch_.clear();
    filters_[c].push(x[c], scratch_);
    filtered_[c].insert(filtered_[c].end(), scratch_.begin(), scratch_.end());
  }
  while (!filtered_[0].empty() && !filtered_[1].empty() && !filtered_[2].empty()) {
    std::array<double, 3> v{};
    bool missing = false;
    for (std::size_t c = 0; c < 3; ++c) {
      v[c] = filtered_[c].front().value;
      missing = missing || filtered_[c].front().missing;
      filtered_[c].pop_front();
    }
    if (missing)
      v = {0, 0, 0};
    emit(next_++, v, missing, out);
  }
}

void OnsetStream::emit(std::uint64_t pos, std::array<double, 3> x, bool missing,
                       std::vector<Token> &out) {
  ring_[ring_head_] = x;
  ring_missing_[ring_head_] = missing;
  ring_head_ = (ring_head_ + 1) % ring_.size();
  ring_count_ = std::min(ring_count_ + 1, ring_.size());

  const auto &c = model_.config();
  auto one = [&](std::uint64_t p, const std::array<double, 3> &xx, bool m) {
    condition(xx, m);
    if (pending_.size() < c.stride * 4)
      return;
    const auto o = model_.step(pending_);
    pending_.clear();
    ++tokens_;
    recent_p_.push_back(o.p);
    const std::size_t ctx_tok = ring_.size() / c.stride;
    if (recent_p_.size() > ctx_tok)
      recent_p_.pop_front();
    refreshed_ = false;
    const auto block = static_cast<std::uint64_t>(kBlockSeconds * 100.0) / c.stride;
    if (tokens_ % block == 0)
      maybe_refresh();
    out.push_back({p, o.p, o.dt, refreshed_, o.geo});
  };

  if (scaled_) {
    one(pos, x, missing);
    return;
  }
  // The first scale is the RMS of the stream's first samples, as for a station
  // with no baseline in training.
  startup_.push_back(x);
  startup_missing_.push_back(missing);
  const auto n0 = static_cast<std::size_t>(c.fallback_scale_seconds * 100.0);
  if (startup_.size() < n0)
    return;
  std::vector<std::uint8_t> m(startup_missing_.begin(), startup_missing_.end());
  scale_ = rms(startup_, m);
  scaled_ = true;
  const std::uint64_t first = pos + 1 - startup_.size();
  for (std::size_t i = 0; i < startup_.size(); ++i)
    one(first + i, startup_[i], startup_missing_[i]);
  startup_.clear();
  startup_missing_.clear();
}

void OnsetStream::condition(const std::array<double, 3> &x, bool missing) {
  for (std::size_t c = 0; c < 3; ++c)
    pending_.push_back(missing ? 0.0f
                               : static_cast<float>(std::asinh(x[c] / scale_[c])));
  pending_.push_back(missing ? 1.0f : 0.0f);
}

void OnsetStream::maybe_refresh() {
  // After each minute: if the last minute is complete and the detector stayed
  // quiet through it, that minute becomes the context and its RMS the scale.
  if (ring_count_ < ring_.size())
    return;
  for (float p : recent_p_)
    if (p >= kQuiet)
      return;
  std::size_t miss = 0;
  for (auto m : ring_missing_)
    miss += m;
  if (static_cast<double>(miss) >= 0.05 * static_cast<double>(ring_.size()))
    return;
  // Oldest first.
  std::vector<std::array<double, 3>> x(ring_.size());
  std::vector<std::uint8_t> m(ring_.size());
  for (std::size_t i = 0; i < ring_.size(); ++i) {
    const std::size_t j = (ring_head_ + i) % ring_.size();
    x[i] = ring_[j];
    m[i] = ring_missing_[j];
  }
  scale_ = rms(x, m);
  std::vector<float> ctx;
  ctx.reserve(x.size() * 4);
  for (std::size_t i = 0; i < x.size(); ++i) {
    for (std::size_t c = 0; c < 3; ++c)
      ctx.push_back(m[i] ? 0.0f : static_cast<float>(std::asinh(x[i][c] / scale_[c])));
    ctx.push_back(m[i] ? 1.0f : 0.0f);
  }
  model_.set_context(ctx);
  ++refreshes_;
  refreshed_ = true;
}

} // namespace ayzek
