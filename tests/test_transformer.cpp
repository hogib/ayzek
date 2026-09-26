// Compares the streaming transformer with PyTorch and scipy on DEMI data
// (tools/export_transformer.py): the causal band-pass across a gap, the network
// one token at a time with and without a station context, and the whole front
// end from raw counts to token outputs.

#include "fixture.hpp"
#include "transformer.hpp"

#include <optional>
#include <vector>

using namespace ayzek;

namespace {

void filter(const Weights &model, const Weights &fx) {
  const auto sos = model.at("filter.sos").f64(), zi = model.at("filter.zi").f64();
  CausalBandpass bp({sos.begin(), sos.end()}, {zi.begin(), zi.end()}, 100);
  const auto raw = fx.at("filter.raw").f64();
  const auto *miss = reinterpret_cast<const std::uint8_t *>(fx.at("filter.missing").bytes.data());
  std::vector<CausalBandpass::Sample> out;
  for (std::size_t i = 0; i < raw.size(); ++i)
    bp.push(miss[i] ? std::nullopt : std::optional<double>(raw[i]), out);
  const auto ref = fx.at("filter.out").f64();
  const auto *ref_m =
      reinterpret_cast<const std::uint8_t *>(fx.at("filter.out_missing").bytes.data());
  // A run still open at the end is released once it is a second long; the
  // fixture's last run is longer than that, so every position has an output.
  CHECK(out.size() == raw.size());
  double worst = 0, scale = 0;
  std::size_t mismatched = 0;
  for (std::size_t i = 0; i < out.size(); ++i) {
    mismatched += out[i].missing != static_cast<bool>(ref_m[i]);
    worst = std::max(worst, std::abs(out[i].value - ref[i]));
    scale = std::max(scale, std::abs(ref[i]));
  }
  std::println("  band-pass: {} samples, max abs err {:.2e} (signal max {:.2e}), "
               "{} missing flags differ",
               out.size(), worst, scale, mismatched);
  CHECK(mismatched == 0);
  // onset keeps filtered samples as float32, so the reference is rounded to
  // about 6e-8 of the signal; the filter itself runs in double on both sides.
  CHECK(worst <= 2e-7 * scale);
}

void network(const Weights &model, const Weights &fx) {
  OnsetTransformer net(model);
  const auto stride = net.config().stride;
  const auto x = fx.at("net.x").f32();
  for (const bool with_ctx : {true, false}) {
    net.reset();
    net.set_context(with_ctx ? fx.at("net.ctx").f32() : std::span<const float>{});
    const auto p = fx.at(with_ctx ? "net.p_ctx" : "net.p_null").f32();
    const auto dt = fx.at(with_ctx ? "net.dt_ctx" : "net.dt_null").f32();
    double ep = 0, ed = 0;
    for (std::size_t t = 0; t < p.size(); ++t) {
      const auto o = net.step(x.subspan(t * stride * 4, stride * 4));
      ep = std::max(ep, static_cast<double>(std::abs(o.p - p[t])));
      ed = std::max(ed, static_cast<double>(std::abs(o.dt - dt[t])));
    }
    std::println("  network, {:<10} {} tokens: max |p| err {:.2e}, max |dt| err {:.2e} s",
                 with_ctx ? "context" : "no context", p.size(), ep, ed);
    CHECK(ep < 1e-4 && ed < 1e-3);
  }
}

void stream(const Weights &model, const Weights &fx) {
  OnsetStream s(model);
  const auto raw = fx.at("stream.raw").f64();
  const auto *miss = reinterpret_cast<const std::uint8_t *>(fx.at("stream.missing").bytes.data());
  const auto p = fx.at("stream.p").f32(), dt = fx.at("stream.dt").f32();
  std::vector<OnsetStream::Token> out;
  const std::size_t n = raw.size() / 3;
  for (std::size_t i = 0; i < n; ++i) {
    std::array<std::optional<double>, 3> x;
    for (std::size_t c = 0; c < 3; ++c)
      if (!miss[i * 3 + c])
        x[c] = raw[i * 3 + c];
    s.push(1000 + i, x, out);
  }
  CHECK(out.size() == p.size());
  double ep = 0, ed = 0;
  for (std::size_t t = 0; t < out.size(); ++t) {
    CHECK(out[t].end == 1000 + t * 10 + 9);
    ep = std::max(ep, static_cast<double>(std::abs(out[t].p - p[t])));
    ed = std::max(ed, static_cast<double>(std::abs(out[t].dt - dt[t])));
  }
  std::println("  raw counts -> tokens, {} tokens: max |p| err {:.2e}, max |dt| err {:.2e} s",
               out.size(), ep, ed);
  CHECK(ep < 1e-4 && ed < 1e-3);
}

} // namespace

int main(int argc, char **argv) {
  const auto root = fixture::root(argc, argv);
  const auto fx = fixture::load_or_skip(root + "/data/fixtures/transformer.ayzw");
  const auto model = fixture::load_or_skip(root + "/models/transformer.ayzw");
  std::println("transformer");
  filter(model, fx);
  network(model, fx);
  stream(model, fx);
  std::println("transformer matches PyTorch");
}
