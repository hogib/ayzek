// Compares the magnitude chain with seismic_cli and PyTorch: noise baseline,
// spectrogram, input normalisation with and without a baseline, and the model.

#include "fixture.hpp"
#include "magnitude.hpp"

#include <format>
#include <vector>

using namespace ayzek;

int main(int argc, char **argv) {
  const auto root = fixture::root(argc, argv);
  auto fx = fixture::load_or_skip(root + "/data/fixtures/magnitude.ayzw");
  const auto bp = dsp::Bandpass::load(
      fixture::load_or_skip(root + "/models/bandpass.ayzw"));
  std::vector<Weights> parts;
  for (int p = 0; p < 3; ++p)
    parts.push_back(fixture::load_or_skip(
        std::format("{}/models/magnitude_p{}.ayzw", root, p)));

  // --- noise baseline from six quiet DEMI windows
  // --------------------------------------
  NoiseBaseline nb(bp);
  const auto &noise_raw = fx.at("noise_raw");
  const std::size_t per = kMagWindow * 3;
  for (std::size_t w = 0; w < noise_raw.shape[0]; ++w)
    nb.add(noise_raw.f64().subspan(w * per, per));
  const auto &noise = nb.noise();
  CHECK(noise.ready);
  fixture::Diff mu, sigma, prof;
  for (std::size_t c = 0; c < 3; ++c) {
    sigma.rel = std::max(
        sigma.rel, std::abs(noise.sigma[c] - fx.at("noise_sigma").f64()[c]) /
                       fx.at("noise_sigma").f64()[c]);
    mu.abs =
        std::max(mu.abs, std::abs(noise.mu[c] - fx.at("noise_mu").f64()[c]));
    auto d = fixture::compare<float, float>(
        noise.profile[c],
        fx.at("noise_profile").f32().subspan(c * kMagBins, kMagBins));
    prof.abs = std::max(prof.abs, d.abs);
  }
  std::println("noise baseline: sigma rel err {:.1e}, mu abs err {:.1e}, "
               "profile max abs err {:.1e} dB",
               sigma.rel, mu.abs, prof.abs);
  CHECK(sigma.rel < 1e-7 && mu.abs < 1e-6);
  // torchaudio computes the STFT in float32 and this code in double; the
  // resulting differences are about 1e-3 dB.
  CHECK(prof.abs < 1e-2);

  // --- inputs for two real windows, with and without a baseline
  // ------------------------
  MagnitudeEstimator est(parts, bp);
  const auto &raw = fx.at("raw");
  const std::size_t plane = 3 * kMagBins * kMagFrames;
  fixture::Diff spec, img, fallback, seq;
  for (std::size_t w = 0; w < raw.shape[0]; ++w) {
    (void)est.estimate(raw.f64().subspan(w * per, per), StationNoise{});
    auto df = fixture::compare<float, float>(
        est.img, fx.at("img_fallback").f32().subspan(w * plane, plane));
    fallback.abs = std::max(fallback.abs, df.abs);

    const float m = est.estimate(raw.f64().subspan(w * per, per), noise);
    spec.abs =
        std::max(spec.abs, fixture::compare<float, float>(
                               est.spec_db,
                               fx.at("spec_db").f32().subspan(w * plane, plane))
                               .abs);
    img.abs = std::max(
        img.abs, fixture::compare<float, float>(
                     est.img, fx.at("img").f32().subspan(w * plane, plane))
                     .abs);
    seq.abs =
        std::max(seq.abs, fixture::compare<float, float>(
                              est.seq, fx.at("seq").f32().subspan(w * per, per))
                              .abs);
    std::println("  window {}: magnitude {:.3f}  (torch p0-p2: {:.3f} {:.3f} "
                 "{:.3f}; ours {:.3f} {:.3f} {:.3f})",
                 w, m, fx.at("mag_p0").f32()[w], fx.at("mag_p1").f32()[w],
                 fx.at("mag_p2").f32()[w], est.per_model[0], est.per_model[1],
                 est.per_model[2]);
  }
  std::println("spectrogram dB max abs err {:.1e}; normalised {:.1e}; fallback "
               "{:.1e}; waveform {:.1e}",
               spec.abs, img.abs, fallback.abs, seq.abs);
  CHECK(spec.abs < 1e-2 && img.abs < 1e-2 && fallback.abs < 1e-3 &&
        seq.abs < 1e-3);

  // --- the model itself, on the encoder's own tensors
  // ------------------------------------
  std::vector<MagnitudeRegressor> regs;
  for (const auto &w : parts)
    regs.emplace_back(w);
  const std::size_t n_demi = raw.shape[0],
                    n_corpus = fx.at("corpus_label").numel();
  fixture::Diff out, b1, b2, lstm;
  for (std::size_t i = 0; i < n_demi + n_corpus; ++i) {
    const bool demi = i < n_demi;
    const std::size_t j = demi ? i : i - n_demi;
    auto s =
        (demi ? fx.at("seq") : fx.at("corpus_seq")).f32().subspan(j * per, per);
    auto im = (demi ? fx.at("img") : fx.at("corpus_img"))
                  .f32()
                  .subspan(j * plane, plane);
    for (std::size_t p = 0; p < regs.size(); ++p) {
      const float got = regs[p].predict(s, im);
      out.abs = std::max(out.abs,
                         std::abs(static_cast<double>(got) -
                                  fx.at(std::format("mag_p{}", p)).f32()[i]));
    }
    b1.abs = std::max(b1.abs, fixture::compare<float, float>(
                                  regs[0].branch1d,
                                  fx.at("p0_b1").f32().subspan(i * 128, 128))
                                  .abs);
    b2.abs = std::max(b2.abs, fixture::compare<float, float>(
                                  regs[0].branch2d,
                                  fx.at("p0_b2").f32().subspan(i * 128, 128))
                                  .abs);
    if (i == 0) {
      (void)regs[0].predict(s, im);
      lstm = fixture::compare<float, float>(regs[0].lstm_out,
                                            fx.at("p0_lstm").f32());
    }
  }
  std::println("regressor: bilstm {:.1e}, 1D branch {:.1e}, 2D branch {:.1e}, "
               "output {:.1e}",
               lstm.abs, b1.abs, b2.abs, out.abs);
  CHECK(lstm.abs < 1e-4 && b1.abs < 1e-4 && b2.abs < 1e-3);
  CHECK(out.abs < 1e-3);
  std::println("magnitude chain matches seismic_cli and PyTorch");
}
