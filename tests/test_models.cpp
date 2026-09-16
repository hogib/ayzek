// Detector and picker against PyTorch on real windows, layer by layer.

#include "fixture.hpp"
#include "models.hpp"
#include "simd.hpp"

#include <format>
#include <vector>

using namespace ayzek;

namespace {

void report(const char* what, fixture::Diff d) { std::println("  {:<28} max abs err {:.2e}", what, d.abs); }

void detector(const std::string& root) {
    auto fx = fixture::load_or_skip(root + "/data/fixtures/detector.ayzw");
    std::vector<Weights> seeds;
    for (int s : {42, 43, 44}) seeds.push_back(fixture::load_or_skip(std::format("{}/models/detector_s{}.ayzw", root, s)));

    const auto& input = fx.at("input");             // (k, 600, 3), asinh already applied
    const std::size_t k = input.shape[0], per = input.shape[1] * input.shape[2];

    std::println("detector ({} windows, {} backend)", k, simd::kBackend);
    Detector first(seeds[0]);
    fixture::Diff conv, lstm, attn, pooled, logit;
    for (std::size_t w = 0; w < k; ++w) {
        const float got = first.logit(input.f32().subspan(w * per, per));
        const auto T = fx.at("conv").shape[1], D = fx.at("conv").shape[2];
        auto cmp = [&](fixture::Diff& acc, const std::vector<float>& ours, const char* name, std::size_t width) {
            auto d = fixture::compare<float, float>(ours, fx.at(name).f32().subspan(w * T * width, T * width));
            acc.abs = std::max(acc.abs, d.abs);
        };
        cmp(conv, first.conv_out, "conv", D);
        cmp(lstm, first.lstm_out, "lstm", 96);
        cmp(attn, first.attn_out, "attn", 96);
        auto dp = fixture::compare<float, float>(first.pooled, fx.at("pooled").f32().subspan(w * 96, 96));
        pooled.abs = std::max(pooled.abs, dp.abs);
        logit.abs = std::max(logit.abs, std::abs(static_cast<double>(got) - fx.at("logit_s42").f32()[w]));
    }
    report("conv stack", conv);
    report("bilstm", lstm);
    report("attention", attn);
    report("norm + mean pool", pooled);
    report("logit (seed 42)", logit);
    CHECK(conv.abs < 1e-4 && lstm.abs < 1e-4 && attn.abs < 1e-4 && pooled.abs < 1e-4);
    CHECK(logit.abs < 1e-3);

    // The ensemble works from standardised input and applies asinh itself.
    DetectorEnsemble ens(seeds);
    std::vector<float> standardized(per);
    double worst = 0;
    for (std::size_t w = 0; w < k; ++w) {
        auto x = input.f32().subspan(w * per, per);
        for (std::size_t i = 0; i < per; ++i) standardized[i] = std::sinh(x[i]);
        const float p = ens.probability(standardized);
        const float ref = fx.at("prob").f32()[w];
        worst = std::max(worst, static_cast<double>(std::abs(p - ref)));
        std::println("  window {}: ensemble p = {:.5f}  (torch {:.5f})", w, p, ref);
    }
    CHECK(worst < 1e-4);
}

void picker(const std::string& root) {
    auto fx = fixture::load_or_skip(root + "/data/fixtures/picker.ayzw");
    Picker pk(fixture::load_or_skip(root + "/models/spicker.ayzw"));
    const auto& input = fx.at("input");             // (k, 3, 6000)
    const std::size_t k = input.shape[0], per = input.shape[1] * input.shape[2];

    std::println("picker ({} windows)", k);
    fixture::Diff stem, pooled, lstm, wave, logits;
    for (std::size_t w = 0; w < k; ++w) {
        auto lg = pk.logits(input.f32().subspan(w * per, per));
        auto upd = [&](fixture::Diff& acc, std::span<const float> ours, const char* name) {
            const auto& ref = fx.at(name);
            const std::size_t each = ref.numel() / k;
            auto d = fixture::compare<float, float>(ours, ref.f32().subspan(w * each, each));
            acc.abs = std::max(acc.abs, d.abs);
        };
        upd(stem, pk.stem_out, "stem");
        upd(pooled, pk.pooled_out, "pooled");
        upd(lstm, pk.lstm_out, "lstm");
        upd(wave, pk.wave_out, "wave");
        upd(logits, lg, "logits");

        auto p = pk.pick(input.f32().subspan(w * per, per));
        std::println("  window {}: P {:.2f} s ({:.2f})  S {:.2f} s ({:.2f})", w, p.p_seconds, p.p_prob, p.s_seconds, p.s_prob);
    }
    report("stem", stem);
    report("pooled", pooled);
    report("bilstm", lstm);
    report("attention + norm", wave);
    report("logits", logits);
    CHECK(stem.abs < 1e-3 && pooled.abs < 1e-3 && lstm.abs < 1e-4 && wave.abs < 1e-4);
    CHECK(logits.abs < 1e-3);
}

}  // namespace

int main(int argc, char** argv) {
    const auto root = fixture::root(argc, argv);
    detector(root);
    picker(root);
    std::println("models match PyTorch");
}
