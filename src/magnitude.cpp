#include "magnitude.hpp"

#include "simd.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ayzek {

// --- regressor ---------------------------------------------------------------------

MagnitudeRegressor::MagnitudeRegressor(const Weights& w)
    : lstm_(w, "b1.lstm"), attn_(w, "b1.attn", 4), norm_(w, "b1.norm"), head_norm_(w, "head.0"),
      p1_(w, "p1"), p2_(w, "p2"), head1_(w, "head.2"), head2_(w, "head.5"),
      c1_(w, "b2.net.0", 1, 1, false), c2_(w, "b2.net.3", 2, 1, false), c3_(w, "b2.net.7", 2, 1, false),
      n1_(w, "b2.net.1"), n2_(w, "b2.net.4"), n3_(w, "b2.net.8") {
    w1_ = w.at("w1").f32()[0];
    w2_ = w.at("w2").f32()[0];
    if (lstm_.out_dim() != 128 || c3_.cout != 128 || head2_.out != 1) throw std::runtime_error(w.path() + ": not a magnitude regressor");
}

float MagnitudeRegressor::predict(std::span<const float> seq, std::span<const float> img) {
    if (seq.size() != kMagWindow * 3 || img.size() != 3 * kMagBins * kMagFrames)
        throw std::invalid_argument("magnitude regressor wants (1000, 3) and (3, 65, 32)");
    const std::size_t T = kMagWindow, E = 128;

    // 1D branch: BiLSTM over raw samples, self-attention, residual + norm, mean.
    lstm_out.resize(T * E);
    lstm_.forward(seq.data(), T, lstm_out.data());
    attn_out_.resize(T * E);
    attn_.forward(lstm_out.data(), T, attn_out_.data());
    simd::add(lstm_out.data(), attn_out_.data(), attn_out_.size());
    norm_.forward(attn_out_.data(), T);
    branch1d.fill(0);
    for (std::size_t t = 0; t < T; ++t) simd::add(attn_out_.data() + t * E, branch1d.data(), E);
    for (auto& v : branch1d) v /= static_cast<float>(T);

    // 2D branch: three conv-BN-GELU stages over the spectrogram, global mean.
    std::size_t H = kMagBins, W = kMagFrames;
    b1_.resize(c1_.cout * c1_.out_len(H) * c1_.out_len(W));
    c1_.forward(img.data(), H, W, b1_.data());
    H = c1_.out_len(H);
    W = c1_.out_len(W);
    n1_.forward(b1_.data(), H * W);
    nn::gelu(b1_.data(), b1_.size());
    b2_.resize(c2_.cout * c2_.out_len(H) * c2_.out_len(W));
    c2_.forward(b1_.data(), H, W, b2_.data());
    H = c2_.out_len(H);
    W = c2_.out_len(W);
    n2_.forward(b2_.data(), H * W);
    nn::gelu(b2_.data(), b2_.size());
    b3_.resize(c3_.cout * c3_.out_len(H) * c3_.out_len(W));
    c3_.forward(b2_.data(), H, W, b3_.data());
    H = c3_.out_len(H);
    W = c3_.out_len(W);
    n3_.forward(b3_.data(), H * W);
    nn::gelu(b3_.data(), b3_.size());
    for (std::size_t c = 0; c < 128; ++c) {
        double s = 0;
        for (std::size_t i = 0; i < H * W; ++i) s += b3_[c * H * W + i];
        branch2d[c] = static_cast<float>(s / static_cast<double>(H * W));
    }

    // a * F1 + b * F2, then the head.
    std::array<float, 128> f1{}, f2{}, g{};
    p1_.forward(branch1d.data(), 1, f1.data());
    p2_.forward(branch2d.data(), 1, f2.data());
    for (std::size_t i = 0; i < 128; ++i) f1[i] = w1_ * f1[i] + w2_ * f2[i];
    head_norm_.forward(f1.data(), 1);
    head1_.forward(f1.data(), 1, g.data());
    nn::gelu(g.data(), g.size());
    float out = 0;
    head2_.forward(g.data(), 1, &out);
    return out;
}

// --- noise baseline ------------------------------------------------------------------

NoiseBaseline::NoiseBaseline(const dsp::Bandpass& bp, std::size_t keep_windows, std::size_t min_windows)
    : cond_(bp, kMagWindow), spec_(128, 32, kMagWindow), keep_(keep_windows), min_(min_windows) {}

void NoiseBaseline::add(std::span<const double> raw) {
    thread_local std::vector<float> unused(kMagWindow * 3);
    cond_.condition(raw, 3, unused);                  // for `cleaned`, as clean_and_filter_1d on each trace

    Entry e;
    std::vector<double> x(kMagWindow), p(kMagBins * kMagFrames);
    for (std::size_t c = 0; c < 3; ++c) {
        for (std::size_t t = 0; t < kMagWindow; ++t) x[t] = cond_.cleaned[t * 3 + c];
        e.sum[c] = simd::sum(x.data(), x.size());
        e.sumsq[c] = simd::dot(x.data(), x.data(), x.size());
        spec_.power(x, p);
        dsp::Spectrogram::db(p);                      // one trace at a time, as the profile builder does
        e.frames[c].resize(kMagFrames * kMagBins);
        for (std::size_t f = 0; f < kMagFrames; ++f)
            for (std::size_t k = 0; k < kMagBins; ++k) e.frames[c][f * kMagBins + k] = p[k * kMagFrames + f];
    }
    entries_.push_back(std::move(e));
    if (entries_.size() > keep_) entries_.pop_front();

    const double n = static_cast<double>(entries_.size() * kMagWindow);
    for (std::size_t c = 0; c < 3; ++c) {
        double s = 0, ss = 0;
        std::vector<double> all;
        all.reserve(entries_.size() * kMagFrames * kMagBins);
        for (const auto& en : entries_) {
            s += en.sum[c];
            ss += en.sumsq[c];
            all.insert(all.end(), en.frames[c].begin(), en.frames[c].end());
        }
        noise_.mu[c] = s / n;
        noise_.sigma[c] = std::sqrt(std::max(ss / n - noise_.mu[c] * noise_.mu[c], 0.0));
        noise_.profile[c] = dsp::median_profile(all, kMagBins);
    }
    noise_.windows = entries_.size();
    noise_.ready = entries_.size() >= min_ && std::ranges::all_of(noise_.sigma, [](double s) { return s > 1e-12; });
}

// --- estimator ------------------------------------------------------------------------

MagnitudeEstimator::MagnitudeEstimator(const std::vector<Weights>& partitions, const dsp::Bandpass& bp)
    : cond_(bp, kMagWindow), spec_(128, 32, kMagWindow) {
    for (const auto& w : partitions) models_.emplace_back(w);
    if (models_.empty()) throw std::runtime_error("no magnitude models");
    seq.resize(kMagWindow * 3);
    img.resize(3 * kMagBins * kMagFrames);
    spec_db.resize(img.size());
    scratch_.resize(kMagWindow * 3);
    channel_.resize(kMagWindow);
    power_.resize(3 * kMagBins * kMagFrames);
}

float MagnitudeEstimator::estimate(std::span<const double> raw, const StationNoise& noise) {
    cond_.condition(raw, 3, scratch_);
    cleaned = cond_.cleaned;

    // Waveform: the station's noise sigma, or the window's own moments without one.
    for (std::size_t c = 0; c < 3; ++c) {
        for (std::size_t t = 0; t < kMagWindow; ++t) channel_[t] = cleaned[t * 3 + c];
        double mu = 0, sd = 0;
        if (noise.ready) {
            mu = noise.mu[c];
            sd = noise.sigma[c];
        } else {
            mu = simd::sum(channel_.data(), kMagWindow) / kMagWindow;
            sd = std::sqrt(std::max(simd::dot(channel_.data(), channel_.data(), kMagWindow) / kMagWindow - mu * mu, 0.0));
        }
        sd = std::max(sd, 1e-12);
        for (std::size_t t = 0; t < kMagWindow; ++t) seq[t * 3 + c] = static_cast<float>((channel_[t] - mu) / sd);
    }

    // Spectrogram of all three components together, so one top_db floor spans them.
    const std::size_t plane = kMagBins * kMagFrames;
    for (std::size_t c = 0; c < 3; ++c) {
        for (std::size_t t = 0; t < kMagWindow; ++t) channel_[t] = cleaned[t * 3 + c];
        spec_.power(channel_, std::span<double>(power_.data() + c * plane, plane));
    }
    dsp::Spectrogram::db(power_);
    for (std::size_t i = 0; i < power_.size(); ++i) spec_db[i] = static_cast<float>(power_[i]);

    if (noise.ready) {
        for (std::size_t c = 0; c < 3; ++c)
            for (std::size_t k = 0; k < kMagBins; ++k)
                for (std::size_t f = 0; f < kMagFrames; ++f) {
                    const std::size_t i = c * plane + k * kMagFrames + f;
                    img[i] = static_cast<float>(power_[i] - noise.profile[c][k]);
                }
    } else {
        // torch's (s - s.mean()) / (s.std() + 1e-6), std with Bessel's correction.
        const double n = static_cast<double>(power_.size());
        const double mean = simd::sum(power_.data(), power_.size()) / n;
        double var = 0;
        for (double v : power_) var += (v - mean) * (v - mean);
        const double sd = std::sqrt(var / (n - 1)) + 1e-6;
        for (std::size_t i = 0; i < power_.size(); ++i) img[i] = static_cast<float>((power_[i] - mean) / sd);
    }

    per_model.clear();
    double acc = 0;
    for (auto& m : models_) {
        per_model.push_back(m.predict(seq, img));
        acc += per_model.back();
    }
    return static_cast<float>(acc / static_cast<double>(models_.size()));
}

}  // namespace ayzek
