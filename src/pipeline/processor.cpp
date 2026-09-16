#include "processor.hpp"

#include <algorithm>
#include <chrono>
#include <thread>

namespace ayzek::pipeline {

namespace {

using Clock = std::chrono::steady_clock;

double ms_since(Clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

std::uint64_t align_up(std::uint64_t pos, std::size_t step) { return (pos + step - 1) / step * step; }

}  // namespace

Processor::Processor(Station& st, const std::vector<Weights>& detector, const Weights& picker, const dsp::Bandpass& bp,
                     ProcessorConfig cfg, Bus& bus, const std::string& scores_path)
    : st_(st), cfg_(cfg), bus_(bus), detector_(detector),
      cond6_(bp, Detector::kWindow), cond60_(bp, Picker::kWindow) {
    if (cfg_.pick) picker_ = std::make_unique<Picker>(picker);
    raw_.resize(Picker::kWindow * 3);
    standardized_.resize(Detector::kWindow * 3);
    planar60_.resize(Picker::kWindow * 3);
    tmp_.resize(Picker::kWindow);
    stats_.window_ms.reserve(1 << 16);
    if (!scores_path.empty()) {
        scores_.open(scores_path);
        scores_ << "station,window_start,probability\n";
    }
}

// Copies [start, start + n) of all three components into `out`, (n, 3)
// interleaved. False if any sample is a gap or already reclaimed; `resume` is
// then the first position a window could start without that gap.
bool Processor::extract(std::uint64_t start, std::size_t n, std::vector<double>& out, std::uint64_t& resume) {
    for (std::size_t c = 0; c < 3; ++c) {
        auto& cs = st_.comp[c];
        const auto base = cs.base.load(std::memory_order_acquire);
        if (base == kUnset || start < base) {
            resume = base == kUnset ? start + 1 : base;
            return false;
        }
        auto got = cs.ring.copy_out(start - base, std::span<std::int32_t>(tmp_.data(), n));
        if (!got) {
            resume = start + 1;
            return false;
        }
        for (std::size_t t = n; t-- > 0;) {
            if (tmp_[t] == kGap) {
                resume = start + t + 1;
                return false;
            }
        }
        for (std::size_t t = 0; t < n; ++t) out[t * 3 + c] = tmp_[t];
    }
    return true;
}

void Processor::score_window(std::uint64_t start) {
    std::uint64_t resume = 0;
    if (!extract(start, Detector::kWindow, raw_, resume)) {
        ++stats_.gap_windows;
        next_ = align_up(std::max(resume, start + 1), cfg_.step);
        return;
    }
    const auto t0 = Clock::now();
    cond6_.condition(std::span<const double>(raw_.data(), Detector::kWindow * 3), 3, standardized_);
    const float p = detector_.probability(standardized_);
    const double ms = ms_since(t0);
    ++stats_.windows;
    stats_.window_ms.push_back(static_cast<float>(ms));

    const double t_start = pos_to_epoch(start);
    if (scores_.is_open()) scores_ << st_.code << ',' << std::format("{:.2f}", t_start) << ',' << p << '\n';

    if (!active_) {
        if (p >= cfg_.threshold && t_start - last_trigger_ >= cfg_.retrigger_seconds) {
            active_ = true;
            below_ = 0;
            last_trigger_ = t_start;
            ++stats_.detections;
            bus_.send(Detection{st_.code, t_start, t_start + Detector::kWindow / kFs, p, ms});
            if (cfg_.pick && pending_pick_ == kUnset) {
                pending_pick_ = start - static_cast<std::uint64_t>(cfg_.picker_lead * kFs);
                pending_trigger_ = t_start;
            }
        }
    } else if (p < cfg_.release) {
        if (++below_ >= cfg_.release_windows) active_ = false;
    } else {
        below_ = 0;
    }
    next_ = start + cfg_.step;
}

void Processor::try_pick(std::uint64_t available) {
    if (pending_pick_ == kUnset || pending_pick_ + Picker::kWindow > available) return;
    const std::uint64_t start = pending_pick_;
    pending_pick_ = kUnset;
    std::uint64_t resume = 0;
    if (!extract(start, Picker::kWindow, raw_, resume)) {
        ++stats_.abandoned_picks;
        Log::get().line("pick", "33", "{:<5} picker window at {} crosses a gap, skipped", st_.code, hms(pos_to_epoch(start)));
        return;
    }
    const auto t0 = Clock::now();
    cond60_.condition_planar(raw_, 3, planar60_);
    const auto pk = picker_->pick(planar60_);
    const double ms = ms_since(t0);
    ++stats_.picks;
    stats_.pick_ms.push_back(static_cast<float>(ms));
    const double t = pos_to_epoch(start);
    bus_.send(Pick{st_.code, pending_trigger_, t + pk.p_seconds, t + pk.s_seconds, pk.p_prob, pk.s_prob, t + Picker::kWindow / kFs, ms});
}

void Processor::run(const std::atomic<bool>& stop) {
    const auto lead = static_cast<std::uint64_t>(cfg_.picker_lead * kFs);
    while (!stop.load(std::memory_order_relaxed)) {
        const bool done = st_.ingest_done.load(std::memory_order_acquire);

        std::uint64_t available = ~std::uint64_t{0}, latest_base = 0;
        bool have_all = true;
        for (auto& cs : st_.comp) {
            const auto base = cs.base.load(std::memory_order_acquire);
            if (base == kUnset) { have_all = false; break; }
            latest_base = std::max(latest_base, base);
            available = std::min(available, cs.end());
        }
        if (!have_all) {
            if (done) break;
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
            continue;
        }
        if (next_ == kUnset) {
            std::uint64_t first = latest_base;
            if (cfg_.from > 0) first = std::max(first, epoch_to_pos(cfg_.from));
            next_ = align_up(first, cfg_.step);
        }

        bool progressed = false;
        while (next_ + Detector::kWindow <= available && !stop.load(std::memory_order_relaxed)) {
            score_window(next_);
            progressed = true;
            // A picker window is due as soon as its data is in, even mid-backlog.
            if (pending_pick_ != kUnset && pending_pick_ + Picker::kWindow <= available) try_pick(available);
        }
        try_pick(available);

        // Reclaim what no window, and no pending or future picker window, can need.
        std::uint64_t keep = next_ > lead ? next_ - lead : 0;
        if (pending_pick_ != kUnset) keep = std::min(keep, pending_pick_);
        for (auto& cs : st_.comp) {
            const auto base = cs.base.load(std::memory_order_acquire);
            if (keep > base) cs.ring.set_floor(std::min<std::uint64_t>(keep - base, cs.ring.written()));
        }

        // Nothing sent from here on is declared before the next window ends or,
        // if a pick is pending, before its window ends.
        double until = pos_to_epoch(next_ + Detector::kWindow);
        if (pending_pick_ != kUnset) until = std::min(until, pos_to_epoch(pending_pick_ + Picker::kWindow));
        if (until - last_progress_ >= 0.5) {
            bus_.send(Progress{st_.code, until});
            last_progress_ = until;
        }

        if (!progressed) {
            // `done` was read before `available`, so everything ingest will ever
            // write was already counted: nothing more can arrive.
            if (done) {
                if (pending_pick_ != kUnset) {
                    ++stats_.abandoned_picks;
                    pending_pick_ = kUnset;
                }
                break;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(2));
        }
    }
    bus_.send(StationDone{st_.code});
}

}  // namespace ayzek::pipeline
