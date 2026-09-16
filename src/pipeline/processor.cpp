#include "processor.hpp"

#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <thread>

namespace ayzek::pipeline {

namespace {

using Clock = std::chrono::steady_clock;

double ms_since(Clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

std::uint64_t align_up(std::uint64_t pos, std::size_t step) { return (pos + step - 1) / step * step; }

}  // namespace

Processor::Processor(Station& st, const std::vector<Weights>& detector, const Weights& picker,
                     const std::vector<Weights>& magnitude, const dsp::Bandpass& bp, ProcessorConfig cfg, Bus& bus,
                     const std::string& scores_path, const std::string& scores_in_path)
    : st_(st), cfg_(cfg), bus_(bus), detector_(detector), noise_(bp),
      cond6_(bp, Detector::kWindow), cond60_(bp, Picker::kWindow) {
    if (cfg_.pick) picker_ = std::make_unique<Picker>(picker);
    if (cfg_.magnitude && !magnitude.empty()) magnitude_ = std::make_unique<MagnitudeEstimator>(magnitude, bp);
    else cfg_.magnitude = false;
    raw_.resize(Picker::kWindow * 3);
    noise_raw_.resize(kMagWindow * 3);
    standardized_.resize(Detector::kWindow * 3);
    planar60_.resize(Picker::kWindow * 3);
    tmp_.resize(Picker::kWindow);
    stats_.window_ms.reserve(1 << 16);
    if (!scores_path.empty()) {
        scores_.open(scores_path);
        scores_ << "station,window_start,probability\n";
    }
    // Detector probabilities from an earlier run with --scores. They do not
    // depend on trigger settings, so re-using them makes a run with different
    // settings give the same result as full inference, much faster.
    if (!scores_in_path.empty()) {
        std::ifstream in(scores_in_path);
        if (!in) throw std::runtime_error("cannot open " + scores_in_path);
        std::string line;
        std::getline(in, line);
        while (std::getline(in, line)) {
            const auto a = line.find(','), b = line.rfind(',');
            if (a == std::string::npos || b == a) continue;
            cached_[epoch_to_pos(std::stod(line.substr(a + 1, b - a - 1)))] = std::stof(line.substr(b + 1));
        }
    }
}

// Copies positions [start, start + n) of the three components into `out`,
// (n, 3) interleaved. Returns false if the range contains a gap or is no longer
// in the ring; `resume` is then the first position after the problem.
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
        above_ = 0;                          // a gap interrupts a run of windows
        next_ = align_up(std::max(resume, start + 1), cfg_.step);
        return;
    }
    const auto t0 = Clock::now();
    float p = 0;
    if (auto it = cached_.find(start); it != cached_.end()) {
        p = it->second;
    } else {
        cond6_.condition(std::span<const double>(raw_.data(), Detector::kWindow * 3), 3, standardized_);
        p = detector_.probability(standardized_);
    }
    const double ms = ms_since(t0);
    ++stats_.windows;
    stats_.window_ms.push_back(static_cast<float>(ms));

    const double t_start = pos_to_epoch(start);
    recent_.emplace_back(start, p);
    while (recent_.size() > 64) recent_.pop_front();
    // Shortest representation that reads back as the same float.
    if (scores_.is_open()) scores_ << std::format("{},{:.2f},{}\n", st_.code, t_start, p);

    // Before `from`, windows only feed the noise baseline; no triggers.
    const bool warmup = cfg_.from > 0 && t_start < cfg_.from;
    if (warmup) {
        active_ = p >= cfg_.release;
    } else if (!active_) {
        // Trigger when `trigger_windows` consecutive windows reach `threshold`,
        // or at once when a window of the run reaches `instant_threshold`. The
        // detection is dated by the first window of the run and declared at the
        // end of the last. Windows within `retrigger_seconds` of the previous
        // trigger do not start a run.
        if (p >= cfg_.threshold && t_start - last_trigger_ >= cfg_.retrigger_seconds) {
            if (above_++ == 0) run_start_ = start;
        } else {
            above_ = 0;
        }
        const double run_t = pos_to_epoch(run_start_);
        if (above_ > 0 && (above_ >= cfg_.trigger_windows || p >= cfg_.instant_threshold)) {
            active_ = true;
            below_ = 0;
            above_ = 0;
            last_trigger_ = run_t;
            ++stats_.detections;
            bus_.send(Detection{st_.code, run_t, t_start + Detector::kWindow / kFs, p, ms});
            if (cfg_.pick && pending_pick_ == kUnset) {
                pending_pick_ = run_start_ - static_cast<std::uint64_t>(cfg_.picker_lead * kFs);
                pending_trigger_ = run_t;
            }
            // Early magnitude window: P is assumed 3.5 s after the start of the
            // first window above the threshold; the regressor's window starts 2 s
            // before P.
            if (cfg_.magnitude && early_mag_ == kUnset) {
                early_mag_ = run_start_ + 150;
                early_trigger_ = run_t;
            }
        }
    } else if (p < cfg_.release) {
        if (++below_ >= cfg_.release_windows) active_ = false;
    } else {
        below_ = 0;
    }
    if (cfg_.magnitude) maybe_add_noise(start);
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
    const double p_time = t + pk.p_seconds, declared = t + Picker::kWindow / kFs;
    bus_.send(Pick{st_.code, pending_trigger_, p_time, t + pk.s_seconds, pk.p_prob, pk.s_prob, declared, ms});

    // If the P pick is confident and near the trigger, estimate the magnitude
    // from the 10 s window starting 2 s before it (the training alignment). The
    // data is already available at this point.
    if (cfg_.magnitude && pk.p_prob >= 0.5 && p_time >= pending_trigger_ - 3 && p_time <= pending_trigger_ + 9)
        run_magnitude(epoch_to_pos(p_time) - 200, pending_trigger_, true, declared);
}

void Processor::try_early_magnitude(std::uint64_t available) {
    if (early_mag_ == kUnset || early_mag_ + kMagWindow > available) return;
    const std::uint64_t start = early_mag_;
    early_mag_ = kUnset;
    run_magnitude(start, early_trigger_, false, pos_to_epoch(start + kMagWindow));
}

void Processor::run_magnitude(std::uint64_t start, double trigger, bool at_pick, double declared_at) {
    std::uint64_t resume = 0;
    if (!extract(start, kMagWindow, raw_, resume)) return;
    const auto t0 = Clock::now();
    const float m = magnitude_->estimate(std::span<const double>(raw_.data(), kMagWindow * 3), noise_.noise());
    const double ms = ms_since(t0);
    ++stats_.magnitudes;
    stats_.magnitude_ms.push_back(static_cast<float>(ms));
    const auto& nz = noise_.noise();
    bus_.send(MagnitudeEstimate{st_.code, trigger, at_pick, pos_to_epoch(start), m, nz.ready ? nz.windows : 0,
                                declared_at, ms});
}

// At most every `noise_every` seconds, adds the 10 s ending at the newest scored
// window to the noise baseline, provided every detector window overlapping it
// scored below `noise_below`. Signal included in the baseline would inflate the
// noise sigma and lower later magnitude estimates.
void Processor::maybe_add_noise(std::uint64_t window_start) {
    const std::uint64_t end = window_start + Detector::kWindow;
    if (end < kMagWindow + Detector::kWindow || active_) return;
    const std::uint64_t begin = end - kMagWindow;
    if (end < last_noise_end_ + static_cast<std::uint64_t>(cfg_.noise_every * kFs)) return;
    if (recent_.empty() || recent_.front().first + Detector::kWindow > begin) return;   // not enough scored history
    for (const auto& [ws, p] : recent_)
        if (ws + Detector::kWindow > begin && p >= cfg_.noise_below) return;
    std::uint64_t resume = 0;
    if (!extract(begin, kMagWindow, noise_raw_, resume)) return;
    noise_.add(std::span<const double>(noise_raw_.data(), kMagWindow * 3));
    stats_.noise_windows = noise_.noise().windows;
    last_noise_end_ = end;
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
            // Start at the beginning of the data even with `from`, so that the
            // earlier windows build the noise baseline.
            std::uint64_t first = latest_base;
            if (cfg_.from > 0 && !cfg_.magnitude) first = std::max(first, epoch_to_pos(cfg_.from));
            next_ = align_up(first, cfg_.step);
        }

        bool progressed = false;
        while (next_ + Detector::kWindow <= available && !stop.load(std::memory_order_relaxed)) {
            score_window(next_);
            progressed = true;
            // Run picker and magnitude jobs as soon as their data is available,
            // so their messages stay in stream-time order with the detections.
            try_pick(available);
            if (cfg_.magnitude) try_early_magnitude(available);
        }
        try_pick(available);
        if (cfg_.magnitude) try_early_magnitude(available);

        // Raise the ring floor, keeping the picker lead, one noise window before
        // the next window, and the windows of scheduled jobs.
        const std::uint64_t history = std::max<std::uint64_t>(lead, kMagWindow);
        std::uint64_t keep = next_ > history ? next_ - history : 0;
        if (pending_pick_ != kUnset) keep = std::min(keep, pending_pick_);
        if (early_mag_ != kUnset) keep = std::min(keep, early_mag_);
        for (auto& cs : st_.comp) {
            const auto base = cs.base.load(std::memory_order_acquire);
            if (keep > base) cs.ring.set_floor(std::min<std::uint64_t>(keep - base, cs.ring.written()));
        }

        // Lower bound on declared_at of any later message: the end of the next
        // detector window or of a scheduled job's window.
        double until = pos_to_epoch(next_ + Detector::kWindow);
        if (pending_pick_ != kUnset) until = std::min(until, pos_to_epoch(pending_pick_ + Picker::kWindow));
        if (early_mag_ != kUnset) until = std::min(until, pos_to_epoch(early_mag_ + kMagWindow));
        if (until - last_progress_ >= 0.5) {
            bus_.send(Progress{st_.code, until});
            last_progress_ = until;
        }

        if (!progressed) {
            // `done` was read before `available`, so all ingested data has been
            // processed.
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
