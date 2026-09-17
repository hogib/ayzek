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
    if (cfg_.detector == DetectorKind::StaLta) stalta_ = std::make_unique<dsp::StaLta>(cfg_.stalta);
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
        stalta_next_ = kUnset;               // and restarts the STA/LTA
        next_ = align_up(std::max(resume, start + 1), cfg_.step);
        return;
    }
    double ms = 0;
    const float p = stalta_ ? score_stalta(start, ms) : score_model(start, ms);
    ++stats_.windows;
    stats_.window_ms.push_back(static_cast<float>(ms));
    recent_.emplace_back(start, p);
    while (recent_.size() > 64) recent_.pop_front();
    // Shortest representation that reads back as the same float.
    if (scores_.is_open()) scores_ << std::format("{},{:.2f},{}\n", st_.code, pos_to_epoch(start), p);
    if (cfg_.magnitude) maybe_add_noise(start);
    next_ = start + cfg_.step;
}

// Detector ensemble on the window; returns its probability after applying the
// trigger rule.
float Processor::score_model(std::uint64_t start, double& ms) {
    const auto t0 = Clock::now();
    float p = 0;
    if (auto it = cached_.find(start); it != cached_.end()) {
        p = it->second;
    } else {
        cond6_.condition(std::span<const double>(raw_.data(), Detector::kWindow * 3), 3, standardized_);
        p = detector_.probability(standardized_);
    }
    ms = ms_since(t0);

    const double t_start = pos_to_epoch(start);
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
        if (above_ > 0 && (above_ >= cfg_.trigger_windows || p >= cfg_.instant_threshold)) {
            active_ = true;
            below_ = 0;
            above_ = 0;
            last_trigger_ = pos_to_epoch(run_start_);
            trigger(run_start_, t_start + Detector::kWindow / kFs, p, ms);
        }
    } else if (p < cfg_.release) {
        if (++below_ >= cfg_.release_windows) active_ = false;
    } else {
        below_ = 0;
    }
    return p;
}

// Feeds the samples of the window not yet seen to the STA/LTA, one at a time,
// and triggers at the first sample whose ratio reaches `stalta_on`. Returns the
// largest ratio among those samples. Within a gap-free stretch each call adds
// the last `step` samples of the window.
float Processor::score_stalta(std::uint64_t start, double& ms) {
    const std::uint64_t end = start + Detector::kWindow;
    std::uint64_t from = stalta_next_;
    if (from == kUnset || from < start) {
        stalta_->reset();
        from = start;
    }
    // Detections carry a window start 3.5 s before P, the convention of the
    // detector windows, so that the picker and magnitude windows and the
    // catalogue matching are placed as for the model.
    const auto onset = static_cast<std::uint64_t>(3.5 * kFs);
    const auto t0 = Clock::now();
    double peak = 0;
    stalta_triggers_.clear();
    for (std::uint64_t pos = from; pos < end; ++pos) {
        const double r = stalta_->step(std::span<const double, 3>(raw_.data() + (pos - start) * 3, 3));
        peak = std::max(peak, r);
        if (cfg_.from > 0 && pos_to_epoch(pos) < cfg_.from) {
            active_ = r >= cfg_.stalta_off;
        } else if (!active_) {
            if (r >= cfg_.stalta_on && pos_to_epoch(pos - onset) - last_trigger_ >= cfg_.retrigger_seconds) {
                stalta_triggers_.emplace_back(pos, r);
                active_ = true;
                last_trigger_ = pos_to_epoch(pos - onset);
            }
        } else if (r < cfg_.stalta_off) {
            active_ = false;
        }
    }
    ms = ms_since(t0);
    stalta_next_ = end;
    for (const auto& [pos, r] : stalta_triggers_) trigger(pos - onset, pos_to_epoch(pos + 1), static_cast<float>(r), ms);
    return static_cast<float>(peak);
}

// Sends a detection dated `run_start` (P assumed 3.5 s later) and schedules the
// picker window and the early magnitude window for it. The caller updates the
// trigger state.
void Processor::trigger(std::uint64_t run_start, double declared_at, float p, double ms) {
    const double run_t = pos_to_epoch(run_start);
    ++stats_.detections;
    bus_.send(Detection{st_.code, run_t, declared_at, p, ms});
    if (cfg_.pick) picks_.emplace_back(run_start - static_cast<std::uint64_t>(cfg_.picker_lead * kFs), run_t);
    // Early magnitude window: the regressor's window starts 2 s before the
    // assumed P.
    if (cfg_.magnitude) early_.emplace_back(run_start + 150, run_t);
}

// Runs the scheduled picker and early magnitude jobs whose windows end at or
// before `limit`. The run loop passes the end of the last scored window, not
// the end of the data ingest has written: in a replay faster than real time
// ingest can be far ahead, and the order of the messages would then depend on
// thread timing. With live data the two are the same.
void Processor::run_jobs(std::uint64_t limit) {
    while (!picks_.empty() && picks_.front().first + Picker::kWindow <= limit) {
        const auto [start, trigger] = picks_.front();
        picks_.pop_front();
        run_pick(start, trigger);
    }
    while (!early_.empty() && early_.front().first + kMagWindow <= limit) {
        const auto [start, trigger] = early_.front();
        early_.pop_front();
        run_early_magnitude(start, trigger);
    }
}

void Processor::run_pick(std::uint64_t start, double trigger) {
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
    Pick pick{st_.code, trigger, p_time, t + pk.s_seconds, pk.p_prob, pk.s_prob, declared, ms, nullptr};

    // If the P pick is confident and near the trigger, the pick carries the 10 s
    // window starting 2 s before it (the training alignment) and the current
    // noise baseline. The network stage estimates the magnitude from them only
    // if the trigger belongs to a declared event; most triggers do not (single-
    // station detections, S and coda re-triggers).
    if (cfg_.magnitude && pk.p_prob >= 0.5 && p_time >= trigger - 3 && p_time <= trigger + 9) {
        const std::uint64_t mstart = epoch_to_pos(p_time) - 200;
        auto w = std::make_shared<MagnitudeWindow>();
        w->raw.resize(kMagWindow * 3);
        if (extract(mstart, kMagWindow, w->raw, resume)) {
            w->start = pos_to_epoch(mstart);
            w->noise = noise_.noise();
            pick.magnitude_window = std::move(w);
        }
    }
    bus_.send(std::move(pick));
}

// Early magnitude estimate for every trigger, so that an alarm can carry a
// magnitude before any pick exists.
void Processor::run_early_magnitude(std::uint64_t start, double trigger) {
    std::uint64_t resume = 0;
    if (!extract(start, kMagWindow, raw_, resume)) return;
    const auto t0 = Clock::now();
    const float m = magnitude_->estimate(std::span<const double>(raw_.data(), kMagWindow * 3), noise_.noise());
    const double ms = ms_since(t0);
    ++stats_.magnitudes;
    stats_.magnitude_ms.push_back(static_cast<float>(ms));
    const auto& nz = noise_.noise();
    bus_.send(MagnitudeEstimate{st_.code, trigger, false, pos_to_epoch(start), m, nz.ready ? nz.windows : 0,
                                pos_to_epoch(start + kMagWindow), ms});
}

// At most every `noise_every` seconds, adds the 10 s ending at the newest scored
// window to the noise baseline, provided every detector window overlapping it
// scored below `noise_below` (`stalta_noise_below` for STA/LTA). Signal included in the baseline would inflate the
// noise sigma and lower later magnitude estimates.
void Processor::maybe_add_noise(std::uint64_t window_start) {
    const std::uint64_t end = window_start + Detector::kWindow;
    if (end < kMagWindow + Detector::kWindow || active_) return;
    const std::uint64_t begin = end - kMagWindow;
    if (end < last_noise_end_ + static_cast<std::uint64_t>(cfg_.noise_every * kFs)) return;
    if (recent_.empty() || recent_.front().first + Detector::kWindow > begin) return;   // not enough scored history
    const double below = stalta_ ? cfg_.stalta_noise_below : cfg_.noise_below;
    for (const auto& [ws, p] : recent_)
        if (ws + Detector::kWindow > begin && p >= below) return;
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
            // STA/LTA also needs its long-term average before `from`.
            if (cfg_.from > 0 && !cfg_.magnitude && !stalta_) first = std::max(first, epoch_to_pos(cfg_.from));
            next_ = align_up(first, cfg_.step);
        }

        // At most 256 windows (about 2 minutes of data) per pass, so that the
        // ring floor below is raised regularly. Ingest may hold a full ring
        // (22 minutes) ahead; scoring all of it before raising the floor would,
        // for a processor slower than ~20x real time, block ingest past its
        // backpressure limit and drop data.
        bool progressed = false;
        for (int n = 0; n < 256 && next_ + Detector::kWindow <= available && !stop.load(std::memory_order_relaxed); ++n) {
            const std::uint64_t end = next_ + Detector::kWindow;
            score_window(next_);
            progressed = true;
            run_jobs(end);
        }

        // Raise the ring floor, keeping the picker lead, one noise window before
        // the next window, and the windows of scheduled jobs.
        const std::uint64_t history = std::max<std::uint64_t>(lead, kMagWindow);
        std::uint64_t keep = next_ > history ? next_ - history : 0;
        if (!picks_.empty()) keep = std::min(keep, picks_.front().first);
        if (!early_.empty()) keep = std::min(keep, early_.front().first);
        for (auto& cs : st_.comp) {
            const auto base = cs.base.load(std::memory_order_acquire);
            if (keep > base) cs.ring.set_floor(std::min<std::uint64_t>(keep - base, cs.ring.written()));
        }

        // Lower bound on declared_at of any later message: the end of the next
        // detector window (for STA/LTA, the next sample to be fed) or of a
        // scheduled job's window. After a gap, the STA/LTA restarts at a later
        // window start.
        double until = pos_to_epoch(next_ + Detector::kWindow);
        if (stalta_) until = pos_to_epoch(stalta_next_ != kUnset ? stalta_next_ : next_);
        if (!picks_.empty()) until = std::min(until, pos_to_epoch(picks_.front().first + Picker::kWindow));
        if (!early_.empty()) until = std::min(until, pos_to_epoch(early_.front().first + kMagWindow));
        if (until - last_progress_ >= 0.5) {
            bus_.send(Progress{st_.code, until});
            last_progress_ = until;
        }

        if (!progressed) {
            // `done` was read before `available`, so all ingested data has been
            // processed.
            if (done) {
                // End of the data: run the jobs whose windows fit in it.
                run_jobs(available);
                stats_.abandoned_picks += picks_.size();
                picks_.clear();
                early_.clear();
                break;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(2));
        }
    }
    bus_.send(StationDone{st_.code});
}

}  // namespace ayzek::pipeline
