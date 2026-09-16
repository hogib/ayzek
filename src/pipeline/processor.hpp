#pragma once

// Processing stage, one thread per station: sliding windows -> preprocessing ->
// detector ensemble -> trigger. A trigger schedules an early magnitude estimate,
// a 60 s picker window, and a magnitude estimate at the picked P. Windows scored
// as noise update the station's noise baseline used by the magnitude regressor.

#include "common.hpp"
#include "dsp.hpp"
#include "magnitude.hpp"
#include "models.hpp"
#include "station.hpp"

#include <atomic>
#include <cstdint>
#include <deque>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace ayzek::pipeline {

struct ProcessorConfig {
    std::size_t step = 50;             // samples between window starts (0.5 s)
    // Outputs saturate at 0.90-0.91 (label smoothing 0.1/0.9), so a single-window
    // threshold at 0.9 depends on the third decimal. A lower threshold held for
    // several windows, or one window at the plateau, triggers (05-pipeline.md).
    float threshold = 0.8f;            // threshold for the persistence rule
    std::size_t trigger_windows = 8;   // consecutive windows at or above `threshold` required to trigger (3.5 s after the first)
    float instant_threshold = 0.9f;    // a single window at or above this triggers (> 1 disables)
    float release = 0.3f;              // trigger resets after `release_windows` windows below this
    std::size_t release_windows = 2;
    double retrigger_seconds = 15.0;   // minimum time between triggers
    double picker_lead = 2.0;          // picker window start, seconds before the triggering window
    bool pick = true;
    bool magnitude = true;
    double noise_every = 30.0;         // seconds between noise-baseline windows
    float noise_below = 0.3f;          // max detector probability over a window accepted as noise
    double from = 0;                   // no triggers before this epoch; earlier windows update the noise baseline only
};

struct ProcessorStats {
    std::uint64_t windows = 0, gap_windows = 0, detections = 0, picks = 0, abandoned_picks = 0, magnitudes = 0;
    std::size_t noise_windows = 0;
    std::vector<float> window_ms;      // conditioning + ensemble, per window
    std::vector<float> pick_ms, magnitude_ms;
};

class Processor {
public:
    Processor(Station& st, const std::vector<Weights>& detector, const Weights& picker,
              const std::vector<Weights>& magnitude, const dsp::Bandpass& bp, ProcessorConfig cfg, Bus& bus,
              const std::string& scores_path, const std::string& scores_in_path = "");
    void run(const std::atomic<bool>& stop);
    [[nodiscard]] const ProcessorStats& stats() const noexcept { return stats_; }

private:
    bool extract(std::uint64_t start, std::size_t n, std::vector<double>& out, std::uint64_t& resume);
    void score_window(std::uint64_t start);
    void try_pick(std::uint64_t available);
    void try_early_magnitude(std::uint64_t available);
    void run_magnitude(std::uint64_t start, double trigger, bool at_pick, double declared_at);
    void maybe_add_noise(std::uint64_t window_start);

    Station& st_;
    ProcessorConfig cfg_;
    Bus& bus_;
    std::ofstream scores_;
    std::unordered_map<std::uint64_t, float> cached_;   // window start -> probability, from --scores-in
    DetectorEnsemble detector_;
    std::unique_ptr<Picker> picker_;
    std::unique_ptr<MagnitudeEstimator> magnitude_;
    NoiseBaseline noise_;
    dsp::Conditioner cond6_, cond60_;
    std::vector<double> raw_, noise_raw_;
    std::vector<float> standardized_, planar60_;
    std::vector<std::int32_t> tmp_;

    std::uint64_t next_ = kUnset;
    bool active_ = false;
    std::size_t below_ = 0;
    std::size_t above_ = 0;                  // consecutive windows at or above the threshold
    std::uint64_t run_start_ = 0;            // start of the first of those windows
    double last_trigger_ = -1e18;
    std::uint64_t pending_pick_ = kUnset;   // start position of a scheduled picker window
    double pending_trigger_ = 0;             // window start of the detection that scheduled it
    std::uint64_t early_mag_ = kUnset;      // start position of a scheduled early magnitude window
    double early_trigger_ = 0;
    std::deque<std::pair<std::uint64_t, float>> recent_;   // (window start, probability)
    std::uint64_t last_noise_end_ = 0;
    double last_progress_ = 0;
    ProcessorStats stats_;
};

}  // namespace ayzek::pipeline
