#pragma once

// Stage 2 per station: sliding windows -> conditioning -> detector ensemble ->
// trigger -> on trigger, an early magnitude, the 60 s picker window, and a
// magnitude at the picked P. Quiet windows feed the station's noise baseline,
// which the magnitude regressor normalises against.

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
#include <vector>

namespace ayzek::pipeline {

struct ProcessorConfig {
    std::size_t step = 50;             // samples between window starts (0.5 s)
    float threshold = 0.9f;            // trigger on; probabilities top out near 0.91 (label smoothing)
    float release = 0.3f;              // trigger off, held for `release_windows`
    std::size_t release_windows = 2;
    double retrigger_seconds = 15.0;   // ignore new triggers this soon after one
    double picker_lead = 2.0;          // picker window starts this long before the triggering window
    bool pick = true;
    bool magnitude = true;
    double noise_every = 30.0;         // seconds between noise-baseline windows
    float noise_below = 0.3f;          // every detector window over a noise window must score under this
    double from = 0;                   // no detections before this epoch; earlier windows only warm the noise baseline
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
              const std::string& scores_path);
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
    double last_trigger_ = -1e18;
    std::uint64_t pending_pick_ = kUnset;   // absolute start of a picker window still waiting for data
    double pending_trigger_ = 0;             // the detection window it belongs to
    std::uint64_t early_mag_ = kUnset;      // absolute start of an early magnitude window still waiting
    double early_trigger_ = 0;
    std::deque<std::pair<std::uint64_t, float>> recent_;   // (window start, probability)
    std::uint64_t last_noise_end_ = 0;
    double last_progress_ = 0;
    ProcessorStats stats_;
};

}  // namespace ayzek::pipeline
