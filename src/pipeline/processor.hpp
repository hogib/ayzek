#pragma once

// Processing stage, one thread per station: sliding windows -> preprocessing ->
// detector ensemble (or STA/LTA) -> trigger. A trigger schedules an early
// magnitude estimate and a 60 s picker window. A confident P pick carries the
// 10 s window for a magnitude estimate at the picked P, which the network stage
// runs only for declared events. Windows scored as noise update the station's
// noise baseline used by the magnitude regressor.

#include "common.hpp"
#include "dsp.hpp"
#include "magnitude.hpp"
#include "models.hpp"
#include "stalta.hpp"
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

enum class DetectorKind { Model, StaLta };

struct ProcessorConfig {
  DetectorKind detector = DetectorKind::Model;
  std::size_t step = 50; // samples between window starts (0.5 s)
  // Outputs saturate at 0.90-0.91 (label smoothing 0.1/0.9), so a single-window
  // threshold at 0.9 depends on the third decimal. A lower threshold held for
  // several windows, or one window at the plateau, triggers (05-pipeline.md).
  float threshold = 0.8f; // threshold for the persistence rule
  std::size_t trigger_windows =
      8; // consecutive windows at or above `threshold` required to trigger (3.5
         // s after the first)
  float instant_threshold =
      0.9f; // a single window at or above this triggers (> 1 disables)
  float release =
      0.3f; // trigger resets after `release_windows` windows below this
  std::size_t release_windows = 2;
  double retrigger_seconds = 15.0; // minimum time between triggers
  double picker_lead =
      2.0; // picker window start, seconds before the triggering window
  bool pick = true;
  bool magnitude = true;
  double noise_every = 30.0; // seconds between noise-baseline windows
  float noise_below =
      0.3f;        // max detector probability over a window accepted as noise
  double from = 0; // no triggers before this epoch; earlier windows update the
                   // noise baseline only

  // STA/LTA (detector = StaLta): trigger when the ratio reaches `stalta_on`,
  // re-arm when it falls below `stalta_off`. The P onset is taken to be the
  // trigger sample. Defaults: the setting with the most detections at the
  // model's number of unmatched alarms on the tuning data (09-sta-lta.md).
  dsp::StaLtaConfig stalta{.sta_seconds = 1.0,
                           .lta_seconds = 30.0,
                           .f_lo = 2.0,
                           .f_hi = 20.0,
                           .three_component = false,
                           .fs = 100.0};
  double stalta_on = 8.0;
  double stalta_off = 1.5;
  double stalta_noise_below = 2.0; // max ratio over a window accepted as noise

  // Anchoring (detector = Model): an STA/LTA runs alongside the detector and
  // its onsets give the P time of a trigger, in place of the assumption that P
  // lies 3.5 s into the first window of the run. It never triggers anything:
  // the detector alone decides (09-sta-lta.md).
  bool anchor = true;
  double anchor_on = 3.0;      // ratio taken as an onset
  double anchor_off = 1.5;     // ratio below which the next onset can be taken
  bool require_onset = false;  // a trigger also needs an onset in that range
  double magnitude_lead = 2.0; // magnitude window start, seconds before P
  // The picker and early magnitude windows are placed from the detector's run
  // start even when the trigger is anchored, unless these are set. Both models
  // were trained on windows placed that way, and moving a window changes what
  // they see (09-sta-lta.md).
  bool anchor_picker = false;
  bool anchor_magnitude = false;
};

// The 10 s window starting 2 s before a picked P and the station's noise
// baseline at that time, sent with the Pick (common.hpp).
struct MagnitudeWindow {
  double start = 0;        // epoch of the first sample
  std::vector<double> raw; // (1000, 3) interleaved counts
  StationNoise noise;
};

struct ProcessorStats {
  std::uint64_t windows = 0, gap_windows = 0, detections = 0, picks = 0,
                abandoned_picks = 0, magnitudes = 0;
  std::uint64_t anchored =
      0; // detections whose P time came from an STA/LTA onset
  std::uint64_t unconfirmed =
      0; // runs of windows dropped for want of an onset (require_onset)
  std::size_t noise_windows = 0;
  std::vector<float>
      window_ms; // conditioning + ensemble (or STA/LTA update), per window
  std::vector<float> pick_ms,
      magnitude_ms; // magnitude_ms: early estimates only
};

class Processor {
public:
  Processor(Station &st, const std::vector<Weights> &detector,
            const Weights &picker, const std::vector<Weights> &magnitude,
            const dsp::Bandpass &bp, ProcessorConfig cfg, Bus &bus,
            const std::string &scores_path,
            const std::string &scores_in_path = "");
  void run(const std::atomic<bool> &stop);
  [[nodiscard]] const ProcessorStats &stats() const noexcept { return stats_; }

private:
  bool extract(std::uint64_t start, std::size_t n, std::vector<double> &out,
               std::uint64_t &resume);
  void score_window(std::uint64_t start);
  float score_model(std::uint64_t start, double &ms);
  float score_stalta(std::uint64_t start, double &ms);
  float update_stalta(std::uint64_t start, double on, double off, double &ms);
  // The onset position for a trigger, or kUnset if there is none in range.
  [[nodiscard]] std::uint64_t onset_for(std::uint64_t run_start) const;
  void trigger(std::uint64_t run_start, std::uint64_t p_pos,
               std::uint64_t model_start, double declared_at, float p,
               double ms);
  void run_jobs(std::uint64_t limit);
  void run_pick(std::uint64_t start, double trigger);
  void run_early_magnitude(std::uint64_t start, double trigger);
  void maybe_add_noise(std::uint64_t window_start);

  Station &st_;
  ProcessorConfig cfg_;
  Bus &bus_;
  std::ofstream scores_;
  std::unordered_map<std::uint64_t, float>
      cached_; // window start -> probability, from --scores-in
  DetectorEnsemble detector_;
  std::unique_ptr<Picker> picker_;
  std::unique_ptr<MagnitudeEstimator> magnitude_;
  std::unique_ptr<dsp::StaLta> stalta_;
  std::uint64_t stalta_next_ =
      kUnset; // next sample to feed to the STA/LTA; kUnset after a gap
  bool stalta_armed_ =
      true; // the ratio has been below `off` since the last onset
  std::deque<std::pair<std::uint64_t, double>>
      onsets_;               // recent (sample, ratio) of STA/LTA onsets
  std::deque<float> ratios_; // STA/LTA ratio of the last 60 s, for anchoring
  std::uint64_t ratios_end_ = 0; // position after the last of them
  NoiseBaseline noise_;
  dsp::Conditioner cond6_, cond60_;
  std::vector<double> raw_, noise_raw_;
  std::vector<float> standardized_, planar60_;
  std::vector<std::int32_t> tmp_;

  std::uint64_t next_ = kUnset;
  bool active_ = false;
  std::size_t below_ = 0;
  std::size_t above_ = 0;       // consecutive windows at or above the threshold
  std::uint64_t run_start_ = 0; // start of the first of those windows
  double last_trigger_ = -1e18;
  // Scheduled picker and early magnitude windows: (start position, window start
  // of the detection that scheduled it), in start order.
  std::deque<std::pair<std::uint64_t, double>> picks_, early_;
  std::deque<std::pair<std::uint64_t, float>>
      recent_; // (window start, probability)
  std::uint64_t last_noise_end_ = 0;
  double last_progress_ = 0;
  ProcessorStats stats_;
};

} // namespace ayzek::pipeline
