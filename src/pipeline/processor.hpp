#pragma once

// Processing stage, one thread per station: sliding windows -> preprocessing ->
// detector ensemble (or STA/LTA, FilterPicker, or the streaming transformer) ->
// trigger. A trigger schedules a 60 s picker window and, with a Pd model, the
// peak-displacement measurements at its onset (pd_magnitude.hpp).

#include "common.hpp"
#include "dsp.hpp"
#include "filterpicker.hpp"
#include "models.hpp"
#include "pd_magnitude.hpp"
#include "stalta.hpp"
#include "station.hpp"
#include "transformer.hpp"
#include "trigger.hpp"

#include <atomic>
#include <cstdint>
#include <deque>
#include <fstream>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace ayzek::pipeline {

// Model: the 3-seed 6 s window detector. Transformer: the streaming onset
// transformer (transformer.hpp), one output per 0.1 s. StaLta: the reference
// STA/LTA trigger. FilterPicker: Lomax's FilterPicker5 (filterpicker.hpp).
enum class DetectorKind { Model, StaLta, Transformer, FilterPicker };

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
  // Magnitude (the Pd estimate, `pd` below). Off, a replay with `from` also
  // starts reading at `from` where no STA/LTA or FilterPicker needs the earlier
  // data, as it did when this switched off the former regressor.
  bool magnitude = true;
  double from = 0; // no triggers before this epoch

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

  // FilterPicker (detector = FilterPicker) on the vertical component: every
  // pick triggers, dated at the pick, and is declared at the sample that
  // accepted it. Defaults: Lomax et al. (2012) at 100 Hz (18-filterpicker.md).
  dsp::FilterPickerConfig fp{};

  // Anchoring (detector = Model): an STA/LTA runs alongside the detector and
  // its onsets give the P time of a trigger, in place of the assumption that P
  // lies 3.5 s into the first window of the run. It never triggers anything:
  // the detector alone decides (09-sta-lta.md).
  bool anchor = true;
  double anchor_on = 3.0;      // ratio taken as an onset
  double anchor_off = 1.5;     // ratio below which the next onset can be taken
  bool require_onset = false;  // a trigger also needs an onset in that range
  // The picker window is placed from the detector's run start even when the
  // trigger is anchored, unless this is set. The picker was trained on
  // windows placed that way, and moving the window changes what it sees
  // (09-sta-lta.md).
  bool anchor_picker = false;
  // 6 s detector: a FilterPicker at its default settings runs alongside on
  // the vertical, and a trigger's onset (Detection::p_time, which dates the
  // P-time location and the Pd measurement) is its pick nearest the
  // trigger's P within onset_fp_window s, or that P without one. It never
  // triggers anything and leaves the detection's window start alone
  // (17-pd-magnitude.md).
  bool onset_fp = true;
  double onset_fp_window = 3.0;

  // Transformer (detector = Transformer): trigger when a token's probability
  // reaches `threshold`, re-arm when it falls below `release`. The detection
  // dates P at the token's end minus its dt output. main() sets both from the
  // model's validation operating point unless given on the command line.
  // With `dt_reset`, a token with p still at or above the threshold also
  // triggers when dt restarts: back to `dt_reset_below` s or less for
  // `dt_reset_tokens` tokens after reaching `dt_reset_from` s, i.e. a new
  // onset inside the coda of the last one (trigger.hpp). Other detectors
  // ignore these.
  // Transformer: the picker looks for P within `pick_p_tol` s of the
  // transformer's own P, and for S only before this station's next trigger,
  // so a later event in the 60 s window is not picked (models.hpp,
  // Picker::Search). Other detectors search the whole window.
  bool pick_on_onset = true;
  double pick_p_tol = 3.0;
  bool dt_reset = true;
  float dt_reset_below = 2.0f;
  float dt_reset_from = 5.0f;
  std::size_t dt_reset_tokens = 2;
  //
  // With the geometry head and `geometry` set, each trigger is followed by a
  // StationGeometry at the trigger token and every `geometry_every` seconds
  // after it, until `geometry_seconds` after P or until p falls below
  // `release`. The estimate sharpens once S is inside the model's ~32 s
  // lookback, which covers S-P times up to ~20 s (~160 km).
  bool geometry = false;
  double geometry_every = 1.0;
  double geometry_seconds = 20.0;
  // With a Pd model (pd_magnitude.hpp) and a gain for the station, the raw
  // vertical component runs through the displacement chain, and each trigger
  // is followed by a PdEstimate as each of the model's windows after the
  // trigger's P ends: the transformer's dated P, the anchoring onset (or 3.5 s
  // into the window) for the window detector, the onset or pick for STA/LTA
  // and FilterPicker. A value needs the chain to have run without a gap from
  // kPdLeadSeconds before P to the end of its window.
  const PdModel *pd = nullptr;
};

// Displacement run needed before P: the filters' run-in and the 10 s noise
// window ending 1 s before P, as in onset's pd.py.
inline constexpr double kPdLeadSeconds = 41.0;

struct ProcessorStats {
  std::uint64_t windows = 0, gap_windows = 0, detections = 0, picks = 0,
                abandoned_picks = 0;
  std::uint64_t anchored =
      0; // detections whose P time came from an STA/LTA onset
  std::uint64_t unconfirmed =
      0; // runs of windows dropped for want of an onset (require_onset)
  std::uint64_t context_refreshes = 0; // transformer: station context updates
  std::uint64_t geometry = 0; // transformer: StationGeometry messages sent
  std::uint64_t pd_values = 0; // PdEstimate messages sent
  std::uint64_t onset_picks = 0; // onset_fp: detections dated by a pick
  std::uint64_t dt_resets = 0; // transformer: detections from a dt restart
  std::vector<float>
      window_ms; // conditioning + ensemble (or STA/LTA update), per window
  std::vector<float> pick_ms;
};

class Processor {
public:
  Processor(Station &st, const std::vector<Weights> &detector,
            const Weights &picker, const dsp::Bandpass &bp,
            ProcessorConfig cfg, Bus &bus,
            const std::string &scores_path,
            const std::string &scores_in_path = "",
            const Weights *transformer = nullptr);
  void run(const std::atomic<bool> &stop);
  [[nodiscard]] const ProcessorStats &stats() const noexcept { return stats_; }

private:
  bool extract(std::uint64_t start, std::size_t n, std::vector<double> &out,
               std::uint64_t &resume);
  void score_window(std::uint64_t start);
  // Transformer: feeds positions [next_, limit) to the stream, at most `max`.
  bool feed_stream(std::uint64_t limit, std::size_t max);
  void on_token(const OnsetStream::Token &tok, std::uint64_t fed, double ms);
  float score_model(std::uint64_t start, double &ms);
  float score_stalta(std::uint64_t start, double &ms);
  float update_stalta(std::uint64_t start, double on, double off, double &ms);
  float score_filterpicker(std::uint64_t start, double &ms);
  // The onset position for a trigger, or kUnset if there is none in range.
  [[nodiscard]] std::uint64_t onset_for(std::uint64_t run_start) const;
  void trigger(std::uint64_t run_start, std::uint64_t p_pos,
               std::uint64_t model_start, double declared_at, float p,
               double ms, bool restart = false, std::uint64_t onset = kUnset);
  // onset_fp: the window's new vertical samples into the onset FilterPicker,
  // and the pick nearest `p_pos` within the window, or kUnset.
  void update_onset_fp(std::uint64_t start);
  [[nodiscard]] std::uint64_t onset_pick(std::uint64_t p_pos) const;
  void run_jobs(std::uint64_t limit);
  void run_pick(std::uint64_t start, double trigger, std::uint64_t p_pos);
  // Pd: one vertical sample into the displacement chain, and the
  // measurements whose windows have ended by `pos`.
  void pd_push(std::uint64_t pos, std::optional<double> z);
  void pd_measure(std::uint64_t pos);
  // Window detectors: the window's new vertical samples into the chain.
  void pd_feed_window(std::uint64_t start);

  Station &st_;
  ProcessorConfig cfg_;
  Bus &bus_;
  std::ofstream scores_;
  std::unordered_map<std::uint64_t, float>
      cached_; // window start -> probability, from --scores-in
  DetectorEnsemble detector_;
  std::unique_ptr<Picker> picker_;
  std::unique_ptr<dsp::StaLta> stalta_;
  std::unique_ptr<dsp::FilterPicker> fp_;
  std::uint64_t fp_next_ = kUnset; // next sample to feed; kUnset after a gap
  std::uint64_t fp_origin_ = 0;    // position of the FilterPicker's sample 0
  // onset_fp: the FilterPicker alongside the 6 s detector, as fp_ above, and
  // the positions of its picks of the last 120 s.
  std::unique_ptr<dsp::FilterPicker> fpo_;
  std::uint64_t fpo_next_ = kUnset, fpo_origin_ = 0;
  std::deque<std::uint64_t> fpo_picks_;
  std::unique_ptr<OnsetStream> onset_;
  std::vector<OnsetStream::Token> tokens_;
  std::optional<TokenTrigger> token_trigger_; // transformer
  // Geometry after the last trigger: next token end to report, the last
  // position to report at (0: not tracking), P and the detection's window.
  std::uint64_t geo_next_ = 0, geo_until_ = 0, geo_p_ = 0;
  double geo_trigger_ = 0;
  // Pd: the chain, |displacement| (counts) of the last kPdRing positions,
  // where its current run began, and the measurements still to send: P, the
  // detection's window, and the index of the next window.
  static constexpr std::size_t kPdRing = 6400;
  std::optional<Displacement> pd_chain_;
  std::vector<double> pd_ring_;
  std::uint64_t pd_run_start_ = 0;
  std::uint64_t pd_next_ = kUnset; // window detectors: next sample to feed
  struct PdPending {
    std::uint64_t p;
    double trigger_window;
    std::size_t next;
  };
  std::deque<PdPending> pd_pending_;
  std::uint64_t stalta_next_ =
      kUnset; // next sample to feed to the STA/LTA; kUnset after a gap
  bool stalta_armed_ =
      true; // the ratio has been below `off` since the last onset
  std::deque<std::pair<std::uint64_t, double>>
      onsets_;               // recent (sample, ratio) of STA/LTA onsets
  std::deque<float> ratios_; // STA/LTA ratio of the last 60 s, for anchoring
  std::uint64_t ratios_end_ = 0; // position after the last of them
  dsp::Conditioner cond6_, cond60_;
  std::vector<double> raw_;
  std::vector<float> standardized_, planar60_;
  std::vector<std::int32_t> tmp_;

  std::uint64_t next_ = kUnset;
  bool active_ = false;
  std::size_t below_ = 0;
  std::size_t above_ = 0;       // consecutive windows at or above the threshold
  std::uint64_t run_start_ = 0; // start of the first of those windows
  double last_trigger_ = -1e18;
  // Scheduled picker windows, in start order.
  struct PickJob {
    std::uint64_t start;  // window start position
    double trigger;       // window_start of the detection it follows
    std::uint64_t p_pos;  // the detector's P, kUnset: search the whole window
  };
  std::deque<PickJob> picks_;
  // Transformer: P positions of this station's accepted triggers, newest
  // last, for ending the S search of an earlier pick at the next onset.
  std::deque<std::uint64_t> onsets_p_;
  double last_progress_ = 0;
  ProcessorStats stats_;
};

} // namespace ayzek::pipeline
