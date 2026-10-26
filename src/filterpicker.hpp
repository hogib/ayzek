#pragma once

// FilterPicker (Lomax, Satriano & Vassallo 2012, SRL 83(3), 531-540), updated
// one sample at a time. A port of FilterPicker5.c v5.0.3 (Anthony Lomax,
// GPL-2.0-or-later; picker_20130212.tgz from alomax.free.fr/FilterPicker,
// sha256 3e0943e3...1cc1222), checked against it by tests/test_filterpicker.cpp
// (docs/impl/18-filterpicker.md).
//
// The first difference of the signal runs through a bank of N + 1 bands, N the
// power of two at or above filter_window / dt: band k is two one-pole
// high-passes and a one-pole low-pass with time constant 2^k dt / (2 pi). The
// characteristic function (CF) of a band is (y^2 - mean) / sd, with mean and
// sd of y^2 decaying over long_term_window, and the CF of a sample is the
// maximum over the bands. A pick is triggered when the CF reaches threshold1
// and declared once the mean CF (clipped at 5 threshold1) over t_up from the
// trigger reaches threshold2; the pick time is where the CF of the triggering
// band last rose past its uncertainty threshold. No pick is triggered before
// long_term_window has passed, and after a pick none until the CF of every
// band has fallen to 2 or below.
//
// One difference from the C function: it starts its first difference from the
// mean of the first long_term_window of samples of the first packet, which a
// stream does not yet have. Here it starts from the first sample, as the C
// function does when its first packet holds one sample.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace ayzek::dsp {

struct FilterPickerConfig {
  // Defaults: Lomax et al. (2012) Table 1 at 100 Hz (300, 500 and 20 sample
  // intervals). Earthworm's pick_FP doubles the three windows at 100 Hz.
  double filter_window = 3.0;    // s
  double long_term_window = 5.0; // s
  double threshold1 = 10.0;
  double threshold2 = 10.0;
  double t_up = 0.2; // s
  double fs = 100.0;
};

// Sample indices count from the first sample after a reset.
struct FilterPick {
  std::int64_t time;       // the pick: middle of [begin, end]
  std::int64_t begin, end; // uncertainty interval, at least period / 20 wide
  std::int64_t trigger;    // the CF reached threshold1
  int polarity;            // +1, -1, 0 unknown
  double polarity_weight;  // 0-1
  double amplitude;        // CF at the trigger
  double period;           // s, of the band that triggered
};

class FilterPicker {
public:
  explicit FilterPicker(const FilterPickerConfig &cfg);
  // One vertical-component sample. Returns the pick declared at it, if any:
  // between t_up after its trigger and later, never before.
  std::optional<FilterPick> step(double x);
  void reset();
  // CF of the last sample, clipped at 5 threshold1.
  [[nodiscard]] double cf() const noexcept { return cf_clipped_; }
  // A new pick can be triggered (the CF has fallen back since the last one).
  [[nodiscard]] bool armed() const noexcept { return allow_.has_value(); }
  [[nodiscard]] std::size_t bands() const noexcept { return bands_.size(); }
  [[nodiscard]] const FilterPickerConfig &config() const noexcept {
    return cfg_;
  }

private:
  struct Band {
    double hp = 0, lp = 0, period = 0; // filter constants, period in s
    double y0 = 0, y1 = 0, y2 = 0, y2_last = 0;
    double e = 0, mean = 0, var = 0,
           sd = 0;              // y2^2 and its long-term statistics
    double cf = 0, cf_last = 0; // clipped CF, for the uncertainty rule
    double uncertainty_threshold = 0;
  };

  FilterPickerConfig cfg_;
  double dt_, decay_factor_, decay_const_, max_cf_, critical_integral_;
  std::int64_t enable_after_;
  std::size_t nup_; // samples in the t_up ring buffers

  std::vector<Band> bands_;
  // Per band, nup_ entries each: where the CF last rose past the uncertainty
  // threshold, and the polarity sums since then.
  std::vector<std::int64_t> index_uncertainty_;
  std::vector<double> pol_sum_, pol_abs_;
  // Ring buffers over the last nup_ samples.
  std::vector<double> integral_, cf_clipped_ring_, cf_ring_;
  std::vector<std::size_t> cf_band_;
  std::size_t ptr_ = 0;

  std::int64_t n_ = 0, total_ = -1;
  bool enabled_ = false;
  std::optional<std::int64_t> allow_; // first sample at which a pick may begin
  double last_sample_ = 0, cf_clipped_ = 0;
};

} // namespace ayzek::dsp
