// Port of FilterPicker5.c v5.0.3 and FilterPicker5_Memory.c
// Copyright (C) 2006-2012 Anthony Lomax <anthony@alomax.net www.alomax.net>
// Port: Copyright (C) 2026 the ayzek authors
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version. The original is GPL-2.0-or-later.
//
// The arithmetic follows the C function statement by statement, so that the
// two agree to rounding (tests/test_filterpicker.cpp). Comments give the C
// names where they differ.

#include "filterpicker.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace ayzek::dsp {

namespace {

constexpr double kMinSd = 1.0e-30; // DOUBLE_MIN_VALUE: a dead trace
constexpr double kCriticalPolarityRatio = 0.667;
constexpr double kMinUncertaintyThreshold = 0.5;
constexpr double kMaxAllowNewPick = 2.0; // maxAllowNewPickThreshold

} // namespace

FilterPicker::FilterPicker(const FilterPickerConfig &cfg)
    : cfg_(cfg), dt_(1.0 / cfg.fs), decay_factor_(dt_ / cfg.long_term_window),
      decay_const_(1.0 - decay_factor_), max_cf_(5.0 * cfg.threshold1) {
  // nLongTermWindow: triggering starts after this many samples.
  enable_after_ = 1 + static_cast<std::int64_t>(cfg.long_term_window / dt_);
  // numRecursive: bands up to the power of two at or above filter_window / dt.
  std::size_t nb = 1;
  {
    const auto previous = static_cast<std::int64_t>(cfg.filter_window / dt_);
    std::int64_t n = 1;
    while (n < previous) {
      ++nb;
      n *= 2;
    }
  }
  bands_.resize(nb);
  double window = dt_ / (2.0 * std::numbers::pi);
  for (auto &b : bands_) {
    b.period = window * 2.0 * std::numbers::pi;
    b.lp = dt_ / (window + dt_);
    b.hp = window / (window + dt_);
    window *= 2.0;
  }
  // nTUpEvent, and criticalIntegralCharFunct = nTUpEvent * threshold2.
  nup_ = static_cast<std::size_t>(
      std::max(1, static_cast<int>(0.5 + cfg.t_up / dt_) + 1));
  critical_integral_ = static_cast<double>(nup_) * cfg.threshold2;
  index_uncertainty_.resize(nb * nup_);
  pol_sum_.resize(nb * nup_);
  pol_abs_.resize(nb * nup_);
  integral_.resize(nup_);
  cf_clipped_ring_.resize(nup_);
  cf_ring_.resize(nup_);
  cf_band_.resize(nup_);
  reset();
}

void FilterPicker::reset() {
  for (auto &b : bands_) {
    b.y0 = b.y1 = b.y2 = b.y2_last = 0;
    b.e = b.mean = b.var = b.sd = 0;
    b.cf = b.cf_last = 0;
    b.uncertainty_threshold = cfg_.threshold1 / 2.0;
  }
  std::fill(index_uncertainty_.begin(), index_uncertainty_.end(), 0);
  std::fill(pol_sum_.begin(), pol_sum_.end(), 0.0);
  std::fill(pol_abs_.begin(), pol_abs_.end(), 0.0);
  std::fill(integral_.begin(), integral_.end(), 0.0);
  std::fill(cf_clipped_ring_.begin(), cf_clipped_ring_.end(), 0.0);
  std::fill(cf_ring_.begin(), cf_ring_.end(), 0.0);
  std::fill(cf_band_.begin(), cf_band_.end(), 0);
  ptr_ = 0;
  n_ = 0;
  total_ = -1;
  enabled_ = false;
  allow_.reset();
  last_sample_ = 0;
  cf_clipped_ = 0;
}

std::optional<FilterPick> FilterPicker::step(double x) {
  const std::int64_t n = n_++;
  if (n == 0)
    last_sample_ = x;
  const std::size_t last = ptr_;
  ptr_ = (ptr_ + 1) % nup_;

  double cf = 0, cf_clipped = 0;
  const double diff = x - last_sample_;
  for (std::size_t k = bands_.size(); k-- > 0;) {
    auto &b = bands_[k];
    // Two one-pole high-passes, then a one-pole low-pass.
    double y = b.hp * (b.y0 + diff);
    const double diff2 = y - b.y0;
    b.y0 = y;
    y = b.hp * (b.y1 + diff2);
    b.y1 = y;
    y = b.y2 + b.lp * (y - b.y2);
    b.y2_last = b.y2;
    b.y2 = y;
    b.e = y * y;
    double clipped = 0;
    if (b.sd > kMinSd) {
      const double c = (b.e - b.mean) / b.sd;
      clipped = c;
      // Clipped, and the energy with it, to recover quickly after a large
      // event.
      if (clipped > max_cf_) {
        clipped = max_cf_;
        b.e = max_cf_ * b.sd + b.mean;
      }
      if (c >= cf) {
        cf = c;
        cf_clipped = clipped;
        cf_band_[ptr_] = k;
      }
      // The band of a trigger is the highest frequency at threshold1.
      if (c >= cfg_.threshold1)
        cf_band_[ptr_] = k;
    }
    // Uncertainty and polarity: restart where the CF rises past the band's
    // uncertainty threshold.
    b.cf = clipped;
    const bool up =
        b.cf_last < b.uncertainty_threshold && b.cf >= b.uncertainty_threshold;
    b.cf_last = b.cf;
    const std::size_t at = k * nup_ + ptr_, prev = k * nup_ + last;
    index_uncertainty_[at] = up ? n - 1 : index_uncertainty_[prev];
    if (up) {
      pol_sum_[at] = 0.0;
      pol_abs_[at] = 0.0;
    } else {
      pol_sum_[at] = pol_sum_[prev];
      pol_abs_[at] = pol_abs_[prev];
    }
    const double inc = b.y2 - b.y2_last;
    pol_sum_[at] += inc;
    pol_abs_[at] += std::fabs(inc);
  }
  cf_clipped_ = cf_clipped;

  bool accepted = false;
  std::int64_t trigger = 0, uncertainty_index = 0;
  std::size_t band = 0;
  double amplitude = 0, polarity_weight = 0;
  int polarity = 0;
  // After the long-term window has passed. As in C, the count stops once
  // triggering is enabled.
  if (enabled_ || total_++ > enable_after_) {
    enabled_ = true;
    // The clipped CF summed over the last nup_ samples.
    integral_[ptr_] = integral_[last] - cf_clipped_ring_[ptr_] + cf_clipped;
    cf_clipped_ring_[ptr_] = cf_clipped;
    cf_ring_[ptr_] = cf;

    // Find the last point in the window where the CF rose past threshold1
    // with the clipped CF from there to now reaching the critical integral.
    if (allow_ && integral_[ptr_] >= critical_integral_) {
      std::size_t m = ptr_;
      double sum = cf_clipped_ring_[m];
      for (std::size_t k = 1;
           k < nup_ && n - static_cast<std::int64_t>(k) > *allow_; ++k) {
        m = m == 0 ? nup_ - 1 : m - 1;
        sum += cf_clipped_ring_[m];
        if (cf_ring_[m] < cfg_.threshold1)
          continue;
        const std::size_t l = m == 0 ? nup_ - 1 : m - 1;
        if (cf_ring_[l] >= cfg_.threshold1 || sum < critical_integral_)
          continue;
        accepted = true;
        amplitude = cf_ring_[m];
        band = cf_band_[m];
        trigger = n - static_cast<std::int64_t>(k);
        uncertainty_index = index_uncertainty_[band * nup_ + m];
        // Polarity from the derivative summed since the uncertainty point,
        // one sample past the trigger.
        const std::size_t ip = (m + 1) % nup_;
        const double s = pol_sum_[band * nup_ + ip],
                     a = pol_abs_[band * nup_ + ip], ratio = s / a;
        if (s > 0.0 && ratio > kCriticalPolarityRatio) {
          polarity = 1;
          polarity_weight =
              (ratio - kCriticalPolarityRatio) / (1.0 - kCriticalPolarityRatio);
        } else if (s < 0.0 && -ratio > kCriticalPolarityRatio) {
          polarity = -1;
          polarity_weight = (-ratio - kCriticalPolarityRatio) /
                            (1.0 - kCriticalPolarityRatio);
        }
        allow_.reset();
        break;
      }
    }

    // New picks once the CF of every band is back at kMaxAllowNewPick.
    if (!accepted && !allow_) {
      bool quiet = true;
      for (const auto &b : bands_)
        if (b.cf > kMaxAllowNewPick) {
          quiet = false;
          break;
        }
      if (quiet)
        allow_ = n;
    }
  }

  // Long-term statistics of every band.
  for (auto &b : bands_) {
    b.mean = b.mean * decay_const_ + b.e * decay_factor_;
    const double dev = b.e - b.mean;
    b.var = b.var * decay_const_ + dev * dev * decay_factor_;
    b.sd = std::sqrt(b.var);
    b.uncertainty_threshold =
        b.uncertainty_threshold * decay_const_ + b.cf * decay_factor_;
    if (b.uncertainty_threshold > cfg_.threshold1 / 2.0)
      b.uncertainty_threshold = cfg_.threshold1 / 2.0;
    else if (b.uncertainty_threshold < kMinUncertaintyThreshold)
      b.uncertainty_threshold = kMinUncertaintyThreshold;
  }
  last_sample_ = x;

  if (!accepted)
    return std::nullopt;
  // The pick is the uncertainty point, the interval symmetric about it,
  // widened to at least a twentieth of the band's period.
  FilterPick p{};
  p.begin = uncertainty_index - (trigger - uncertainty_index);
  p.end = trigger;
  p.period = bands_[band].period;
  const double width = dt_ * static_cast<double>(p.end - p.begin);
  if (width < p.period / 20.0) {
    const auto shift =
        static_cast<std::int64_t>(0.5 * (p.period / 20.0 - width) / dt_);
    p.begin -= shift;
    p.end += shift;
  }
  p.time = (p.begin + p.end) / 2;
  p.trigger = trigger;
  p.polarity = polarity;
  p.polarity_weight = polarity_weight;
  p.amplitude = amplitude;
  return p;
}

} // namespace ayzek::dsp
