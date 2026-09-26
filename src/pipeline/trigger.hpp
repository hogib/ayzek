#pragma once

// The trigger rule of the streaming transformer (--detector transformer), one
// decision per 0.1 s token from its two outputs:
//
//   p   an event is under way
//   dt  seconds since it began
//
// Rising edge: the first token with p >= threshold after p last fell below
// `release`. On its own that misses most of an aftershock sequence: p stays
// high between events, never falls below `release`, and the next onset has no
// rising edge.
//
// dt restart: while p stays at or above the threshold, a token whose dt is
// back at or below `dt_reset_below` s, after dt had reached `dt_reset_from` s
// since the last trigger, is a new onset rather than the coda of the old one.
// It fires after `dt_reset_tokens` such tokens in a row, so a single-token dip
// in dt does not. This is what the dt head exists for (onset DESIGN.md §2).
//
// The processor applies the minimum time between triggers on top of this.

#include <algorithm>
#include <cstddef>

namespace ayzek::pipeline {

struct TokenTriggerConfig {
  float threshold = 0.5f;
  float release = 0.25f;
  bool dt_reset = true;
  float dt_reset_below = 2.0f; // dt, seconds, that counts as a fresh onset
  float dt_reset_from = 5.0f;  // dt the event must have reached before
  std::size_t dt_reset_tokens = 2;
};

class TokenTrigger {
public:
  enum class Fire { None, Rising, DtReset };

  explicit TokenTrigger(TokenTriggerConfig cfg) : cfg_(cfg) {}

  // One token. During `warmup` nothing fires, but the state advances as if it
  // had, so the first trigger after it is a real rising edge or restart.
  Fire step(float p, float dt, bool warmup) {
    if (p < cfg_.release) {
      armed_ = true;
      peak_ = 0;
      low_ = 0;
      return Fire::None;
    }
    if (p < cfg_.threshold) {
      peak_ = std::max(peak_, dt);
      low_ = 0;
      return Fire::None;
    }
    Fire f = Fire::None;
    if (armed_) {
      f = Fire::Rising;
    } else if (cfg_.dt_reset && peak_ >= cfg_.dt_reset_from &&
               dt <= cfg_.dt_reset_below) {
      if (++low_ >= cfg_.dt_reset_tokens)
        f = Fire::DtReset;
    } else {
      low_ = 0;
    }
    armed_ = false;
    if (f == Fire::None) {
      peak_ = std::max(peak_, dt);
      return f;
    }
    peak_ = dt; // a new event starts here
    low_ = 0;
    return warmup ? Fire::None : f;
  }

  // Not re-armed since the last crossing: an event is under way.
  [[nodiscard]] bool active() const noexcept { return !armed_; }

private:
  TokenTriggerConfig cfg_;
  bool armed_ = true;
  float peak_ = 0;       // largest dt since the last trigger
  std::size_t low_ = 0;  // consecutive tokens with a restarted dt
};

} // namespace ayzek::pipeline
