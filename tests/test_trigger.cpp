// Tests for the transformer's token trigger (pipeline/trigger.hpp): rising
// edges as before, and a dt restart inside a long event firing once, as an
// aftershock would, while the coda, a one-token dip in dt and a disabled rule
// do not.

#include "pipeline/trigger.hpp"

#include <cstdlib>
#include <print>
#include <vector>

#define CHECK(cond)                                                            \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::println(stderr, "CHECK failed: {}  ({}:{})", #cond, __FILE__,       \
                   __LINE__);                                                  \
      std::abort();                                                            \
    }                                                                          \
  } while (0)

using namespace ayzek::pipeline;
using Fire = TokenTrigger::Fire;

namespace {

struct Tok {
  float p, dt;
};

// Noise, then an event from token `onset`, whose dt grows 0.1 s per token up
// to 10 s, with p held at 0.99.
std::vector<Tok> event(std::size_t quiet, std::size_t length) {
  std::vector<Tok> v(quiet, {0.01f, 0.0f});
  for (std::size_t i = 0; i < length; ++i)
    v.push_back({0.99f, std::min(10.0f, 0.1f * static_cast<float>(i))});
  return v;
}

std::vector<std::pair<std::size_t, Fire>> run(TokenTrigger &t,
                                              const std::vector<Tok> &v,
                                              std::size_t warmup = 0) {
  std::vector<std::pair<std::size_t, Fire>> out;
  for (std::size_t i = 0; i < v.size(); ++i)
    if (auto f = t.step(v[i].p, v[i].dt, i < warmup); f != Fire::None)
      out.emplace_back(i, f);
  return out;
}

const TokenTriggerConfig kCfg{.threshold = 0.9f, .release = 0.45f};

void rising_edge_and_coda() {
  TokenTrigger t(kCfg);
  const auto f = run(t, event(50, 400)); // 40 s of coda with dt at 10 s
  CHECK(f.size() == 1 && f[0].first == 50 && f[0].second == Fire::Rising);
  CHECK(t.active());
}

void aftershock_in_the_coda() {
  // A second onset 30 s into the first event: p never drops, dt restarts.
  auto v = event(50, 300);
  for (std::size_t i = 0; i < 200; ++i)
    v.push_back({0.99f, std::min(10.0f, 0.1f * static_cast<float>(i))});
  TokenTrigger t(kCfg);
  const auto f = run(t, v);
  CHECK(f.size() == 2);
  CHECK(f[0].second == Fire::Rising && f[0].first == 50);
  // dt 0.0 then 0.1: fires on the second restarted token.
  CHECK(f[1].second == Fire::DtReset && f[1].first == 351);

  TokenTrigger off({.threshold = 0.9f, .release = 0.45f, .dt_reset = false});
  CHECK(run(off, v).size() == 1);
}

void one_token_dip_does_not_fire() {
  auto v = event(50, 300);
  v[200].dt = 0.3f;
  TokenTrigger t(kCfg);
  CHECK(run(t, v).size() == 1);
}

void early_restart_does_not_fire() {
  // dt restarting before it reached 5 s is the same onset re-dated, not a new
  // event.
  auto v = event(50, 30); // dt reaches 2.9 s
  for (std::size_t i = 0; i < 100; ++i)
    v.push_back({0.99f, 0.1f * static_cast<float>(i)});
  TokenTrigger t(kCfg);
  CHECK(run(t, v).size() == 1);
}

void restart_needs_p_above_threshold() {
  auto v = event(50, 300);
  for (std::size_t i = 0; i < 50; ++i)
    v.push_back({0.8f, 0.1f * static_cast<float>(i)}); // under way, below threshold
  TokenTrigger t(kCfg);
  CHECK(run(t, v).size() == 1);
}

void release_rearms() {
  auto v = event(50, 100);
  const auto second = event(20, 100);
  v.insert(v.end(), second.begin(), second.end());
  TokenTrigger t(kCfg);
  const auto f = run(t, v);
  CHECK(f.size() == 2 && f[1].first == 170 && f[1].second == Fire::Rising);
}

void warmup_consumes_without_firing() {
  TokenTrigger t(kCfg);
  CHECK(run(t, event(50, 300), 60).empty()); // crossing at 50 is in warm-up
}

} // namespace

int main() {
  rising_edge_and_coda();
  aftershock_in_the_coda();
  one_token_dip_does_not_fire();
  early_restart_does_not_fire();
  restart_needs_p_above_threshold();
  release_rearms();
  warmup_consumes_without_firing();
  std::println("token trigger: rising edges and dt restarts as specified");
}
