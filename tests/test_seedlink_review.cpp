// Review tests for the live input (branch seedlink-review). Each shows a
// defect found in review; the test fails until the defect is fixed.
//
// 1. A component outage longer than about a ring and a half stops the station
//    for good. While a station holds more than kMaxQueuedSamples back,
//    LiveInput::on_record drops *every* record of the station, whichever
//    component it belongs to. During an outage of one component (here HHE)
//    the processor cannot advance (it scores the minimum of the three ends),
//    so Z and N fill their rings and then the queue. When HHE comes back its
//    records are dropped too, before they reach the reorderer that would
//    declare the outage a gap; without that gap the processor never advances,
//    the queue never drains, and every later record of the station is dropped
//    for the rest of the run. A replay of the same records recovers (run_ingest
//    gives up holding the source after 60 s and offers HHE).
#include "pipeline/live.hpp"
#include "records.hpp"
#include "seedlink.hpp"

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <memory>
#include <print>

using namespace ayzek;
using namespace ayzek::pipeline;

namespace {

// The start-time fields of a record `e4` units of 1e-4 s after 18:00:00.
testrec::Spec at(testrec::Spec s, std::uint64_t e4) {
  s.hour = static_cast<std::uint8_t>(18 + e4 / 36'000'000);
  s.minute = static_cast<std::uint8_t>(e4 / 600'000 % 60);
  s.second = static_cast<std::uint8_t>(e4 / 10'000 % 60);
  s.fract = static_cast<std::uint16_t>(e4 % 10'000);
  return s;
}

int component_outage() {
  auto st = std::make_unique<Station>();
  st->code = "TST";
  LiveConfig cfg;
  cfg.client.streams = seedlink::parse_streams("KO.TST");
  LiveInput live(cfg, {st.get()});

  std::array<std::int32_t, 112> samples{};
  for (std::size_t i = 0; i < samples.size(); ++i)
    samples[i] = static_cast<std::int32_t>(i);
  // Record k of a channel: 112 samples from 18:00:00 + 1.12 k s, contiguous.
  auto offer = [&](const char *cha, std::uint64_t k) {
    testrec::Spec r;
    r.cha = cha;
    live.on_record(0, testrec::int32_record(at(r, 11'200 * k), samples), 0.0);
  };
  // A processor that keeps up: it frees each ring up to the minimum of the
  // three components' ends, as Processor::run does (less its history).
  auto process = [&] {
    std::uint64_t avail = std::numeric_limits<std::uint64_t>::max();
    for (const auto &c : st->comp) {
      if (c.base.load() == kUnset)
        return;
      avail = std::min(avail, c.end());
    }
    for (auto &c : st->comp) {
      const auto base = c.base.load();
      if (avail > base)
        c.ring.set_floor(std::min<std::uint64_t>(avail - base, c.ring.written()));
    }
    live.drain();
  };

  std::uint64_t k = 0;
  for (; k < 10; ++k) { // all three components
    for (const char *cha : {"HHZ", "HHN", "HHE"})
      offer(cha, k);
    process();
  }
  const std::uint64_t e_before = st->comp[E].end();
  // HHE stops; Z and N go on until records are being dropped (about 33 min).
  while (live.stats(0).dropped < 20) {
    offer("HHZ", k);
    offer("HHN", k);
    process();
    ++k;
  }
  const std::uint64_t outage_end = k;
  const std::uint64_t dropped_in_outage = live.stats(0).dropped;
  // HHE is back. Ten more minutes of all three components.
  for (std::uint64_t i = 0; i < 540; ++i, ++k) {
    for (const char *cha : {"HHZ", "HHN", "HHE"})
      offer(cha, k);
    for (int j = 0; j < 4; ++j)
      process();
  }
  const std::uint64_t e_after = st->comp[E].end();
  const std::uint64_t dropped_after = live.stats(0).dropped - dropped_in_outage;
  std::println("  component outage: HHE end {} -> {} (outage ended at {}); "
               "{} records dropped during the outage, {} in the 10 min after "
               "it ({} offered)",
               e_before, e_after,
               st->comp[E].base.load() + outage_end * 112 - 10 * 112,
               dropped_in_outage, dropped_after, 540 * 3);
  // After the outage HHE must be written again (its outage as a gap), and the
  // station must not go on dropping everything.
  if (e_after <= e_before || dropped_after >= 540 * 3 / 2) {
    std::println(stderr,
                 "FAIL: after an HHE outage the station drops every record "
                 "for good (HHE never written again)");
    return 1;
  }
  return 0;
}

} // namespace

int main() {
  int failed = 0;
  failed += component_outage();
  if (failed)
    return 1;
  std::println("all review tests passed");
}
