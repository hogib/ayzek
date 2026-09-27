// Tests for the alarm assessment (pipeline/assess.hpp): picks from one source
// are an earthquake, picks no source explains are a misfire, too few picks
// are unclassified, and one bad station among good ones is left out.

#include "models.hpp"
#include "pipeline/network.hpp"

#include <cmath>
#include <cstdlib>
#include <print>

#define CHECK(cond)                                                            \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::println(stderr, "CHECK failed: {}  ({}:{})", #cond, __FILE__,       \
                   __LINE__);                                                  \
      std::abort();                                                            \
    }                                                                          \
  } while (0)

using namespace ayzek::pipeline;

namespace {

constexpr double kLat = 39.20, kLon = 28.10, kOrigin = 1000.0;
struct Sta {
  const char *code;
  double lat, lon;
};
const Sta kStations[] = {{"A", 39.60, 28.00},
                         {"B", 38.90, 27.60},
                         {"C", 39.10, 28.70},
                         {"D", 39.45, 28.55}};

std::map<std::string, StationInfo> stations() {
  std::map<std::string, StationInfo> m;
  for (const auto &s : kStations)
    m[s.code] = {s.lat, s.lon};
  return m;
}

double travel(const Sta &s, double v) {
  return std::hypot(distance_km(s.lat, s.lon, kLat, kLon), 10.0) / v;
}

// Declares an event from all four stations; `picked` stations get a pick with
// S shifted by `s_shift` seconds (0: the true S). Each station triggered once
// in `span` seconds, which sets the chance-coincidence rate.
Assessment run(const std::vector<std::pair<int, double>> &picked,
               bool magnitudes = false, double prob = 0.9, double span = 1200) {
  NetworkConfig cfg;
  cfg.assess = true;
  Network net(stations(), cfg);
  for (const auto &s : kStations) {
    const double p = kOrigin + travel(s, 6.0);
    net.on(Detection{s.code, p - 3.5, p + 1.0, 0.99f, 0.0});
  }
  for (const auto &[i, shift] : picked) {
    const auto &s = kStations[i];
    const double p = kOrigin + travel(s, 6.0), sv = kOrigin + travel(s, 3.5) + shift;
    net.on(Pick{s.code, p - 3.5, p, sv, prob, prob, p + 56.5, 1.0, nullptr});
  }
  if (magnitudes)
    for (const auto &s : kStations) {
      const double p = kOrigin + travel(s, 6.0);
      net.on(MagnitudeEstimate{s.code, p - 3.5, false, p - 2, 3.0f + 0.1f * (s.code[0] - 'A'),
                               0, p + 8, 1.0});
    }
  CHECK(net.declared_count() == 1);
  return net.assess(net.events().front(), kOrigin - span / 2, kOrigin + span / 2);
}

void four_consistent_stations_are_an_earthquake() {
  const auto a = run({{0, 0}, {1, 0}, {2, 0}, {3, 0}});
  CHECK(a.verdict == Verdict::Earthquake && a.consistent);
  CHECK(distance_km(a.lat, a.lon, kLat, kLon) < 2.0 && std::abs(a.origin - kOrigin) < 0.3);
  CHECK(a.stations.size() == 4 && a.rms < 0.1);
}

void two_stations_need_support() {
  // Stations triggering every few seconds: a coincidence is no evidence.
  CHECK(run({{0, 0}, {2, 0}}, false, 0.9, 10).verdict == Verdict::Possible);
  // Rare triggers: it is.
  CHECK(run({{0, 0}, {2, 0}}, false, 0.9, 1e6).verdict == Verdict::Earthquake);
  // Station magnitudes that agree support it too, with nothing against.
  CHECK(run({{0, 0}, {2, 0}}, true, 0.9, 400).verdict == Verdict::Earthquake);
}

void s_p_that_no_source_explains_is_a_misfire() {
  // C's S 14 s late: ~120 km too far. Two stations cannot say which pick is
  // wrong, so it takes an objection besides: frequent triggers.
  CHECK(run({{0, 0}, {2, 14.0}}, false, 0.9, 10).verdict == Verdict::Misfire);
  CHECK(run({{0, 0}, {2, 14.0}}, false, 0.9, 1e6).verdict == Verdict::Unclassified);
  // Three stations that no single source reconciles, even leaving one out.
  CHECK(run({{0, 0}, {1, 14.0}, {2, 10.0}}, false, 0.9, 1e6).verdict ==
        Verdict::Misfire);
}

void an_implausible_s_p_is_not_used() {
  const auto a = run({{0, 0}, {1, 0}, {2, 0}, {3, 40.0}}); // D's S-P ~50 s
  CHECK(a.verdict == Verdict::Earthquake && a.stations.size() == 3);
}

void one_bad_station_among_good_ones_is_left_out() {
  const auto a = run({{0, 0}, {1, 0}, {2, 0}, {3, 7.0}});
  CHECK(a.verdict == Verdict::Earthquake);
  CHECK(!a.stations[3].used && a.stations[0].used);
}

void too_few_picks_are_unclassified() {
  CHECK(run({{0, 0}}).verdict == Verdict::Unclassified);
  const auto a = run({{0, 0}, {1, 0}}, false, 0.3); // both fail the pick QC
  CHECK(a.verdict == Verdict::Unclassified && a.picks_rejected == 2);
}

// A station that did not trigger is asked for a pick at its predicted P, and
// that pick joins the assessment.
void pick_at_a_station_that_did_not_trigger() {
  NetworkConfig cfg;
  cfg.assess = true;
  Network net(stations(), cfg);
  for (int i = 0; i < 3; ++i) {
    const auto &s = kStations[i];
    const double p = kOrigin + travel(s, 6.0);
    net.on(Detection{s.code, p - 3.5, p + 1.0, 0.99f, 0.0});
    net.on(Pick{s.code, p - 3.5, p, kOrigin + travel(s, 3.5), 0.9, 0.9,
                p + 56.5, 1.0, nullptr});
  }
  const auto req = net.pick_requests(kOrigin - 600, kOrigin + 600);
  CHECK(req.size() == 1 && req[0].station == "D");
  const auto &d = kStations[3];
  CHECK(std::abs(req[0].p_pred - (kOrigin + travel(d, 6.0))) < 0.5);
  net.on(Pick{"D", req[0].p_pred - 3.5, kOrigin + travel(d, 6.0),
              kOrigin + travel(d, 3.5), 0.8, 0.8, req[0].p_pred + 54.5, 1.0,
              nullptr, req[0].event_id, req[0].event_alarm});
  const auto a = net.assess(net.events().front(), kOrigin - 600, kOrigin + 600);
  CHECK(a.stations.size() == 4 && a.stations[3].requested);
  CHECK(a.verdict == Verdict::Earthquake);
}

// The picker's search: P only near the expected P, S only after it and before
// the next onset.
void picker_search_window() {
  using ayzek::Picker;
  std::vector<float> lg(Picker::kChunks * Picker::kClasses, 0.0f);
  auto set = [&](double seconds, int cls, float v) {
    lg[static_cast<std::size_t>(seconds / 0.24) * Picker::kClasses + cls] = v;
  };
  set(5.5, 1, 4.0f);  // the event's P
  set(40.0, 1, 6.0f); // a later event's, stronger
  set(12.0, 2, 3.0f); // the event's S
  set(47.0, 2, 5.0f); // the later event's
  set(3.0, 2, 6.0f);  // "S" before P
  const auto free = Picker::select(lg, nullptr);
  CHECK(std::abs(free.p_seconds - 40.0) < 0.3 && std::abs(free.s_seconds - 3.0) < 0.3);
  const Picker::Search search{3.0, 8.0, 36.0, 0.5};
  const auto on = Picker::select(lg, &search);
  CHECK(std::abs(on.p_seconds - 5.5) < 0.3 && std::abs(on.s_seconds - 12.0) < 0.3);
}

} // namespace

int main() {
  Log::get().quiet = true;
  four_consistent_stations_are_an_earthquake();
  two_stations_need_support();
  s_p_that_no_source_explains_is_a_misfire();
  one_bad_station_among_good_ones_is_left_out();
  an_implausible_s_p_is_not_used();
  too_few_picks_are_unclassified();
  pick_at_a_station_that_did_not_trigger();
  picker_search_window();
  std::println("alarm assessment: verdicts as specified");
}
