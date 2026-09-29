// Tests for the geometry locator (pipeline/locate.hpp): the scenarios of
// onset's tests/test_locate.py, whose `onset.locate` this transcribes, with the
// same expected answers, and the network stage locating a declared event from
// StationGeometry messages alone.

#include "pipeline/locate.hpp"
#include "pipeline/network.hpp"

#include <cmath>
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

namespace {

constexpr double kEpiLat = 39.20, kEpiLon = 28.10, kOrigin = 1000.0;
const GeometryLocatorConfig kCfg{}; // vp 6, depth 10, sigma_p 0.5

struct Sta {
  const char *code;
  double lat, lon;
};
const Sta kStations[] = {{"A", 39.60, 28.00},
                         {"B", 38.90, 27.60},
                         {"C", 39.10, 28.70},
                         {"D", 39.45, 28.55}};

GeometryObs observe(const Sta &s, double sd = 0.1, double dp = 0.0,
                    double dist_scale = 1.0) {
  const double d = distance_km(s.lat, s.lon, kEpiLat, kEpiLon);
  return {s.code, s.lat, s.lon, kOrigin + std::hypot(d, 10.0) / 6.0 + dp,
          std::log(d * dist_scale), sd};
}

double off(const GeometryFit &f) {
  return distance_km(f.lat, f.lon, kEpiLat, kEpiLon);
}

void geodesy() {
  CHECK(std::abs(distance_km(39.0, 28.0, 40.0, 28.0) - 111.19) < 0.01);
}

void exact() {
  std::vector<GeometryObs> obs;
  for (const auto &s : kStations)
    obs.push_back(observe(s));
  const auto f = locate_geometry(obs, kCfg);
  CHECK(f && off(*f) < 1.0);
  CHECK(std::abs(f->origin - kOrigin) < 0.1 && f->rms < 0.05);
  CHECK(f->n_stations == 4);
}

// Two rings cross in two mirror points: two stations are not enough.
void too_few_stations() {
  CHECK(!locate_geometry({observe(kStations[0])}, kCfg));
  CHECK(!locate_geometry({observe(kStations[0]), observe(kStations[1])}, kCfg));
}

void three_stations() {
  std::vector<GeometryObs> obs;
  for (int i = 0; i < 3; ++i)
    obs.push_back(observe(kStations[i]));
  const auto f = locate_geometry(obs, kCfg);
  CHECK(f && off(*f) < 1.0);
}

void uncertainty_weights() {
  std::vector<GeometryObs> good;
  for (int i = 0; i < 3; ++i)
    good.push_back(observe(kStations[i]));
  auto unsure = good, sure = good;
  unsure.push_back(observe(kStations[3], 1.0, 0.0, 1.5));
  sure.push_back(observe(kStations[3], 0.02, 0.0, 1.5));
  const auto u = locate_geometry(unsure, kCfg), s = locate_geometry(sure, kCfg);
  CHECK(u && s);
  CHECK(off(*u) < 2.0);
  CHECK(off(*s) > off(*u) + 2.0);
  // onset.locate gives 0.35 and 15.40 km.
  CHECK(std::abs(off(*u) - 0.35) < 0.3 && std::abs(off(*s) - 15.40) < 0.5);
  std::println("  uncertain station: {:.2f} km off; overconfident: {:.2f} km",
               off(*u), off(*s));
}

// With loose P times, the rings' widths set the error radius.
void error_radius() {
  const GeometryLocatorConfig loose_p{.sigma_p = 10.0};
  std::vector<GeometryObs> tight, loose;
  for (int i = 0; i < 3; ++i) {
    tight.push_back(observe(kStations[i], 0.05));
    loose.push_back(observe(kStations[i], 0.4));
  }
  const auto t = locate_geometry(tight, loose_p);
  const auto l = locate_geometry(loose, loose_p);
  CHECK(t && l && l->err_km > 2 * t->err_km);
  // onset.locate gives 2.78 (half a coarse step, the floor) and 36.19 km.
  CHECK(std::abs(t->err_km - 2.78) < 0.05 && std::abs(l->err_km - 36.19) < 0.05);
}

// The network stage: two detections declare an event, geometry estimates
// locate it, and a station with a bad P time is left out.
void network() {
  std::map<std::string, StationInfo> st;
  for (const auto &s : kStations)
    st[s.code] = {s.lat, s.lon};
  NetworkConfig cfg;
  cfg.locator = Locator::Geometry;
  cfg.max_rms = 1.0;
  Network net(st, cfg);
  double t = kOrigin + 20;
  for (const auto &s : kStations) {
    const auto o = observe(s, 0.1, s.code[0] == 'D' ? 8.0 : 0.0);
    const double win = o.p_time - 3.5;
    net.on(Detection{s.code, win, t, 0.99f, 0.0});
    net.on(StationGeometry{s.code, win, o.p_time, 1.0, o.log_dist,
                           o.log_dist_sd, t});
    t += 0.1;
  }
  CHECK(net.declared_count() == 1);
  const auto &e = net.events().front();
  CHECK(e.location);
  CHECK(distance_km(e.location->lat, e.location->lon, kEpiLat, kEpiLon) < 1.0);
  CHECK(e.location->dropped.size() == 1 && e.location->dropped[0] == "D");
  CHECK(e.location->n_stations == 3 && !std::isnan(e.location->err_km));
  CHECK(e.first_located_at > 0 && e.first_located_at <= kOrigin + 20.2);
}

// Two stations declare the event but cannot locate it: the event keeps their
// distances as ranges and has no location until a third station reports.
void ranges_until_three() {
  std::map<std::string, StationInfo> st;
  for (const auto &s : kStations)
    st[s.code] = {s.lat, s.lon};
  NetworkConfig cfg;
  cfg.locator = Locator::Geometry;
  cfg.min_stations = 2;
  Network net(st, cfg);
  double t = kOrigin + 20;
  auto report = [&](const Sta &s) {
    const auto o = observe(s);
    const double win = o.p_time - 3.5;
    net.on(Detection{s.code, win, t, 0.99f, 0.0});
    net.on(StationGeometry{s.code, win, o.p_time, 1.0, o.log_dist,
                           o.log_dist_sd, t});
    t += 0.1;
  };
  report(kStations[0]);
  report(kStations[1]);
  CHECK(net.declared_count() == 1);
  const auto &e = net.events().front();
  CHECK(!e.location);
  CHECK(e.ranged.size() == 2 && e.ranged.contains("A") && e.ranged.contains("B"));
  report(kStations[2]);
  CHECK(e.location && e.location->n_stations == 3);
  CHECK(distance_km(e.location->lat, e.location->lon, kEpiLat, kEpiLon) < 1.0);
}

// A re-trigger within the coda window is absorbed as coda, unless it is a
// transformer dt restart: two of those declare a second event.
void restart_is_not_coda() {
  std::map<std::string, StationInfo> st;
  for (const auto &s : kStations)
    st[s.code] = {s.lat, s.lon};
  for (const bool restart : {false, true}) {
    Network net(st, NetworkConfig{});
    net.on(Detection{"A", 1000.0, 1004.0, 0.99f, 0.0});
    net.on(Detection{"B", 1001.0, 1005.0, 0.99f, 0.0});
    net.on(Detection{"A", 1020.0, 1024.0, 0.99f, 0.0, restart});
    net.on(Detection{"B", 1021.0, 1025.0, 0.99f, 0.0, restart});
    CHECK(net.declared_count() == (restart ? 2u : 1u));
    CHECK(net.events().front().coda == (restart ? 0u : 2u));
  }
}

} // namespace

int main() {
  Log::get().quiet = true;
  geodesy();
  exact();
  too_few_stations();
  three_stations();
  uncertainty_weights();
  error_radius();
  network();
  ranges_until_three();
  restart_is_not_coda();
  std::println("geometry locator agrees with onset.locate");
}
