// Tests for the geometry locator (pipeline/locate.hpp): the scenarios of
// onset's tests/test_locate.py, whose `onset.locate` this transcribes, with the
// same expected answers, and the network stage locating a declared event from
// StationGeometry messages alone.

#include "pipeline/locate.hpp"
#include "pipeline/network.hpp"

#include <cmath>
#include <cstdlib>
#include <numbers>
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
constexpr double kDeg = std::numbers::pi / 180.0;
const GeometryLocatorConfig kCfg{}; // vp 6, depth 10, sigma_p 0.5

struct Sta {
  const char *code;
  double lat, lon;
};
const Sta kStations[] = {{"A", 39.60, 28.00},
                         {"B", 38.90, 27.60},
                         {"C", 39.10, 28.70},
                         {"D", 39.45, 28.55}};

GeometryObs observe(const Sta &s, double sd = 0.1, double kappa = 20.0,
                    double dp = 0.0, double dist_scale = 1.0,
                    double baz_off_deg = 0.0) {
  const double d = distance_km(s.lat, s.lon, kEpiLat, kEpiLon);
  return {s.code,
          s.lat,
          s.lon,
          kOrigin + std::hypot(d, 10.0) / 6.0 + dp,
          std::log(d * dist_scale),
          sd,
          azimuth_rad(s.lat, s.lon, kEpiLat, kEpiLon) + baz_off_deg * kDeg,
          kappa};
}

double off(const GeometryFit &f) {
  return distance_km(f.lat, f.lon, kEpiLat, kEpiLon);
}

void geodesy() {
  double lat = 0, lon = 0;
  destination(39.0, 28.0, 60.0 * kDeg, 50.0, lat, lon);
  CHECK(std::abs(distance_km(39.0, 28.0, lat, lon) - 50.0) < 1e-6 * 50);
  CHECK(std::abs(azimuth_rad(39.0, 28.0, lat, lon) / kDeg - 60.0) < 1e-6);
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

void one_station() {
  const auto f = locate_geometry({observe(kStations[0])}, kCfg);
  CHECK(f && off(*f) < 1.0);
  auto o = observe(kStations[0]);
  o.kappa = 0;
  CHECK(!locate_geometry({o}, kCfg));
}

void distance_and_time_only() {
  std::vector<GeometryObs> obs;
  for (int i = 0; i < 3; ++i)
    obs.push_back(observe(kStations[i], 0.1, 0.0));
  const auto f = locate_geometry(obs, kCfg);
  CHECK(f && off(*f) < 2.0);
}

void uncertainty_weights() {
  std::vector<GeometryObs> good;
  for (int i = 0; i < 3; ++i)
    good.push_back(observe(kStations[i]));
  auto unsure = good, sure = good;
  unsure.push_back(observe(kStations[3], 1.0, 0.5, 0.0, 1.3, 40.0));
  sure.push_back(observe(kStations[3], 0.02, 200.0, 0.0, 1.3, 40.0));
  const auto u = locate_geometry(unsure, kCfg), s = locate_geometry(sure, kCfg);
  CHECK(u && s);
  CHECK(off(*u) < 2.0);
  CHECK(off(*s) > off(*u) + 2.0);
  // onset.locate gives 0.37 and 15.78 km.
  CHECK(std::abs(off(*u) - 0.37) < 0.3 && std::abs(off(*s) - 15.78) < 0.5);
  std::println("  uncertain station: {:.2f} km off; overconfident: {:.2f} km",
               off(*u), off(*s));
}

void error_radius() {
  const auto t = locate_geometry({observe(kStations[0], 0.05, 100.0)}, kCfg);
  const auto l = locate_geometry({observe(kStations[0], 0.4, 5.0)}, kCfg);
  CHECK(t && l && l->err_km > 2 * t->err_km);
  // onset.locate gives 4.31 and 37.56 km.
  CHECK(std::abs(t->err_km - 4.31) < 0.05 && std::abs(l->err_km - 37.56) < 0.05);
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
    const auto o = observe(s, 0.1, 20.0, s.code[0] == 'D' ? 8.0 : 0.0);
    const double win = o.p_time - 3.5;
    net.on(Detection{s.code, win, t, 0.99f, 0.0});
    net.on(StationGeometry{s.code, win, o.p_time, 1.0, o.log_dist,
                           o.log_dist_sd, o.baz, o.kappa, t});
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
  one_station();
  distance_and_time_only();
  uncertainty_weights();
  error_radius();
  network();
  restart_is_not_coda();
  std::println("geometry locator agrees with onset.locate");
}
