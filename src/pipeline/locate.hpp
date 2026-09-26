#pragma once

// Locating an event from the transformer's geometry head, in place of the
// S-P picks. A transcription of onset's `onset/locate.py`; tests/test_locate.cpp
// holds it to the same scenarios as onset's tests/test_locate.py.
//
// Each station contributes a likelihood over the epicentre: its distance
// (Gaussian in log km) and back-azimuth (von Mises), exactly the terms the head
// was trained on, and its P time (dated by the detector from dt). The location
// minimises
//
//   J(x) = sum_i (log max(D_i(x), 1) - log_dist_i)^2 / (2 sd_i^2)
//        + sum_i kappa_i (1 - cos(AZ_i(x) - baz_i))
//        + sum_i (r_i(x) - origin(x))^2 / (2 sigma_p^2)
//
// where r_i is the P time minus the travel time from x in a uniform
// half-space and origin(x) their mean. Grid search at fixed depth: +-3 deg at
// 0.05 deg around the mean of the single-station epicentres, then +-0.1 deg at
// 0.005 deg. One station with a back-azimuth is enough; without one, two.

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace ayzek::pipeline {

struct GeometryObs {
  std::string station;
  double lat = 0, lon = 0;
  double p_time = 0;      // epoch
  double log_dist = 0;    // log km
  double log_dist_sd = 1; // standard deviation of log_dist
  double baz = 0;         // station -> event, radians clockwise from north
  double kappa = 0;       // 0: no back-azimuth
};

struct GeometryFit {
  double lat, lon, origin;
  double rms;    // P-time residual rms, seconds
  double misfit; // J at the solution
  double err_km; // radius of the 68% region (J <= min + 1.15) on the coarse grid
  std::size_t n_stations;
  // The largest |log D_i - log_dist_i| / sd_i at the solution: how far the
  // station that disagrees most is from its own distance, in its own sds.
  double dist_z_max;
};

struct GeometryLocatorConfig {
  double vp = 6.0;       // km/s
  double depth_km = 10.0;
  double sigma_p = 0.5;  // P-time uncertainty, seconds
};

double azimuth_rad(double lat1, double lon1, double lat2, double lon2);
// The point `km` from (lat, lon) along bearing `az`, radians.
void destination(double lat, double lon, double az, double km, double &lat2,
                 double &lon2);

// nullopt if the stations cannot fix an epicentre. With `worst`, the station
// contributing most to J at the solution.
std::optional<GeometryFit> locate_geometry(const std::vector<GeometryObs> &obs,
                                           const GeometryLocatorConfig &cfg,
                                           std::string *worst = nullptr);

} // namespace ayzek::pipeline
