#include "locate.hpp"

#include "network.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace ayzek::pipeline {

namespace {
constexpr double kOneSigmaDJ = 1.15; // chi-square(2 dof, 68%) / 2

struct Eval {
  double j, origin, rms;
};

// J at (lat, lon); `parts` receives each station's share of it.
Eval evaluate(const std::vector<GeometryObs> &obs,
              const GeometryLocatorConfig &cfg, double lat, double lon,
              std::vector<double> &parts, std::vector<double> &resid) {
  const std::size_t n = obs.size();
  parts.resize(n);
  resid.resize(n);
  double sum = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const auto &o = obs[i];
    const double d = distance_km(o.lat, o.lon, lat, lon);
    const double e = std::log(std::max(d, 1.0)) - o.log_dist;
    parts[i] = e * e / (2.0 * o.log_dist_sd * o.log_dist_sd);
    resid[i] = o.p_time - std::hypot(d, cfg.depth_km) / cfg.vp;
    sum += resid[i];
  }
  const double origin = sum / static_cast<double>(n);
  double j = 0, sq = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double r = resid[i] - origin;
    parts[i] += r * r / (2.0 * cfg.sigma_p * cfg.sigma_p);
    sq += r * r;
    j += parts[i];
  }
  return {j, origin, std::sqrt(sq / static_cast<double>(n))};
}
} // namespace

std::optional<GeometryFit> locate_geometry(const std::vector<GeometryObs> &obs,
                                           const GeometryLocatorConfig &cfg,
                                           std::string *worst) {
  if (obs.size() < kMinGeometryStations)
    return std::nullopt;

  // Start: the station centroid.
  double lat0 = 0, lon0 = 0;
  for (const auto &o : obs) {
    lat0 += o.lat;
    lon0 += o.lon;
  }
  lat0 /= static_cast<double>(obs.size());
  lon0 /= static_cast<double>(obs.size());

  std::vector<double> parts, resid;
  struct Point {
    double lat, lon, j;
  };
  std::vector<Point> coarse;
  Point best{lat0, lon0, std::numeric_limits<double>::infinity()};
  auto search = [&](double clat, double clon, double half, double step,
                    std::vector<Point> *keep) {
    const auto n = static_cast<long>(std::floor(2 * half / step + 1e-9)) + 1;
    Point b = best;
    for (long a = 0; a < n; ++a)
      for (long c = 0; c < n; ++c) {
        const double la = clat + static_cast<double>(a) * step - half;
        const double lo = clon + static_cast<double>(c) * step - half;
        const double j = evaluate(obs, cfg, la, lo, parts, resid).j;
        if (keep)
          keep->push_back({la, lo, j});
        if (j < b.j)
          b = {la, lo, j};
      }
    best = b;
  };
  search(lat0, lon0, 3.0, 0.05, &coarse);
  search(best.lat, best.lon, 0.1, 0.005, nullptr);

  const auto e = evaluate(obs, cfg, best.lat, best.lon, parts, resid);
  double err = distance_km(0, 0, 0, 0.05) / 2;
  for (const auto &p : coarse)
    if (p.j <= e.j + kOneSigmaDJ)
      err = std::max(err, distance_km(best.lat, best.lon, p.lat, p.lon));
  if (worst)
    *worst = obs[static_cast<std::size_t>(
                     std::ranges::max_element(parts) - parts.begin())]
                 .station;
  double z = 0;
  for (const auto &o : obs)
    z = std::max(z, std::abs(std::log(std::max(
                                 distance_km(o.lat, o.lon, best.lat, best.lon),
                                 1.0)) -
                             o.log_dist) /
                        o.log_dist_sd);
  return GeometryFit{best.lat, best.lon, e.origin, e.rms, e.j,
                     err,      obs.size(), z};
}

} // namespace ayzek::pipeline
