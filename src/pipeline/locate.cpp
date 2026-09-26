#include "locate.hpp"

#include "network.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace ayzek::pipeline {

namespace {
constexpr double kRad = std::numbers::pi / 180.0;
constexpr double kEarthKm = 6371.0;
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
    if (o.kappa > 0)
      parts[i] +=
          o.kappa * (1.0 - std::cos(azimuth_rad(o.lat, o.lon, lat, lon) - o.baz));
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

double azimuth_rad(double lat1, double lon1, double lat2, double lon2) {
  const double p1 = lat1 * kRad, p2 = lat2 * kRad, dl = (lon2 - lon1) * kRad;
  return std::atan2(std::sin(dl) * std::cos(p2),
                    std::cos(p1) * std::sin(p2) -
                        std::sin(p1) * std::cos(p2) * std::cos(dl));
}

void destination(double lat, double lon, double az, double km, double &lat2,
                 double &lon2) {
  const double d = km / kEarthKm, p1 = lat * kRad, l1 = lon * kRad;
  const double p2 = std::asin(std::sin(p1) * std::cos(d) +
                              std::cos(p1) * std::sin(d) * std::cos(az));
  const double l2 =
      l1 + std::atan2(std::sin(az) * std::sin(d) * std::cos(p1),
                      std::cos(d) - std::sin(p1) * std::sin(p2));
  lat2 = p2 / kRad;
  lon2 = l2 / kRad;
}

std::optional<GeometryFit> locate_geometry(const std::vector<GeometryObs> &obs,
                                           const GeometryLocatorConfig &cfg,
                                           std::string *worst) {
  const bool any_baz =
      std::ranges::any_of(obs, [](const GeometryObs &o) { return o.kappa > 0; });
  if (obs.empty() || (obs.size() < 2 && !any_baz))
    return std::nullopt;

  // Start: the mean of the single-station epicentres, or the centroid.
  double lat0 = 0, lon0 = 0;
  std::size_t n0 = 0;
  for (const auto &o : obs) {
    if (any_baz && !(o.kappa > 0))
      continue;
    double la = o.lat, lo = o.lon;
    if (any_baz)
      destination(o.lat, o.lon, o.baz, std::exp(o.log_dist), la, lo);
    lat0 += la;
    lon0 += lo;
    ++n0;
  }
  lat0 /= static_cast<double>(n0);
  lon0 /= static_cast<double>(n0);

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
  return GeometryFit{best.lat, best.lon, e.origin, e.rms,
                     e.j,      err,      obs.size()};
}

} // namespace ayzek::pipeline
