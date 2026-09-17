// Evaluates station subsets from one ayzek recording.
//
//   network_subsets RECORDING.tsv --catalog CSV [options] > subsets.csv
//
// For every subset of the recorded stations with a size in [min-k, max-k], the
// network stage is run on that subset's messages only, and the result is
// compared with the catalogue. Because station processing does not depend on
// the other stations, this equals running ayzek on each subset.
//
// Output: one CSV row per subset. The main event is the largest catalogue
// event within the evaluation radius; its columns give the alert delay, the
// S-wave blind-zone radius at the alert, the first and final magnitude, and
// the location error. A summary of the best subsets per size and of the
// marginal gain from adding each station is written to stderr.
//
// Options:
//   --catalog CSV       AFAD catalogue (required)
//   --min-k N           smallest subset size (default 2)
//   --max-k N           largest subset size (default: all stations)
//   --min-mag M         only catalogue events with magnitude >= M (default 0)
//   --radius KM         evaluation radius around the full network's centroid
//   (default 250)
//   --min-stations N    detections required to declare an event (default 2)

#include "pipeline/network.hpp"
#include "pipeline/recording.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <map>
#include <print>
#include <string>
#include <vector>

using namespace ayzek::pipeline;

namespace {

constexpr double kVs = 3.5; // km/s, for the blind-zone radius

struct Row {
  std::vector<std::string> stations;
  std::size_t catalogue = 0, declared = 0, located = 0, unmatched = 0;
  double mean_alert = NAN, median_location_error = NAN,
         mean_abs_magnitude_error = NAN;
  // main event
  bool main_declared = false;
  double main_alert = NAN, main_blind_km = NAN, main_first_mag = NAN,
         main_first_mag_delay = NAN, main_final_mag = NAN,
         main_location_error = NAN;
};

double median(std::vector<double> v) {
  if (v.empty())
    return NAN;
  std::ranges::sort(v);
  return v.size() % 2 ? v[v.size() / 2]
                      : 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]);
}

double mean(const std::vector<double> &v) {
  if (v.empty())
    return NAN;
  double s = 0;
  for (double x : v)
    s += x;
  return s / static_cast<double>(v.size());
}

std::string fmt(double x, int prec = 2) {
  return std::isnan(x) ? "" : std::format("{:.{}f}", x, prec);
}

} // namespace

int main(int argc, char **argv) try {
  if (argc < 2) {
    std::println(stderr, "usage: network_subsets RECORDING.tsv --catalog CSV "
                         "[--min-k N] [--max-k N] [--min-mag M] "
                         "[--radius KM] [--min-stations N]");
    return 2;
  }
  const auto rec = Recording::load(argv[1]);
  std::string catalog_path;
  std::size_t min_k = 2, max_k = rec.stations.size();
  double min_mag = 0;
  NetworkConfig base;
  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&] {
      if (i + 1 >= argc)
        throw std::runtime_error(a + " needs a value");
      return std::string(argv[++i]);
    };
    if (a == "--catalog")
      catalog_path = next();
    else if (a == "--min-k")
      min_k = std::stoul(next());
    else if (a == "--max-k")
      max_k = std::stoul(next());
    else if (a == "--min-mag")
      min_mag = std::stod(next());
    else if (a == "--radius")
      base.catalog_radius_km = std::stod(next());
    else if (a == "--min-stations")
      base.min_stations = std::stoul(next());
    else
      throw std::runtime_error("unknown option " + a);
  }
  if (catalog_path.empty())
    throw std::runtime_error("--catalog is required");

  // Catalogue events within the radius of the full network's centroid, so
  // that every subset is scored on the same events.
  double clat = 0, clon = 0;
  for (const auto &[_, s] : rec.stations) {
    clat += s.lat;
    clon += s.lon;
  }
  clat /= static_cast<double>(rec.stations.size());
  clon /= static_cast<double>(rec.stations.size());
  for (auto &c :
       load_afad_catalog(catalog_path, rec.t_first - 30, rec.t_last - 30))
    if (c.magnitude >= min_mag &&
        distance_km(clat, clon, c.lat, c.lon) <= base.catalog_radius_km)
      base.catalog.push_back(c);
  base.catalog_radius_km = 1e9; // selection done above
  if (base.catalog.empty())
    throw std::runtime_error(
        "no catalogue events in the recording's span and radius");
  const auto main_event =
      std::ranges::max_element(base.catalog, {}, &CatalogEvent::magnitude);

  std::vector<std::string> codes;
  for (const auto &[code, _] : rec.stations)
    codes.push_back(code);
  const std::size_t n = codes.size();
  if (n > 20)
    throw std::runtime_error("more than 20 stations: too many subsets");
  std::println(stderr,
               "{} stations, {} catalogue events (main: {}{:.1f} at {}), "
               "subset sizes {}..{}",
               n, base.catalog.size(), main_event->type, main_event->magnitude,
               hms(main_event->time), min_k, max_k);

  Log::get().quiet = true;
  std::vector<Row> rows;
  for (std::uint32_t mask = 1; mask < (std::uint32_t{1} << n); ++mask) {
    const auto k = static_cast<std::size_t>(std::popcount(mask));
    if (k < min_k || k > max_k)
      continue;
    Row row;
    std::map<std::string, StationInfo> subset;
    for (std::size_t i = 0; i < n; ++i)
      if (mask & (std::uint32_t{1} << i)) {
        subset[codes[i]] = rec.stations.at(codes[i]);
        row.stations.push_back(codes[i]);
      }

    Network net(subset, base);
    for (const auto &m : rec.messages) {
      if (!subset.contains(station_of(m)))
        continue;
      if (auto *d = std::get_if<Detection>(&m))
        net.on(*d);
      else if (auto *p = std::get_if<Pick>(&m))
        net.on(*p);
      else if (auto *g = std::get_if<MagnitudeEstimate>(&m))
        net.on(*g);
    }

    std::vector<double> alerts, loc_errors, mag_errors;
    for (const auto &sc : net.score()) {
      ++row.catalogue;
      const auto &c = *sc.event;
      const Event *e = sc.detected;
      if (!e)
        continue;
      ++row.declared;
      const double alert = e->declared_at - c.time;
      alerts.push_back(alert);
      double loc_error = NAN;
      if (e->location) {
        ++row.located;
        loc_error =
            distance_km(e->location->lat, e->location->lon, c.lat, c.lon);
        loc_errors.push_back(loc_error);
      }
      if (e->magnitude)
        mag_errors.push_back(std::abs(*e->magnitude - c.magnitude));
      // The network holds its own copy of the catalogue; compare by value.
      if (c.time == main_event->time && c.lat == main_event->lat &&
          c.lon == main_event->lon) {
        row.main_declared = true;
        row.main_alert = alert;
        row.main_blind_km = std::sqrt(
            std::max(kVs * kVs * alert * alert - c.depth * c.depth, 0.0));
        if (e->first_magnitude) {
          row.main_first_mag = *e->first_magnitude;
          row.main_first_mag_delay = e->first_magnitude_at - c.time;
        }
        if (e->magnitude)
          row.main_final_mag = *e->magnitude;
        row.main_location_error = loc_error;
      }
    }
    row.unmatched = net.unmatched_count();
    row.mean_alert = mean(alerts);
    row.median_location_error = median(loc_errors);
    row.mean_abs_magnitude_error = mean(mag_errors);
    rows.push_back(std::move(row));
  }

  std::println(
      "k,stations,catalogue,declared,located,unmatched,mean_alert_s,median_"
      "location_error_km,"
      "mean_abs_magnitude_error,main_declared,main_alert_s,main_blind_zone_km,"
      "main_first_magnitude,"
      "main_first_magnitude_s,main_final_magnitude,main_location_error_km");
  for (const auto &r : rows) {
    std::string names;
    for (const auto &s : r.stations)
      names += (names.empty() ? "" : "+") + s;
    std::println("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                 r.stations.size(), names, r.catalogue, r.declared, r.located,
                 r.unmatched, fmt(r.mean_alert), fmt(r.median_location_error),
                 fmt(r.mean_abs_magnitude_error), r.main_declared ? 1 : 0,
                 fmt(r.main_alert), fmt(r.main_blind_km, 1),
                 fmt(r.main_first_mag), fmt(r.main_first_mag_delay),
                 fmt(r.main_final_mag), fmt(r.main_location_error));
  }

  // stderr summary: per size, the subset with the earliest main-event alert
  // (ties broken by more catalogue events declared, then fewer unmatched).
  std::println(stderr,
               "\nbest subset per size (earliest alert for the main event):");
  std::println(stderr, "  {:>2}  {:<40} {:>8} {:>9} {:>7} {:>9} {:>9} {:>8}",
               "k", "stations", "alert s", "blind km", "first M", "final M",
               "declared", "extra");
  for (std::size_t k = min_k; k <= max_k; ++k) {
    const Row *best = nullptr;
    std::size_t with_alert = 0, total = 0;
    for (const auto &r : rows) {
      if (r.stations.size() != k)
        continue;
      ++total;
      if (!r.main_declared)
        continue;
      ++with_alert;
      auto key = [](const Row *x) {
        return std::tuple(x->main_alert, -static_cast<double>(x->declared),
                          x->unmatched);
      };
      if (!best || key(&r) < key(best))
        best = &r;
    }
    if (!best) {
      std::println(stderr,
                   "  {:>2}  no subset declared the main event ({} subsets)", k,
                   total);
      continue;
    }
    std::string names;
    for (const auto &s : best->stations)
      names += (names.empty() ? "" : "+") + s;
    std::println(stderr,
                 "  {:>2}  {:<40} {:>8} {:>9} {:>7} {:>9} {:>6}/{:<2} {:>8}   "
                 "({} of {} subsets declare it)",
                 k, names, fmt(best->main_alert, 1),
                 fmt(best->main_blind_km, 0), fmt(best->main_first_mag, 1),
                 fmt(best->main_final_mag, 1), best->declared, best->catalogue,
                 best->unmatched, with_alert, total);
  }

  // Marginal gain of each station: for every evaluated subset S without the
  // station whose extension S + station was also evaluated, the change from S
  // to S + station, averaged over all such S.
  std::map<std::uint32_t, const Row *> by_mask;
  std::vector<std::uint32_t> masks;
  {
    std::size_t r = 0;
    for (std::uint32_t mask = 1; mask < (std::uint32_t{1} << n); ++mask) {
      const auto k = static_cast<std::size_t>(std::popcount(mask));
      if (k < min_k || k > max_k)
        continue;
      by_mask[mask] = &rows[r++];
    }
  }
  std::println(stderr,
               "\nmarginal gain from adding each station to a subset without "
               "it (mean over {}..{}-station subsets):",
               min_k, max_k - 1);
  std::println(stderr, "  {:<6} {:>9} {:>12} {:>14} {:>16} {:>18}", "",
               "subsets", "+declared", "+not in cat.", "main alert s",
               "enables main alert");
  for (std::size_t i = 0; i < n; ++i) {
    const std::uint32_t bit = std::uint32_t{1} << i;
    double d_declared = 0, d_unmatched = 0, d_alert = 0;
    std::size_t count = 0, both = 0, enables = 0, without_alert = 0;
    for (const auto &[mask, row] : by_mask) {
      if (mask & bit)
        continue;
      auto it = by_mask.find(mask | bit);
      if (it == by_mask.end())
        continue;
      const Row &a = *row;
      const Row &b = *it->second;
      ++count;
      d_declared +=
          static_cast<double>(b.declared) - static_cast<double>(a.declared);
      d_unmatched +=
          static_cast<double>(b.unmatched) - static_cast<double>(a.unmatched);
      if (a.main_declared && b.main_declared) {
        d_alert += b.main_alert - a.main_alert;
        ++both;
      }
      if (!a.main_declared) {
        ++without_alert;
        enables += b.main_declared;
      }
    }
    if (!count)
      continue;
    std::println(stderr, "  {:<6} {:>9} {:>12} {:>14} {:>16} {:>18}", codes[i],
                 count, fmt(d_declared / count, 2), fmt(d_unmatched / count, 2),
                 both ? fmt(d_alert / both, 2) : "",
                 without_alert ? std::format("{} of {}", enables, without_alert)
                               : "");
  }
  return 0;
} catch (const std::exception &e) {
  std::println(stderr, "network_subsets: {}", e.what());
  return 1;
}
