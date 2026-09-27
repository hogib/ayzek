// Network::assess and its output (assess.hpp, --assess).

#include "network.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <ranges>

namespace ayzek::pipeline {

namespace {

double median(std::vector<double> v) {
  if (v.empty())
    return NAN;
  std::ranges::sort(v);
  return v.size() % 2 ? v[v.size() / 2]
                      : 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]);
}

std::string join(const std::vector<std::string> &v, const char *sep) {
  std::string out;
  for (const auto &s : v)
    out += (out.empty() ? "" : sep) + s;
  return out;
}

// Replaces the CSV separator and quotes in free text.
std::string csv_text(std::string s) {
  std::ranges::replace(s, ',', ';');
  std::ranges::replace(s, '"', '\'');
  return s;
}

} // namespace

Assessment Network::assess(const Event &e, double t0, double t1) const {
  Assessment a;
  a.picks_rejected = e.picks_rejected;
  const double km_per_s = cfg_.vp * cfg_.vs / (cfg_.vp - cfg_.vs);
  std::vector<std::string> support, against;

  // Stations with confident P and S picks (Network::on(Pick) keeps only
  // those) and known coordinates.
  // An S-P beyond `assess_max_sp_km` is a pick on the wrong phase (in an
  // aftershock sequence the 60 s picker window often holds another event),
  // not evidence about this one.
  std::map<std::string, Pick> picks;
  std::string implausible;
  for (const auto &[code, pk] : e.picks) {
    if (!stations_.contains(code))
      continue;
    const double km = (pk.s_time - pk.p_time) * km_per_s;
    if (km > cfg_.assess_max_sp_km)
      implausible += std::format("{}{} {:.0f} s", implausible.empty() ? "" : ", ",
                                 code, pk.s_time - pk.p_time);
    else
      picks.emplace(code, pk);
  }
  for (const auto &[code, pk] : picks) {
    StationSP s{code, pk.s_time - pk.p_time, (pk.s_time - pk.p_time) * km_per_s};
    s.requested = pk.event_id != 0;
    if (auto g = e.geometry.find(code); g != e.geometry.end()) {
      s.geo_km = std::exp(g->second.log_dist);
      s.geo_z = (std::log(std::max(s.km, 1.0)) - g->second.log_dist) /
                (g->second.log_dist_sd * cfg_.geo_sd_scale);
    }
    a.stations.push_back(s);
  }
  auto list = [&] {
    std::string out;
    for (const auto &s : a.stations)
      if (s.used)
        out += std::format("{}{} {:.1f} s ({:.0f} km{})", out.empty() ? "" : ", ",
                           s.station, s.sp, s.km,
                           s.requested ? ", not triggered" : "");
    return out;
  };

  // The S-P test: one source for every station's P and S. Two stations'
  // circles must be able to meet, the grid-search fit must be tight overall,
  // and no station may be far off it.
  std::string why_not;
  auto test = [&](const std::map<std::string, Pick> &pk,
                  std::string *worst_station) -> std::optional<Location> {
    for (auto i = pk.begin(); i != pk.end(); ++i)
      for (auto j = std::next(i); j != pk.end(); ++j) {
        const auto &si = stations_.at(i->first), &sj = stations_.at(j->first);
        const double sep = distance_km(si.lat, si.lon, sj.lat, sj.lon);
        const double di = (i->second.s_time - i->second.p_time) * km_per_s;
        const double dj = (j->second.s_time - j->second.p_time) * km_per_s;
        if (std::abs(di - dj) > sep + 5.0 || sep > di + dj + 5.0) {
          why_not = std::format("{} and {} are {:.0f} km apart, but S-P puts "
                                "the source {:.0f} and {:.0f} km from them",
                                i->first, j->first, sep, di, dj);
          *worst_station = di > dj ? i->first : j->first;
          return std::nullopt;
        }
      }
    auto fit = locate(pk, nullptr);
    if (!fit) {
      why_not = "no fit";
      return std::nullopt;
    }
    double largest = 0;
    for (const auto &[code, p] : pk) {
      const auto &st = stations_.at(code);
      const double d = std::hypot(distance_km(fit->lat, fit->lon, st.lat, st.lon),
                                  cfg_.depth_km);
      const double r = std::max(std::abs(p.p_time - fit->origin - d / cfg_.vp),
                                std::abs(p.s_time - fit->origin - d / cfg_.vs));
      if (r > largest) {
        largest = r;
        *worst_station = code;
      }
    }
    if (fit->rms > cfg_.assess_max_rms) {
      why_not = std::format("the best single source misses the P and S times "
                            "by {:.1f} s rms", fit->rms);
      return std::nullopt;
    }
    if (largest > cfg_.assess_station_tol) {
      why_not = std::format("{} is {:.1f} s off the best single source",
                            *worst_station, largest);
      return std::nullopt;
    }
    return fit;
  };

  if (picks.empty()) {
    a.reasons.push_back(
        a.picks_rejected
            ? std::format("no confident P and S at the {} station{} picked",
                          a.picks_rejected, a.picks_rejected == 1 ? "" : "s")
            : std::string("no picks"));
  } else if (picks.size() == 1) {
    const auto &s = a.stations.front();
    a.reasons.push_back(std::format(
        "one station with P and S: {} {:.1f} s ({:.0f} km){}", s.station, s.sp,
        s.km,
        std::isnan(s.geo_z)
            ? std::string()
            : std::format(", the transformer says {:.0f} km", s.geo_km)));
  } else {
    std::string worst;
    auto fit = test(picks, &worst);
    std::string dropped;
    if (!fit && picks.size() >= 3) {
      // One bad pick should not condemn an event the others agree on.
      const std::string first_why = why_not;
      auto rest = picks;
      rest.erase(worst);
      std::string w2;
      fit = test(rest, &w2);
      if (fit) {
        dropped = worst;
        for (auto &s : a.stations)
          s.used = s.station != worst;
      } else {
        why_not = first_why;
      }
    }
    a.consistent = fit.has_value();
    if (!fit) {
      a.reasons.push_back("S-P does not fit one source: " + why_not);
    } else {
      a.lat = fit->lat;
      a.lon = fit->lon;
      a.origin = fit->origin;
      a.rms = fit->rms;
      for (auto &s : a.stations) {
        const auto &st = stations_.at(s.station);
        const auto &p = picks.at(s.station);
        const double d = std::hypot(distance_km(a.lat, a.lon, st.lat, st.lon),
                                    cfg_.depth_km);
        s.p_resid = p.p_time - a.origin - d / cfg_.vp;
        s.s_resid = p.s_time - a.origin - d / cfg_.vs;
      }
      a.reasons.push_back(std::format(
          "S-P fits one source at {} stations, rms {:.2f} s: {}",
          picks.size() - !dropped.empty(), a.rms, list()));
      if (!dropped.empty())
        a.reasons.push_back(dropped + "'s picks do not fit the others (left out)");
    }
  }

  // Geometry: the transformer's distances against the picker's.
  std::vector<double> z;
  double far = 0;
  for (const auto &s : a.stations)
    if (s.used && !std::isnan(s.geo_z)) {
      z.push_back(std::abs(s.geo_z));
      far = std::max(far, s.km);
    }
  if (!z.empty()) {
    a.geo_z_median = median(z);
    if (a.consistent && e.location && cfg_.locator == Locator::Geometry)
      a.geo_offset_km = distance_km(a.lat, a.lon, e.location->lat, e.location->lon);
    const double near_km =
        std::max(30.0, e.location && !std::isnan(e.location->err_km)
                           ? 2 * e.location->err_km
                           : 0.0);
    const bool agree =
        a.geo_z_median <= cfg_.assess_geo_z &&
        (std::isnan(a.geo_offset_km) || a.geo_offset_km <= near_km);
    const std::string where =
        std::isnan(a.geo_offset_km)
            ? std::string()
            : std::format(", epicentres {:.0f} km apart", a.geo_offset_km);
    if (agree)
      support.push_back(std::format("the transformer's distances agree "
                                    "(median {:.1f} sd{})",
                                    a.geo_z_median, where));
    else if (far > 60.0)
      a.reasons.push_back(std::format(
          "the transformer's distances differ (median {:.1f} sd{}), but S-P "
          "puts a station {:.0f} km out, beyond the distances it was trained on",
          a.geo_z_median, where, far));
    else if (a.geo_z_median > 2 * cfg_.assess_geo_z)
      against.push_back(std::format("the transformer's distances disagree "
                                    "(median {:.1f} sd{})",
                                    a.geo_z_median, where));
  }

  // Station magnitudes.
  if (e.magnitudes.size() >= 2) {
    std::vector<double> m;
    for (const auto &[_, sm] : e.magnitudes)
      m.push_back(sm.magnitude);
    const double mean = std::accumulate(m.begin(), m.end(), 0.0) / m.size();
    double sq = 0;
    for (double x : m)
      sq += (x - mean) * (x - mean);
    a.mag_spread = std::sqrt(sq / m.size());
    const auto [lo, hi] = std::ranges::minmax(m);
    if (a.mag_spread <= cfg_.assess_mag_spread)
      support.push_back(std::format("station magnitudes agree (M{:.1f}-{:.1f})",
                                    lo, hi));
    else if (a.mag_spread > 2 * cfg_.assess_mag_spread)
      against.push_back(std::format(
          "station magnitudes disagree (M{:.1f}-{:.1f})", lo, hi));
  }

  // Chance coincidence: the other declaring stations each triggering inside
  // their association window with the first, at their own trigger rates.
  std::vector<const Detection *> declaring;
  for (const auto &[_, d] : e.detections)
    if (d.declared_at <= e.declared_at && stations_.contains(d.station))
      declaring.push_back(&d);
  const double span = t1 - t0;
  if (declaring.size() >= 2 && span > 0) {
    std::ranges::sort(declaring, {}, &Detection::window_start);
    const auto &first = stations_.at(declaring.front()->station);
    a.p_chance = 1.0;
    for (std::size_t i = 1; i < declaring.size(); ++i) {
      const auto &st = stations_.at(declaring[i]->station);
      const double w =
          2 * (distance_km(first.lat, first.lon, st.lat, st.lon) / cfg_.vp +
               cfg_.slack_seconds);
      const auto n = triggers_.contains(declaring[i]->station)
                         ? triggers_.at(declaring[i]->station)
                         : 0;
      a.p_chance *= 1.0 - std::exp(-static_cast<double>(n) / span * w);
    }
    if (a.p_chance <= cfg_.assess_p_chance)
      support.push_back(
          std::format("unlikely to be a chance coincidence (p {:.1e})", a.p_chance));
    else if (a.p_chance > 0.1)
      against.push_back(std::format(
          "a chance coincidence is plausible (p {:.2f})", a.p_chance));
  }

  const std::size_t used = static_cast<std::size_t>(
      std::ranges::count_if(a.stations, &StationSP::used));
  if (!implausible.empty())
    a.reasons.push_back(std::format("S-P too long to be this event, not used: {}",
                                    implausible));
  // Two stations that disagree cannot say which pick is wrong: that is a
  // misfire only with an objection besides. Three or more that cannot be
  // reconciled, even leaving one out, are one.
  if (picks.size() < 2)
    a.verdict = Verdict::Unclassified;
  else if (!a.consistent)
    a.verdict = picks.size() >= 3 || !against.empty() ? Verdict::Misfire
                                                      : Verdict::Unclassified;
  else if (used >= 3 || (!support.empty() && against.empty()))
    a.verdict = Verdict::Earthquake;
  else
    a.verdict = Verdict::Possible;
  for (const auto &s : support)
    a.reasons.push_back("+ " + s);
  for (const auto &s : against)
    a.reasons.push_back("- " + s);

  // Catalogue context.
  const auto &ctx = cfg_.context_catalog.empty() ? cfg_.catalog : cfg_.context_catalog;
  if (const auto *c = match(e)) {
    a.context = std::format("catalogued: AFAD {} {:.1f} at {}", c->type,
                            c->magnitude, hms(c->time));
  } else {
    double lat = a.lat, lon = a.lon, origin = a.origin;
    if (std::isnan(lat) && e.location) {
      lat = e.location->lat;
      lon = e.location->lon;
      origin = e.location->origin;
    }
    if (!std::isnan(lat)) {
      const CatalogEvent *near = nullptr, *zone = nullptr;
      for (const auto &c : ctx) {
        const double km = distance_km(lat, lon, c.lat, c.lon);
        if (std::abs(c.time - origin) <= 60.0 && km <= 100.0 &&
            (!near || std::abs(c.time - origin) < std::abs(near->time - origin)))
          near = &c;
        if (c.time < origin - 60.0 && origin - c.time <= 30 * 86400.0 &&
            km <= 50.0 && c.magnitude >= 3.0 &&
            (!zone || c.magnitude > zone->magnitude))
          zone = &c;
      }
      if (near)
        a.context = std::format(
            "AFAD {} {:.1f} at {}, {:.0f} km away, origin {:+.1f} s off: "
            "probably this event, missed by the catalogue matching",
            near->type, near->magnitude, hms(near->time),
            distance_km(lat, lon, near->lat, near->lon), origin - near->time);
      else if (zone) {
        const double h = (origin - zone->time) / 3600.0;
        a.context = std::format(
            "in the aftershock zone of AFAD {} {:.1f} ({} UTC), {:.0f} km "
            "away, {}",
            zone->type, zone->magnitude, ymd_hms(zone->time),
            distance_km(lat, lon, zone->lat, zone->lon),
            h < 48 ? std::format("{:.1f} h earlier", h)
                   : std::format("{:.0f} days earlier", h / 24));
      }
    }
  }
  return a;
}

void Network::print_assessment(Event &e, double now) {
  const double t0 = std::isnan(first_seen_) ? now : first_seen_;
  const auto a = assess(e, t0, now);
  const std::size_t n = a.stations.size();
  if (static_cast<int>(a.verdict) == e.assessed_verdict && n == e.assessed_stations)
    return;
  e.assessed_verdict = static_cast<int>(a.verdict);
  e.assessed_stations = n;
  Log::get().line("ASSESS", a.verdict == Verdict::Misfire ? "1;33" : "1;36",
                  "#{:<3} {}  {}: {}", e.id, hms(now), verdict_name(a.verdict),
                  join(a.reasons, "; "));
}

bool Network::write_assessments(const std::string &path, double t0, double t1,
                                const std::string &record) const {
  std::ofstream f(path);
  if (!f)
    return false;
  // time_utc, station_a, station_b, dt_s and record are what
  // tools/plot_candidates.py reads (the columns of tools/misfire_review.py).
  f << "event,alarm_utc,matched,verdict,picked_stations,sp,rms_s,lat,lon,"
       "origin_utc,geo_z_median,geo_offset_km,mag_spread,p_chance,context,"
       "reasons,record,time_utc,station_a,station_b,dt_s\n";
  for (const auto &e : events_) {
    if (!e.declared)
      continue;
    const auto a = assess(e, t0, t1);
    std::string sp;
    for (const auto &s : a.stations)
      sp += std::format("{}{}:{:.2f}s:{:.0f}km{}", sp.empty() ? "" : " ",
                        s.station, s.sp, s.km, s.used ? "" : ":out");
    std::vector<const Detection *> d;
    for (const auto &[_, x] : e.detections)
      d.push_back(&x);
    std::ranges::sort(d, {}, &Detection::window_start);
    const auto *c = match(e);
    auto num = [](double v, int prec) {
      return std::isnan(v) ? std::string() : std::format("{:.{}f}", v, prec);
    };
    f << std::format(
        "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{:.2f},{},{},{}\n",
        e.id, ymd_hms(e.declared_at),
        c ? std::format("{} {:.1f} {}", c->type, c->magnitude, hms(c->time))
          : std::string(),
        verdict_name(a.verdict), a.stations.size(), sp, num(a.rms, 2),
        num(a.lat, 4), num(a.lon, 4),
        std::isnan(a.origin) ? std::string() : ymd_hms(a.origin),
        num(a.geo_z_median, 2), num(a.geo_offset_km, 1), num(a.mag_spread, 2),
        std::isnan(a.p_chance) ? std::string() : std::format("{:.3g}", a.p_chance),
        csv_text(a.context), csv_text(join(a.reasons, " | ")), record,
        d.front()->window_start, d.front()->station,
        d.size() > 1 ? d[1]->station : std::string(),
        d.size() > 1 ? num(d[1]->window_start - d.front()->window_start, 2)
                     : std::string());
  }
  return true;
}

} // namespace ayzek::pipeline
