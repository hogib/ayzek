#include "network.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <numbers>
#include <numeric>
#include <ranges>
#include <sstream>
#include <stdexcept>

namespace ayzek::pipeline {

double distance_km(double lat1, double lon1, double lat2, double lon2) {
  const double r = std::numbers::pi / 180.0;
  const double a = std::pow(std::sin((lat2 - lat1) * r / 2), 2) +
                   std::cos(lat1 * r) * std::cos(lat2 * r) *
                       std::pow(std::sin((lon2 - lon1) * r / 2), 2);
  return 2 * 6371.0 * std::asin(std::sqrt(a));
}

std::vector<CatalogEvent> load_afad_catalog(const std::string &path, double t0,
                                            double t1) {
  std::ifstream f(path);
  if (!f)
    throw std::runtime_error("cannot open " + path);
  std::vector<CatalogEvent> out;
  std::string line;
  std::getline(f, line);
  while (std::getline(f, line)) {
    std::stringstream ss(line);
    std::string date, lon, lat, depth, rms, type, mag;
    std::getline(ss, date, ',');
    std::getline(ss, lon, ',');
    std::getline(ss, lat, ',');
    std::getline(ss, depth, ',');
    std::getline(ss, rms, ',');
    std::getline(ss, type, ',');
    std::getline(ss, mag, ',');
    const auto parsed = parse_utc(date, "dmyhms");
    if (!parsed)
      continue;
    const double epoch = *parsed;
    if (epoch < t0 || epoch > t1)
      continue;
    try {
      out.push_back({epoch, std::stod(lat), std::stod(lon), std::stod(depth),
                     std::stod(mag), type});
    } catch (const std::exception &) {
    }
  }
  std::ranges::sort(out, {}, &CatalogEvent::time);
  return out;
}

Network::Network(std::map<std::string, StationInfo> stations, NetworkConfig cfg)
    : stations_(std::move(stations)), cfg_(std::move(cfg)) {}

// Two detections are compatible if their window times differ by no more than
// the P travel time between the two stations plus `slack_seconds`, which allows
// for the position of P within each 6 s window.
bool Network::compatible(const Detection &a, const Detection &b) const {
  auto ia = stations_.find(a.station), ib = stations_.find(b.station);
  if (ia == stations_.end() || ib == stations_.end())
    return std::abs(a.window_start - b.window_start) <= 30.0;
  const double km = distance_km(ia->second.lat, ia->second.lon, ib->second.lat,
                                ib->second.lon);
  return std::abs(a.window_start - b.window_start) <=
         km / cfg_.vp + cfg_.slack_seconds;
}

double Network::one_source_rms(
    const std::vector<std::pair<std::string, double>> &obs) const {
  // P dates (window start + 3.5 s) against one epicentre at fixed depth,
  // within merge_radius_deg of the stations, and its least-squares origin;
  // grid search as for the picks locator.
  struct Xy {
    double lat, lon, t;
  };
  std::vector<Xy> x;
  double clat = 0, clon = 0;
  for (const auto &[code, t] : obs) {
    const auto &s = stations_.at(code);
    x.push_back({s.lat, s.lon, t});
    clat += s.lat;
    clon += s.lon;
  }
  clat /= static_cast<double>(x.size());
  clon /= static_cast<double>(x.size());
  std::vector<double> r(x.size());
  auto rms = [&](double lat, double lon) {
    double sum = 0, sq = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
      r[i] = x[i].t - std::hypot(distance_km(lat, lon, x[i].lat, x[i].lon),
                                 cfg_.depth_km) /
                          cfg_.vp;
      sum += r[i];
    }
    const double origin = sum / static_cast<double>(r.size());
    for (double v : r)
      sq += (v - origin) * (v - origin);
    return std::sqrt(sq / static_cast<double>(r.size()));
  };
  double best = INFINITY, blat = clat, blon = clon;
  auto search = [&](double lat0, double lon0, double half, double step) {
    const auto n = static_cast<long>(std::floor(2 * half / step + 1e-9)) + 1;
    for (long a = 0; a < n; ++a)
      for (long b = 0; b < n; ++b) {
        const double la = lat0 - half + static_cast<double>(a) * step;
        const double lo = lon0 - half + static_cast<double>(b) * step;
        if (const double m = rms(la, lo); m < best)
          best = m, blat = la, blon = lo;
      }
  };
  search(clat, clon, cfg_.merge_radius_deg, 0.05);
  search(blat, blon, 0.075, 0.005);
  return best;
}

// `e2` is about to be declared. If its P times and those of a recently
// declared event fit one source, it is that event seen again: a station that
// triggered twice within one earthquake (once early or on its coda, once on
// its P) splits the stations between two events. Then e2's detections join
// the declared event and e2 is not declared.
//
// One detection per station: where both events have one, each choice is
// tried. At most one detection may be left out (the stray trigger), and the
// fit needs merge_min_fit stations, at least two of them e2's own and one the
// declared event's. A declared event whose own four or more P times fit one
// source loses at most one detection to the merge; so two real earthquakes
// seconds apart, each consistent on its own, are not merged.
bool Network::merge_duplicate(Event &e2) {
  auto p_date = [](const Detection &d) { return d.window_start + 3.5; };
  for (const auto &[code, _] : e2.detections)
    if (!stations_.contains(code))
      return false;
  const double t2 =
      std::ranges::min(e2.detections | std::views::values, {}, p_date)
          .window_start;
  for (auto &e1 : events_ | std::views::reverse) {
    if (!e1.declared || &e1 == &e2)
      continue;
    const double t1 =
        std::ranges::min(e1.detections | std::views::values, {}, p_date)
            .window_start;
    if (std::abs(t2 - t1) > cfg_.merge_window)
      continue;
    if (std::ranges::any_of(e1.detections, [&](const auto &kv) {
          return !stations_.contains(kv.first);
        }))
      continue;
    std::vector<std::string> shared;
    for (const auto &[code, _] : e2.detections)
      if (e1.detections.contains(code))
        shared.push_back(code);
    // A declared event whose own four or more P times fit one source is an
    // earthquake in its own right: the merge may drop at most one of its
    // detections (the stray). Three always fit, so they say nothing.
    std::size_t max_removed = e1.detections.size();
    if (e1.detections.size() >= 4) {
      std::vector<std::pair<std::string, double>> own1;
      for (const auto &[code, d] : e1.detections)
        own1.emplace_back(code, p_date(d));
      if (one_source_rms(own1) <= cfg_.merge_rms)
        max_removed = 1;
    }
    if (shared.size() > 4)
      continue;
    struct Obs {
      std::string code;
      double t;
      bool own; // e2's
    };
    for (unsigned combo = 0; combo < (1u << shared.size()); ++combo) {
      std::vector<Obs> obs;
      for (const auto &[code, d] : e1.detections) {
        const auto k = std::ranges::find(shared, code) - shared.begin();
        if (k == static_cast<long>(shared.size()) || !(combo >> k & 1u))
          obs.push_back({code, p_date(d), false});
      }
      for (const auto &[code, d] : e2.detections) {
        const auto k = std::ranges::find(shared, code) - shared.begin();
        if (k == static_cast<long>(shared.size()) || (combo >> k & 1u))
          obs.push_back({code, p_date(d), true});
      }
      // -1: all of them; i: all but obs[i].
      for (long out = -1; out < static_cast<long>(obs.size()); ++out) {
        std::vector<std::pair<std::string, double>> fit;
        std::size_t own = 0, theirs = 0;
        const std::size_t removed =
            static_cast<std::size_t>(std::ranges::count_if(
                shared, [&](const std::string &c) {
                  return combo >> (std::ranges::find(shared, c) -
                                   shared.begin()) & 1u;
                })) +
            (out >= 0 && !obs[out].own ? 1 : 0);
        if (removed > max_removed)
          continue;
        for (long i = 0; i < static_cast<long>(obs.size()); ++i) {
          if (i == out)
            continue;
          fit.emplace_back(obs[i].code, obs[i].t);
          (obs[i].own ? own : theirs) += 1;
        }
        if (fit.size() < cfg_.merge_min_fit || own < 2 || theirs < 1)
          continue;
        const double rms = one_source_rms(fit);
        if (rms > cfg_.merge_rms)
          continue;
        // Merge: e1 keeps the fitted detections, with whatever it knew from
        // their stations (geometry, magnitude, picks); the rest is dropped.
        std::string moved;
        for (long i = 0; i < static_cast<long>(obs.size()); ++i) {
          const auto &o = obs[i];
          if (i == out) {
            if (!o.own) {
              e1.detections.erase(o.code);
              e1.geometry.erase(o.code);
              e1.magnitudes.erase(o.code);
              e1.picks.erase(o.code);
              e1.pd.erase(o.code);
            }
            continue;
          }
          if (!o.own)
            continue;
          e1.detections.insert_or_assign(o.code, e2.detections.at(o.code));
          if (auto it = e2.geometry.find(o.code); it != e2.geometry.end())
            e1.geometry.insert_or_assign(o.code, it->second);
          else
            e1.geometry.erase(o.code);
          if (auto it = e2.magnitudes.find(o.code); it != e2.magnitudes.end())
            e1.magnitudes.insert_or_assign(o.code, it->second);
          else
            e1.magnitudes.erase(o.code);
          if (auto it = e2.picks.find(o.code); it != e2.picks.end())
            e1.picks.insert_or_assign(o.code, it->second);
          else
            e1.picks.erase(o.code);
          if (auto it = e2.pd.find(o.code); it != e2.pd.end())
            e1.pd.insert_or_assign(o.code, it->second);
          else
            e1.pd.erase(o.code);
          moved += (moved.empty() ? "" : ", ") + o.code;
        }
        const double now =
            std::ranges::max(e2.detections | std::views::values, {},
                             &Detection::declared_at)
                .declared_at;
        Log::get().line("MERGE", "33",
                        "#{:<3} {}  {} joined: one source with its P times, "
                        "rms {:.2f} s{}",
                        e1.id, hms(now), moved, rms,
                        out >= 0 ? std::format(", {} left out", obs[out].code)
                                 : std::string());
        Event &kept = e1;
        std::erase_if(events_, [&](const Event &e) { return &e == &e2; });
        report_magnitude(kept, now);
        if (cfg_.locator == Locator::Geometry)
          report_location(kept, now);
        report_pd(kept, now);
        return true;
      }
    }
  }
  return false;
}

void Network::on(const Detection &d) {
  ++triggers_[d.station];
  if (std::isnan(first_seen_))
    first_seen_ = d.window_start;
  // The detector can trigger again on the S wave and coda of an event. A later
  // detection at a station within `coda_seconds` of that station's detection
  // of a declared event is assigned to that event. A separate event in that
  // interval is therefore not declared from this station, unless the
  // detection is a transformer dt restart: the model itself says a new onset
  // began inside the coda, which is what an aftershock looks like.
  for (auto &e : events_) {
    if (d.restart)
      break;
    auto it = e.detections.find(d.station);
    if (e.declared && it != e.detections.end() &&
        d.window_start > it->second.window_start &&
        d.window_start - it->second.window_start <= cfg_.coda_seconds) {
      ++e.coda;
      if (cfg_.verbose)
        Log::get().line("detect", "2",
                        "{:<5} p={:.2f}  window {}  -- coda of #{}", d.station,
                        d.probability, hms(d.window_start), e.id);
      return;
    }
  }
  if (cfg_.verbose)
    Log::get().line("detect", "36", "{:<5} p={:.2f}  window {}  ({:.1f} ms){}",
                    d.station, d.probability, hms(d.window_start),
                    d.compute_ms, d.restart ? "  new onset in a coda" : "");

  // The first compatible event, or with merge_duplicates the one with the
  // most stations (declared first): a stray early trigger at one station must
  // not keep the rest of an earthquake's stations in an event of its own.
  Event *ev = nullptr;
  for (auto &e : events_) {
    if (e.detections.contains(d.station))
      continue;
    if (std::ranges::all_of(
            e.detections | std::views::values,
            [&](const Detection &o) { return compatible(o, d); })) {
      if (!cfg_.merge_duplicates) {
        ev = &e;
        break;
      }
      if (!ev || std::pair(e.declared, e.detections.size()) >
                     std::pair(ev->declared, ev->detections.size()))
        ev = &e;
    }
  }
  if (!ev) {
    events_.push_back(Event{});
    ev = &events_.back();
  }
  ev->detections.emplace(d.station, d);

  if (cfg_.merge_duplicates && !ev->declared &&
      ev->detections.size() >= cfg_.min_stations && merge_duplicate(*ev))
    return;
  if (!ev->declared && ev->detections.size() >= cfg_.min_stations) {
    ev->declared = true;
    ev->id = next_id_++;
    ev->declared_at = std::ranges::max(ev->detections | std::views::values, {},
                                       &Detection::declared_at)
                          .declared_at;
    std::string names;
    for (const auto &[code, _] : ev->detections)
      names += (names.empty() ? "" : ", ") + code;
    const auto *c = match(*ev);
    Log::get().line("ALARM", "1;31", "#{:<3} {}  earthquake detected by {}{}",
                    ev->id, hms(ev->declared_at), names,
                    c ? std::format("  (AFAD {} {:.1f} at {}, {:.1f} s ago)",
                                    c->type, c->magnitude, hms(c->time),
                                    ev->declared_at - c->time)
                      : "");
    report_magnitude(
        *ev,
        ev->declared_at); // estimates received before the event was declared
    if (cfg_.locator == Locator::Geometry)
      report_location(*ev, ev->declared_at); // likewise
    report_pd(*ev, ev->declared_at);         // likewise
  } else if (ev->declared && cfg_.verbose) {
    Log::get().line("EVENT", "31", "#{} joined by {}", ev->id, d.station);
  }
}

void Network::on(const Pick &p) {
  const double sp = p.s_time - p.p_time;
  const double km =
      sp * cfg_.vp * cfg_.vs / (cfg_.vp - cfg_.vs); // S-P lag to distance

  // Pick quality control. The picker returns the most probable chunk even when
  // there is no phase, so low-probability picks, S before P, and P far from the
  // trigger window are not used.
  std::string reject;
  if (p.p_prob < cfg_.min_pick_prob || p.s_prob < cfg_.min_pick_prob)
    reject = "low confidence";
  else if (sp <= 0 || sp > 60)
    reject = "S not after P";
  else if (p.p_time < p.trigger_window - 3 || p.p_time > p.trigger_window + 9)
    reject = "P outside the trigger window";
  // The event a pick belongs to: the one with the detection it follows, or,
  // for an --assess pick at a station that did not trigger, the event it was
  // made for (by id and alarm time; failing that, as when a recording is
  // replayed for a subset of its stations, the declared event with a
  // detection nearest its predicted window).
  Event *target = nullptr;
  if (p.event_id) {
    double nearest = 40.0;
    for (auto &e : events_) {
      if (!e.declared || e.picks.contains(p.station))
        continue;
      if (e.id == p.event_id && e.declared_at == p.event_alarm) {
        target = &e;
        break;
      }
      for (const auto &[_, d] : e.detections)
        if (std::abs(d.window_start - p.trigger_window) < nearest) {
          nearest = std::abs(d.window_start - p.trigger_window);
          target = &e;
        }
    }
  } else {
    for (auto &e : events_) {
      auto it = e.detections.find(p.station);
      if (it != e.detections.end() &&
          it->second.window_start == p.trigger_window) {
        target = &e;
        break;
      }
    }
  }
  if (!reject.empty()) {
    if (cfg_.verbose)
      Log::get().line(
          "pick", "2", "{:<5} P {} ({:.2f})  S {} ({:.2f})  -- not used: {}",
          p.station, hms(p.p_time), p.p_prob, hms(p.s_time), p.s_prob, reject);
    if (target) {
      ++target->picks_rejected;
      if (target->declared && cfg_.assess)
        print_assessment(*target, p.declared_at);
    }
    return;
  }
  if (cfg_.verbose)
    Log::get().line("pick", "35",
                    "{:<5} P {} ({:.2f})  S {} ({:.2f})  S-P {:.2f} s ~ {:.0f} "
                    "km  ({:.0f} ms){}",
                    p.station, hms(p.p_time), p.p_prob, hms(p.s_time), p.s_prob,
                    sp, km, p.compute_ms,
                    p.event_id ? std::format("  for #{}, not triggered", p.event_id)
                               : std::string());
  if (!target)
    return;
  target->picks.emplace(p.station, p);
  // A pick made for the assessment only informs the assessment: the event
  // keeps the location the run gave it.
  if (target->declared && cfg_.locator == Locator::Picks && !p.event_id)
    report_location(*target, p.declared_at);
  if (target->declared && cfg_.assess)
    print_assessment(*target, p.declared_at);
}

std::vector<PickRequest> Network::pick_requests(double t0, double t1) const {
  std::vector<PickRequest> out;
  for (const auto &e : events_) {
    if (!e.declared)
      continue;
    // The source: the S-P one if the picks agree on it, else the location.
    const auto a = assess(e, t0, t1);
    double lat = a.lat, lon = a.lon, origin = a.origin;
    if (!a.consistent) {
      if (!e.location)
        continue;
      lat = e.location->lat;
      lon = e.location->lon;
      origin = e.location->origin;
    }
    for (const auto &[code, st] : stations_) {
      if (e.detections.contains(code) || e.picks.contains(code))
        continue;
      const double km = distance_km(lat, lon, st.lat, st.lon);
      if (km > cfg_.assess_max_sp_km)
        continue;
      const double p = origin + std::hypot(km, cfg_.depth_km) / cfg_.vp;
      // S no later than this station's next detected onset.
      double s_cap = INFINITY;
      for (const auto &o : events_)
        if (auto it = o.detections.find(code); it != o.detections.end()) {
          const double q = it->second.window_start + 3.5;
          if (q > p + 1.0)
            s_cap = std::min(s_cap, q - 0.5);
        }
      out.push_back({e.id, e.declared_at, code, p, s_cap});
    }
  }
  return out;
}

void Network::on(const StationGeometry &g) {
  if (cfg_.verbose)
    Log::get().line(
        "geo", "35",
        "{:<5} P {} +{:.1f} s  {:.0f} km (x/{:.2f})", g.station,
        hms(g.p_time), g.since_p, std::exp(g.log_dist),
        std::exp(g.log_dist_sd));
  for (auto &e : events_) {
    auto it = e.detections.find(g.station);
    if (it == e.detections.end() || it->second.window_start != g.trigger_window)
      continue;
    e.geometry.insert_or_assign(g.station, g);
    if (e.declared && cfg_.locator == Locator::Geometry)
      report_location(e, g.declared_at);
    return;
  }
}

// Locates the event; while rms exceeds max_rms and more stations remain than
// the locator needs (two with picks, kMinGeometryStations with geometry), removes the station with the largest residual and relocates. The
// uniform velocity model underestimates P speed beyond ~150 km, where the first
// arrival travels through the upper mantle, so distant stations are the ones
// typically removed.
void Network::on(const MagnitudeEstimate &m) {
  if (cfg_.verbose)
    Log::get().line(
        "mag", "34", "{:<5} M{:.2f}  {}  noise baseline {}  ({:.0f} ms)",
        m.station, m.magnitude, m.at_pick ? "at picked P" : "early      ",
        m.noise_windows ? std::format("{} windows", m.noise_windows)
                        : std::string("not ready, per-window"),
        m.compute_ms);
  for (auto &e : events_) {
    auto it = e.detections.find(m.station);
    if (it == e.detections.end() || it->second.window_start != m.trigger_window)
      continue;
    // Estimates at the picked P are made for declared events only
    // (app/ayzek.cpp). A recording (--record) contains them for every pick;
    // ignoring the others here makes a replay of the recording agree with a
    // run on the same stations.
    if (m.at_pick && !e.declared)
      return;
    auto have = e.magnitudes.find(m.station);
    if (have != e.magnitudes.end() && have->second.at_pick && !m.at_pick)
      return;
    e.magnitudes.insert_or_assign(m.station, m);
    if (e.declared)
      report_magnitude(e, m.declared_at);
    return;
  }
}

void Network::on(const PdEstimate &q) {
  if (cfg_.verbose)
    Log::get().line("pd", "34", "{:<5} P {} +{:g} s  Pd {:.2e} m  noise {:.2e} m",
                    q.station, hms(q.p_time), q.tau, q.pd, q.pd_noise);
  for (auto &e : events_) {
    auto it = e.detections.find(q.station);
    if (it == e.detections.end() || it->second.window_start != q.trigger_window)
      continue;
    auto have = e.pd.find(q.station);
    if (have == e.pd.end() || q.tau >= have->second.tau)
      e.pd.insert_or_assign(q.station, q);
    if (e.declared)
      report_pd(e, q.declared_at);
    return;
  }
}

void Network::report_pd(Event &e, double now) {
  if (!cfg_.pd || e.pd.empty())
    return;
  std::vector<PdObservation> obs;
  double tau = 0;
  const bool located = e.location && !(e.location->err_km > cfg_.pd_max_err_km);
  for (const auto &[code, q] : e.pd) {
    auto st = stations_.find(code);
    if (st == stations_.end())
      continue;
    double d = NAN;
    if (located)
      d = distance_km(e.location->lat, e.location->lon, st->second.lat,
                      st->second.lon);
    else if (auto g = e.geometry.find(code); g != e.geometry.end())
      d = std::exp(g->second.log_dist);
    if (std::isnan(d))
      continue;
    obs.push_back({code, q.tau, q.pd, q.pd_noise, d});
    tau = std::max(tau, q.tau);
  }
  const auto m = pd_magnitude(*cfg_.pd, obs);
  if (!m)
    return;
  const bool changed = !e.pd_magnitude || std::isnan(e.pd_reported) ||
                       std::abs(m->mean - e.pd_reported) >= 0.1 ||
                       m->stations != e.pd_magnitude->stations;
  e.pd_magnitude = m;
  if (!e.pd_first) {
    e.pd_first = m->mean;
    e.pd_first_at = now;
  }
  if (!changed)
    return;
  e.pd_reported = m->mean;
  Log::get().line("PDMAG", "1;34",
                  "#{:<3} {}  Pd magnitude M{:.1f} ± {:.1f} from {} station{} "
                  "({} above noise), windows up to {:g} s, distances from {}",
                  e.id, hms(now), m->mean, m->sd, m->stations,
                  m->stations == 1 ? "" : "s", m->above, tau,
                  located ? "the location" : "the stations' geometry");
}

bool Network::declared(const std::string &station,
                       double trigger_window) const {
  return std::ranges::any_of(events_, [&](const Event &e) {
    auto it = e.detections.find(station);
    return e.declared && it != e.detections.end() &&
           it->second.window_start == trigger_window;
  });
}

// Event magnitude: median of the station estimates. Printed when the value
// changes by 0.05 or more, or when the number of stations changes.
void Network::report_magnitude(Event &e, double now) {
  if (e.magnitudes.empty())
    return;
  std::vector<double> v;
  std::string parts;
  for (const auto &[code, sm] : e.magnitudes) {
    v.push_back(sm.magnitude);
    parts += std::format("{}{} {:.1f}", parts.empty() ? "" : ", ", code,
                         sm.magnitude);
  }
  std::ranges::sort(v);
  const double med = v.size() % 2
                         ? v[v.size() / 2]
                         : 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]);
  const bool changed = !e.magnitude || std::abs(*e.magnitude - med) >= 0.05 ||
                       e.magnitudes.size() != e.magnitude_stations;
  e.magnitude = med;
  e.magnitude_stations = e.magnitudes.size();
  if (!e.first_magnitude) {
    e.first_magnitude = med;
    e.first_magnitude_at = now;
  }
  if (!changed)
    return;
  Log::get().line("MAG", "1;34",
                  "#{:<3} {}  magnitude M{:.1f} from {} station{} ({})", e.id,
                  hms(now), med, v.size(), v.size() == 1 ? "" : "s", parts);
}

std::optional<Location> Network::locate(const Event &e) const {
  const std::size_t min_stations =
      cfg_.locator == Locator::Geometry ? kMinGeometryStations : 2;
  auto relocate = [&](auto obs) -> std::optional<Location> {
    std::vector<std::string> dropped;
    for (;;) {
      std::string worst;
      auto loc = locate(obs, &worst);
      if (!loc || !failing(*loc) || obs.size() <= min_stations) {
        if (loc)
          loc->dropped = dropped;
        return loc;
      }
      dropped.push_back(worst);
      obs.erase(worst);
    }
  };
  return cfg_.locator == Locator::Geometry ? relocate(e.geometry)
                                           : relocate(e.picks);
}

// Geometry locator (locate.hpp): each station's latest distance and P time. `worst` is the station contributing most to the misfit.
std::optional<Location>
Network::locate(const std::map<std::string, StationGeometry> &geo,
                std::string *worst) const {
  std::vector<GeometryObs> obs;
  for (const auto &[code, g] : geo) {
    auto it = stations_.find(code);
    if (it == stations_.end())
      continue;
    obs.push_back({code, it->second.lat, it->second.lon, g.p_time, g.log_dist,
                   g.log_dist_sd * cfg_.geo_sd_scale});
  }
  const auto fit = locate_geometry(
      obs, {.vp = cfg_.vp, .depth_km = cfg_.depth_km, .sigma_p = cfg_.sigma_p},
      worst);
  if (!fit)
    return std::nullopt;
  Location loc{fit->lat, fit->lon, fit->origin, fit->rms, fit->n_stations, {}};
  loc.err_km = fit->err_km;
  loc.dist_z = fit->dist_z_max;
  return loc;
}

std::optional<Location>
Network::locate(const std::map<std::string, Pick> &picks,
                std::string *worst) const {
  struct Obs {
    double lat, lon, t, v;
    const std::string *station;
  };
  std::vector<Obs> obs;
  double clat = 0, clon = 0;
  for (const auto &[code, pk] : picks) {
    auto it = stations_.find(code);
    if (it == stations_.end())
      continue;
    obs.push_back({it->second.lat, it->second.lon, pk.p_time, cfg_.vp, &code});
    obs.push_back({it->second.lat, it->second.lon, pk.s_time, cfg_.vs, &code});
    clat += it->second.lat;
    clon += it->second.lon;
  }
  const std::size_t n_st = obs.size() / 2;
  if (n_st < 2)
    return std::nullopt; // 2 stations x (P, S) = 4 times for 3 unknowns
  clat /= static_cast<double>(n_st);
  clon /= static_cast<double>(n_st);

  // Grid search over epicentre; depth fixed. For a given epicentre the least-
  // squares origin time is the mean of (observed - predicted travel time).
  std::vector<double> r(obs.size());
  auto misfit = [&](double lat, double lon, double &origin) {
    double sum = 0, sq = 0;
    for (std::size_t i = 0; i < obs.size(); ++i) {
      r[i] =
          obs[i].t - std::hypot(distance_km(lat, lon, obs[i].lat, obs[i].lon),
                                cfg_.depth_km) /
                         obs[i].v;
      sum += r[i];
    }
    origin = sum / static_cast<double>(r.size());
    for (double x : r)
      sq += (x - origin) * (x - origin);
    return std::sqrt(sq / static_cast<double>(r.size()));
  };
  Location best{clat, clon, 0, 1e18, n_st, {}};
  auto search = [&](double lat0, double lon0, double half, double step) {
    for (double la = lat0 - half; la <= lat0 + half; la += step)
      for (double lo = lon0 - half; lo <= lon0 + half; lo += step) {
        double origin = 0;
        const double m = misfit(la, lo, origin);
        if (m < best.rms)
          best = {la, lo, origin, m, n_st, {}};
      }
  };
  search(clat, clon, 3.0, 0.05);
  search(best.lat, best.lon, 0.1, 0.005);

  if (worst) {
    double origin = 0, largest = -1;
    misfit(best.lat, best.lon, origin);
    for (std::size_t i = 0; i < obs.size(); ++i)
      if (std::abs(r[i] - origin) > largest) {
        largest = std::abs(r[i] - origin);
        *worst = *obs[i].station;
      }
  }
  return best;
}

std::size_t Network::geometry_stations(const Event &e) const {
  return static_cast<std::size_t>(std::ranges::count_if(
      e.geometry, [&](const auto &kv) { return stations_.contains(kv.first); }));
}

std::string Network::ranges(const Event &e) const {
  std::string out;
  for (const auto &[code, g] : e.geometry) {
    if (!stations_.contains(code))
      continue;
    const double sd = g.log_dist_sd * cfg_.geo_sd_scale;
    out += std::format("{}{:.0f} km from {} ({:.0f}-{:.0f} km)",
                       out.empty() ? "" : ", ", std::exp(g.log_dist), code,
                       std::exp(g.log_dist - sd), std::exp(g.log_dist + sd));
  }
  return out;
}

// Too few stations to locate: report how far the event is from each one
// that has an estimate. Prints when a station joins or a distance moves by
// more than 10%, as estimates arrive every second per station.
void Network::report_ranges(Event &e, double now) {
  std::map<std::string, double> now_ranged;
  for (const auto &[code, g] : e.geometry)
    if (stations_.contains(code))
      now_ranged[code] = g.log_dist;
  if (now_ranged.empty())
    return;
  bool same = now_ranged.size() == e.ranged.size();
  for (const auto &[code, ld] : now_ranged) {
    auto it = e.ranged.find(code);
    same = same && it != e.ranged.end() && std::abs(it->second - ld) <= 0.1;
  }
  if (same)
    return;
  e.ranged = std::move(now_ranged);
  Log::get().line("RANGE", "1;36", "#{:<3} {}  {}; {} of {} stations to locate",
                  e.id, hms(now), ranges(e), e.ranged.size(),
                  kMinGeometryStations);
}

void Network::report_location(Event &e, double now) {
  auto loc = locate(e);
  if (!loc) {
    // The geometry locator refuses only for too few stations.
    if (cfg_.locator == Locator::Geometry && !e.location)
      report_ranges(e, now);
    return;
  }
  if (failing(*loc) || loc->err_km > cfg_.geo_max_err_km) {
    if (cfg_.verbose)
      Log::get().line("LOCATE", "33",
                      "#{} rejected: best fit rms {:.1f} s{} from {} stations",
                      e.id, loc->rms,
                      std::isnan(loc->dist_z)
                          ? std::string()
                          : std::format(", a station's distance {:.1f} sd "
                                        "off, +-{:.0f} km",
                                        loc->dist_z, loc->err_km),
                      loc->n_stations);
    return;
  }
  // Do not print an unchanged solution.
  // Geometry estimates arrive every second per station: print a solution
  // only when it moves by 2 km or more, or its station count changes.
  const bool same =
      e.location &&
      (cfg_.locator == Locator::Geometry
           ? distance_km(e.location->lat, e.location->lon, loc->lat,
                         loc->lon) < 2.0 &&
                 e.location->n_stations == loc->n_stations
           : e.location->lat == loc->lat && e.location->lon == loc->lon &&
                 e.location->origin == loc->origin);
  if (!e.location)
    e.first_located_at = now;
  e.location = loc;
  report_pd(e, now); // the distances changed
  if (same)
    return;
  std::string left_out;
  for (const auto &d : loc->dropped)
    left_out += (left_out.empty() ? ", left out " : " ") + d;
  Log::get().line("LOCATE", "1;32",
                  "#{:<3} {}  located at {:.3f}N {:.3f}E{}, origin {}, rms "
                  "{:.2f} s from {} station{}{}",
                  e.id, hms(now), loc->lat, loc->lon,
                  std::isnan(loc->err_km)
                      ? std::string()
                      : std::format(" +-{:.0f} km", loc->err_km),
                  hms(loc->origin), loc->rms, loc->n_stations,
                  loc->n_stations == 1 ? "" : "s", left_out);
}

// Matching catalogue event, if any. A catalogue event matches if the location
// agrees (origin within 5 s, epicentre within 30 km), or if its predicted P
// arrival at a station that raised the alarm is within 6 s of that detection's
// window start plus 3.5 s. The second criterion does not require a location.
// Detections that joined after the alarm are not used for it: they show that
// the event was recorded, not that the alarm was raised by it.
const CatalogEvent *Network::match(const Event &e) const {
  const CatalogEvent *best = nullptr;
  double best_miss = 1e18;
  for (const auto &c : cfg_.catalog) {
    double miss = 1e18;
    if (e.location &&
        distance_km(e.location->lat, e.location->lon, c.lat, c.lon) <= 30.0 &&
        std::abs(e.location->origin - c.time) <= 5.0)
      miss = std::abs(e.location->origin - c.time) /
             100.0; // location matches rank first
    for (const auto &[code, d] : e.detections) {
      auto it = stations_.find(code);
      if (it == stations_.end() || d.declared_at > e.declared_at)
        continue;
      const double p =
          c.time +
          std::hypot(distance_km(c.lat, c.lon, it->second.lat, it->second.lon),
                     c.depth) /
              cfg_.vp;
      miss = std::min(miss, std::abs(p - (d.window_start + 3.5)));
    }
    if (miss <= 6.0 && miss < best_miss) {
      best = &c;
      best_miss = miss;
    }
  }
  return best;
}

std::vector<CatalogScore> Network::score() const {
  std::vector<CatalogScore> out;
  if (cfg_.catalog.empty() || stations_.empty())
    return out;
  double clat = 0, clon = 0;
  for (const auto &[_, s] : stations_) {
    clat += s.lat;
    clon += s.lon;
  }
  clat /= static_cast<double>(stations_.size());
  clon /= static_cast<double>(stations_.size());
  for (const auto &c : cfg_.catalog) {
    const double dist = distance_km(clat, clon, c.lat, c.lon);
    if (dist > cfg_.catalog_radius_km)
      continue;
    const Event *hit = nullptr;
    for (const auto &e : events_)
      if (e.declared && match(e) == &c) {
        hit = &e;
        break;
      }
    out.push_back({&c, dist, hit});
  }
  return out;
}

std::size_t Network::declared_count() const {
  return static_cast<std::size_t>(
      std::ranges::count_if(events_, &Event::declared));
}

std::size_t Network::unmatched_count() const {
  const auto scores = score();
  std::size_t n = 0;
  for (const auto &e : events_)
    if (e.declared && std::ranges::none_of(scores, [&](const CatalogScore &s) {
          return s.detected == &e;
        }))
      ++n;
  return n;
}

std::vector<Network::Warning> Network::warnings(const Event &e,
                                                const CatalogEvent *c) const {
  std::vector<Warning> out;
  double lat = 0, lon = 0, depth = cfg_.depth_km, origin = 0;
  if (c) {
    lat = c->lat;
    lon = c->lon;
    depth = c->depth;
    origin = c->time;
  } else if (e.location) {
    lat = e.location->lat;
    lon = e.location->lon;
    origin = e.location->origin;
  } else {
    return out;
  }
  auto add = [&](const std::string &name, double slat, double slon, bool site) {
    const double km = distance_km(lat, lon, slat, slon);
    Warning w{name, km, origin + std::hypot(km, depth) / cfg_.vs, false, site};
    if (!site)
      if (auto it = e.picks.find(name); it != e.picks.end()) {
        w.s_time = it->second.s_time;
        w.picked = true;
      }
    out.push_back(w);
  };
  for (const auto &[code, s] : stations_)
    add(code, s.lat, s.lon, false);
  for (const auto &site : cfg_.sites)
    add(site.name, site.lat, site.lon, true);
  std::ranges::sort(out, {}, &Warning::distance_km);
  return out;
}

void Network::report(double t_first, double t_last) const {
  auto &log = Log::get();
  const std::string rule(78, '-');
  log.plain(false, "");
  log.plain(true, "{}", rule);
  log.plain(true, "REPORT  {} stations, {} to {} UTC", stations_.size(),
            ymd_hms(t_first), hms(t_last));
  log.plain(true, "{}", rule);

  const auto scores = score();
  std::vector<const Event *> declared;
  for (const auto &e : events_)
    if (e.declared)
      declared.push_back(&e);
  std::ranges::sort(declared, {}, &Event::declared_at);

  std::vector<double> alarm_delays, magnitude_errors, positive_warnings,
      pd_errors; // signed
  std::size_t arrivals = 0, warned = 0;
  for (const Event *e : declared) {
    const CatalogEvent *c = nullptr;
    for (const auto &sc : scores)
      if (sc.detected == e)
        c = sc.event;

    // Headline.
    std::string headline = std::format("Event #{}: detected ", e->id);
    headline += e->magnitude ? std::format("M{:.1f} earthquake", *e->magnitude)
                             : std::string("earthquake (no magnitude)");
    if (e->location)
      headline +=
          std::format(" at {:.3f}N {:.3f}E, origin {} UTC", e->location->lat,
                      e->location->lon, hms(e->location->origin));
    else if (cfg_.locator == Locator::Geometry && geometry_stations(*e) > 0)
      headline += std::format(", not located ({} of {} stations)",
                              geometry_stations(*e), kMinGeometryStations);
    else
      headline += ", not located";
    log.plain(false, "");
    log.plain(true, "{}", headline);

    std::string first, later;
    for (const auto &[code, d] : e->detections)
      (d.declared_at <= e->declared_at ? first : later) += std::format(
          "{}{}",
          (d.declared_at <= e->declared_at ? first : later).empty() ? "" : ", ",
          code);
    log.plain(false, "  alarm      {} UTC by {}{}", hms(e->declared_at), first,
              later.empty() ? "" : std::format("; later {}", later));
    if (e->magnitude) {
      const double lag = e->first_magnitude_at - e->declared_at;
      log.plain(false,
                "  magnitude  M{:.1f} {}, M{:.1f} final from {} station{}",
                *e->first_magnitude,
                lag < 0.05 ? std::string("at the alarm")
                           : std::format("{:.1f} s after the alarm", lag),
                *e->magnitude, e->magnitude_stations,
                e->magnitude_stations == 1 ? "" : "s");
    }
    if (e->pd_magnitude) {
      const double lag = e->pd_first_at - e->declared_at;
      log.plain(false,
                "  Pd magn.   M{:.1f} {}, M{:.1f} ± {:.1f} final from {} "
                "station{} ({} above noise)",
                *e->pd_first,
                lag < 0.05 ? std::string("at the alarm")
                           : std::format("{:.1f} s after the alarm", lag),
                e->pd_magnitude->mean, e->pd_magnitude->sd,
                e->pd_magnitude->stations,
                e->pd_magnitude->stations == 1 ? "" : "s",
                e->pd_magnitude->above);
    }
    if (e->location)
      log.plain(false, "  location   rms {:.2f} s from {} station{}{}, first {:.1f} s "
                "after the alarm ({})",
                e->location->rms, e->location->n_stations,
                e->location->n_stations == 1 ? "" : "s",
                std::isnan(e->location->err_km)
                    ? std::string()
                    : std::format(", +-{:.0f} km", e->location->err_km),
                e->first_located_at - e->declared_at,
                cfg_.locator == Locator::Geometry ? "transformer geometry"
                                                  : "P and S picks");
    else if (cfg_.locator == Locator::Geometry && geometry_stations(*e) > 0)
      log.plain(false, "  distance   {}", ranges(*e));

    const double origin =
        c ? c->time : (e->location ? e->location->origin : NAN);
    if (c) {
      std::string cmp = std::format(
          "  AFAD       {} {:.1f} at {} UTC: alarm {:.1f} s after origin",
          c->type, c->magnitude, hms(c->time), e->declared_at - c->time);
      // Rounded before formatting, so that a difference below 0.05 prints as
      // +0.0.
      if (e->magnitude)
        cmp += std::format(
            ", magnitude {:+.1f}",
            std::round((*e->magnitude - c->magnitude) * 10) / 10 + 0.0);
      if (e->pd_magnitude)
        cmp += std::format(
            ", Pd magnitude {:+.1f}",
            std::round((e->pd_magnitude->mean - c->magnitude) * 10) / 10 + 0.0);
      if (e->location)
        cmp += std::format(
            ", epicentre {:.1f} km off, origin {:+.1f} s",
            distance_km(e->location->lat, e->location->lon, c->lat, c->lon),
            e->location->origin - c->time);
      log.plain(false, "{}", cmp);
      alarm_delays.push_back(e->declared_at - c->time);
      if (e->magnitude)
        magnitude_errors.push_back(std::abs(*e->magnitude - c->magnitude));
      if (e->pd_magnitude)
        pd_errors.push_back(e->pd_magnitude->mean - c->magnitude);
    } else {
      log.plain(false, "  AFAD       no matching catalogue event");
    }
    if (cfg_.assess) {
      const auto a = assess(*e, t_first, t_last);
      log.plain(false, "  assessed   {}: {}", verdict_name(a.verdict),
                a.reasons.empty() ? std::string("-") : a.reasons.front());
      for (std::size_t i = 1; i < a.reasons.size(); ++i)
        log.plain(false, "             {}", a.reasons[i]);
      if (!a.context.empty() && !c)
        log.plain(false, "             context: {}", a.context);
    }

    const auto ws = warnings(*e, c);
    if (ws.empty())
      continue;
    const double depth = c ? c->depth : cfg_.depth_km;
    const double blind =
        std::sqrt(std::max(cfg_.vs * cfg_.vs * (e->declared_at - origin) *
                                   (e->declared_at - origin) -
                               depth * depth,
                           0.0));
    log.plain(false,
              "  S wave     had travelled {:.0f} km from the epicentre when "
              "the alarm sounded ({})",
              blind, c ? "AFAD hypocentre" : "ayzek location");
    log.plain(false, "  warning    {:<10} {:>8}   {:<24} {:>8}", "", "distance",
              "S wave arrives", "warning");
    for (const auto &w : ws) {
      const double lead = w.s_time - e->declared_at;
      log.plain(
          false,
          "             {:<10} {:>5.0f} km   {} UTC {:<10} {:>+7.1f} s  {}",
          w.name, w.distance_km, hms(w.s_time),
          w.picked ? "(picked)" : "(predicted)", lead,
          lead > 0 ? "before S" : "after S");
      if (c && !w.site) {
        ++arrivals;
        if (lead > 0) {
          ++warned;
          positive_warnings.push_back(lead);
        }
      }
    }
  }

  auto median = [](std::vector<double> v) -> double {
    if (v.empty())
      return NAN;
    std::ranges::sort(v);
    return v.size() % 2 ? v[v.size() / 2]
                        : 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]);
  };
  auto pct = [](std::size_t a, std::size_t b) {
    return b ? std::format("{:.0f}%", 100.0 * static_cast<double>(a) /
                                          static_cast<double>(b))
             : std::string("-");
  };

  std::size_t located = 0;
  std::vector<const Event *> unmatched;
  for (const Event *e : declared) {
    located += e->location.has_value();
    if (std::ranges::none_of(
            scores, [&](const CatalogScore &sc) { return sc.detected == e; }))
      unmatched.push_back(e);
  }
  std::vector<const CatalogScore *> missed;
  for (const auto &sc : scores)
    if (!sc.detected)
      missed.push_back(&sc);

  log.plain(false, "");
  log.plain(true, "{}", rule);
  log.plain(true, "SUMMARY");
  log.plain(true, "{}", rule);
  log.plain(false, "  Alarms raised                        {:>4}",
            declared.size());
  log.plain(false, "  single-station detections, no alarm  {:>4}",
            events_.size() - declared.size());
  // --assess: verdict counts, catalogued alarms against the rest. The
  // catalogued ones measure how often a real event passes the assessment.
  auto assessment_summary = [&] {
    if (!cfg_.assess || declared.empty())
      return;
    std::map<Verdict, std::size_t> on_matched, on_unmatched;
    for (const Event *e : declared) {
      const bool matched = std::ranges::any_of(
          scores, [&](const CatalogScore &s) { return s.detected == e; });
      ++(matched || cfg_.catalog.empty() ? on_matched : on_unmatched)
          [assess(*e, t_first, t_last).verdict];
    }
    auto row = [&](const std::map<Verdict, std::size_t> &m) {
      std::string out;
      for (auto v : {Verdict::Earthquake, Verdict::Possible, Verdict::Misfire,
                     Verdict::Unclassified})
        out += std::format("{}{} {}", out.empty() ? "" : ", ",
                           m.contains(v) ? m.at(v) : 0, verdict_name(v));
      return out;
    };
    log.plain(false, "");
    log.plain(false, "  Alarm assessment (S-P test; "
                     "docs/impl/16-alarm-assessment.md)");
    if (cfg_.catalog.empty()) {
      log.plain(false, "    all alarms                          {}",
                row(on_matched));
    } else {
      log.plain(false, "    catalogued alarms                   {}",
                row(on_matched));
      log.plain(false, "    unmatched alarms                    {}",
                row(on_unmatched));
    }
  };
  if (cfg_.catalog.empty()) {
    log.plain(false, "  (no catalogue given: correct, false and missed alarms "
                     "not evaluated)");
    assessment_summary();
    log.plain(true, "{}", rule);
    return;
  }
  const std::size_t correct = declared.size() - unmatched.size();
  log.plain(false, "");
  log.plain(false, "  Alarms compared with the AFAD catalogue");
  log.plain(false, "    correct  (match an AFAD event)     {:>4}  {:>4}",
            correct, pct(correct, declared.size()));
  log.plain(false, "    false    (no AFAD event)           {:>4}  {:>4}",
            unmatched.size(), pct(unmatched.size(), declared.size()));
  log.plain(false, "");
  log.plain(false, "  AFAD events within {:.0f} km             {:>4}",
            cfg_.catalog_radius_km, scores.size());
  log.plain(false, "    detected                           {:>4}  {:>4}",
            scores.size() - missed.size(),
            pct(scores.size() - missed.size(), scores.size()));
  log.plain(false, "    missed                             {:>4}  {:>4}",
            missed.size(), pct(missed.size(), scores.size()));

  // Detection rate by magnitude.
  log.plain(false, "");
  log.plain(false, "  By magnitude      AFAD events   detected   missed");
  const std::array<std::pair<double, double>, 5> bands{
      {{-10, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 10}}};
  for (auto [lo, hi] : bands) {
    std::size_t n = 0, hit = 0;
    for (const auto &sc : scores)
      if (sc.event->magnitude >= lo && sc.event->magnitude < hi) {
        ++n;
        hit += sc.detected != nullptr;
      }
    if (!n)
      continue;
    const std::string label = lo < 0   ? "M < 2"
                              : hi > 9 ? std::format("M {:.0f}+", lo)
                                       : std::format("M {:.0f}-{:.0f}", lo, hi);
    log.plain(false, "    {:<14} {:>10} {:>10} {:>8}", label, n, hit, n - hit);
  }

  if (!alarm_delays.empty()) {
    log.plain(false, "");
    log.plain(false, "  Correct alarms");
    log.plain(false, "    alarm after origin time             median {:.1f} s",
              median(alarm_delays));
    if (!magnitude_errors.empty())
      log.plain(false,
                "    magnitude error                     mean {:.2f} units",
                std::accumulate(magnitude_errors.begin(),
                                magnitude_errors.end(), 0.0) /
                    static_cast<double>(magnitude_errors.size()));
    if (!pd_errors.empty()) {
      double abs_sum = 0, sum = 0;
      for (double x : pd_errors) {
        abs_sum += std::abs(x);
        sum += x;
      }
      const auto n = static_cast<double>(pd_errors.size());
      log.plain(false,
                "    Pd magnitude error                  mean {:.2f} units, bias "
                "{:+.2f}, {} events",
                abs_sum / n, sum / n, pd_errors.size());
    }
    log.plain(false, "    located                             {} of {}",
              located -
                  static_cast<std::size_t>(std::ranges::count_if(
                      unmatched,
                      [](const Event *e) { return e->location.has_value(); })),
              correct);
    if (arrivals)
      log.plain(false,
                "    alarm before the S wave             at {} of {} station "
                "arrivals, median warning {:.1f} s, max {:.1f} s",
                warned, arrivals, median(positive_warnings),
                positive_warnings.empty()
                    ? NAN
                    : std::ranges::max(positive_warnings));
  }

  if (!missed.empty()) {
    log.plain(false, "");
    log.plain(false, "  Missed AFAD events");
    for (const auto *sc : missed)
      log.plain(false,
                "    {} UTC  {} {:.1f}  {:>4.0f} km from the network centre",
                hms(sc->event->time), sc->event->type, sc->event->magnitude,
                sc->distance_km);
  }
  if (!unmatched.empty()) {
    log.plain(false, "");
    log.plain(false, "  False alarms (no AFAD event; small uncatalogued "
                     "earthquakes are possible)");
    for (const Event *e : unmatched)
      log.plain(false, "    {} UTC  #{}  {}{}{}", hms(e->declared_at), e->id,
                e->magnitude ? std::format("M{:.1f}", *e->magnitude)
                             : std::string("no magnitude"),
                e->location ? std::format(", located {:.2f}N {:.2f}E",
                                          e->location->lat, e->location->lon)
                            : std::string(", not located"),
                cfg_.assess ? std::format(", assessed {}",
                                          verdict_name(assess(*e, t_first, t_last).verdict))
                            : std::string());
  }
  assessment_summary();
  log.plain(true, "{}", rule);
}

} // namespace ayzek::pipeline
