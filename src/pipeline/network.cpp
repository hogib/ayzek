#include "network.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <numeric>
#include <numbers>
#include <ranges>
#include <sstream>
#include <stdexcept>

namespace ayzek::pipeline {

double distance_km(double lat1, double lon1, double lat2, double lon2) {
    const double r = std::numbers::pi / 180.0;
    const double a = std::pow(std::sin((lat2 - lat1) * r / 2), 2) +
                     std::cos(lat1 * r) * std::cos(lat2 * r) * std::pow(std::sin((lon2 - lon1) * r / 2), 2);
    return 2 * 6371.0 * std::asin(std::sqrt(a));
}

std::vector<CatalogEvent> load_afad_catalog(const std::string& path, double t0, double t1) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
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
        if (!parsed) continue;
        const double epoch = *parsed;
        if (epoch < t0 || epoch > t1) continue;
        try {
            out.push_back({epoch, std::stod(lat), std::stod(lon), std::stod(depth), std::stod(mag), type});
        } catch (const std::exception&) {
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
bool Network::compatible(const Detection& a, const Detection& b) const {
    auto ia = stations_.find(a.station), ib = stations_.find(b.station);
    if (ia == stations_.end() || ib == stations_.end()) return std::abs(a.window_start - b.window_start) <= 30.0;
    const double km = distance_km(ia->second.lat, ia->second.lon, ib->second.lat, ib->second.lon);
    return std::abs(a.window_start - b.window_start) <= km / cfg_.vp + cfg_.slack_seconds;
}

void Network::on(const Detection& d) {
    // The detector can trigger again on the S wave and coda of an event. A later
    // detection at a station within `coda_seconds` of that station's detection
    // of a declared event is assigned to that event. A separate event in that
    // interval is therefore not declared from this station.
    for (auto& e : events_) {
        auto it = e.detections.find(d.station);
        if (e.declared && it != e.detections.end() && d.window_start > it->second.window_start &&
            d.window_start - it->second.window_start <= cfg_.coda_seconds) {
            ++e.coda;
            if (cfg_.verbose)
                Log::get().line("detect", "2", "{:<5} p={:.2f}  window {}  -- coda of #{}", d.station, d.probability,
                                hms(d.window_start), e.id);
            return;
        }
    }
    if (cfg_.verbose)
        Log::get().line("detect", "36", "{:<5} p={:.2f}  window {}  ({:.1f} ms)", d.station, d.probability,
                        hms(d.window_start), d.compute_ms);

    Event* ev = nullptr;
    for (auto& e : events_) {
        if (e.detections.contains(d.station)) continue;
        if (std::ranges::all_of(e.detections | std::views::values, [&](const Detection& o) { return compatible(o, d); })) {
            ev = &e;
            break;
        }
    }
    if (!ev) {
        events_.push_back(Event{});
        ev = &events_.back();
    }
    ev->detections.emplace(d.station, d);

    if (!ev->declared && ev->detections.size() >= cfg_.min_stations) {
        ev->declared = true;
        ev->id = next_id_++;
        ev->declared_at = std::ranges::max(ev->detections | std::views::values, {}, &Detection::declared_at).declared_at;
        std::string names;
        for (const auto& [code, _] : ev->detections) names += (names.empty() ? "" : ", ") + code;
        const auto* c = match(*ev);
        Log::get().line("ALARM", "1;31", "#{:<3} {}  earthquake detected by {}{}", ev->id, hms(ev->declared_at), names,
                        c ? std::format("  (AFAD {} {:.1f} at {}, {:.1f} s ago)", c->type, c->magnitude, hms(c->time),
                                        ev->declared_at - c->time)
                          : "");
        report_magnitude(*ev, ev->declared_at);    // estimates received before the event was declared
    } else if (ev->declared && cfg_.verbose) {
        Log::get().line("EVENT", "31", "#{} joined by {}", ev->id, d.station);
    }
}

void Network::on(const Pick& p) {
    const double sp = p.s_time - p.p_time;
    const double km = sp * cfg_.vp * cfg_.vs / (cfg_.vp - cfg_.vs);   // S-P lag to distance

    // Pick quality control. The picker returns the most probable chunk even when
    // there is no phase, so low-probability picks, S before P, and P far from the
    // trigger window are not used.
    std::string reject;
    if (p.p_prob < cfg_.min_pick_prob || p.s_prob < cfg_.min_pick_prob) reject = "low confidence";
    else if (sp <= 0 || sp > 60) reject = "S not after P";
    else if (p.p_time < p.trigger_window - 3 || p.p_time > p.trigger_window + 9) reject = "P outside the trigger window";
    if (!reject.empty()) {
        if (cfg_.verbose) Log::get().line("pick", "2", "{:<5} P {} ({:.2f})  S {} ({:.2f})  -- not used: {}", p.station, hms(p.p_time), p.p_prob,
                        hms(p.s_time), p.s_prob, reject);
        return;
    }
    if (cfg_.verbose)
        Log::get().line("pick", "35", "{:<5} P {} ({:.2f})  S {} ({:.2f})  S-P {:.2f} s ~ {:.0f} km  ({:.0f} ms)", p.station,
                        hms(p.p_time), p.p_prob, hms(p.s_time), p.s_prob, sp, km, p.compute_ms);
    for (auto& e : events_) {
        auto it = e.detections.find(p.station);
        if (it == e.detections.end() || it->second.window_start != p.trigger_window) continue;
        e.picks.emplace(p.station, p);
        if (e.declared) report_location(e);
        return;
    }
}

// Locates the event; while rms exceeds max_rms and more than two stations
// remain, removes the station with the largest residual and relocates. The
// uniform velocity model underestimates P speed beyond ~150 km, where the first
// arrival travels through the upper mantle, so distant stations are the ones
// typically removed.
void Network::on(const MagnitudeEstimate& m) {
    if (cfg_.verbose) Log::get().line("mag", "34", "{:<5} M{:.2f}  {}  noise baseline {}  ({:.0f} ms)", m.station, m.magnitude,
                    m.at_pick ? "at picked P" : "early      ",
                    m.noise_windows ? std::format("{} windows", m.noise_windows) : std::string("not ready, per-window"),
                    m.compute_ms);
    for (auto& e : events_) {
        auto it = e.detections.find(m.station);
        if (it == e.detections.end() || it->second.window_start != m.trigger_window) continue;
        auto have = e.magnitudes.find(m.station);
        if (have != e.magnitudes.end() && have->second.at_pick && !m.at_pick) return;
        e.magnitudes.insert_or_assign(m.station, m);
        if (e.declared) report_magnitude(e, m.declared_at);
        return;
    }
}

// Event magnitude: median of the station estimates. Printed when the value
// changes by 0.05 or more, or when the number of stations changes.
void Network::report_magnitude(Event& e, double now) {
    if (e.magnitudes.empty()) return;
    std::vector<double> v;
    std::string parts;
    for (const auto& [code, sm] : e.magnitudes) {
        v.push_back(sm.magnitude);
        parts += std::format("{}{} {:.1f}", parts.empty() ? "" : ", ", code, sm.magnitude);
    }
    std::ranges::sort(v);
    const double med = v.size() % 2 ? v[v.size() / 2] : 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]);
    const bool changed = !e.magnitude || std::abs(*e.magnitude - med) >= 0.05 || e.magnitudes.size() != e.magnitude_stations;
    e.magnitude = med;
    e.magnitude_stations = e.magnitudes.size();
    if (!e.first_magnitude) {
        e.first_magnitude = med;
        e.first_magnitude_at = now;
    }
    if (!changed) return;
    Log::get().line("MAG", "1;34", "#{:<3} {}  magnitude M{:.1f} from {} station{} ({})", e.id, hms(now), med, v.size(),
                    v.size() == 1 ? "" : "s", parts);
}

std::optional<Location> Network::locate(const Event& e) const {
    auto picks = e.picks;
    std::vector<std::string> dropped;
    for (;;) {
        std::string worst;
        auto loc = locate(picks, &worst);
        if (!loc || loc->rms <= cfg_.max_rms || picks.size() <= 2) {
            if (loc) loc->dropped = dropped;
            return loc;
        }
        dropped.push_back(worst);
        picks.erase(worst);
    }
}

std::optional<Location> Network::locate(const std::map<std::string, Pick>& picks, std::string* worst) const {
    struct Obs {
        double lat, lon, t, v;
        const std::string* station;
    };
    std::vector<Obs> obs;
    double clat = 0, clon = 0;
    for (const auto& [code, pk] : picks) {
        auto it = stations_.find(code);
        if (it == stations_.end()) continue;
        obs.push_back({it->second.lat, it->second.lon, pk.p_time, cfg_.vp, &code});
        obs.push_back({it->second.lat, it->second.lon, pk.s_time, cfg_.vs, &code});
        clat += it->second.lat;
        clon += it->second.lon;
    }
    const std::size_t n_st = obs.size() / 2;
    if (n_st < 2) return std::nullopt;        // 2 stations x (P, S) = 4 times for 3 unknowns
    clat /= static_cast<double>(n_st);
    clon /= static_cast<double>(n_st);

    // Grid search over epicentre; depth fixed. For a given epicentre the least-
    // squares origin time is the mean of (observed - predicted travel time).
    std::vector<double> r(obs.size());
    auto misfit = [&](double lat, double lon, double& origin) {
        double sum = 0, sq = 0;
        for (std::size_t i = 0; i < obs.size(); ++i) {
            r[i] = obs[i].t - std::hypot(distance_km(lat, lon, obs[i].lat, obs[i].lon), cfg_.depth_km) / obs[i].v;
            sum += r[i];
        }
        origin = sum / static_cast<double>(r.size());
        for (double x : r) sq += (x - origin) * (x - origin);
        return std::sqrt(sq / static_cast<double>(r.size()));
    };
    Location best{clat, clon, 0, 1e18, n_st, {}};
    auto search = [&](double lat0, double lon0, double half, double step) {
        for (double la = lat0 - half; la <= lat0 + half; la += step)
            for (double lo = lon0 - half; lo <= lon0 + half; lo += step) {
                double origin = 0;
                const double m = misfit(la, lo, origin);
                if (m < best.rms) best = {la, lo, origin, m, n_st, {}};
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

void Network::report_location(Event& e) {
    auto loc = locate(e);
    if (!loc) return;
    if (loc->rms > cfg_.max_rms) {
        if (cfg_.verbose)
            Log::get().line("LOCATE", "33", "#{} rejected: best fit rms {:.1f} s from {} stations", e.id, loc->rms, loc->n_stations);
        return;
    }
    // Do not print an unchanged solution.
    const bool same = e.location && e.location->lat == loc->lat && e.location->lon == loc->lon &&
                      e.location->origin == loc->origin;
    e.location = loc;
    if (same) return;
    std::string left_out;
    for (const auto& d : loc->dropped) left_out += (left_out.empty() ? ", left out " : " ") + d;
    const double now = std::ranges::max(e.picks | std::views::values, {}, &Pick::declared_at).declared_at;
    Log::get().line("LOCATE", "1;32", "#{:<3} {}  located at {:.3f}N {:.3f}E, origin {}, rms {:.2f} s from {} stations{}", e.id,
                    hms(now), loc->lat, loc->lon, hms(loc->origin), loc->rms, loc->n_stations, left_out);
}

// Matching catalogue event, if any. A catalogue event matches if the location
// agrees (origin within 5 s, epicentre within 30 km), or if its predicted P
// arrival at a detecting station is within 6 s of that detection's window start
// plus 3.5 s. The second criterion does not require a location.
const CatalogEvent* Network::match(const Event& e) const {
    const CatalogEvent* best = nullptr;
    double best_miss = 1e18;
    for (const auto& c : cfg_.catalog) {
        double miss = 1e18;
        if (e.location && distance_km(e.location->lat, e.location->lon, c.lat, c.lon) <= 30.0 &&
            std::abs(e.location->origin - c.time) <= 5.0)
            miss = std::abs(e.location->origin - c.time) / 100.0;   // location matches rank first
        for (const auto& [code, d] : e.detections) {
            auto it = stations_.find(code);
            if (it == stations_.end()) continue;
            const double p = c.time + std::hypot(distance_km(c.lat, c.lon, it->second.lat, it->second.lon), c.depth) / cfg_.vp;
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
    if (cfg_.catalog.empty() || stations_.empty()) return out;
    double clat = 0, clon = 0;
    for (const auto& [_, s] : stations_) {
        clat += s.lat;
        clon += s.lon;
    }
    clat /= static_cast<double>(stations_.size());
    clon /= static_cast<double>(stations_.size());
    for (const auto& c : cfg_.catalog) {
        const double dist = distance_km(clat, clon, c.lat, c.lon);
        if (dist > cfg_.catalog_radius_km) continue;
        const Event* hit = nullptr;
        for (const auto& e : events_)
            if (e.declared && match(e) == &c) {
                hit = &e;
                break;
            }
        out.push_back({&c, dist, hit});
    }
    return out;
}

std::size_t Network::declared_count() const {
    return static_cast<std::size_t>(std::ranges::count_if(events_, &Event::declared));
}

std::size_t Network::unmatched_count() const {
    const auto scores = score();
    std::size_t n = 0;
    for (const auto& e : events_)
        if (e.declared && std::ranges::none_of(scores, [&](const CatalogScore& s) { return s.detected == &e; })) ++n;
    return n;
}

std::vector<Network::Warning> Network::warnings(const Event& e, const CatalogEvent* c) const {
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
    auto add = [&](const std::string& name, double slat, double slon, bool site) {
        const double km = distance_km(lat, lon, slat, slon);
        Warning w{name, km, origin + std::hypot(km, depth) / cfg_.vs, false, site};
        if (!site)
            if (auto it = e.picks.find(name); it != e.picks.end()) {
                w.s_time = it->second.s_time;
                w.picked = true;
            }
        out.push_back(w);
    };
    for (const auto& [code, s] : stations_) add(code, s.lat, s.lon, false);
    for (const auto& site : cfg_.sites) add(site.name, site.lat, site.lon, true);
    std::ranges::sort(out, {}, &Warning::distance_km);
    return out;
}

void Network::report(double t_first, double t_last) const {
    auto& log = Log::get();
    const std::string rule(78, '-');
    log.plain(false, "");
    log.plain(true, "{}", rule);
    log.plain(true, "REPORT  {} stations, {} to {} UTC", stations_.size(), ymd_hms(t_first), hms(t_last));
    log.plain(true, "{}", rule);

    const auto scores = score();
    std::vector<const Event*> declared;
    for (const auto& e : events_)
        if (e.declared) declared.push_back(&e);
    std::ranges::sort(declared, {}, &Event::declared_at);

    std::vector<double> alarm_delays, magnitude_errors, positive_warnings;
    std::size_t arrivals = 0, warned = 0;
    for (const Event* e : declared) {
        const CatalogEvent* c = nullptr;
        for (const auto& sc : scores)
            if (sc.detected == e) c = sc.event;

        // Headline.
        std::string headline = std::format("Event #{}: detected ", e->id);
        headline += e->magnitude ? std::format("M{:.1f} earthquake", *e->magnitude) : std::string("earthquake (no magnitude)");
        if (e->location)
            headline += std::format(" at {:.3f}N {:.3f}E, origin {} UTC", e->location->lat, e->location->lon, hms(e->location->origin));
        else
            headline += ", not located";
        log.plain(false, "");
        log.plain(true, "{}", headline);

        std::string first, later;
        for (const auto& [code, d] : e->detections)
            (d.declared_at <= e->declared_at ? first : later) += std::format("{}{}", (d.declared_at <= e->declared_at ? first : later).empty() ? "" : ", ", code);
        log.plain(false, "  alarm      {} UTC by {}{}", hms(e->declared_at), first, later.empty() ? "" : std::format("; later {}", later));
        if (e->magnitude) {
            const double lag = e->first_magnitude_at - e->declared_at;
            log.plain(false, "  magnitude  M{:.1f} {}, M{:.1f} final from {} station{}", *e->first_magnitude,
                      lag < 0.05 ? std::string("at the alarm") : std::format("{:.1f} s after the alarm", lag),
                      *e->magnitude, e->magnitude_stations, e->magnitude_stations == 1 ? "" : "s");
        }
        if (e->location)
            log.plain(false, "  location   rms {:.2f} s from {} stations", e->location->rms, e->location->n_stations);

        const double origin = c ? c->time : (e->location ? e->location->origin : NAN);
        if (c) {
            std::string cmp = std::format("  AFAD       {} {:.1f} at {} UTC: alarm {:.1f} s after origin", c->type, c->magnitude,
                                          hms(c->time), e->declared_at - c->time);
            // Rounded before formatting, so that a difference below 0.05 prints as +0.0.
            if (e->magnitude) cmp += std::format(", magnitude {:+.1f}", std::round((*e->magnitude - c->magnitude) * 10) / 10 + 0.0);
            if (e->location)
                cmp += std::format(", epicentre {:.1f} km off, origin {:+.1f} s",
                                   distance_km(e->location->lat, e->location->lon, c->lat, c->lon), e->location->origin - c->time);
            log.plain(false, "{}", cmp);
            alarm_delays.push_back(e->declared_at - c->time);
            if (e->magnitude) magnitude_errors.push_back(std::abs(*e->magnitude - c->magnitude));
        } else {
            log.plain(false, "  AFAD       no matching catalogue event");
        }

        const auto ws = warnings(*e, c);
        if (ws.empty()) continue;
        const double depth = c ? c->depth : cfg_.depth_km;
        const double blind = std::sqrt(std::max(cfg_.vs * cfg_.vs * (e->declared_at - origin) * (e->declared_at - origin) - depth * depth, 0.0));
        log.plain(false, "  S wave     had travelled {:.0f} km from the epicentre when the alarm sounded ({})", blind,
                  c ? "AFAD hypocentre" : "ayzek location");
        log.plain(false, "  warning    {:<10} {:>8}   {:<24} {:>8}", "", "distance", "S wave arrives", "warning");
        for (const auto& w : ws) {
            const double lead = w.s_time - e->declared_at;
            log.plain(false, "             {:<10} {:>5.0f} km   {} UTC {:<10} {:>+7.1f} s  {}", w.name, w.distance_km, hms(w.s_time),
                      w.picked ? "(picked)" : "(predicted)", lead, lead > 0 ? "before S" : "after S");
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
        if (v.empty()) return NAN;
        std::ranges::sort(v);
        return v.size() % 2 ? v[v.size() / 2] : 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]);
    };
    auto pct = [](std::size_t a, std::size_t b) { return b ? std::format("{:.0f}%", 100.0 * static_cast<double>(a) / static_cast<double>(b)) : std::string("-"); };

    std::size_t located = 0;
    std::vector<const Event*> unmatched;
    for (const Event* e : declared) {
        located += e->location.has_value();
        if (std::ranges::none_of(scores, [&](const CatalogScore& sc) { return sc.detected == e; })) unmatched.push_back(e);
    }
    std::vector<const CatalogScore*> missed;
    for (const auto& sc : scores)
        if (!sc.detected) missed.push_back(&sc);

    log.plain(false, "");
    log.plain(true, "{}", rule);
    log.plain(true, "SUMMARY");
    log.plain(true, "{}", rule);
    log.plain(false, "  Alarms raised                        {:>4}", declared.size());
    log.plain(false, "  single-station detections, no alarm  {:>4}", events_.size() - declared.size());
    if (cfg_.catalog.empty()) {
        log.plain(false, "  (no catalogue given: correct, false and missed alarms not evaluated)");
        log.plain(true, "{}", rule);
        return;
    }
    const std::size_t correct = declared.size() - unmatched.size();
    log.plain(false, "");
    log.plain(false, "  Alarms compared with the AFAD catalogue");
    log.plain(false, "    correct  (match an AFAD event)     {:>4}  {:>4}", correct, pct(correct, declared.size()));
    log.plain(false, "    false    (no AFAD event)           {:>4}  {:>4}", unmatched.size(), pct(unmatched.size(), declared.size()));
    log.plain(false, "");
    log.plain(false, "  AFAD events within {:.0f} km             {:>4}", cfg_.catalog_radius_km, scores.size());
    log.plain(false, "    detected                           {:>4}  {:>4}", scores.size() - missed.size(), pct(scores.size() - missed.size(), scores.size()));
    log.plain(false, "    missed                             {:>4}  {:>4}", missed.size(), pct(missed.size(), scores.size()));

    // Detection rate by magnitude.
    log.plain(false, "");
    log.plain(false, "  By magnitude      AFAD events   detected   missed");
    const std::array<std::pair<double, double>, 5> bands{{{-10, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 10}}};
    for (auto [lo, hi] : bands) {
        std::size_t n = 0, hit = 0;
        for (const auto& sc : scores)
            if (sc.event->magnitude >= lo && sc.event->magnitude < hi) {
                ++n;
                hit += sc.detected != nullptr;
            }
        if (!n) continue;
        const std::string label = lo < 0 ? "M < 2" : hi > 9 ? std::format("M {:.0f}+", lo) : std::format("M {:.0f}-{:.0f}", lo, hi);
        log.plain(false, "    {:<14} {:>10} {:>10} {:>8}", label, n, hit, n - hit);
    }

    if (!alarm_delays.empty()) {
        log.plain(false, "");
        log.plain(false, "  Correct alarms");
        log.plain(false, "    alarm after origin time             median {:.1f} s", median(alarm_delays));
        if (!magnitude_errors.empty())
            log.plain(false, "    magnitude error                     mean {:.2f} units",
                      std::accumulate(magnitude_errors.begin(), magnitude_errors.end(), 0.0) / static_cast<double>(magnitude_errors.size()));
        log.plain(false, "    located                             {} of {}", located - static_cast<std::size_t>(std::ranges::count_if(unmatched, [](const Event* e) { return e->location.has_value(); })), correct);
        if (arrivals)
            log.plain(false, "    alarm before the S wave             at {} of {} station arrivals, median warning {:.1f} s, max {:.1f} s",
                      warned, arrivals, median(positive_warnings), positive_warnings.empty() ? NAN : std::ranges::max(positive_warnings));
    }

    if (!missed.empty()) {
        log.plain(false, "");
        log.plain(false, "  Missed AFAD events");
        for (const auto* sc : missed)
            log.plain(false, "    {} UTC  {} {:.1f}  {:>4.0f} km from the network centre", hms(sc->event->time), sc->event->type,
                      sc->event->magnitude, sc->distance_km);
    }
    if (!unmatched.empty()) {
        log.plain(false, "");
        log.plain(false, "  False alarms (no AFAD event; small uncatalogued earthquakes are possible)");
        for (const Event* e : unmatched)
            log.plain(false, "    {} UTC  #{}  {}{}", hms(e->declared_at), e->id,
                      e->magnitude ? std::format("M{:.1f}", *e->magnitude) : std::string("no magnitude"),
                      e->location ? std::format(", located {:.2f}N {:.2f}E", e->location->lat, e->location->lon) : std::string(", not located"));
    }
    log.plain(true, "{}", rule);
}

}  // namespace ayzek::pipeline
