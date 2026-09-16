#include "network.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
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
            Log::get().line("detect", "2", "{:<5} p={:.2f}  window {}  -- coda of {}", d.station, d.probability,
                            hms(d.window_start), std::format("#{}", e.id));
            return;
        }
    }
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
        Log::get().line("EVENT", "1;31", "#{} declared at {} by {}", ev->id, hms(ev->declared_at), names);
        if (const auto* c = match(*ev)) {
            Log::get().line("EVENT", "31", "#{} catalogue: {}{:.1f} at {}, alert {:.1f} s after origin", ev->id, c->type,
                            c->magnitude, hms(c->time), ev->declared_at - c->time);
        }
        report_magnitude(*ev, ev->declared_at);    // estimates received before the event was declared
    } else if (ev->declared) {
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
        Log::get().line("pick", "2", "{:<5} P {} ({:.2f})  S {} ({:.2f})  -- not used: {}", p.station, hms(p.p_time), p.p_prob,
                        hms(p.s_time), p.s_prob, reject);
        return;
    }
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
    Log::get().line("mag", "34", "{:<5} M{:.2f}  {}  noise baseline {}  ({:.0f} ms)", m.station, m.magnitude,
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
    if (!changed) return;
    const auto* c = match(e);
    Log::get().line("MAG", "1;34", "#{} M{:.1f} at {} from {} station{} ({}){}", e.id, med, hms(now), v.size(),
                    v.size() == 1 ? "" : "s", parts,
                    c ? std::format("  catalogue {}{:.1f}: {:+.1f}, {:.0f} s after origin", c->type, c->magnitude,
                                    med - c->magnitude, now - c->time)
                      : "");
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
    Log::get().line("LOCATE", "1;32", "#{} {:.3f}N {:.3f}E  origin {}  rms {:.2f} s  ({} stations, P+S{})", e.id, loc->lat,
                    loc->lon, hms(loc->origin), loc->rms, loc->n_stations, left_out);
    if (const auto* c = match(e)) {
        Log::get().line("LOCATE", "32", "#{} catalogue {}{:.1f} {:.3f}N {:.3f}E: {:.1f} km off, origin {:+.1f} s", e.id,
                        c->type, c->magnitude, c->lat, c->lon, distance_km(loc->lat, loc->lon, c->lat, c->lon),
                        loc->origin - c->time);
    }
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

void Network::summary() const {
    std::size_t declared = 0;
    for (const auto& e : events_) declared += e.declared;
    Log::get().line("summary", "1", "{} events declared, {} single-station detections not confirmed", declared,
                    events_.size() - declared);
    if (cfg_.catalog.empty()) return;

    double clat = 0, clon = 0;
    for (const auto& [_, s] : stations_) { clat += s.lat; clon += s.lon; }
    clat /= static_cast<double>(stations_.size());
    clon /= static_cast<double>(stations_.size());

    Log::get().line("catalog", "1", "{:<11} {:>6} {:>7}   ayzek", "origin", "mag", "dist");
    std::size_t n = 0, found = 0, located = 0;
    std::vector<const Event*> matched;
    for (const auto& c : cfg_.catalog) {
        const double dist = distance_km(clat, clon, c.lat, c.lon);
        if (dist > cfg_.catalog_radius_km) continue;
        ++n;
        const Event* hit = nullptr;
        for (const auto& e : events_)
            if (e.declared && match(e) == &c) { hit = &e; break; }
        std::string what = "missed";
        if (hit) {
            ++found;
            matched.push_back(hit);
            what = std::format("#{} alert +{:.1f} s", hit->id, hit->declared_at - c.time);
            if (hit->magnitude) what += std::format(", M{:.1f} ({:+.1f})", *hit->magnitude, *hit->magnitude - c.magnitude);
            if (hit->location) {
                ++located;
                what += std::format(", located {:.1f} km off, origin {:+.1f} s", distance_km(hit->location->lat, hit->location->lon, c.lat, c.lon),
                                    hit->location->origin - c.time);
            }
        }
        Log::get().line("catalog", hit ? "32" : "33", "{:<11} {:>2}{:>4.1f} {:>4.0f} km   {}", hms(c.time), c.type, c.magnitude, dist, what);
    }
    std::size_t unmatched = 0;
    for (const auto& e : events_)
        if (e.declared && std::ranges::find(matched, &e) == matched.end()) ++unmatched;
    Log::get().line("catalog", "1", "{} of {} catalogue events within {:.0f} km declared, {} located; {} declared events not in the catalogue",
                    found, n, cfg_.catalog_radius_km, located, unmatched);
}

}  // namespace ayzek::pipeline
