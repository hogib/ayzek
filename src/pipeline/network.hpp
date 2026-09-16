#pragma once

// Network stage (one instance): associates station detections into events,
// declares an event when enough stations detect it, locates it from the picks,
// combines station magnitudes, and optionally compares events with a catalogue.

#include "common.hpp"

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace ayzek::pipeline {

struct StationInfo {
    double lat = 0, lon = 0;
};

struct CatalogEvent {
    double time, lat, lon, depth, magnitude;
    std::string type;
};

// AFAD's CSV export: Date (DD/MM/YYYY HH:MM:SS, UTC), Longitude, Latitude,
// Depth, Rms, Type, Magnitude, ... Only events in [t0, t1] are kept.
std::vector<CatalogEvent> load_afad_catalog(const std::string& path, double t0, double t1);

struct NetworkConfig {
    std::size_t min_stations = 2;      // detections needed to declare an event
    double slack_seconds = 3.0;        // tolerance added to the inter-station P travel time
    double coda_seconds = 40.0;        // later detections at a station within this time belong to the same event
    double vp = 6.0, vs = 3.5;         // km/s, uniform half-space
    double depth_km = 10.0;            // fixed during location
    double max_rms = 2.0;              // maximum accepted location rms, seconds
    double min_pick_prob = 0.5;        // minimum P and S pick probability
    std::vector<CatalogEvent> catalog; // reference events for evaluation (optional)
    double catalog_radius_km = 250;    // evaluation radius around the station centroid
};

struct Location {
    double lat, lon, origin, rms;
    std::size_t n_stations;
    std::vector<std::string> dropped;   // stations removed during relocation
};

struct Event {
    int id = 0;
    bool declared = false;
    double declared_at = 0;
    std::map<std::string, Detection> detections;
    std::map<std::string, Pick> picks;
    std::size_t coda = 0;               // later detections assigned to this event
    std::map<std::string, MagnitudeEstimate> magnitudes;   // per station; pick-based replaces early
    std::optional<double> magnitude;    // median over stations
    std::size_t magnitude_stations = 0;
    std::optional<Location> location;
};

class Network {
public:
    Network(std::map<std::string, StationInfo> stations, NetworkConfig cfg);
    void on(const Detection& d);
    void on(const Pick& p);
    void on(const MagnitudeEstimate& m);
    void summary() const;

private:
    [[nodiscard]] bool compatible(const Detection& a, const Detection& b) const;
    [[nodiscard]] std::optional<Location> locate(const Event& e) const;
    [[nodiscard]] std::optional<Location> locate(const std::map<std::string, Pick>& picks, std::string* worst) const;
    [[nodiscard]] const CatalogEvent* match(const Event& e) const;
    void report_location(Event& e);
    void report_magnitude(Event& e, double now);

    std::map<std::string, StationInfo> stations_;
    NetworkConfig cfg_;
    std::vector<Event> events_;
    int next_id_ = 1;
};

double distance_km(double lat1, double lon1, double lat2, double lon2);

}  // namespace ayzek::pipeline
