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

// A place for which the S-wave warning time is reported, e.g. a city.
struct Site {
    std::string name;
    double lat = 0, lon = 0;
};

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
    std::vector<Site> sites;           // additional places for the warning-time report
    bool verbose = false;              // per-station detection, pick and magnitude lines
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
    std::optional<double> first_magnitude;   // first reported event magnitude
    double first_magnitude_at = 0;           // and the stream time it was reported
    std::optional<Location> location;
};

// One catalogue event within the evaluation radius and the declared event
// matched to it, if any.
struct CatalogScore {
    const CatalogEvent* event;
    double distance_km;                 // from the station centroid
    const Event* detected;              // nullptr if missed
};

class Network {
public:
    Network(std::map<std::string, StationInfo> stations, NetworkConfig cfg);
    void on(const Detection& d);
    void on(const Pick& p);
    void on(const MagnitudeEstimate& m);
    // True if the detection at `station` dated `trigger_window` belongs to a
    // declared event. Magnitude estimates at the picked P are computed only then.
    [[nodiscard]] bool declared(const std::string& station, double trigger_window) const;
    // Final report: one block per declared event (magnitude, alarm, location,
    // catalogue comparison, S-wave warning time at each station and site), the
    // catalogue events not detected, and totals.
    void report(double t_first, double t_last) const;
    [[nodiscard]] std::vector<CatalogScore> score() const;
    [[nodiscard]] const std::vector<Event>& events() const noexcept { return events_; }
    [[nodiscard]] std::size_t declared_count() const;
    [[nodiscard]] std::size_t unmatched_count() const;   // declared events with no catalogue match

private:
    [[nodiscard]] bool compatible(const Detection& a, const Detection& b) const;
    [[nodiscard]] std::optional<Location> locate(const Event& e) const;
    [[nodiscard]] std::optional<Location> locate(const std::map<std::string, Pick>& picks, std::string* worst) const;
    [[nodiscard]] const CatalogEvent* match(const Event& e) const;
    void report_location(Event& e);
    void report_magnitude(Event& e, double now);

    struct Warning {
        std::string name;
        double distance_km;
        double s_time;
        bool picked;          // S time from this event's S pick, otherwise predicted
        bool site;
    };
    // S arrival at each station and site for event `e`, using the catalogue
    // hypocentre if `c` is given, otherwise ayzek's location. Empty if neither
    // is available. Sorted by distance.
    [[nodiscard]] std::vector<Warning> warnings(const Event& e, const CatalogEvent* c) const;

    std::map<std::string, StationInfo> stations_;
    NetworkConfig cfg_;
    std::vector<Event> events_;
    int next_id_ = 1;
};

double distance_km(double lat1, double lon1, double lat2, double lon2);

}  // namespace ayzek::pipeline
