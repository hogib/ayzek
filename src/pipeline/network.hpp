#pragma once

// Network stage (one instance): associates station detections into events,
// declares an event when enough stations detect it, locates it (from the
// transformer's per-station geometry, or from P and S picks), combines station
// magnitudes, and optionally compares events with a catalogue.

#include "assess.hpp"
#include "common.hpp"
#include "locate.hpp"

#include <cmath>
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
std::vector<CatalogEvent> load_afad_catalog(const std::string &path, double t0,
                                            double t1);

// A place for which the S-wave warning time is reported, e.g. a city.
struct Site {
  std::string name;
  double lat = 0, lon = 0;
};

// Picks: grid search on P and S picks (the 60 s picker window). Geometry: the
// transformer's per-station distance and back-azimuth with its P times
// (locate.hpp), available from the trigger on.
enum class Locator { Picks, Geometry };

struct NetworkConfig {
  Locator locator = Locator::Picks;
  // Geometry locator (locate.hpp): the P-time uncertainty; a factor on every
  // station's stated distance sd (the head's is overconfident, 15-geometry-
  // location.md); and the quality gates. While a station's distance is more
  // than `geo_max_dist_z` of its (scaled) sds from the solution, or the P rms
  // exceeds `max_rms`, and more than two stations remain, the station
  // contributing most to the misfit is left out; a solution still failing
  // either, or with a 68% radius over `geo_max_err_km`, is not reported.
  double sigma_p = 0.5;
  double geo_sd_scale = 1.0;
  double geo_max_dist_z = 1e9;
  double geo_max_err_km = 1e9;
  std::size_t min_stations = 2; // detections needed to declare an event
  double slack_seconds =
      3.0; // tolerance added to the inter-station P travel time
  double coda_seconds = 40.0; // later detections at a station within this time
                              // belong to the same event
  double vp = 6.0, vs = 3.5;  // km/s, uniform half-space
  double depth_km = 10.0;     // fixed during location
  double max_rms = 2.0;       // maximum accepted location rms, seconds
  double min_pick_prob = 0.5; // minimum P and S pick probability
  std::vector<CatalogEvent>
      catalog; // reference events for evaluation (optional)
  double catalog_radius_km =
      250;                 // evaluation radius around the station centroid
  std::vector<Site> sites; // additional places for the warning-time report
  bool verbose = false;    // per-station detection, pick and magnitude lines
  // --assess (assess.hpp): an ASSESS line while running and a section in the
  // report. `context_catalog` reaches further back than `catalog`, for the
  // aftershock-zone context; empty means `catalog`.
  bool assess = false;
  std::vector<CatalogEvent> context_catalog;
  double assess_max_rms = 2.5;     // S-P fit rms, seconds
  double assess_station_tol = 3.5; // largest P or S residual at one station, s
  double assess_max_sp_km = 200;   // an S-P distance beyond this is a bad pick
  double assess_mag_spread = 0.5; // station magnitudes agree within this sd
  double assess_p_chance = 0.01;  // below this, not a chance coincidence
  double assess_geo_z = 2.0;      // geometry agrees within this many sds
};

// --assess: a pick wanted at a station that did not trigger for event
// `event_id`, P predicted at `p_pred` from the event's source, and S no later
// than `s_cap` (the station's next detected onset).
struct PickRequest {
  int event_id;
  double event_alarm;
  std::string station;
  double p_pred, s_cap;
};

struct Location {
  double lat, lon, origin, rms;
  std::size_t n_stations;
  std::vector<std::string> dropped; // stations removed during relocation
  double err_km = NAN; // geometry: radius of the 68% region; picks: not known
  double dist_z = NAN; // geometry: largest station distance residual, in sds
};

struct Event {
  int id = 0;
  bool declared = false;
  double declared_at = 0;
  std::map<std::string, Detection> detections;
  std::map<std::string, Pick> picks;
  std::size_t coda = 0; // later detections assigned to this event
  std::map<std::string, MagnitudeEstimate>
      magnitudes;                  // per station; pick-based replaces early
  std::optional<double> magnitude; // median over stations
  std::size_t magnitude_stations = 0;
  std::optional<double> first_magnitude; // first reported event magnitude
  double first_magnitude_at = 0;         // and the stream time it was reported
  std::optional<Location> location;
  double first_located_at = 0; // stream time of the first accepted location
  // Geometry locator: the latest estimate per station (common.hpp).
  std::map<std::string, StationGeometry> geometry;
  std::size_t picks_rejected = 0; // picks for this event that failed QC
  // --assess: what the last ASSESS line said, so it prints on changes only.
  int assessed_verdict = -1;
  std::size_t assessed_stations = 0;
};

// One catalogue event within the evaluation radius and the declared event
// matched to it, if any.
struct CatalogScore {
  const CatalogEvent *event;
  double distance_km;    // from the station centroid
  const Event *detected; // nullptr if missed
};

class Network {
public:
  Network(std::map<std::string, StationInfo> stations, NetworkConfig cfg);
  void on(const Detection &d);
  void on(const Pick &p);
  void on(const MagnitudeEstimate &m);
  void on(const StationGeometry &g);
  // True if the detection at `station` dated `trigger_window` belongs to a
  // declared event. Magnitude estimates at the picked P are computed only then.
  [[nodiscard]] bool declared(const std::string &station,
                              double trigger_window) const;
  // Final report: one block per declared event (magnitude, alarm, location,
  // catalogue comparison, S-wave warning time at each station and site), the
  // catalogue events not detected, and totals.
  void report(double t_first, double t_last) const;
  [[nodiscard]] std::vector<CatalogScore> score() const;
  [[nodiscard]] const std::vector<Event> &events() const noexcept {
    return events_;
  }
  [[nodiscard]] std::size_t declared_count() const;
  [[nodiscard]] std::size_t
  unmatched_count() const; // declared events with no catalogue match
  // --assess: after a replay, the picks worth making at stations that did not
  // trigger, for every declared event with a source estimate within
  // `assess_max_sp_km` of them. The caller makes them and passes them to
  // on(Pick).
  [[nodiscard]] std::vector<PickRequest> pick_requests(double t0,
                                                       double t1) const;
  // --assess: event `e` over stream times [t0, t1] (for trigger rates).
  [[nodiscard]] Assessment assess(const Event &e, double t0, double t1) const;
  // One CSV row per declared event; `record` is the --record file's name, for
  // tools/plot_candidates.py. Returns false if the file cannot be written.
  bool write_assessments(const std::string &path, double t0, double t1,
                         const std::string &record) const;

private:
  [[nodiscard]] bool compatible(const Detection &a, const Detection &b) const;
  [[nodiscard]] std::optional<Location> locate(const Event &e) const;
  // P rms over max_rms, or (geometry) a station's distance over geo_max_dist_z.
  [[nodiscard]] bool failing(const Location &l) const {
    return l.rms > cfg_.max_rms ||
           (!std::isnan(l.dist_z) && l.dist_z > cfg_.geo_max_dist_z);
  }
  [[nodiscard]] std::optional<Location>
  locate(const std::map<std::string, Pick> &picks, std::string *worst) const;
  [[nodiscard]] std::optional<Location>
  locate(const std::map<std::string, StationGeometry> &geo,
         std::string *worst) const;
  [[nodiscard]] const CatalogEvent *match(const Event &e) const;
  void report_location(Event &e, double now);
  void report_magnitude(Event &e, double now);

  struct Warning {
    std::string name;
    double distance_km;
    double s_time;
    bool picked; // S time from this event's S pick, otherwise predicted
    bool site;
  };
  // S arrival at each station and site for event `e`, using the catalogue
  // hypocentre if `c` is given, otherwise ayzek's location. Empty if neither
  // is available. Sorted by distance.
  [[nodiscard]] std::vector<Warning> warnings(const Event &e,
                                              const CatalogEvent *c) const;

  std::map<std::string, StationInfo> stations_;
  NetworkConfig cfg_;
  std::vector<Event> events_;
  int next_id_ = 1;
  // Every detection per station, coda included, for trigger rates.
  std::map<std::string, std::size_t> triggers_;
  double first_seen_ = NAN; // window start of the first detection
  void print_assessment(Event &e, double now);
};

double distance_km(double lat1, double lon1, double lat2, double lon2);

} // namespace ayzek::pipeline
