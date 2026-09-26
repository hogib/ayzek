#pragma once

// Alarm assessment (--assess): is a declared event an earthquake, or a
// misfire? Mostly for the alarms the catalogue does not list, since the
// catalogue misses small events, but without a catalogue it scores every
// alarm. docs/impl/16-alarm-assessment.md.
//
// The core test is S-P. A station's S-P interval gives its distance to the
// source, d = (tS - tP) vp vs / (vp - vs), a circle around it; the P and S
// times of all stations with confident picks must then fit one source. Noise
// that happens to trigger two stations together does not produce a confident
// S phase at each at a consistent S-P, and arrival times that match the
// implied distances. Independent support, when available:
//
//   geometry    the transformer's own distance (and epicentre) agrees with
//               the picker's S-P distances
//   magnitude   the station magnitudes agree with each other
//   chance      the other declaring stations were unlikely to trigger inside
//               their association windows by chance, given their trigger rates
//
// and, with a catalogue, where the event sits relative to catalogued ones:
// near a catalogued event at nearly the same time (mis-associated), or in the
// aftershock zone of an earlier one.

#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace ayzek::pipeline {

enum class Verdict { Earthquake, Possible, Misfire, Unclassified };

inline const char *verdict_name(Verdict v) {
  switch (v) {
  case Verdict::Earthquake:
    return "earthquake";
  case Verdict::Possible:
    return "possible";
  case Verdict::Misfire:
    return "misfire";
  case Verdict::Unclassified:
    break;
  }
  return "unclassified";
}

// One station with a confident P and S pick.
struct StationSP {
  std::string station;
  double sp = 0;       // S - P, seconds
  double km = 0;       // distance implied by S - P
  double p_resid = NAN, s_resid = NAN; // at the S-P source, seconds
  double geo_km = NAN; // the geometry head's distance, if any
  double geo_z = NAN;  // (log km - log geo_km) / the head's sd
  bool used = true;    // false: left out as inconsistent
};

struct Assessment {
  Verdict verdict = Verdict::Unclassified;
  std::vector<std::string> reasons; // for, against, and context, in order
  std::vector<StationSP> stations;
  std::size_t picks_rejected = 0; // picks received but not confident
  // The S-P source, if at least two stations have picks.
  double lat = NAN, lon = NAN, origin = NAN, rms = NAN;
  bool consistent = false;
  double geo_z_median = NAN;   // median |geo_z| over stations
  double geo_offset_km = NAN;  // geometry epicentre to S-P source
  double mag_spread = NAN;     // standard deviation of station magnitudes
  double p_chance = NAN;       // all other declaring stations by chance
  std::string context;         // catalogue context, if any
};

} // namespace ayzek::pipeline
