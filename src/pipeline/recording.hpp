#pragma once

// Recording of the station processors' outputs (detections, picks, geometry
// estimates, Pd values) in the order the network stage processed them.
//
// Station processing is independent per station; only the network stage
// combines stations. A recording of all stations therefore allows the network
// stage to be re-run for any subset of stations without reprocessing the
// waveforms, with the same result as running ayzek on that subset.
//
// Format: tab-separated text, one record per line. Numbers are written in the
// shortest form that reads back to the same value. Lines starting with '#'
// are comments, except `# min_stations N`: the stations the run needed to
// declare an event (absent in older recordings), the default for a replay.
//   S  code  lat  lon
//   T  t_first  t_last
//   D  station  window_start  declared_at  probability  compute_ms  [restart
//      [p_time]]  (restart 1 for a transformer dt restart; p_time the
//      detection's onset, Detection::p_time; both absent in older recordings)
//   P  station  trigger_window  p_time  s_time  p_prob  s_prob  declared_at
//      compute_ms  [event_id  event_alarm]  (--assess picks at stations that
//      did not trigger; absent otherwise)
//   M  a magnitude-regressor estimate, in recordings made before it was
//      retired (v0.1.0-beta.6 and earlier); skipped
//   G  station  trigger_window  p_time  since_p  log_dist  log_dist_sd
//      declared_at  (older recordings have baz and kappa before declared_at;
//      they are skipped)
//   Q  station  trigger_window  p_time  tau  pd  pd_noise  declared_at
// PdTrace messages (two per second per station) are not recorded, so a replay
// of a recording estimates Pd without the bounds of silent stations.

#include "common.hpp"
#include "network.hpp"

#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace ayzek::pipeline {

class Recorder {
public:
  Recorder(const std::string &path,
           const std::map<std::string, StationInfo> &stations, double t_first,
           double t_last, std::size_t min_stations);
  void write(const Message &m);

private:
  std::ofstream out_;
};

struct Recording {
  std::map<std::string, StationInfo> stations;
  double t_first = 0, t_last = 0;
  std::size_t min_stations = 0; // 0: not recorded
  std::vector<Message> messages; // Detection, Pick, StationGeometry,
                                 // PdEstimate

  static Recording load(const std::string &path);
};

// Station code of a Detection, Pick, StationGeometry or PdEstimate.
const std::string &station_of(const Message &m);

} // namespace ayzek::pipeline
