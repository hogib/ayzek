#pragma once

// Recording of the station processors' outputs (detections, picks, magnitude
// estimates) in the order the network stage processed them.
//
// Station processing is independent per station; only the network stage
// combines stations. A recording of all stations therefore allows the network
// stage to be re-run for any subset of stations without reprocessing the
// waveforms, with the same result as running ayzek on that subset.
//
// Format: tab-separated text, one record per line. Numbers are written in the
// shortest form that reads back to the same value.
//   S  code  lat  lon
//   T  t_first  t_last
//   D  station  window_start  declared_at  probability  compute_ms
//   P  station  trigger_window  p_time  s_time  p_prob  s_prob  declared_at  compute_ms
//   M  station  trigger_window  at_pick  window_start  magnitude  noise_windows  declared_at  compute_ms

#include "common.hpp"
#include "network.hpp"

#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace ayzek::pipeline {

class Recorder {
public:
    Recorder(const std::string& path, const std::map<std::string, StationInfo>& stations, double t_first, double t_last);
    void write(const Message& m);

private:
    std::ofstream out_;
};

struct Recording {
    std::map<std::string, StationInfo> stations;
    double t_first = 0, t_last = 0;
    std::vector<Message> messages;       // Detection, Pick, MagnitudeEstimate

    static Recording load(const std::string& path);
};

// Station code of a Detection, Pick or MagnitudeEstimate.
const std::string& station_of(const Message& m);

}  // namespace ayzek::pipeline
