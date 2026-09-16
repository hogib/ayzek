#pragma once

// Ingest stage, one thread per station: miniSEED records -> Steim2 decoding ->
// reordering -> component rings.
//
// ReplaySource reads records from a file. Each record is released when the
// stream clock reaches the time of its last sample, the earliest time a
// datalogger could transmit it. A real-time client (e.g. SeedLink) would take
// the place of ReplaySource.

#include "common.hpp"
#include "station.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace ayzek::pipeline {

class ReplaySource {
public:
    // Reads and indexes every 512-byte record. Throws on an unreadable file.
    explicit ReplaySource(const std::string& path);

    [[nodiscard]] const std::string& station() const noexcept { return station_; }
    [[nodiscard]] double first_time() const noexcept { return first_; }
    [[nodiscard]] double last_time() const noexcept { return last_; }

    struct Record {
        double end_time;
        std::size_t comp;
        std::array<std::byte, 512> bytes;
    };
    // Records sorted by end time, i.e. the delivery order of a stream without
    // late records. Late arrivals present in the archive are not reproduced.
    [[nodiscard]] const std::vector<Record>& records() const noexcept { return records_; }

private:
    std::string station_;
    double first_ = 0, last_ = 0;
    std::vector<Record> records_;
};

// Feeds all records of `src` into the station's rings, then sets
// st.ingest_done. When a ring is full it waits for the processor (replay only;
// see ingest.cpp).
void run_ingest(const ReplaySource& src, Station& st, const StreamClock& clock, const std::atomic<bool>& stop);

}  // namespace ayzek::pipeline
