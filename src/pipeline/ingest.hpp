#pragma once

// Stage 1 per station: miniSEED records -> decode -> reorder -> component rings.
//
// The source is a file replayed against the stream clock: each record is
// released when stream time reaches its last sample, which is the earliest a
// datalogger could have sent it. A live SeedLink source would replace
// `ReplaySource` and keep everything downstream.

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
    // In release order: by end time, which is delivery order for an on-time
    // stream. The archive's own late records are not reproduced by default.
    [[nodiscard]] const std::vector<Record>& records() const noexcept { return records_; }

private:
    std::string station_;
    double first_ = 0, last_ = 0;
    std::vector<Record> records_;
};

// Runs until the source is exhausted; sets st.ingest_done. Blocks on a full
// ring rather than dropping, which is right for replay and would be wrong live.
void run_ingest(const ReplaySource& src, Station& st, const StreamClock& clock, const std::atomic<bool>& stop);

}  // namespace ayzek::pipeline
