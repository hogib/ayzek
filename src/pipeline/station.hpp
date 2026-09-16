#pragma once

// A station's live state: three component rings on one absolute sample clock.
//
// Ring index i of a component holds absolute position `base + i`. Gaps are
// written as kGap samples rather than skipped, so that mapping stays exact and
// a window over a gap is recognisable by content (DESIGN.md, "gaps").

#include "common.hpp"
#include "ring.hpp"

#include <array>
#include <atomic>
#include <cstdint>
#include <limits>
#include <string>

namespace ayzek::pipeline {

inline constexpr std::int32_t kGap = std::numeric_limits<std::int32_t>::min();
inline constexpr std::size_t kRingCapacity = std::size_t{1} << 17;   // ~21.8 min at 100 Hz

enum Comp : std::size_t { Z = 0, N = 1, E = 2 };

struct ComponentStream {
    SpscRing<std::int32_t, kRingCapacity> ring;
    std::atomic<std::uint64_t> base{kUnset};      // absolute position of ring index 0
    std::atomic<std::uint64_t> gaps{0}, gap_samples{0}, stale{0};

    // Absolute position one past the newest sample, or 0 before any data.
    [[nodiscard]] std::uint64_t end() const noexcept {
        const auto b = base.load(std::memory_order_acquire);
        return b == kUnset ? 0 : b + ring.written();
    }
};

struct Station {
    std::string code;
    double lat = 0, lon = 0;
    std::array<ComponentStream, 3> comp;
    std::atomic<bool> ingest_done{false};
    std::atomic<std::uint64_t> records{0};
};

}  // namespace ayzek::pipeline
