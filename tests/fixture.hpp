#pragma once

// Shared by the tests that compare against tools/export_models.py outputs.
// Fixtures are generated, not committed; a test without them exits 77, which
// Meson reports as skipped rather than passed.

#include "weights.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <print>
#include <span>
#include <string>

#define CHECK(cond)                                                                      \
    do {                                                                                 \
        if (!(cond)) {                                                                   \
            std::println(stderr, "CHECK failed: {}  ({}:{})", #cond, __FILE__, __LINE__); \
            std::abort();                                                                \
        }                                                                                \
    } while (0)

namespace fixture {

inline std::string root(int argc, char** argv) { return argc > 1 ? argv[1] : "."; }

inline ayzek::Weights load_or_skip(const std::string& path) {
    if (!std::filesystem::exists(path)) {
        std::println("skipped: {} not found (run tools/export_models.py)", path);
        std::exit(77);
    }
    return ayzek::Weights::load(path);
}

// Largest |a - b| and largest |a - b| / max(|b|, floor).
struct Diff {
    double abs = 0, rel = 0;
};

template <typename A, typename B>
Diff compare(std::span<const A> a, std::span<const B> b, double floor = 1e-6) {
    CHECK(a.size() == b.size());
    Diff d;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double e = std::abs(static_cast<double>(a[i]) - static_cast<double>(b[i]));
        d.abs = std::max(d.abs, e);
        d.rel = std::max(d.rel, e / std::max(std::abs(static_cast<double>(b[i])), floor));
    }
    return d;
}

}  // namespace fixture
