// The two bit-twiddling primitives Steim2 rests on, checked in isolation. The
// real test is bit-exact agreement with ObsPy on archive data; these exist so
// that when that fails, the failure can be localised.
#include "mseed.hpp"

#include <array>
#include <print>
#include <cstdlib>

// CHECK() vanishes under NDEBUG; a test that stops checking in release builds
// is not a test.
#define CHECK(cond)                                                                      \
    do {                                                                                 \
        if (!(cond)) {                                                                   \
            std::println(stderr, "CHECK failed: {}  ({}:{})", #cond, __FILE__, __LINE__); \
            std::abort();                                                                \
        }                                                                                \
    } while (0)

using namespace ayzek::mseed;
using detail::sign_extend;
using detail::unpack;

int main() {
    // Sign extension at every width Steim2 uses.
    CHECK(sign_extend(0x7, 4) == 7);
    CHECK(sign_extend(0x8, 4) == -8);
    CHECK(sign_extend(0xF, 4) == -1);
    CHECK(sign_extend(0x1F, 6) == 31);
    CHECK(sign_extend(0x20, 6) == -32);
    CHECK(sign_extend(0x1FFFFFFF, 30) == 536870911);
    CHECK(sign_extend(0x20000000, 30) == -536870912);
    CHECK(sign_extend(0x3FFFFFFF, 30) == -1);
    std::println("  sign extension           ok");

    // Seven 4-bit values right-aligned under a dnib of 10:
    //   1 -1 7 -8 0 3 -2  ->  0x1 0xF 0x7 0x8 0x0 0x3 0xE
    const std::uint32_t w7 = (0b10u << 30) | 0x1F7803Eu;
    std::array<std::int32_t, 7> o7{};
    unpack(w7, 7, 4, o7);
    CHECK((o7 == std::array<std::int32_t, 7>{1, -1, 7, -8, 0, 3, -2}));
    std::println("  seven 4-bit              ok");

    // Two 15-bit values: 16383 and -16384, the extremes.
    const std::uint32_t w2 = (0b10u << 30) | (0x3FFFu << 15) | 0x4000u;
    std::array<std::int32_t, 2> o2{};
    unpack(w2, 2, 15, o2);
    CHECK((o2 == std::array<std::int32_t, 2>{16383, -16384}));
    std::println("  two 15-bit extremes      ok");

    // Day-of-year arithmetic across a leap year boundary.
    auto a = btime_to_ns(2024, 60, 0, 0, 0, 0, 0);    // 29 Feb 2024
    auto b = btime_to_ns(2024, 61, 0, 0, 0, 0, 0);    //  1 Mar 2024
    CHECK(a && b && *b - *a == 86'400'000'000'000LL);
    auto c = btime_to_ns(2025, 1, 0, 0, 0, 1, 5);     // +100 us +5 us
    auto d = btime_to_ns(2025, 1, 0, 0, 0, 0, 0);
    CHECK(c && d && *c - *d == 105'000);
    std::println("  btime + microseconds     ok");

    std::println("all steim tests passed");
}
