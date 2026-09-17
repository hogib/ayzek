// Unit tests for the Steim2 bit operations, integer record decoding and BTIME
// conversion. Full decoding is validated against ObsPy by
// tools/validate_mseed.py.
#include "mseed.hpp"

#include <array>
#include <cstdlib>
#include <print>

// Aborting check that, unlike assert(), stays active in release builds.
#define CHECK(cond)                                                            \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::println(stderr, "CHECK failed: {}  ({}:{})", #cond, __FILE__,       \
                   __LINE__);                                                  \
      std::abort();                                                            \
    }                                                                          \
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

  // Two 15-bit values at the limits of the range: 16383 and -16384.
  const std::uint32_t w2 = (0b10u << 30) | (0x3FFFu << 15) | 0x4000u;
  std::array<std::int32_t, 2> o2{};
  unpack(w2, 2, 15, o2);
  CHECK((o2 == std::array<std::int32_t, 2>{16383, -16384}));
  std::println("  two 15-bit extremes      ok");

  // Day of year across 29 February in a leap year.
  auto a = btime_to_ns(2024, 60, 0, 0, 0, 0, 0); // 29 Feb 2024
  auto b = btime_to_ns(2024, 61, 0, 0, 0, 0, 0); //  1 Mar 2024
  CHECK(a && b && *b - *a == 86'400'000'000'000LL);
  auto c = btime_to_ns(2025, 1, 0, 0, 0, 1, 5); // +100 us +5 us
  auto d = btime_to_ns(2025, 1, 0, 0, 0, 0, 0);
  CHECK(c && d && *c - *d == 105'000);
  std::println("  btime + microseconds     ok");

  // Uncompressed integer data (encodings 3 and 1), in both byte orders.
  {
    std::array<std::byte, kRecordLength> rec{};
    Header h;
    h.data_offset = 56;
    h.num_samples = 3;
    const std::array<std::uint8_t, 12> be32{0x00, 0x00, 0x00, 0x01, 0xFF, 0xFF,
                                            0xFF, 0xFE, 0x7F, 0xFF, 0xFF, 0xFF};
    std::array<std::int32_t, 3> out{};

    h.encoding = kInt32;
    h.word_order = 1;
    for (std::size_t i = 0; i < be32.size(); ++i)
      rec[56 + i] = std::byte{be32[i]};
    CHECK(decode(rec, h, out) &&
          (out == std::array<std::int32_t, 3>{1, -2, 2147483647}));

    h.word_order = 0; // the same values, little-endian
    for (std::size_t w = 0; w < 3; ++w)
      for (std::size_t b = 0; b < 4; ++b)
        rec[56 + w * 4 + b] = std::byte{be32[w * 4 + 3 - b]};
    out = {};
    CHECK(decode(rec, h, out) &&
          (out == std::array<std::int32_t, 3>{1, -2, 2147483647}));

    h.encoding = kInt16;
    h.word_order = 1;
    h.num_samples = 2;
    const std::array<std::uint8_t, 4> be16{0xFF, 0xFF, 0x7F, 0xFF};
    for (std::size_t i = 0; i < be16.size(); ++i)
      rec[56 + i] = std::byte{be16[i]};
    std::array<std::int32_t, 2> out16{};
    CHECK(decode(rec, h, out16) &&
          (out16 == std::array<std::int32_t, 2>{-1, 32767}));

    h.encoding = kInt32;
    h.num_samples = 115; // (512 - 56) / 4 = 114 fit
    std::array<std::int32_t, 200> big{};
    CHECK(!decode(rec, h, big) &&
          decode(rec, h, big).error() == Error::TooManySamples);
  }
  std::println("  int16 / int32 records    ok");

  std::println("all steim tests passed");
}
