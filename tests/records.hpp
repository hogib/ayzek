#pragma once

// Synthetic 512-byte miniSEED 2 records for the tests: a big-endian fixed
// header, blockette 1000 at byte 48 and the data at byte 64.

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

namespace testrec {

using Bytes = std::array<std::byte, 512>;

struct Spec {
  std::string_view net = "KO", sta = "TST", loc = "", cha = "HHZ";
  // Start time: year, day of year, hour, minute, second, 1e-4 s.
  std::uint16_t year = 2025, doy = 314;
  std::uint8_t hour = 18, minute = 0, second = 0;
  std::uint16_t fract = 0;
  std::int16_t rate_factor = 100, rate_mult = 1;
  std::uint8_t encoding = 3, word_order = 1;
};

inline void put16(Bytes &r, std::size_t off, std::uint16_t v) {
  r[off] = std::byte(v >> 8);
  r[off + 1] = std::byte(v & 0xFF);
}

inline void put32(Bytes &r, std::size_t off, std::uint32_t v) {
  for (std::size_t i = 0; i < 4; ++i)
    r[off + i] = std::byte((v >> (24 - 8 * i)) & 0xFF);
}

inline void put_text(Bytes &r, std::size_t off, std::size_t len,
                     std::string_view s) {
  for (std::size_t i = 0; i < len; ++i)
    r[off + i] = std::byte(i < s.size() ? s[i] : ' ');
}

// Header and blockette 1000 for `n` samples; the data section is left zero.
inline Bytes header(const Spec &s, std::uint16_t n) {
  Bytes r{};
  put_text(r, 0, 6, "000001");
  r[6] = std::byte{'D'};
  r[7] = std::byte{' '};
  put_text(r, 8, 5, s.sta);
  put_text(r, 13, 2, s.loc);
  put_text(r, 15, 3, s.cha);
  put_text(r, 18, 2, s.net);
  put16(r, 20, s.year);
  put16(r, 22, s.doy);
  r[24] = std::byte(s.hour);
  r[25] = std::byte(s.minute);
  r[26] = std::byte(s.second);
  put16(r, 28, s.fract);
  put16(r, 30, n);
  put16(r, 32, static_cast<std::uint16_t>(s.rate_factor));
  put16(r, 34, static_cast<std::uint16_t>(s.rate_mult));
  r[39] = std::byte{1}; // blockettes that follow
  put16(r, 44, 64);     // data offset
  put16(r, 46, 48);     // first blockette
  put16(r, 48, 1000);
  put16(r, 50, 0);
  r[52] = std::byte(s.encoding);
  r[53] = std::byte(s.word_order);
  r[54] = std::byte{9}; // 2^9 = 512 bytes
  return r;
}

// An uncompressed big-endian 32-bit record; at most 112 samples fit.
inline Bytes int32_record(const Spec &s, std::span<const std::int32_t> v) {
  Bytes r = header(s, static_cast<std::uint16_t>(v.size()));
  for (std::size_t i = 0; i < v.size(); ++i)
    put32(r, 64 + 4 * i, static_cast<std::uint32_t>(v[i]));
  return r;
}

} // namespace testrec
