#pragma once

// miniSEED 2 record parsing and Steim2 decompression.
//
// Scope is deliberately narrow: exactly what the TDVMS archive contains, and an
// explicit error for anything else. A scan of all 1,484,074 records in a
// 21-day DEMI chunk found every one identical in structure -- quality D, 512
// bytes, blockette 1000 at offset 48 and 1001 at 56, data at 64, Steim2,
// big-endian, no time correction. A lenient decoder that guessed at the
// unusual cases would convert a malformed record into plausible wrong samples,
// which is the one failure a detector cannot recover from.
//
// Steim2 is lossless integer compression, so the acceptance test is exact:
// every decoded sample must equal what ObsPy (libmseed) produces, bit for bit.

#include <array>
#include <bit>
#include <chrono>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <expected>
#include <span>
#include <string_view>

namespace ayzek::mseed {

inline constexpr std::size_t kRecordLength = 512;
inline constexpr std::size_t kFrameLength = 64;

enum class Error : std::uint8_t {
    BadRecordLength,
    UnsupportedByteOrder,
    BadTime,
    MissingBlockette1000,
    UnsupportedEncoding,
    UnsupportedRecordLength,
    BadDataOffset,
    TooManySamples,
    InvalidSteimCode,
    IntegrityMismatch,   // last decoded sample != the record's stored Xn
};

[[nodiscard]] constexpr std::string_view to_string(Error e) noexcept {
    switch (e) {
        case Error::BadRecordLength:         return "bad record length";
        case Error::UnsupportedByteOrder:    return "unsupported byte order";
        case Error::BadTime:                 return "bad start time";
        case Error::MissingBlockette1000:    return "missing blockette 1000";
        case Error::UnsupportedEncoding:     return "unsupported encoding";
        case Error::UnsupportedRecordLength: return "unsupported record length";
        case Error::BadDataOffset:           return "bad data offset";
        case Error::TooManySamples:          return "more samples than the record can hold";
        case Error::InvalidSteimCode:        return "invalid Steim2 code";
        case Error::IntegrityMismatch:       return "Steim2 integrity check failed";
    }
    return "unknown";
}

// ---------------------------------------------------------------------------
// Trap 1: reading integers out of a byte buffer.
//
// In C you might cast `(uint32_t*)buf` or pun through a union. In C++ both are
// undefined behaviour (strict aliasing; C++ does not bless union punning the
// way C does). `std::memcpy` into a local is the sanctioned route, and
// compilers turn it into a single load -- it costs nothing. Rust's
// `u32::from_be_bytes` is the same idea with nicer spelling.
//
// `std::byteswap` is new in C++23. The comparison against `std::endian::native`
// is resolved at compile time, so on a little-endian machine this is a load and
// a bswap, and on a big-endian one just the load.
template <std::integral T>
[[nodiscard]] inline T read_be(std::span<const std::byte> s, std::size_t off) noexcept {
    T v;
    std::memcpy(&v, s.data() + off, sizeof(T));
    return std::endian::native == std::endian::big ? v : std::byteswap(v);
}

struct Header {
    std::array<char, 3> network{};
    std::array<char, 6> station{};
    std::array<char, 3> location{};
    std::array<char, 4> channel{};

    std::int64_t start_ns{};     // nanoseconds since the Unix epoch, UTC
    double sample_rate{};
    std::uint16_t num_samples{};
    std::uint16_t data_offset{};

    [[nodiscard]] std::string_view net() const noexcept { return trimmed(network); }
    [[nodiscard]] std::string_view sta() const noexcept { return trimmed(station); }
    [[nodiscard]] std::string_view loc() const noexcept { return trimmed(location); }
    [[nodiscard]] std::string_view cha() const noexcept { return trimmed(channel); }

    // SEED pads identifiers with spaces; the arrays carry one extra byte so a
    // full-width field still has room, and trimming stops at either.
    template <std::size_t N>
    [[nodiscard]] static std::string_view trimmed(const std::array<char, N>& a) noexcept {
        std::size_t n = 0;
        while (n < N && a[n] != ' ' && a[n] != '\0') ++n;
        return {a.data(), n};
    }
};

// SEED BTIME plus blockette 1001's microseconds -> epoch nanoseconds.
//
// std::chrono does the calendar arithmetic, including day-of-year, so there is
// no hand-written days-from-civil to get wrong. A leap second (sec == 60) rolls
// into the next minute, which matches ObsPy.
[[nodiscard]] inline std::expected<std::int64_t, Error>
btime_to_ns(std::uint16_t year, std::uint16_t doy, std::uint8_t hour, std::uint8_t min,
            std::uint8_t sec, std::uint16_t fract_1e4, std::int8_t micro) noexcept {
    using namespace std::chrono;
    if (year < 1900 || year > 2500 || doy < 1 || doy > 366 || hour > 23 || min > 59 ||
        sec > 60 || fract_1e4 > 9999) {
        return std::unexpected(Error::BadTime);
    }
    const sys_days day0 = sys_days{std::chrono::year{year} / January / 1} + days{doy - 1};
    const auto t = day0 + hours{hour} + minutes{min} + seconds{sec} +
                   microseconds{fract_1e4 * 100LL + micro};
    return duration_cast<nanoseconds>(t.time_since_epoch()).count();
}

[[nodiscard]] inline double sample_rate(std::int16_t factor, std::int16_t mult) noexcept {
    if (factor == 0 || mult == 0) return 0.0;
    const double f = factor, m = mult;
    if (factor > 0 && mult > 0) return f * m;
    if (factor > 0 && mult < 0) return -f / m;
    if (factor < 0 && mult > 0) return -m / f;
    return 1.0 / (f * m);
}

[[nodiscard]] inline std::expected<Header, Error>
parse_header(std::span<const std::byte> rec) noexcept {
    if (rec.size() != kRecordLength) return std::unexpected(Error::BadRecordLength);

    // Byte order is not flagged in the fixed header itself; the convention is
    // to read the year and see whether it is sane. This archive is big-endian
    // throughout, and anything else is refused rather than guessed at.
    const auto year = read_be<std::uint16_t>(rec, 20);
    const auto doy = read_be<std::uint16_t>(rec, 22);
    if (year < 1900 || year > 2500 || doy < 1 || doy > 366) {
        return std::unexpected(Error::UnsupportedByteOrder);
    }

    Header h;
    std::memcpy(h.station.data(), rec.data() + 8, 5);
    std::memcpy(h.location.data(), rec.data() + 13, 2);
    std::memcpy(h.channel.data(), rec.data() + 15, 3);
    std::memcpy(h.network.data(), rec.data() + 18, 2);

    h.num_samples = read_be<std::uint16_t>(rec, 30);
    h.sample_rate = sample_rate(read_be<std::int16_t>(rec, 32), read_be<std::int16_t>(rec, 34));
    h.data_offset = read_be<std::uint16_t>(rec, 44);

    const auto activity = static_cast<std::uint8_t>(rec[36]);
    const auto correction = read_be<std::int32_t>(rec, 40);   // units of 1e-4 s
    const auto first_blk = read_be<std::uint16_t>(rec, 46);

    bool have_1000 = false;
    std::int8_t micro = 0;
    // Walk the blockette chain. Each carries the offset of the next, 0 ends it.
    // The bound on iterations stops a corrupt self-referencing chain from
    // looping forever.
    std::uint16_t off = first_blk;
    for (int guard = 0; off != 0 && guard < 16; ++guard) {
        if (off + 8u > kRecordLength) return std::unexpected(Error::MissingBlockette1000);
        const auto type = read_be<std::uint16_t>(rec, off);
        const auto next = read_be<std::uint16_t>(rec, off + 2u);
        if (type == 1000) {
            have_1000 = true;
            const auto encoding = static_cast<std::uint8_t>(rec[off + 4u]);
            const auto word_order = static_cast<std::uint8_t>(rec[off + 5u]);
            const auto len_exp = static_cast<std::uint8_t>(rec[off + 6u]);
            if (encoding != 11 || word_order != 1) {
                return std::unexpected(Error::UnsupportedEncoding);
            }
            if (len_exp >= 16 || (std::size_t{1} << len_exp) != kRecordLength) {
                return std::unexpected(Error::UnsupportedRecordLength);
            }
        } else if (type == 1001) {
            micro = static_cast<std::int8_t>(rec[off + 5u]);
        }
        off = next;
    }
    if (!have_1000) return std::unexpected(Error::MissingBlockette1000);

    if (h.data_offset < 48 || h.data_offset >= kRecordLength ||
        (kRecordLength - h.data_offset) % kFrameLength != 0) {
        return std::unexpected(Error::BadDataOffset);
    }

    auto t = btime_to_ns(year, doy, static_cast<std::uint8_t>(rec[24]),
                         static_cast<std::uint8_t>(rec[25]), static_cast<std::uint8_t>(rec[26]),
                         read_be<std::uint16_t>(rec, 28), micro);
    if (!t) return std::unexpected(t.error());
    h.start_ns = *t;
    // Activity flag 0x02 means the time correction is already folded into the
    // start time. Otherwise it still has to be applied.
    if ((activity & 0x02) == 0) h.start_ns += std::int64_t{correction} * 100'000;

    return h;
}

// ---------------------------------------------------------------------------
// Steim2.
//
// Data is a run of 64-byte frames. Each frame is one control word holding 16
// two-bit codes -- code i describes word i -- followed by 15 data words. The
// first frame's words 1 and 2 are special: X0, the first sample, and Xn, the
// last, stored so the decoder can check itself.
//
// Everything else is first differences, packed at a width chosen per word:
//
//   code 01              four  8-bit
//   code 10  dnib 01     one  30-bit
//            dnib 10     two  15-bit
//            dnib 11     three 10-bit
//   code 11  dnib 00     five  6-bit
//            dnib 01     six   5-bit
//            dnib 10     seven 4-bit
//
// with dnib the word's top two bits and the values right-aligned below them.

namespace detail {

// Trap 2: sign extension from an arbitrary bit width.
//
// `(v ^ m) - m` with m the sign bit is branchless and exact. The subtraction
// is done in *unsigned* arithmetic and only then converted, because signed
// overflow in C++ is undefined behaviour -- not wrapping, as in Rust release
// builds, and not a panic, as in Rust debug builds. Undefined. The conversion
// from unsigned to signed, by contrast, has been defined as two's complement
// since C++20.
[[nodiscard]] constexpr std::int32_t sign_extend(std::uint32_t v, unsigned bits) noexcept {
    const std::uint32_t m = std::uint32_t{1} << (bits - 1);
    return static_cast<std::int32_t>((v ^ m) - m);
}

// Unpack `n` values of `bits` width, most significant first, from the low
// `n * bits` bits of `word`.
constexpr void unpack(std::uint32_t word, unsigned n, unsigned bits,
                      std::span<std::int32_t> out) noexcept {
    const std::uint32_t mask = (std::uint32_t{1} << bits) - 1;
    for (unsigned i = 0; i < n; ++i) {
        out[i] = sign_extend((word >> ((n - 1 - i) * bits)) & mask, bits);
    }
}

}  // namespace detail

// Worst case: 7 frames, 15 words, 7 differences per word.
inline constexpr std::size_t kMaxDiffs = 7 * 15 * 7;

// Decode `h.num_samples` samples from `rec` into `out`, which must hold that
// many. No allocation: differences go to a fixed scratch array on the stack.
[[nodiscard]] inline std::expected<std::size_t, Error>
decode_steim2(std::span<const std::byte> rec, const Header& h,
              std::span<std::int32_t> out) noexcept {
    const std::size_t n = h.num_samples;
    if (n > kMaxDiffs || n > out.size()) return std::unexpected(Error::TooManySamples);
    if (n == 0) return std::size_t{0};

    std::array<std::int32_t, kMaxDiffs> diffs{};
    std::size_t nd = 0;
    std::int32_t x0 = 0, xn = 0;

    const std::size_t frames = (kRecordLength - h.data_offset) / kFrameLength;
    for (std::size_t f = 0; f < frames && nd < n; ++f) {
        const std::size_t base = h.data_offset + f * kFrameLength;
        const auto ctrl = read_be<std::uint32_t>(rec, base);

        for (unsigned w = 1; w < 16 && nd < n; ++w) {
            const auto word = read_be<std::uint32_t>(rec, base + 4u * w);
            const unsigned code = (ctrl >> (30 - 2 * w)) & 0b11;
            const unsigned dnib = word >> 30;
            std::span<std::int32_t> dst(diffs.data() + nd, kMaxDiffs - nd);

            switch (code) {
                case 0b00:
                    if (f == 0 && w == 1) x0 = static_cast<std::int32_t>(word);
                    if (f == 0 && w == 2) xn = static_cast<std::int32_t>(word);
                    break;
                case 0b01:
                    detail::unpack(word, 4, 8, dst);
                    nd += 4;
                    break;
                case 0b10:
                    switch (dnib) {
                        case 0b01: detail::unpack(word, 1, 30, dst); nd += 1; break;
                        case 0b10: detail::unpack(word, 2, 15, dst); nd += 2; break;
                        case 0b11: detail::unpack(word, 3, 10, dst); nd += 3; break;
                        default: return std::unexpected(Error::InvalidSteimCode);
                    }
                    break;
                case 0b11:
                    switch (dnib) {
                        case 0b00: detail::unpack(word, 5, 6, dst); nd += 5; break;
                        case 0b01: detail::unpack(word, 6, 5, dst); nd += 6; break;
                        case 0b10: detail::unpack(word, 7, 4, dst); nd += 7; break;
                        default: return std::unexpected(Error::InvalidSteimCode);
                    }
                    break;
            }
        }
    }
    if (nd < n) return std::unexpected(Error::TooManySamples);

    // Integrate. The first difference bridges from the *previous* record's last
    // sample, so it is discarded here and X0 stands in for sample 0.
    //
    // Trap 3: accumulate in uint32. For valid data the running sum never
    // overflows, but a decoder must not have undefined behaviour on malformed
    // input, and signed overflow is exactly that.
    std::uint32_t acc = static_cast<std::uint32_t>(x0);
    out[0] = x0;
    for (std::size_t i = 1; i < n; ++i) {
        acc += static_cast<std::uint32_t>(diffs[i]);
        out[i] = static_cast<std::int32_t>(acc);
    }

    // Trap 4, and the reason Xn exists: without this check a mis-decoded frame
    // produces a smooth, plausible, entirely wrong waveform.
    if (out[n - 1] != xn) return std::unexpected(Error::IntegrityMismatch);
    return n;
}

}  // namespace ayzek::mseed
