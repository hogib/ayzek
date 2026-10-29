// Robustness of the miniSEED decoder (src/mseed.hpp) on records it was not
// written for: every truncated length, single damaged fields with the error
// each must give, random damage to valid records, and random records. Each
// must give an error, or the record's own sample count within the output,
// never a crash or an access outside the record or the output; and a record
// whose decoding fields are intact must decode to the samples it holds,
// whatever else in it is damaged. The valid records are built here (a Steim2
// encoder after the frame layout in mseed.hpp, and the two integer
// encodings), and first decoded undamaged against the samples they were built
// from. Records and outputs are heap buffers of their exact size, so that
// AddressSanitizer reports any access past either: under
// -Db_sanitize=address,undefined "never outside" is a checked statement.
// Deterministic: splitmix64 from a fixed seed. The callers (ingest.cpp,
// mseed_dump.cpp, replay_check.cpp) read whole 512-byte records and pass each
// to parse_header and then to decode, as run() does here.
//
// A behaviour the decoder's own checks prescribe otherwise is reported as a
// KNOWN ISSUE line, not a failure (see KNOWN_ISSUE); there is one, the start
// time of years after 2262, which -fsanitize=undefined also reports as a
// signed overflow in bits/chrono.h when this test runs under it.

#include "mseed.hpp"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <optional>
#include <print>
#include <string>
#include <vector>

#define CHECK(cond)                                                            \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::println(stderr, "CHECK failed: {}  ({}:{})", #cond, __FILE__,       \
                   __LINE__);                                                  \
      std::abort();                                                            \
    }                                                                          \
  } while (0)

using namespace ayzek::mseed;

namespace {

// As in test_network.cpp: `holds` is what the decoder's checks prescribe.
// While it does not hold, the test prints it and goes on; once it holds, the
// line says so, and the statement can become a CHECK.
int g_known_issues = 0;
void known_issue(bool holds, int line, std::string_view what) {
  if (holds) {
    std::println("NOW HOLDS    test_mseed_fuzz.cpp:{}: {} -- make it a CHECK",
                 line, what);
    return;
  }
  ++g_known_issues;
  std::println("KNOWN ISSUE  test_mseed_fuzz.cpp:{}: {}", line, what);
}
#define KNOWN_ISSUE(holds, what) known_issue((holds), __LINE__, (what))

// splitmix64 (Steele, Lea and Flood 2014): the same sequence everywhere,
// unlike the standard library's distributions.
struct Rng {
  std::uint64_t s;
  std::uint64_t next() {
    std::uint64_t z = (s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }
  unsigned below(unsigned n) { return static_cast<unsigned>(next() % n); }
};

using Record = std::vector<std::byte>;

void put8(Record &r, std::size_t off, unsigned v) {
  r[off] = static_cast<std::byte>(v & 0xFF);
}
void put16(Record &r, std::size_t off, unsigned v) {
  put8(r, off, v >> 8);
  put8(r, off + 1, v);
}
void put32(Record &r, std::size_t off, std::uint32_t v) {
  put16(r, off, v >> 16);
  put16(r, off + 2, v & 0xFFFF);
}
std::uint32_t get32(const Record &r, std::size_t off) {
  std::uint32_t v = 0;
  for (std::size_t i = 0; i < 4; ++i)
    v = v << 8 | static_cast<std::uint32_t>(r[off + i]);
  return v;
}
void text(Record &r, std::size_t off, const char *s) {
  for (std::size_t i = 0; s[i]; ++i)
    r[off + i] = static_cast<std::byte>(s[i]);
}

// 2025-11-10 18:20:51 UTC is 1762798851 s after the epoch (20402 days of
// 86400 s, plus 66051 s), and day 314 of 2025 (304 days to the end of
// October, plus 10); the header's fraction adds 0.1234 s.
constexpr std::int64_t kStartNs = 1762798851'1234'00000;

// A 512-byte record's fixed header and blockette 1000 (SEED 2.4): station
// FUZZ, channel HHZ, network XX, 2025 day 314 18:20:51.1234, 100 Hz, `n`
// samples from byte 64 in encoding `enc` and byte order `order`, record
// length 2^9; blockette 1000 at 48, the only one. Bytes 56-63 are zero.
Record header(unsigned enc, unsigned order, unsigned n) {
  Record r(kRecordLength, std::byte{0});
  text(r, 0, "000001D ");
  text(r, 8, "FUZZ ");
  text(r, 13, "  ");
  text(r, 15, "HHZ");
  text(r, 18, "XX");
  put16(r, 20, 2025);
  put16(r, 22, 314);
  put8(r, 24, 18);
  put8(r, 25, 20);
  put8(r, 26, 51);
  put16(r, 28, 1234);
  put16(r, 30, n);
  put16(r, 32, 100); // sample rate factor
  put16(r, 34, 1);   // and multiplier
  put8(r, 39, 1);    // number of blockettes
  put16(r, 44, 64);  // data offset
  put16(r, 46, 48);  // first blockette
  put16(r, 48, 1000);
  put16(r, 50, 0); // next blockette: none
  put8(r, 52, enc);
  put8(r, 53, order);
  put8(r, 54, 9);
  return r;
}

// Steim2, as tabled in mseed.hpp: a data word holds `count` differences of
// `bits` bits each, under the 2-bit code of the frame's control word and, for
// codes 10 and 11, the 2-bit dnib in its top bits.
struct Packing {
  unsigned count, bits, code, dnib;
};
constexpr Packing kPackings[] = {{7, 4, 0b11, 0b10},  {6, 5, 0b11, 0b01},
                                 {5, 6, 0b11, 0b00},  {4, 8, 0b01, 0},
                                 {3, 10, 0b10, 0b11}, {2, 15, 0b10, 0b10},
                                 {1, 30, 0b10, 0b01}};

bool fits(std::int64_t v, unsigned bits) {
  const std::int64_t lim = std::int64_t{1} << (bits - 1);
  return v >= -lim && v < lim;
}

// Packs the longest prefix of `x` that fits the seven frames from byte 64,
// its first difference taken from `prev`, each word with as many differences
// as it can hold; X0 and Xn in the first frame's words 1 and 2. Sets the
// sample count and returns it; `used[k]` counts the words of kPackings[k].
std::size_t encode_steim2(Record &r, const std::vector<std::int32_t> &x,
                          std::int32_t prev, std::array<unsigned, 7> &used) {
  std::vector<std::int64_t> d(x.size());
  for (std::size_t i = 0; i < x.size(); ++i)
    d[i] = std::int64_t{x[i]} - (i ? x[i - 1] : prev);
  std::size_t n = 0;
  for (std::size_t f = 0; f < 7; ++f) {
    const std::size_t base = 64 + 64 * f;
    std::uint32_t ctrl = 0;
    for (unsigned w = 1; w < 16; ++w) {
      if (f == 0 && w <= 2)
        continue;
      for (std::size_t k = 0; k < std::size(kPackings); ++k) {
        const auto &pk = kPackings[k];
        bool ok = n + pk.count <= d.size();
        for (unsigned i = 0; ok && i < pk.count; ++i)
          ok = fits(d[n + i], pk.bits);
        if (!ok)
          continue;
        std::uint32_t word = pk.code == 0b01 ? 0 : pk.dnib << 30;
        const std::uint32_t mask = (std::uint32_t{1} << pk.bits) - 1;
        for (unsigned i = 0; i < pk.count; ++i)
          word |= (static_cast<std::uint32_t>(d[n + i]) & mask)
                  << ((pk.count - 1 - i) * pk.bits);
        put32(r, base + 4 * w, word);
        ctrl |= pk.code << (30 - 2 * w);
        n += pk.count;
        ++used[k];
        break;
      }
    }
    put32(r, base, ctrl);
  }
  put32(r, 64 + 4, static_cast<std::uint32_t>(x[0]));
  put32(r, 64 + 8, static_cast<std::uint32_t>(x[n - 1]));
  put16(r, 30, static_cast<unsigned>(n));
  return n;
}

// Samples whose differences come in runs of eight of one width, through
// every Steim2 width (4 to 30 bits), alternating in sign so that the signal
// stays within +-2^29.
std::vector<std::int32_t> varied(Rng &rng, std::size_t n) {
  constexpr unsigned kBits[] = {4, 5, 6, 8, 10, 15, 30};
  std::vector<std::int32_t> x(n);
  std::int64_t v = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const unsigned bits = kBits[(i / 8) % 7];
    const auto mag = static_cast<std::int64_t>(rng.below(1u << (bits - 2)));
    v += i % 2 ? mag : -mag;
    x[i] = static_cast<std::int32_t>(v);
  }
  return x;
}

// Integer data from byte 64, most significant byte first if `order` is 1.
Record int_record(unsigned enc, unsigned order,
                  const std::vector<std::int32_t> &x) {
  Record r = header(enc, order, static_cast<unsigned>(x.size()));
  const std::size_t w = enc == kInt16 ? 2 : 4;
  for (std::size_t i = 0; i < x.size(); ++i)
    for (std::size_t b = 0; b < w; ++b)
      put8(r, 64 + i * w + b,
           static_cast<std::uint32_t>(x[i]) >>
               (8 * (order == 1 ? w - 1 - b : b)));
  return r;
}

struct Result {
  std::optional<Error> error;
  Header h;
  std::vector<std::int32_t> samples;
};

// parse_header, then decode into a heap output of `out` samples (kMaxDiffs,
// what the callers use, by default).
Result run(const Record &r, std::size_t out = kMaxDiffs) {
  Result res;
  const auto h = parse_header(r);
  if (!h) {
    res.error = h.error();
    return res;
  }
  res.h = *h;
  std::vector<std::int32_t> buf(out);
  const auto n = decode(r, *h, buf);
  if (!n) {
    res.error = n.error();
    return res;
  }
  CHECK(*n == h->num_samples && *n <= buf.size());
  res.samples.assign(buf.begin(), buf.begin() + static_cast<long>(*n));
  return res;
}

// Both output sizes: what the callers pass, and exactly the record's count
// (ASan then sees a write past the count). Returns the first.
Result run_both(const Record &r) {
  const Result a = run(r);
  if (const auto h = parse_header(r); h && h->num_samples <= kMaxDiffs) {
    const Result b = run(r, h->num_samples);
    CHECK(b.error == a.error && b.samples == a.samples);
  }
  return a;
}

// The bytes decoding does not read: sequence number, quality and the
// identifiers (0-19); the activity, I/O and quality flags, the number of
// blockettes and the time correction (36-43), which change the start time
// only; and 56-63, outside the blockette chain and before the data.
bool inert(std::size_t i) {
  return i < 20 || (i >= 36 && i < 44) || (i >= 56 && i < 64);
}

std::map<std::string, unsigned> g_outcomes;
void tally(const char *phase, const Result &r) {
  ++g_outcomes[std::format("{:<16} {}", phase,
                           r.error ? to_string(*r.error) : "decoded")];
}

} // namespace

int main() {
  Rng rng{20251110};

  // --- valid records, undamaged ---------------------------------------------
  const auto x = varied(rng, 900);
  Record steim = header(kSteim2, 1, 0);
  std::array<unsigned, 7> used{};
  const std::size_t n = encode_steim2(steim, x, 17, used);
  const std::vector<std::int32_t> xs(x.begin(),
                                     x.begin() + static_cast<long>(n));
  CHECK(n > 100 && xs[n - 2] != xs[n - 1]);
  for (const unsigned u : used)
    CHECK(u > 0); // every Steim2 packing is exercised
  const auto s = run_both(steim);
  CHECK(!s.error && s.samples == xs);
  CHECK(s.h.sta() == "FUZZ" && s.h.cha() == "HHZ" && s.h.net() == "XX" &&
        s.h.loc().empty() && s.h.sample_rate == 100.0 &&
        s.h.start_ns == kStartNs);
  std::println("  Steim2 record: {} samples in 103 words, every packing used",
               n);

  std::vector<std::int32_t> x32(112), x16(224);
  for (auto &v : x32)
    v = static_cast<std::int32_t>(static_cast<std::uint32_t>(rng.next()));
  for (auto &v : x16)
    v = static_cast<std::int16_t>(static_cast<std::uint16_t>(rng.next()));
  struct Base {
    const char *name;
    Record r;
    std::vector<std::int32_t> samples;
  };
  std::vector<Base> bases = {{"Steim2", steim, xs},
                             {"int32 BE", int_record(kInt32, 1, x32), x32},
                             {"int32 LE", int_record(kInt32, 0, x32), x32},
                             {"int16 BE", int_record(kInt16, 1, x16), x16},
                             {"int16 LE", int_record(kInt16, 0, x16), x16}};
  for (const auto &b : bases) {
    const auto r = run_both(b.r);
    CHECK(!r.error && r.samples == b.samples);
  }
  std::println("  int32 (112 samples, the most from byte 64) and int16 (224), "
               "both byte orders");

  // --- truncated -------------------------------------------------------------
  for (const auto &b : bases) {
    for (std::size_t len = 0; len < kRecordLength; ++len) {
      const Record t(b.r.begin(), b.r.begin() + static_cast<long>(len));
      CHECK(run(t).error == Error::BadRecordLength);
    }
    Record longer = b.r;
    longer.push_back(std::byte{0});
    CHECK(run(longer).error == Error::BadRecordLength);
  }
  std::println(
      "  every truncation (0-511 bytes) and 513 bytes: bad record length");

  // --- one damaged field -----------------------------------------------------
  // Each on a copy of the Steim2 record (n samples, X0 at byte 68, Xn at 72,
  // the first data word at 76, word 3 of frame 0, whose code is bits 24-25 of
  // the control word at 64), with the outcome the format and mseed.hpp's
  // checks prescribe: an error, or a decoded record (nullopt).
  const unsigned nn = static_cast<unsigned>(n);
  const auto x0 = static_cast<std::uint32_t>(xs[0]),
             xn = static_cast<std::uint32_t>(xs[n - 1]);
  struct Damage {
    const char *what;
    std::optional<Error> expect;
    void (*apply)(Record &, unsigned, std::uint32_t, std::uint32_t);
  };
  const Damage damages[] = {
      {"year 1899", Error::UnsupportedByteOrder,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 20, 1899);
       }},
      {"year 2501", Error::UnsupportedByteOrder,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 20, 2501);
       }},
      {"day 0", Error::UnsupportedByteOrder,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 22, 0);
       }},
      {"day 367", Error::UnsupportedByteOrder,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 22, 367);
       }},
      {"year little-endian (0xE907 = 59655)", Error::UnsupportedByteOrder,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 20, 0xE907);
       }},
      {"hour 24", Error::BadTime,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put8(r, 24, 24);
       }},
      {"minute 60", Error::BadTime,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put8(r, 25, 60);
       }},
      {"second 61", Error::BadTime,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put8(r, 26, 61);
       }},
      {"fraction 10000", Error::BadTime,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 28, 10000);
       }},
      {"no blockette", Error::MissingBlockette1000,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 46, 0);
       }},
      {"first blockette at 505 (505 + 8 > 512)", Error::MissingBlockette1000,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 46, 505);
       }},
      {"blockette 1001, not 1000", Error::MissingBlockette1000,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 48, 1001);
       }},
      {"chain into the zeros at 56", Error::MissingBlockette1000,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 46, 56);
       }},
      {"a cycle without 1000", Error::MissingBlockette1000,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 48, 2000);
         put16(r, 50, 48);
       }},
      {"Steim1 (10)", Error::UnsupportedEncoding,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put8(r, 52, 10);
       }},
      {"encoding 0", Error::UnsupportedEncoding,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put8(r, 52, 0);
       }},
      {"Steim2 little-endian", Error::UnsupportedEncoding,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put8(r, 53, 0);
       }},
      {"byte order 2", Error::UnsupportedEncoding,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put8(r, 53, 2);
       }},
      {"256-byte record length", Error::UnsupportedRecordLength,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put8(r, 54, 8);
       }},
      {"4096-byte record length", Error::UnsupportedRecordLength,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put8(r, 54, 12);
       }},
      {"record length 2^255", Error::UnsupportedRecordLength,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put8(r, 54, 255);
       }},
      {"data offset 0", Error::BadDataOffset,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 44, 0);
       }},
      {"data offset 47", Error::BadDataOffset,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 44, 47);
       }},
      {"data offset 100 (412 bytes, not whole frames)", Error::BadDataOffset,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 44, 100);
       }},
      {"data offset 512", Error::BadDataOffset,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 44, 512);
       }},
      {"data offset 65535", Error::BadDataOffset,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 44, 65535);
       }},
      {"one sample more than encoded", Error::TooManySamples,
       [](Record &r, unsigned m, std::uint32_t, std::uint32_t) {
         put16(r, 30, m + 1);
       }},
      {"736 samples (7 x 15 x 7 + 1)", Error::TooManySamples,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 30, 736);
       }},
      {"65535 samples", Error::TooManySamples,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 30, 65535);
       }},
      {"one sample fewer (the last but one is not Xn)",
       Error::IntegrityMismatch,
       [](Record &r, unsigned m, std::uint32_t, std::uint32_t) {
         put16(r, 30, m - 1);
       }},
      {"Xn + 1", Error::IntegrityMismatch,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t v) {
         put32(r, 72, v + 1);
       }},
      {"X0 + 1 (every sample + 1)", Error::IntegrityMismatch,
       [](Record &r, unsigned, std::uint32_t v, std::uint32_t) {
         put32(r, 68, v + 1);
       }},
      {"code 10, dnib 00", Error::InvalidSteimCode,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put32(r, 64, (get32(r, 64) & ~(3u << 24)) | 2u << 24);
         put32(r, 76, get32(r, 76) & 0x3FFFFFFF);
       }},
      {"code 11, dnib 11", Error::InvalidSteimCode,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put32(r, 64, get32(r, 64) | 3u << 24);
         put32(r, 76, get32(r, 76) | 3u << 30);
       }},
      {"zero samples", std::nullopt,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 30, 0);
       }},
      {"leap second (60)", std::nullopt,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put8(r, 26, 60);
       }},
      {"blockette 1000 pointing to itself", std::nullopt,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 50, 48);
       }},
      {"time correction +1 s, not applied", std::nullopt,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put32(r, 40, 10000);
       }},
      {"time correction +1 s, applied (activity bit 0x02)", std::nullopt,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put32(r, 40, 10000);
         put8(r, 36, 0x02);
       }},
      {"blockette 1001 at 56: +50 us", std::nullopt,
       [](Record &r, unsigned, std::uint32_t, std::uint32_t) {
         put16(r, 50, 56);
         put16(r, 56, 1001);
         put8(r, 61, 50);
       }},
  };
  // The start times of the six that decode, in that order: 18:20:51.1234;
  // 9 s later (18:21:00.1234); unchanged; + 1 s; unchanged; + 50 us.
  const std::int64_t starts[] = {kStartNs, kStartNs + 9'000'000'000,
                                 kStartNs, kStartNs + 1'000'000'000,
                                 kStartNs, kStartNs + 50'000};
  std::size_t decoded = 0;
  for (const auto &d : damages) {
    Record r = steim;
    d.apply(r, nn, x0, xn);
    const auto res = run_both(r);
    if (res.error != d.expect)
      std::println(stderr, "{}: {}", d.what,
                   res.error ? to_string(*res.error) : "decoded");
    CHECK(res.error == d.expect);
    if (!d.expect) {
      CHECK(res.h.start_ns == starts[decoded++]);
      CHECK(res.samples.size() == res.h.num_samples);
      if (res.h.num_samples)
        CHECK(res.samples == xs);
    }
  }
  CHECK(decoded == std::size(starts));
  // Integer data: the count must fit between the offset and the end, and the
  // output must hold it.
  {
    Record r = bases[1].r;
    put16(r, 30, 113); // (512 - 64) / 4 = 112
    CHECK(run(r).error == Error::TooManySamples);
    r = bases[3].r;
    put16(r, 30, 225); // (512 - 64) / 2 = 224
    CHECK(run(r).error == Error::TooManySamples);
    r = bases[1].r;
    put16(r, 44, 48); // integers need no whole frames: 116 fit from 48
    put16(r, 30, 116);
    CHECK(!run(r).error && run(r).samples.size() == 116);
    put16(r, 30, 117);
    CHECK(run(r).error == Error::TooManySamples);
    for (const auto &b : bases)
      CHECK(run(b.r, b.samples.size() - 1).error == Error::TooManySamples);
  }
  std::println("  {} single damages: the error each must give",
               std::size(damages));

  // Years 2263 to 2500 pass parse_header's range check (mseed.hpp:135, year
  // <= 2500), but nanoseconds since the epoch exceed int64 from 2^63 ns =
  // 9223372036.854775807 s, 2262-04-11 23:47:16.85 UTC, on: the duration_cast
  // to nanoseconds in btime_to_ns (mseed.hpp:143) then overflows a signed
  // integer, undefined behaviour, which in practice wraps to a start before
  // 1970. 1 January 2263 00:00:00 is 293 years after 1970, 71 of them leap
  // (1972 to 2260), so 107016 days = 9246182400 s, past the limit. Its record
  // should be refused as a bad start time; it parses.
  {
    Record r = steim;
    put16(r, 20, 2263);
    put16(r, 22, 1);
    put8(r, 24, 0);
    put8(r, 25, 0);
    put8(r, 26, 0);
    put16(r, 28, 0);
    const auto res = run_both(r);
    KNOWN_ISSUE(res.error == Error::BadTime,
                std::format("a record dated 2263-01-01 parses with start_ns {} "
                            "(nanoseconds since 1970 overflow int64 after "
                            "2262-04-11; mseed.hpp:135 accepts years to 2500, "
                            "mseed.hpp:143 overflows for them)",
                            res.error ? -1 : res.h.start_ns));
  }

  // --- random damage to valid records ----------------------------------------
  // 1 to 8 random bytes of a record replaced, or one bit flipped. Where every
  // damaged byte is inert, the record must decode as before.
  std::size_t inert_cases = 0;
  for (const auto &b : bases)
    for (int i = 0; i < 20000; ++i) {
      Record r = b.r;
      bool all_inert = true;
      const unsigned k = 1 + rng.below(8);
      for (unsigned j = 0; j < k; ++j) {
        const std::size_t at = rng.below(kRecordLength);
        all_inert = all_inert && inert(at);
        if (rng.below(2))
          r[at] ^= static_cast<std::byte>(1u << rng.below(8));
        else
          r[at] = static_cast<std::byte>(rng.next() & 0xFF);
      }
      const auto res = run_both(r);
      tally(b.name, res);
      if (all_inert) {
        ++inert_cases;
        CHECK(!res.error && res.samples == b.samples);
      }
    }
  std::println("  100000 randomly damaged records, {} damaged in inert bytes "
               "only, which decode unchanged",
               inert_cases);

  // --- random records with a plausible header --------------------------------
  // A valid date and blockette 1000 at 48; the encoding, byte order, record
  // length, sample count and data offset mostly plausible, sometimes anything;
  // the rest, data included, random.
  for (int i = 0; i < 100000; ++i) {
    Record r(kRecordLength);
    for (auto &c : r)
      c = static_cast<std::byte>(rng.next() & 0xFF);
    put16(r, 20, 1900 + rng.below(601));
    put16(r, 22, 1 + rng.below(366));
    put8(r, 24, rng.below(24));
    put8(r, 25, rng.below(60));
    put8(r, 26, rng.below(61));
    put16(r, 28, rng.below(10000));
    put16(r, 46, 48);
    put16(r, 48, 1000);
    if (rng.below(4))
      put16(r, 50, 0);
    const unsigned encs[] = {kInt16, kInt32, kSteim2, kSteim2};
    if (rng.below(8))
      put8(r, 52, encs[rng.below(4)]);
    if (rng.below(8))
      put8(r, 53, 1);
    if (rng.below(8))
      put8(r, 54, 9);
    if (rng.below(4))
      put16(r, 44,
            rng.below(2) ? 64 * (1 + rng.below(7)) : 48 + rng.below(464));
    if (rng.below(4))
      put16(r, 30, rng.below(800));
    tally("plausible", run_both(r));
  }
  // --- random records --------------------------------------------------------
  for (int i = 0; i < 100000; ++i) {
    Record r(kRecordLength);
    for (auto &c : r)
      c = static_cast<std::byte>(rng.next() & 0xFF);
    tally("random", run_both(r));
  }
  std::println(
      "  100000 random records with a plausible header, 100000 random records");
  for (const auto &[what, count] : g_outcomes)
    std::println("    {:<56} {:>6}", what, count);
  std::println("mseed decoder: errors, never a crash, on damaged and random "
               "records; {} known issue{} reported above",
               g_known_issues, g_known_issues == 1 ? "" : "s");
}
