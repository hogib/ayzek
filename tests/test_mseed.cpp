// Decoding against ObsPy: every record of the small files in
// tests/reference/mseed, one per supported encoding and byte order, must give
// exactly the samples ObsPy (libmseed) reads from them, NAME.i32. The files
// are written by tools/make_mseed_fixtures.py, which records how in
// fixtures.json. Skipped (77) if they are missing.
#include "fixture.hpp"
#include "mseed.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <print>
#include <string>
#include <vector>

using namespace ayzek;

int main(int argc, char **argv) {
  const std::filesystem::path dir =
      std::filesystem::path(fixture::root(argc, argv)) / "tests/reference/mseed";
  struct Case {
    const char *name;
    std::uint8_t encoding, word_order;
  };
  for (const Case c : {Case{"steim1", mseed::kSteim1, 1},
                       Case{"steim2", mseed::kSteim2, 1},
                       Case{"int16_be", mseed::kInt16, 1},
                       Case{"int16_le", mseed::kInt16, 0},
                       Case{"int32_be", mseed::kInt32, 1},
                       Case{"int32_le", mseed::kInt32, 0}}) {
    const auto path = dir / (std::string(c.name) + ".mseed");
    const auto ref_path = dir / (std::string(c.name) + ".i32");
    if (!std::filesystem::exists(path) || !std::filesystem::exists(ref_path)) {
      std::println("skipped: {} not found (tools/make_mseed_fixtures.py)",
                   path.string());
      return 77;
    }

    std::vector<std::int32_t> ours;
    std::size_t records = 0;
    std::ifstream f(path, std::ios::binary);
    std::array<std::byte, mseed::kRecordLength> rec{};
    std::array<std::int32_t, mseed::kMaxDiffs> buf{};
    while (f.read(reinterpret_cast<char *>(rec.data()), rec.size())) {
      const auto h = mseed::parse_header(rec);
      CHECK(h.has_value());
      CHECK(h->encoding == c.encoding && h->word_order == c.word_order);
      CHECK(h->sample_rate == 100.0 && h->cha() == "HHZ");
      if (records == 0)
        CHECK(h->start_ns == 1'762'797'600'000'000'000); // 2025-11-10T18:00:00
      const auto n = mseed::decode(rec, *h, buf);
      CHECK(n.has_value());
      ours.insert(ours.end(), buf.begin(),
                  buf.begin() + static_cast<std::ptrdiff_t>(*n));
      ++records;
    }

    // The reference, little-endian whatever the host.
    std::ifstream r(ref_path, std::ios::binary);
    std::vector<std::int32_t> ref;
    for (std::array<unsigned char, 4> b{};
         r.read(reinterpret_cast<char *>(b.data()), 4);)
      ref.push_back(static_cast<std::int32_t>(
          std::uint32_t{b[0]} | std::uint32_t{b[1]} << 8 |
          std::uint32_t{b[2]} << 16 | std::uint32_t{b[3]} << 24));

    CHECK(records > 0 && !ref.empty());
    CHECK(ours == ref);
    std::println("  {:<9} {} records, {:>4} samples, identical to ObsPy",
                 c.name, records, ours.size());
  }
  std::println("all mseed reference tests passed");
}
