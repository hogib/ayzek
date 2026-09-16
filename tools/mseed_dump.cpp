// Decodes a miniSEED file record by record.
//
//   mseed_dump FILE.mseed OUT_DIR
//
// Writes OUT_DIR/<channel>.i32 (all samples in file order, native int32) and
// OUT_DIR/index.csv (one row per record), for comparison with ObsPy by
// tools/validate_mseed.py. Records that fail to parse or decode are reported;
// the exit status is non-zero if there are any.

#include "mseed.hpp"

#include <array>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <print>
#include <string>

using namespace ayzek;

int main(int argc, char** argv) {
    if (argc != 3) {
        std::println(stderr, "usage: {} FILE.mseed OUT_DIR", argv[0]);
        return 2;
    }
    const std::filesystem::path in = argv[1], out = argv[2];
    std::filesystem::create_directories(out);

    std::ifstream f(in, std::ios::binary);
    if (!f) { std::println(stderr, "cannot open {}", in.string()); return 2; }

    std::ofstream index(out / "index.csv");
    index << "record,channel,start_ns,num_samples,sample_rate\n";
    std::map<std::string, std::ofstream> dumps;
    std::map<std::string, std::uint64_t> samples_per_channel;

    std::array<std::byte, mseed::kRecordLength> rec{};
    std::array<std::int32_t, mseed::kMaxDiffs> buf{};
    std::uint64_t n_rec = 0, n_bad = 0, n_samp = 0;

    const auto t0 = std::chrono::steady_clock::now();
    while (f.read(reinterpret_cast<char*>(rec.data()), rec.size())) {
        const std::uint64_t idx = n_rec++;

        auto h = mseed::parse_header(rec);
        if (!h) {
            if (n_bad++ < 10) std::println(stderr, "record {}: {}", idx, mseed::to_string(h.error()));
            continue;
        }
        auto n = mseed::decode(rec, *h, buf);
        if (!n) {
            if (n_bad++ < 10)
                std::println(stderr, "record {} ({}): {}", idx, h->cha(), mseed::to_string(n.error()));
            continue;
        }

        const std::string cha(h->cha());
        auto [it, fresh] = dumps.try_emplace(cha);
        if (fresh) it->second.open(out / (cha + ".i32"), std::ios::binary);
        it->second.write(reinterpret_cast<const char*>(buf.data()),
                         static_cast<std::streamsize>(*n * sizeof(std::int32_t)));

        index << idx << ',' << cha << ',' << h->start_ns << ',' << *n << ',' << h->sample_rate << '\n';
        samples_per_channel[cha] += *n;
        n_samp += *n;
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    const double mb = static_cast<double>(n_rec * mseed::kRecordLength) / 1e6;
    std::println("{} records ({:.0f} MB), {} samples in {:.2f} s  ->  {:.0f} MB/s, {:.1f} M samples/s",
                 n_rec, mb, n_samp, secs, mb / secs, static_cast<double>(n_samp) / secs / 1e6);
    for (const auto& [cha, cnt] : samples_per_channel) std::println("  {}  {} samples", cha, cnt);
    if (n_bad) {
        std::println(stderr, "FAILED: {} of {} records could not be decoded", n_bad, n_rec);
        return 1;
    }
    std::println("all records decoded, integrity checks passed");
    return 0;
}
