#include "ingest.hpp"

#include "mseed.hpp"
#include "reorder.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <stdexcept>
#include <thread>

namespace ayzek::pipeline {

namespace {

std::size_t component_of(std::string_view cha) {
    switch (cha.empty() ? '?' : cha.back()) {
        case 'Z': return Z;
        case 'N': case '1': return N;
        case 'E': case '2': return E;
        default: return 3;
    }
}

}  // namespace

ReplaySource::ReplaySource(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::array<std::byte, 512> raw{};
    std::map<std::string, std::size_t> skipped;     // reason -> records
    while (f.read(reinterpret_cast<char*>(raw.data()), raw.size())) {
        auto h = mseed::parse_header(raw);
        if (!h) { ++skipped[std::string(mseed::to_string(h.error()))]; continue; }
        const std::size_t comp = component_of(h->cha());
        if (comp > E || h->cha().substr(0, 2) != "HH" || h->sample_rate != 100.0) { ++skipped["not a 100 Hz HH? channel"]; continue; }
        if (station_.empty()) station_ = std::string(h->sta());
        const double start = static_cast<double>(h->start_ns) / 1e9;
        const double end = start + h->num_samples / h->sample_rate;
        records_.push_back({end, comp, raw});
    }
    if (records_.empty()) throw std::runtime_error(path + ": no usable 100 Hz HH? records");
    std::stable_sort(records_.begin(), records_.end(), [](const Record& a, const Record& b) { return a.end_time < b.end_time; });
    first_ = records_.front().end_time;
    last_ = records_.back().end_time;
    for (const auto& [reason, count] : skipped)
        Log::get().line("replay", "33", "{}: skipped {} records ({})", station_, count, reason);
}

void run_ingest(const ReplaySource& src, Station& st, const StreamClock& clock, const std::atomic<bool>& stop) {
    // One reorderer per component with max_lateness 0: holes become gaps at once
    // and late records are discarded as stale (docs/DESIGN.md).
    std::array<std::unique_ptr<Reorderer<8, mseed::kMaxDiffs>>, 3> reorder;
    for (auto& r : reorder) r = std::make_unique<Reorderer<8, mseed::kMaxDiffs>>(0);
    std::array<std::int32_t, mseed::kMaxDiffs> buf{};
    std::array<std::int32_t, 4096> gapfill{};
    gapfill.fill(kGap);

    for (const auto& rec : src.records()) {
        if (stop.load(std::memory_order_relaxed)) break;
        while (rec.end_time > clock.now() && !stop.load(std::memory_order_relaxed))
            std::this_thread::sleep_for(std::chrono::milliseconds(2));

        auto h = mseed::parse_header(rec.bytes);
        auto n = mseed::decode(rec.bytes, *h, buf);
        if (!n) continue;
        const auto pos = static_cast<std::uint64_t>(h->start_ns) / 10'000'000;
        auto& cs = st.comp[rec.comp];

        auto push = [&](std::span<const std::int32_t> s) {
            int waited_ms = 0;
            while (!s.empty()) {
                // Ring full: wait for the processor to raise its floor. The wait
                // is limited to 60 s of wall time, because the processor cannot
                // advance while another component of the station has no data.
                // After the limit the samples are dropped, counted and reported.
                const auto room = kRingCapacity - (cs.ring.written() - cs.ring.floor());
                if (room == 0 && waited_ms < 60'000) {
                    if (stop.load(std::memory_order_relaxed)) return;
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                    ++waited_ms;
                    continue;
                }
                if (room == 0) {
                    Log::get().line("warn", "33", "{} component {}: processor stalled for 60 s, dropping {} samples",
                                    st.code, rec.comp, s.size());
                    cs.ring.push_bulk(s);
                    return;
                }
                const auto take = std::min<std::size_t>(room, s.size());
                cs.ring.push_bulk(s.first(take));
                s = s.subspan(take);
            }
        };
        auto sink = [&](std::uint64_t p, std::span<const std::int32_t> s) {
            if (cs.base.load(std::memory_order_relaxed) == kUnset) cs.base.store(p, std::memory_order_release);
            push(s);
        };
        auto gap = [&](std::uint64_t a, std::uint64_t b) {
            cs.gaps.fetch_add(1, std::memory_order_relaxed);
            cs.gap_samples.fetch_add(b - a, std::memory_order_relaxed);
            for (std::uint64_t left = b - a; left > 0;) {
                const auto k = std::min<std::uint64_t>(left, gapfill.size());
                push(std::span<const std::int32_t>(gapfill.data(), k));
                left -= k;
            }
        };
        reorder[rec.comp]->offer(pos, std::span<const std::int32_t>(buf.data(), *n), sink, gap);
        cs.stale.store(reorder[rec.comp]->stats().stale, std::memory_order_relaxed);
        st.records.fetch_add(1, std::memory_order_relaxed);
    }
    st.ingest_done.store(true, std::memory_order_release);
}

}  // namespace ayzek::pipeline
