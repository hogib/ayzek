// A reorderer that is subtly wrong emits a plausible stream with one record in
// the wrong place, so the checks here verify every sample against its position.

#include "reorder.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <print>
#include <random>
#include <span>
#include <utility>
#include <vector>

using namespace ayzek;

// assert() vanishes under NDEBUG; a test that stops checking in release builds
// is not a test.
#define CHECK(cond)                                                                      \
    do {                                                                                 \
        if (!(cond)) {                                                                   \
            std::println(stderr, "CHECK failed: {}  ({}:{})", #cond, __FILE__, __LINE__); \
            std::abort();                                                                \
        }                                                                                \
    } while (0)

namespace {

constexpr std::size_t kMaxRec = 735;          // the Steim2 decoder's per-record bound

// Every sample encodes its own position, so a misplaced one cannot pass.
std::int32_t value_at(std::uint64_t p) noexcept {
    return static_cast<std::int32_t>(static_cast<std::uint32_t>(p * 2654435761u));
}

struct Record {
    std::uint64_t pos;
    std::vector<std::int32_t> data;
};

Record make(std::uint64_t pos, std::size_t n) {
    Record r{pos, std::vector<std::int32_t>(n)};
    for (std::size_t i = 0; i < n; ++i) r.data[i] = value_at(pos + i);
    return r;
}

// Checks, on every call, that output is contiguous: each run or gap must begin
// exactly where the last one ended, and every emitted sample must be the one
// that belongs at its position.
struct Collector {
    static constexpr std::uint64_t kUnset = ~std::uint64_t{0};
    std::uint64_t first = kUnset, cursor = kUnset, emitted = 0;
    std::vector<std::pair<std::uint64_t, std::uint64_t>> gaps;
    bool ok = true;

    void run(std::uint64_t pos, std::span<const std::int32_t> s) {
        if (cursor == kUnset) first = cursor = pos;
        if (pos != cursor) ok = false;
        for (std::size_t i = 0; i < s.size(); ++i)
            if (s[i] != value_at(pos + i)) ok = false;
        cursor = pos + s.size();
        emitted += s.size();
    }
    void gap(std::uint64_t from, std::uint64_t to) {
        if (cursor == kUnset || from != cursor || to <= from) ok = false;
        gaps.emplace_back(from, to);
        cursor = to;
    }
};

template <std::size_t Held>
Collector feed(Reorderer<Held, kMaxRec>& r, const std::vector<Record>& arrivals) {
    Collector c;
    auto sink = [&](std::uint64_t p, std::span<const std::int32_t> s) { c.run(p, s); };
    auto gap = [&](std::uint64_t a, std::uint64_t b) { c.gap(a, b); };
    for (const Record& rec : arrivals) r.offer(rec.pos, rec.data, sink, gap);
    r.flush(sink, gap);
    return c;
}

std::vector<Record> contiguous(std::uint64_t base, std::initializer_list<std::size_t> lens) {
    std::vector<Record> v;
    std::uint64_t p = base;
    for (std::size_t n : lens) { v.push_back(make(p, n)); p += n; }
    return v;
}

// --- unit cases ------------------------------------------------------------

void in_order_passes_through() {
    auto recs = contiguous(1000, {300, 300, 300, 300});
    Reorderer<16, kMaxRec> r(500);
    auto c = feed(r, recs);
    CHECK(c.ok);
    CHECK(c.first == 1000 && c.cursor == 2200 && c.emitted == 1200);
    CHECK(c.gaps.empty());
    CHECK(r.stats().stale == 0);
}

void swapped_pair_is_reordered() {
    auto recs = contiguous(0, {200, 200, 200, 200, 200});
    std::swap(recs[1], recs[2]);
    Reorderer<16, kMaxRec> r(1000);
    auto c = feed(r, recs);
    CHECK(c.ok && c.gaps.empty() && c.emitted == 1000);
    CHECK(r.stats().held > 0);
}

void a_hole_never_filled_becomes_one_gap() {
    auto recs = contiguous(0, {200, 200, 200, 200, 200, 200, 200});
    recs.erase(recs.begin() + 2);                              // [400, 600) never arrives
    Reorderer<16, kMaxRec> r(300);
    auto c = feed(r, recs);
    CHECK(c.ok);
    CHECK(c.gaps.size() == 1);
    // Extra parentheses: the preprocessor splits macro arguments on commas and
    // does not treat <> as brackets, so the comma inside pair<...> would.
    CHECK((c.gaps[0] == std::pair<std::uint64_t, std::uint64_t>(400, 600)));
    CHECK(c.emitted == 1200);
}

void a_record_later_than_the_bound_is_stale_not_spliced() {
    auto recs = contiguous(0, {200, 200, 200, 200, 200, 200, 200});
    Record late = recs[2];
    recs.erase(recs.begin() + 2);
    recs.push_back(late);                                       // arrives after everything
    Reorderer<16, kMaxRec> r(300);
    auto c = feed(r, recs);
    CHECK(c.ok);
    CHECK(c.gaps.size() == 1);                                  // given up on...
    CHECK(r.stats().stale == 1);                                // ...and counted when it finally came
    CHECK(c.emitted == 1200);                                   // not spliced back in
}

void duplicates_are_dropped_and_counted() {
    auto recs = contiguous(0, {200, 200, 200});
    recs.push_back(recs[1]);
    Reorderer<16, kMaxRec> r(1000);
    auto c = feed(r, recs);
    CHECK(c.ok && c.emitted == 600);
    CHECK(r.stats().stale == 1);
}

void overlap_is_trimmed_not_duplicated() {
    std::vector<Record> recs{make(0, 300), make(250, 300)};     // 50 samples overlap
    Reorderer<16, kMaxRec> r(0);
    auto c = feed(r, recs);
    CHECK(c.ok && c.emitted == 550 && c.cursor == 550);
    CHECK(r.stats().trimmed_samples == 50);
}

void vertical_policy_never_waits() {
    Reorderer<16, kMaxRec> r(0);
    Collector c;
    auto sink = [&](std::uint64_t p, std::span<const std::int32_t> s) { c.run(p, s); };
    auto gap = [&](std::uint64_t a, std::uint64_t b) { c.gap(a, b); };

    auto a = make(0, 200), b = make(200, 200), d = make(600, 200), late = make(400, 200);
    r.offer(a.pos, a.data, sink, gap);
    r.offer(b.pos, b.data, sink, gap);
    CHECK(c.emitted == 400);                                    // no holding back
    r.offer(d.pos, d.data, sink, gap);
    CHECK(c.gaps.size() == 1 && c.emitted == 600);              // hole is a gap at once
    r.offer(late.pos, late.data, sink, gap);
    CHECK(r.stats().stale == 1 && c.emitted == 600);            // late Z is an anomaly
    CHECK(r.held() == 0 && c.ok);
}

void a_full_pool_forces_the_oldest_hole() {
    auto recs = contiguous(0, {100, 100, 100, 100, 100, 100, 100, 100, 100});
    recs.erase(recs.begin() + 1);                               // hole at [100, 200)
    Reorderer<4, kMaxRec> r(100'000);                           // lateness alone never gives up
    auto c = feed(r, recs);
    CHECK(c.ok);
    CHECK(c.gaps.size() == 1);
    CHECK(r.stats().forced_gaps == 1);
    CHECK(c.emitted == 800);
}

void a_late_opening_record_is_not_lost() {
    auto recs = contiguous(5000, {200, 200, 200, 200, 200, 200});
    std::swap(recs[0], recs[1]);                                // the true first arrives second
    Reorderer<16, kMaxRec> r(500);
    auto c = feed(r, recs);
    CHECK(c.ok && c.first == 5000 && c.emitted == 1200 && c.gaps.empty());
}

// --- property test -----------------------------------------------------------

// Arrival order is a sort on (end + jitter), jitter uniform in [0, J]. A record
// can then only follow records whose end is at most J beyond its own, so its
// lateness is bounded by J *by construction*, and a reorderer given that bound
// must reproduce the stream exactly.
void shuffled_streams_come_back_exactly(bool with_drops) {
    std::mt19937_64 rng(with_drops ? 0xA12EC : 0x5EED);
    const int trials = 2000;
    std::uint64_t total_samples = 0, total_gaps = 0, total_held = 0;

    for (int t = 0; t < trials; ++t) {
        std::uniform_int_distribution<std::size_t> len(50, 700), count(5, 60);
        std::uniform_int_distribution<std::uint64_t> jit_max(0, 2500);
        const std::uint64_t J = jit_max(rng);

        std::vector<Record> stream;
        std::uint64_t p = 1'000'000 + rng() % 1'000'000;
        const std::size_t nrec = count(rng);
        for (std::size_t i = 0; i < nrec; ++i) { stream.push_back(make(p, len(rng))); p += stream.back().data.size(); }

        // Drop a random interior subset. First and last are kept, so every
        // dropped run lies strictly inside the stream and must become a gap.
        std::vector<bool> drop(nrec, false);
        if (with_drops)
            for (std::size_t i = 1; i + 1 < nrec; ++i) drop[i] = (rng() % 5 == 0);

        std::vector<std::pair<std::uint64_t, std::size_t>> keyed;
        std::uniform_int_distribution<std::uint64_t> jit(0, J);
        for (std::size_t i = 0; i < nrec; ++i)
            if (!drop[i]) keyed.emplace_back(stream[i].pos + stream[i].data.size() + jit(rng), i);
        std::sort(keyed.begin(), keyed.end());

        std::vector<Record> arrivals;
        for (auto& [k, i] : keyed) arrivals.push_back(stream[i]);

        Reorderer<128, kMaxRec> r(J);
        auto c = feed(r, arrivals);

        // Expected gaps: maximal runs of dropped records.
        std::vector<std::pair<std::uint64_t, std::uint64_t>> expect;
        for (std::size_t i = 0; i < nrec; ++i) {
            if (!drop[i]) continue;
            const std::uint64_t a = stream[i].pos, b = a + stream[i].data.size();
            if (!expect.empty() && expect.back().second == a) expect.back().second = b;
            else expect.emplace_back(a, b);
        }
        std::uint64_t expect_emitted = 0;
        for (std::size_t i = 0; i < nrec; ++i) if (!drop[i]) expect_emitted += stream[i].data.size();

        const bool pass = c.ok && c.first == stream.front().pos && c.cursor == p &&
                          c.emitted == expect_emitted && c.gaps == expect &&
                          r.stats().stale == 0 && r.stats().forced_gaps == 0 &&
                          r.stats().max_lateness <= J;
        if (!pass) {
            std::println(stderr, "trial {} (drops={}) failed: J={} nrec={} ok={} emitted={}/{} gaps={}/{} "
                                 "stale={} forced={} max_lateness={}",
                         t, with_drops, J, nrec, c.ok, c.emitted, expect_emitted, c.gaps.size(),
                         expect.size(), r.stats().stale, r.stats().forced_gaps, r.stats().max_lateness);
            std::abort();
        }
        total_samples += c.emitted;
        total_gaps += c.gaps.size();
        total_held += r.stats().held;
    }
    std::println("  property ({})  {} trials, {} samples, {} held, {} gaps -- all exact",
                 with_drops ? "with drops" : "no drops  ", trials, total_samples, total_held, total_gaps);
}

}  // namespace

int main() {
    in_order_passes_through();                          std::println("  in order passes through           ok");
    swapped_pair_is_reordered();                        std::println("  swapped pair reordered            ok");
    a_hole_never_filled_becomes_one_gap();              std::println("  unfilled hole -> one gap          ok");
    a_record_later_than_the_bound_is_stale_not_spliced(); std::println("  too-late record stale, not spliced ok");
    duplicates_are_dropped_and_counted();               std::println("  duplicates dropped and counted    ok");
    overlap_is_trimmed_not_duplicated();                std::println("  overlap trimmed                   ok");
    vertical_policy_never_waits();                      std::println("  vertical policy never waits       ok");
    a_full_pool_forces_the_oldest_hole();               std::println("  full pool forces oldest hole      ok");
    a_late_opening_record_is_not_lost();                std::println("  late opening record kept          ok");
    shuffled_streams_come_back_exactly(false);
    shuffled_streams_come_back_exactly(true);
    std::println("all reorder tests passed");
}
