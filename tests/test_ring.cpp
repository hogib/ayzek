// Tests for SpscRing: views, floors, wrap-around, bulk writes, and concurrent use.
#include "ring.hpp"

#include <print>
#include <thread>
#include <vector>
#include <cstdlib>

// Aborting check that, unlike assert(), stays active in release builds.
#define CHECK(cond)                                                                      \
    do {                                                                                 \
        if (!(cond)) {                                                                   \
            std::println(stderr, "CHECK failed: {}  ({}:{})", #cond, __FILE__, __LINE__); \
            std::abort();                                                                \
        }                                                                                \
    } while (0)

using namespace ayzek;

namespace {

Sample mk(int i) { return {float(i), float(i) + 0.5f, float(i) + 0.25f}; }

void basic_append_and_view() {
    SpscRing<Sample, 1024> r;
    for (int i = 0; i < 100; ++i) r.push(mk(i));
    CHECK(r.written() == 100);

    auto v = r.view(10, 20);
    CHECK(v.has_value());
    CHECK(v->size() == 20);
    CHECK((*v)[0].z == 10.0f);
    CHECK((*v)[19].z == 29.0f);
}

void reads_past_the_end_are_refused() {
    SpscRing<Sample, 1024> r;
    for (int i = 0; i < 10; ++i) r.push(mk(i));
    auto v = r.view(5, 20);                       // only 10 written
    CHECK(!v.has_value());
    CHECK(v.error() == RingError::NotYetWritten);
}

void writer_drops_rather_than_overwriting() {
    SpscRing<Sample, 64> r;
    // With the floor at 0 and 64 samples written, further writes must be dropped.
    for (int i = 0; i < 200; ++i) r.push(mk(i));
    CHECK(r.written() == 64);
    CHECK(r.dropped() == 136);

    // Raising the floor allows writes again.
    r.set_floor(64);
    CHECK(r.push(mk(1000)));
    CHECK(r.written() == 65);

    // Positions below the floor are no longer readable.
    auto gone = r.view(0, 8);
    CHECK(!gone.has_value());
    CHECK(gone.error() == RingError::Expired);
}

void copy_out_crosses_the_wrap_seam() {
    SpscRing<Sample, 64> r;
    // Raise the floor so that writing can continue past the end of the array.
    for (int i = 0; i < 100; ++i) {
        if (i >= 64) r.set_floor(std::uint64_t(i) - 63);
        r.push(mk(i));
    }

    // Positions 60..68 cross the end of the 64-element array, so view() fails.
    auto v = r.view(60, 8);
    CHECK(!v.has_value());
    CHECK(v.error() == RingError::Wrapped);

    // copy_out copies both parts.
    std::array<Sample, 8> dst{};
    auto c = r.copy_out(60, dst);
    CHECK(c.has_value());
    for (int i = 0; i < 8; ++i) CHECK(dst[size_t(i)].z == float(60 + i));

    // A range that does not cross the end is available as a view.
    auto flat = r.view(40, 8);
    CHECK(flat.has_value());
    CHECK((*flat)[0].z == 40.0f);
}

void bulk_matches_single() {
    SpscRing<Sample, 1024> a, b;
    std::vector<Sample> batch;
    for (int i = 0; i < 300; ++i) batch.push_back(mk(i));

    a.push_bulk(batch);
    for (const auto& s : batch) b.push(s);

    CHECK(a.written() == b.written());
    auto va = a.view(0, 300);
    auto vb = b.view(0, 300);
    CHECK(va && vb);
    for (size_t i = 0; i < 300; ++i) CHECK((*va)[i].z == (*vb)[i].z);
}

// Concurrent producer and consumer. Intended to be run under ThreadSanitizer
// (build-tsan) as well as normally.
void concurrent_producer_and_consumer() {
    constexpr int kTotal = 2'000'000;
    SpscRing<Sample, 8192> r;

    std::thread producer([&] {
        for (int i = 0; i < kTotal; ++i) {
            // Retry until there is room, so that every sample is written and a
            // corrupted read cannot be mistaken for a drop.
            while (!r.push(mk(i))) std::this_thread::yield();
        }
    });

    long reads = 0, torn = 0;
    std::thread consumer([&] {
        std::array<Sample, 64> dst{};
        std::uint64_t pos = 0;
        while (pos + 64 <= kTotal) {
            // Set the floor before reading, so the producer cannot overwrite
            // [pos, pos + 64) during the copy.
            r.set_floor(pos);
            auto c = r.copy_out(pos, dst);
            if (!c) { std::this_thread::yield(); continue; }
            // Each sample must equal the value written at its position.
            for (size_t i = 0; i < dst.size(); ++i) {
                if (dst[i].z != float(pos + i)) { ++torn; break; }
            }
            ++reads;
            pos += 64;
        }
        r.set_floor(kTotal);   // allow the producer to finish
    });

    producer.join();
    consumer.join();
    std::println("  concurrent: {} reads, {} torn, {} dropped", reads, torn, r.dropped());
    CHECK(torn == 0);
}

}  // namespace

int main() {
    basic_append_and_view();            std::println("  basic view                ok");
    reads_past_the_end_are_refused();   std::println("  refuse unwritten          ok");
    writer_drops_rather_than_overwriting();
    std::println("  drop, never overwrite     ok");
    copy_out_crosses_the_wrap_seam();   std::println("  copy across wrap          ok");
    bulk_matches_single();              std::println("  bulk == single            ok");
    concurrent_producer_and_consumer();
    std::println("all ring tests passed");
}
