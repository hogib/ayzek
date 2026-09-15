// A ring that is subtly wrong drops samples silently, so these are the point.
#include "ring.hpp"

#include <cassert>
#include <print>
#include <thread>
#include <vector>

using namespace ayzek;

namespace {

Sample mk(int i) { return {float(i), float(i) + 0.5f, float(i) + 0.25f}; }

void basic_append_and_view() {
    SpscRing<Sample, 1024> r;
    for (int i = 0; i < 100; ++i) r.push(mk(i));
    assert(r.written() == 100);

    auto v = r.view(10, 20);
    assert(v.has_value());
    assert(v->size() == 20);
    assert((*v)[0].z == 10.0f);
    assert((*v)[19].z == 29.0f);
}

void reads_past_the_end_are_refused() {
    SpscRing<Sample, 1024> r;
    for (int i = 0; i < 10; ++i) r.push(mk(i));
    auto v = r.view(5, 20);                       // only 10 written
    assert(!v.has_value());
    assert(v.error() == RingError::NotYetWritten);
}

void writer_drops_rather_than_overwriting() {
    SpscRing<Sample, 64> r;
    // Floor at 0 means the reader still wants everything, so once 64 samples
    // are in, the writer must refuse rather than overwrite.
    for (int i = 0; i < 200; ++i) r.push(mk(i));
    assert(r.written() == 64);
    assert(r.dropped() == 136);

    // Releasing the floor lets the writer proceed again.
    r.set_floor(64);
    assert(r.push(mk(1000)));
    assert(r.written() == 65);

    // And the released range is now refused to the reader, as it may have
    // been reclaimed.
    auto gone = r.view(0, 8);
    assert(!gone.has_value());
    assert(gone.error() == RingError::Expired);
}

void copy_out_crosses_the_wrap_seam() {
    SpscRing<Sample, 64> r;
    // Keep the floor moving so the writer has room to pass the seam.
    for (int i = 0; i < 100; ++i) {
        if (i >= 64) r.set_floor(std::uint64_t(i) - 63);
        r.push(mk(i));
    }

    // Position 60 sits at storage index 60, so 60..68 runs off the physical
    // end of the array. A single span cannot describe two disjoint pieces.
    auto v = r.view(60, 8);
    assert(!v.has_value());
    assert(v.error() == RingError::Wrapped);

    // copy_out stitches the two pieces together.
    std::array<Sample, 8> dst{};
    auto c = r.copy_out(60, dst);
    assert(c.has_value());
    for (int i = 0; i < 8; ++i) assert(dst[size_t(i)].z == float(60 + i));

    // A range that happens not to straddle the seam still works as a view.
    auto flat = r.view(40, 8);
    assert(flat.has_value());
    assert((*flat)[0].z == 40.0f);
}

void bulk_matches_single() {
    SpscRing<Sample, 1024> a, b;
    std::vector<Sample> batch;
    for (int i = 0; i < 300; ++i) batch.push_back(mk(i));

    a.push_bulk(batch);
    for (const auto& s : batch) b.push(s);

    assert(a.written() == b.written());
    auto va = a.view(0, 300);
    auto vb = b.view(0, 300);
    assert(va && vb);
    for (size_t i = 0; i < 300; ++i) assert((*va)[i].z == (*vb)[i].z);
}

// The real check: a producer and a consumer hammering it at once. Run under
// ThreadSanitizer to make it meaningful.
void concurrent_producer_and_consumer() {
    constexpr int kTotal = 2'000'000;
    SpscRing<Sample, 8192> r;

    std::thread producer([&] {
        for (int i = 0; i < kTotal; ++i) {
            // Spin until the reader releases room. Dropping here would be
            // correct in production; for the test we want every sample, so
            // that a torn read cannot hide as a drop.
            while (!r.push(mk(i))) std::this_thread::yield();
        }
    });

    long reads = 0, torn = 0;
    std::thread consumer([&] {
        std::array<Sample, 64> dst{};
        std::uint64_t pos = 0;
        while (pos + 64 <= kTotal) {
            // Publish the floor before reading: this is what forbids the
            // writer from touching [pos, pos+64) while the copy runs.
            r.set_floor(pos);
            auto c = r.copy_out(pos, dst);
            if (!c) { std::this_thread::yield(); continue; }
            // Every sample must be exactly the one written at that position.
            // A torn read shows up here as a discontinuity.
            for (size_t i = 0; i < dst.size(); ++i) {
                if (dst[i].z != float(pos + i)) { ++torn; break; }
            }
            ++reads;
            pos += 64;
        }
        r.set_floor(kTotal);   // release everything so the producer can finish
    });

    producer.join();
    consumer.join();
    std::println("  concurrent: {} reads, {} torn, {} dropped", reads, torn, r.dropped());
    assert(torn == 0);
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
