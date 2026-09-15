#pragma once

// A single-producer / single-consumer ring of samples, with history.
//
// This is deliberately NOT a queue. A queue pops: each item is consumed once
// and gone. Seismic windowing reads overlapping windows, so sample N appears
// in many of them, and the reader trails the writer at an arbitrary distance
// rather than draining it.
//
// The reader publishes a floor -- the oldest position it still needs -- and the
// writer will not overwrite at or above it, dropping new samples and counting
// them instead. That ordering is deliberate: it makes concurrent
// read-while-overwrite *impossible* rather than merely detectable.
//
// An earlier version let the writer always overwrite and re-checked afterwards
// whether the reader had been lapped. ThreadSanitizer correctly rejected it:
// reading a non-atomic object while another thread writes it is undefined
// behaviour, and no after-the-fact check repairs that. The test passed and
// reported zero torn reads, which is precisely how UB presents itself until it
// stops.
//
// Dropping is the right failure for this domain. Size the ring so a healthy
// reader is never near the floor; a non-zero `dropped()` then means the
// pipeline is not keeping up, which is an alarm rather than something to
// paper over.

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>

namespace ayzek {

// One instant of ground motion. Three components share a clock and every stage
// consumes all three, so they live together and travel under one index.
struct Sample {
    float z{}, n{}, e{};
};

enum class RingError : std::uint8_t {
    Expired,        // older than the ring retains
    NotYetWritten,
    Wrapped,        // the range straddles the physical end of storage
};

// `Capacity` must be a power of two so index arithmetic is a mask rather than
// a division. `std::has_single_bit` is C++20's is-power-of-two.
template <typename T, std::size_t Capacity>
    requires std::is_trivially_copyable_v<T>
class SpscRing {
    static_assert(std::has_single_bit(Capacity),
                  "Capacity must be a power of two");

    // `alignas(64)` puts each counter on its own cache line. Without it the
    // writer's store to `write_` and the reader's load of `read_` land in the
    // same line, and the two cores ping the line back and forth for no reason
    // -- false sharing. Same idea as Rust's `CachePadded`.
    alignas(64) std::atomic<std::uint64_t> write_{0};
    alignas(64) std::atomic<std::uint64_t> read_{0};      // reader's floor
    alignas(64) std::atomic<std::uint64_t> dropped_{0};
    alignas(64) std::array<T, Capacity> buf_{};

    static constexpr std::uint64_t kMask = Capacity - 1;

public:
    static constexpr std::size_t capacity() noexcept { return Capacity; }

    // --- producer side. Called from exactly one thread. ----------------------

    // Append one sample. Never blocks. Drops, and counts the drop, rather than
    // overwriting a slot the reader has reserved.
    bool push(const T& v) noexcept {
        const std::uint64_t w = write_.load(std::memory_order_relaxed);
        // `acquire` here pairs with the reader's `release` in `set_floor`, so
        // the floor we see is at least as new as the reader's last publish.
        const std::uint64_t floor = read_.load(std::memory_order_acquire);
        if (w - floor >= Capacity) {
            dropped_.fetch_add(1, std::memory_order_relaxed);
            return false;
        }
        buf_[w & kMask] = v;
        // `release` publishes the *data write above* to whichever thread does a
        // matching `acquire` on this counter. Rust spells it
        // `Ordering::Release`; C11 spells it `memory_order_release`. Without
        // it the compiler or CPU may make the counter visible before the
        // sample, and the reader sees uninitialised memory.
        write_.store(w + 1, std::memory_order_release);
        return true;
    }

    // Append many at once -- how packets actually arrive. Publishing the
    // counter once at the end costs one release instead of `n` of them, and
    // makes a packet appear atomically to the reader.
    std::size_t push_bulk(std::span<const T> src) noexcept {
        std::uint64_t w = write_.load(std::memory_order_relaxed);
        const std::uint64_t floor = read_.load(std::memory_order_acquire);
        const std::uint64_t room = Capacity - (w - floor);
        const std::size_t take = std::min<std::size_t>(src.size(), room);
        for (std::size_t i = 0; i < take; ++i) {
            buf_[(w + i) & kMask] = src[i];
        }
        if (take < src.size()) {
            dropped_.fetch_add(src.size() - take, std::memory_order_relaxed);
        }
        write_.store(w + take, std::memory_order_release);
        return take;
    }

    [[nodiscard]] std::uint64_t dropped() const noexcept {
        return dropped_.load(std::memory_order_relaxed);
    }

    // --- consumer side. Called from exactly one thread. ----------------------

    // Total samples ever written. Monotonic: it is a position in the stream,
    // not an index into storage. Storage indices come from masking it.
    //
    // Counting forever rather than wrapping is what makes "is this range still
    // live?" answerable. With wrapped indices, full and empty look identical.
    [[nodiscard]] std::uint64_t written() const noexcept {
        return write_.load(std::memory_order_acquire);
    }

    [[nodiscard]] std::uint64_t oldest_live() const noexcept {
        const std::uint64_t w = written();
        return w > Capacity ? w - Capacity : 0;
    }

    // Publish the oldest position still needed. The writer will not overwrite
    // at or above this, so anything below it may be reclaimed.
    //
    // This must be called *before* reading a range and must not be raised
    // while a view into that range is still held -- raising it is exactly what
    // authorises the writer to overwrite.
    void set_floor(std::uint64_t pos) noexcept {
        read_.store(pos, std::memory_order_release);
    }

    [[nodiscard]] std::uint64_t floor() const noexcept {
        return read_.load(std::memory_order_relaxed);
    }

    // A contiguous view of [from, from + count). No copy: this is a pointer
    // into the ring's own storage, like Rust's `&[T]`.
    //
    // Returns `std::expected`, which is C++23's `Result<T, E>` -- same idea,
    // same ergonomics, no exceptions involved.
    //
    // Fails with `Wrapped` when the range straddles the end of storage, since a
    // single span cannot describe two disjoint pieces. Callers who need to
    // cross that seam use `copy_out` instead.
    [[nodiscard]] std::expected<std::span<const T>, RingError>
    view(std::uint64_t from, std::size_t count) const noexcept {
        const std::uint64_t w = written();
        if (from + count > w) return std::unexpected(RingError::NotYetWritten);
        if (from < floor()) return std::unexpected(RingError::Expired);

        const std::size_t begin = static_cast<std::size_t>(from & kMask);
        if (begin + count > Capacity) return std::unexpected(RingError::Wrapped);
        return std::span<const T>(buf_.data() + begin, count);
    }

    // Same range, copied into caller storage, and therefore able to cross the
    // wrap seam. Costs a memcpy; use `view` when the range happens to be
    // contiguous and this when it is not.
    [[nodiscard]] std::expected<std::span<T>, RingError>
    copy_out(std::uint64_t from, std::span<T> dst) const noexcept {
        const std::size_t count = dst.size();
        const std::uint64_t w = written();
        if (from + count > w) return std::unexpected(RingError::NotYetWritten);
        // No post-copy re-check is needed, and none would have been sufficient:
        // the floor is what stops the writer touching these bytes at all.
        if (from < floor()) return std::unexpected(RingError::Expired);

        const std::size_t begin = static_cast<std::size_t>(from & kMask);
        const std::size_t first = std::min(count, Capacity - begin);
        std::copy_n(buf_.data() + begin, first, dst.data());
        if (first < count) {
            std::copy_n(buf_.data(), count - first, dst.data() + first);
        }
        return dst;
    }
};

}  // namespace ayzek
