#pragma once

// Single-producer, single-consumer ring buffer that retains history.
//
// The consumer reads arbitrary ranges by absolute position instead of popping
// items, because consecutive analysis windows overlap and each sample is read
// several times.
//
// Overwrite protection: the consumer publishes a floor, the oldest position it
// may still read. The producer never writes over positions at or above the
// floor. If the ring is full, new samples are dropped and counted. A range the
// consumer is reading therefore cannot be modified concurrently. A non-zero
// `dropped()` means the consumer is not keeping up.

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

// Three-component sample. The pipeline keeps one ring per component; this type
// is used by the ring tests.
struct Sample {
    float z{}, n{}, e{};
};

enum class RingError : std::uint8_t {
    Expired,        // below the consumer's floor
    NotYetWritten,
    Wrapped,        // the range crosses the end of the storage array
};

// `Capacity` must be a power of two, so that a position maps to a storage index
// with a bit mask.
template <typename T, std::size_t Capacity>
    requires std::is_trivially_copyable_v<T>
class SpscRing {
    static_assert(std::has_single_bit(Capacity),
                  "Capacity must be a power of two");

    // Each counter is aligned to its own 64-byte cache line, so producer and
    // consumer updates do not invalidate each other's cache lines.
    alignas(64) std::atomic<std::uint64_t> write_{0};
    alignas(64) std::atomic<std::uint64_t> read_{0};      // consumer's floor
    alignas(64) std::atomic<std::uint64_t> dropped_{0};
    alignas(64) std::array<T, Capacity> buf_{};

    static constexpr std::uint64_t kMask = Capacity - 1;

public:
    static constexpr std::size_t capacity() noexcept { return Capacity; }

    // --- producer side (one thread) --------------------------------------------

    // Appends one sample. Returns false, and counts a drop, if the ring is full up
    // to the consumer's floor. Never blocks.
    bool push(const T& v) noexcept {
        const std::uint64_t w = write_.load(std::memory_order_relaxed);
        const std::uint64_t floor = read_.load(std::memory_order_acquire);
        if (w - floor >= Capacity) {
            dropped_.fetch_add(1, std::memory_order_relaxed);
            return false;
        }
        buf_[w & kMask] = v;
        // The counter is published after the sample is stored, so the consumer
        // never sees a position whose data has not been written.
        write_.store(w + 1, std::memory_order_release);
        return true;
    }

    // Appends as many samples as fit and returns that number. The counter is
    // published once, so the consumer sees the block appear as a whole.
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

    // --- consumer side (one thread) --------------------------------------------

    // Number of samples written since construction. Positions increase
    // monotonically and are not wrapped; the storage index is `position & kMask`.
    // An unwrapped count distinguishes a full ring from an empty one.
    [[nodiscard]] std::uint64_t written() const noexcept {
        return write_.load(std::memory_order_acquire);
    }

    // Oldest position still physically present in storage.
    [[nodiscard]] std::uint64_t oldest_live() const noexcept {
        const std::uint64_t w = written();
        return w > Capacity ? w - Capacity : 0;
    }

    // Publishes the oldest position the consumer may still read. Positions below
    // it may be overwritten. Raise it only after finishing with any view into the
    // range it releases.
    void set_floor(std::uint64_t pos) noexcept {
        read_.store(pos, std::memory_order_release);
    }

    [[nodiscard]] std::uint64_t floor() const noexcept {
        return read_.load(std::memory_order_relaxed);
    }

    // Zero-copy view of [from, from + count) into the ring's storage. Fails with
    // `Wrapped` if the range crosses the end of the array, since one span cannot
    // describe two pieces; `copy_out` handles that case.
    [[nodiscard]] std::expected<std::span<const T>, RingError>
    view(std::uint64_t from, std::size_t count) const noexcept {
        const std::uint64_t w = written();
        if (from + count > w) return std::unexpected(RingError::NotYetWritten);
        if (from < floor()) return std::unexpected(RingError::Expired);

        const std::size_t begin = static_cast<std::size_t>(from & kMask);
        if (begin + count > Capacity) return std::unexpected(RingError::Wrapped);
        return std::span<const T>(buf_.data() + begin, count);
    }

    // Copies [from, from + dst.size()) into `dst`, crossing the end of the array
    // if necessary.
    [[nodiscard]] std::expected<std::span<T>, RingError>
    copy_out(std::uint64_t from, std::span<T> dst) const noexcept {
        const std::size_t count = dst.size();
        const std::uint64_t w = written();
        if (from + count > w) return std::unexpected(RingError::NotYetWritten);
        // Positions at or above the floor are not written by the producer, so the
        // copy cannot race with a write.
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
