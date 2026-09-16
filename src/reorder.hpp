#pragma once

// Bounded reordering of out-of-order records into a contiguous sample stream.
//
// About 0.3% of horizontal records in the TDVMS archive arrive late -- by 7.7 s
// to ~10 min, median ~40 s -- and the continuous IIR filter downstream needs
// samples in time order. This sits
// between decode and the ring: records go in in arrival order, contiguous runs
// come out in strictly increasing position.
//
// Positions are exact. Every record in the archive starts on a 10 ms boundary,
// so position = start_ns / 10'000'000 with no rounding, and a hole is simply a
// jump in position.
//
// Policy, all of it loud:
//   - a record at the frontier is emitted immediately, without copying
//   - a record ahead of the frontier is held until the hole before it fills
//   - a hole left open for more than `max_lateness` samples of later data is
//     declared a gap and the frontier moves past it
//   - a record that arrives after its hole was given up on is counted as stale,
//     never spliced in after the fact -- splicing would rewrite samples the
//     filter has already consumed
//
// Every channel uses max_lateness = 0: never wait, any hole is an immediate
// gap, any late record is counted as stale. Data late enough to be late here is
// too late for early warning (docs/DESIGN.md), so the class earns its place by
// declaring gaps, trimming overlaps and catching duplicates. A nonzero bound is
// supported for offline replay.

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>

namespace ayzek {

struct ReorderStats {
    std::uint64_t records = 0;
    std::uint64_t in_order = 0;          // arrived exactly at the frontier
    std::uint64_t held = 0;              // arrived early and waited for a hole to fill
    std::uint64_t stale = 0;             // wholly behind the frontier, dropped
    std::uint64_t trimmed_samples = 0;   // leading samples behind the frontier, dropped
    std::uint64_t gaps = 0;
    std::uint64_t gap_samples = 0;
    std::uint64_t forced_gaps = 0;       // declared early because the pool was full
    std::uint64_t rejected = 0;          // longer than a slot can hold
    // How far a record's end trailed the newest data already received, in
    // samples. The same quantity the constructor's max_lateness bounds.
    std::uint64_t max_lateness = 0;
};

// `MaxHeld` records can wait at once; `MaxRecordSamples` bounds one record.
// Storage is a fixed array, so offer() never allocates.
template <std::size_t MaxHeld, std::size_t MaxRecordSamples>
class Reorderer {
    struct Slot {
        std::uint64_t pos = 0;
        std::uint32_t n = 0;
        bool used = false;
        std::array<std::int32_t, MaxRecordSamples> data{};
    };

    std::array<Slot, MaxHeld> pool_{};
    std::size_t count_ = 0;
    std::uint64_t max_lateness_;

    // The frontier is unset until the stream has shown enough data to be sure
    // which record is really first. Starting it at the first *arrival* would
    // drop the start of any stream whose opening records were themselves late.
    bool started_ = false;
    std::uint64_t next_ = 0;             // position of the next sample to emit
    std::uint64_t seen_end_ = 0;         // highest record end ever offered

    ReorderStats stats_{};

public:
    // `max_lateness_samples`: a record whose end trails the newest data already
    // received by more than this is given up on, and its hole becomes a gap.
    // Zero means never wait, which is the real-time policy on every channel.
    explicit Reorderer(std::uint64_t max_lateness_samples) noexcept
        : max_lateness_(max_lateness_samples) {}

    [[nodiscard]] const ReorderStats& stats() const noexcept { return stats_; }
    [[nodiscard]] std::size_t held() const noexcept { return count_; }
    [[nodiscard]] std::uint64_t frontier() const noexcept { return next_; }

    // `sink(pos, samples)` receives contiguous runs in increasing position;
    // `gap(from, to)` marks [from, to) as missing.
    //
    // These are templates constrained with `std::invocable`, which is how C++
    // spells Rust's `impl FnMut(u64, &[i32])`: the call is monomorphised and
    // inlined, with no allocation and no indirection. `std::function` would be
    // the `Box<dyn Fn>` equivalent, and it can allocate.
    template <typename Sink, typename Gap>
        requires std::invocable<Sink&, std::uint64_t, std::span<const std::int32_t>> &&
                 std::invocable<Gap&, std::uint64_t, std::uint64_t>
    void offer(std::uint64_t pos, std::span<const std::int32_t> samples, Sink&& sink, Gap&& gap) {
        const std::uint64_t n = samples.size();
        if (n == 0) return;
        ++stats_.records;

        const std::uint64_t end = pos + n;
        if (seen_end_ > end) stats_.max_lateness = std::max(stats_.max_lateness, seen_end_ - end);
        seen_end_ = std::max(seen_end_, end);

        if (started_) {
            if (end <= next_) {                              // entirely behind: too late
                ++stats_.stale;
                return;
            }
            if (pos <= next_) {                              // at, or straddling, the frontier
                const std::uint64_t skip = next_ - pos;
                stats_.trimmed_samples += skip;
                if (skip == 0) ++stats_.in_order;
                sink(next_, samples.subspan(static_cast<std::size_t>(skip)));
                next_ = end;
                drain(sink);
                resolve(sink, gap, false);
                return;
            }
        }

        // Ahead of the frontier, or frontier not yet established: hold it.
        if (n > MaxRecordSamples) {
            // Unreachable while MaxRecordSamples matches the decoder's bound,
            // but a record that cannot be stored is counted, not lost silently.
            ++stats_.rejected;
            return;
        }
        if (count_ == MaxHeld) resolve(sink, gap, true);
        Slot& s = free_slot();
        s.pos = pos;
        s.n = static_cast<std::uint32_t>(n);
        s.used = true;
        std::copy(samples.begin(), samples.end(), s.data.begin());
        ++count_;
        ++stats_.held;

        resolve(sink, gap, false);
    }

    // End of stream: release everything, declaring remaining holes as gaps.
    template <typename Sink, typename Gap>
    void flush(Sink&& sink, Gap&& gap) {
        while (count_ > 0) {
            if (!started_) start_at_earliest();
            Slot* m = earliest();
            if (m->pos > next_) {
                declare_gap(gap, m->pos);
            }
            drain(sink);
        }
    }

private:
    Slot& free_slot() noexcept {
        for (Slot& s : pool_) if (!s.used) return s;
        return pool_[0];   // guarded by the caller's count_ check
    }

    // Linear scan. MaxHeld is tens of slots, and a scan over contiguous memory
    // beats a heap's pointer-chasing at that size.
    Slot* earliest() noexcept {
        Slot* best = nullptr;
        for (Slot& s : pool_)
            if (s.used && (best == nullptr || s.pos < best->pos)) best = &s;
        return best;
    }

    void release(Slot& s) noexcept {
        s.used = false;
        --count_;
    }

    void start_at_earliest() noexcept {
        next_ = earliest()->pos;
        started_ = true;
    }

    template <typename Gap>
    void declare_gap(Gap& gap, std::uint64_t to) {
        gap(next_, to);
        ++stats_.gaps;
        stats_.gap_samples += to - next_;
        next_ = to;
    }

    // Emit every held record now reachable from the frontier, in order.
    template <typename Sink>
    void drain(Sink& sink) {
        for (;;) {
            Slot* m = earliest();
            if (m == nullptr || m->pos > next_) return;
            const std::uint64_t end = m->pos + m->n;
            if (end <= next_) {
                ++stats_.stale;                              // duplicate of emitted data
            } else {
                const std::uint64_t skip = next_ - m->pos;
                stats_.trimmed_samples += skip;
                sink(next_, std::span<const std::int32_t>(m->data.data() + skip, m->n - skip));
                next_ = end;
            }
            release(*m);
        }
    }

    // Give up on holes that have been open too long, or on the oldest hole when
    // the pool has no room.
    template <typename Sink, typename Gap>
    void resolve(Sink& sink, Gap& gap, bool pool_full) {
        while (count_ > 0) {
            if (!started_) {
                const std::uint64_t first = earliest()->pos;
                if (!pool_full && seen_end_ - first <= max_lateness_) return;
                start_at_earliest();
                drain(sink);
                continue;
            }
            Slot* m = earliest();
            // Measured from where the hole *ends* (the next held record), not
            // where it starts. From the start, a hole's own length counts as
            // lateness: a 3.6 s record arriving 8 s late would read as ~11.6 s
            // and be declared a gap despite arriving in time. From the end, this
            // is the same quantity as stats().max_lateness -- how far the missing
            // data's end trails the newest data received -- so the measured
            // maximum on real data is directly the setting to use.
            const bool waited_out = seen_end_ - m->pos > max_lateness_;
            if (!waited_out && !pool_full) return;
            if (m->pos > next_) {
                if (pool_full && !waited_out) ++stats_.forced_gaps;
                declare_gap(gap, m->pos);
            }
            drain(sink);
            pool_full = false;                               // one slot is enough
        }
    }
};

}  // namespace ayzek
