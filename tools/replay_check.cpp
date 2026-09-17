// Feeds an archive file through Reorderer and evaluates the output.
//
//   replay_check FILE.mseed
//
// Records are offered per channel in file order, which contains the archive's
// out-of-order records. For each channel and max_lateness setting it reports:
//
//   exact     whether the output equals the time-sorted records sample for
//             sample, with gaps exactly at the holes in the data
//   stale     records discarded as too late
//   delay     how far the stream advanced between a record's arrival and its
//             emission, i.e. the latency added by waiting
//   gaps      the length distribution of the gaps in the output
//
// The archive has no arrival times. Arrival is modelled on a stream clock (the
// newest sample position seen on the channel), assuming file order is arrival
// order. A final table estimates detector windows lost when horizontal gaps up
// to a given length are bridged.

#include "mseed.hpp"
#include "reorder.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <map>
#include <print>
#include <string>
#include <vector>

using namespace ayzek;

namespace {

constexpr std::uint64_t kPeriodNs = 10'000'000; // 100 Hz

struct Rec {
  std::uint64_t pos;
  std::vector<std::int32_t> s;
};

struct Holes {
  std::vector<std::pair<std::uint64_t, std::uint64_t>> gaps;
  std::uint64_t samples = 0;
};

// Holes between records sorted by position.
Holes reference_holes(const std::vector<const Rec *> &sorted) {
  Holes h;
  std::uint64_t end = sorted.front()->pos;
  for (const Rec *r : sorted) {
    if (r->pos > end)
      h.gaps.emplace_back(end, r->pos);
    end = std::max(end, r->pos + r->s.size());
    h.samples += r->s.size();
  }
  return h;
}

double pct(std::vector<std::uint64_t> &v, double q) {
  if (v.empty())
    return 0.0;
  const auto k =
      static_cast<std::size_t>(q * static_cast<double>(v.size() - 1));
  std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(k),
                   v.end());
  return static_cast<double>(v[k]) / 100.0; // samples -> seconds
}

using Gaps = std::vector<std::pair<std::uint64_t, std::uint64_t>>;

// Fraction of detector window start positions lost to gaps when horizontal gaps
// up to a given length are bridged and longer ones are not. A window is lost if
// it overlaps an unbridged gap on any channel, so blocked intervals are merged
// across channels. Vertical-component gaps are never bridged.
void fill_policy(const std::map<std::string, Gaps> &gaps,
                 std::uint64_t span_begin, std::uint64_t span_end) {
  constexpr std::uint64_t kWindow = 600; // the detector's 6 s at 100 Hz
  const double starts = static_cast<double>(span_end - span_begin - kWindow);

  std::println("bridging horizontal gaps up to a length, skipping longer ones "
               "({} s windows):",
               kWindow / 100);
  for (std::uint64_t fill :
       {0, 100, 200, 300, 400, 500, 800, 1000, 1500, 9000}) {
    Gaps blocked; // window starts that would overlap a skipped gap
    std::uint64_t bridged = 0, skipped = 0, fabricated = 0;
    for (const auto &[cha, gs] : gaps) {
      for (auto [a, b] : gs) {
        if (cha != "HHZ" && b - a <= fill) {
          ++bridged;
          fabricated += b - a;
          continue;
        }
        ++skipped;
        blocked.emplace_back(
            std::max(span_begin, (a + 1 > kWindow) ? a + 1 - kWindow : 0), b);
      }
    }
    std::ranges::sort(blocked);
    std::uint64_t lost = 0, cur_a = 0, cur_b = 0;
    for (auto [a, b] : blocked) {
      if (a > cur_b) {
        lost += cur_b - cur_a;
        cur_a = a;
      }
      cur_b = std::max(cur_b, b);
    }
    lost += cur_b - cur_a;

    std::println("  up to {:>4.1f}s   bridged {:>4}  skipped {:>4}   windows "
                 "lost {:>5.3f}% ({:>5.0f} s)"
                 "   fabricated {:>5.0f} s",
                 static_cast<double>(fill) / 100.0, bridged, skipped,
                 100.0 * static_cast<double>(lost) / starts,
                 static_cast<double>(lost) / 100.0,
                 static_cast<double>(fabricated) / 100.0);
  }
}

Gaps run_channel(const std::string &cha, const std::vector<Rec> &recs,
                 std::uint64_t lateness) {
  // Record pointers sorted by position, used for the reference holes and to
  // look up the expected value of every emitted sample without copying data.
  std::vector<const Rec *> sorted;
  sorted.reserve(recs.size());
  for (const Rec &r : recs)
    sorted.push_back(&r);
  std::sort(sorted.begin(), sorted.end(),
            [](auto *a, auto *b) { return a->pos < b->pos; });
  const Holes ref = reference_holes(sorted);

  Reorderer<256, mseed::kMaxDiffs> r(lateness);
  std::uint64_t clock = 0, emitted = 0, mismatches = 0,
                cursor = ~std::uint64_t{0};
  std::size_t lookup = 0;
  std::vector<std::pair<std::uint64_t, std::uint64_t>> gaps;
  std::map<std::uint64_t, std::uint64_t>
      arrived_at;                    // record start -> clock at arrival
  std::vector<std::uint64_t> delays; // per released run, in samples
  std::uint64_t delayed_samples = 0;
  bool contiguous = true;

  auto sink = [&](std::uint64_t pos, std::span<const std::int32_t> s) {
    if (cursor != ~std::uint64_t{0} && pos != cursor)
      contiguous = false;
    cursor = pos + s.size();

    auto it = arrived_at.upper_bound(pos);
    const std::uint64_t arrival =
        (it == arrived_at.begin()) ? clock : std::prev(it)->second;
    const std::uint64_t delay = clock - arrival;
    delays.push_back(delay);
    if (delay > 0)
      delayed_samples += s.size();

    for (std::size_t i = 0; i < s.size(); ++i) {
      const std::uint64_t p = pos + i;
      while (lookup < sorted.size() &&
             sorted[lookup]->pos + sorted[lookup]->s.size() <= p)
        ++lookup;
      if (lookup == sorted.size()) {
        ++mismatches;
        continue;
      } // emitted past the last record
      const Rec *home = sorted[lookup];
      if (home->pos > p || home->s[p - home->pos] != s[i])
        ++mismatches;
    }
    emitted += s.size();
  };
  auto gap = [&](std::uint64_t a, std::uint64_t b) {
    if (cursor != a)
      contiguous = false;
    gaps.emplace_back(a, b);
    cursor = b;
  };

  for (const Rec &rec : recs) {
    clock = std::max(clock, rec.pos + rec.s.size());
    arrived_at.emplace(rec.pos, clock);
    r.offer(rec.pos, rec.s, sink, gap);
  }
  r.flush(sink, gap);

  // Lengths of all gaps in the output, whether caused by missing or late data.
  std::vector<std::uint64_t> gap_len;
  for (auto [a, b] : gaps)
    gap_len.push_back(b - a);
  const auto under = [&](std::uint64_t n) {
    return std::ranges::count_if(gap_len, [n](auto g) { return g <= n; });
  };

  const auto &st = r.stats();
  const bool exact = contiguous && mismatches == 0 && emitted == ref.samples &&
                     gaps == ref.gaps;
  std::uint64_t total = emitted ? emitted : 1;
  std::println(
      "{}  L={:>5.1f}s  rec {:>7}  held {:>5}  stale {:>4}  gaps {:>3}/{:<3}  "
      "delayed {:>5.2f}%  delay p50 {:>4.1f}s p99 {:>4.1f}s max {:>4.1f}s  "
      "max_late {:>4.1f}s  {}",
      cha, static_cast<double>(lateness) / 100.0, st.records, st.held, st.stale,
      gaps.size(), ref.gaps.size(),
      100.0 * static_cast<double>(delayed_samples) / static_cast<double>(total),
      pct(delays, 0.50), pct(delays, 0.99), pct(delays, 1.0),
      static_cast<double>(st.max_lateness) / 100.0,
      exact ? "EXACT" : (st.stale ? "lost data" : "MISMATCH"));
  std::println("     gap length min {:>5.2f}s p50 {:>5.1f}s p90 {:>5.1f}s max "
               "{:>6.1f}s   <=1s {}  <=5s {}  <=10s {}  <=30s {}",
               pct(gap_len, 0.0), pct(gap_len, 0.50), pct(gap_len, 0.90),
               pct(gap_len, 1.0), under(100), under(500), under(1000),
               under(3000));
  return gaps;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::println(stderr, "usage: {} FILE.mseed", argv[0]);
    return 2;
  }

  std::ifstream f(argv[1], std::ios::binary);
  if (!f) {
    std::println(stderr, "cannot open {}", argv[1]);
    return 2;
  }

  std::map<std::string, std::vector<Rec>> by_channel;
  std::array<std::byte, mseed::kRecordLength> raw{};
  std::array<std::int32_t, mseed::kMaxDiffs> buf{};
  std::uint64_t misaligned = 0;

  while (f.read(reinterpret_cast<char *>(raw.data()), raw.size())) {
    auto h = mseed::parse_header(raw);
    if (!h) {
      std::println(stderr, "header: {}", mseed::to_string(h.error()));
      return 1;
    }
    auto n = mseed::decode(raw, *h, buf);
    if (!n) {
      std::println(stderr, "decode: {}", mseed::to_string(n.error()));
      return 1;
    }
    if (h->start_ns < 0 || h->start_ns % kPeriodNs != 0)
      ++misaligned;
    const auto pos = static_cast<std::uint64_t>(h->start_ns) / kPeriodNs;
    by_channel[std::string(h->cha())].push_back(
        {pos, {buf.begin(), buf.begin() + *n}});
  }
  std::println("decoded; {} records off the 10 ms grid\n", misaligned);

  std::map<std::string, Gaps> realtime; // gaps per channel at max_lateness = 0
  std::uint64_t span_begin = ~std::uint64_t{0}, span_end = 0;
  for (const auto &[cha, recs] : by_channel) {
    for (const Rec &rec : recs) {
      span_begin = std::min(span_begin, rec.pos);
      span_end = std::max(span_end, rec.pos + rec.s.size());
    }
  }

  for (const auto &[cha, recs] : by_channel) {
    const std::array<std::uint64_t, 6> settings =
        cha == "HHZ"
            ? std::array<std::uint64_t, 6>{0, 0, 0, 0, 0, 0}
            : std::array<std::uint64_t, 6>{0, 200, 400, 600, 800, 1000};
    std::uint64_t last = ~std::uint64_t{0};
    for (std::uint64_t L : settings) {
      if (L == last)
        continue;
      last = L;
      Gaps g = run_channel(cha, recs, L);
      if (L == 0)
        realtime[cha] = std::move(g);
    }
    std::println("");
  }
  fill_policy(realtime, span_begin, span_end);
}
