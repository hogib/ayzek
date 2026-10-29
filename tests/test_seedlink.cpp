// Unit tests for the SeedLink client's parsing (packets, handshake replies,
// stream lists, sequence numbers), its duplicate rule, and the live input's
// handling of each record: what it accepts, what it rejects and why, the
// archive, and the drop of records when a station is a ring behind. No
// network: the client's socket loop is tested end to end against
// tools/seedlink_mock.py by tests/seedlink_check.sh.
#include "pipeline/live.hpp"
#include "records.hpp"
#include "seedlink.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <print>
#include <string>
#include <vector>

#define CHECK(cond)                                                            \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::println(stderr, "CHECK failed: {}  ({}:{})", #cond, __FILE__,       \
                   __LINE__);                                                  \
      std::abort();                                                            \
    }                                                                          \
  } while (0)

using namespace ayzek;
using namespace ayzek::seedlink;

namespace {

template <typename F> bool throws(F f) {
  try {
    f();
  } catch (const std::exception &) {
    return true;
  }
  return false;
}

std::vector<std::byte> bytes(std::string_view s) {
  std::vector<std::byte> out;
  for (char c : s)
    out.push_back(static_cast<std::byte>(c));
  return out;
}

// A packet: the 8-byte header and the record.
std::vector<std::byte> packet(std::string_view head, const testrec::Bytes &r) {
  auto out = bytes(head);
  out.insert(out.end(), r.begin(), r.end());
  return out;
}

// The start-time fields of a record `e4` units of 1e-4 s after 18:00:00.
testrec::Spec at(testrec::Spec s, std::uint32_t e4) {
  s.hour = static_cast<std::uint8_t>(18 + e4 / 36'000'000);
  s.minute = static_cast<std::uint8_t>(e4 / 600'000 % 60);
  s.second = static_cast<std::uint8_t>(e4 / 10'000 % 60);
  s.fract = static_cast<std::uint16_t>(e4 % 10'000);
  return s;
}

Frame data_frame(std::uint32_t seq, const testrec::Bytes &r) {
  Frame f;
  f.seq = seq;
  f.record = r;
  return f;
}

} // namespace

int main() {
  // Sequence numbers: 24 bit, wrapping.
  CHECK(seq_next(0x00ABCC) == 0x00ABCD);
  CHECK(seq_next(0xFFFFFF) == 0);
  CHECK(seq_after(6, 5) && !seq_after(5, 5) && !seq_after(4, 5));
  CHECK(seq_after(0, 0xFFFFFF) && seq_after(3, 0xFFFFF0));
  CHECK(seq_after(0x800004, 5) && !seq_after(0x800005, 5));
  std::println("  sequence numbers         ok");

  // Stream lists.
  const auto s = parse_streams("KO.DKL,ko.kula,IU.ANMO.00,GE.WLF.--");
  CHECK(s.size() == 4);
  CHECK(s[0].net == "KO" && s[0].sta == "DKL" && !s[0].loc);
  CHECK(s[1].sta == "KULA" && s[1].id() == "KO.KULA");
  CHECK(s[2].loc == "00" && s[2].id() == "IU.ANMO.00");
  CHECK(s[3].loc && s[3].loc->empty() && s[3].id() == "GE.WLF.--");
  CHECK(s[0].selector() == "HH?.D" && s[2].selector() == "00HH?.D" &&
        s[3].selector() == "HH?.D");
  for (const char *bad : {"", "DKL", "KO.", ".DKL", "KO.TOOLONG", "KO.DKL.000",
                          "KO.D-L", "KO.DKL,TU.DKL", "KO.DKL,"})
    CHECK(throws([&] { parse_streams(bad); }));
  std::println("  stream lists             ok");

  // Server addresses.
  CHECK(parse_endpoint("localhost:18001").host == "localhost" &&
        parse_endpoint("localhost:18001").port == "18001");
  CHECK(parse_endpoint("geofon.gfz.de").port == "18000");
  CHECK(parse_endpoint("[::1]:18002").host == "::1" &&
        parse_endpoint("[::1]:18002").port == "18002");
  for (const char *bad : {"", ":18000", "host:", "host:18a", "[::1", "[::1]x"})
    CHECK(throws([&] { parse_endpoint(bad); }));
  std::println("  server addresses         ok");

  // Handshake replies and commands.
  CHECK(parse_reply("OK") && parse_reply("OK")->ok);
  CHECK(parse_reply("OK\r") && parse_reply("OK\r")->ok);
  CHECK(parse_reply("ERROR") && !parse_reply("ERROR")->ok);
  CHECK(parse_reply("ERROR station not found")->text == "station not found");
  CHECK(!parse_reply("SeedLink v3.1") && !parse_reply(""));
  CHECK(station_command(s[0]) == "STATION DKL KO");
  CHECK(select_command(s[2]) == "SELECT 00HH?.D");
  CHECK(data_command(std::nullopt) == "DATA");
  CHECK(data_command(0x00ABCC) == "DATA 00ABCD");
  CHECK(data_command(0xFFFFFF) == "DATA 000000");
  std::println("  replies and commands     ok");

  // Record identifiers and start time.
  testrec::Spec spec;
  spec.loc = "00";
  spec.cha = "HHN";
  const std::array<std::int32_t, 3> v{1, -2, 3};
  const auto rec = testrec::int32_record(at(spec, 12'3456), v);
  const RecordId id = record_id(rec);
  CHECK(id.net == "KO" && id.sta == "TST" && id.loc == "00" && id.cha == "HHN");
  CHECK(id.start_ns && *id.start_ns % 1'000'000'000 == 345'600'000);
  std::println("  record identifiers       ok");

  // Frames, however the bytes are split.
  {
    std::vector<std::byte> stream;
    for (const auto &part :
         {packet("SL00A1B2", rec), packet("SLINFO *", rec),
          packet("SLINFO  ", rec), bytes("ERROR\r\n"), packet("SLffffff", rec),
          bytes("END")})
      stream.insert(stream.end(), part.begin(), part.end());
    for (std::size_t step : {std::size_t{1}, std::size_t{7}, std::size_t{520},
                             stream.size()}) {
      FrameReader reader;
      std::vector<Frame> got;
      for (std::size_t i = 0; i < stream.size(); i += step) {
        reader.feed(std::span(stream).subspan(i, std::min(step, stream.size() - i)));
        while (auto f = reader.next())
          got.push_back(*f);
      }
      CHECK(got.size() == 6 && reader.buffered() == 0);
      CHECK(got[0].kind == Frame::Kind::Data && got[0].seq == 0x00A1B2 &&
            got[0].record == rec);
      CHECK(got[1].kind == Frame::Kind::Info && got[1].more);
      CHECK(got[2].kind == Frame::Kind::Info && !got[2].more);
      CHECK(got[3].kind == Frame::Kind::Error);
      CHECK(got[4].kind == Frame::Kind::Data && got[4].seq == 0xFFFFFF);
      CHECK(got[5].kind == Frame::Kind::End);
    }
    for (const auto &bad : {packet("SL00G1B2", rec), packet("SLINFO x", rec),
                            bytes("HELLO\r\n"), bytes("OK\r\n"),
                            bytes("ERROR" + std::string(300, ' '))}) {
      FrameReader reader;
      reader.feed(bad);
      CHECK(throws([&] { reader.next(); }));
    }
    FrameReader partial; // a packet cut short waits for the rest
    const auto p = packet("SL000001", rec);
    partial.feed(std::span(p).first(300));
    CHECK(!partial.next() && partial.buffered() == 300);
  }
  std::println("  frames                   ok");

  // The client's duplicate rule, without a connection.
  {
    ClientConfig cfg;
    cfg.streams = parse_streams("KO.TST");
    std::vector<std::uint32_t> got;
    std::uint32_t seq = 0;
    auto make = [&] {
      return Client(
          cfg,
          [&](std::size_t i, Record) {
            CHECK(i == 0);
            got.push_back(seq);
          },
          [] {});
    };
    Client client = make();
    auto deliver = [&](Client &c, std::uint32_t q, const char *cha,
                       std::uint32_t e4) {
      testrec::Spec r;
      r.cha = cha;
      seq = q;
      return c.deliver(data_frame(q, testrec::int32_record(at(r, e4), v)));
    };
    CHECK(deliver(client, 10, "HHZ", 0) && deliver(client, 11, "HHN", 0));
    CHECK(!deliver(client, 11, "HHN", 0)); // sent again
    CHECK(!deliver(client, 10, "HHZ", 0)); // an earlier one sent again
    CHECK(deliver(client, 12, "HHZ", 40'000));
    CHECK(client.stats().duplicates == 2 && client.stats().seq_restarts == 0);
    // The numbering goes back with new data: a restarted server.
    CHECK(deliver(client, 3, "HHZ", 80'000) && deliver(client, 4, "HHN", 80'000));
    CHECK(client.stats().seq_restarts == 1);
    CHECK(!deliver(client, 4, "HHN", 80'000));
    CHECK((got == std::vector<std::uint32_t>{10, 11, 12, 3, 4}));
    testrec::Spec other;
    other.sta = "XYZ";
    CHECK(!client.deliver(data_frame(1, testrec::int32_record(other, v))));
    CHECK(client.stats().unknown == 1 && client.stats().packets == 5 &&
          client.stats().duplicates == 3);

    // Across the wrap.
    got.clear();
    Client wrap = make();
    CHECK(deliver(wrap, 0xFFFFFE, "HHE", 0) && deliver(wrap, 0xFFFFFF, "HHE", 10'000));
    CHECK(deliver(wrap, 0, "HHE", 20'000) && deliver(wrap, 1, "HHE", 30'000));
    CHECK(!deliver(wrap, 0xFFFFFF, "HHE", 10'000));
    CHECK((got == std::vector<std::uint32_t>{0xFFFFFE, 0xFFFFFF, 0, 1}));
    CHECK(wrap.stats().seq_restarts == 0 && wrap.stats().duplicates == 1);
  }
  std::println("  duplicates and restarts  ok");

  // The live input's handling of each record.
  {
    const auto dir = std::filesystem::temp_directory_path() /
                     std::format("ayzek_test_seedlink_{}",
                                 std::chrono::steady_clock::now()
                                     .time_since_epoch()
                                     .count());
    pipeline::Station st;
    st.code = "TST";
    std::uint64_t received = 0;
    {
      pipeline::LiveConfig cfg;
      cfg.client.streams = parse_streams("KO.TST");
      cfg.archive_dir = dir.string();
      pipeline::LiveInput live(cfg, {&st});
      std::array<std::int32_t, 112> samples{};
      for (std::size_t i = 0; i < samples.size(); ++i)
        samples[i] = static_cast<std::int32_t>(i) - 50;
      auto offer = [&](testrec::Spec r, std::uint32_t e4) {
        live.on_record(0, testrec::int32_record(at(r, e4), samples), 0.0);
      };
      testrec::Spec r;
      for (const char *cha : {"HHZ", "HHN", "HHE"}) {
        r.cha = cha;
        offer(r, 0);
        offer(r, 11'200); // 112 samples = 1.12 s later: contiguous
      }
      CHECK(live.stats(0).accepted == 6 && live.stats(0).rejected.empty());
      for (const auto &c : st.comp) {
        CHECK(c.ring.written() == 224 && c.gaps.load() == 0);
        std::array<std::int32_t, 224> got{};
        CHECK(c.ring.copy_out(0, got));
        CHECK(got[0] == -50 && got[111] == 61 && got[112] == -50);
      }
      // 2025 day 314 = 2025-11-10, 18:00:00 UTC, plus 2 x 1.12 s.
      CHECK(std::abs(live.first_time() - (1762797600.0 + 1.12)) < 1e-6);
      CHECK(std::abs(live.last_time() - (1762797600.0 + 2.24)) < 1e-6);

      // Rejections, by reason.
      r.cha = "HH1";
      offer(r, 22'400);
      r.cha = "BHZ";
      offer(r, 22'400);
      r.cha = "HHZ";
      r.rate_factor = 200;
      offer(r, 22'400);
      r.rate_factor = 100;
      r.loc = "10"; // the first records were location ""
      offer(r, 22'400);
      r.loc = "";
      r.encoding = 4; // IEEE floats
      offer(r, 22'400);
      r.encoding = 3;
      const auto &rej = live.stats(0).rejected;
      CHECK(rej.size() == 5);
      std::uint64_t total = 0;
      for (const auto &[reason, n] : rej) {
        total += n;
        CHECK(n == 1);
      }
      CHECK(total == 5 && live.stats(0).accepted == 6);
      CHECK(rej.contains("unsupported encoding"));
      CHECK(std::ranges::any_of(rej, [](const auto &kv) {
        return kv.first.contains("HH1") && kv.first.contains("rotate");
      }));
      CHECK(std::ranges::any_of(
          rej, [](const auto &kv) { return kv.first.contains("200 Hz"); }));
      CHECK(std::ranges::any_of(rej, [](const auto &kv) {
        return kv.first.contains("location '10'");
      }));

      // A record whose data does not decode (Steim2 with a wrong Xn).
      testrec::Spec steim = at(r, 22'400);
      steim.encoding = 11;
      auto bad = testrec::header(steim, 4);
      testrec::put32(bad, 64, 0b01u << 24);   // word 3: four 8-bit diffs
      testrec::put32(bad, 68, 10);            // X0
      testrec::put32(bad, 72, 99);            // Xn, should be 16
      testrec::put32(bad, 76, 0x00010203);    // 0 1 2 3
      live.on_record(0, bad, 0.0);
      CHECK(st.undecodable.load() == 1 && live.stats(0).accepted == 6);

      // A station a ring behind: Z alone fills its ring (2^17 samples), then
      // the queue (kMaxQueuedSamples), then records are dropped.
      r.cha = "HHZ";
      std::uint32_t k = 2;
      while (live.stats(0).dropped < 10)
        offer(r, 11'200 * k++);
      const std::uint64_t dropped = live.stats(0).dropped;
      auto &z = st.comp[pipeline::Z];
      CHECK(z.ring.written() == pipeline::kRingCapacity);
      // The processor catches up; the next record follows a gap exactly as
      // long as the dropped records.
      for (int i = 0; i < 4; ++i) {
        z.ring.set_floor(z.ring.written());
        live.drain();
      }
      offer(r, 11'200 * k++);
      CHECK(z.gaps.load() == 1 && z.gap_samples.load() == dropped * 112);
      CHECK(z.end() == z.base.load() + k * 112);
      received = live.stats(0).received;
      CHECK(received == live.stats(0).accepted + dropped + total + 1);
    }
    // Every record received is in the archive, as received.
    CHECK(std::filesystem::file_size(dir / "TST.mseed") == 512 * received);
    std::filesystem::remove_all(dir);
  }
  std::println("  live records             ok");

  std::println("all seedlink tests passed");
}
