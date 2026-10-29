#pragma once

// SeedLink 3 client: the subset of the protocol a real-time client needs,
// written from the SeisComP description of it (docs/impl/19-seedlink.md).
//
// A session in multi-station mode:
//
//   HELLO                 two lines: server version, organisation
//   STATION sta net       OK or ERROR       } once per station
//   SELECT HH?.D          OK or ERROR       }
//   DATA [seq]            OK or ERROR       }
//   END                   no reply; the server starts sending packets
//
// A data packet is "SL", six hexadecimal digits of sequence number and one
// 512-byte miniSEED record. Sequence numbers are per station, 24 bit, and
// wrap after FFFFFF; DATA with a number asks for the packets from that number
// on, so a client resumes after packet s with DATA s+1. INFO requests, which
// are the only command allowed once data flows, are answered by packets
// headed "SLINFO" whose 8th byte is '*' while more of the response follows,
// or by "ERROR\r\n" in the packet stream if the server refuses them.
//
// The parsing is platform-independent. Client needs POSIX sockets and is not
// available in the Windows build.

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace ayzek::seedlink {

inline constexpr std::size_t kHeaderLength = 8;
inline constexpr std::size_t kRecordLength = 512;
inline constexpr std::size_t kPacketLength = kHeaderLength + kRecordLength;
inline constexpr std::uint32_t kSeqMask = 0xFFFFFF;

using Record = std::span<const std::byte, kRecordLength>;

// The sequence number after `s`.
[[nodiscard]] constexpr std::uint32_t seq_next(std::uint32_t s) noexcept {
  return (s + 1) & kSeqMask;
}

// True if `s` comes after `last` in the wrapping sequence: less than half the
// sequence space ahead of it.
[[nodiscard]] constexpr bool seq_after(std::uint32_t s,
                                       std::uint32_t last) noexcept {
  const std::uint32_t d = (s - last) & kSeqMask;
  return d != 0 && d <= kSeqMask / 2;
}

// One requested station: NET.STA, or NET.STA.LOC for one location code
// ("--" for the empty one).
struct StreamSpec {
  std::string net, sta;
  std::optional<std::string> loc; // unset: any location
  [[nodiscard]] std::string id() const; // "KO.DKL", "IU.ANMO.00"
  // The SELECT pattern: the HH? data records of the location.
  [[nodiscard]] std::string selector() const; // "HH?.D", "00HH?.D"
};

// Parses "NET.STA[.LOC][,NET.STA[.LOC]...]". Throws std::invalid_argument on
// a malformed entry or a station given twice.
std::vector<StreamSpec> parse_streams(std::string_view s);

struct Endpoint {
  std::string host, port;
};

// "HOST:PORT", "HOST" (port 18000) or "[IPV6]:PORT". Throws
// std::invalid_argument.
Endpoint parse_endpoint(std::string_view s);

// --- data-transfer phase ----------------------------------------------------

struct Frame {
  enum class Kind { Data, Info, Error, End };
  Kind kind = Kind::Data;
  std::uint32_t seq = 0; // Data
  bool more = false;     // Info: more packets of this response follow
  std::array<std::byte, kRecordLength> record{}; // Data and Info
};

class ProtocolError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

// Splits the bytes of the data-transfer phase into frames, however they are
// split by the network.
class FrameReader {
public:
  void feed(std::span<const std::byte> bytes);
  // The next complete frame, or nullopt until more bytes arrive. Throws
  // ProtocolError on bytes that begin no frame; the stream cannot be
  // resynchronised after that, so the connection has to be dropped.
  std::optional<Frame> next();
  void clear() noexcept;
  [[nodiscard]] std::size_t buffered() const noexcept {
    return buf_.size() - head_;
  }

private:
  std::vector<std::byte> buf_;
  std::size_t head_ = 0;
};

// --- handshake --------------------------------------------------------------

// The reply to a modifier command, without its CR LF: OK, or ERROR with
// whatever text the server added. nullopt for any other line.
struct Reply {
  bool ok = false;
  std::string text;
};
std::optional<Reply> parse_reply(std::string_view line);

std::string station_command(const StreamSpec &s); // "STATION DKL KO"
std::string select_command(const StreamSpec &s);  // "SELECT HH?.D"
// "DATA" (real time, from the next packet), or "DATA 00ABCD" to resume after
// packet 00ABCC.
std::string data_command(std::optional<std::uint32_t> last);

// Fields of a record's fixed header that do not depend on its byte order: the
// identifiers, trimmed, and the start time (BTIME, without blockette 1001's
// microseconds or the time correction) if it reads as a big-endian time.
struct RecordId {
  std::string net, sta, loc, cha;
  std::optional<std::int64_t> start_ns;
};
RecordId record_id(Record rec);

// --- client -----------------------------------------------------------------

struct ClientConfig {
  Endpoint server;
  std::vector<StreamSpec> streams;
  double timeout_s = 120;  // nothing received for this long: reconnect
  double keepalive_s = 60; // nothing received for this long: INFO ID; 0: never
  double first_delay_s = 1; // reconnection delay, doubling from this ...
  double max_delay_s = 30;  // ... up to this; back to the first after data
};

struct ClientStats {
  std::uint64_t connections = 0; // TCP connections made
  std::uint64_t sessions = 0;    // handshakes completed
  std::uint64_t packets = 0;     // data packets delivered
  std::uint64_t duplicates = 0;  // data packets discarded as already delivered
  std::uint64_t unknown = 0;     // data packets of no requested station
  std::uint64_t info_packets = 0;
  std::uint64_t keepalives = 0;    // INFO ID requests sent
  std::uint64_t server_errors = 0; // ERROR in the packet stream
  std::uint64_t seq_restarts = 0;  // sequence numbers that went back
};

class Client {
public:
  // `on_record(i, record)`: a data record of station `streams[i]`, in arrival
  // order and once only across reconnections. `on_idle()`: called at least
  // every 20 ms, also while connecting and while waiting to reconnect.
  using RecordFn = std::function<void(std::size_t, Record)>;
  using IdleFn = std::function<void()>;
  Client(ClientConfig cfg, RecordFn on_record, IdleFn on_idle);

  // Connects and receives until `stop()` returns true. A lost connection is
  // reconnected after the configured delay, and each station resumes after
  // its last packet. Throws std::runtime_error if the server refuses a
  // station, its selector or its DATA command in the first session (a wrong
  // stream list), and in the Windows build.
  void run(const std::function<bool()> &stop);

  [[nodiscard]] const ClientStats &stats() const noexcept { return stats_; }
  // False in the Windows build, which has no POSIX sockets.
  [[nodiscard]] static bool available() noexcept;

  // A data frame as run() delivers it: discarded if it is a duplicate, else
  // passed to on_record. False if it was discarded. Public for the tests.
  bool deliver(const Frame &f);

private:
  struct Station {
    std::optional<std::uint32_t> last_seq;
    std::map<std::string, std::int64_t> latest_start; // channel -> ns
    bool warned_restart = false;
  };

  void session(const std::function<bool()> &stop, bool &got_data);

  ClientConfig cfg_;
  RecordFn on_record_;
  IdleFn on_idle_;
  ClientStats stats_;
  std::vector<Station> stations_;
  std::map<std::string, std::size_t> index_; // "NET.STA" -> stream
  // The socket state, unused in the Windows build.
  [[maybe_unused]] bool handshaken_ = false; // a session was established before
  [[maybe_unused]] int fd_ = -1;
};

} // namespace ayzek::seedlink
