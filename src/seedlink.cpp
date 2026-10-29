#include "seedlink.hpp"

#include "mseed.hpp"
#include "pipeline/common.hpp"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <format>
#include <thread>

#ifndef _WIN32
#include <cerrno>
#include <csignal>
#include <fcntl.h>
#include <netdb.h>
#include <netinet/in.h>
#include <poll.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>
#endif

namespace ayzek::seedlink {

using pipeline::Log;

namespace {

bool is_code(std::string_view s, std::size_t max) {
  return !s.empty() && s.size() <= max &&
         std::ranges::all_of(s, [](char c) {
           return (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9');
         });
}

std::string upper(std::string_view s) {
  std::string out(s);
  for (char &c : out)
    if (c >= 'a' && c <= 'z')
      c = static_cast<char>(c - 'a' + 'A');
  return out;
}

int hex_digit(std::byte b) {
  const auto c = static_cast<char>(b);
  if (c >= '0' && c <= '9')
    return c - '0';
  if (c >= 'A' && c <= 'F')
    return c - 'A' + 10;
  if (c >= 'a' && c <= 'f')
    return c - 'a' + 10;
  return -1;
}

std::string trimmed(Record rec, std::size_t off, std::size_t len) {
  std::string s(reinterpret_cast<const char *>(rec.data()) + off, len);
  while (!s.empty() && (s.back() == ' ' || s.back() == '\0'))
    s.pop_back();
  return s;
}

} // namespace

std::string StreamSpec::id() const {
  return loc ? std::format("{}.{}.{}", net, sta, loc->empty() ? "--" : *loc)
             : std::format("{}.{}", net, sta);
}

std::string StreamSpec::selector() const {
  // The empty location has no selector spelling that every server reads the
  // same way, so it is selected client-side (live.cpp).
  return loc && !loc->empty() ? *loc + "HH?.D" : std::string("HH?.D");
}

std::vector<StreamSpec> parse_streams(std::string_view s) {
  if (s.empty())
    throw std::invalid_argument("no streams given");
  std::vector<StreamSpec> out;
  for (std::size_t from = 0; from <= s.size();) {
    const auto comma = std::min(s.find(',', from), s.size());
    const std::string item = upper(s.substr(from, comma - from));
    from = comma + 1;
    std::vector<std::string> parts;
    for (std::size_t at = 0;;) {
      const auto dot = item.find('.', at);
      parts.push_back(item.substr(at, dot - at));
      if (dot == std::string::npos)
        break;
      at = dot + 1;
    }
    if (parts.size() < 2 || parts.size() > 3 || !is_code(parts[0], 2) ||
        !is_code(parts[1], 5) ||
        (parts.size() == 3 && parts[2] != "--" && !is_code(parts[2], 2)))
      throw std::invalid_argument("bad stream '" + item +
                                  "': want NET.STA or NET.STA.LOC");
    StreamSpec spec{parts[0], parts[1], std::nullopt};
    if (parts.size() == 3)
      spec.loc = parts[2] == "--" ? std::string() : parts[2];
    for (const auto &o : out)
      if (o.sta == spec.sta)
        throw std::invalid_argument(
            std::format("{} and {}: a station code may be given only once",
                        o.id(), spec.id()));
    out.push_back(std::move(spec));
  }
  return out;
}

Endpoint parse_endpoint(std::string_view s) {
  const auto bad = [&] {
    return std::invalid_argument("bad server address '" + std::string(s) +
                                 "': want HOST:PORT");
  };
  Endpoint e;
  std::string_view port;
  if (s.starts_with('[')) { // [IPv6]:port
    const auto close = s.find(']');
    if (close == std::string_view::npos)
      throw bad();
    e.host = s.substr(1, close - 1);
    const auto rest = s.substr(close + 1);
    if (!rest.empty() && !rest.starts_with(':'))
      throw bad();
    port = rest.empty() ? "18000" : rest.substr(1);
  } else {
    const auto colon = s.rfind(':');
    e.host = s.substr(0, colon);
    port = colon == std::string_view::npos ? "18000" : s.substr(colon + 1);
  }
  e.port = port;
  if (e.host.empty() || e.port.empty() ||
      !std::ranges::all_of(e.port, [](char c) { return c >= '0' && c <= '9'; }))
    throw bad();
  return e;
}

void FrameReader::feed(std::span<const std::byte> bytes) {
  // Drop what has been consumed once it is the larger part of the buffer.
  if (head_ > 0 && head_ * 2 >= buf_.size()) {
    buf_.erase(buf_.begin(), buf_.begin() + static_cast<std::ptrdiff_t>(head_));
    head_ = 0;
  }
  buf_.insert(buf_.end(), bytes.begin(), bytes.end());
}

void FrameReader::clear() noexcept {
  buf_.clear();
  head_ = 0;
}

std::optional<Frame> FrameReader::next() {
  const std::size_t avail = buf_.size() - head_;
  if (avail == 0)
    return std::nullopt;
  const std::byte *p = buf_.data() + head_;
  auto at = [&](std::size_t i) { return static_cast<char>(p[i]); };
  // True if the buffered bytes agree with `s` as far as they go.
  auto agrees = [&](std::string_view s) {
    for (std::size_t i = 0; i < std::min(avail, s.size()); ++i)
      if (at(i) != s[i])
        return false;
    return true;
  };
  auto bad = [&] {
    std::string shown;
    for (std::size_t i = 0; i < std::min<std::size_t>(avail, 8); ++i)
      shown += at(i) >= ' ' && at(i) <= '~'
                   ? std::string(1, at(i))
                   : std::format("\\x{:02x}", static_cast<unsigned char>(at(i)));
    return ProtocolError("unexpected bytes in the packet stream: \"" + shown +
                         "\"");
  };

  Frame f;
  if (agrees("SL")) {
    if (avail < kHeaderLength)
      return std::nullopt;
    if (agrees("SLINFO ")) {
      if (at(7) != '*' && at(7) != ' ')
        throw bad();
      f.kind = Frame::Kind::Info;
      f.more = at(7) == '*';
    } else {
      std::uint32_t seq = 0;
      for (std::size_t i = 2; i < kHeaderLength; ++i) {
        const int d = hex_digit(p[i]);
        if (d < 0)
          throw bad();
        seq = seq << 4 | static_cast<std::uint32_t>(d);
      }
      f.seq = seq;
    }
    if (avail < kPacketLength)
      return std::nullopt;
    std::memcpy(f.record.data(), p + kHeaderLength, kRecordLength);
    head_ += kPacketLength;
    return f;
  }
  if (agrees("ERROR")) {
    // "ERROR\r\n", possibly with a description before the line end.
    for (std::size_t i = 0; i < std::min<std::size_t>(avail, 256); ++i)
      if (at(i) == '\n') {
        f.kind = Frame::Kind::Error;
        head_ += i + 1;
        return f;
      }
    if (avail >= 256)
      throw bad();
    return std::nullopt;
  }
  if (agrees("END")) { // the end of a dial-up transfer, never sent after DATA
    if (avail < 3)
      return std::nullopt;
    f.kind = Frame::Kind::End;
    head_ += 3;
    return f;
  }
  throw bad();
}

std::optional<Reply> parse_reply(std::string_view line) {
  while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
    line.remove_suffix(1);
  if (line == "OK")
    return Reply{true, {}};
  if (line.starts_with("ERROR")) {
    line.remove_prefix(5);
    while (!line.empty() && line.front() == ' ')
      line.remove_prefix(1);
    return Reply{false, std::string(line)};
  }
  return std::nullopt;
}

std::string station_command(const StreamSpec &s) {
  return std::format("STATION {} {}", s.sta, s.net);
}

std::string select_command(const StreamSpec &s) {
  return "SELECT " + s.selector();
}

std::string data_command(std::optional<std::uint32_t> last) {
  return last ? std::format("DATA {:06X}", seq_next(*last)) : "DATA";
}

RecordId record_id(Record rec) {
  RecordId id{trimmed(rec, 18, 2), trimmed(rec, 8, 5), trimmed(rec, 13, 2),
              trimmed(rec, 15, 3), std::nullopt};
  const auto year = mseed::read_be<std::uint16_t>(rec, 20);
  const auto doy = mseed::read_be<std::uint16_t>(rec, 22);
  if (auto t = mseed::btime_to_ns(year, doy, static_cast<std::uint8_t>(rec[24]),
                                  static_cast<std::uint8_t>(rec[25]),
                                  static_cast<std::uint8_t>(rec[26]),
                                  mseed::read_be<std::uint16_t>(rec, 28), 0))
    id.start_ns = *t;
  return id;
}

Client::Client(ClientConfig cfg, RecordFn on_record, IdleFn on_idle)
    : cfg_(std::move(cfg)), on_record_(std::move(on_record)),
      on_idle_(std::move(on_idle)), stations_(cfg_.streams.size()) {
  for (std::size_t i = 0; i < cfg_.streams.size(); ++i)
    index_[cfg_.streams[i].net + "." + cfg_.streams[i].sta] = i;
}

bool Client::deliver(const Frame &f) {
  const RecordId id = record_id(f.record);
  const auto it = index_.find(id.net + "." + id.sta);
  if (it == index_.end()) {
    ++stats_.unknown;
    return false;
  }
  Station &s = stations_[it->second];
  // A packet at or before the last one delivered is one the server has sent
  // again, unless its data is newer than anything of its channel delivered:
  // then the server's numbering has started again (a restarted server), and
  // it is followed from here.
  if (s.last_seq && !seq_after(f.seq, *s.last_seq)) {
    const auto latest = s.latest_start.find(id.cha);
    const bool newer = id.start_ns && (latest == s.latest_start.end() ||
                                       *id.start_ns > latest->second);
    if (!newer) {
      ++stats_.duplicates;
      return false;
    }
    ++stats_.seq_restarts;
    if (!s.warned_restart) {
      s.warned_restart = true;
      Log::get().line("warn", "33",
                      "{}: sequence number {:06X} after {:06X} with new data; "
                      "taking it as a restarted numbering",
                      cfg_.streams[it->second].id(), f.seq, *s.last_seq);
    }
  }
  s.last_seq = f.seq;
  if (id.start_ns) {
    auto &latest = s.latest_start[id.cha];
    latest = std::max(latest, *id.start_ns);
  }
  ++stats_.packets;
  on_record_(it->second, f.record);
  return true;
}

#ifdef _WIN32

bool Client::available() noexcept { return false; }

void Client::run(const std::function<bool()> &) {
  throw std::runtime_error(
      "live SeedLink input is not available in the Windows build");
}

void Client::session(const std::function<bool()> &, bool &) {}

#else

bool Client::available() noexcept { return true; }

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point t) {
  return std::chrono::duration<double>(Clock::now() - t).count();
}

// A failure of one connection: logged, and the client reconnects.
class Dropped : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

// A failure that reconnecting cannot mend: the stream list is wrong.
class Fatal : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

std::string errno_text(int e) { return std::strerror(e); }

#ifdef MSG_NOSIGNAL
constexpr int kSendFlags = MSG_NOSIGNAL;
#else
constexpr int kSendFlags = 0;
#endif

// A connected non-blocking socket. Every wait is a poll of at most 20 ms,
// after which `idle` runs and `stop` is checked.
class Socket {
public:
  Socket(int &fd, const std::function<bool()> &stop,
         const std::function<void()> &idle)
      : fd_(fd), stop_(stop), idle_(idle) {}

  // Tries every address of host:port until one connects within `timeout_s`.
  // False if `stop` was set meanwhile.
  bool connect(const Endpoint &ep, double timeout_s) {
    addrinfo hints{};
    hints.ai_family = AF_UNSPEC;
    hints.ai_socktype = SOCK_STREAM;
    addrinfo *res = nullptr;
    if (const int rc = getaddrinfo(ep.host.c_str(), ep.port.c_str(), &hints, &res);
        rc != 0)
      throw Dropped(std::format("cannot resolve {}: {}", ep.host,
                                gai_strerror(rc)));
    std::string why = "no address";
    for (addrinfo *a = res; a != nullptr; a = a->ai_next) {
      const int fd = ::socket(a->ai_family, a->ai_socktype, a->ai_protocol);
      if (fd < 0) {
        why = errno_text(errno);
        continue;
      }
      ::fcntl(fd, F_SETFL, ::fcntl(fd, F_GETFL, 0) | O_NONBLOCK);
      const int one = 1;
      ::setsockopt(fd, SOL_SOCKET, SO_KEEPALIVE, &one, sizeof one);
      int rc = ::connect(fd, a->ai_addr, a->ai_addrlen);
      if (rc != 0 && errno == EINPROGRESS) {
        const auto t0 = Clock::now();
        rc = -1;
        errno = ETIMEDOUT;
        while (seconds_since(t0) < timeout_s) {
          if (stop_()) {
            ::close(fd);
            ::freeaddrinfo(res);
            return false;
          }
          pollfd p{fd, POLLOUT, 0};
          if (::poll(&p, 1, 20) > 0) {
            int err = 0;
            socklen_t len = sizeof err;
            ::getsockopt(fd, SOL_SOCKET, SO_ERROR, &err, &len);
            rc = err == 0 ? 0 : -1;
            errno = err;
            break;
          }
          idle_();
        }
      }
      if (rc == 0) {
        fd_ = fd;
        ::freeaddrinfo(res);
        return true;
      }
      why = errno_text(errno);
      ::close(fd);
    }
    ::freeaddrinfo(res);
    throw Dropped(std::format("cannot connect to {}:{}: {}", ep.host, ep.port,
                              why));
  }

  void send_line(std::string_view cmd, double timeout_s) {
    const std::string line = std::string(cmd) + "\r\n";
    std::size_t sent = 0;
    const auto t0 = Clock::now();
    while (sent < line.size()) {
      const ssize_t n =
          ::send(fd_, line.data() + sent, line.size() - sent, kSendFlags);
      if (n > 0) {
        sent += static_cast<std::size_t>(n);
        continue;
      }
      if (n < 0 && errno != EAGAIN && errno != EWOULDBLOCK && errno != EINTR)
        throw Dropped("send: " + errno_text(errno));
      if (seconds_since(t0) > timeout_s)
        throw Dropped("send: timed out");
      pollfd p{fd_, POLLOUT, 0};
      ::poll(&p, 1, 20);
      idle_();
    }
  }

  // Receives what has arrived, waiting up to 20 ms for something. Throws
  // Dropped when the server has closed the connection.
  std::span<const std::byte> receive() {
    pollfd p{fd_, POLLIN, 0};
    if (::poll(&p, 1, 20) <= 0)
      return {};
    const ssize_t n = ::recv(fd_, buf_.data(), buf_.size(), 0);
    if (n == 0)
      throw Dropped("connection closed by the server");
    if (n < 0) {
      if (errno == EAGAIN || errno == EWOULDBLOCK || errno == EINTR)
        return {};
      throw Dropped("receive: " + errno_text(errno));
    }
    return std::span<const std::byte>(buf_.data(), static_cast<std::size_t>(n));
  }

  // One line of the handshake, without its line end. Bytes after it stay in
  // `pending` for the next call or for the packet stream.
  std::string read_line(double timeout_s, std::vector<std::byte> &pending) {
    const auto t0 = Clock::now();
    for (;;) {
      const auto nl = std::ranges::find(pending, std::byte{'\n'});
      if (nl != pending.end()) {
        std::string line(reinterpret_cast<const char *>(pending.data()),
                         static_cast<std::size_t>(nl - pending.begin()));
        pending.erase(pending.begin(), nl + 1);
        if (!line.empty() && line.back() == '\r')
          line.pop_back();
        return line;
      }
      if (pending.size() > 4096)
        throw Dropped("handshake: reply too long, not a SeedLink server?");
      if (stop_())
        throw Dropped("stopped");
      if (seconds_since(t0) > timeout_s)
        throw Dropped(std::format("handshake: no reply in {:g} s", timeout_s));
      const auto got = receive();
      pending.insert(pending.end(), got.begin(), got.end());
      idle_();
    }
  }

private:
  int &fd_;
  const std::function<bool()> &stop_;
  const std::function<void()> &idle_;
  // On the heap: musl gives a thread 128 KiB of stack.
  std::vector<std::byte> buf_ = std::vector<std::byte>(1 << 16);
};

} // namespace

void Client::session(const std::function<bool()> &stop, bool &got_data) {
  Socket sock(fd_, stop, on_idle_);
  if (!sock.connect(cfg_.server, std::min(cfg_.timeout_s, 10.0)))
    return;
  ++stats_.connections;
  const std::string where =
      std::format("{}:{}", cfg_.server.host, cfg_.server.port);
  std::vector<std::byte> pending;

  sock.send_line("HELLO", cfg_.timeout_s);
  const std::string version = sock.read_line(cfg_.timeout_s, pending);
  const std::string organisation = sock.read_line(cfg_.timeout_s, pending);
  if (!handshaken_)
    Log::get().line("seedlink", "36", "{}: {} ({})", where, version,
                    organisation);

  // A refusal in the first session means a wrong stream list and ends the
  // run; later, the station is left out of this session only.
  std::size_t accepted = 0;
  auto ask = [&](const std::string &cmd, const StreamSpec &spec) {
    sock.send_line(cmd, cfg_.timeout_s);
    const std::string line = sock.read_line(cfg_.timeout_s, pending);
    const auto reply = parse_reply(line);
    if (!reply)
      throw Dropped(std::format("handshake: unexpected reply \"{}\" to {}",
                                line, cmd));
    if (reply->ok)
      return true;
    const std::string what =
        std::format("{} refused {}{}", where, cmd,
                    reply->text.empty() ? "" : " (" + reply->text + ")");
    if (!handshaken_)
      throw Fatal(std::format("{}: {}", spec.id(), what));
    Log::get().line("warn", "33", "{}: {}; left out until the next connection",
                    spec.id(), what);
    return false;
  };
  for (std::size_t i = 0; i < cfg_.streams.size(); ++i) {
    const auto &spec = cfg_.streams[i];
    if (ask(station_command(spec), spec) && ask(select_command(spec), spec) &&
        ask(data_command(stations_[i].last_seq), spec))
      ++accepted;
  }
  if (accepted == 0)
    throw Dropped("the server accepted none of the stations");
  sock.send_line("END", cfg_.timeout_s);
  ++stats_.sessions;
  Log::get().line("seedlink", "36", "{}: {} streaming {} of {} stations", where,
                  handshaken_ ? "reconnected," : "connected,", accepted,
                  cfg_.streams.size());
  handshaken_ = true;

  FrameReader reader;
  reader.feed(pending);
  auto last_rx = Clock::now();
  bool keepalive_pending = false;
  for (;;) {
    while (auto f = reader.next()) {
      switch (f->kind) {
      case Frame::Kind::Data:
        got_data = deliver(*f) || got_data;
        break;
      case Frame::Kind::Info:
        ++stats_.info_packets;
        keepalive_pending = keepalive_pending && f->more;
        break;
      case Frame::Kind::Error:
        // The only command sent in this phase is INFO ID.
        if (stats_.server_errors++ == 0)
          Log::get().line("warn", "33",
                          "{}: the server answered INFO with ERROR (INFO "
                          "disabled?); keepalives go on regardless",
                          where);
        keepalive_pending = false;
        break;
      case Frame::Kind::End:
        throw Dropped("the server ended the transfer");
      }
    }
    if (stop()) {
      // Best effort: the server closes its side on BYE.
      const std::string bye = "BYE\r\n";
      [[maybe_unused]] const auto n =
          ::send(fd_, bye.data(), bye.size(), kSendFlags);
      return;
    }
    const auto got = sock.receive();
    if (!got.empty()) {
      last_rx = Clock::now();
      reader.feed(got);
    }
    on_idle_();
    const double quiet = seconds_since(last_rx);
    if (quiet > cfg_.timeout_s)
      throw Dropped(std::format("nothing received for {:g} s", cfg_.timeout_s));
    if (cfg_.keepalive_s > 0 && quiet > cfg_.keepalive_s && !keepalive_pending) {
      sock.send_line("INFO ID", cfg_.timeout_s);
      ++stats_.keepalives;
      keepalive_pending = true;
    }
  }
}

void Client::run(const std::function<bool()> &stop) {
  // A write to a connection the server has reset must not kill the process.
  std::signal(SIGPIPE, SIG_IGN);
  const std::string where =
      std::format("{}:{}", cfg_.server.host, cfg_.server.port);
  double delay = cfg_.first_delay_s;
  while (!stop()) {
    bool got_data = false;
    try {
      session(stop, got_data);
    } catch (const Fatal &e) {
      if (fd_ >= 0)
        ::close(fd_);
      fd_ = -1;
      throw std::runtime_error(e.what());
    } catch (const std::exception &e) { // Dropped, ProtocolError
      if (!stop())
        Log::get().line("warn", "33", "{}: {}", where, e.what());
    }
    if (fd_ >= 0)
      ::close(fd_);
    fd_ = -1;
    if (stop())
      break;
    if (got_data)
      delay = cfg_.first_delay_s;
    Log::get().line("seedlink", "36", "{}: reconnecting in {:g} s", where,
                    delay);
    for (const auto t0 = Clock::now(); seconds_since(t0) < delay && !stop();) {
      std::this_thread::sleep_for(std::chrono::milliseconds(20));
      on_idle_();
    }
    delay = std::min(delay * 2, cfg_.max_delay_s);
  }
}

#endif

} // namespace ayzek::seedlink
