#pragma once

// Types shared by the pipeline threads: time helpers, the replay clock, the
// messages sent from station processors to the network stage, the message
// queue, and the logger.

#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <format>
#include <mutex>
#include <optional>
#include <print>
#include <string>
#include <variant>

namespace ayzek::pipeline {

inline constexpr double kFs = 100.0;
inline constexpr std::uint64_t kUnset = ~std::uint64_t{0};

inline double pos_to_epoch(std::uint64_t pos) { return static_cast<double>(pos) / kFs; }
inline std::uint64_t epoch_to_pos(double t) { return static_cast<std::uint64_t>(t * kFs + 0.5); }

// "18:21:02.50" in UTC.
inline std::string hms(double epoch) {
    using namespace std::chrono;
    const auto t = sys_time<milliseconds>(milliseconds(static_cast<std::int64_t>(epoch * 1000.0 + 0.5)));
    return std::format("{:%H:%M:%S}", t).substr(0, 11);
}

inline std::string ymd_hms(double epoch) {
    using namespace std::chrono;
    const auto t = sys_time<milliseconds>(milliseconds(static_cast<std::int64_t>(epoch * 1000.0 + 0.5)));
    return std::format("{:%Y-%m-%d %H:%M:%S}", t).substr(0, 22);
}

// Parses a UTC date-time into epoch seconds, taking the first six integers in
// `s` in the given order: "ymdhms" for 2025-11-10T18:20:51, "dmyhms" for the
// AFAD format 10/11/2025 18:20:51. Written by hand because libc++ (used by the
// Raspberry Pi build) does not provide std::chrono::parse.
inline std::optional<double> parse_utc(std::string_view s, std::string_view order) {
    int v[6] = {};
    std::size_t n = 0, i = 0;
    while (n < 6 && i < s.size()) {
        if (s[i] < '0' || s[i] > '9') { ++i; continue; }
        int x = 0;
        while (i < s.size() && s[i] >= '0' && s[i] <= '9') x = x * 10 + (s[i++] - '0');
        v[n++] = x;
    }
    if (n < 6) return std::nullopt;
    int y = 0, mo = 0, d = 0;
    for (std::size_t k = 0; k < 3; ++k) {
        if (order[k] == 'y') y = v[k];
        else if (order[k] == 'm') mo = v[k];
        else d = v[k];
    }
    const std::chrono::year_month_day ymd{std::chrono::year{y}, std::chrono::month{static_cast<unsigned>(mo)},
                                          std::chrono::day{static_cast<unsigned>(d)}};
    if (!ymd.ok()) return std::nullopt;
    const auto days = std::chrono::sys_days{ymd}.time_since_epoch().count();
    return static_cast<double>(days) * 86400.0 + v[3] * 3600.0 + v[4] * 60.0 + v[5];
}

// Stream time for replay: `start` plus `speed` times the elapsed wall time.
// With speed 0, now() returns a time beyond any data, so records are released
// without delay.
class StreamClock {
public:
    StreamClock(double start_epoch, double speed)
        : start_(start_epoch), speed_(speed), wall0_(std::chrono::steady_clock::now()) {}
    [[nodiscard]] double now() const {
        if (speed_ <= 0.0) return 1e18;
        const std::chrono::duration<double> dt = std::chrono::steady_clock::now() - wall0_;
        return start_ + dt.count() * speed_;
    }
    [[nodiscard]] double speed() const noexcept { return speed_; }
    [[nodiscard]] double wall_seconds() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0_).count();
    }
private:
    double start_, speed_;
    std::chrono::steady_clock::time_point wall0_;
};

// --- messages ---------------------------------------------------------------------

struct Detection {
    std::string station;
    double window_start;      // epoch of the first window above threshold
    double declared_at;       // stream time at which the decision is available (window end)
    float probability;
    double compute_ms;        // conditioning + ensemble for that window
};

struct Pick {
    std::string station;
    double trigger_window;    // window_start of the Detection this pick follows
    double p_time, s_time;    // epoch
    double p_prob, s_prob;
    double declared_at;       // stream time at the end of the 60 s picker window
    double compute_ms;
};

struct MagnitudeEstimate {
    std::string station;
    double trigger_window;    // window_start of the Detection this follows
    bool at_pick;             // false: window placed from the trigger time; true: from the picked P
    double window_start;
    float magnitude;
    std::size_t noise_windows; // noise windows in the station baseline; 0 = per-window normalisation
    double declared_at;
    double compute_ms;
};

struct StationDone {
    std::string station;
};

// Sent by a station processor to state that none of its later messages has
// declared_at earlier than `until`. The network stage processes messages in
// declared_at order up to the minimum `until` over all stations, so the result
// does not depend on thread scheduling or replay speed.
struct Progress {
    std::string station;
    double until;
};

using Message = std::variant<Detection, Pick, MagnitudeEstimate, StationDone, Progress>;

inline double declared_at(const Message& m) {
    if (auto* d = std::get_if<Detection>(&m)) return d->declared_at;
    if (auto* p = std::get_if<Pick>(&m)) return p->declared_at;
    if (auto* g = std::get_if<MagnitudeEstimate>(&m)) return g->declared_at;
    return 0;
}

// Mutex-protected message queue. Message rates are low (a few per event); the
// sample data itself goes through the lock-free rings.
class Bus {
public:
    void send(Message m) {
        {
            std::lock_guard lk(mu_);
            q_.push_back(std::move(m));
        }
        cv_.notify_one();
    }
    std::optional<Message> receive(std::chrono::milliseconds timeout) {
        std::unique_lock lk(mu_);
        if (!cv_.wait_for(lk, timeout, [&] { return !q_.empty(); })) return std::nullopt;
        Message m = std::move(q_.front());
        q_.pop_front();
        return m;
    }
private:
    std::mutex mu_;
    std::condition_variable cv_;
    std::deque<Message> q_;
};

// Serialised line output, so that lines from different threads do not mix.
class Log {
public:
    static Log& get() {
        static Log l;
        return l;
    }
    bool color = true;
    bool quiet = false;

    template <typename... Args>
    void line(std::string_view tag, std::string_view tint, std::format_string<Args...> fmt, Args&&... args) {
        std::string body = std::format(fmt, std::forward<Args>(args)...);
        std::lock_guard lk(mu_);
        if (color) std::println("\x1b[{}m{:<8}\x1b[0m {}", tint, tag, body);
        else std::println("{:<8} {}", tag, body);
        std::fflush(stdout);
    }
private:
    std::mutex mu_;
};

}  // namespace ayzek::pipeline
