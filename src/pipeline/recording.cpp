#include "recording.hpp"

#include <sstream>
#include <stdexcept>
#include <vector>

namespace ayzek::pipeline {

Recorder::Recorder(const std::string &path,
                   const std::map<std::string, StationInfo> &stations,
                   double t_first, double t_last,
                   std::size_t min_stations, double from)
    : out_(path) {
  if (!out_)
    throw std::runtime_error("cannot write " + path);
  out_ << "# ayzek recording v1\n";
  out_ << std::format("# min_stations {}\n", min_stations);
  out_ << std::format("# from {}\n", from);
  for (const auto &[code, s] : stations)
    out_ << std::format("S\t{}\t{}\t{}\n", code, s.lat, s.lon);
  out_ << std::format("T\t{}\t{}\n", t_first, t_last);
}

void Recorder::write(const Message &m) {
  if (auto *d = std::get_if<Detection>(&m)) {
    out_ << std::format("D\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", d->station,
                        d->window_start, d->declared_at, d->probability,
                        d->compute_ms, d->restart ? 1 : 0, p_of(*d));
  } else if (auto *p = std::get_if<Pick>(&m)) {
    out_ << std::format("P\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", p->station,
                        p->trigger_window, p->p_time, p->s_time, p->p_prob,
                        p->s_prob, p->declared_at, p->compute_ms);
    if (p->event_id)
      out_ << std::format("\t{}\t{}", p->event_id, p->event_alarm);
    out_ << '\n';
  } else if (auto *o = std::get_if<StationGeometry>(&m)) {
    out_ << std::format("G\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", o->station,
                        o->trigger_window, o->p_time, o->since_p, o->log_dist,
                        o->log_dist_sd, o->declared_at);
  } else if (auto *q = std::get_if<PdEstimate>(&m)) {
    out_ << std::format("Q\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", q->station,
                        q->trigger_window, q->p_time, q->tau, q->pd, q->pd_noise,
                        q->declared_at);
  }
}

Recording Recording::load(const std::string &path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("cannot open " + path);
  Recording r;
  std::string line;
  std::size_t lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    if (line.starts_with("# min_stations ")) {
      r.min_stations = std::stoul(line.substr(15));
      continue;
    }
    if (line.starts_with("# from ")) {
      r.from = std::stod(line.substr(7));
      continue;
    }
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream f(line);
    std::string tag;
    f >> tag;
    if (tag == "S") {
      std::string code;
      StationInfo s;
      f >> code >> s.lat >> s.lon;
      r.stations[code] = s;
    } else if (tag == "T") {
      f >> r.t_first >> r.t_last;
    } else if (tag == "D") {
      Detection d;
      f >> d.station >> d.window_start >> d.declared_at >> d.probability >>
          d.compute_ms;
      // The restart flag and the P time are absent in recordings made
      // before they existed; P is then the window start + 3.5 s.
      d.p_time = d.window_start + 3.5;
      if (!f.fail()) {
        int restart = 0;
        if (f >> restart) {
          d.restart = restart != 0;
          double p = 0;
          if (f >> p)
            d.p_time = p;
          else if (!f.eof())
            throw std::runtime_error(
                std::format("{}:{}: malformed record", path, lineno));
        } else if (!f.eof()) {
          throw std::runtime_error(
              std::format("{}:{}: malformed record", path, lineno));
        }
        f.clear();
      }
      r.messages.emplace_back(std::move(d));
    } else if (tag == "P") {
      Pick p;
      f >> p.station >> p.trigger_window >> p.p_time >> p.s_time >> p.p_prob >>
          p.s_prob >> p.declared_at >> p.compute_ms;
      if (!f.fail()) {
        if (f >> p.event_id >> p.event_alarm) {
        } else if (!f.eof()) {
          throw std::runtime_error(
              std::format("{}:{}: malformed record", path, lineno));
        }
        f.clear();
      }
      r.messages.emplace_back(std::move(p));
    } else if (tag == "M") {
      continue; // the retired magnitude regressor (recording.hpp)
    } else if (tag == "G") {
      StationGeometry o;
      f >> o.station >> o.trigger_window >> o.p_time >> o.since_p >>
          o.log_dist >> o.log_dist_sd;
      // declared_at is last: the 7th field, or the 9th in recordings made
      // when the head also gave a back-azimuth and its concentration.
      std::vector<double> rest;
      for (double v; f >> v;)
        rest.push_back(v);
      if (!f.eof())
        throw std::runtime_error(
            std::format("{}:{}: malformed record", path, lineno));
      f.clear(); // the loop ends by failing at the end of the line
      if (rest.size() != 1 && rest.size() != 3)
        throw std::runtime_error(
            std::format("{}:{}: malformed record", path, lineno));
      o.declared_at = rest.back();
      r.messages.emplace_back(std::move(o));
    } else if (tag == "Q") {
      PdEstimate q;
      f >> q.station >> q.trigger_window >> q.p_time >> q.tau >> q.pd >>
          q.pd_noise >> q.declared_at;
      r.messages.emplace_back(std::move(q));
    } else {
      throw std::runtime_error(
          std::format("{}:{}: unknown record type '{}'", path, lineno, tag));
    }
    if (f.fail())
      throw std::runtime_error(
          std::format("{}:{}: malformed record", path, lineno));
  }
  return r;
}

const std::string &station_of(const Message &m) {
  if (auto *d = std::get_if<Detection>(&m))
    return d->station;
  if (auto *p = std::get_if<Pick>(&m))
    return p->station;
  if (auto *o = std::get_if<StationGeometry>(&m))
    return o->station;
  if (auto *q = std::get_if<PdEstimate>(&m))
    return q->station;
  if (auto *q = std::get_if<PdTrace>(&m))
    return q->station;
  throw std::invalid_argument("message has no station");
}

} // namespace ayzek::pipeline
