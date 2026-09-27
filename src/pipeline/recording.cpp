#include "recording.hpp"

#include <sstream>
#include <stdexcept>

namespace ayzek::pipeline {

Recorder::Recorder(const std::string &path,
                   const std::map<std::string, StationInfo> &stations,
                   double t_first, double t_last)
    : out_(path) {
  if (!out_)
    throw std::runtime_error("cannot write " + path);
  out_ << "# ayzek recording v1\n";
  for (const auto &[code, s] : stations)
    out_ << std::format("S\t{}\t{}\t{}\n", code, s.lat, s.lon);
  out_ << std::format("T\t{}\t{}\n", t_first, t_last);
}

void Recorder::write(const Message &m) {
  if (auto *d = std::get_if<Detection>(&m)) {
    out_ << std::format("D\t{}\t{}\t{}\t{}\t{}\t{}\n", d->station,
                        d->window_start, d->declared_at, d->probability,
                        d->compute_ms, d->restart ? 1 : 0);
  } else if (auto *p = std::get_if<Pick>(&m)) {
    out_ << std::format("P\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", p->station,
                        p->trigger_window, p->p_time, p->s_time, p->p_prob,
                        p->s_prob, p->declared_at, p->compute_ms);
    if (p->event_id)
      out_ << std::format("\t{}\t{}", p->event_id, p->event_alarm);
    out_ << '\n';
  } else if (auto *g = std::get_if<MagnitudeEstimate>(&m)) {
    out_ << std::format("M\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", g->station,
                        g->trigger_window, g->at_pick ? 1 : 0, g->window_start,
                        g->magnitude, g->noise_windows, g->declared_at,
                        g->compute_ms);
  } else if (auto *o = std::get_if<StationGeometry>(&m)) {
    out_ << std::format("G\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        o->station, o->trigger_window, o->p_time, o->since_p,
                        o->log_dist, o->log_dist_sd, o->baz, o->kappa,
                        o->declared_at);
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
      // The restart flag is absent in recordings made before it existed.
      if (!f.fail()) {
        int restart = 0;
        if (f >> restart)
          d.restart = restart != 0;
        else if (!f.eof())
          throw std::runtime_error(
              std::format("{}:{}: malformed record", path, lineno));
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
      MagnitudeEstimate g;
      int at_pick = 0;
      f >> g.station >> g.trigger_window >> at_pick >> g.window_start >>
          g.magnitude >> g.noise_windows >> g.declared_at >> g.compute_ms;
      g.at_pick = at_pick != 0;
      r.messages.emplace_back(std::move(g));
    } else if (tag == "G") {
      StationGeometry o;
      f >> o.station >> o.trigger_window >> o.p_time >> o.since_p >>
          o.log_dist >> o.log_dist_sd >> o.baz >> o.kappa >> o.declared_at;
      r.messages.emplace_back(std::move(o));
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
  if (auto *g = std::get_if<MagnitudeEstimate>(&m))
    return g->station;
  if (auto *o = std::get_if<StationGeometry>(&m))
    return o->station;
  throw std::invalid_argument("message has no station");
}

} // namespace ayzek::pipeline
