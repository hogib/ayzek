#include "pd_magnitude.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace ayzek {

namespace {

// Rows of a CSV file with a header, as {column: value} maps.
std::vector<std::map<std::string, std::string>> read_csv(const std::string &path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("cannot open " + path);
  auto split = [](const std::string &line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ','))
      out.push_back(cell);
    if (!line.empty() && line.back() == ',')
      out.emplace_back();
    return out;
  };
  std::string line;
  auto next = [&] {
    if (!std::getline(in, line))
      return false;
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    return true;
  };
  next();
  const auto header = split(line);
  std::vector<std::map<std::string, std::string>> rows;
  while (next()) {
    if (line.empty())
      continue;
    const auto cells = split(line);
    std::map<std::string, std::string> r;
    for (std::size_t i = 0; i < header.size() && i < cells.size(); ++i)
      r[header[i]] = cells[i];
    rows.push_back(std::move(r));
  }
  return rows;
}

double num(const std::map<std::string, std::string> &r, const std::string &k) {
  auto it = r.find(k);
  if (it == r.end() || it->second.empty())
    throw std::runtime_error("pd model: missing column " + k);
  return std::stod(it->second);
}

// log of the standard normal CDF, accurate in the far left tail.
double log_ndtr(double z) {
  if (z > -20.0)
    return std::log(0.5 * std::erfc(-z / std::numbers::sqrt2));
  // Asymptotic series: Phi(z) ~ phi(z)/(-z) (1 - 1/z^2 + 3/z^4).
  const double z2 = z * z;
  return -0.5 * z2 - std::log(-z) - 0.5 * std::log(2 * std::numbers::pi) +
         std::log1p(-1.0 / z2 + 3.0 / (z2 * z2));
}

} // namespace

double PdRelation::mean(double m, double d_km) const {
  const double lr = std::log10(std::hypot(d_km, depth_km));
  return alpha + beta_m * m + beta_above * std::max(0.0, m - m_knot) +
         gamma_r * lr + gamma_k1 * std::max(0.0, lr - std::log10(k1_km)) +
         gamma_k2 * std::max(0.0, lr - std::log10(k2_km));
}

std::optional<PdModel> PdModel::load(const std::string &dir) {
  const std::string rel = dir + "/pd_relation.csv";
  if (!std::ifstream(rel))
    return std::nullopt;
  PdModel m;
  for (const auto &r : read_csv(rel)) {
    m.windows.push_back({num(r, "tau_s"), num(r, "alpha"), num(r, "beta_m"),
                         num(r, "beta_above"), num(r, "m_knot"),
                         num(r, "gamma_r"), num(r, "gamma_k1"),
                         num(r, "gamma_k2"), num(r, "k1_km"), num(r, "k2_km"),
                         num(r, "sigma"), num(r, "b_value"), num(r, "min_snr"),
                         num(r, "depth_km"), num(r, "flag_log10"),
                         num(r, "dead_noise_factor")});
    m.highpass = {num(r, "hp_b0"), num(r, "hp_b1"), num(r, "hp_b2"),
                  num(r, "hp_a1"), num(r, "hp_a2")};
  }
  if (m.windows.empty())
    throw std::runtime_error(rel + ": no windows");
  std::ranges::sort(m.windows, {}, &PdRelation::tau_s);
  for (const auto &r : read_csv(dir + "/pd_gain.csv")) {
    const auto &end = r.at("end");
    m.gains[r.at("station")].push_back(
        {num(r, "start"),
         end.empty() ? std::numeric_limits<double>::infinity() : std::stod(end),
         num(r, "sensitivity")});
  }
  for (const auto &r : read_csv(dir + "/pd_terms.csv"))
    m.terms[r.at("station")][num(r, "tau_s")] = num(r, "term");
  for (const auto &r : read_csv(dir + "/pd_noise.csv"))
    m.noise[r.at("station")] = num(r, "noise_m");
  return m;
}

std::optional<double> PdModel::sensitivity(const std::string &station,
                                           double t) const {
  auto it = gains.find(station);
  if (it == gains.end())
    return std::nullopt;
  for (const auto &e : it->second)
    if (t >= e.start && t < e.end)
      return e.sensitivity;
  return std::nullopt;
}

double PdModel::term(const std::string &station, double tau) const {
  auto it = terms.find(station);
  if (it == terms.end())
    return 0.0;
  auto jt = it->second.find(tau);
  return jt == it->second.end() ? 0.0 : jt->second;
}

bool PdModel::flagged(const std::string &station, double tau) const {
  const auto *w = window(tau);
  return w && std::abs(term(station, tau)) > w->flag_log10;
}

bool PdModel::not_recording(const std::string &station, double pd_noise) const {
  auto it = noise.find(station);
  return it != noise.end() && !windows.empty() &&
         pd_noise * windows.front().dead_noise_factor < it->second;
}

const PdRelation *PdModel::window(double tau) const {
  for (const auto &w : windows)
    if (std::abs(w.tau_s - tau) < 1e-9)
      return &w;
  return nullptr;
}

double Displacement::section(double x, std::array<double, 2> &z) const {
  // Transposed direct form II, as scipy's sosfilt.
  const double y = hp_[0] * x + z[0];
  z[0] = hp_[1] * x - hp_[3] * y + z[1];
  z[1] = hp_[2] * x - hp_[4] * y;
  return y;
}

std::optional<double> Displacement::push(std::optional<double> counts) {
  if (!counts) {
    run_ = 0;
    return std::nullopt;
  }
  if (run_ == 0) {
    offset_ = *counts;
    z1_ = {};
    z2_ = {};
    integral_ = 0;
  }
  ++run_;
  // onset's pd.displacement: high-pass, cumulative sum / fs, high-pass.
  integral_ += section(*counts - offset_, z1_) / fs_;
  return section(integral_, z2_);
}

std::optional<PdMagnitude> pd_magnitude(const PdModel &model,
                                        const std::vector<PdObservation> &obs) {
  struct Row {
    const PdRelation *rel;
    double y, base;
    bool censored;
  };
  std::vector<Row> rows;
  std::size_t above = 0;
  for (const auto &o : obs) {
    const auto *rel = model.window(o.tau_s);
    if (!rel || !(o.pd_noise > 0) || model.flagged(o.station, o.tau_s) ||
        model.not_recording(o.station, o.pd_noise))
      continue;
    const bool up = o.pd > 0 && o.pd >= rel->min_snr * o.pd_noise;
    above += up;
    rows.push_back({rel, std::log10(up ? o.pd : rel->min_snr * o.pd_noise),
                    rel->mean(0.0, o.distance_km) + model.term(o.station, o.tau_s),
                    !up});
  }
  if (above == 0)
    return std::nullopt;
  const double b = rows.front().rel->b_value;
  constexpr int kN = 751; // M = 1.00 ... 8.50
  std::array<double, kN> lp{};
  double peak = -std::numeric_limits<double>::infinity();
  for (int k = 0; k < kN; ++k) {
    const double m = 1.0 + 0.01 * k;
    double s = -b * std::numbers::ln10 * m;
    for (const auto &r : rows) {
      const auto &w = *r.rel;
      const double mu =
          r.base + w.beta_m * m + w.beta_above * std::max(0.0, m - w.m_knot);
      const double z = (r.y - mu) / w.sigma;
      s += r.censored ? log_ndtr(z)
                      : -0.5 * z * z - std::log(w.sigma) -
                            0.5 * std::log(2 * std::numbers::pi);
    }
    lp[static_cast<std::size_t>(k)] = s;
    peak = std::max(peak, s);
  }
  double z = 0, mean = 0, sq = 0;
  for (int k = 0; k < kN; ++k) {
    const double p = std::exp(lp[static_cast<std::size_t>(k)] - peak);
    const double m = 1.0 + 0.01 * k;
    z += p;
    mean += p * m;
    sq += p * m * m;
  }
  mean /= z;
  return PdMagnitude{mean, std::sqrt(std::max(0.0, sq / z - mean * mean)),
                     rows.size(), above};
}

} // namespace ayzek
