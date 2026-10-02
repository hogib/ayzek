#pragma once

// Magnitude from peak P-wave displacement (Pd), the estimator calibrated in
// onset (onset/docs/MAGNITUDE.md; docs/impl/17-pd-magnitude.md here).
//
// A station's raw vertical counts pass through a causal chain, a two-pole
// high-pass at 0.075 Hz, integration and the same high-pass again
// (Displacement), and Pd(tau) is the largest |displacement| in the tau
// seconds after the P that the transformer dates, divided by the channel's
// sensitivity. An event's magnitude is the posterior mean, under a
// Gutenberg-Richter prior, of the censored likelihood of its stations' values
// under the relation
//
//   log10 Pd = alpha + beta_m M + beta_above max(0, M - m_knot)
//            + gamma_r lr + gamma_k1 max(0, lr - log10 k1) + gamma_k2 max(0, lr - log10 k2)
//            + s_station,     lr = log10 sqrt(D^2 + depth^2),
//
// with residual sd sigma, where a value below min_snr times the pre-P noise
// level is censored. The coefficients, the station terms and the stations'
// gains are exported by tools/export_pd.py; tests/test_pd.cpp holds the chain
// and the estimator to onset's pd.py and pd_fit.py.

#include <array>
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace ayzek {

struct PdRelation {
  double tau_s = 0, alpha = 0, beta_m = 0, beta_above = 0, m_knot = 0;
  double gamma_r = 0, gamma_k1 = 0, gamma_k2 = 0, k1_km = 0, k2_km = 0;
  double sigma = 1, b_value = 0, min_snr = 3, depth_km = 10, flag_log10 = 1;
  double dead_noise_factor = 10;
  double coda_noise_factor = 10;
  // log10 Pd predicted at magnitude m and epicentral distance d_km, without
  // the station term.
  [[nodiscard]] double mean(double m, double d_km) const;
};

struct PdModel {
  std::vector<PdRelation> windows; // by increasing tau_s
  // Second-order section of the high-pass: b0 b1 b2 a1 a2 (a0 = 1).
  std::array<double, 5> highpass{};
  struct Epoch {
    double start, end; // POSIX seconds; end infinite while open
    double sensitivity; // counts per m/s
  };
  std::map<std::string, std::vector<Epoch>> gains;
  std::map<std::string, std::map<double, double>> terms; // station -> tau -> term
  std::map<std::string, double> noise; // station -> median pre-P noise level, m

  // models/pd_relation.csv, pd_gain.csv and pd_terms.csv in `dir`; nullopt
  // when pd_relation.csv is absent.
  static std::optional<PdModel> load(const std::string &dir);
  [[nodiscard]] std::optional<double> sensitivity(const std::string &station,
                                                  double t) const;
  [[nodiscard]] double term(const std::string &station, double tau) const;
  // A station whose term exceeds flag_log10 in magnitude: a probable
  // response error, left out of the estimates (as in onset).
  [[nodiscard]] bool flagged(const std::string &station, double tau) const;
  // A value whose pre-P noise lies outside the station's usual band, as in
  // onset: more than dead_noise_factor below its median (the channel is not
  // recording, and its spurious upper bound would pull the estimate down) or
  // more than coda_noise_factor above it (measured inside another event's
  // coda, whose energy inflates the peak after P).
  [[nodiscard]] bool outside_noise_band(const std::string &station,
                                        double pd_noise) const;
  [[nodiscard]] const PdRelation *window(double tau) const;
};

// The causal displacement chain on raw counts. Each contiguous run starts
// from rest at its first sample, which is subtracted, as onset's pd.py starts
// its window; a gap ends the run. The output is in counts (divide by the
// sensitivity for metres).
class Displacement {
public:
  explicit Displacement(const std::array<double, 5> &highpass,
                        double fs = 100.0)
      : hp_(highpass), fs_(fs) {}
  // nullopt for a missing sample.
  std::optional<double> push(std::optional<double> counts);
  [[nodiscard]] std::size_t run_length() const { return run_; }

private:
  double section(double x, std::array<double, 2> &z) const;
  std::array<double, 5> hp_;
  double fs_;
  std::array<double, 2> z1_{}, z2_{};
  double offset_ = 0, integral_ = 0;
  std::size_t run_ = 0;
};

// One station's values for an event.
struct PdObservation {
  std::string station;
  double tau_s;       // window
  double pd;          // metres
  double pd_noise;    // metres
  double distance_km; // epicentral
};

struct PdMagnitude {
  double mean, sd;
  std::size_t stations; // used
  std::size_t above;    // of them, above the threshold
};

// The event's posterior over M on [1, 8.5] in steps of 0.01, as onset's
// pd_fit.estimate_censored. nullopt when no usable station records the event
// above its threshold.
std::optional<PdMagnitude> pd_magnitude(const PdModel &model,
                                        const std::vector<PdObservation> &obs);

} // namespace ayzek
