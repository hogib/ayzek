// Compares the Pd magnitude chain with onset's pd.py and pd_fit.py, through
// the reference outputs of tools/export_pd.py: the displacement of a
// synthetic record sample by sample, Pd and the noise level from it, and the
// event estimator on three cases; then the network stage's estimate for a
// declared event from PdEstimate messages.

#include "fixture.hpp"
#include "pd_magnitude.hpp"
#include "pipeline/network.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <vector>

using namespace ayzek;

namespace {

std::vector<std::vector<std::string>> csv_rows(const std::string &path) {
  if (!std::filesystem::exists(path)) {
    std::println("skipped: {} not found (run tools/export_pd.py)", path);
    std::exit(77);
  }
  std::ifstream in(path);
  std::vector<std::vector<std::string>> out;
  std::string line;
  std::getline(in, line); // header
  while (std::getline(in, line)) {
    std::vector<std::string> cells;
    std::stringstream ss(line);
    std::string c;
    while (std::getline(ss, c, ','))
      cells.push_back(c);
    out.push_back(cells);
  }
  return out;
}

} // namespace

int main(int argc, char **argv) {
  const auto root = fixture::root(argc, argv);
  const auto model = PdModel::load(root + "/models");
  if (!model) {
    std::println("skipped: models/pd_relation.csv not found (run tools/export_pd.py)");
    return 77;
  }

  // --- the displacement chain, sample by sample ---------------------------
  const auto chain = csv_rows(root + "/data/fixtures/pd_chain.csv");
  const auto peaks = csv_rows(root + "/data/fixtures/pd_chain_peaks.csv").at(0);
  const double sens = std::stod(peaks[0]);
  const auto p = static_cast<std::size_t>(std::stod(peaks[1]));
  Displacement disp(model->highpass);
  std::vector<double> d;
  double worst = 0, scale = 0;
  for (const auto &r : chain) {
    d.push_back(*disp.push(std::stod(r[0])) / sens);
    const double ref = std::stod(r[1]);
    worst = std::max(worst, std::abs(d.back() - ref));
    scale = std::max(scale, std::abs(ref));
  }
  std::println("  displacement: {} samples, max abs err {:.2e} (signal max {:.2e})",
               d.size(), worst, scale);
  CHECK(worst <= 1e-9 * scale);

  // Pd(tau) and the noise level, as pd.peak_displacement.
  const double windows[] = {1, 2, 3, 4, 5, 7, 10};
  for (std::size_t k = 0; k < std::size(windows); ++k) {
    double pd = 0;
    for (std::size_t i = p; i < p + static_cast<std::size_t>(windows[k] * 100); ++i)
      pd = std::max(pd, std::abs(d[i]));
    CHECK(std::abs(pd / std::stod(peaks[2 + k]) - 1) < 1e-9);
  }
  double noise = 0;
  for (std::size_t i = p - 1100; i < p - 100; ++i)
    noise = std::max(noise, std::abs(d[i]));
  CHECK(std::abs(noise / std::stod(peaks[2 + std::size(windows)]) - 1) < 1e-9);

  // A gap restarts the chain from rest.
  Displacement g(model->highpass);
  CHECK(g.push(5.0) == 0.0 && g.push(7.0).has_value());
  CHECK(!g.push(std::nullopt).has_value() && g.run_length() == 0);
  CHECK(*g.push(1e6) == 0.0);

  // --- the event estimator ------------------------------------------------
  std::map<int, std::vector<PdObservation>> cases;
  std::map<int, std::pair<double, double>> expected;
  for (const auto &r : csv_rows(root + "/data/fixtures/pd_estimate.csv")) {
    const int c = std::stoi(r[0]);
    cases[c].push_back({r[1], 5.0, std::stod(r[3]), std::stod(r[4]), std::stod(r[2])});
    expected[c] = {std::stod(r[5]), std::stod(r[6])};
  }
  for (const auto &[c, obs] : cases) {
    const auto got = pd_magnitude(*model, obs);
    const auto [mean, sd] = expected[c];
    CHECK(got.has_value() == !std::isnan(mean));
    if (!got)
      continue;
    std::println("  case {}: M{:.3f} ± {:.3f} (onset M{:.3f} ± {:.3f}), {} of {} above noise",
                 c, got->mean, got->sd, mean, sd, got->above, got->stations);
    CHECK(std::abs(got->mean - mean) < 1e-6 && std::abs(got->sd - sd) < 1e-6);
    const auto &cp = got->class_p;
    CHECK(std::abs(cp[0] + cp[1] + cp[2] - 1) < 1e-12);
    CHECK(std::ranges::all_of(cp, [](double p) { return p >= 0 && p <= 1; }));
  }
  // Size classes: values at the relation's own prediction for M 2.5 and M 6.5,
  // at three stations 30 km away and at their median noise, fall almost
  // wholly in the small and the large class.
  for (const double m : {2.5, 6.5}) {
    std::vector<PdObservation> obs;
    for (const auto *st : {"ADVT", "AFSR", "AGRB"})
      obs.push_back({st, 5.0,
                     std::pow(10.0, model->window(5.0)->mean(m, 30.0) +
                                        model->term(st, 5.0)),
                     model->noise.at(st), 30.0});
    const auto got = pd_magnitude(*model, obs);
    CHECK(got.has_value());
    std::println("  M{:.1f} at 30 km: M{:.2f}, classes {:.3f} {:.3f} {:.3f}", m,
                 got->mean, got->class_p[0], got->class_p[1], got->class_p[2]);
    CHECK(got->class_p[m < 4 ? 0 : 2] > 0.99);
  }
  // No station above the threshold: no estimate.
  CHECK(!pd_magnitude(*model, {{"ADVT", 5.0, 1e-7, 1e-7, 50.0}}).has_value());

  // --- the network stage --------------------------------------------------
  // Three stations declare an event; their PdEstimates, at the distances of
  // their geometry estimates (the event is not located), give the estimator's
  // value for those distances. A value for another trigger is not counted.
  {
    using namespace ayzek::pipeline;
    Log::get().quiet = true;
    const auto &c0 = cases.at(0);
    std::map<std::string, StationInfo> st;
    for (std::size_t i = 0; i < c0.size(); ++i)
      st[c0[i].station] = {39.0 + 0.3 * static_cast<double>(i), 28.0};
    NetworkConfig cfg;
    cfg.min_stations = 3;
    cfg.pd = &*model;
    Network net(st, cfg);
    double t = 1000.0;
    for (const auto &o : c0) {
      net.on(Detection{o.station, t, t + 4.0, 0.99f, 0.0});
      net.on(StationGeometry{o.station, t, t + 3.5, 5.0, std::log(o.distance_km), 0.2,
                             t + 4.5});
      net.on(PdEstimate{o.station, t, t + 3.5, 5.0, o.pd, o.pd_noise, t + 9.0});
      net.on(PdEstimate{o.station, t + 100.0, t + 103.5, 5.0, 1.0, 1e-7, t + 9.0});
      t += 0.5;
    }
    CHECK(net.declared_count() == 1);
    const auto &e = net.events().front();
    CHECK(e.pd_magnitude.has_value() && e.pd.size() == 3);
    const auto direct = pd_magnitude(*model, c0);
    std::println("  network: M{:.3f} ± {:.3f} from {} stations (estimator M{:.3f})",
                 e.pd_magnitude->mean, e.pd_magnitude->sd, e.pd_magnitude->stations,
                 direct->mean);
    CHECK(std::abs(e.pd_magnitude->mean - direct->mean) < 1e-9);
    CHECK(e.pd_first.has_value());
  }
  // --- silent stations (NetworkConfig::pd_silent_bounds) ---------------------
  // Five stations around a source trigger and report Pd as the relation
  // predicts for M 3.5; a sixth records only noise (PdTraces at its median
  // noise level). Once the P times locate the event, the sixth enters as an
  // upper bound and lowers the estimate; without its traces it does not.
  {
    using namespace ayzek::pipeline;
    Log::get().quiet = true;
    const double lat0 = 39.0, lon0 = 28.0, origin = 1000.0, m_true = 3.5;
    const std::vector<std::pair<std::string, std::pair<double, double>>> sites = {
        {"ADVT", {39.3, 28.0}}, {"AFSR", {38.7, 28.1}}, {"AGRB", {39.0, 28.4}},
        {"AKDG", {39.1, 27.6}}, {"AKDM", {38.8, 27.7}}, {"ANTB", {39.2, 28.3}}};
    auto run = [&](bool traces) {
      std::map<std::string, StationInfo> st;
      for (const auto &[c, ll] : sites)
        st[c] = {ll.first, ll.second};
      NetworkConfig cfg;
      cfg.pd = &*model;
      cfg.p_sigma = 0.2; // exact P times: a radius within pd_max_err_km
      Network net(st, cfg);
      // The silent station's traces, from 60 s before the origin to 30 s after.
      if (traces)
        for (double t = origin - 60; t < origin + 30; t += 0.5)
          net.on(PdTrace{"ANTB", t, t + 0.5, model->noise.at("ANTB"), origin - 300, t + 0.5});
      const auto *w = model->window(5.0);
      for (std::size_t i = 0; i + 1 < sites.size(); ++i) {
        const auto &[c, ll] = sites[i];
        const double d = distance_km(lat0, lon0, ll.first, ll.second);
        const double p = origin + std::hypot(d, 10.0) / 6.0;
        Detection det{c, p - 3.5, p + 4.0, 0.99f, 0.0};
        det.p_time = p;
        net.on(det);
        net.on(PdEstimate{c, p - 3.5, p, 5.0,
                          std::pow(10.0, w->mean(m_true, d) + model->term(c, 5.0)),
                          model->noise.at(c), p + 5.0});
      }
      // Reported after every value is in: the traces arrived first here.
      const auto &e = net.events().front();
      CHECK(e.declared && e.p_location && e.p_location->err_km <= cfg.pd_max_err_km);
      return *e.pd_magnitude;
    };
    const auto without = run(false), with = run(true);
    std::println("  silent station: M{:.2f} from {} stations without its bound, "
                 "M{:.2f} from {} with it",
                 without.mean, without.stations, with.mean, with.stations);
    CHECK(without.stations == 5 && with.stations == 6);
    CHECK(with.mean < without.mean);
  }
  std::println("Pd chain and estimator agree with onset");
}
