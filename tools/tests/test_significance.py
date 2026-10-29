"""tools/significance.py: exact McNemar, Clopper-Pearson, the paired bootstrap
of chance-corrected differences (event and block), the equivalence verdict,
--outcome, and the quiet-window comparison.

Expected values: scipy.stats (binomtest, its exact proportion interval,
binom quantiles, special.betainc) as independent implementations; hand
arithmetic; bootstrap distributions that are known in closed form for the
constructed cases (a binomial count of the draws of one event); invariants.
"""
import csv
import json
import math
import random

import pytest
from scipy import special, stats

import significance as sg
from conftest import run_tool


def origin(i):
    """Origin time of the i-th event: one per minute, 5 s past it."""
    return f"{i // 60:02d}:{i % 60:02d}:05.00"


def events(ds, n, start=0):
    return [(ds, origin(start + i), "ML", "2.0") for i in range(n)]


def outcomes(evs, detected):
    """{event: 1 if its index is in `detected` else 0}."""
    return {e: int(i in detected) for i, e in enumerate(evs)}


# --- exact tests and intervals -----------------------------------------------------

@pytest.mark.parametrize("b, c", [(0, 1), (1, 0), (0, 5), (5, 0), (3, 3), (1, 9), (10, 20),
                                  (40, 25), (0, 50), (150, 170), (300, 280)])
def test_mcnemar_matches_scipy(b, c):
    # Exact McNemar: two-sided binomial test of b discordant pairs out of b + c
    # at p = 1/2.
    assert sg.mcnemar(b, c) == pytest.approx(stats.binomtest(b, b + c, 0.5).pvalue,
                                             rel=1e-9, abs=1e-300)


def test_mcnemar_no_discordant_pairs():
    assert sg.mcnemar(0, 0) == 1.0


@pytest.mark.parametrize("b, c", [(510, 530), (954, 1046), (0, 1040), (1040, 0), (2400, 2600),
                                  (2500, 2500), (0, 5000), (2300, 2700)])
def test_mcnemar_large_counts(b, c):
    # Alarm totals of test networks exceed 1000 (development: 954 for one
    # detector); scipy's exact binomial test is the reference, to n = 5000.
    assert sg.mcnemar(b, c) == pytest.approx(stats.binomtest(b, b + c, 0.5).pvalue, rel=1e-9, abs=1e-300)


@pytest.mark.parametrize("k, n, p", [(3, 10, 0.3), (0, 50, 0.01), (50, 50, 0.99), (2000, 5000, 0.4),
                                     (4999, 5000, 0.999), (1, 3000, 0.0002)])
def test_binom_cdf_general_p_matches_scipy(k, n, p):
    assert sg.binom_cdf(k, n, p) == pytest.approx(stats.binom.cdf(k, n, p), rel=1e-8, abs=1e-300)


def test_binom_cdf_exact_cases_by_hand():
    # P(X <= 1) for Binomial(3, 1/2) = (1 + 3) / 8; P(X <= n) = 1; P(X <= 0) = 2^-n.
    assert sg.binom_cdf(1, 3) == 0.5
    assert sg.binom_cdf(3, 3) == 1.0 and sg.binom_cdf(0, 10) == 2 ** -10
    assert sg.binom_cdf(5000, 5000) == 1.0 and sg.binom_cdf(0, 5000) == 2.0 ** -5000


@pytest.mark.parametrize("k, n", [(0, 1), (1, 1), (0, 10), (3, 10), (10, 10), (8, 9), (4, 12),
                                  (1, 300), (50, 197), (128, 197), (196, 197)])
@pytest.mark.parametrize("alpha", [0.05, 0.10])
def test_clopper_pearson_matches_scipy(k, n, alpha):
    ci = stats.binomtest(k, n).proportion_ci(confidence_level=1 - alpha, method="exact")
    lo, hi = sg.clopper_pearson(k, n, alpha)
    assert lo == pytest.approx(ci.low, abs=1e-9)
    assert hi == pytest.approx(ci.high, abs=1e-9)


@pytest.mark.parametrize("a", [0.5, 1, 2, 5, 30, 150])
@pytest.mark.parametrize("b", [0.5, 1, 3, 40, 200])
@pytest.mark.parametrize("x", [1e-4, 0.05, 0.3, 0.5, 0.7, 0.95, 0.9999])
def test_betainc_matches_scipy(a, b, x):
    assert sg.betainc(a, b, x) == pytest.approx(special.betainc(a, b, x), abs=1e-10)


def test_pct_by_hand():
    assert sg.pct([5, 1, 3], 0.5) == 3
    assert sg.pct([7], 0.025) == sg.pct([7], 0.975) == 7
    v = list(range(100))
    # rank round(q * 99): 0.025 * 99 = 2.475 -> 2; 0.975 * 99 = 96.525 -> 97
    assert (sg.pct(v, 0.025), sg.pct(v, 0.975), sg.pct(v, 0), sg.pct(v, 1)) == (2, 97, 0, 99)
    # With B = 10000 resamples the 95% interval leaves 250 values strictly
    # below and 250 strictly above it, the 90% one 500 and 500.
    v = list(range(sg.B))
    lo, hi = sg.pct(v, 0.025), sg.pct(v, 0.975)
    assert (sum(x < lo for x in v), sum(x > hi for x in v)) == (250, 250)
    lo, hi = sg.pct(v, 0.05), sg.pct(v, 0.95)
    assert (sum(x < lo for x in v), sum(x > hi for x in v)) == (500, 500)


# --- chance-corrected detections -----------------------------------------------------

def test_corrected_by_hand():
    # Replay s1: 10 events, 6 detected; its shifted catalogues 1/10, 2/10,
    # 0/10, 1/10 -> c = 0.4 / 4 = 0.1, corrected (0.6 - 0.1) / 0.9 * 10 = 50/9.
    # Replay s2: 5 events, 4 detected; shifted 1/5, 0/5, 0/5, 1/5 -> c = 0.1,
    # corrected (0.8 - 0.1) / 0.9 * 5 = 35/9. Summed per replay: 85/9; raw 10.
    e1, e2 = events("s1", 10), events("s2", 5)
    obs = {**outcomes(e1, range(6)), **outcomes(e2, range(4))}
    sh = {127: {**outcomes(e1, {0}), **outcomes(e2, {0})},
          -127: {**outcomes(e1, {0, 1}), **outcomes(e2, set())},
          263: {**outcomes(e1, set()), **outcomes(e2, set())},
          -263: {**outcomes(e1, {3}), **outcomes(e2, {2})}}
    total, raw = sg.corrected(obs, sh, sorted(obs), {s: sorted(v) for s, v in sh.items()})
    assert total == pytest.approx(85 / 9, abs=1e-12)
    assert raw == 10


# A chance level of 1 (every shifted-catalogue event of a replay detected, in
# the data of a short replay or in one resample): tools/sweep_trigger.py:181
# reports 0.0 corrected detections for that replay, and significance.py:141
# follows it since 943e031 (before, it divided by 1 - c).


def test_chance_level_one_as_sweep_trigger():
    e = events("s1", 2)
    obs = outcomes(e, {0, 1})
    sh = {127: outcomes(e, {0, 1}), -127: outcomes(e, {0, 1})}
    assert sg.corrected(obs, sh, e, {s: e for s in sh})[0] == 0.0


def test_resample_without_shifted_events_is_chance_level_zero():
    # 943e031: a resample that drew no shifted event of a replay has c = 0, so
    # the corrected count is the raw one: 3 of 4 detected -> 0.75 * 4 = 3.
    e = events("s1", 4)
    sh = {127: outcomes(e, {0}), -127: outcomes(e, {1})}
    total, raw = sg.corrected(outcomes(e, {0, 1, 2}), sh, e, {127: [], -127: []})
    assert (total, raw) == (3.0, 3)


def test_compare_with_a_replay_at_chance_level_one():
    # Group of two replays. s1 (3 events) is short: every event is detected at
    # every shift, c = 1 in the data and in every resample; it contributes 0
    # to both detectors by sweep_trigger's rule. s2 (4 events) has c = 0:
    # learned 3, other 2. So corrected 3 vs 2, difference 1, and the paired
    # bootstrap difference is X draws of s2's extra event, X ~ Binomial(4, 1/4):
    # P(X=0) = 0.3164, P(X<=2) = 0.9492, P(X<=3) = 0.9961 -> 95% interval (0, 3).
    e1, e2 = events("s1", 3), events("s2", 4, start=100)
    all1 = outcomes(e1, range(3))
    det = {"learned": {0: {**all1, **outcomes(e2, range(3))},
                       127: {**all1, **outcomes(e2, set())},
                       -127: {**all1, **outcomes(e2, set())}},
           "other": {0: {**all1, **outcomes(e2, range(2))},
                     127: {**all1, **outcomes(e2, set())},
                     -127: {**all1, **outcomes(e2, set())}}}
    r = sg.compare(det, "learned", "other", {"s1", "s2"}, random.Random(5))
    assert (r["a_corrected"], r["b_corrected"], r["difference"]) == (3.0, 2.0, 1.0)
    assert (r["a_detected"], r["b_detected"], r["events"]) == (6, 5, 7)
    assert r["diff_ci"] == pytest.approx((0, 3))
    assert r["p_boot_le0"] == pytest.approx(0.75 ** 4, abs=0.02)
    assert stats.binom(4, 0.25).ppf(0.975) == 3


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=(
    "A shifted catalogue with no scorable event in a replay has no rows, and "
    "significance averages c over the other shifts only (significance.py:132-140, the `if cs` fallback to 0 covering only a resample with no shifted event at all; "
    "compare() takes the shifts from the rows), while tools/sweep_trigger.py:179-180 "
    "counts it as a chance level of 0 (detected / max(1, 0)). Here c = 0.05 against "
    "sweep_trigger's 0.025, corrected 9.47 against 9.74. Which rule is meant is for the "
    "user; unlikely to arise on a 3 h sequence replay, which keeps events at every shift."))
def test_shift_without_events_as_sweep_trigger():
    # Replay of 20 events, 10 detected; +127: 1 of 20; -127: no event at all
    # (the shifted catalogue is empty for this replay, so the table has no
    # rows for it). sweep_trigger: c = (1/20 + 0/max(1, 0)) / 2 = 0.025,
    # corrected = (0.5 - 0.025) / 0.975 * 20 = 9.7436.
    e = events("s1", 20)
    det = {"learned": {0: outcomes(e, range(10)), 127: outcomes(e, {0})},
           "other": {0: outcomes(e, range(10)), 127: outcomes(e, {0})}}
    r = sg.compare(det, "learned", "other", {"s1"}, random.Random(1))
    assert r["a_corrected"] == pytest.approx((0.5 - 0.025) / 0.975 * 20)


# --- the paired bootstrap ---------------------------------------------------------------

def test_identical_detectors_every_resample_differs_by_zero():
    # Invariant of a paired bootstrap: with the same draws for both detectors,
    # two identical detectors differ by exactly 0 in every resample, whatever
    # their outcomes and chance levels.
    e1, e2 = events("s1", 15), events("s2", 8, start=100)
    obs = {**outcomes(e1, {0, 2, 3, 7, 11, 12}), **outcomes(e2, {1, 4, 5})}
    sh = {s: {**outcomes(e1, {s % 7, 9}), **outcomes(e2, {s % 3})} for s in (127, -127, 263, -263)}
    det = {"learned": {0: obs, **sh}, "other": {0: dict(obs), **{s: dict(v) for s, v in sh.items()}}}
    r = sg.compare(det, "learned", "other", {"s1", "s2"}, random.Random(7))
    assert r["difference"] == 0
    assert r["diff_ci"] == (0.0, 0.0) and r["diff_ci90"] == (0.0, 0.0)
    assert r["p_boot_le0"] == 1.0
    assert (r["a_only"], r["b_only"], r["mcnemar_p"]) == (0, 0, 1.0)
    assert r["a_cp"] == r["b_cp"]


@pytest.mark.parametrize("chance", [0.0, 0.5])
def test_one_extra_detection_bootstrap_is_binomial(chance):
    # One replay of n = 20 events; learned detects events 0-10, the other 0-9.
    # Shifted catalogues: with chance 0, two shifts with nothing detected; with
    # chance 0.5, three shifts with everything detected and three with nothing,
    # so that every resample has c = 0.5 exactly. Then corrected =
    # (k/n - c) / (1 - c) * n = (k - c n) / (1 - c), and the paired difference
    # in a resample is X / (1 - c), X the number of draws of event 10:
    # X ~ Binomial(20, 1/20). Point estimate: 1 / (1 - c).
    n = 20
    e = events("s1", n)
    if chance == 0.0:
        sh = {127: outcomes(e, set()), -127: outcomes(e, set())}
    else:
        sh = {s: outcomes(e, range(n) if s > 0 else set()) for s in (127, -127, 263, -263, 391, -391)}
    det = {"learned": {0: outcomes(e, range(11)), **sh},
           "other": {0: outcomes(e, range(10)), **{s: dict(v) for s, v in sh.items()}}}
    r = sg.compare(det, "learned", "other", {"s1"}, random.Random(20261022))
    scale = 1 / (1 - chance)
    x = stats.binom(n, 1 / n)
    # Binomial(20, 0.05): P(X<=2) = 0.9245, P(X<=3) = 0.9841, so the 95th and
    # 97.5th percentiles are 3 and the 2.5th and 5th are 0, each more than 7
    # standard errors from the next value at B = 10000.
    assert (x.ppf(0.025), x.ppf(0.975), x.ppf(0.05), x.ppf(0.95)) == (0, 3, 0, 3)
    assert r["diff_ci"] == pytest.approx((0, 3 * scale))
    assert r["diff_ci90"] == pytest.approx((0, 3 * scale))
    # P(diff <= 0) = P(X = 0) = 0.95^20 = 0.3585; binomial s.e. at B = 10000
    # is 0.0048.
    assert r["p_boot_le0"] == pytest.approx(x.pmf(0), abs=0.02)
    assert r["difference"] == pytest.approx(scale)
    assert (r["a_detected"], r["b_detected"], r["events"]) == (11, 10, 20)
    assert (r["a_only"], r["b_only"]) == (1, 0)
    assert r["mcnemar_p"] == pytest.approx(stats.binomtest(1, 1, 0.5).pvalue)


def test_same_shift_draws_for_both_detectors():
    # The two detectors share the shifted outcomes (events 0-3 at +127, 0-1 at
    # -127), and the learned one detects one more observed event. If both use
    # the same draws of the shifted events, c is the same for both in every
    # resample and the difference is X / (1 - c) >= 0, zero exactly when the
    # extra event is not drawn: P(diff <= 0) = 0.95^20 = 0.3585, and the
    # interval starts at 0. Separate draws would give negative differences.
    n = 20
    e = events("s1", n)
    sh = {127: outcomes(e, range(4)), -127: outcomes(e, range(2))}
    det = {"learned": {0: outcomes(e, range(1, 9)), **sh},
           "other": {0: outcomes(e, range(1, 8)), **{s: dict(v) for s, v in sh.items()}}}
    r = sg.compare(det, "learned", "other", {"s1"}, random.Random(3))
    assert r["diff_ci"][0] == 0.0
    assert r["p_boot_le0"] == pytest.approx(0.95 ** 20, abs=0.02)


def test_resample_within_each_replay():
    # Every replay keeps its size, and draws come only from that replay.
    evs = events("s1", 7) + events("s2", 3, start=50)
    rng = random.Random(5)
    for _ in range(200):
        out = sg.resample(evs, rng, 0)
        assert len(out) == 10
        assert sum(e[0] == "s1" for e in out) == 7


def test_block_resample_draws_whole_blocks():
    # Ten one-minute blocks of two events (origins mm:05 and mm:35): with
    # 60 s blocks a resample is ten whole blocks, so both events of a block
    # are drawn equally often.
    evs = [("s1", f"00:{m:02d}:{s}.00", "ML", "2.0") for m in range(10) for s in ("05", "35")]
    rng = random.Random(11)
    for _ in range(200):
        out = sg.resample(evs, rng, 60)
        assert len(out) == 20
        for m in range(10):
            assert out.count(evs[2 * m]) == out.count(evs[2 * m + 1])


def test_block_boundaries_by_hand():
    assert sg.seconds("01:02:03.50") == 3723.5
    # 00:00:59.99 is in the first 60 s block, 00:01:00.00 in the second.
    evs = [("s1", "00:00:59.99", "ML", "1.0"), ("s1", "00:01:00.00", "ML", "1.0")]
    rng = random.Random(2)
    singles = [sg.resample(evs, rng, 60) for _ in range(100)]
    assert all(len(o) == 2 for o in singles)
    # With 120 s blocks both are one block, always drawn together.
    assert all(sg.resample(evs, rng, 120) == evs for _ in range(100))


@pytest.mark.parametrize("block, want_hi, want_p0", [
    # Blocks (60 s): the learned detector alone detects both events of block 0;
    # diff = 2 * (draws of block 0), Binomial(10, 0.1): P(<=2) = 0.9298,
    # P(<=3) = 0.9872 -> 97.5th percentile 2 * 3 = 6; P(0) = 0.9^10 = 0.3487.
    (60, 6, 0.9 ** 10),
    # Single events: diff = draws of either event, Binomial(20, 0.1):
    # P(<=4) = 0.9568, P(<=5) = 0.9887 -> 5; P(0) = 0.9^20 = 0.1216.
    (0, 5, 0.9 ** 20),
])
def test_block_bootstrap_distribution(block, want_hi, want_p0):
    evs = [("s1", f"00:{m:02d}:{s}.00", "ML", "2.0") for m in range(10) for s in ("05", "35")]
    base = {e: int(i % 4 == 2) for i, e in enumerate(evs)}       # shared: events 2, 6, 10, 14, 18
    extra = {**base, evs[0]: 1, evs[1]: 1}                       # plus both of block 0
    sh = {127: {e: 0 for e in evs}, -127: {e: 0 for e in evs}}
    det = {"learned": {0: extra, **sh}, "other": {0: base, **{s: dict(v) for s, v in sh.items()}}}
    r = sg.compare(det, "learned", "other", {"s1"}, random.Random(99), block)
    assert r["diff_ci"] == pytest.approx((0, want_hi))
    assert r["p_boot_le0"] == pytest.approx(want_p0, abs=0.02)
    assert stats.binom(10 if block else 20, 0.1).ppf(0.975) * (2 if block else 1) == want_hi


# --- load_outcomes, main ---------------------------------------------------------------------

FIELDS = ["dataset", "detector", "shift", "origin", "type", "magnitude", "detected", "detected_pre_s"]


def write_table(path, detector, rows):
    """rows: (dataset, shift, origin, detected, detected_pre_s) as event_outcomes writes them."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(FIELDS)
        for ds, shift, o, d, p in rows:
            w.writerow([ds, detector, shift, o, "ML", "2.0", d, p])


def extra_detection_tables(d, pre_s_learned=None):
    """The one-extra-detection case of test_one_extra_detection_bootstrap_is_binomial
    (chance 0) as event_outcomes tables: learned 0-10 of 20, other 0-9."""
    pre = pre_s_learned or (lambda i: int(i < 11))
    write_table(d / "learned.csv", "learned",
                [("s1", 0, origin(i), int(i < 11), pre(i)) for i in range(20)]
                + [("s1", s, origin(i), 0, 0) for s in (127, -127) for i in range(20)])
    write_table(d / "other.csv", "other",
                [("s1", 0, origin(i), int(i < 10), int(i < 10)) for i in range(20)]
                + [("s1", s, origin(i), 0, 0) for s in (127, -127) for i in range(20)])


def test_load_outcomes(tmp_path):
    write_table(tmp_path / "a.csv", "learned", [("s1", 0, "00:00:05.00", 1, ""),
                                                ("s1", 127, "00:00:05.00", 0, 0)])
    det = sg.load_outcomes(tmp_path)
    assert det["learned"][0] == {("s1", "00:00:05.00", "ML", "2.0"): 1}
    assert det["learned"][127] == {("s1", "00:00:05.00", "ML", "2.0"): 0}
    # A table from a build before the pre-S line has no detected_pre_s values.
    with pytest.raises(SystemExit):
        sg.load_outcomes(tmp_path, "detected_pre_s")


def run_main(d, *args):
    p = run_tool("significance", "--outcomes", d, "--against", "other", "--groups", "g=s1",
                 "--json", d / "out.json", *args)
    assert p.returncode == 0, p.stderr
    return p, json.loads((d / "out.json").read_text())["summary"]


def test_main_summary(tmp_path):
    extra_detection_tables(tmp_path)
    p, s = run_main(tmp_path)
    # By hand (c = 0): corrected = detections; 11 vs 10, difference 1; the 95%
    # interval (0, 3) as in test_one_extra_detection_bootstrap_is_binomial;
    # McNemar on 1 discordant pair: p = 1.
    assert s == {"g: learned corrected": 11.0, "g: other detected": 10, "g: other corrected": 10.0,
                 "g: other difference": 1.0, "g: other difference low": 0.0,
                 "g: other difference high": 3.0, "g: other mcnemar p": 1.0}
    assert "learned vs other: 20 events" in p.stdout


@pytest.mark.parametrize("fraction, equivalent", [
    # 90% interval (0, 3) (Binomial(20, 0.05): 5th percentile 0, 95th 3);
    # margin m = fraction * 20 events; equivalent iff -m < 0 and 3 < m.
    (0.0, False),     # m = 0: -0 < 0 fails
    (0.15, False),    # m = 3.0: 3 < 3 fails (the interval reaches the margin)
    (0.16, True),     # m = 3.2
    (0.5, True),      # m = 10
])
def test_equivalence_verdict_flips_at_the_margin(tmp_path, fraction, equivalent):
    extra_detection_tables(tmp_path)
    p, s = run_main(tmp_path, "--equivalence", f"other={fraction}")
    assert s["g: other difference 90% low"] == 0.0
    assert s["g: other difference 90% high"] == 3.0
    assert s["g: other equivalence margin"] == round(fraction * 20, 1)
    assert s["g: other equivalent"] is equivalent
    assert ("equivalent (both one-sided tests p < 0.05)" in p.stdout) is equivalent


def test_outcome_detected_pre_s(tmp_path):
    # The learned detector detects events 0-10 but only 0-7 before S; the other
    # detects 0-9, all before S. --outcome detected_pre_s compares 8 with 10.
    extra_detection_tables(tmp_path, pre_s_learned=lambda i: int(i < 8))
    _, s = run_main(tmp_path, "--outcome", "detected_pre_s")
    assert (s["g: learned corrected"], s["g: other corrected"], s["g: other difference"]) == (8.0, 10.0, -2.0)
    _, s = run_main(tmp_path)
    assert s["g: other difference"] == 1.0


def test_outcome_detected_pre_s_needs_the_column(tmp_path):
    extra_detection_tables(tmp_path, pre_s_learned=lambda i: "")
    p = run_tool("significance", "--outcomes", tmp_path, "--against", "other", "--groups", "g=s1",
                 "--outcome", "detected_pre_s")
    assert p.returncode != 0 and "no detected_pre_s column" in p.stderr


def test_seeded_runs_reproduce(tmp_path):
    # Invariant: the seed is fixed, so two runs give the same figures, also
    # under different string-hash seeds (set iteration order). Three replays in
    # one group, varied outcomes and chance levels.
    rows_l, rows_o = [], []
    rng = random.Random(4)
    for ds in ("s1", "s2", "s3"):
        for s in (0, 127, -127, 263, -263):
            for i in range(rng.randint(8, 15)):
                rows_l.append((ds, s, origin(i), int(rng.random() < (0.7 if s == 0 else 0.1)), 0))
                rows_o.append((ds, s, origin(i), int(rng.random() < (0.6 if s == 0 else 0.1)), 0))
    # Both detectors need the same events per replay and shift.
    keys = {}
    for r in rows_l:
        keys[(r[0], r[1], r[2])] = r
    rows_o = [r for r in rows_o if (r[0], r[1], r[2]) in keys]
    seen = {(r[0], r[1], r[2]) for r in rows_o}
    rows_l = [r for r in rows_l if (r[0], r[1], r[2]) in seen]
    write_table(tmp_path / "learned.csv", "learned", rows_l)
    write_table(tmp_path / "other.csv", "other", rows_o)
    outs = []
    for hashseed in ("0", "1", "2", "3"):
        p = run_tool("significance", "--outcomes", tmp_path, "--against", "other",
                     "--groups", "g=s1,s2,s3", "--json", tmp_path / f"{hashseed}.json",
                     "--equivalence", "other=0.1", env={"PYTHONHASHSEED": hashseed})
        assert p.returncode == 0, p.stderr
        outs.append(((tmp_path / f"{hashseed}.json").read_text(), p.stdout))
    assert all(o == outs[0] for o in outs)


# --- quiet windows ------------------------------------------------------------------------------

def test_quiet_by_hand():
    # Two windows of 5.5 h; learned 3 + 1 = 4 alarms, stalta 6 + 2 = 8.
    hours = {"w1": 5.5, "w2": 5.5}
    q = sg.quiet(hours, {"learned": {"w1": 3, "w2": 1}, "stalta": {"w1": 6, "w2": 2}},
                 random.Random(1))
    assert q["hours"] == 11.0
    assert (q["learned"], q["stalta"]) == (4, 8)
    assert q["learned_per_day"] == pytest.approx(24 * 4 / 11)
    assert q["stalta_per_day"] == pytest.approx(24 * 8 / 11)
    assert q["ratio"] == 0.5
    # Conditional on 12 alarms, 4 ~ Binomial(12, 1/2) under equal rates.
    assert q["p_equal_rates"] == pytest.approx(stats.binomtest(4, 12, 0.5).pvalue, rel=1e-9)
    ci = stats.binomtest(4, 12).proportion_ci(method="exact")
    assert q["ratio_exact_ci"] == pytest.approx((ci.low / (1 - ci.low), ci.high / (1 - ci.high)))
    # Bootstrap over the two windows: a resample is (w1, w1), (w1, w2) or
    # (w2, w2), with probabilities 1/4, 1/2, 1/4: learned per day
    # 24*6/11, 24*4/11, 24*2/11, so the 2.5th and 97.5th percentiles are the
    # extremes; stalta is twice learned in every window, so the ratio is 0.5
    # in every resample.
    assert q["learned_per_day_boot_ci"] == pytest.approx((24 * 2 / 11, 24 * 6 / 11))
    assert q["stalta_per_day_boot_ci"] == pytest.approx((24 * 4 / 11, 24 * 12 / 11))
    assert q["ratio_boot_ci"] == pytest.approx((0.5, 0.5))


def test_quiet_second_detector_without_alarms():
    # stalta raises no alarm: conditional on the 3 alarms, all are learned's,
    # p = P(X in {0, 3}) for X ~ Binomial(3, 1/2) = 2/8 = 0.25; the ratio of
    # rates is undefined and reported as inf, the exact interval's upper bound
    # inf and its lower bound lo/(1 - lo) with lo the Clopper-Pearson lower
    # bound of 3 of 3 (scipy); per-day rates 24*3/11 and 0; the bootstrap's
    # ratio is inf in every resample.
    q = sg.quiet({"w1": 5.5, "w2": 5.5}, {"learned": {"w1": 2, "w2": 1}, "stalta": {"w1": 0, "w2": 0}},
                 random.Random(1))
    assert q["p_equal_rates"] == pytest.approx(0.25) == pytest.approx(stats.binomtest(3, 3, 0.5).pvalue)
    assert q["ratio"] == math.inf
    lo = stats.binomtest(3, 3).proportion_ci(method="exact").low
    assert q["ratio_exact_ci"] == pytest.approx((lo / (1 - lo), math.inf))
    assert (q["learned_per_day"], q["stalta_per_day"]) == (pytest.approx(24 * 3 / 11), 0.0)
    assert q["ratio_boot_ci"] == (math.inf, math.inf)
    assert q["stalta_per_day_boot_ci"] == (0.0, 0.0)


def test_quiet_first_detector_without_alarms():
    # learned raises none of the 4: p = 2/16 = 0.125, ratio 0, exact interval
    # (0, hi/(1 - hi)) with hi the upper bound of 0 of 4.
    q = sg.quiet({"w1": 5.5, "w2": 5.5}, {"learned": {"w1": 0, "w2": 0}, "stalta": {"w1": 3, "w2": 1}},
                 random.Random(1))
    assert q["p_equal_rates"] == pytest.approx(0.125)
    hi = stats.binomtest(0, 4).proportion_ci(method="exact").high
    assert q["ratio"] == 0.0 and q["ratio_exact_ci"] == pytest.approx((0.0, hi / (1 - hi)))


def test_quiet_no_alarms_at_all():
    # No alarm by either: p = 1 (no discordant count), ratio nan, interval
    # (0, inf), nothing to bootstrap but zeros.
    q = sg.quiet({"w1": 5.5}, {"learned": {"w1": 0}, "stalta": {"w1": 0}}, random.Random(1))
    assert q["p_equal_rates"] == 1.0 and math.isnan(q["ratio"])
    assert q["ratio_exact_ci"] == (0.0, math.inf)
    assert q["learned_per_day_boot_ci"] == (0.0, 0.0)


def test_main_quiet_without_alarms_reports_inf(tmp_path):
    # Through main(): the printed ratio and the JSON stay well-formed (Python's
    # json writes Infinity, which json.load reads back).
    index, run = quiet_dirs(tmp_path, {"w1": "  Alarms raised                           3\n"})
    other = tmp_path / "other"
    other.mkdir()
    (other / "w1.txt").write_text("  Alarms raised                           0\n")
    outcomes = tmp_path / "outcomes"
    outcomes.mkdir()
    p = run_tool("significance", "--outcomes", outcomes, "--groups", "g=none", "--against", "stalta_eqfa",
                 "--quiet", index, "--quiet-runs", f"learned={run}", f"stalta={other}",
                 "--json", tmp_path / "s.json")
    assert p.returncode == 0, p.stderr
    assert "  ratio                        inf\n" in p.stdout
    assert "  ratio_exact_ci               " in p.stdout and " to inf\n" in p.stdout
    js = json.loads((tmp_path / "s.json").read_text())["summary"]
    assert js["quiet: stalta ratio"] == math.inf and js["quiet: stalta p"] == 0.25
    assert (js["quiet: learned alarms"], js["quiet: stalta alarms"]) == (3, 0)


def quiet_dirs(tmp_path, reports):
    index = tmp_path / "quiet_index.csv"
    index.write_text("name,scored_hours\n" + "".join(f"{n},5.5\n" for n in reports))
    run = tmp_path / "run"
    run.mkdir()
    for n, text in reports.items():
        (run / f"{n}.txt").write_text(text)
    return index, run


def test_quiet_counts(tmp_path):
    index, run = quiet_dirs(tmp_path, {"w1": "  Alarms raised                           3\n",
                                       "w2": "  Alarms raised                           0\n"})
    hours, counts = sg.quiet_counts(index, {"learned": run})
    assert hours == {"w1": 5.5, "w2": 5.5}
    assert counts == {"learned": {"w1": 3, "w2": 0}}


@pytest.mark.parametrize("text", ["ayzek    error: cannot open w2/KULA.mseed\n", "",
                                  "  Alarms compared with the AFAD catalogue\n"])
def test_quiet_counts_report_without_summary_stops(tmp_path, text):
    # A report without the "Alarms raised" line (a run that crashed or was cut
    # short) stops with the file named, rather than counting as 0 alarms.
    index, run = quiet_dirs(tmp_path, {"w1": "  Alarms raised                           3\n", "w2": text})
    with pytest.raises(SystemExit) as e:
        sg.quiet_counts(index, {"learned": run})
    assert str(run / "w2.txt") in str(e.value) and "Alarms raised" in str(e.value)


def test_quiet_counts_missing_report_stops(tmp_path):
    index, run = quiet_dirs(tmp_path, {"w1": "  Alarms raised                           3\n"})
    index.write_text("name,scored_hours\nw1,5.5\nw2,5.5\n")
    with pytest.raises(FileNotFoundError):
        sg.quiet_counts(index, {"learned": run})
