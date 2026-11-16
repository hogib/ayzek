"""tools/heldout_grid_matched.py (exploratory): the fixed rule that picks baseline
settings for a learned grid setting, on small synthetic summaries. No
significance.py run.

Expected values: the rule as the tool's RULE states it (ratio = count(L) /
count(B), matched when 0.8 <= ratio <= 1.25, closest on a log scale, first in
grid order on a tie, last in grid order for bracket_above), arithmetic by hand.
"""
import heldout_grid_matched as gm


def row(setting, unmatched, quiet=0):
    return {"setting": setting, "unmatched_alarms": str(unmatched), "quiet_alarms": str(quiet)}


def names(picks):
    return [(rule, r["setting"], flag, tie) for rule, r, flag, tie in picks]


def test_the_rule_constants():
    assert (gm.LO, gm.HI) == (0.8, 1.25)
    assert gm.FAMILIES == ("stalta", "phasenet", "eqt")
    assert gm.NETS == ("greece_now", "nz_now", "ko_now")


def test_ratio_bounds_are_inclusive():
    # 8 / 10 = 0.8 and 10 / 8 = 1.25: both still matched
    assert gm.pick(row("L", 8), [row("b", 10)], "unmatched_alarms")[1:] == ("matched", False)
    assert gm.pick(row("L", 10), [row("b", 8)], "unmatched_alarms")[1:] == ("matched", False)
    assert gm.pick(row("L", 10), [row("b", 10)], "unmatched_alarms")[1:] == ("matched", False)


def test_outside_the_bounds_the_closest_is_still_returned_with_its_ratio():
    # 79 / 100 = 0.79, 126 / 100 = 1.26
    b = row("b", 100)
    assert gm.pick(row("L", 79), [b], "unmatched_alarms") == (b, "not matched (ratio 0.79)", False)
    assert gm.pick(row("L", 126), [b], "unmatched_alarms") == (b, "not matched (ratio 1.26)", False)


def test_closest_is_on_a_log_scale_not_by_difference():
    # L = 10: 4 is nearer by difference (6 against 12), 22 by ratio (2.2 against 2.5)
    cands = [row("b4", 4), row("b22", 22)]
    got, flag, tie = gm.pick(row("L", 10), cands, "unmatched_alarms")
    assert got["setting"] == "b22" and flag == "not matched (ratio 0.45)" and not tie


def test_tie_on_the_log_scale_takes_the_first_in_grid_order():
    # 10 / 5 = 2 and 10 / 20 = 1/2: the same distance
    cands = [row("b5", 5), row("b20", 20)]
    got, flag, tie = gm.pick(row("L", 10), cands, "unmatched_alarms")
    assert got["setting"] == "b5" and flag == "not matched (ratio 2.00)" and tie
    got, _, tie = gm.pick(row("L", 10), list(reversed(cands)), "unmatched_alarms")
    assert got["setting"] == "b20" and tie


def test_equal_counts_tie_takes_the_first_in_grid_order():
    cands = [row("b3", 3), row("first", 9), row("second", 9)]
    got, flag, tie = gm.pick(row("L", 10), cands, "unmatched_alarms")
    assert got["setting"] == "first" and flag == "matched" and tie       # 10 / 9 = 1.11


def test_baseline_settings_with_a_zero_count_are_not_candidates_for_the_ratio():
    # 5 / 4 = 1.25; the setting with no alarm has no ratio and is passed over, without a tie
    cands = [row("zero", 0), row("b4", 4)]
    got, flag, tie = gm.pick(row("L", 5), cands, "unmatched_alarms")
    assert (got["setting"], flag, tie) == ("b4", "matched", False)


def test_zero_count_rule_when_L_is_zero():
    # closest by absolute difference among all, first in grid order, flagged
    cands = [row("b3", 3), row("b1", 1), row("b1_again", 1)]
    got, flag, tie = gm.pick(row("L", 0), cands, "unmatched_alarms")
    assert (got["setting"], flag, tie) == ("b1", "zero count", True)
    got, flag, tie = gm.pick(row("L", 0), [row("b0", 0), row("b2", 2)], "unmatched_alarms")
    assert (got["setting"], flag, tie) == ("b0", "zero count", False)


def test_zero_count_rule_when_no_baseline_setting_has_a_count():
    cands = [row("a", 0), row("b", 0)]
    got, flag, tie = gm.pick(row("L", 5), cands, "unmatched_alarms")
    assert (got["setting"], flag, tie) == ("a", "zero count", True)


def test_matched_quiet_uses_the_quiet_column():
    L = row("L", 100, quiet=4)
    cands = [row("u100", 100, quiet=40), row("q5", 1, quiet=5)]
    assert gm.pick(L, cands, "unmatched_alarms")[0]["setting"] == "u100"
    got, flag, tie = gm.pick(L, cands, "quiet_alarms")      # 4 / 5 = 0.8
    assert (got["setting"], flag, tie) == ("q5", "matched", False)


def test_bracket_takes_the_nearest_counts_strictly_below_and_above():
    cands = [row("b3", 3), row("b7", 7), row("b12", 12), row("b20", 20)]
    assert names(gm.bracket(row("L", 10), cands)) == [("bracket_below", "b7", "bracket", False),
                                                      ("bracket_above", "b12", "bracket", False)]


def test_bracket_ties_first_below_last_above():
    cands = [row("b7a", 7), row("b7b", 7), row("b12a", 12), row("b12b", 12), row("b20", 20)]
    assert names(gm.bracket(row("L", 10), cands)) == [("bracket_below", "b7a", "bracket", True),
                                                      ("bracket_above", "b12b", "bracket", True)]


def test_bracket_needs_both_sides_and_an_equal_count_is_neither():
    L = row("L", 10)
    assert gm.bracket(L, [row("b12", 12), row("b20", 20)]) == []
    assert gm.bracket(L, [row("b3", 3), row("b7", 7)]) == []
    assert gm.bracket(L, [row("b10", 10), row("b12", 12)]) == []      # 10 is not strictly below
    assert names(gm.bracket(L, [row("b9", 9), row("b10", 10), row("b11", 11)])) == [
        ("bracket_below", "b9", "bracket", False), ("bracket_above", "b11", "bracket", False)]


def test_significance_output_is_read_into_the_named_figures():
    # as tools/significance.py prints a comparison (its print statements; figures made up)
    text = ("learned vs g_stalta_8: 120 events\n"
            "  detected            100 vs 61  (95% CP 0.75-0.89 vs 0.42-0.60)\n"
            "  chance-corrected    98.5 vs 57.0, difference +41.5 (95% bootstrap +32.0 to +50.9, P(diff <= 0) = 0.0000)\n"
            "  discordant          learned only 42, g_stalta_8 only 3, exact McNemar p = 1.1e-09\n")
    (m,) = gm.SIG.finditer(text)
    assert m.group(1) == "g_stalta_8"
    assert dict(zip(gm.KEYS, m.groups()[1:])) == {
        "events": "120", "detected_L": "100", "detected_B": "61", "corrected_L": "98.5", "corrected_B": "57.0",
        "difference": "+41.5", "ci_low": "+32.0", "ci_high": "+50.9", "p_diff_le0": "0.0000",
        "L_only": "42", "B_only": "3", "mcnemar_p": "1.1e-09"}
