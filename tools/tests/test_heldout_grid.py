"""tools/heldout_grid.py and tools/heldout_grid_summary.py (exploratory grids): the
settings a grid is run with, and the summary's reading of significance.py's output
and of the reports. No run: nothing here starts ayzek or run_heldout.

Expected values: the grid tables of the former drivers (logs/grid_run.py,
logs/instant_run.py) read by hand; the settings files those drivers wrote
(data/runs/heldout_grid, data/runs/heldout_instant; skipped where absent);
significance.py's print statements; the fixture report read by hand.
"""
import json

import pytest

import heldout_grid as hg
import heldout_grid_summary as hs
from conftest import REPO, fixture_text

FROZEN = {
    "_comment": "ignored",
    "learned": {"detector": "learned", "args": []},
    "stalta_eqfa": {"detector": "stalta", "args": ["--stalta-band", "1,10", "--stalta-on", "8", "--stalta-off", "1.5"]},
    "phasenet_causal_eqfa": {"detector": "learned", "model": "phasenet", "weights": "instance", "mode": "causal",
                             "args": ["--threshold", "0.8", "--trigger-windows", "1"]},
    "eqt_causal_eqfa": {"detector": "learned", "model": "eqtransformer", "args": ["--threshold", "0.7"]},
    "fp_eqfa": {"detector": "filterpicker", "args": ["--fp-windows", "6,10,0.4", "--fp-thresholds", "20,10"]},
    "transformer": {"detector": "transformer", "args": []},
}


def test_tradeoff_grid_has_the_five_families_in_order_with_filterpicker_last():
    s = hg.settings_for(FROZEN, "tradeoff")
    keys = list(s)
    assert keys[0] == "_comment"
    assert keys[1:] == (
        [f"g_learned_{v}" for v in ("0p4", "0p5", "0p6", "0p7", "0p75", "0p8", "0p85", "0p9")]
        + [f"g_stalta_{v}" for v in ("4", "5", "6", "7", "8", "10", "12", "15", "20")]
        + [f"g_phasenet_{v}" for v in ("0p3", "0p4", "0p5", "0p6", "0p7", "0p8", "0p9")]
        + [f"g_eqt_{v}" for v in ("0p3", "0p4", "0p5", "0p6", "0p7", "0p8", "0p9")]
        + [f"g_fp_{v}" for v in ("10", "12", "15", "20", "25", "30", "40")])
    assert len(keys) == 1 + 8 + 9 + 7 + 7 + 7


def test_a_grid_point_replaces_the_option_in_place_or_appends_it():
    s = hg.settings_for(FROZEN, "tradeoff")
    # the frozen learned entry carries no --threshold: appended
    assert s["g_learned_0p75"] == {"detector": "learned", "args": ["--threshold", "0.75"]}
    # replaced where it stands, the other options and keys untouched
    assert s["g_stalta_10"]["args"] == ["--stalta-band", "1,10", "--stalta-on", "10", "--stalta-off", "1.5"]
    assert s["g_phasenet_0p3"] == {"detector": "learned", "model": "phasenet", "weights": "instance",
                                   "mode": "causal", "args": ["--threshold", "0.3", "--trigger-windows", "1"]}
    # FilterPicker: the first threshold varies, the second stays 10; the label has the first only
    assert s["g_fp_25"]["args"] == ["--fp-windows", "6,10,0.4", "--fp-thresholds", "25,10"]
    # the frozen value is a point of its own curve, unchanged
    assert s["g_stalta_8"] == FROZEN["stalta_eqfa"]
    assert s["g_fp_20"] == FROZEN["fp_eqfa"]


def test_instant_grid_is_one_family_on_the_learned_detector():
    s = hg.settings_for(FROZEN, "instant")
    assert list(s) == ["_comment"] + [f"g_instant_{v}" for v in ("0p85", "0p88", "0p895", "0p9", "0p903", "0p905", "2")]
    assert s["g_instant_0p895"] == {"detector": "learned", "args": ["--instant-threshold", "0.895"]}
    assert s["g_instant_2"]["args"] == ["--instant-threshold", "2"]


def test_the_frozen_settings_are_not_modified():
    before = json.dumps(FROZEN)
    s = hg.settings_for(FROZEN, "tradeoff")
    s["g_stalta_4"]["args"].append("x")
    assert json.dumps(FROZEN) == before


def test_the_comment_says_exploratory_and_names_the_tool_and_grid():
    for grid in ("tradeoff", "instant"):
        (c,) = hg.settings_for(FROZEN, grid)["_comment"]
        assert "outside the protocol" in c and c.endswith(f"tools/heldout_grid.py --grid {grid}")


def test_default_output_directories_are_the_former_drivers():
    assert hg.DEFAULT_OUT == {"tradeoff": "data/runs/heldout_grid", "instant": "data/runs/heldout_instant"}
    assert hs.DEFAULT_OUT == hg.DEFAULT_OUT


@pytest.mark.parametrize("grid, driver", [("tradeoff", "logs/grid_run.py"), ("instant", "logs/instant_run.py")])
def test_settings_equal_the_file_the_grid_was_run_with(grid, driver):
    """The moved tool against the settings file its former driver wrote: the same
    text, entry for entry and in the same order, once the driver's name in the
    comment is replaced by the tool's (the one intended difference)."""
    ran = REPO / hg.DEFAULT_OUT[grid] / "grid_settings.json"
    if not ran.exists():
        pytest.skip(f"{ran} is not on this machine")
    frozen = json.loads((REPO / "tools/heldout_settings.json").read_text())
    got = hg.settings_for(frozen, grid)
    old = json.loads(ran.read_text())
    assert list(got) == list(old)
    for k in old:
        if k != "_comment":
            assert got[k] == old[k] and list(got[k]) == list(old[k]), k
    assert driver in ran.read_text()
    assert json.dumps(got, indent=1) == ran.read_text().replace(driver, f"tools/heldout_grid.py --grid {grid}")


# ---- the summary ---------------------------------------------------------------

# Two comparisons as tools/significance.py prints them (its print statements, figures made up).
SIGNIFICANCE = """\
group greece_now: 3 sets
learned vs g_learned_0p4: 120 events
  detected            100 vs 104  (95% CP 0.75-0.89 vs 0.79-0.92)
  chance-corrected    98.5 vs 101.2, difference -2.7 (95% bootstrap -6.1 to +0.4, P(diff <= 0) = 0.9512)
  discordant          learned only 1, g_learned_0p4 only 5, exact McNemar p = 0.22
learned vs g_stalta_8: 120 events
  detected            100 vs 61  (95% CP 0.75-0.89 vs 0.42-0.60)
  chance-corrected    98.5 vs 57.0, difference +41.5 (95% bootstrap +32.0 to +50.9, P(diff <= 0) = 0.0000)
  discordant          learned only 42, g_stalta_8 only 3, exact McNemar p = 1.1e-09
"""


def test_parse_significance_keeps_the_printed_figures_per_setting():
    sig = hs.parse_significance(SIGNIFICANCE)
    assert list(sig) == ["g_learned_0p4", "g_stalta_8"]
    assert sig["g_learned_0p4"] == ("g_learned_0p4", "120", "100", "104", "98.5", "101.2", "-2.7", "-6.1", "+0.4")
    assert sig["g_stalta_8"] == ("g_stalta_8", "120", "100", "61", "98.5", "57.0", "+41.5", "+32.0", "+50.9")


def test_report_counts_sum_alarms_and_unmatched_over_the_reports():
    # the fixture's summary: "Alarms raised 10", "false (no AFAD event) 2"
    rep = fixture_text("report_new_learned_demo_ko.txt")
    quiet = "SUMMARY\n  Alarms raised                           3\n"
    none = "SUMMARY\n  single-station detections, no alarm     1\n"
    assert hs.report_counts([rep, rep], [quiet, none, quiet]) == {"alarms": 20, "unmatched_alarms": 4, "quiet_alarms": 6}
    assert hs.report_counts([], []) == {"alarms": 0, "unmatched_alarms": 0, "quiet_alarms": 0}
