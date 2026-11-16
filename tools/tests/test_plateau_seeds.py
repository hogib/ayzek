"""tools/plateau_seeds.py (exploratory) and tools/heldout_quiet_scores.py
(exploratory): the row of plateau.csv from a list of window probabilities, and
the rewrite of a release command for the quiet-window scores. No run.

Expected values: counting by hand on made-up lists; the percentile is the
sorted value at position round(q * n), counted from 1, as the tool computes it.
"""
import pytest

import heldout_quiet_scores as qs
import plateau_seeds as ps


def test_percentiles_and_counts_on_a_thousand_evenly_spaced_values():
    p = [i / 1000 for i in range(1, 1001)]                 # 0.001 ... 1.000
    assert ps.plateau_row("ensemble", "simav26", p) == {
        "run": "ensemble", "set": "simav26", "windows": 1000, "max": "1.000000",
        "p99": "0.990000",            # position 990
        "p99.9": "0.999000",          # position 999
        "p99.99": "1.000000",         # position round(999.9) = 1000
        "ge_0.895": 106,              # 0.895 ... 1.000
        "ge_0.9": 101, "ge_0.905": 96}


def test_the_input_need_not_be_sorted_and_the_thresholds_are_inclusive():
    row = ps.plateau_row("s42", "x", [0.905, 0.1, 0.9, 0.895, 0.8949])
    assert (row["windows"], row["max"]) == (5, "0.905000")
    assert (row["ge_0.895"], row["ge_0.9"], row["ge_0.905"]) == (3, 2, 1)
    # n = 5: positions round(4.95) = 5, round(4.995) = 5, round(4.9995) = 5: the maximum
    assert (row["p99"], row["p99.9"], row["p99.99"]) == ("0.905000",) * 3


def test_the_percentile_position_is_rounded_not_the_nearest_rank():
    """n = 1051, values 1 ... 1051: 0.99 n = 1040.49, taken at position 1040 (the
    nearest-rank definition, the smallest value with at least 99 % at or below
    it, is position 1041). Kept as the tool has computed it."""
    row = ps.plateau_row("s43", "x", [float(i) for i in range(1, 1052)])
    assert row["p99"] == "1040.000000"
    assert row["p99.9"] == "1050.000000"       # 1049.949 -> 1050
    assert row["p99.99"] == "1051.000000"      # 1050.89 -> 1051


def test_one_window():
    row = ps.plateau_row("s44", "x", [0.42])
    assert (row["windows"], row["max"], row["p99"], row["p99.99"], row["ge_0.9"]) == (1, "0.420000", "0.420000",
                                                                                        "0.420000", 0)


def test_the_runs_are_the_ensemble_and_its_three_members():
    assert ps.RUNS == {"ensemble": None, "s42": 42, "s43": 43, "s44": 44}


# ---- heldout_quiet_scores.plan -------------------------------------------------

LINE = ("build-linux/app/ayzek --speed 0 --no-color --catalog '../ayzek_heldout/data/gr now/catalog.csv' "
        "--record data/runs/heldout/greece_now/quiet/learned/gr_now_quiet07.rec --from 1800 a.mseed b.mseed")


def test_a_release_command_is_repeated_with_the_recording_moved_and_scores_added(tmp_path):
    name, argv, tmp, final, report = qs.plan(LINE, tmp_path)
    base = tmp_path / "greece_now"
    assert name == "gr_now_quiet07"
    assert argv == ["build-linux/app/ayzek", "--speed", "0", "--no-color", "--catalog",
                    "../ayzek_heldout/data/gr now/catalog.csv",
                    "--record", str(base / "rec/gr_now_quiet07.rec"),
                    "--scores", str(base / "scores/learned/gr_now_quiet07.tmp"),
                    "--from", "1800", "a.mseed", "b.mseed"]
    assert (tmp, final) == (base / "scores/learned/gr_now_quiet07.tmp", base / "scores/learned/gr_now_quiet07")
    assert report == base / "reports/gr_now_quiet07.txt"
    assert all((base / d).is_dir() for d in ("rec", "reports", "scores/learned"))
    assert not final.exists()


@pytest.mark.parametrize("line", [
    LINE + " --scores somewhere",                                           # already writes scores
    LINE.replace("quiet/learned", "quiet/stalta_eqfa"),                     # another detector's run
    LINE.replace("data/runs/heldout/", "data/runs/heldout_min3/"),          # not a release run
    LINE.replace("quiet/learned", "runs/learned"),                          # a sequence set
])
def test_an_unexpected_command_stops_the_script(tmp_path, line):
    with pytest.raises(SystemExit):
        qs.plan(line, tmp_path)


def test_a_window_marked_done_is_kept_without_running(tmp_path):
    job = qs.plan(LINE, tmp_path)
    job[3].mkdir()
    (job[3] / ".done").touch()
    assert qs.run(job) == ("gr_now_quiet07", 0.0, "kept")
    assert not job[4].exists()                  # no report written: nothing ran
