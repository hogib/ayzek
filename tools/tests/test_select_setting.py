"""tools/select_setting.py: the three selection rules and their documented
tie-breaks, the sums over replays, the quiet tables, the boundaries.

Expected values: small tables built by hand, the rules as the docstring (and
docs/protocol-heldout.md section 6) states them, decimal arithmetic.
"""
import csv
import re
from decimal import Decimal

import pytest

import select_setting as ss
from conftest import run_tool

RULES = ("equal unmatched", "equal detections", "equal quiet alarms")


def write_sweeps(root, source, datasets, rows):
    """rows: (setting dict, {dataset: (detected, catalogue, corrected, false)}); a
    setting absent from a dataset's table is left out of that table."""
    d = root / source
    d.mkdir()
    for ds in datasets:
        keys = sorted({k for s, _ in rows for k in s})
        with open(d / f"sweep_{ds}.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, ["dataset", *keys, "alarms", "correct", "false", "catalogue",
                                    "detected", "missed", "alarm_delay_s", "chance", "corrected"])
            w.writeheader()
            for setting, per in rows:
                if ds not in per:
                    continue
                det, cat, cor, false = per[ds]
                w.writerow({"dataset": ds, **setting, "alarms": "", "correct": "", "false": false,
                            "catalogue": cat, "detected": det, "missed": cat - det,
                            "alarm_delay_s": "", "chance": "", "corrected": cor})


def write_quiet(path, keys, rows):
    """A tools/quiet_sweep.py table: rows (setting dict, alarms)."""
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, [*keys, "alarms", "hours", "per_day"])
        w.writeheader()
        for setting, alarms in rows:
            w.writerow({**setting, "alarms": alarms, "hours": "100.00", "per_day": ""})


def select(root, pattern, datasets, unmatched, corrected, *extra, env=None):
    p = run_tool("select_setting", "--sweeps", f"{root}/{pattern}", "--datasets", datasets,
                 "--unmatched", unmatched, "--corrected", corrected, *extra, env=env)
    assert p.returncode == 0, p.stderr
    out = {"stdout": p.stdout}
    for line in p.stdout.splitlines():
        for rule in RULES:
            if line.startswith(f"{rule}: "):
                out[rule] = line[len(rule) + 2:]
    return out


def setting_of(text):
    """'SRC  k=v k=v  detected ...' -> 'k=v k=v'."""
    return text.split("  ")[1]


# --- threshold(): the first rule's last tie-break --------------------------------------------

def test_threshold_of_each_detector():
    assert ss.threshold((("threshold", "0.8"), ("trigger-windows", "8"))) == (0.8,)
    assert ss.threshold((("stalta-band", "1,10"), ("stalta-on", "6"))) == (6.0,)
    assert ss.threshold((("fp-windows", "3,5,0.2"), ("fp-thresholds", "15,8"))) == (15.0, 8.0)
    # FilterPicker's pair compares in order: (10, 15) is above (10, 10) and
    # below (15, 8).
    assert ss.threshold((("fp-thresholds", "10,10"),)) < ss.threshold((("fp-thresholds", "10,15"),)) \
        < ss.threshold((("fp-thresholds", "15,8"),))


# --- rule 1: equal unmatched alarms ------------------------------------------------------------

T = [{"threshold": "0.5"}, {"threshold": "0.6"}, {"threshold": "0.7"}, {"threshold": "0.8"},
     {"threshold": "0.9"}]


def test_rule1_most_corrected_then_fewer_unmatched(tmp_path):
    # corrected, unmatched: 0.5 (10, 3); 0.6 (10, 4); 0.7 (10, 4); 0.8 (10, 4); 0.9 (12, 6).
    rows = list(zip(T, [{"d": (10, 20, "10.0", 3)}, {"d": (10, 20, "10.0", 4)},
                        {"d": (10, 20, "10.0", 4)}, {"d": (10, 20, "10.0", 4)},
                        {"d": (12, 20, "12.0", 6)}]))
    write_sweeps(tmp_path, "s", ["d"], rows)
    # <= 5 unmatched: 0.9 is out; the four at 10.0 tie, 0.5 has fewer unmatched.
    assert setting_of(select(tmp_path, "s", "d", 5, 99)["equal unmatched"]) == "threshold=0.5"
    # <= 6: 0.9 qualifies and has the most.
    assert setting_of(select(tmp_path, "s", "d", 6, 99)["equal unmatched"]) == "threshold=0.9"
    # <= 3: only 0.5. <= 2: none.
    assert setting_of(select(tmp_path, "s", "d", 3, 99)["equal unmatched"]) == "threshold=0.5"
    assert select(tmp_path, "s", "d", 2, 99)["equal unmatched"] == "none satisfies the rule"


def test_rule1_then_higher_threshold(tmp_path):
    # 0.6, 0.7, 0.8 all (10, 4): the highest threshold.
    rows = list(zip(T[1:4], [{"d": (10, 20, "10.0", 4)}] * 3))
    write_sweeps(tmp_path, "s", ["d"], rows)
    assert setting_of(select(tmp_path, "s", "d", 4, 99)["equal unmatched"]) == "threshold=0.8"


def test_rule1_higher_ratio_for_stalta(tmp_path):
    rows = [({"stalta-band": "1,10", "stalta-on": "4"}, {"d": (10, 20, "10.0", 4)}),
            ({"stalta-band": "1,10", "stalta-on": "10"}, {"d": (10, 20, "10.0", 4)}),
            ({"stalta-band": "2,10", "stalta-on": "6"}, {"d": (10, 20, "10.0", 4)})]
    write_sweeps(tmp_path, "s", ["d"], rows)
    assert setting_of(select(tmp_path, "s", "d", 4, 99)["equal unmatched"]) == "stalta-band=1,10 stalta-on=10"


def test_rule1_boundary_at_most(tmp_path):
    # "at most --unmatched": a setting with exactly that many qualifies.
    rows = [({"threshold": "0.8"}, {"d": (10, 20, "10.0", 26)}),
            ({"threshold": "0.7"}, {"d": (11, 20, "11.0", 27)})]
    write_sweeps(tmp_path, "s", ["d"], rows)
    assert setting_of(select(tmp_path, "s", "d", 26, 99)["equal unmatched"]) == "threshold=0.8"
    assert setting_of(select(tmp_path, "s", "d", 27, 99)["equal unmatched"]) == "threshold=0.7"


# --- rule 2: equal detections --------------------------------------------------------------------

def test_rule2_fewest_unmatched_then_more_detections(tmp_path):
    # corrected, unmatched: p (12, 6); q (13, 6); r (11, 5); s (9, 2).
    rows = [({"threshold": "0.9"}, {"d": (12, 20, "12.0", 6)}),
            ({"threshold": "0.8"}, {"d": (13, 20, "13.0", 6)}),
            ({"threshold": "0.7"}, {"d": (11, 20, "11.0", 5)}),
            ({"threshold": "0.6"}, {"d": (9, 20, "9.0", 2)})]
    write_sweeps(tmp_path, "s", ["d"], rows)
    # >= 11: p, q, r; r has the fewest unmatched.
    assert setting_of(select(tmp_path, "s", "d", 99, 11)["equal detections"]) == "threshold=0.7"
    # >= 12: p and q tie on unmatched; q has more detections.
    assert setting_of(select(tmp_path, "s", "d", 99, 12)["equal detections"]) == "threshold=0.8"
    # >= 13.0: q alone (at least); >= 13.1: none.
    assert setting_of(select(tmp_path, "s", "d", 99, 13.0)["equal detections"]) == "threshold=0.8"
    assert select(tmp_path, "s", "d", 99, 13.1)["equal detections"] == "none satisfies the rule"
    # >= 9: s, with 2 unmatched.
    assert setting_of(select(tmp_path, "s", "d", 99, 9)["equal detections"]) == "threshold=0.6"


# --- rule 3: equal quiet alarms --------------------------------------------------------------------

def test_rule3_most_corrected_then_fewer_quiet_then_fewer_unmatched(tmp_path):
    # corrected, unmatched, quiet: x (10, 7, 5); y (10, 9, 4); z (10, 8, 4);
    # w (11, 1, 9); v (20, 0, no quiet row).
    rows = [({"threshold": "0.5"}, {"d": (10, 20, "10.0", 7)}),
            ({"threshold": "0.6"}, {"d": (10, 20, "10.0", 9)}),
            ({"threshold": "0.7"}, {"d": (10, 20, "10.0", 8)}),
            ({"threshold": "0.8"}, {"d": (11, 20, "11.0", 1)}),
            ({"threshold": "0.9"}, {"d": (20, 20, "20.0", 0)})]
    write_sweeps(tmp_path, "s", ["d"], rows)
    q = tmp_path / "quiet.csv"
    write_quiet(q, ["threshold"], [({"threshold": "0.5"}, 5), ({"threshold": "0.6"}, 4),
                                   ({"threshold": "0.7"}, 4), ({"threshold": "0.8"}, 9)])
    # <= 8 quiet alarms: x, y, z tie at 10.0; y and z have fewer quiet alarms;
    # z fewer unmatched. v has no quiet row and is never a candidate.
    r = select(tmp_path, "s", "d", 99, 99, "--quiet-alarms", 8, "--quiet-table", f"s={q}")
    assert setting_of(r["equal quiet alarms"]) == "threshold=0.7"
    assert "  (quiet windows: 4 alarms, at most 8)" in r["stdout"]
    # <= 9: w qualifies with the most. <= 3: none.
    r = select(tmp_path, "s", "d", 99, 99, "--quiet-alarms", 9, "--quiet-table", f"s={q}")
    assert setting_of(r["equal quiet alarms"]) == "threshold=0.8"
    r = select(tmp_path, "s", "d", 99, 99, "--quiet-alarms", 3, "--quiet-table", f"s={q}")
    assert r["equal quiet alarms"] == "none satisfies the rule"
    # Without --quiet-alarms the rule is not applied.
    assert "equal quiet alarms" not in select(tmp_path, "s", "d", 99, 99)


def test_rule3_quiet_table_applies_to_its_source_only(tmp_path):
    rows = [({"threshold": "0.8"}, {"d": (10, 20, "10.0", 1)})]
    write_sweeps(tmp_path, "a", ["d"], rows)
    write_sweeps(tmp_path, "b", ["d"], [({"threshold": "0.8"}, {"d": (15, 20, "15.0", 1)})])
    q = tmp_path / "quiet.csv"
    write_quiet(q, ["threshold"], [({"threshold": "0.8"}, 2)])
    # Source b has the better setting but no quiet table: a's is selected.
    # (The glob names the two source directories: a path that is not a sweep
    # directory, such as quiet.csv here, stops the tool with "missing".)
    r = select(tmp_path, "[ab]", "d", 99, 99, "--quiet-alarms", 5, "--quiet-table", f"a={q}")
    assert r["equal quiet alarms"].startswith("a  threshold=0.8")


# --- sums over replays, sources, missing tables -----------------------------------------------------

def test_sums_over_datasets_and_drops_a_setting_missing_in_one(tmp_path):
    rows = [({"threshold": "0.8", "trigger-windows": "8"},
             {"d1": (5, 9, "4.5", 2), "d2": (30, 50, "25.5", 10)}),
            ({"threshold": "0.9", "trigger-windows": "8"}, {"d1": (9, 9, "9.0", 0)})]  # d2 missing
    write_sweeps(tmp_path, "s", ["d1", "d2"], rows)
    r = select(tmp_path, "s", "d1,d2", 99, 0)
    assert r["equal unmatched"] == ("s  threshold=0.8 trigger-windows=8  detected 35/59, "
                                    "corrected 30.0, unmatched 12")
    assert r["equal detections"] == r["equal unmatched"]
    assert "1 settings over 1 sources; best per source at <= 99 unmatched:" in r["stdout"]


def test_best_per_source_list(tmp_path):
    write_sweeps(tmp_path, "a", ["d"], [({"threshold": "0.8"}, {"d": (10, 20, "10.0", 3)}),
                                        ({"threshold": "0.7"}, {"d": (12, 20, "12.0", 8)})])
    write_sweeps(tmp_path, "b", ["d"], [({"threshold": "0.8"}, {"d": (11, 20, "11.0", 9)})])
    r = select(tmp_path, "*", "d", 5, 99)
    assert r["equal unmatched"].startswith("a  threshold=0.8")
    assert "3 settings over 2 sources; best per source at <= 5 unmatched:" in r["stdout"]
    assert "\n  a: a  threshold=0.8  detected 10/20, corrected 10.0, unmatched 3\n" in r["stdout"]
    assert "\n  b: none at <= 5 unmatched\n" in r["stdout"]


def test_missing_table_stops(tmp_path):
    write_sweeps(tmp_path, "s", ["d1"], [({"threshold": "0.8"}, {"d1": (1, 1, "1.0", 0)})])
    p = run_tool("select_setting", "--sweeps", f"{tmp_path}/s", "--datasets", "d1,d2",
                 "--unmatched", 1, "--corrected", 1)
    assert p.returncode != 0 and "missing" in p.stderr and "sweep_d2.csv" in p.stderr


def test_no_settings_stops(tmp_path):
    p = run_tool("select_setting", "--sweeps", f"{tmp_path}/nothing*", "--unmatched", 1,
                 "--corrected", 1)
    assert p.returncode != 0 and "no settings" in p.stderr


# --- boundaries: decimal sums, full ties ----------------------------------------------------------------

FLOAT_REASON = (
    "select_setting sums the tables' one-decimal `corrected` values as floats "
    "(select_setting.py:73): 0.7 + 0.1 = 0.7999999999999999, so a setting whose sum is "
    "exactly the --corrected target (rule 2, 'at least') can be excluded, and two settings "
    "with equal sums (rule 1 tie-break) can compare unequal. Decimal sums were exact for "
    "the 2026-10-29 tuning sweeps; about 10% of five-term sums of one-decimal values are "
    "not.")


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=FLOAT_REASON)
def test_rule2_at_least_holds_at_the_exact_decimal_sum(tmp_path):
    # A: corrected 0.7 + 0.1 = 0.8, unmatched 2; B: 0.4 + 0.4 = 0.8, unmatched 3.
    rows = [({"threshold": "0.7"}, {"d1": (1, 9, "0.7", 1), "d2": (1, 9, "0.1", 1)}),
            ({"threshold": "0.8"}, {"d1": (1, 9, "0.4", 3), "d2": (1, 9, "0.4", 0)})]
    write_sweeps(tmp_path, "s", ["d1", "d2"], rows)
    assert sum(Decimal(x) for x in ("0.7", "0.1")) == Decimal("0.8")
    # Both have at least 0.8; A has fewer unmatched alarms.
    assert setting_of(select(tmp_path, "s", "d1,d2", 99, 0.8)["equal detections"]) == "threshold=0.7"


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=FLOAT_REASON)
def test_rule1_tie_on_equal_decimal_sums(tmp_path):
    rows = [({"threshold": "0.7"}, {"d1": (1, 9, "0.7", 1), "d2": (1, 9, "0.1", 1)}),
            ({"threshold": "0.8"}, {"d1": (1, 9, "0.4", 3), "d2": (1, 9, "0.4", 0)})]
    write_sweeps(tmp_path, "s", ["d1", "d2"], rows)
    # Equal corrected sums: the tie goes to fewer unmatched alarms, A.
    assert setting_of(select(tmp_path, "s", "d1,d2", 99, 0)["equal unmatched"]) == "threshold=0.7"


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=(
    "Among settings equal in corrected detections, unmatched alarms and threshold (here "
    "four trigger-windows values), max() takes the first in the iteration order of the "
    "set built at select_setting.py:67, which varies with PYTHONHASHSEED: a full tie is "
    "resolved differently from one run to the next. The 2026-10-29 tuning sweeps had no "
    "such tie (12 hash seeds gave the same selections)."))
def test_full_tie_is_resolved_the_same_way_in_every_run(tmp_path):
    rows = [({"threshold": "0.8", "trigger-windows": w}, {"d": (8, 9, "7.0", 2)})
            for w in ("4", "6", "8", "10")]
    write_sweeps(tmp_path, "s", ["d"], rows)
    picks = {select(tmp_path, "s", "d", 5, 1, env={"PYTHONHASHSEED": str(seed)})["equal unmatched"]
             for seed in range(12)}
    assert len(picks) == 1, sorted(picks)


def test_output_format(tmp_path):
    write_sweeps(tmp_path, "src1", ["d"], [({"stalta-band": "1,10", "sta": "0.5", "lta": "10",
                                             "stalta-on": "8", "stalta-3c": "off"},
                                            {"d": (118, 197, "110.3", 47)})])
    r = select(tmp_path, "src1", "d", 47, 110.3)
    want = ("src1  stalta-band=1,10 sta=0.5 lta=10 stalta-on=8 stalta-3c=off  "
            "detected 118/197, corrected 110.3, unmatched 47")
    assert r["equal unmatched"] == want and r["equal detections"] == want
    assert re.search(r"^1 settings over 1 sources; best per source at <= 47 unmatched:$",
                     r["stdout"], re.M)
