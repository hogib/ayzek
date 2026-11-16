"""tools/heldout_explore_common.py (helpers of the exploratory held-out reruns):
the readers and the rule that finds the station detections an alarm was raised
by. Synthetic recordings; the fixture report for the report reader.

Expected values: geometry by hand (one degree of a 6371 km sphere is 111.195
km), the recording format of src/pipeline/recording.hpp as the reader's
docstring gives it, the fixture report read by hand, calendar.timegm for its
epochs.
"""
import calendar

import pytest

import heldout_explore_common as hc
from conftest import FIXTURES

DEG_KM = 111.19492664455873        # 6371 * pi / 180


def test_distance_is_the_haversine_on_a_6371_km_sphere():
    assert hc.distance_km(0, 0, 1, 0) == pytest.approx(DEG_KM)
    assert hc.distance_km(0, 0, 0, 90) == pytest.approx(90 * DEG_KM)
    assert hc.distance_km(60, 0, 60, 1) == pytest.approx(DEG_KM / 2, rel=1e-4)    # cos 60 = 1/2
    assert hc.km_to_deg(DEG_KM) == pytest.approx(1.0)


def test_hms_cuts_to_hundredths_as_ayzek_prints():
    assert hc.hms(3661.239) == "01:01:01.23"
    assert hc.hms(86400.0 * 3 + 0.0) == "00:00:00.00"
    assert hc.hms(86399.995) == "23:59:59.99"       # cut, not rounded up into the next day


def test_epoch_near_takes_the_day_that_puts_the_time_nearest():
    day10 = 86400.0 * 10
    assert hc.epoch_near("00:01:40.00", day10 + 60) == day10 + 100
    assert hc.epoch_near("23:59:50.00", day10 + 100) == day10 - 10          # the evening before
    assert hc.epoch_near("00:00:05.50", day10 + 86390) == day10 + 86400 + 5.5   # the next morning


def recording(tmp_path, lines):
    f = tmp_path / "x.rec"
    f.write_text("\n".join(lines) + "\n")
    return f


def test_read_recording(tmp_path):
    rec = hc.read_recording(recording(tmp_path, [
        "# ayzek recording", "# from 1800.5",
        "S\tAAA\t39.5\t28.25", "S\tBBB\t-41.0\t174.0", "T\t1000.0\t9000.0",
        "D\tAAA\t1500.0\t1503.5", "D\tBBB\t1501.0\t1504.0", "D\tAAA\t2500.0\t2503.5"]))
    assert rec == {"stations": {"AAA": (39.5, 28.25), "BBB": (-41.0, 174.0)}, "t_first": 1000.0, "t_last": 9000.0,
                   "from": 1800.5, "detections": {"AAA": [(1500.0, 1503.5), (2500.0, 2503.5)], "BBB": [(1501.0, 1504.0)]}}


def test_read_recording_without_a_T_line_stops(tmp_path):
    with pytest.raises(SystemExit):
        hc.read_recording(recording(tmp_path, ["S\tAAA\t39.5\t28.25"]))


def test_read_report_on_the_fixture():
    # REPORT  6 stations, 2025-11-10 18:04:02.73 to 18:35:00.01 UTC; 10 alarms
    day = calendar.timegm((2025, 11, 10, 0, 0, 0))
    rec = {"t_first": day + 18 * 3600 + 4 * 60 + 2.73, "t_last": day + 18 * 3600 + 35 * 60 + 0.01}
    alarms = hc.read_report(FIXTURES / "report_new_learned_demo_ko.txt", rec)
    assert [a["number"] for a in alarms] == list(range(1, 11))
    a3 = alarms[2]          # Event #3: detected M4.8 earthquake at 39.206N 28.137E, origin 18:20:52.66 UTC
    assert abs(a3["time"] - (day + 18 * 3600 + 21 * 60 + 9.0)) < 1e-3      # alarm 18:21:09.00 UTC
    assert a3["hms"] == "18:21:09.00"
    assert a3["by"] == ["KULA", "SIMA"] and a3["later"] == ["DKL", "KCTX", "TVSB", "YAYO"]
    assert a3["location"][:2] == (39.206, 28.137)
    assert abs(a3["location"][2] - (day + 18 * 3600 + 20 * 60 + 52.66)) < 1e-3
    assert a3["magnitude"] == 4.8
    assert a3["credit"] == ("18:20:51.00", "MW", "4.9")
    a8 = alarms[7]          # Event #8: no magnitude, not located, no matching catalogue event
    assert (a8["location"], a8["magnitude"], a8["credit"]) == (None, None, None)
    assert a8["by"] == ["KCTX", "KULA"] and a8["later"] == ["SIMA"]
    assert alarms[0]["magnitude"] is None and alarms[0]["credit"] == ("18:07:17.00", "ML", "1.4")


def test_read_catalogue_takes_the_rows_ayzeks_reader_takes(tmp_path):
    f = tmp_path / "catalog.csv"
    f.write_text("﻿Date,Longitude,Latitude,Depth,Rms,Type,Magnitude,Location,EventID\n"
                 "10/11/2025 18:20:51,28.1,39.2,7.5,0.3,MW,4.9,SINDIRGI,E2\n"
                 "10/11/2025 18:07:17,28.2,39.1,10,0.2,ML,1.4,SINDIRGI,E1\n"
                 "not a date,28.2,39.1,10,0.2,ML,1.4,X,E3\n"
                 "short,row\n", encoding="utf-8")
    day = calendar.timegm((2025, 11, 10, 0, 0, 0))
    cat = hc.read_catalogue(f)
    assert [e["id"] for e in cat] == ["E1", "E2"]                 # time order, the unreadable rows left out
    assert cat[1] == {"time": day + 18 * 3600 + 20 * 60 + 51, "lon": 28.1, "lat": 39.2, "depth": 7.5, "mag": 4.9,
                      "type": "MW", "depth_text": "7.5", "id": "E2"}


def test_guard_out_refuses_the_release_outputs(tmp_path):
    assert hc.guard_out(tmp_path / "anywhere") == (tmp_path / "anywhere").resolve()
    for d in ("heldout", "heldout/greece_now", "heldout_scores", "heldout_quiet_scores/x"):
        with pytest.raises(SystemExit):
            hc.guard_out(hc.ROOT / "data/runs" / d)
    assert hc.guard_out(hc.ROOT / "data/runs/heldout_min3") == (hc.ROOT / "data/runs/heldout_min3").resolve()


def test_tree_digest_depends_on_names_and_contents(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "a/x.txt").write_text("1")
    (tmp_path / "a/y.txt").write_text("2")
    d1 = hc.tree_digest(tmp_path / "a")
    assert d1["files"] == 2
    (tmp_path / "a/y.txt").write_text("3")
    assert hc.tree_digest(tmp_path / "a")["sha256"] != d1["sha256"]
    (tmp_path / "a/y.txt").write_text("2")
    assert hc.tree_digest(tmp_path / "a") == d1


# ---- alarm_detections ----------------------------------------------------------
# Four stations on the equator. A, B and C are close (A-B 11 km, B-C 22 km, A-C 33 km:
# window starts compatible within 4.9, 6.7 and 8.6 s); X is 122 km from C (23.4 s).
STATIONS = {"A": (0.0, 0.0), "B": (0.0, 0.1), "C": (0.0, 0.3), "X": (0.0, 1.4)}


def rec_of(detections, stations=STATIONS):
    return {"stations": stations, "detections": detections, "t_first": 0.0, "t_last": 1000.0, "from": 0.0}


def alarm(number, time, by, later=()):
    return {"number": number, "time": time, "by": list(by), "later": list(later)}


def test_the_detection_declared_at_the_alarm_time_anchors_and_the_others_take_their_earliest_compatible():
    rec = rec_of({"A": [(50.0, 53.0), (100.0, 103.0), (102.0, 103.9)], "B": [(101.0, 104.0)]})
    dets, problems = hc.alarm_detections(alarm(1, 104.0, ["A", "B"]), rec)
    # B is declared at the alarm time; A's detection at 50 s is not compatible with it
    # (51 s apart, 4.9 s allowed), its next two are, and the earlier one is taken
    assert dets == {"A": (100.0, 103.0), "B": (101.0, 104.0)} and problems == []


def test_the_alarm_time_is_matched_within_the_hundredth_the_report_prints():
    rec = rec_of({"A": [(100.0, 103.0)], "B": [(101.0, 104.009)]})
    dets, problems = hc.alarm_detections(alarm(1, 104.0, ["A", "B"]), rec)
    assert dets["B"] == (101.0, 104.009) and problems == []


def test_detections_declared_after_the_alarm_are_not_its_own():
    rec = rec_of({"A": [(100.0, 104.5)], "B": [(101.0, 104.0)]})
    dets, problems = hc.alarm_detections(alarm(1, 104.0, ["A", "B"]), rec)
    assert dets == {"B": (101.0, 104.0)} and problems == ["A: no compatible detection"]


def test_no_anchor_is_a_named_problem_and_the_latest_compatible_are_taken():
    # nothing is declared at 110 s (the transformer's merge replaced the detection)
    rec = rec_of({"A": [(90.0, 93.0), (100.0, 103.0)], "B": [(101.0, 104.0)]})
    dets, problems = hc.alarm_detections(alarm(1, 110.0, ["A", "B"]), rec)
    assert problems == ["no detection declared at the alarm time"]
    assert dets == {"A": (100.0, 103.0), "B": (101.0, 104.0)}


def test_an_unknown_station_is_compatible_within_30_s():
    stations = {"A": (0.0, 0.0)}                    # B has no coordinates
    rec = rec_of({"A": [(100.0, 104.0)], "B": [(60.0, 63.0), (71.0, 74.0)]}, stations)
    dets, _ = hc.alarm_detections(alarm(1, 104.0, ["A", "B"]), rec)
    assert dets["B"] == (71.0, 74.0)                # 40 s apart is not, 29 s is


def two_alarms():
    """Alarm #1 by A and B, joined later by C (c1, declared after the alarm). Alarm #2 by
    X and C: C's own detection for it is c2; c1 is earlier and also compatible with X's
    (15 s apart, 23.4 s allowed)."""
    c1, c2 = (105.0, 108.0), (118.0, 121.0)
    rec = rec_of({"A": [(100.0, 103.0)], "B": [(101.0, 104.0)], "C": [c1, c2], "X": [(120.0, 123.0)]})
    return rec, alarm(1, 104.0, ["A", "B"], later=["C"]), alarm(2, 123.0, ["X", "C"]), c1, c2


def test_regression_the_later_stations_detections_are_reserved_for_their_alarm():
    """The bug fixed on 2026-11-14: a station that joined alarm #1 after it was raised
    (`later`) had its detection handed to alarm #2 (greece_now_seq020, transformer,
    alarm #11 given HORT's detection of alarm #10). With the alarms taken in number
    order and one `used` set, alarm #2 gets C's own detection."""
    rec, a1, a2, c1, c2 = two_alarms()
    used = set()
    d1, p1 = hc.alarm_detections(a1, rec, used=used)
    assert d1 == {"A": (100.0, 103.0), "B": (101.0, 104.0)} and p1 == []     # `later` is reserved, not returned
    assert ("C", c1[0]) in used
    d2, p2 = hc.alarm_detections(a2, rec, used=used)
    assert d2 == {"X": (120.0, 123.0), "C": c2} and p2 == []
    assert used == {("A", 100.0), ("B", 101.0), ("C", 105.0), ("X", 120.0), ("C", 118.0)}


def test_without_a_used_set_nothing_is_reserved():
    """The call without `used` (as tools/heldout_globalcat.py makes it) reserves
    nothing: the same alarm #2 is then given C's earliest compatible detection,
    the one of alarm #1."""
    rec, _, a2, c1, _ = two_alarms()
    d2, _ = hc.alarm_detections(a2, rec)
    assert d2["C"] == c1


def test_a_reserved_detection_is_still_used_when_the_station_has_no_other():
    rec, a1, a2, c1, _ = two_alarms()
    rec["detections"]["C"] = [c1]
    used = set()
    hc.alarm_detections(a1, rec, used=used)
    d2, p2 = hc.alarm_detections(a2, rec, used=used)
    assert d2["C"] == c1 and p2 == []               # unused candidates come first; they do not exclude


def test_only_a_compatible_detection_of_a_later_station_is_reserved_one_per_station():
    # C: the window at 80 s is declared after the alarm but 20 s from A's window (8.6 s
    # allowed); those at 105 s and 106 s are compatible, and the first declared is reserved
    rec = rec_of({"A": [(100.0, 103.0)], "B": [(101.0, 104.0)], "C": [(80.0, 105.0), (105.0, 108.0), (106.0, 109.0)]})
    used = set()
    hc.alarm_detections(alarm(1, 104.0, ["A", "B"], later=["C"]), rec, used=used)
    assert used == {("A", 100.0), ("B", 101.0), ("C", 105.0)}
