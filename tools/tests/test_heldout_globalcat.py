"""tools/heldout_globalcat.py (exploratory): the rule that decides whether a
catalogue event explains an alarm, on synthetic detections and events, and the
reader of the FDSN text listings. No request is made; the arrival times are
prescribed by a stand-in for the TauP lookup, so obspy is not needed (one test
with the real lookup is skipped where obspy is not importable).

Expected values: the rule as the tool's docstring states it (rules 1 to 5),
arithmetic by hand: a detection's time is its window start + 3.5 s, a residual
is that time minus the predicted arrival.
"""
import calendar

import pytest

import heldout_globalcat as gc

TIERS = [(4.5, 180.0), (4.0, 30.0), (2.5, 10.0)]          # the default --tiers
STATIONS = {"A": (0.0, 0.0), "B": (0.0, 0.1), "C": (0.0, 0.2)}
# window starts 1000, 1001, 1002: detection times 1003.5, 1004.5, 1005.5
DETS = {"A": (1000.0, 1003.0), "B": (1001.0, 1004.0), "C": (1002.0, 1005.0)}
NEVER = 1e12


def event(eid, time, p=None, s=None, mag=5.0, lat=0.0, lon=50.0):
    """An event whose predicted arrivals are prescribed: p and s map a station to
    the epoch of its first P-type and S-type arrival (absent: no such arrival)."""
    return {"id": eid, "time": time, "lat": lat, "lon": lon, "mag": mag, "depth": 10.0, "p": p or {}, "s": s or {}}


def arrivals(ev, lat, lon):
    st = next(k for k, v in STATIONS.items() if v == (lat, lon))
    return ev["p"].get(st), ev["s"].get(st)


def explain(events, dets=DETS, tiers=TIERS, tol=6.0, need=2, shift=0):
    return gc.explain(dets, STATIONS, events, arrivals, tiers, tol, need, shift)


def test_an_alarm_is_explained_when_two_of_its_stations_fit_one_event_and_phase_type():
    ev = event("e", 900.0, p={"A": 1005.5, "B": 1001.5, "C": 1100.0})
    best, near = explain([ev])
    assert best[0] is ev and best[1] == "P"
    assert best[2] == {"A": -2.0, "B": 3.0}            # 1003.5 - 1005.5, 1004.5 - 1001.5; C is 94.5 s off
    assert best[3] == 2.5                              # mean absolute residual
    assert near[0] is ev and near[3] == 3.0            # the 2nd smallest |residual|


def test_one_station_within_the_tolerance_does_not_explain_and_the_nearest_miss_is_listed():
    ev = event("e", 900.0, p={"A": 1005.5, "B": 1014.5, "C": 1100.0})
    best, near = explain([ev])
    assert best is None
    assert near[0] is ev and near[1] == "P" and near[3] == 10.0      # B: 1004.5 - 1014.5
    assert near[2] == {"A": -2.0, "B": -10.0, "C": -94.5}


def test_the_tolerance_is_inclusive():
    at = event("at", 900.0, p={"A": 997.5, "B": 1010.5})           # +6.0 and -6.0
    assert explain([at])[0][2] == {"A": 6.0, "B": -6.0}
    past = event("past", 900.0, p={"A": 997.5, "B": 1010.75})      # +6.0 and -6.25
    assert explain([past])[0] is None
    assert explain([past], tol=6.25)[0] is not None


def test_p_and_s_fits_of_the_same_event_do_not_add_up():
    ev = event("e", 900.0, p={"A": 1003.5, "B": 1100.0}, s={"A": 1100.0, "B": 1004.5})
    best, _ = explain([ev], dets={"A": DETS["A"], "B": DETS["B"]})
    assert best is None
    ev_s = event("s", 900.0, s={"A": 1004.5, "B": 1004.5})
    best, _ = explain([ev_s])
    assert best[1] == "S" and best[2] == {"A": -1.0, "B": 0.0}


def test_the_event_with_the_smallest_mean_residual_is_listed():
    loose = event("loose", 900.0, p={"A": 1007.5, "B": 1008.5})      # -4, -4: mean 4
    tight = event("tight", 950.0, p={"A": 1004.5, "B": 1004.0})      # -1, +0.5: mean 0.75
    best, _ = explain([loose, tight])
    assert best[0] is tight and best[3] == 0.75
    best, _ = explain([tight, loose])
    assert best[0] is tight


def test_min_explained_three_needs_three_stations():
    ev = event("e", 900.0, p={"A": 1003.5, "B": 1004.5, "C": 1020.0})
    assert explain([ev], need=2)[0] is not None
    assert explain([ev], need=3)[0] is None
    ev3 = event("e3", 900.0, p={"A": 1003.5, "B": 1004.5, "C": 1006.5})
    assert explain([ev3], need=3)[0][2] == {"A": 0.0, "B": 0.0, "C": -1.0}


def test_rule_1_an_event_is_a_candidate_by_magnitude_and_distance_to_the_station():
    fit = {"A": 1003.5, "B": 1004.5}
    # 50 degrees away: M >= 4.5 only
    assert explain([event("m45", 900.0, p=fit, mag=4.5, lon=50.0)])[0] is not None
    assert explain([event("m44", 900.0, p=fit, mag=4.4, lon=50.0)]) == (None, None)
    # 20 degrees away: M >= 4.0
    assert explain([event("m40", 900.0, p=fit, mag=4.0, lon=20.0)])[0] is not None
    assert explain([event("m39", 900.0, p=fit, mag=3.9, lon=20.0)]) == (None, None)
    # 5 degrees away: M >= 2.5
    assert explain([event("m25", 900.0, p=fit, mag=2.5, lon=5.0)])[0] is not None
    assert explain([event("m24", 900.0, p=fit, mag=2.4, lon=5.0)]) == (None, None)


def test_rule_1_is_per_station():
    # M 2.6 at 9.95 degrees from C (0.2 E), 10.05 from B, 10.15 from A: a candidate at C alone
    ev = event("edge", 900.0, p={"A": 1003.5, "B": 1004.5, "C": 1005.5}, mag=2.6, lon=10.15)
    assert explain([ev]) == (None, None)
    ev = event("in", 900.0, p={"A": 1003.5, "B": 1004.5, "C": 1005.5}, mag=2.6, lon=10.05)
    best, _ = explain([ev])                             # now B and C
    assert best[2] == {"B": 0.0, "C": 0.0}


def test_an_event_after_the_detections_or_long_before_them_is_not_looked_at():
    fit = {"A": 1003.5, "B": 1004.5, "C": 1005.5}
    # origin later than the last detection time + tolerance (1005.5 + 6)
    assert explain([event("late", 1011.5, p=fit)])[0] is not None
    assert explain([event("later", 1011.6, p=fit)]) == (None, None)
    # origin more than 2400 s before the first detection time (1003.5)
    assert explain([event("old", 1003.5 - 2400.0, p=fit)])[0] is not None
    assert explain([event("older", 1003.4 - 2400.0, p=fit)]) == (None, None)


def test_the_chance_shift_moves_every_detection_time():
    ev = event("e", 900.0, p={"A": 1003.5, "B": 1004.5})
    assert explain([ev], shift=0)[0] is not None
    assert explain([ev], shift=600)[0] is None
    moved = event("m", 900.0, p={"A": 1603.5, "B": 1604.5})
    assert explain([moved], shift=0)[0] is None
    assert explain([moved], shift=600)[0][2] == {"A": 0.0, "B": 0.0}
    assert explain([event("b", 300.0, p={"A": 403.5, "B": 404.5})], shift=-600)[0] is not None


def test_stations_without_coordinates_do_not_count_and_no_detection_is_no_answer():
    ev = event("e", 900.0, p={"A": 1003.5})
    dets = {"A": DETS["A"], "ZZZ": (1001.0, 1004.0)}
    assert explain([ev], dets=dets)[0] is None
    assert explain([ev], dets={"ZZZ": (1001.0, 1004.0)}) == (None, None)
    assert explain([ev], dets={}) == (None, None)


def test_the_shifts_and_the_request_span_of_the_rule():
    assert gc.SHIFTS == [0, 600, -600, 1200, -1200, 1800, -1800]
    assert (gc.BEFORE_S, gc.AFTER_S) == (7200, 1800)
    assert (gc.GLOBAL_MIN_MAG, gc.REGIONAL_MIN_MAG, gc.REGIONAL_RADIUS_DEG) == (4.0, 2.5, 10.0)


def test_parse_fdsn_text():
    body = (b"#EventID | Time | Latitude | Longitude | Depth/km | Author | Catalog | Contributor | ContributorID | "
            b"MagType | Magnitude | MagAuthor | EventLocationName\n"
            b"us7000aaaa|2025-01-01T00:00:10.500|10.5|-20.25|33.0|us|us|us|x|mww|5.1|us|SOMEWHERE\n"
            b"us7000bbbb|2025-01-01T01:00:00Z|1.0|2.0||us|us|us|x|mb|4.2|us|NO DEPTH\n"
            b"us7000cccc|2025-01-01T02:00:00|1.0|2.0|5.0|us|us|us|x|||us|NO MAGNITUDE\n"
            b"\n")
    t0 = calendar.timegm((2025, 1, 1, 0, 0, 0))
    evs = gc.parse_fdsn_text(body)
    assert evs == [
        {"id": "us7000aaaa", "time": t0 + 10.5, "lat": 10.5, "lon": -20.25, "mag": 5.1, "magtype": "mww",
         "place": "SOMEWHERE", "depth": 33.0},
        {"id": "us7000bbbb", "time": t0 + 3600.0, "lat": 1.0, "lon": 2.0, "mag": 4.2, "magtype": "mb",
         "place": "NO DEPTH", "depth": None}]
    assert gc.parse_fdsn_text(b"") == []
    assert gc.parse_fdsn_text(b"us7000aaaa|2025-01-01T00:00:10|1|2|3|a|b|c|d|mb|4.0|e|no header\n") == []


def test_real_iasp91_arrivals_where_obspy_is_installed():
    pytest.importorskip("obspy.taup")
    arr = gc.Arrivals()
    ev = {"id": "e", "time": 1000.0, "lat": 0.0, "lon": 0.0, "mag": 6.0, "depth": None}     # no depth: 10 km
    p, s = arr(ev, 0.0, 30.0)
    # P takes about 6.2 min to 30 degrees and S about 11 min (travel-time tables); S after P
    assert 1000.0 + 340 < p < 1000.0 + 400
    assert 1000.0 + 620 < s < 1000.0 + 720
    assert arr(ev, 0.0, 30.0) == (p, s)
    assert arr({**ev, "depth": 10.0, "id": "same depth"}, 0.0, 30.0) == (p, s)
