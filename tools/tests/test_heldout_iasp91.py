"""tools/heldout_iasp91.py (exploratory): the repeated catalogue match with the
half-space travel times (ayzek's own, the tool's control), on synthetic alarms
and events. No obspy, no release output read.

Expected values: Network::match and Network::score as the tool's docstring
states them (P within 6 s of window start + 3.5 s, or the location within 30 km
and 5 s with misfit |origin difference| / 100; smallest misfit credited; scored
within 250 km of the stations' centroid; first credited alarm counts; before S
when every crediting station's declaration precedes S), with P = hypot(D,
depth) / 6 and S = hypot(D, depth) / 3.5 by hand.
"""
import math

import pytest

import heldout_iasp91 as hi

# A and B 0.1 degree (11.12 km) apart on the equator
STATIONS = {"A": (0.0, 0.0), "B": (0.0, 0.1)}
HALF = hi.TravelTimes("halfspace", 40.0)


def ev(time, lat=0.0, lon=0.0, depth=0.0, mag=3.0, typ="ML"):
    return {"time": time, "lat": lat, "lon": lon, "depth": depth, "mag": mag, "type": typ}


def alarm(number, location=None):
    return {"number": number, "location": location}


def test_halfspace_travel_times():
    p, s = HALF(30.0, 40.0)                      # hypot = 50 km
    assert p == pytest.approx(50 / 6) and s == pytest.approx(50 / 3.5)
    h = math.hypot(30.0, 400.0)                  # the half-space model has no other rule for deep events
    assert HALF(30.0, 400.0) == pytest.approx((h / 6, h / 3.5))


def test_the_constants_are_the_matchers():
    assert (hi.VP, hi.VS, hi.TOLERANCE_S, hi.RADIUS_KM) == (6.0, 3.5, 6.0, 250.0)
    assert hi.SHIFTS == [0, 127, -127, 263, -263, 391, -391]


def test_an_alarm_is_credited_by_the_predicted_p_and_scored_before_s():
    # 35 km below A: P at A 1000 + 35/6 = 1005.83, S at A 1010.0; S at B (hypot 36.72 km) 1010.49
    events = [ev(1000.0, depth=35.0)]
    dets = [{"A": (1002.0, 1005.0), "B": (1003.0, 1009.0)}]       # A: 1005.5 against 1005.83
    first, scored, unmatched = hi.rematch([alarm(7)], dets, events, STATIONS, HALF)
    assert first == {0: 0} and unmatched == 0
    assert scored == {0: (1, 1, 7)}                               # both declared before S
    dets = [{"A": (1002.0, 1005.0), "B": (1003.0, 1012.0)}]       # B declared after S at B
    assert hi.rematch([alarm(7)], dets, events, STATIONS, HALF)[1] == {0: (1, 0, 7)}


def test_the_p_tolerance_is_six_seconds_inclusive():
    events = [ev(1000.0)]                                         # at A, surface: P at A = 1000.0
    at = [{"A": (1002.5, 1005.5)}]                                # 1006.0: 6.0 s off
    assert hi.rematch([alarm(1)], at, events, STATIONS, HALF)[1] == {0: (1, 0, 1)}
    past = [{"A": (1002.75, 1005.5)}]                             # 6.25 s off
    first, scored, unmatched = hi.rematch([alarm(1)], past, events, STATIONS, HALF)
    assert (first, scored, unmatched) == ({}, {0: (0, 0, "")}, 1)


def test_the_location_criterion_gives_the_smaller_misfit():
    # e0 fits by P with misfit 1.0 s; e1 by location: 2 s from the alarm's origin, misfit 0.02
    e0 = ev(1000.0)
    e1 = ev(990.0, lat=1.0, lon=0.0)                              # 111 km north: its P at A is 18.5 s after 990
    events = [e0, e1]
    dets = [{"A": (997.5, 1000.5)}]                               # 1001.0 against e0's 1000.0; e1's is 1008.5
    assert hi.rematch([alarm(1)], dets, events, STATIONS, HALF)[0] == {0: 0}
    located = alarm(1, location=(1.0, 0.1, 992.0))                # 11 km from e1, 2 s
    assert hi.rematch([located], dets, events, STATIONS, HALF)[0] == {0: 1}
    far = alarm(1, location=(1.0, 0.5, 992.0))                    # 56 km from e1: not within 30 km
    assert hi.rematch([far], dets, events, STATIONS, HALF)[0] == {0: 0}
    late = alarm(1, location=(1.0, 0.1, 996.0))                   # 6 s: not within 5 s
    assert hi.rematch([late], dets, events, STATIONS, HALF)[0] == {0: 0}


def test_the_first_credited_alarm_counts_and_the_second_is_unmatched():
    events = [ev(1000.0)]
    dets = [{"A": (996.5, 1001.0)}, {"A": (997.5, 1002.0)}]       # 1000.0 and 1001.0: both within 6 s
    first, scored, unmatched = hi.rematch([alarm(4), alarm(5)], dets, events, STATIONS, HALF)
    assert first == {0: 0} and scored == {0: (1, 0, 4)} and unmatched == 1      # declared after S (1000.0)


def test_an_event_beyond_250_km_of_the_centroid_is_credited_but_not_scored():
    # 2.3 degrees (255.7 km) east of the centroid at 0.05 E: lon 2.35; P at A after 261.3 / 6 = 43.55 s
    events = [ev(1000.0, lon=2.35)]
    dets = [{"A": (1040.0, 1043.0)}]                              # 1043.5 against 1043.55
    first, scored, unmatched = hi.rematch([alarm(1)], dets, events, STATIONS, HALF)
    assert scored == {} and first == {} and unmatched == 1
    events = [ev(1000.0, lon=2.25)]                               # 2.2 degrees, 244.6 km: scored; P 41.70 s
    dets = [{"A": (1038.0, 1041.0)}]
    assert hi.rematch([alarm(1)], dets, events, STATIONS, HALF)[1] == {0: (1, 1, 1)}


def test_a_detection_of_a_station_without_coordinates_is_ignored():
    events = [ev(1000.0)]
    dets = [{"ZZZ": (996.5, 999.0)}]
    assert hi.rematch([alarm(1)], dets, events, STATIONS, HALF) == ({}, {0: (0, 0, "")}, 1)


def test_event_key_is_the_unshifted_origin_type_and_magnitude():
    e = ev(3661.239 + 127, mag=4.96, typ="MW")
    assert hi.event_key(e, 127) == ("01:01:01.23", "MW", "5.0")
    assert hi.centroid(STATIONS) == (0.0, 0.05)
