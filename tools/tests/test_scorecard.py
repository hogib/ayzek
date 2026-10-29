"""tools/scorecard.py's report parsers: the EVENT block regex on reports with
and without the pre-S line, score_dataset, assessed, summarise.

Expected values: read off the fixture reports by hand; arithmetic by hand on
a small hand-built scorecard.
"""
import pytest

import scorecard as sc
from conftest import fixture_text

# The matched events of the demo_ko learned reports (new build and frozen
# build; the same alarms, the new one adds a pre-S line after each AFAD line):
# (type, M, origin, alarm delay, Pd error, epicentre km, stations, Pd first, Pd final),
# as the report prints them.
LEARNED_EVENTS = [
    ("ML", "1.4", "18:07:17.00", "16.5", None, "2.7", "KULA, SIMA; later KCTX, YAYO", None, None),
    ("ML", "1.9", "18:16:19.00", "17.5", None, "2.3", "KULA, SIMA; later KCTX, TVSB, YAYO", None, None),
    ("MW", "4.9", "18:20:51.00", "18.0", "-0.1", "3.7", "KULA, SIMA; later DKL, KCTX, TVSB, YAYO", "4.7", "4.8"),
    ("ML", "3.1", "18:23:21.00", "22.5", None, "158.2", "KCTX, KULA; later TVSB", None, None),
    ("ML", "3.6", "18:25:51.00", "18.0", None, "3.5", "KULA, SIMA; later DKL, KCTX, TVSB, YAYO", None, None),
    ("ML", "2.4", "18:29:11.00", "6.0", None, "5.2", "KULA, YAYO; later KCTX, TVSB", None, None),
    ("ML", "2.7", "18:30:02.00", "18.0", None, "3.7", "KULA, SIMA; later KCTX, TVSB, YAYO", None, None),
    ("ML", "2.4", "18:32:40.00", "12.5", None, None, "KCTX, SIMA; later KULA, YAYO", None, None),
]


def groups(m):
    return (m.group(6), m.group(7), m.group(8), m.group(9), m.group(10), m.group(11), m.group(3),
            m.group(4), m.group(5))


@pytest.mark.parametrize("name", ["report_new_learned_demo_ko.txt", "report_frozen_learned_demo_ko.txt"])
def test_event_regex_with_and_without_pre_s_line(name):
    ms = list(sc.EVENT.finditer(fixture_text(name)))
    assert [groups(m) for m in ms] == LEARNED_EVENTS
    # Event #3's headline magnitude and alarm time.
    assert (ms[2].group(1), ms[2].group(2)) == ("4.8", "18:21:09.00")
    # The two unmatched alarms (Events #8 and #10) have no AFAD match line.
    assert all(m.group(1) is None for m in ms if m is not ms[2])


def test_event_regex_on_the_stalta_report():
    ms = list(sc.EVENT.finditer(fixture_text("report_new_stalta3_demo_ko.txt")))
    # 8 matched events; the unmatched M4.8 alarm (#20) a few lines before the
    # matched MW 4.9 (#21) is not taken for it.
    assert len(ms) == 8
    assert [m.group(8) for m in ms] == ["18:07:17.00", "18:16:19.00", "18:20:51.00", "18:23:21.00",
                                        "18:25:51.00", "18:29:11.00", "18:30:02.00", "18:32:40.00"]
    mw = ms[2]
    assert groups(mw) == ("MW", "4.9", "18:20:51.00", "17.9", "-0.3", "4.3", "SIMA, YAYO; later DKL",
                          "4.6", "4.6")
    assert (mw.group(1), mw.group(2)) == ("4.6", "18:21:08.92")
    assert ms[3].group(11) is None and ms[5].group(11) is None and ms[7].group(11) is None


def test_event_regex_geometry_block():
    # The geometry locator's block: "not located (2 of 3 stations)" headline,
    # a distance line, no location line, the pre-S line after AFAD.
    out = ("Event #4: detected M4.0 earthquake, not located (2 of 3 stations)\n"
           "  alarm      10:00:12.00 UTC by AA, BB\n"
           "  Pd magn.   M3.9 at the alarm, M4.0 ± 0.3 (M<4 40%, 4-5.5 60%, M>=5.5 0%) final from 2 stations (2 above noise)\n"
           "  distance   AA 40-60 km, BB 70-90 km\n"
           "  AFAD       ML 4.2 at 10:00:00.00 UTC: alarm 12.0 s after origin, Pd magnitude -0.2\n"
           "  pre-S      2 of 2 alarm detections declared before the predicted S\n")
    (m,) = sc.EVENT.finditer(out)
    assert groups(m) == ("ML", "4.2", "10:00:00.00", "12.0", "-0.2", None, "AA, BB", "3.9", "4.0")
    assert m.group(1) == "4.0"


# --- score_dataset -------------------------------------------------------------------------------

def without(d, *keys):
    return {k: v for k, v in d.items() if k not in keys}


def test_score_dataset_new_learned():
    # From the report's summary: 10 alarms, 8 correct, 2 false, 9 events,
    # 8 detected, median 17.8 s, Pd error mean 0.12 bias -0.12 over 1 event,
    # alarm before the S wave at 47 of 48 arrivals, median warning 15.0 s.
    s = sc.score_dataset("demo_ko", fixture_text("report_new_learned_demo_ko.txt"))
    assert without(s, "stations", "events") == {
        "role": "public", "alarms": 10, "matched": 8, "unmatched": 2, "catalogue": 9,
        "detected": 8, "alarm_delay_median_s": 17.8,
        "bands": {"M < 2": "2/3", "M 2-3": "3/3", "M 3-4": "2/2", "M 4-5": "1/1"},
        "pd_magnitude_mae": 0.12, "pd_magnitude_bias": -0.12, "pd_magnitude_events": 1,
        "before_s": "47/48", "warning_median_s": 15.0, "assessed": {}, "dropped_warnings": 0}
    assert s["events"][2] == {"afad_type": "MW", "afad_M": 4.9, "origin": "18:20:51.00",
                              "alarm_delay_s": 18.0, "epicentre_km": 3.7,
                              "stations": "KULA, SIMA + DKL, KCTX, TVSB, YAYO",
                              "Pd_alarm": 4.7, "Pd_final": 4.8, "Pd_error": -0.1}
    assert s["events"][7] == {"afad_type": "ML", "afad_M": 2.4, "origin": "18:32:40.00",
                              "alarm_delay_s": 12.5, "epicentre_km": None,
                              "stations": "KCTX, SIMA + KULA, YAYO"}
    assert len(s["events"]) == 8


def test_score_dataset_frozen_and_stalta():
    f = sc.score_dataset("demo_ko", fixture_text("report_frozen_learned_demo_ko.txt"))
    n = sc.score_dataset("demo_ko", fixture_text("report_new_learned_demo_ko.txt"))
    # The two reports differ only by the pre-S lines (diff of the fixtures):
    # the same card, the warning summary read off line 220 of the frozen one
    # ("at 47 of 48 station arrivals, median warning 15.0 s").
    assert (f["before_s"], f["warning_median_s"]) == ("47/48", 15.0)
    assert f == n
    s = sc.score_dataset("demo_ko", fixture_text("report_new_stalta3_demo_ko.txt"))
    assert (s["alarms"], s["matched"], s["unmatched"], s["catalogue"], s["detected"]) == (33, 8, 25, 9, 8)
    assert (s["alarm_delay_median_s"], s["pd_magnitude_mae"], s["pd_magnitude_bias"]) == (15.0, 0.35, -0.35)
    assert s["events"][2]["Pd_error"] == -0.3 and s["events"][2]["epicentre_km"] == 4.3


NEGATIVE_BLOCK = """\
Event #1: detected M-0.2 earthquake at 39.000N 28.000E, origin 10:01:00.00 UTC
  alarm      10:00:58.50 UTC by PA, PB
  Pd magn.   M-0.3 at the alarm, M-0.2 ± 0.1 (M<4 100%, 4-5.5 0%, M>=5.5 0%) final from 2 stations (2 above noise)
  AFAD       ML -0.3 at 10:01:00.00 UTC: alarm -1.5 s after origin, Pd magnitude +0.1, epicentre 1.0 km off, origin +0.5 s
  pre-S      2 of 2 alarm detections declared before the predicted S
"""


def test_event_regex_negative_magnitudes():
    # Magnitudes printed with {:.1f} (src/pipeline/network.cpp:1090, 1134):
    # a GeoNet ML -0.3 and a Pd magnitude below zero keep their sign.
    m = sc.EVENT.search(NEGATIVE_BLOCK)
    assert m is not None
    assert groups(m) == ("ML", "-0.3", "10:01:00.00", "-1.5", "+0.1", "1.0", "PA, PB", "-0.3", "-0.2")
    assert m.group(1) == "-0.2"
    s = sc.score_dataset("demo_ko", NEGATIVE_BLOCK + "  AFAD events within 250 km                1\n")
    assert s["events"] == [{"afad_type": "ML", "afad_M": -0.3, "origin": "10:01:00.00", "alarm_delay_s": -1.5,
                            "epicentre_km": 1.0, "stations": "PA, PB", "Pd_alarm": -0.3, "Pd_final": -0.2,
                            "Pd_error": 0.1}]


def test_score_dataset_noise_rate():
    # noise_a1 is scored from 01:20 to 07:16: 5 h 56 min = 5.9333 h;
    # 3 unmatched alarms -> 24 * 3 / 5.9333 = 12.135 per day.
    out = ("  Alarms raised                           3\n"
           "    correct  (match an AFAD event)        0    0%\n"
           "    false    (no AFAD event)              3  100%\n"
           "  AFAD events within 250 km                0\n")
    s = sc.score_dataset("noise_a1", out)
    assert s["role"] == "noise"
    assert s["hours"] == pytest.approx(5 + 56 / 60)
    assert s["false_per_day"] == pytest.approx(72 / (5 + 56 / 60))
    assert (s["alarms"], s["unmatched"], s["catalogue"], s["detected"]) == (3, 3, 0, 0)
    assert s["alarm_delay_median_s"] is None and s["events"] == []


def test_assessed():
    out = ("  Alarm assessment (S-P test; docs/impl/16-alarm-assessment.md)\n"
           "    catalogued alarms                   3 earthquake, 1 possible, 4 misfire, 2 unclassified\n"
           "    unmatched alarms                    0 earthquake, 2 possible, 5 misfire, 1 unclassified\n")
    assert sc.assessed(out) == {
        "catalogued": {"earthquake": 3, "possible": 1, "misfire": 4, "unclassified": 2},
        "unmatched": {"earthquake": 0, "possible": 2, "misfire": 5, "unclassified": 1}}
    assert sc.assessed("no assessment here") == {}


# --- summarise ---------------------------------------------------------------------------------------

def test_summarise_by_hand():
    sets = {
        "demo": {"role": "tuning", "detected": 7, "catalogue": 9, "unmatched": 2,
                 "events": [{"afad_M": 4.9, "Pd_error": -0.1, "epicentre_km": 3.7},
                            {"afad_M": 1.4, "epicentre_km": 2.7}],
                 "assessed": {"catalogued": {"earthquake": 5, "possible": 1, "misfire": 1, "unclassified": 0},
                              "unmatched": {"earthquake": 1, "possible": 0, "misfire": 1, "unclassified": 0}}},
        "sindirgi61": {"role": "held out", "detected": 39, "catalogue": 66, "unmatched": 7,
                       "events": [{"afad_M": 6.1, "Pd_error": 0.3, "epicentre_km": 10.0}]},
        "demo_ko": {"role": "public", "detected": 8, "catalogue": 9, "unmatched": 2,
                    "events": [{"afad_M": 4.9, "Pd_error": -0.3, "epicentre_km": 4.3}]},
        "noise_a1": {"role": "noise", "unmatched": 3, "hours": 5 + 56 / 60, "events": []},
        "noise_b": {"role": "noise", "unmatched": 1, "hours": 12 + 8 / 60, "events": []},
        "quiet_ko": {"role": "public noise", "unmatched": 2, "hours": 5 + 50 / 60, "events": []},
    }
    h = sc.summarise(sets)
    assert h == {
        "detected (tuning)": "7/9", "unmatched alarms (tuning)": 2,
        "detected (held out)": "39/66", "unmatched alarms (held out)": 7,
        "detected (public)": "8/9", "unmatched alarms (public)": 2,
        # unmatched alarms judged real: 1 + 0 of 1 + 0 + 1 + 0; catalogued: 6 of 7
        "unmatched alarms assessed real": "1/2 (1 earthquake)",
        "catalogued alarms assessed real": "6/7",
        # noise: 4 alarms in 5.9333 + 12.1333 = 18.0667 h: 24*4/18.0667 = 5.3137
        "false alarms per day (noise)": 5.31, "hours (noise)": 18.1,
        # public noise: 2 in 5.8333 h: 8.2286
        "false alarms per day (public noise)": 8.23, "hours (public noise)": 5.8,
        # Pd errors -0.1, 0.3, -0.3 (all M >= 4): MAE 0.7/3 = 0.2333, bias -0.1/3
        "Pd magnitude MAE, all matched": 0.233, "Pd magnitude MAE, AFAD M>=4": 0.233,
        "Pd magnitude bias, AFAD M>=4": -0.033,
        # epicentres 3.7, 2.7, 10.0, 4.3: median (3.7 + 4.3)/2; AFAD replays
        # (tuning, held out) 3.7, 2.7, 10.0: median 3.7
        "epicentre error median km": 4.0, "epicentre error median km, AFAD replays": 3.7}


def test_summarise_without_events_or_noise():
    sets = {"demo": {"role": "tuning", "detected": 0, "catalogue": 9, "unmatched": 0, "events": []}}
    assert sc.summarise(sets) == {"detected (tuning)": "0/9", "unmatched alarms (tuning)": 0}
