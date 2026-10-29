"""tools/sweep_trigger.py: report parsing, catalogue copies, chance correction,
the CSV.

Expected values: the fixtures' own summaries and station tables (read off the
report text by hand), arithmetic by hand in the comments, and invariants.
"""
import calendar
import csv
import json
import math
import re

import pytest

import sweep_trigger as st
from conftest import fixture_text, run_tool, write_fake_binary, write_mseed_stub

FIELDS_EMPTY = ["alarms", "correct", "false", "catalogue", "detected", "missed", "alarm_delay_s"]


def epoch(dmy_hms):
    """'DD/MM/YYYY HH:MM:SS[.f]' -> epoch seconds, parsed by hand (not datetime)."""
    day, clock = dmy_hms.split(" ")
    d, m, y = (int(x) for x in day.split("/"))
    hh, mm, ss = clock.split(":")
    return calendar.timegm((y, m, d, int(hh), int(mm), 0, 0, 0, 0)) + float(ss)


# --- parse_report on real reports --------------------------------------------

# From each fixture's SUMMARY block and station table (by hand):
#   new learned / frozen learned: Alarms raised 10, correct 8, false 2, AFAD
#   events within 250 km 9, detected 8, missed 1, median 17.8 s; bands M<2 3/2,
#   M2-3 3/3, M3-4 2/2, M4-5 1/1 (no M5+ row); largest matched event MW 4.9,
#   alarm 18.0 s after origin, "by KULA, SIMA; later DKL, KCTX, TVSB, YAYO";
#   triggers 2+13+18+9+11+8 = 61 (the report's own "for 21 of 61 triggers").
#   STA/LTA: Alarms raised 33, correct 8, false 25, 9 events, detected 8,
#   missed 1, median 15.0 s, same bands; MW 4.9 matched by Event #21, alarm
#   17.9 s, "by SIMA, YAYO; later DKL"; triggers 33+30+31+24+25+26 = 169.
LEARNED = {"alarms": "10", "correct": "8", "false": "2", "catalogue": "9", "detected": "8",
           "missed": "1", "alarm_delay_s": "17.8", "M < 2": "2/3", "M 2-3": "3/3",
           "M 3-4": "2/2", "M 4-5": "1/1", "M 5+": "", "main_M": 4.9, "main_alarm_s": "18.0",
           "main_stations": "KULA SIMA + DKL KCTX TVSB YAYO", "triggers": 61,
           "dropped_warnings": 0}
STALTA = {**LEARNED, "alarms": "33", "false": "25", "alarm_delay_s": "15.0",
          "main_alarm_s": "17.9", "main_stations": "SIMA YAYO + DKL", "triggers": 169}


@pytest.mark.parametrize("name, want", [
    ("report_new_learned_demo_ko.txt", LEARNED),
    ("report_frozen_learned_demo_ko.txt", LEARNED),
    ("report_new_stalta3_demo_ko.txt", STALTA),
])
def test_parse_report_real(name, want):
    assert st.parse_report(fixture_text(name)) == want


def test_detected_before_s_line_is_not_the_detected_count():
    # The pre-S build adds "detected before S" right after "detected"; the
    # detected count must stay the first line's (7), not the second's (3).
    out = ("  AFAD events within 250 km               10\n"
           "    detected                              7   70%\n"
           "    detected before S                     3   30%\n"
           "    missed                                3   30%\n")
    row = st.parse_report(out)
    assert (row["catalogue"], row["detected"], row["missed"]) == ("10", "7", "3")


SYNTHETIC = """\
Event #1: detected earthquake (no magnitude), not located
  alarm      10:00:05.00 UTC by XA, XB
  AFAD       no matching catalogue event

Event #2: detected M4.4 earthquake at 39.000N 28.000E, origin 10:01:00.00 UTC
  alarm      10:01:12.50 UTC by PA, PB; later PC
  Pd magn.   M4.3 at the alarm, M4.4 ± 0.2 (M<4 10%, 4-5.5 90%, M>=5.5 0%) final from 3 stations (3 above noise)
  location   rms 0.10 s from 3 stations, first 50.0 s after the alarm (P and S picks)
  P loc.     39.010N 28.010E +-9 km from 3 stations, first 2.0 s after the alarm, epicentre 1.4 km off
  AFAD       ML 4.5 at 10:01:00.00 UTC: alarm 12.5 s after origin, Pd magnitude -0.1, epicentre 1.0 km off, origin +0.5 s
  pre-S      3 of 3 alarm detections declared before the predicted S

Event #3: detected earthquake (no magnitude), not located
  alarm      10:05:20.00 UTC by QA, QB
  AFAD       ML 2.0 at 10:05:00.00 UTC: alarm 20.0 s after origin
  pre-S      1 of 2 alarm detections declared before the predicted S
"""


def test_largest_event_its_own_alarm_line():
    # The largest matched event is Event #2 (ML 4.5 > ML 2.0), three lines
    # between its alarm and AFAD lines (the most the report has). The unmatched
    # Event #1 just before it must not lend it its stations (XA XB).
    row = st.parse_report(SYNTHETIC)
    assert (row["main_M"], row["main_alarm_s"], row["main_stations"]) == (4.5, "12.5", "PA PB + PC")


def test_empty_report():
    row = st.parse_report("")
    assert all(row[k] == "" for k in FIELDS_EMPTY + st.BANDS + ["main_M", "main_alarm_s",
                                                                 "main_stations"])
    assert (row["triggers"], row["dropped_warnings"]) == (0, 0)


# --- record_epoch, data_span ---------------------------------------------------

def test_data_span(tmp_path):
    # 2025-11-10 is day 314 of 2025 (31+28+31+30+31+30+31+31+30+31 = 304, +10).
    t = calendar.timegm((2025, 11, 10, 18, 4, 2, 0, 0, 0))
    a = write_mseed_stub(tmp_path / "A.mseed", t + 0.73, t + 1856.0)
    b = write_mseed_stub(tmp_path / "B.mseed", t + 3.0, t + 1857.01)
    assert st.record_epoch(a.read_bytes()[:512]) == pytest.approx(t + 0.73, abs=1e-6)
    lo, hi = st.data_span([str(a), str(b)])
    # The earliest first record and the latest last record over the files.
    assert lo == pytest.approx(t + 0.73, abs=1e-6)
    assert hi == pytest.approx(t + 1857.01, abs=1e-6)


# --- write_catalogue -------------------------------------------------------------

HEADER = "Date,Longitude,Latitude,Depth,Rms,Type,Magnitude,Location,EventID\n"


def test_write_catalogue_window_and_shift(tmp_path):
    span = (epoch("31/12/2025 12:00:00"), epoch("01/01/2026 06:00:00"))
    # The window is the span plus a day either side: 30/12/2025 12:00:00 to
    # 02/01/2026 06:00:00, both ends kept.
    src = tmp_path / "src.csv"
    src.write_text(HEADER
                   + "30/12/2025 11:59:59,28.1,39.2,7.0,0.3,ML,1.0,a,1\n"   # 1 s before: out
                   + "30/12/2025 12:00:00,28.1,39.2,7.0,0.3,ML,1.1,b,2\n"   # first second: in
                   + "31/12/2025 23:59:00,28.1,39.2,7.0,0.3,MW,4.9,c,3\n"
                   + "not a date,28.1,39.2,7.0,0.3,ML,1.2,d,4\n"
                   + "02/01/2026 06:00:00,28.1,39.2,7.0,0.3,ML,1.3,e,5\n"   # last second: in
                   + "02/01/2026 06:00:01,28.1,39.2,7.0,0.3,ML,1.4,f,6\n")  # 1 s after: out
    dst = tmp_path / "dst.csv"
    st.write_catalogue(src, dst, span, 127)
    # +127 s = +2 min 7 s; 23:59:00 + 127 s crosses into the new year.
    assert dst.read_text() == (HEADER
                               + "30/12/2025 12:02:07,28.1,39.2,7.0,0.3,ML,1.1,b,2\n"
                               + "01/01/2026 00:01:07,28.1,39.2,7.0,0.3,MW,4.9,c,3\n"
                               + "02/01/2026 06:02:07,28.1,39.2,7.0,0.3,ML,1.3,e,5\n")
    st.write_catalogue(src, dst, span, -391)
    # -391 s = -6 min 31 s: 12:00:00 -> 11:53:29, 23:59:00 -> 23:52:29,
    # 06:00:00 -> 05:53:29. The window applies to the unshifted times.
    assert dst.read_text() == (HEADER
                               + "30/12/2025 11:53:29,28.1,39.2,7.0,0.3,ML,1.1,b,2\n"
                               + "31/12/2025 23:52:29,28.1,39.2,7.0,0.3,MW,4.9,c,3\n"
                               + "02/01/2026 05:53:29,28.1,39.2,7.0,0.3,ML,1.3,e,5\n")


@pytest.mark.parametrize("shift", [127, 263, 391])
def test_write_catalogue_shift_round_trip(tmp_path, shift):
    # Invariant: shifting the real catalogue excerpt by +S then -S gives it
    # back line for line; and every +S line is the original's time + S (by the
    # hand parser) with the rest of the line unchanged.
    from conftest import REPO
    src = REPO / "tests/catalogs/demo_ko.csv"
    wide = (0.0, 4.0e9)
    plus, back = tmp_path / "plus.csv", tmp_path / "back.csv"
    st.write_catalogue(src, plus, wide, shift)
    st.write_catalogue(plus, back, wide, -shift)
    orig = src.read_text(encoding="utf-8").splitlines()
    assert back.read_text().splitlines() == orig
    shifted = plus.read_text().splitlines()
    assert len(shifted) == len(orig) == 803
    for o, s in zip(orig[1:], shifted[1:]):
        (od, orest), (sd, srest) = o.split(",", 1), s.split(",", 1)
        assert srest == orest
        assert epoch(sd) - epoch(od) == shift


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=(
    "write_catalogue parses '%d/%m/%Y %H:%M:%S' and silently skips a line whose "
    "time has fractional seconds (sweep_trigger.py:110-112), which ayzek reads since "
    "6942e49 ('GeoNet's and NOA's origin times carry them'); the shifted and unshifted "
    "copies then lack the event. Latent while ayzek_heldout writes whole seconds "
    "(heldout/catalog.py format_date)."))
def test_write_catalogue_keeps_fractional_seconds(tmp_path):
    t = epoch("10/11/2025 18:20:51")
    src = tmp_path / "src.csv"
    src.write_text(HEADER + "10/11/2025 18:20:51.25,28.145,39.239,9.6,0.6,MW,4.9,x,1\n")
    dst = tmp_path / "dst.csv"
    st.write_catalogue(src, dst, (t, t), 127)
    lines = dst.read_text().splitlines()
    assert len(lines) == 2
    assert abs(epoch(lines[1].split(",")[0]) - (t + 0.25 + 127)) < 1.0


# --- main: the grid, the shifted runs, the chance correction, the CSV -------------

FAKE_AYZEK = """
import calendar, json, os, sys
argv = sys.argv[1:]
with open(os.environ["FAKE_LOG"], "a") as f:
    f.write(json.dumps(argv) + "\\n")
sc = json.load(open(os.environ["FAKE_SCENARIO"]))
if sc.get("exit"):
    sys.exit(sc["exit"])
if sc.get("no_summary"):
    print("ayzek    something went wrong")
    sys.exit(0)
# The shift of this run: the catalogue's first event against the scenario's t0.
line = open(argv[argv.index("--catalog") + 1]).read().splitlines()[1]
day, clock = line.split(",")[0].split(" ")
d, m, y = (int(x) for x in day.split("/"))
hh, mm, ss = (int(x) for x in clock.split(":"))
shift = calendar.timegm((y, m, d, hh, mm, ss, 0, 0, 0)) - sc["t0"]
key = argv[argv.index("--stalta-on") + 1] if "--stalta-on" in argv else ""
c = sc["counts"][key][str(shift)]
print(f"  Alarms raised                        {c['alarms']:>4}")
print()
print("  Alarms compared with the AFAD catalogue")
print(f"    correct  (match an AFAD event)     {c['alarms'] - c['false']:>4}   -")
print(f"    false    (no AFAD event)           {c['false']:>4}   -")
print()
print(f"  AFAD events within 250 km             {c['catalogue']:>4}")
print(f"    detected                           {c['detected']:>4}   -")
print(f"    missed                             {c['catalogue'] - c['detected']:>4}   -")
"""


def counts(alarms, false, catalogue, detected):
    return {"alarms": alarms, "false": false, "catalogue": catalogue, "detected": detected}


@pytest.fixture
def sweep_env(tmp_path):
    """A fake ayzek, two stub data files spanning 18:04-18:34 on 2025-11-10,
    and a one-event catalogue (its event at t0 tells the fake the shift)."""
    t0 = calendar.timegm((2025, 11, 10, 18, 20, 51, 0, 0, 0))
    write_mseed_stub(tmp_path / "S1.mseed", t0 - 1011, t0 + 789)
    write_mseed_stub(tmp_path / "S2.mseed", t0 - 1000, t0 + 780)
    cat = tmp_path / "cat.csv"
    cat.write_text(HEADER + "10/11/2025 18:20:51,28.145,39.239,9.6,0.6,MW,4.9,x,1\n")
    fake = write_fake_binary(tmp_path / "ayzek", FAKE_AYZEK)
    log = tmp_path / "argv.jsonl"

    def run(scenario, *args):
        (tmp_path / "scenario.json").write_text(json.dumps({"t0": t0, **scenario}))
        p = run_tool("sweep_trigger", "--binary", fake, "--catalog", cat, "--jobs", 1,
                     "--out", tmp_path / "out.csv", *args,
                     env={"FAKE_LOG": str(log), "FAKE_SCENARIO": str(tmp_path / "scenario.json")})
        rows = (list(csv.DictReader(open(tmp_path / "out.csv")))
                if p.returncode == 0 else None)
        argvs = [json.loads(x) for x in log.read_text().splitlines()] if log.exists() else []
        return p, rows, argvs

    run.files = f"{tmp_path}/S1.mseed,{tmp_path}/S2.mseed"
    return run


def test_chance_correction(sweep_env):
    # Shifts 127,263 -> runs at 0, +127, -127, +263, -263.
    #
    # stalta-on 4: observed 6 of 10; shifted 1/10, 2/10, 0/9, 1/10.
    #   c = (0.1 + 0.2 + 0 + 0.1) / 4 = 0.1;
    #   corrected = (0.6 - 0.1) / (1 - 0.1) * 10 = 0.5 / 0.9 * 10 = 5.556 -> "5.6".
    # stalta-on 8: observed 4 of 10; shifted 0 of 0 events (+127), 1/10, 1/10,
    #   0/10; a shift without events counts 0 / max(1, 0) = 0, so
    #   c = (0 + 0.1 + 0.1 + 0) / 4 = 0.05;
    #   corrected = (0.4 - 0.05) / 0.95 * 10 = 3.684 -> "3.7".
    # stalta-on 10: every event detected at every shift: c = 1; (o - c)/(1 - c)
    #   is undefined, and a detector chance explains entirely shows no skill:
    #   the corrected count must be finite and not positive.
    sc = {"counts": {
        "4": {"0": counts(9, 2, 10, 6), "127": counts(5, 4, 10, 1), "-127": counts(5, 3, 10, 2),
              "263": counts(5, 5, 9, 0), "-263": counts(5, 4, 10, 1)},
        "8": {"0": counts(5, 1, 10, 4), "127": counts(3, 3, 0, 0), "-127": counts(3, 2, 10, 1),
              "263": counts(3, 2, 10, 1), "-263": counts(3, 3, 10, 0)},
        "10": {s: counts(12, 0, 10, 10) for s in ("0", "127", "-127", "263", "-263")}}}
    p, rows, argvs = sweep_env(sc, "--dataset", f"d1::{sweep_env.files}",
                               "--set", "stalta-on=4 8 10", "--shifts", "127,263")
    assert p.returncode == 0, p.stderr
    assert len(argvs) == 3 * 5
    by = {r["stalta-on"]: r for r in rows}
    assert set(by) == {"4", "8", "10"}
    assert (by["4"]["chance"], by["4"]["corrected"]) == ("0.100", "5.6")
    assert (by["8"]["chance"], by["8"]["corrected"]) == ("0.050", "3.7")
    assert by["10"]["chance"] == "1.000"
    assert math.isfinite(float(by["10"]["corrected"])) and float(by["10"]["corrected"]) <= 0
    # The other columns are the unshifted run's.
    assert {k: by["4"][k] for k in ("alarms", "correct", "false", "catalogue", "detected",
                                    "missed")} == {"alarms": "9", "correct": "7", "false": "2",
                                                   "catalogue": "10", "detected": "6",
                                                   "missed": "4"}
    assert list(rows[0])[:2] == ["dataset", "stalta-on"]
    assert all(r["dataset"] == "d1" for r in rows)


def test_without_shifts_no_chance_columns(sweep_env):
    sc = {"counts": {"4": {"0": counts(9, 2, 10, 6)}}}
    p, rows, argvs = sweep_env(sc, "--dataset", f"d1::{sweep_env.files}", "--set", "stalta-on=4")
    assert p.returncode == 0, p.stderr
    assert len(argvs) == 1
    assert (rows[0]["chance"], rows[0]["corrected"]) == ("", "")


def test_grid_and_options_passed_to_ayzek(sweep_env):
    # 'on' passes the bare flag, 'off' leaves it out, a value passes OPTION VALUE;
    # every run gets --speed 0 --no-color --no-magnitude, its catalogue copy,
    # --args, and --scores-in CACHE when the dataset names one; files last.
    sc = {"counts": {"4": {"0": counts(1, 0, 1, 1)}}}
    p, rows, argvs = sweep_env(sc, "--dataset", f"d1:CACHEDIR:{sweep_env.files}",
                               "--args", "--detector stalta --stalta-off 1.5",
                               "--set", "stalta-on=4", "--set", "stalta-3c=off on",
                               "--set", "sta=0.5")
    assert p.returncode == 0, p.stderr
    assert len(argvs) == 2 and len(rows) == 2
    for argv in argvs:
        pairs = {(a, b) for a, b in zip(argv, argv[1:])}
        assert {("--speed", "0"), ("--scores-in", "CACHEDIR"), ("--detector", "stalta"),
                ("--stalta-off", "1.5"), ("--stalta-on", "4"), ("--sta", "0.5")} <= pairs
        assert "--no-color" in argv and "--no-magnitude" in argv and "--catalog" in argv
        assert argv[-2:] == sweep_env.files.split(",")
    assert sorted("--stalta-3c" in argv for argv in argvs) == [False, True]
    assert sorted(r["stalta-3c"] for r in rows) == ["off", "on"]


@pytest.mark.parametrize("scenario", [{"exit": 3}, {"no_summary": True}])
def test_stops_when_ayzek_fails(sweep_env, scenario):
    # a24f610: a run that died or printed no catalogue summary stops the sweep.
    p, rows, _ = sweep_env(scenario, "--dataset", f"d1::{sweep_env.files}", "--set", "stalta-on=4")
    assert p.returncode != 0
    assert "without a report" in p.stderr
    assert rows is None
