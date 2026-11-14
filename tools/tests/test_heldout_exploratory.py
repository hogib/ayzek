"""tools/heldout_exploratory.py: the exploratory figures read from a held-out
run's outputs.

Expected values: arithmetic by hand on small constructed tables and reports
(written in each test), great-circle distances along the equator and a
meridian (one degree = 6371 pi / 180 km), and counts read off the constructed
report text. The reports are written here in the layout of ayzek's (the lines
the tool reads: ALARM, the event blocks' "Pd magn.", the "False alarms" list,
the "Station processing" table; for item 9 the event blocks' "alarm", "AFAD",
"pre-S" and warning-table lines).
"""

import csv
import json
import math

import pytest
from conftest import run_tool

import heldout_exploratory as he

DEG_KM = 6371.0 * math.pi / 180


# --- paired latency --------------------------------------------------------------------


def latency_rows():
    """Alarm delays (learned, STA/LTA): e1 3.0, 2.0; e2 1.5, 1.5 (a tie); e3
    1.0, 1.5; e4 4.0, 2.0; e5 learned only; e6 STA/LTA only; and in another set
    an event with e1's origin, 9.0 against 9.0 for a third detector only."""
    rows = []
    for i, (a, b) in enumerate(
        [(3.0, 2.0), (1.5, 1.5), (1.0, 1.5), (4.0, 2.0), (5.0, None), (None, 7.0)], 1
    ):
        for det, delay in (("learned", a), ("stalta_eqfa", b)):
            if delay is not None:
                rows.append(
                    {
                        "detector": det,
                        "set": "s1",
                        "origin": f"00:0{i}:00.00",
                        "type": "ML",
                        "magnitude": "2.0",
                        "alarm_delay_s": str(delay),
                    }
                )
    rows.append(
        {
            "detector": "learned",
            "set": "s2",
            "origin": "00:01:00.00",
            "type": "ML",
            "magnitude": "2.0",
            "alarm_delay_s": "9.0",
        }
    )
    rows.append(
        {
            "detector": "phasenet_causal_eqfa",
            "set": "s2",
            "origin": "00:01:00.00",
            "type": "ML",
            "magnitude": "2.0",
            "alarm_delay_s": "9.0",
        }
    )
    return rows


def test_quantile_linear_interpolation():
    # Positions q (n - 1) in 1, 2, 3, 4: 0.75, 1.5, 2.25.
    v = [4, 1, 3, 2]
    assert [he.quantile(v, q) for q in (0.25, 0.5, 0.75)] == pytest.approx(
        [1.75, 2.5, 3.25]
    )
    assert he.quantile([7.0], 0.25) == 7.0
    assert (he.quantile(v, 0), he.quantile(v, 1)) == (1, 4)


def test_paired_latency_with_a_tie():
    # Both detected e1-e4: differences +1, 0, -0.5, +2; sorted -0.5, 0, 1, 2:
    # median 0.5, quartiles -0.5 + 0.75 * 0.5 = -0.125 and 1 + 0.25 * 1 = 1.25.
    # STA/LTA first on e1 and e4, learned first on e3, e2 a tie. The event of
    # s2 has e1's origin but is another set's, and STA/LTA has no row for it.
    r = he.paired_latency(latency_rows(), "learned", "stalta_eqfa")
    assert r == {
        "events": 4,
        "learned_first": 1,
        "stalta_eqfa_first": 2,
        "ties": 1,
        "median": pytest.approx(0.5),
        "q1": pytest.approx(-0.125),
        "q3": pytest.approx(1.25),
    }


def test_paired_latency_without_common_events():
    r = he.paired_latency(latency_rows(), "stalta_eqfa", "phasenet_causal_eqfa")
    assert r == {
        "events": 0,
        "stalta_eqfa_first": 0,
        "phasenet_causal_eqfa_first": 0,
        "ties": 0,
    }


# --- sequences and chance levels -------------------------------------------------------


def write_outcomes(d, label, s, detected, shifted, magnitudes=None):
    """OUTCOMES/LABEL_SET.csv: len(detected) events (1/0 each), and per shift
    the same events with `shifted[shift]` giving those detected."""
    d.mkdir(parents=True, exist_ok=True)
    mags = magnitudes or ["2.0"] * len(detected)
    with open(d / f"{label}_{s}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "dataset",
                "detector",
                "shift",
                "origin",
                "type",
                "magnitude",
                "detected",
                "detected_pre_s",
                "depth_km",
            ]
        )
        for shift, hit in [(0, detected), *shifted.items()]:
            for i, m in enumerate(mags):
                w.writerow(
                    [
                        s,
                        label,
                        shift,
                        f"00:{i:02d}:05.00",
                        "ML",
                        m,
                        hit[i] if shift == 0 else int(i in hit),
                        0,
                        "10",
                    ]
                )


def three_sequences(d):
    """s1: 4 events, learned 3, one of its 8 shifted events detected; STA/LTA 2.
    s2: 2 events, both detect one (not the same one). s3: 2 events, learned 0,
    STA/LTA 1."""
    none = {127: set(), -127: set()}
    write_outcomes(d, "learned", "s1", [1, 1, 1, 0], {127: {0}, -127: set()})
    write_outcomes(d, "stalta_eqfa", "s1", [1, 1, 0, 0], none)
    write_outcomes(d, "learned", "s2", [1, 0], none)
    write_outcomes(d, "stalta_eqfa", "s2", [0, 1], none)
    write_outcomes(d, "learned", "s3", [0, 0], none)
    write_outcomes(d, "stalta_eqfa", "s3", [1, 0], none)
    return he.load_outcomes(d, ["learned", "stalta_eqfa"], ["s1", "s2", "s3"])


def test_sequence_row_and_chance_level_by_hand(tmp_path):
    # s1, learned: o = 3/4, c = (1/4 + 0/4) / 2 = 1/8, corrected
    # 4 (3/4 - 1/8) / (7/8) = 20/7. STA/LTA: c = 0, corrected = detected = 2.
    det = three_sequences(tmp_path)
    r = he.sequence_row(det["learned"], "s1")
    assert (r["scored"], r["detected"], r["chance"]) == (4, 3, 0.125)
    assert r["corrected"] == pytest.approx(20 / 7)
    assert he.sequence_row(det["stalta_eqfa"], "s1") == {
        "scored": 4,
        "detected": 2,
        "chance": 0.0,
        "corrected": 2.0,
    }


def test_more_same_fewer_and_network_chance(tmp_path):
    # Detected, learned against STA/LTA: s1 3 > 2, s2 1 = 1 (different events:
    # the count is of detections, not of pairs), s3 0 < 1.
    # Network, learned: n = 8, raw 4, corrected 20/7 + 1 + 0 = 27/7, effective
    # c = (4 - 27/7) / (8 - 27/7) = (1/7) / (29/7) = 1/29. STA/LTA: c = 0.
    det = three_sequences(tmp_path)
    table = {
        s: {label: he.sequence_row(det[label], s) for label in det}
        for s in ("s1", "s2", "s3")
    }
    assert he.more_same_fewer(table, "learned", "stalta_eqfa") == {
        "more": 1,
        "same": 1,
        "fewer": 1,
    }
    assert he.more_same_fewer(table, "stalta_eqfa", "learned") == {
        "more": 1,
        "same": 1,
        "fewer": 1,
    }
    n = he.network_chance(table, "learned")
    assert (n["scored"], n["detected"]) == (8, 4)
    assert n["corrected"] == pytest.approx(27 / 7)
    assert n["effective_chance"] == pytest.approx(1 / 29)
    assert he.network_chance(table, "stalta_eqfa") == {
        "scored": 8,
        "detected": 4,
        "corrected": 4.0,
        "effective_chance": 0.0,
    }


# --- magnitude bands ----------------------------------------------------------------------


def test_magnitude_band_edges(tmp_path):
    # The lower edge belongs to the band: 2.0 is 2-3, 3.0 is 3-4, 4.0 is >=4;
    # 1.9, 2.9 and 3.9 stay below; a negative magnitude is <2.
    mags = ["1.9", "2.0", "2.9", "3.0", "3.9", "4.0", "-0.3", "5.6"]
    none = {127: set(), -127: set()}
    write_outcomes(tmp_path, "learned", "s1", [1, 1, 1, 1, 1, 1, 0, 1], none, mags)
    write_outcomes(tmp_path, "stalta_eqfa", "s1", [0, 1, 0, 1, 0, 1, 0, 1], none, mags)
    det = he.load_outcomes(tmp_path, ["learned", "stalta_eqfa"], ["s1"])
    assert he.magnitude_bands(det, ["learned", "stalta_eqfa"]) == [
        {"band": "<2", "scored": 2, "learned": 1, "stalta_eqfa": 0},
        {"band": "2-3", "scored": 2, "learned": 2, "stalta_eqfa": 1},
        {"band": "3-4", "scored": 2, "learned": 2, "stalta_eqfa": 1},
        {"band": ">=4", "scored": 2, "learned": 2, "stalta_eqfa": 2},
    ]


# --- closest pair --------------------------------------------------------------------------

COORDS = {"AAA": (0.0, 0.0), "BBB": (0.0, 0.3), "CCC": (1.0, 0.0), "DDD": (0.0, 0.5)}


def test_haversine_along_the_equator_and_a_meridian():
    assert he.haversine_km(0, 0, 0, 1) == pytest.approx(DEG_KM)
    assert he.haversine_km(10, 25, 12, 25) == pytest.approx(2 * DEG_KM)
    assert he.haversine_km(-41, 172, -41, 172) == 0


def test_closest_pair():
    # On the equator AAA 0.0, BBB 0.3, DDD 0.5 degrees east, CCC one degree
    # north of AAA: BBB-DDD, 0.2 degrees, is the nearest pair (AAA-BBB 0.3).
    pair, km = he.closest_pair(COORDS)
    assert pair == ("BBB", "DDD") and km == pytest.approx(0.2 * DEG_KM)
    # Without DDD: AAA-BBB at 0.3 degrees; the pair in alphabetical order.
    pair, km = he.closest_pair({k: COORDS[k] for k in ("CCC", "BBB", "AAA")})
    assert pair == ("AAA", "BBB") and km == pytest.approx(0.3 * DEG_KM)
    assert he.closest_pair({"AAA": (0, 0)}) is None


# --- reports -----------------------------------------------------------------------------------


def event_block(number, pd=None, afad=None):
    head = (
        f"detected M{pd[1]} earthquake" if pd else "detected earthquake (no magnitude)"
    )
    lines = [
        f"Event #{number}: {head}, not located",
        "  alarm      00:00:00.00 UTC by AAA, BBB",
    ]
    if pd:
        when = "at the alarm" if pd[2] is None else f"{pd[2]} s after the alarm"
        lines.append(
            f"  Pd magn.   M{pd[0]} {when}, M{pd[1]} ± 0.2 (M<4 50%, 4-5.5 50%, M>=5.5 0%) "
            "final from 3 stations (2 above noise)"
        )
    lines.append(
        f"  AFAD       {afad}" if afad else "  AFAD       no matching catalogue event"
    )
    return "\n".join(lines) + "\n"


def make_report(alarms, unmatched, stations=("AAA", "BBB", "CCC")):
    """alarms: (number, time, stations, AFAD note or None, (Pd first, final,
    seconds after the alarm or None) or None); unmatched: the numbers listed
    under "False alarms"."""
    out = [f"ayzek    {len(stations)} stations, scalar backend, replay"]
    for n, t, by, note, _ in alarms:
        out.append(
            f"ALARM    #{n:<3d} {t}  earthquake detected by {by}"
            + (f"  ({note})" if note else "")
        )
        out.append(
            f"PDMAG    #{n:<3d} {t}  Pd magnitude M9.9 ± 0.2 (M<4 0%, 4-5.5 0%, M>=5.5 100%) from 3 stations"
        )
    out.append("-" * 78 + "\n")
    out += [
        event_block(
            n, pd, note and "ML 2.0 at 00:00:00.00 UTC: alarm 10.0 s after origin"
        )
        for n, _, _, note, pd in alarms
    ]
    out.append("Summary\n  Alarms raised" + " " * 20 + str(len(alarms)))
    out.append(
        "  Missed AFAD events\n    00:59:00.00 UTC  ML 1.0    19 km from the network centre\n"
    )
    if unmatched:
        out.append(
            "  False alarms (no AFAD event; small uncatalogued earthquakes are possible)"
        )
        times = {n: t for n, t, *_ in alarms}
        out += [
            f"    {times[n]} UTC  #{n}  no magnitude, not located" for n in unmatched
        ]
    out.append("-" * 78)
    out.append(
        "Station processing (scalar backend)\n"
        "         windows  with gap triggers  picks Pd values    ms/window (p99)   ms/pick"
    )
    out += [
        f"  {s:5s}   21588         0      220    219      1540     94.866 (201.063)       207"
        for s in stations
    ]
    out.append("  Pd: 5130 values sent")
    return "\n".join(out) + "\n"


AFAD_NOTE = "AFAD ML 2.0 at 00:00:00.00, 10.0 s ago"
# Alarm 2 is matched; 3 carries the note on its ALARM line but the report's
# final matching lists it as unmatched; 10 has a two-digit number and three stations.
SEQUENCE_ALARMS = [
    (1, "00:01:00.00", "AAA, BBB", None, None),
    (2, "00:02:00.00", "AAA, CCC", AFAD_NOTE, ("4.1", "4.4", None)),
    (3, "00:03:00.00", "BBB, CCC", AFAD_NOTE, ("4.6", "4.0", "1.5")),
    (4, "00:04:00.50", "BBB, AAA", None, ("4.2", "3.9", "2.0")),
    (10, "00:10:00.00", "AAA, BBB, CCC", None, None),
]
SEQUENCE_REPORT = make_report(SEQUENCE_ALARMS, [1, 3, 4, 10])
QUIET_REPORT = make_report(
    [
        (1, "05:00:00.00", "CCC, AAA", None, ("-0.3", "-0.1", None)),
        (2, "05:30:00.00", "BBB, AAA", None, ("3.0", "5.2", "0.5")),
    ],
    [1, 2],
)


def test_parse_report():
    r = he.parse_report(SEQUENCE_REPORT)
    assert sorted(r["alarms"]) == [1, 2, 3, 4, 10]
    assert r["alarms"][3] == {
        "time": "00:03:00.00",
        "pd_final": 4.0,
        "stations": ["BBB", "CCC"],
    }
    assert r["alarms"][10]["stations"] == ["AAA", "BBB", "CCC"]
    # The final value of the "Pd magn." line, not its first nor a PDMAG line's 9.9.
    assert [r["alarms"][n]["pd_final"] for n in (1, 2, 3, 4, 10)] == [
        None,
        4.4,
        4.0,
        3.9,
        None,
    ]
    assert r["unmatched"] == [1, 3, 4, 10]  # not the missed event's line above the list
    assert r["stations"] == ["AAA", "BBB", "CCC"]
    q = he.parse_report(QUIET_REPORT)
    assert [q["alarms"][n]["pd_final"] for n in (1, 2)] == [-0.1, 5.2]


def test_report_without_unmatched_alarms():
    r = he.parse_report(make_report(SEQUENCE_ALARMS[1:2], []))
    assert r["unmatched"] == [] and list(r["alarms"]) == [2]


def test_large_estimates_among_unmatched_alarms():
    # Unmatched 1, 3, 4, 10: estimates for 3 (final 4.0) and 4 (final 3.9,
    # first 4.2): two with an estimate, one at or above 4.0 by its final value.
    r = he.parse_report(SEQUENCE_REPORT)
    m = he.alarm_magnitudes([("s1", r["alarms"][n]) for n in r["unmatched"]], 4.0)
    assert (m["alarms"], m["with_pd"], m["large"]) == (4, 2, 1)
    assert m["large_cases"] == [{"set": "s1", "time": "00:03:00.00", "pd_final": 4.0}]
    assert (
        he.alarm_magnitudes([("s1", r["alarms"][n]) for n in r["unmatched"]], 3.9)[
            "large"
        ]
        == 2
    )
    q = he.parse_report(QUIET_REPORT)
    m = he.alarm_magnitudes([("q1", a) for a in q["alarms"].values()], 4.0)
    assert (m["alarms"], m["with_pd"], m["large"], m["large_cases"][0]["pd_final"]) == (
        2,
        2,
        1,
        5.2,
    )


def test_alarms_declared_by_the_closest_pair():
    # Closest pair AAA-BBB. Sequence: unmatched 1 (AAA, BBB) and 4 (BBB, AAA:
    # the order does not matter) are the pair's, 3 (BBB, CCC) and 10 (three
    # stations) are not; the one matched alarm, 2, is by AAA, CCC. Quiet: both
    # alarms count as quiet whatever the list says; one is by BBB, AAA.
    pairs = {"s1": (("AAA", "BBB"), 33.4), "q1": (("AAA", "BBB"), 33.4), "q2": None}
    reports = [
        ("s1", "sequence", he.parse_report(SEQUENCE_REPORT)),
        ("q1", "quiet", he.parse_report(QUIET_REPORT)),
        ("q2", "quiet", he.parse_report(QUIET_REPORT)),
    ]
    assert he.pair_shares(reports, pairs) == {
        "matched": {"alarms": 1, "by_closest_pair": 0},
        "unmatched": {"alarms": 4, "by_closest_pair": 2},
        "quiet": {"alarms": 4, "by_closest_pair": 1},
    }


# --- quiet windows by station set ----------------------------------------------------------------


def test_quiet_groups_by_station_set():
    # Windows q1 and q3 on set 7, q2 on set 2; 5.5 h each. Set 7: 11 h, learned
    # 3 + 8 = 11 alarms, 24 per day; STA/LTA 0. Set 2: 5.5 h, 1 and 11 alarms,
    # 24 / 5.5 = 4.36 and 48 per day. The larger group first.
    hours = {"q1": 5.5, "q2": 5.5, "q3": 5.5}
    counts = {
        "learned": {"q1": 3, "q2": 1, "q3": 8},
        "stalta_eqfa": {"q1": 0, "q2": 11, "q3": 0},
    }
    g = he.quiet_groups(
        {"q1": 7, "q2": 2, "q3": 7}, hours, counts, ["learned", "stalta_eqfa"]
    )
    assert g == [
        {
            "set": "007",
            "windows": 2,
            "hours": 11.0,
            "alarms": {"learned": 11, "stalta_eqfa": 0},
            "per_day": {"learned": pytest.approx(24.0), "stalta_eqfa": 0.0},
        },
        {
            "set": "002",
            "windows": 1,
            "hours": 5.5,
            "alarms": {"learned": 1, "stalta_eqfa": 11},
            "per_day": {
                "learned": pytest.approx(24 / 5.5),
                "stalta_eqfa": pytest.approx(48.0),
            },
        },
    ]


# --- end to end ------------------------------------------------------------------------------------


def first(lon, hms, magnitude=2.04, kind="ML"):
    """A selection's sequence-initiating event on the equator."""
    return {
        "time": f"2024-07-06T{hms}",
        "lat": 0.0,
        "lon": lon,
        "depth_km": 11.08,
        "magnitude": magnitude,
        "type": kind,
        "id": "x",
    }


def make_network(tmp_path):
    """A run directory and a selection: sequence sets xx_now_seq001 (stations
    AAA, BBB, CCC) and xx_now_seq003 (BBB, CCC, DDD; sequence 2 is not
    selected), quiet windows 01 and 03 on set 1 and 02 on set 3."""
    run, sel = tmp_path / "run", tmp_path / "selection"
    sel.mkdir()
    (sel / "selection.json").write_text(
        json.dumps(
            {
                "group": "xx",
                "name": "now",
                "stations": {
                    "usable": [
                        {"network": "XX", "station": k, "lat": lat, "lon": lon}
                        for k, (lat, lon) in COORDS.items()
                    ]
                },
                "sequences": [
                    {"index": 1, "selected": True, "first": first(0.3, "00:01:05.30")},
                    {"index": 2, "selected": False},
                    {"index": 3, "selected": True, "first": first(0.4, "00:09:00.00")},
                ],
                "quiet": {
                    "windows": [{"sequence": 1}, {"sequence": 3}, {"sequence": 1}]
                },
            }
        )
    )
    none = {127: set(), -127: set()}
    write_outcomes(
        run / "outcomes",
        "learned",
        "xx_now_seq001",
        [1, 1, 1, 0],
        {127: {0}, -127: set()},
        ["1.9", "2.0", "4.0", "3.9"],
    )
    write_outcomes(
        run / "outcomes",
        "stalta_eqfa",
        "xx_now_seq001",
        [1, 1, 0, 0],
        none,
        ["1.9", "2.0", "4.0", "3.9"],
    )
    write_outcomes(run / "outcomes", "learned", "xx_now_seq003", [1, 0], none)
    write_outcomes(run / "outcomes", "stalta_eqfa", "xx_now_seq003", [0, 1], none)
    (run / "analysis").mkdir()
    with open(run / "analysis" / "magnitude_latency.csv", "w", newline="") as fh:
        rows = [
            {**r, "set": "xx_now_seq001" if r["set"] == "s1" else "xx_now_seq003"}
            for r in latency_rows()
        ]
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    other = ("BBB", "CCC", "DDD")
    quiet_other = make_report([(1, "06:00:00.00", "DDD, BBB", None, None)], [1], other)
    reports = {
        "runs": {
            "xx_now_seq001": SEQUENCE_REPORT,
            "xx_now_seq003": make_report(SEQUENCE_ALARMS[:2], [1], other),
        },
        "quiet": {
            "xx_now_quiet01": QUIET_REPORT,
            "xx_now_quiet02": quiet_other,
            "xx_now_quiet03": make_report([], []),
        },
    }
    for sub, sets in reports.items():
        for label in ("learned", "stalta_eqfa"):
            (run / sub / label).mkdir(parents=True)
            for s, text in sets.items():
                (run / sub / label / f"{s}.txt").write_text(text)
    # STA/LTA's third quiet window has three alarms.
    (run / "quiet" / "stalta_eqfa" / "xx_now_quiet03.txt").write_text(
        make_report(
            [(n, f"07:0{n}:00.00", "AAA, BBB", None, None) for n in (1, 2, 3)],
            [1, 2, 3],
        )
    )
    with open(run / "quiet" / "quiet_index.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "name",
                "start_utc",
                "score_from_utc",
                "end_utc",
                "scored_hours",
                "stations",
            ]
        )
        for j in (1, 2, 3):
            w.writerow([f"xx_now_quiet{j:02d}", "", "", "", "5.5000", 3])
    return run, sel


def test_main_end_to_end(tmp_path):
    run, sel = make_network(tmp_path)
    before = sorted(p for p in run.rglob("*"))
    out = tmp_path / "out.json"
    p = run_tool(
        "heldout_exploratory",
        "--run",
        run,
        "--selection",
        sel,
        "--detectors",
        "learned,stalta_eqfa",
        "--quiet-sets",
        "1",
        "--json",
        out,
        "--mainshocks",
        tmp_path / "mainshocks.csv",
    )
    assert p.returncode == 0, p.stderr
    assert sorted(run.rglob("*")) == before  # nothing written beside the inputs
    r = json.loads(out.read_text())
    assert p.stdout.startswith("== xx_now: exploratory, outside the frozen protocol\n")
    # 1: as test_paired_latency_with_a_tie.
    assert (
        r["paired_latency"]["events"],
        r["paired_latency"]["median"],
        r["paired_latency"]["ties"],
    ) == (4, 0.5, 1)
    assert (
        "median +0.50 s (interquartile range -0.12 to +1.25); stalta_eqfa first 2, learned first 1, ties 1"
        in p.stdout
    )
    # 2: the sequences of the selection only; unmatched from each report's list.
    assert list(r["sequences"]) == ["xx_now_seq001", "xx_now_seq003"]
    row = r["sequences"]["xx_now_seq001"]["learned"]
    assert (row["scored"], row["detected"], row["chance"], row["unmatched"]) == (
        4,
        3,
        0.125,
        4,
    )
    assert row["corrected"] == pytest.approx(20 / 7)
    assert r["sequences"]["xx_now_seq003"]["stalta_eqfa"]["unmatched"] == 1
    assert r["sequence_comparison"] == {"more": 1, "same": 1, "fewer": 0}
    # 3: n = 6, raw 4, corrected 20/7 + 1 = 27/7: c = (1/7) / (15/7) = 1/15.
    assert r["network_chance"]["learned"]["effective_chance"] == pytest.approx(1 / 15)
    # 4: set 1 holds windows 01 and 03 (11 h): learned 2 + 0 alarms, STA/LTA
    # 2 + 3 = 5 (10.9 per day); set 3 window 02: one alarm each.
    assert r["quiet_windows"]["xx_now_quiet03"] == {
        "set": "001",
        "learned": 0,
        "stalta_eqfa": 3,
    }
    assert [
        (g["set"], g["windows"], g["hours"], g["alarms"]) for g in r["quiet_sets"]
    ] == [
        ("001", 2, 11.0, {"learned": 2, "stalta_eqfa": 5}),
        ("003", 1, 5.5, {"learned": 1, "stalta_eqfa": 1}),
    ]
    assert r["quiet_sets"][0]["per_day"]["stalta_eqfa"] == pytest.approx(24 * 5 / 11)
    assert [(g["set"], g["windows"]) for g in r["quiet_named_sets"]] == [
        ("001", 2),
        ("other", 1),
    ]
    # 5: the stations that ran: AAA-BBB (0.3 degrees) on set 1, BBB-DDD (0.2) on set 3.
    assert r["closest_pairs"]["xx_now_seq001"] == {
        "pair": ["AAA", "BBB"],
        "km": pytest.approx(0.3 * DEG_KM),
    }
    assert r["closest_pairs"]["xx_now_quiet02"] == {
        "pair": ["BBB", "DDD"],
        "km": pytest.approx(0.2 * DEG_KM),
    }
    # Learned: matched alarm 2 of each sequence set (AAA, CCC): 0 of 2;
    # unmatched 1 and 4 of set 1 by AAA-BBB, the one of set 3 (AAA, BBB) is not
    # by BBB-DDD: 2 of 5; quiet 01 one of two (BBB, AAA), 02 one of one (DDD, BBB).
    assert r["closest_pair_alarms"]["learned"] == {
        "matched": {"alarms": 2, "by_closest_pair": 0},
        "unmatched": {"alarms": 5, "by_closest_pair": 2},
        "quiet": {"alarms": 3, "by_closest_pair": 2},
    }
    # 6: 1.9 | 2.0 and the two 2.0 of set 3 | 3.9 | 4.0.
    assert [
        (b["band"], b["scored"], b["learned"], b["stalta_eqfa"])
        for b in r["magnitude_bands"]
    ] == [("<2", 1, 1, 1), ("2-3", 3, 2, 2), ("3-4", 1, 0, 0), (">=4", 1, 1, 0)]
    # 7: unmatched 4 + 1 alarms, estimates for 3 and 4 of set 1, one final >= 4.0;
    # quiet 2 + 1 + 0 alarms, two estimates (-0.1 and 5.2).
    m = r["alarm_magnitudes"]["learned"]
    assert (
        m["unmatched"]["alarms"],
        m["unmatched"]["with_pd"],
        m["unmatched"]["large"],
    ) == (5, 2, 1)
    assert (m["quiet"]["alarms"], m["quiet"]["with_pd"], m["quiet"]["large"]) == (
        3,
        2,
        1,
    )
    assert m["quiet"]["large_cases"] == [
        {"set": "xx_now_quiet01", "time": "05:30:00.00", "pd_final": 5.2}
    ]
    assert "      xx_now_quiet01 05:30:00.00 M5.2\n" in p.stdout


def test_missing_outcome_table_stops(tmp_path):
    run, sel = make_network(tmp_path)
    (run / "outcomes" / "stalta_eqfa_xx_now_seq003.csv").unlink()
    p = run_tool(
        "heldout_exploratory",
        "--run",
        run,
        "--selection",
        sel,
        "--detectors",
        "learned,stalta_eqfa",
    )
    assert p.returncode != 0 and "no outcome table" in p.stderr


# --- item 8: another bulletin's events and the unmatched alarms after them -----------------------


def test_followed_events_counts_events_and_distinct_alarms():
    # By hand. Events at 06:45:11 and 06:45:29 lie 18 s apart: alarm 10 at
    # 06:45:40 follows both (29 s and 11 s), so two events are followed by one
    # alarm. The event at 06:30:16 has alarm 3 at 06:31:16.00, exactly 60 s
    # later, which counts; alarm 4 at 06:31:16.50 is 60.5 s later and does not
    # on its own. The event at 08:56:26 has no alarm within 60 s (alarm 9 is
    # 1 s before it). 3 of 4 events followed, by 2 distinct alarms.
    events = ["06:30:16.000", "06:45:11.452", "06:45:29.901", "08:56:26.795"]
    alarms = {3: "06:31:16.00", 4: "06:31:16.50", 10: "06:45:40.00", 9: "08:56:25.79"}
    got = he.followed_events(events, alarms, 60.0)
    assert (got["events"], got["followed"], got["distinct_alarms"]) == (4, 3, 2)
    assert got["followed_events"] == events[:3]


def test_followed_events_across_midnight():
    # 23:59:50 to 00:00:20 is 30 s.
    assert he.followed_events(["23:59:50.0"], {1: "00:00:20.00"}, 60.0)["followed"] == 1


def test_read_other_events_keeps_one_author(tmp_path):
    f = tmp_path / "isc.txt"
    f.write_text(
        "#EventID|Time|Latitude|Longitude|Depth/km|Author|Catalog\n"
        "1|2023-11-03T06:24:27.500|38.7|23.5|14.0|ATH|ATH\n"
        "2|2023-11-03T06:30:16.114|38.7|23.5|7.7|THE|THE\n"
        "# Agencies whose data contributed\n"
    )
    assert he.read_other_events(f) == ["06:24:27.500", "06:30:16.114"]
    assert he.read_other_events(f, "THE") == ["06:30:16.114"]


# --- item 9: the sequence-initiating events --------------------------------------------------------


def mainshock_block(number, time, by, afad=None, pd=None, pre=None, warning=()):
    """An event block in ayzek's layout. afad: (type, magnitude, origin, delay);
    pd: (first, seconds after the alarm or None, final); pre: (before, used);
    warning: (station, km, picked or predicted, signed seconds)."""
    lines = [
        f"Event #{number}: detected earthquake at 37.366N 38.835E, origin {time} UTC",
        f"  alarm      {time} UTC by {by}",
    ]
    if pd:
        when = "at the alarm" if pd[1] is None else f"{pd[1]} s after the alarm"
        lines.append(
            f"  Pd magn.   M{pd[0]} {when}, M{pd[2]} ± 0.1 (M<4 0%, 4-5.5 24%, M>=5.5 76%) "
            "final from 4 stations (4 above noise)"
        )
    lines.append("  location   rms 0.81 s from 6 stations, first 50.0 s after the alarm")
    lines.append(
        f"  AFAD       {afad[0]} {afad[1]} at {afad[2]} UTC: alarm {afad[3]} s after origin, "
        "epicentre 6.8 km off, origin +0.2 s"
        if afad
        else "  AFAD       no matching catalogue event"
    )
    if pre:
        lines.append(
            f"  pre-S      {pre[0]} of {pre[1]} alarm detections declared before the predicted S"
        )
    if warning:
        lines.append(
            "  S wave     had travelled 75 km from the epicentre when the alarm sounded"
        )
        lines.append(
            "  warning               distance   S wave arrives            warning"
        )
        lines += [
            f"             {st:5s}       {km:4d} km   07:41:29.72 UTC ({src})   {lead:+6.1f} s  "
            f"{'before' if lead > 0 else 'after'} S"
            for st, km, src, lead in warning
        ]
    return "\n".join(lines) + "\n"


MAINSHOCK = first(0.0, "17:34:29.746", 4.93)
MAINSHOCK_KEY = ("s1", "17:34:30.00", "ML", "4.9")
# The mainshock's block: two stations at the same printed 34 km, both before S.
OWN_BLOCK = mainshock_block(
    2,
    "17:34:37.70",
    "CCC, DDD; later AAA",
    ("ML", "4.9", "17:34:30.00", "7.7"),
    ("4.6", None, "4.8"),
    (1, 2),
    [
        ("AAA", 20, "picked", -1.0),
        ("BBB", 34, "picked", 3.4),
        ("CCC", 34, "predicted", 2.4),
        ("DDD", 60, "predicted", 9.0),
    ],
)


def mainshock_report(own=OWN_BLOCK):
    """Alarm 1 on an ML 1.9 19 s before the mainshock; 2 the mainshock's (or
    `own`); 3 unmatched, 59.754 s after the selection's origin; 4, 60.254 s
    after, and 5, 30.246 s before, lie outside the 30 s and 60 s."""
    return "ayzek    4 stations\n" + "-" * 78 + "\n\n" + "\n".join(
        [
            mainshock_block(
                1, "17:34:18.50", "AAA, BBB", ("ML", "1.9", "17:34:11.00", "7.5")
            ),
            own,
            mainshock_block(3, "17:35:29.50", "AAA, CCC", pd=("3.0", "1.5", "3.1")),
            mainshock_block(4, "17:35:30.00", "AAA, CCC"),
            mainshock_block(5, "17:33:59.50", "AAA, CCC"),
        ]
    ) + "\n" + "-" * 78 + "\nSUMMARY\n"


def test_seconds_after_within_half_a_day():
    assert he.seconds_after("00:00:10.00", "23:59:50.50") == pytest.approx(19.5)
    assert he.seconds_after("23:59:50.50", "00:00:10.00") == pytest.approx(-19.5)
    assert he.seconds_after("17:34:18.50", "17:34:29.746") == pytest.approx(-11.246)


def test_parse_events_reads_the_block_lines():
    ev = he.parse_events(mainshock_report())
    assert sorted(ev) == [1, 2, 3, 4, 5]
    assert ev[2] == {
        "time": "17:34:37.70",
        "declaring": ["CCC", "DDD"],  # not the stations after "; later"
        "pd_first": 4.6,
        "pd_first_after_alarm_s": 0.0,  # "at the alarm"
        "pd_final": 4.8,
        "afad": ("ML", "4.9", "17:34:30.00"),
        "delay": 7.7,
        "pre_s": (1, 2),
        "warning": [
            ("AAA", 20, "picked", -1.0, "after"),
            ("BBB", 34, "picked", 3.4, "before"),
            ("CCC", 34, "predicted", 2.4, "before"),
            ("DDD", 60, "predicted", 9.0, "before"),
        ],
    }
    assert (ev[1]["afad"], ev[1]["delay"], ev[1]["pd_final"], ev[1]["pre_s"]) == (
        ("ML", "1.9", "17:34:11.00"),
        7.5,
        None,
        None,
    )
    assert (ev[3]["afad"], ev[3]["pd_first"], ev[3]["pd_first_after_alarm_s"]) == (
        None,
        3.0,
        1.5,
    )
    assert ev[3]["warning"] == [] and ev[3]["declaring"] == ["AAA", "CCC"]


def test_scored_mainshock_matches_the_whole_second_origin():
    # The catalogue's origin is the selection's to the second: 29.746 is 30.00
    # (0.254 s off), 4.93 prints as 4.9. An event of the same type and
    # magnitude 5.5 min later is not it, nor the ML 1.9 19 s before.
    events = [
        ("s1", "17:34:11.00", "ML", "1.9"),
        MAINSHOCK_KEY,
        ("s1", "17:40:00.00", "ML", "4.9"),
    ]
    assert he.scored_mainshock(MAINSHOCK, events) == MAINSHOCK_KEY
    # Truncated instead of rounded (0.746 s off) still matches; the type's case
    # does not matter; of two within a second the nearer in time.
    assert he.scored_mainshock(MAINSHOCK, [("s1", "17:34:29.00", "ml", "4.9")])
    assert (
        he.scored_mainshock(
            MAINSHOCK, [("s1", "17:34:29.00", "ML", "4.9"), MAINSHOCK_KEY]
        )
        == MAINSHOCK_KEY
    )
    # 1.254 s off; another type; a magnitude 0.07 off (4.93 against 5.0, and 4.86
    # would print 4.9): no scored event.
    assert he.scored_mainshock(MAINSHOCK, [("s1", "17:34:31.00", "ML", "4.9")]) is None
    assert he.scored_mainshock(MAINSHOCK, [("s1", "17:34:30.00", "Mw", "4.9")]) is None
    assert he.scored_mainshock(MAINSHOCK, [("s1", "17:34:30.00", "ML", "5.0")]) is None
    assert he.scored_mainshock(MAINSHOCK, []) is None


def test_station_distances_from_the_epicentre():
    # Epicentre at AAA (0, 0): BBB 0.3 degrees east, DDD 0.5, CCC 1 north.
    d = he.station_distances(MAINSHOCK, COORDS, ["CCC", "DDD", "BBB"])
    assert [s for _, s in d] == ["BBB", "DDD", "CCC"]
    assert [km for km, _ in d] == pytest.approx([0.3 * DEG_KM, 0.5 * DEG_KM, DEG_KM])


def test_mainshock_detected_with_its_block():
    # The block credited to the mainshock is alarm 2: 7.7 s after the catalogue's
    # 17:34:30.00 as printed, 37.70 - 29.746 = 7.954 s after the selection's
    # origin. Three of four warning rows are before S; the nearest are BBB and
    # CCC at 34 km, of which the smaller warning, CCC's 2.4 s, is taken.
    r = he.mainshock_detector(
        MAINSHOCK, MAINSHOCK_KEY, 1, 0, he.parse_events(mainshock_report()), (30, 60)
    )
    note = r.pop("note")
    assert r == {
        "status": "detected",
        "detected_pre_s": 0,
        "alarm": 2,
        "alarm_time": "17:34:37.70",
        "alarm_delay_s": 7.7,
        "alarm_after_selection_origin_s": pytest.approx(7.95),
        "declaring": ["CCC", "DDD"],
        "pre_s": "1 of 2",
        "pd_first": 4.6,
        "pd_first_after_alarm_s": 0.0,
        "pd_final": 4.8,
        "warning_stations": 4,
        "stations_ahead": 3,
        "nearest_ahead_station": "CCC",
        "nearest_ahead_km": 34,
        "nearest_ahead_warning_s": 2.4,
        "nearest_ahead_s_source": "predicted",
    }
    # Alarms 1 (-11.246 s) and 3 (+59.754 s) are within 30 s before to 60 s
    # after; 4 (+60.254 s) and 5 (-30.246 s) are not; 2 is the row itself.
    assert note == [
        "alarm #1 at 17:34:18.50 (-11.2 s) by AAA, BBB credited to ML 1.9 at 17:34:11.00",
        "alarm #3 at 17:35:29.50 (+59.8 s) by AAA, CCC credited to no catalogue event, "
        "Pd final M3.1",
    ]
    # A wider window takes 4 and 5 in.
    wide = he.mainshock_detector(
        MAINSHOCK, MAINSHOCK_KEY, 1, 0, he.parse_events(mainshock_report()), (31, 61)
    )
    assert [n.split()[1] for n in wide["note"]] == ["#1", "#3", "#4", "#5"]


def test_mainshock_missed_keeps_the_other_alarms_in_the_note():
    # No block names the mainshock: alarm 2, 7.954 s after its origin, is
    # credited to no event. The row has no delay, estimate or warning; the
    # alarms around the origin are in the note only.
    unmatched = mainshock_block(2, "17:34:37.70", "CCC, DDD", pd=("4.6", None, "4.8"))
    events = he.parse_events(mainshock_report(unmatched))
    r = he.mainshock_detector(MAINSHOCK, MAINSHOCK_KEY, 0, 0, events, (30, 60))
    assert r == {
        "status": "missed",
        "detected_pre_s": 0,
        "note": [
            "alarm #1 at 17:34:18.50 (-11.2 s) by AAA, BBB credited to ML 1.9 at 17:34:11.00",
            "alarm #2 at 17:34:37.70 (+8.0 s) by CCC, DDD credited to no catalogue event, "
            "Pd final M4.8",
            "alarm #3 at 17:35:29.50 (+59.8 s) by AAA, CCC credited to no catalogue event, "
            "Pd final M3.1",
        ],
    }
    # The outcome table saying detected against such a report is noted, and so
    # is the reverse; the status stays the table's.
    r = he.mainshock_detector(MAINSHOCK, MAINSHOCK_KEY, 1, 1, events, (30, 60))
    assert r["status"] == "detected" and "alarm" not in r
    assert r["note"][0] == (
        "outcome table detected, the report has no event block credited to it"
    )
    r = he.mainshock_detector(
        MAINSHOCK, MAINSHOCK_KEY, 0, 0, he.parse_events(mainshock_report()), (30, 60)
    )
    assert r["status"] == "missed" and r["alarm"] == 2
    assert r["note"][0] == (
        "outcome table missed, the report has an event block credited to it"
    )


def test_mainshock_not_scored_and_two_credited_alarms():
    # Without an event in the outcome table nothing of the report is the
    # mainshock's, whatever its blocks say.
    events = he.parse_events(mainshock_report())
    r = he.mainshock_detector(MAINSHOCK, None, None, None, events, (30, 60))
    assert (r["status"], r["detected_pre_s"], "alarm" in r) == ("not scored", None, False)
    assert [n.split()[1] for n in r["note"]] == ["#1", "#2", "#3"]
    # Two blocks credited to it: the earlier alarm (6, at 17:34:36.00) is the
    # row, the other is noted.
    second = mainshock_block(
        6, "17:34:36.00", "AAA, DDD", ("ML", "4.9", "17:34:30.00", "6.0")
    )
    events = he.parse_events(mainshock_report(OWN_BLOCK + "\n" + second))
    r = he.mainshock_detector(MAINSHOCK, MAINSHOCK_KEY, 1, 1, events, (30, 60))
    assert (r["alarm"], r["alarm_delay_s"], r["declaring"]) == (6, 6.0, ["AAA", "DDD"])
    assert (r["warning_stations"], r["stations_ahead"]) == (0, 0)
    assert "nearest_ahead_station" not in r
    assert r["note"][0] == "2 alarms credited to it (#6, #2), the earliest taken"
    assert (
        "alarm #2 at 17:34:37.70 (+8.0 s) by CCC, DDD credited to it too, Pd final M4.8"
        in r["note"]
    )


def test_earlier_events_and_joined_events_note():
    events = [
        ("s1", "17:30:00.00", "ML", "2.2"),
        ("s1", "17:34:11.00", "ML", "1.9"),
        MAINSHOCK_KEY,
        ("s1", "17:40:00.00", "ML", "4.9"),
    ]
    # Two scored events before 17:34:29.746; the last 18.746 s before.
    joined = [first(0.0, "17:50:01.655362", 4.926)]
    assert he.earlier_events_note(MAINSHOCK, MAINSHOCK_KEY, events, joined) == [
        "2 scored events earlier in the set, the last ML 1.9 at 17:34:11.00 (18.7 s before)",
        "the selection joined 1 later event: ML 4.9 at 17:50:01.65",
    ]
    assert he.earlier_events_note(MAINSHOCK, MAINSHOCK_KEY, events[2:], []) == []


def test_mainshock_summary_by_hand():
    # Five mainshocks: three detected (delays 7.7, 10.0, 20.0: median 10.0), two
    # of them before S by the outcome table, one with a Pd estimate at the alarm
    # and two with one at all, two with a station before S (warnings 2.4 and
    # 6.0: median 4.2); one missed; one not scored.
    rows = [
        {
            "status": "detected",
            "detected_pre_s": 1,
            "alarm_delay_s": 7.7,
            "pd_first_after_alarm_s": 0.0,
            "pd_final": 4.8,
            "nearest_ahead_warning_s": 2.4,
        },
        {
            "status": "detected",
            "detected_pre_s": 0,
            "alarm_delay_s": 20.0,
            "pd_first_after_alarm_s": 2.5,
            "pd_final": 5.0,
            "nearest_ahead_warning_s": 6.0,
        },
        {
            "status": "detected",
            "detected_pre_s": 1,
            "alarm_delay_s": 10.0,
            "pd_first_after_alarm_s": None,
            "pd_final": None,
        },
        {"status": "missed", "detected_pre_s": 0},
        {"status": "not scored", "detected_pre_s": None},
    ]
    assert he.mainshock_summary(rows) == {
        "mainshocks": 5,
        "scored": 4,
        "detected": 3,
        "detected_pre_s": 2,
        "pd_at_alarm": 1,
        "pd_any": 2,
        "station_ahead": 2,
        "median_alarm_delay_s": 10.0,
        "median_nearest_ahead_warning_s": pytest.approx(4.2),
    }
    none = he.mainshock_summary(rows[3:])
    assert (none["detected"], none["median_alarm_delay_s"]) == (0, None)
    assert none["median_nearest_ahead_warning_s"] is None


def test_mainshock_table_end_to_end(tmp_path):
    # Sequence 1 (stations AAA, BBB, CCC): the selection's first event, ML 2.04
    # at 00:01:05.30 and 0.3 degrees east on the equator (at BBB; AAA 0.3
    # degrees away, CCC further), is the outcome table's ML 2.0 at 00:01:05.00:
    # detected by both detectors there, with one scored event (ML 1.9 at
    # 00:00:05.00, 60.3 s) before it. The reports' blocks name another event, which is noted.
    # Sequence 3: no event of the table is within a second of 00:09:00.00.
    run, sel = make_network(tmp_path)
    out, table = tmp_path / "out.json", tmp_path / "mainshocks.csv"
    p = run_tool(
        "heldout_exploratory",
        "--run",
        run,
        "--selection",
        sel,
        "--detectors",
        "learned,stalta_eqfa",
        "--json",
        out,
        "--mainshocks",
        table,
    )
    assert p.returncode == 0, p.stderr
    m = json.loads(out.read_text())["mainshocks"]
    assert m["alarm_window_s"] == [30.0, 60.0]
    assert [e["set"] for e in m["events"]] == ["xx_now_seq001", "xx_now_seq003"]
    e1, e3 = m["events"]
    assert e1["scored_as"] == "00:01:05.00 ML 2.0"
    assert [x["station"] for x in e1["nearest"]] == ["BBB", "AAA"]
    assert [x["km"] for x in e1["nearest"]] == pytest.approx([0.0, 0.3 * DEG_KM])
    assert e1["note"] == [
        "1 scored event earlier in the set, the last ML 1.9 at 00:00:05.00 (60.3 s before)"
    ]
    assert e1["detectors"]["learned"] == {
        "status": "detected",
        "detected_pre_s": 0,
        "note": ["outcome table detected, the report has no event block credited to it"],
    }
    assert e3["scored_as"] is None
    # Its two scored events (00:00:05.00 and 00:01:05.00) are both earlier.
    assert e3["note"] == [
        "no such event in the outcome table",
        "2 scored events earlier in the set, the last ML 2.0 at 00:01:05.00 (475.0 s before)",
    ]
    assert e3["detectors"]["stalta_eqfa"]["status"] == "not scored"
    assert m["summary"]["learned"]["mainshocks"] == 2
    assert (m["summary"]["learned"]["scored"], m["summary"]["learned"]["detected"]) == (
        1,
        1,
    )
    assert "   learned                2 / 1 / 1 / 0 / 0 (0) / 0; none s; none s\n" in p.stdout
    rows = list(csv.DictReader(open(table)))
    assert list(rows[0]) == he.MAINSHOCK_FIELDS
    assert [(r["set"], r["detector"], r["status"]) for r in rows] == [
        ("xx_now_seq001", "learned", "detected"),
        ("xx_now_seq001", "stalta_eqfa", "detected"),
        ("xx_now_seq003", "learned", "not scored"),
        ("xx_now_seq003", "stalta_eqfa", "not scored"),
    ]
    assert (rows[0]["magnitude"], rows[0]["depth_km"], rows[0]["nearest_km"]) == (
        "2.04",
        "11.1",
        "0.0",
    )
    assert rows[0]["second_station"] == "AAA" and rows[0]["alarm_delay_s"] == ""
    assert rows[0]["note"] == (
        "1 scored event earlier in the set, the last ML 1.9 at 00:00:05.00 (60.3 s before); "
        "outcome table detected, the report has no event block credited to it"
    )
