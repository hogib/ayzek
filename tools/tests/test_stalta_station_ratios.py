"""tools/stalta_station_ratios.py: times, distances, the visibility rule and
the quiet-time mask, rates, the recording parser, and main() (each group's
target from its own tuning windows, and the ratio rule), as rewritten in
943e031.

Expected values: calendar.timegm and the spherical law of cosines as
independent arithmetic, the mask and visibility rule as docs/protocol-heldout.md
section 4 states them, small cases by hand, and a fake ayzek that writes the
recording the test prescribes.
"""
import calendar
import json
import math
import time

import pytest

import stalta_station_ratios as sr
from conftest import run_tool, write_fake_binary


def tg(*ymdhms):
    return calendar.timegm((*ymdhms, 0, 0, 0))


# --- epoch, km, catalogue -----------------------------------------------------------------

@pytest.mark.parametrize("text, want", [
    ("2025-11-10T18:04:00", tg(2025, 11, 10, 18, 4, 0)),
    ("2025-11-10T18:04:00+00:00", tg(2025, 11, 10, 18, 4, 0)),
    ("2025-11-10T20:04:00+02:00", tg(2025, 11, 10, 18, 4, 0)),
    ("10/11/2025 18:04:00", tg(2025, 11, 10, 18, 4, 0)),
    ("10/11/2025 18:04:00.5", tg(2025, 11, 10, 18, 4, 0) + 0.5),
    ("  01/01/2026 00:00:00 \n", tg(2026, 1, 1, 0, 0, 0)),
])
def test_epoch(text, want):
    assert sr.epoch(text) == want


def test_km_by_hand():
    # One degree along a meridian: 6371 km * pi/180 = 111.1949 km; half the
    # earth: pi * 6371 = 20015.09 km; the same point: 0.
    assert sr.km(0, 0, 1, 0) == pytest.approx(6371 * math.pi / 180, rel=1e-9)
    assert sr.km(0, 0, 0, 1) == pytest.approx(6371 * math.pi / 180, rel=1e-9)
    assert sr.km(0, 0, 0, 180) == pytest.approx(math.pi * 6371, rel=1e-9)
    assert sr.km(39.2, 28.1, 39.2, 28.1) == 0
    assert sr.km(39.2, 28.1, 38.1, 27.0) == sr.km(38.1, 27.0, 39.2, 28.1)


def test_km_against_the_law_of_cosines():
    def cosines(lat1, lon1, lat2, lon2):
        p1, p2 = math.radians(lat1), math.radians(lat2)
        return 6371 * math.acos(math.sin(p1) * math.sin(p2)
                                + math.cos(p1) * math.cos(p2) * math.cos(math.radians(lon2 - lon1)))
    for a, b in [((39.0, 28.0), (39.5, 28.5)), ((39.24, 28.15), (38.51, 28.66)),
                 ((-39.37, 175.53), (-41.0, 174.0)), ((37.5, 22.1), (39.75, 22.13))]:
        assert sr.km(*a, *b) == pytest.approx(cosines(*a, *b), rel=1e-9)


def test_catalogue_reads_ayzek_format(tmp_path):
    f = tmp_path / "cat.csv"
    f.write_text("﻿Date,Longitude,Latitude,Depth,Rms,Type,Magnitude,Location,EventID\n"
                 "10/11/2025 18:20:51,28.145,39.23889,9.64,0.56,MW,4.9,Sındırgı (Balıkesir),692187\n"
                 "27/01/2022 05:20:02.5,175.52683,-39.37058,9.54,0.45,MLv,1.29,King Country,x\n",
                 encoding="utf-8")
    assert sr.catalogue(f) == [(tg(2025, 11, 10, 18, 20, 51), 39.23889, 28.145, 4.9),
                               (tg(2022, 1, 27, 5, 20, 2) + 0.5, -39.37058, 175.52683, 1.29)]


# --- arrivals: visibility and the mask ------------------------------------------------------

def north(km):
    """A point `km` north of (0, 0): one degree of latitude is 6371*pi/180 km."""
    return (km / (6371 * math.pi / 180), 0.0)


@pytest.mark.parametrize("km, mag, visible", [
    (49.9, 0.0, True), (50.1, 1.9, False),        # within 50 km at any magnitude
    (149.9, 2.0, True), (150.1, 2.0, False), (149.9, 1.9, False),   # 150 km from M2
    (399.0, 3.0, True), (401.0, 3.0, False), (300.0, 2.9, False),   # 400 km from M3
    (1499.0, 4.5, True), (1501.0, 4.5, False), (500.0, 4.4, False),  # 1500 km from M4.5
])
def test_visibility_rule(km, mag, visible):
    lat, lon = north(km)
    spans = sr.arrivals([(1000.0, lat, lon, mag)], 0.0, 0.0, 0.0, 1e5)
    assert bool(spans) is visible


def test_mask_window_by_hand():
    # Event at t = 1000, 80 km away: [1000 + 80/8 - 10, 1000 + 80/6 + 120]
    # = [1000, 1133.33]. At the station itself (0 km): [990, 1120].
    lat, lon = north(80.0)
    (a, b), = sr.arrivals([(1000.0, lat, lon, 3.0)], 0.0, 0.0, 0.0, 5000.0)
    assert a == pytest.approx(1000.0, abs=1e-6) and b == pytest.approx(1000 + 80 / 6 + 120, abs=1e-6)
    assert sr.arrivals([(1000.0, 0.0, 0.0, 1.0)], 0.0, 0.0, 0.0, 5000.0) == [(990.0, 1120.0)]


def test_mask_clipped_and_merged():
    ev = (1000.0, 0.0, 0.0, 1.0)                      # mask [990, 1120]
    assert sr.arrivals([ev], 0, 0, 1050.0, 5000.0) == [(1050.0, 1120.0)]
    assert sr.arrivals([ev], 0, 0, 0.0, 1100.0) == [(990.0, 1100.0)]
    assert sr.arrivals([ev], 0, 0, 1120.0, 5000.0) == []     # a < b fails
    assert sr.arrivals([ev], 0, 0, 0.0, 990.0) == []
    # A second event 100 s later overlaps ([1090, 1220]) and merges; one
    # 200 s later touches the end exactly (1120) and merges too; one 500 s
    # later is a span of its own.
    assert sr.arrivals([ev, (1100.0, 0.0, 0.0, 1.0)], 0, 0, 0.0, 5000.0) == [(990.0, 1220.0)]
    assert sr.arrivals([ev, (1130.0, 0.0, 0.0, 1.0)], 0, 0, 0.0, 5000.0) == [(990.0, 1250.0)]
    assert sr.arrivals([(1500.0, 0.0, 0.0, 1.0), ev], 0, 0, 0.0, 5000.0) == [(990.0, 1120.0),
                                                                              (1490.0, 1620.0)]


# --- quiet_counts ---------------------------------------------------------------------------

def test_quiet_counts_by_hand():
    # Two hours with one masked event [3590, 3720] (130 s): quiet hours
    # (7200 - 130) / 3600. Triggers: 0 (counted, t0 <= t), 3589.9 (counted),
    # 3590 (masked: a <= t), 3719.9 (masked), 3720 (counted: t < b fails),
    # 7199.9 (counted), 7200 (not in [t0, t1)), -1 (before t0) -> 4.
    stations = {"A": (0.0, 0.0), "B": (0.0, 0.0), "C": (5.0, 5.0)}
    times = {"A": [0, 3589.9, 3590, 3719.9, 3720, 7199.9, 7200, -1], "B": []}
    out = sr.quiet_counts(stations, times, [(3600.0, 0.0, 0.0, 1.0)], 0.0, 7200.0)
    assert out["A"] == (pytest.approx(7070 / 3600), 4)
    assert out["B"] == (pytest.approx(7070 / 3600), 0)
    # C is 786 km away: the M1 event is not visible there, nothing masked.
    assert out["C"] == (2.0, 0)


def test_quiet_counts_without_events():
    assert sr.quiet_counts({"A": (0.0, 0.0)}, {"A": [10.0, 20.0]}, [], 0.0, 3600.0) == {"A": (1.0, 2)}


# --- replay: the recording parser ----------------------------------------------------------------

FAKE_RECORDER = """
import os, sys
argv = sys.argv[1:]
if os.environ.get("FAKE_EXIT"):
    sys.exit(int(os.environ["FAKE_EXIT"]))
if os.environ.get("FAKE_NO_REC"):
    sys.exit(0)
open(argv[argv.index("--record") + 1], "w").write(os.environ["FAKE_REC"])
"""

REC = ("# ayzek recording v1\n# min_stations 2\n"
       "S\tAAA\t39.1\t28.2\nS\tBBB\t39.3\t28.4\nT\t100\t200\n"
       "D\tAAA\t110\t116.5\t0.9\t0.1\t0\t113.5\n"
       "P\tAAA\t110\t113.4\t120.1\t0.8\t0.7\t170\t0.3\n"
       "D\tBBB\t120\t125\t0.8\t0.1\t0\t123.5\n"
       "D\tAAA\t150\t155\t0.95\t0.1\t0\t153.5\n")


def test_replay_parses_stations_and_detection_times(tmp_path, monkeypatch):
    fake = write_fake_binary(tmp_path / "ayzek", FAKE_RECORDER)
    monkeypatch.setenv("FAKE_REC", REC)
    stations, times = sr.replay(str(fake), ["--detector", "stalta"], ["x.mseed"])
    assert stations == {"AAA": (39.1, 28.2), "BBB": (39.3, 28.4)}
    assert dict(times) == {"AAA": [116.5, 155.0], "BBB": [125.0]}


@pytest.mark.parametrize("var", ["FAKE_EXIT", "FAKE_NO_REC"])
def test_replay_stops_when_ayzek_fails(tmp_path, monkeypatch, var):
    fake = write_fake_binary(tmp_path / "ayzek", FAKE_RECORDER)
    monkeypatch.setenv(var, "3")
    with pytest.raises(SystemExit):
        sr.replay(str(fake), [], ["x.mseed"])


# --- main ----------------------------------------------------------------------------------------

FAKE_AYZEK = """
import json, os, sys
argv = sys.argv[1:]
with open(os.environ["FAKE_LOG"], "a") as f:
    f.write(json.dumps(argv) + "\\n")
sc = json.load(open(os.environ["FAKE_SCENARIO"]))
ratio = argv[argv.index("--stalta-on") + 1]
run = sc[os.path.basename(os.path.dirname(argv[argv.index("--stations") + 1]))]
with open(argv[argv.index("--record") + 1], "w") as f:
    f.write("# ayzek recording v1\\n# min_stations 2\\n")
    for code, (lat, lon) in run["stations"].items():
        f.write(f"S\\t{code}\\t{lat}\\t{lon}\\n")
    f.write("T\\t0\\t1\\n")
    for code, times in run["times"].get(ratio, {}).items():
        for t in times:
            f.write(f"D\\t{code}\\t{t - 3.5}\\t{t}\\t0.9\\t0.1\\t0\\t{t}\\n")
"""

HEADER = "Date,Longitude,Latitude,Depth,Rms,Type,Magnitude,Location,EventID\n"
T1, T2, T3 = tg(2022, 1, 15, 0, 0, 0), tg(2022, 2, 1, 0, 0, 0), tg(2026, 7, 4, 0, 0, 0)
SPAN_S = 23400          # 6.5 h sets: a 30 min lead-in, 6 h counted


def iso(t):
    return time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(t))


def ddmmyyyy(t):
    return time.strftime("%d/%m/%Y %H:%M:%S", time.gmtime(t))


def spaced(t0, n):
    """n triggers 100 s apart from t0 + 50 s: all within the first 2000 s."""
    return [t0 + 50 + 100 * i for i in range(n)]


def tuning_set(root, name, start, stations, pd_gain, catalogue=""):
    d = root / name
    d.mkdir()
    (d / "S.mseed").write_bytes(b"")
    (d / "dataset.json").write_text(json.dumps({"span": [iso(start), iso(start + SPAN_S)]}))
    (d / "stations.csv").write_text("code,lat,lon\n" + "".join(
        f"{c},{lat},{lon}\n" for c, (lat, lon) in stations.items()))
    (d / "catalog.csv").write_text(HEADER + catalogue)
    if pd_gain:
        (d / "pd_gain.csv").write_text("station,gain\n")
    return d


@pytest.fixture
def main_env(tmp_path):
    # Grid 4, 6, 8, the network-wide ratio 8. Counted span of a set: from
    # 30 min after its start (LEAD_S) to its end, 6 h.
    #
    # Group gr: X in both sets (12 h), Y and W in set 1 (6 h). Triggers in the
    # counted span per ratio 4/6/8:
    #   X: set 1 10/3/3 (plus one in the lead-in and one at the end of the
    #      span, at every ratio: not counted), set 2 4/3/2 -> 14/6/5 over 12 h:
    #      1.1667, 0.5000, 0.4167
    #   Y: 2/1/0 over 6 h: 0.3333, 0.1667, 0.0000
    #   W: 6/3/3 over 6 h: 1.0000, 0.5000, 0.5000
    # Target gr = median of the rates at 8 (0.4167, 0, 0.5) = 0.4167 (X's), so
    #   X meets it at 8 only, by equality (rate <= target)      -> 8, met
    #   Y meets it at the lowest ratio                           -> 4, met
    #   W is above it at every ratio                             -> 8 (the
    #      highest of the grid), not met.
    # Set 1's catalogue has one M1.0 event 170 km from X: not visible (50 km
    # at M < 2), nothing masked.
    #
    # Group ko: Z alone, 20/10/5 triggers plus one at 3000 s into the span at
    # every ratio, inside the mask of set 3's M1.0 event at Z's own position at
    # that time ([3000 - 10, 3000 + 120] after the origin, 130 s): quiet hours
    # 6 - 130/3600 = 5.9639, rates 20/5.9639 = 3.3535, 1.6768, 0.8384; target
    # = the median of one rate = 0.8384; Z meets it at 8 -> 8, met.
    x, y, w, z = (39.0, 28.0), (39.5, 28.5), (39.1, 28.1), (38.7, 23.4)
    c1, c2, c3 = T1 + 1800, T2 + 1800, T3 + 1800
    sets = [tuning_set(tmp_path, "tune_gr_1", T1, {"X": x, "Y": y, "W": w}, True,
                       f"{ddmmyyyy(c1 + 600)},30.0,39.0,10,0.3,ML,1.0,far,1\n"),
            tuning_set(tmp_path, "tune_gr_2", T2, {"X": x}, False),
            tuning_set(tmp_path, "tune_ko_1", T3, {"Z": z}, False,
                       f"{ddmmyyyy(c3 + 3000)},{z[1]},{z[0]},10,0.3,ML,1.0,here,2\n")]
    scenario = {
        "tune_gr_1": {"stations": {"X": x, "Y": y, "W": w},
                      "times": {r: {"X": [T1 + 600, *spaced(c1, n), T1 + SPAN_S],
                                    "Y": spaced(c1, m), "W": spaced(c1, k)}
                                for r, n, m, k in (("4", 10, 2, 6), ("6", 3, 1, 3), ("8", 3, 0, 3))}},
        "tune_gr_2": {"stations": {"X": x},
                      "times": {r: {"X": spaced(c2, n)} for r, n in (("4", 4), ("6", 3), ("8", 2))}},
        "tune_ko_1": {"stations": {"Z": z},
                      "times": {r: {"Z": [*spaced(c3, n), c3 + 3000]}
                                for r, n in (("4", 20), ("6", 10), ("8", 5))}}}
    (tmp_path / "sc.json").write_text(json.dumps(scenario))
    fake = write_fake_binary(tmp_path / "ayzek", FAKE_AYZEK)
    log = tmp_path / "argv.jsonl"
    out = tmp_path / "out"

    def run(*args, tuning=None, on="8"):
        p = run_tool("stalta_station_ratios", "--binary", fake, "--on", on, "--grid", "4,6,8",
                     "--jobs", 1, "--out", out, *args, "--tuning", *(tuning or sets),
                     env={"FAKE_LOG": str(log), "FAKE_SCENARIO": str(tmp_path / "sc.json")})
        argvs = [json.loads(x) for x in log.read_text().splitlines()] if log.exists() else []
        return p, argvs

    run.out, run.tmp = out, tmp_path
    return run


def test_main_target_rates_and_ratios(main_env):
    p, argvs = main_env("--args", "--detector stalta --stalta-band 1,10")
    assert p.returncode == 0, p.stderr
    out = main_env.out
    assert (out / "targets.csv").read_text().splitlines() == [
        "group,stations,sets,median_per_hour_at_8", "gr,3,2,0.4167", "ko,1,1,0.8384"]
    assert (out / "rates_gr.csv").read_text().splitlines() == [
        "station,quiet_hours,triggers_4,triggers_6,triggers_8,per_hour_4,per_hour_6,per_hour_8,"
        "ratio,meets_target",
        "W,6.00,6,3,3,1.0000,0.5000,0.5000,8,0",
        "X,12.00,14,6,5,1.1667,0.5000,0.4167,8,1",
        "Y,6.00,2,1,0,0.3333,0.1667,0.0000,4,1"]
    assert (out / "rates_ko.csv").read_text().splitlines()[1:] == [
        "Z,5.96,20,10,5,3.3535,1.6768,0.8384,8,1"]
    assert (out / "stalta_on_gr.csv").read_text() == "station,ratio\nW,8\nX,8\nY,4\n"
    assert (out / "stalta_on_ko.csv").read_text() == "station,ratio\nZ,8\n"
    assert p.stdout.splitlines() == [
        "gr: 3 stations over 2 sets, target 0.4167 triggers per quiet hour at ratio 8; "
        "ratios 4: 1, 6: 0, 8: 2; 1 above the target at every ratio",
        "ko: 1 stations over 1 sets, target 0.8384 triggers per quiet hour at ratio 8; "
        "ratios 4: 0, 6: 0, 8: 1; 0 above the target at every ratio"]


def test_main_target_follows_the_network_wide_ratio(main_env):
    # --on 4: target gr = median(1.1667, 0.3333, 1.0) = 1.0 (W's); X meets it
    # at 6 (0.5), Y at 4, W at 4 (equality). ko: target 3.3535, Z -> 4.
    p, _ = main_env("--args", "--detector stalta", on="4")
    assert p.returncode == 0, p.stderr
    assert (main_env.out / "targets.csv").read_text().splitlines()[1:] == ["gr,3,2,1.0000", "ko,1,1,3.3535"]
    assert (main_env.out / "stalta_on_gr.csv").read_text() == "station,ratio\nW,4\nX,6\nY,4\n"
    assert (main_env.out / "stalta_on_ko.csv").read_text() == "station,ratio\nZ,4\n"


def test_main_runs_ayzek_as_documented(main_env):
    p, argvs = main_env("--args", "--detector stalta --stalta-band 1,10")
    assert p.returncode == 0, p.stderr
    # 3 sets x 3 ratios, each a recorded replay of the set's files with its
    # own station table, the --args and the grid's ratio.
    assert len(argvs) == 9
    by_set = {}
    for a in argvs:
        pairs = {(x, y) for x, y in zip(a, a[1:])}
        assert {("--speed", "0"), ("--detector", "stalta"), ("--stalta-band", "1,10")} <= pairs
        assert "--no-color" in a and "--record" in a and "--from" not in a
        assert a[-1].endswith("/S.mseed")
        by_set.setdefault(a[a.index("--stations") + 1].split("/")[-2], []).append(a)
    assert {k: sorted(a[a.index("--stalta-on") + 1] for a in v) for k, v in by_set.items()} == {
        "tune_gr_1": ["4", "6", "8"], "tune_gr_2": ["4", "6", "8"], "tune_ko_1": ["4", "6", "8"]}
    assert all(a[a.index("--stations") + 1].endswith("tune_gr_1/stations.csv") for a in by_set["tune_gr_1"])
    # --pd-gains only where the set has pd_gain.csv.
    assert all(a[a.index("--pd-gains") + 1].endswith("tune_gr_1/pd_gain.csv") for a in by_set["tune_gr_1"])
    assert all("--pd-gains" not in a for a in by_set["tune_gr_2"] + by_set["tune_ko_1"])


@pytest.mark.parametrize("args", ["--detector stalta --stalta-on 8",
                                  "--detector stalta --stalta-on-stations x.csv"])
def test_main_refuses_a_ratio_in_args(main_env, args):
    p, argvs = main_env("--args", args)
    assert p.returncode != 0 and "without its ratio" in p.stderr and not argvs


def test_main_refuses_a_ratio_outside_the_grid(main_env):
    p, argvs = main_env("--args", "--detector stalta", on="5")
    assert p.returncode != 0 and "must be one of the grid" in p.stderr and not argvs


@pytest.mark.parametrize("name", ["quiet_ko", "tune_ko_x", "tune_KO_1"])
def test_main_refuses_a_directory_that_is_not_a_tuning_set(main_env, name):
    bad = main_env.tmp / name
    bad.mkdir()
    (bad / "S.mseed").write_bytes(b"")
    (bad / "dataset.json").write_text("{}")
    p, argvs = main_env("--args", "--detector stalta", tuning=[bad])
    assert p.returncode != 0 and "not a fetched tuning set" in p.stderr and not argvs


def test_main_refuses_a_set_without_data(main_env):
    bad = main_env.tmp / "tune_ko_2"
    bad.mkdir()
    (bad / "dataset.json").write_text("{}")   # no miniSEED
    p, argvs = main_env("--args", "--detector stalta", tuning=[bad])
    assert p.returncode != 0 and "not a fetched tuning set" in p.stderr and not argvs
