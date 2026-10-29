"""tools/quiet_without_stations.py: the ranking of stations by their detections,
the exclusion, and the rates with and without them.

Expected values: hand arithmetic on recordings the test writes, with a fake
network_subsets whose alarm count for a subset is the number of recorded
detections at its stations (a rule the test chooses; the real tool's network
stage is tested in C++).
"""
import re

import pytest

import quiet_without_stations as qw
from conftest import run_tool, write_fake_binary

# The real header (tools/network_subsets.cpp): 22 columns.
HEADER = ("k,stations,catalogue,declared,located,unmatched,mean_alert_s,median_location_error_km,"
          "mean_abs_magnitude_error,main_declared,main_alert_s,main_blind_zone_km,"
          "main_first_magnitude,main_first_magnitude_s,main_final_magnitude,main_location_error_km,"
          "located_within_30km,located_beyond_50km,assessed_real_catalogued,"
          "assessed_real_unmatched,assessed_misfire_catalogued,assessed_misfire_unmatched")

FAKE_SUBSETS = f"""
import itertools, os, sys
argv = sys.argv[1:]
if os.environ.get("FAKE_EXIT"):
    sys.exit(int(os.environ["FAKE_EXIT"]))
rec = argv[0]
assert "--quiet" in argv
k = int(argv[argv.index("--min-k") + 1])
assert argv[argv.index("--max-k") + 1] == str(k)
stations, dets = [], []
for line in open(rec).read().splitlines():
    f = line.split("\\t")
    if f[0] == "S":
        stations.append(f[1])
    elif f[0] == "D":
        dets.append(f[1])
print({HEADER!r})
for sub in itertools.combinations(stations, k):
    n = sum(d in sub for d in dets)
    # Station names in reverse order, so that a match must be by set.
    print(f"{{k}},{{'+'.join(reversed(sub))}},0,{{n}},0,{{n}}" + ",-" * 16)
"""


def rec(stations, detections):
    """A recording with the given station codes and one D line per detection."""
    s = "# ayzek recording v1\n# min_stations 2\n"
    s += "".join(f"S\t{c}\t39.{i}\t28.{i}\n" for i, c in enumerate(stations))
    s += "T\t0\t100\n"
    s += "".join(f"D\t{c}\t{10 * i}\t{10 * i + 5}\t0.9\t0.1\t0\t{10 * i + 3.5}\n"
                 for i, c in enumerate(detections))
    return s


def test_stations_of(tmp_path):
    f = tmp_path / "w.rec"
    f.write_text(rec(["KULA", "SIMA", "DKL"], ["SIMA"]))
    assert qw.stations_of(f) == ["KULA", "SIMA", "DKL"]


def test_alarms_matches_the_subset_by_set(tmp_path):
    fake = write_fake_binary(tmp_path / "network_subsets", FAKE_SUBSETS)
    f = tmp_path / "w.rec"
    f.write_text(rec(["A", "B", "C"], ["A", "A", "B", "C", "A"]))
    # Detections: A 3, B 1, C 1. The fake prints the stations reversed.
    assert qw.alarms(fake, f, ["A", "B", "C"]) == 5
    assert qw.alarms(fake, f, ["C", "A"]) == 4
    assert qw.alarms(fake, f, ["B", "C"]) == 2
    with pytest.raises(SystemExit):
        qw.alarms(fake, f, ["A", "Z"])


def test_alarms_stops_when_network_subsets_fails(tmp_path, monkeypatch):
    fake = write_fake_binary(tmp_path / "network_subsets", FAKE_SUBSETS)
    f = tmp_path / "w.rec"
    f.write_text(rec(["A", "B"], []))
    monkeypatch.setenv("FAKE_EXIT", "2")
    with pytest.raises(SystemExit):
        qw.alarms(fake, f, ["A", "B"])


@pytest.fixture
def env(tmp_path):
    # Windows w1 (5 h) and w2 (7 h): 12 h.
    index = tmp_path / "quiet_index.csv"
    index.write_text("name,scored_hours\nw1,5.0\nw2,7.0\n")
    # stalta: w1 stations A B C with detections A 5, B 2, C 1; w2 stations
    # A B D E with A 3, D 2, B 1, E 1. Totals: A 8, B 3, D 2, C 1, E 1.
    st = tmp_path / "stalta"
    st.mkdir()
    (st / "w1.rec").write_text(rec(["A", "B", "C"], ["A"] * 5 + ["B"] * 2 + ["C"]))
    (st / "w2.rec").write_text(rec(["A", "B", "D", "E"], ["A"] * 3 + ["D"] * 2 + ["B", "E"]))
    # learned: w1 A 1, C 1; w2 E 2.
    le = tmp_path / "learned"
    le.mkdir()
    (le / "w1.rec").write_text(rec(["A", "B", "C"], ["A", "C"]))
    (le / "w2.rec").write_text(rec(["A", "B", "D", "E"], ["E", "E"]))
    fake = write_fake_binary(tmp_path / "network_subsets", FAKE_SUBSETS)

    def run(*args):
        p = run_tool("quiet_without_stations", "--index", index, "--binary", fake,
                     "--runs", f"learned={le}", f"stalta={st}", *args)
        rates = {m.group(1): tuple(float(x) for x in m.groups()[1:]) for m in re.finditer(
            r"^\s+(\S+)\s+all stations\s+(\d+) alarms,\s+([\d.]+) per day; without them\s+(\d+),"
            r"\s+([\d.]+) per day$", p.stdout, re.M)}
        return p, rates

    return run


def test_exclude_top_two_by_stalta(env):
    p, rates = env("--exclude-top", "2", "--by", "stalta")
    assert p.returncode == 0, p.stderr
    assert "station detections of stalta: A 8, B 3, D 2, C 1, E 1" in p.stdout
    assert "left out: A, B; 2 windows, 12.0 h" in p.stdout
    # stalta: all stations 8 + 7 = 15 alarms, 24*15/12 = 30 per day; without
    # A and B: w1 keeps C alone (fewer than 2 stations: 0), w2 keeps D, E:
    # 2 + 1 = 3 alarms, 24*3/12 = 6 per day.
    assert rates["stalta"] == (15, 30.0, 3, 6.0)
    # learned: 2 + 2 = 4 alarms, 8 per day; without A, B: 0 + 2 = 2, 4 per day.
    assert rates["learned"] == (4, 8.0, 2, 4.0)


def test_explicit_exclusion_adds_to_the_top(env):
    p, rates = env("--exclude", "C,D", "--exclude-top", "1", "--by", "stalta")
    assert p.returncode == 0, p.stderr
    assert "left out: A, C, D; 2 windows, 12.0 h" in p.stdout
    # stalta without A, C, D: w1 keeps B (alone: 0); w2 keeps B, E: 1 + 1 = 2.
    assert rates["stalta"] == (15, 30.0, 2, 4.0)
    # learned: w2 keeps B, E: E 2.
    assert rates["learned"] == (4, 8.0, 2, 4.0)


def test_explicit_exclusion_only(env):
    p, rates = env("--exclude", "E")
    assert p.returncode == 0, p.stderr
    assert "station detections of" not in p.stdout
    # stalta without E: w1 all three kept: 8; w2 keeps A, B, D: 6 -> 14, 28/day.
    assert rates["stalta"] == (15, 30.0, 14, 28.0)
    assert rates["learned"] == (4, 8.0, 2, 4.0)


def test_exclude_top_needs_by(env):
    p, _ = env("--exclude-top", "2")
    assert p.returncode != 0 and "needs --by" in p.stderr
    p, _ = env("--exclude-top", "2", "--by", "other")
    assert p.returncode != 0 and "needs --by" in p.stderr
