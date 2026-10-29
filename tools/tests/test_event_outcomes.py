"""tools/event_outcomes.py: per-event outcomes (detected, detected_pre_s) of a
report, at unshifted times, and the table main() writes.

Expected values: the fixtures' summaries (counts read off by hand), the
catalogue tests/catalogs/demo_ko.csv with ayzek's replay window, the pre-S
rule of the held-out protocol (section 7.5) as ayzek's summary applies it
(src/pipeline/network.cpp:1293, `used > 0 && before == used`), hand arithmetic.
"""
import calendar
import csv
import json

import pytest

import event_outcomes as eo
from conftest import REPO, fixture_text, run_tool, write_fake_binary, write_mseed_stub


def hms(t):
    """Seconds of the day -> 'HH:MM:SS.ff', by hand."""
    t %= 86400
    h, r = divmod(round(t * 100), 360000)
    m, r = divmod(r, 6000)
    return f"{h:02d}:{m:02d}:{r // 100:02d}.{r % 100:02d}"


def catalogue_day(path, day):
    """(seconds of day, type, magnitude as printed with one decimal) of the
    catalogue's events on DD/MM/YYYY `day`, parsed by hand."""
    out = []
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        f = line.split(",")
        d, clock = f[0].split(" ")
        if d != day:
            continue
        hh, mm, ss = (int(x) for x in clock.split(":"))
        out.append((hh * 3600 + mm * 60 + ss, f[5], f"{float(f[6]):.1f}"))
    return out


# --- real reports ---------------------------------------------------------------

# demo_ko learned (new build): detected (each "pre-S 2 of 2"): 18:07:17 ML 1.4,
# 18:16:19 ML 1.9, 18:20:51 MW 4.9, 18:23:21 ML 3.1, 18:25:51 ML 3.6,
# 18:29:11 ML 2.4, 18:30:02 ML 2.7, 18:32:40 ML 2.4; missed 18:10:29 ML 1.0.
# The STA/LTA report matches the same eight and misses the same one.
DETECTED = [("18:07:17.00", "ML", "1.4"), ("18:16:19.00", "ML", "1.9"),
            ("18:20:51.00", "MW", "4.9"), ("18:23:21.00", "ML", "3.1"),
            ("18:25:51.00", "ML", "3.6"), ("18:29:11.00", "ML", "2.4"),
            ("18:30:02.00", "ML", "2.7"), ("18:32:40.00", "ML", "2.4")]
MISSED = [("18:10:29.00", "ML", "1.0")]


@pytest.mark.parametrize("name", ["report_new_learned_demo_ko.txt",
                                  "report_new_stalta3_demo_ko.txt"])
def test_new_build_reports(name):
    # Summary of both: AFAD events within 250 km 9, detected 8, detected
    # before S 8, missed 1.
    ev = eo.outcomes(fixture_text(name), 0)
    assert ev == {**{k: (1, 1) for k in DETECTED}, **{k: (0, 0) for k in MISSED}}
    assert (len(ev), sum(d for d, _ in ev.values()), sum(p for _, p in ev.values())) == (9, 8, 8)


def test_old_build_report_has_no_pre_s():
    # The frozen build (9c78e4e) prints no pre-S line: detected_pre_s is ""
    # for a detected event (unknown); a missed event is 0 either way.
    ev = eo.outcomes(fixture_text("report_frozen_learned_demo_ko.txt"), 0)
    assert ev == {**{k: (1, "") for k in DETECTED}, **{k: (0, 0) for k in MISSED}}


def expected_shifted(shift):
    """The events a run against the catalogue shifted by `shift` scores, at
    their unshifted times: ayzek loads the catalogue events whose (shifted)
    time lies in [start - 30 s, end - 30 s] of the replay (app/ayzek.cpp:524),
    the replay being 18:04:02.73 to 18:35:00.01 (the report's header), and all
    events of that hour are within 20 km of the network (the reports' missed
    lists). Derived from the catalogue file, not from the tool."""
    lo, hi = 18 * 3600 + 4 * 60 + 2.73 - 30, 18 * 3600 + 35 * 60 + 0.01 - 30
    return {(hms(t), kind, mag): (0, 0)
            for t, kind, mag in catalogue_day(REPO / "tests/catalogs/demo_ko.csv", "10/11/2025")
            if lo <= t + shift <= hi}


@pytest.mark.parametrize("name, shift", [("report_new_learned_demo_ko_p127.txt", 127),
                                         ("report_new_learned_demo_ko_m263.txt", -263)])
def test_shifted_reports_unshift_to_catalogue_events(name, shift):
    want = expected_shifted(shift)
    # By hand: +127 adds 18:02:17 ML 1.4 (shifted to 18:04:24) and loses
    # 18:32:40 (to 18:34:47, after 18:34:30.01): 9 events; -263 loses 18:07:17
    # (to 18:02:54) and adds 18:34:51 ML 2.6 and 18:37:39 ML 1.8: 10 events.
    # Both reports: detected 0, so every outcome is (0, 0).
    if shift == 127:
        assert len(want) == 9 and ("18:02:17.00", "ML", "1.4") in want
        assert ("18:32:40.00", "ML", "2.4") not in want
    else:
        assert len(want) == 10 and ("18:37:39.00", "ML", "1.8") in want
        assert ("18:07:17.00", "ML", "1.4") not in want
    assert eo.outcomes(fixture_text(name), shift) == want


def test_unshifted_report_matches_catalogue_window():
    # The same derivation at shift 0 gives the 9 events of the unshifted report.
    assert set(eo.outcomes(fixture_text("report_new_learned_demo_ko.txt"), 0)) == set(expected_shifted(0))


# --- synthetic reports ---------------------------------------------------------------

def block(n, alarm, kind, mag, origin, pre_s=None):
    s = (f"\nEvent #{n}: detected earthquake (no magnitude), not located\n"
         f"  alarm      {alarm} UTC by AA, BB\n"
         f"  AFAD       {kind} {mag} at {origin} UTC: alarm 9.0 s after origin\n")
    if pre_s:
        s += f"  pre-S      {pre_s[0]} of {pre_s[1]} alarm detections declared before the predicted S\n"
    return s


def summary(catalogue, missed=()):
    s = ("\n------------------------------------------------------------------------------\n"
         "SUMMARY\n"
         f"  AFAD events within 250 km             {catalogue:>4}\n")
    if missed:
        s += "\n  Missed AFAD events\n" + "".join(
            f"    {t} UTC  {k} {m}    20 km from the network centre\n" for t, k, m in missed)
    return s + ("\n  False alarms (no AFAD event; small uncatalogued earthquakes are possible)\n"
                "    10:40:00.00 UTC  #9  no magnitude, not located\n"
                "------------------------------------------------------------------------------\n")


def test_pre_s_rule():
    # Detected before S when every detection the match used was declared
    # before the predicted S: 2 of 3 -> 0, 3 of 3 -> 1, 0 of 2 -> 0, and
    # 0 of 0 -> 0 (ayzek's summary counts `used > 0 && before == used`).
    out = (block(1, "10:00:05.00", "ML", "2.0", "09:59:56.00", (2, 3))
           + block(2, "10:10:09.00", "ML", "2.5", "10:10:00.00", (3, 3))
           + block(3, "10:20:09.00", "ML", "1.5", "10:20:00.00", (0, 2))
           + block(4, "10:25:09.00", "ML", "1.1", "10:25:00.00", (0, 0))
           + summary(5, [("10:30:00.00", "ML", "1.2")]))
    assert eo.outcomes(out, 0) == {("09:59:56.00", "ML", "2.0"): (1, 0),
                                   ("10:10:00.00", "ML", "2.5"): (1, 1),
                                   ("10:20:00.00", "ML", "1.5"): (1, 0),
                                   ("10:25:00.00", "ML", "1.1"): (1, 0),
                                   ("10:30:00.00", "ML", "1.2"): (0, 0)}


def test_missed_list_at_the_end_of_the_report():
    # Missed list followed directly by the closing rule (no false alarms).
    out = (block(1, "10:00:05.00", "ML", "2.0", "09:59:56.00", (1, 1))
           + "SUMMARY\n  AFAD events within 250 km                2\n"
           + "\n  Missed AFAD events\n    10:30:00.00 UTC  ML 1.2    20 km from the network centre\n"
           + "------------------------------------------------------------------------------\n")
    assert eo.outcomes(out, 0) == {("09:59:56.00", "ML", "2.0"): (1, 1),
                                   ("10:30:00.00", "ML", "1.2"): (0, 0)}


@pytest.mark.parametrize("out", [
    # 2 parsable events, the report counts 3
    block(1, "10:00:05.00", "ML", "2.0", "09:59:56.00", (1, 1)) + summary(3, [("10:30:00.00", "ML", "1.2")]),
    # no catalogue summary at all (a24f610: no longer passes for zero events)
    block(1, "10:00:05.00", "ML", "2.0", "09:59:56.00", (1, 1)),
    "",
])
def test_count_mismatch_or_no_summary_stops(out):
    with pytest.raises(SystemExit):
        eo.outcomes(out, 0)


@pytest.mark.xfail(strict=True, raises=SystemExit, reason=(
    "DETECTED and MISSED (event_outcomes.py:44-46) take the magnitude as [\\d.]+, so "
    "an event of negative magnitude (ML -0.3, possible in GeoNet's catalogue) is not "
    "parsed and the report stops with a count mismatch. Latent: the tuning "
    "catalogues' smallest magnitude is 0.63."))
def test_negative_magnitude():
    out = (block(1, "10:00:05.00", "ML", "-0.3", "09:59:56.00", (2, 2))
           + summary(2, [("10:30:00.00", "ML", "-0.5")]))
    assert eo.outcomes(out, 0) == {("09:59:56.00", "ML", "-0.3"): (1, 1),
                                   ("10:30:00.00", "ML", "-0.5"): (0, 0)}


# --- unshift ----------------------------------------------------------------------------

@pytest.mark.parametrize("shifted, shift, origin", [
    ("18:22:58.00", 127, "18:20:51.00"),     # 18:20:51 + 2:07
    ("00:01:00.00", 127, "23:58:53.00"),     # across midnight, backwards
    ("23:58:00.00", -391, "00:04:31.00"),    # 23:58:00 + 6:31 = 24:04:31
    ("12:00:00.25", 263, "11:55:37.25"),     # 12:00:00.25 - 4:23
])
def test_unshift_by_hand(shifted, shift, origin):
    assert eo.unshift(shifted, shift) == origin


def test_unshift_inverts_the_shift():
    # Invariant: an origin shifted by S and unshifted by S is the origin.
    for t in range(0, 86400, 977):
        for s in (127, -127, 263, -263, 391, -391):
            assert eo.unshift(hms(t + s), s) == hms(t)


# --- main ---------------------------------------------------------------------------------

FAKE_AYZEK = """
import calendar, json, os, sys
argv = sys.argv[1:]
with open(os.environ["FAKE_LOG"], "a") as f:
    f.write(json.dumps(argv) + "\\n")
sc = json.load(open(os.environ["FAKE_SCENARIO"]))
if sc.get("exit"):
    sys.exit(sc["exit"])
events = []
for line in open(argv[argv.index("--catalog") + 1]).read().splitlines()[1:]:
    f = line.split(",")
    day, clock = f[0].split(" ")
    d, m, y = (int(x) for x in day.split("/"))
    hh, mm, ss = (int(x) for x in clock.split(":"))
    events.append((calendar.timegm((y, m, d, hh, mm, ss, 0, 0, 0)), f[5], f[6]))
shift = events[0][0] - sc["t0"]
fmt = lambda t: "%02d:%02d:%02d.00" % (t % 86400 // 3600, t % 3600 // 60, t % 60)
outcome = sc["outcomes"][str(shift)]
missed = []
for i, (t, kind, mag) in enumerate(events):
    o = outcome[i]
    if o is None:
        missed.append((t, kind, mag))
        continue
    print(f"\\nEvent #{i + 1}: detected earthquake (no magnitude), not located")
    print(f"  alarm      {fmt(t + 15)} UTC by AA, BB")
    print(f"  AFAD       {kind} {mag} at {fmt(t)} UTC: alarm 15.0 s after origin")
    print(f"  pre-S      {o[0]} of {o[1]} alarm detections declared before the predicted S")
print("SUMMARY")
print(f"  AFAD events within 250 km             {len(events):>4}")
if missed:
    print("\\n  Missed AFAD events")
    for t, kind, mag in missed:
        print(f"    {fmt(t)} UTC  {kind} {mag}    20 km from the network centre")
print("------------------------------------------------------------------------------")
"""


@pytest.fixture
def outcomes_env(tmp_path):
    t0 = calendar.timegm((2025, 11, 10, 18, 10, 0, 0, 0, 0))
    write_mseed_stub(tmp_path / "S1.mseed", t0 - 360, t0 + 1500)
    cat = tmp_path / "cat.csv"
    cat.write_text("Date,Longitude,Latitude,Depth,Rms,Type,Magnitude,Location,EventID\n"
                   "10/11/2025 18:10:00,28.1,39.2,7.0,0.3,ML,2.0,a,1\n"
                   "10/11/2025 18:20:51,28.1,39.2,7.0,0.3,MW,4.9,b,2\n"
                   "10/11/2025 18:25:00,28.1,39.2,7.0,0.3,ML,1.5,c,3\n")
    fake = write_fake_binary(tmp_path / "ayzek", FAKE_AYZEK)
    log = tmp_path / "argv.jsonl"

    def run(scenario, *args):
        (tmp_path / "sc.json").write_text(json.dumps({"t0": t0, **scenario}))
        p = run_tool("event_outcomes", "--binary", fake, "--catalog", cat, "--jobs", 1,
                     "--out", tmp_path / "out.csv", *args,
                     env={"FAKE_LOG": str(log), "FAKE_SCENARIO": str(tmp_path / "sc.json")})
        rows = list(csv.reader(open(tmp_path / "out.csv"))) if p.returncode == 0 else None
        argvs = [json.loads(x) for x in log.read_text().splitlines()] if log.exists() else []
        return p, rows, argvs

    run.dataset = f"dk:CACHE:{tmp_path}/S1.mseed"
    return run


# Per event (in catalogue order): [before, used] when detected, None when missed.
SCENARIO = {"outcomes": {"0": [[2, 2], [1, 2], None],
                         "127": [[2, 2], None, None],
                         "-127": [None, None, None]}}


def test_main_table(outcomes_env):
    p, rows, argvs = outcomes_env(SCENARIO, "--dataset", outcomes_env.dataset,
                                  "--detector", "learned", "--shifts", "127")
    assert p.returncode == 0, p.stderr
    # By hand: shift 0 detects the first two (pre-S 2 of 2 -> 1, 1 of 2 -> 0)
    # and misses the third; +127 detects the first only; -127 none. Origins
    # are the unshifted catalogue times in every run.
    assert rows == [
        ["dataset", "detector", "shift", "origin", "type", "magnitude", "detected", "detected_pre_s"],
        ["dk", "learned", "0", "18:10:00.00", "ML", "2.0", "1", "1"],
        ["dk", "learned", "0", "18:20:51.00", "MW", "4.9", "1", "0"],
        ["dk", "learned", "0", "18:25:00.00", "ML", "1.5", "0", "0"],
        ["dk", "learned", "127", "18:10:00.00", "ML", "2.0", "1", "1"],
        ["dk", "learned", "127", "18:20:51.00", "MW", "4.9", "0", "0"],
        ["dk", "learned", "127", "18:25:00.00", "ML", "1.5", "0", "0"],
        ["dk", "learned", "-127", "18:10:00.00", "ML", "2.0", "0", "0"],
        ["dk", "learned", "-127", "18:20:51.00", "MW", "4.9", "0", "0"],
        ["dk", "learned", "-127", "18:25:00.00", "ML", "1.5", "0", "0"]]
    assert p.stdout.strip() == "dk learned: 2 of 3 detected, 1 before S; +127: 1/3, -127: 0/3"
    assert len(argvs) == 3
    for argv in argvs:
        assert argv[argv.index("--scores-in") + 1] == "CACHE"
        assert "--no-magnitude" in argv


def test_main_detector_options_and_label(outcomes_env):
    p, rows, argvs = outcomes_env(SCENARIO, "--dataset", outcomes_env.dataset,
                                  "--detector", "stalta_eqfa", "--shifts", "127",
                                  "--label", "stalta_x", "--", "--stalta-on", "6")
    assert p.returncode == 0, p.stderr
    assert {r[1] for r in rows[1:]} == {"stalta_x"}
    for argv in argvs:
        assert argv[argv.index("--detector") + 1] == "stalta"
        assert argv[argv.index("--stalta-on") + 1] == "6"
        # 943e031: the dataset's cache goes to every detector that has one
        # (the docstring: "The dataset's CACHE may be empty for any detector").
        assert argv[argv.index("--scores-in") + 1] == "CACHE"


def test_main_empty_cache_passes_no_scores(outcomes_env):
    # An empty CACHE field: no --scores-in, the learned detector running its
    # own networks (docstring of event_outcomes.py, 943e031).
    p, rows, argvs = outcomes_env(SCENARIO, "--dataset", outcomes_env.dataset.replace("CACHE", ""),
                                  "--detector", "learned", "--shifts", "127")
    assert p.returncode == 0, p.stderr
    assert len(argvs) == 3 and all("--scores-in" not in argv for argv in argvs)
    assert rows[1][:3] == ["dk", "learned", "0"]


def test_main_transformer_detector(outcomes_env):
    p, rows, argvs = outcomes_env(SCENARIO, "--dataset", outcomes_env.dataset,
                                  "--detector", "transformer", "--shifts", "127")
    assert p.returncode == 0, p.stderr
    assert {r[1] for r in rows[1:]} == {"transformer"}
    assert all(argv[argv.index("--detector") + 1] == "transformer" for argv in argvs)


def test_main_stops_when_ayzek_fails(outcomes_env):
    p, rows, _ = outcomes_env({"exit": 2}, "--dataset", outcomes_env.dataset,
                              "--detector", "learned", "--shifts", "127")
    assert p.returncode != 0 and "without a report" in p.stderr and rows is None
