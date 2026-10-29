"""tools/run_heldout.py: the pure parts of the runner, with no stage run: the
set names from a selection, each detector's ayzek options (per-station
ratios, model caches), the refusal of unfinished settings, the quiet index,
the report regexes on the fixtures and the magnitude/latency table built
from them.

Expected values: tools/heldout_settings.json read by hand, the protocol's
option names (docs/protocol-heldout.md section 4: --stalta-on-stations), the
fixture reports and tests/catalogs/demo_ko.csv read by hand, arithmetic.
"""
import argparse
import csv
import json
import re
import shutil

import pytest

import run_heldout as rh
from conftest import FIXTURES, REPO, fixture_text

SETTINGS = {
    "_comment": "ignored",
    "learned": {"detector": "learned", "args": []},
    "learned_neighbours": {"detector": "learned", "grid": {"threshold": ["0.75", "0.85"],
                                                           "trigger-windows": ["6", "10"]}},
    "stalta_eqfa": {"detector": "stalta", "args": ["--stalta-band", "1,10", "--stalta-on", "8"]},
    "stalta_station": {"detector": "stalta", "args": ["--stalta-band", "1,10", "--stalta-off", "1.5"],
                       "on": "8", "ratios": {"gr": "tools/stalta_stations/stalta_on_gr.csv",
                                             "ko": "tools/stalta_stations/stalta_on_ko.csv"}},
    "fp_eqfa": {"detector": "filterpicker", "args": ["--fp-thresholds", "20,10"]},
    "phasenet_causal_eqfa": {"detector": "learned", "model": "phasenet", "weights": "instance",
                             "mode": "causal", "args": ["--threshold", "0.8"]},
    "transformer": {"detector": "transformer", "args": []},
}

SELECTION = {
    "group": "gr", "name": "nm", "tool": {"commit": "abc"},
    "sequences": [{"index": 1, "selected": False, "stations": ["GR.AAA"]},
                  {"index": 3, "selected": True, "stations": ["GR.AAA", "GR.BBB"]},
                  {"index": 12, "selected": True, "stations": ["GR.CCC"]}],
    "quiet": {"windows": [{"start": "2026-08-01T03:00:00", "end": "2026-08-01T09:00:00"},
                          {"start": "2026-08-02T22:10:00", "end": "2026-08-03T01:40:00"}]},
}
SETS = ["gr_nm_seq003", "gr_nm_seq012", "gr_nm_quiet01", "gr_nm_quiet02"]


def make_runner(tmp_path, settings=SETTINGS, selection=SELECTION, sets=SETS, pd_gain=True):
    sel = tmp_path / "sel"
    sel.mkdir(exist_ok=True)
    (sel / "selection.json").write_text(json.dumps(selection))
    (tmp_path / "settings.json").write_text(json.dumps(settings))
    data = tmp_path / "data"
    for s in sets:
        d = data / s
        d.mkdir(parents=True, exist_ok=True)
        (d / "dataset.json").write_text("{}")
        (d / "catalog.csv").write_text("Date,Longitude,Latitude,Depth,Rms,Type,Magnitude,Location,EventID\n")
        (d / "AAA.mseed").write_bytes(b"")
        (d / "BBB.mseed").write_bytes(b"")
        if pd_gain:
            (d / "pd_gain.csv").write_text("station,gain\n")
    a = argparse.Namespace(out=str(tmp_path / "out"), selection=str(sel), settings=str(tmp_path / "settings.json"),
                           data=str(data), scores=str(tmp_path / "scores"), binary=str(tmp_path / "ayzek"),
                           network_subsets=str(tmp_path / "network_subsets"), station_split=None,
                           stages=["analysis"], jobs=1)
    return rh.Runner(a)


# --- the sets -------------------------------------------------------------------------------------

def test_set_names_from_the_selection(tmp_path):
    r = make_runner(tmp_path)
    # Selected sequences by their index, zero-padded to three; quiet windows
    # numbered from 1 in the selection's order, two digits.
    assert r.sequences == ["gr_nm_seq003", "gr_nm_seq012"]
    assert [q for q, _ in r.quiet] == ["gr_nm_quiet01", "gr_nm_quiet02"]
    assert (r.group, r.name) == ("gr", "nm")
    # learned_neighbours (a grid) is not a detector; the comment is not a setting.
    assert list(r.detectors) == ["learned", "stalta_eqfa", "stalta_station", "fp_eqfa",
                                 "phasenet_causal_eqfa", "transformer"]
    assert "_comment" not in r.settings and "learned_neighbours" in r.settings


def test_unfetched_set_is_refused(tmp_path):
    with pytest.raises(SystemExit, match="gr_nm_quiet02: not fetched"):
        make_runner(tmp_path, sets=SETS[:-1])


def test_set_args(tmp_path):
    r = make_runner(tmp_path)
    d = tmp_path / "data" / "gr_nm_seq003"
    assert [str(x) for x in r.set_args("gr_nm_seq003")] == [
        "--stations", str(d / "stations.csv"), "--catalog", str(d / "catalog.csv"),
        "--pd-gains", str(d / "pd_gain.csv")]
    assert [str(x) for x in r.files("gr_nm_seq003")] == [str(d / "AAA.mseed"), str(d / "BBB.mseed")]


def test_set_args_without_gains(tmp_path):
    r = make_runner(tmp_path, pd_gain=False)
    assert [str(x) for x in r.set_args("gr_nm_seq003")][-2:] == ["--catalog", str(tmp_path / "data/gr_nm_seq003/catalog.csv")]
    assert "--pd-gains" not in map(str, r.set_args("gr_nm_seq003"))


# --- detector options -----------------------------------------------------------------------------

def test_detector_args_per_detector(tmp_path):
    r = make_runner(tmp_path)
    s = "gr_nm_seq003"
    args = {k: [str(x) for x in r.detector_args(k, s)] for k in r.detectors}
    # The learned detector: its own networks, nothing added.
    assert args["learned"] == []
    # STA/LTA: --detector and the setting's options as written.
    assert args["stalta_eqfa"] == ["--detector", "stalta", "--stalta-band", "1,10", "--stalta-on", "8"]
    # The per-station baseline: the network-wide ratio and the group's table
    # of per-station ratios, from the repository root (protocol section 4).
    assert args["stalta_station"] == [
        "--detector", "stalta", "--stalta-band", "1,10", "--stalta-off", "1.5",
        "--stalta-on", "8", "--stalta-on-stations", str(REPO / "tools/stalta_stations/stalta_on_gr.csv")]
    assert args["fp_eqfa"] == ["--detector", "filterpicker", "--fp-thresholds", "20,10"]
    assert args["transformer"] == ["--detector", "transformer"]
    # A published model: the learned detector's trigger options with the
    # model's cache SCORES/MODEL_WEIGHTS_MODE/SET.
    assert args["phasenet_causal_eqfa"] == [
        "--threshold", "0.8", "--scores-in", str(tmp_path / "scores/phasenet_instance_causal" / s)]


def test_cache_only_for_published_models(tmp_path):
    r = make_runner(tmp_path)
    assert r.cache("learned", "gr_nm_seq003") is None
    assert r.cache("stalta_station", "gr_nm_seq003") is None
    assert str(r.cache("phasenet_causal_eqfa", "gr_nm_quiet01")) == str(tmp_path / "scores/phasenet_instance_causal/gr_nm_quiet01")


def test_ratios_follow_the_selection_group(tmp_path):
    r = make_runner(tmp_path, selection={**SELECTION, "group": "ko"},
                    sets=[s.replace("gr_", "ko_") for s in SETS])
    args = [str(x) for x in r.detector_args("stalta_station", "ko_nm_seq003")]
    assert args[-1] == str(REPO / "tools/stalta_stations/stalta_on_ko.csv")


# --- settings not final ---------------------------------------------------------------------------

@pytest.mark.parametrize("change, named", [
    ({"fp_eqfa": {"detector": "filterpicker", "args": "TBD"}}, "fp_eqfa"),
    ({"phasenet_causal_eqfa": {**SETTINGS["phasenet_causal_eqfa"], "weights": "TBD"}}, "phasenet_causal_eqfa"),
    ({"stalta_station": {**SETTINGS["stalta_station"], "ratios": {"gr": "TBD", "ko": "x.csv"}}}, "stalta_station"),
    ({"stalta_station": {**SETTINGS["stalta_station"], "ratios": {"ko": "x.csv"}}}, "stalta_station (no ratios for gr)"),
])
def test_tbd_setting_stops_before_anything_runs(tmp_path, change, named):
    with pytest.raises(SystemExit, match=re.escape(f"settings not final: {named}")):
        make_runner(tmp_path, settings={**SETTINGS, **change})
    assert not (tmp_path / "out" / "runs").exists()


def test_committed_settings_are_final_for_every_group(tmp_path):
    # The committed tools/heldout_settings.json names no TBD and has a
    # per-station table for each held-out group (greece, ko, nz), each present.
    committed = json.loads((REPO / "tools/heldout_settings.json").read_text())
    for group in ("greece", "ko", "nz"):
        make_runner(tmp_path, settings=committed, selection={**SELECTION, "group": group},
                    sets=[s.replace("gr_", f"{group}_") for s in SETS])
    for path in committed["stalta_station"]["ratios"].values():
        rows = list(csv.DictReader(open(REPO / path)))
        assert rows and set(rows[0]) == {"station", "ratio"}
        assert all(float(r["ratio"]) in {3, 4, 5, 6, 8, 10, 12, 15} for r in rows)


# --- the quiet windows ----------------------------------------------------------------------------

def test_quiet_from_and_index(tmp_path):
    r = make_runner(tmp_path)
    # Scored from 30 min after the window's start (the lead-in warms the
    # filters): 03:00 -> 03:30, 22:10 -> 22:40.
    assert r.quiet_from(SELECTION["quiet"]["windows"][0]) == "2026-08-01T03:30:00"
    assert r.quiet_from(SELECTION["quiet"]["windows"][1]) == "2026-08-02T22:40:00"
    index = r.write_quiet_index()
    assert str(index) == str(tmp_path / "out/quiet/quiet_index.csv")
    # Scored hours: 6 h - 0.5 h = 5.5; 3 h 30 min - 0.5 h = 3.0. Two files per set.
    assert index.read_text().splitlines() == [
        "name,start_utc,score_from_utc,end_utc,scored_hours,stations",
        "gr_nm_quiet01,2026-08-01T03:00:00+00:00,2026-08-01T03:30:00+00:00,2026-08-01T09:00:00+00:00,5.5000,2",
        "gr_nm_quiet02,2026-08-02T22:10:00+00:00,2026-08-02T22:40:00+00:00,2026-08-03T01:40:00+00:00,3.0000,2"]


# --- the report regexes on real reports ----------------------------------------------------------------

def test_event_regex_on_the_new_report():
    # 8 matched events (the fixture's summary: detected 8), each with a pre-S
    # line "2 of 2"; event 3 (the M4.9) is the only one with a Pd magnitude:
    # M4.7 3.5 s after the alarm, M4.8 final (fixture lines 103-109). The 13
    # groups: final M, alarm time, stations, Pd first, its delay, Pd final,
    # catalogue type, M, origin, alarm delay, the rest of the AFAD line, pre-S
    # before, pre-S used.
    ms = list(rh.EVENT.finditer(fixture_text("report_new_learned_demo_ko.txt")))
    assert len(ms) == 8 and all(len(m.groups()) == 13 for m in ms)
    assert [m.group(9) for m in ms] == ["18:07:17.00", "18:16:19.00", "18:20:51.00", "18:23:21.00",
                                        "18:25:51.00", "18:29:11.00", "18:30:02.00", "18:32:40.00"]
    assert {(m.group(12), m.group(13)) for m in ms} == {("2", "2")}
    assert ms[2].groups() == ("4.8", "18:21:09.00", "KULA, SIMA; later DKL, KCTX, TVSB, YAYO",
                              "4.7", "3.5", "4.8", "MW", "4.9", "18:20:51.00", "18.0",
                              ", Pd magnitude -0.1, epicentre 3.7 km off, origin +1.7 s", "2", "2")
    assert ms[3].group(3) == "KCTX, KULA; later TVSB"
    assert ms[7].groups()[:6] == (None, "18:32:52.50", "KCTX, SIMA; later KULA, YAYO", None, None, None)
    assert ms[7].group(11) == ""        # not located: nothing after the alarm delay


def test_event_regex_on_the_old_report():
    # The frozen build prints no pre-S line: 8 events, the two pre-S groups unset.
    ms = list(rh.EVENT.finditer(fixture_text("report_frozen_learned_demo_ko.txt")))
    assert len(ms) == 8
    assert {(m.group(12), m.group(13)) for m in ms} == {(None, None)}
    assert [m.group(9) for m in ms][2] == "18:20:51.00"


def test_warning_regex_counts_the_predicted_arrivals():
    # 18 "(predicted)" lines in the new report (grep -c), one of them after S:
    # SIMA 74 km -1.2 s for the 18:23:21 ML 3.1 event (fixture line 127).
    found = rh.WARNING.findall(fixture_text("report_new_learned_demo_ko.txt"))
    assert len(found) == 18
    assert [f for f in found if f[3] == "after"] == [("SIMA", "74", "-1.2", "after")]
    assert found[0] == ("TVSB", "110", "+14.9", "before")


# --- magnitude and latency table ------------------------------------------------------------------------

def test_magnitude_and_latency_from_the_fixture(tmp_path):
    r = make_runner(tmp_path, settings={"learned": SETTINGS["learned"], "old": SETTINGS["learned"]})
    for label, name in (("learned", "report_new_learned_demo_ko.txt"), ("old", "report_frozen_learned_demo_ko.txt")):
        d = tmp_path / "out/runs" / label
        d.mkdir(parents=True)
        for s in r.sequences:
            shutil.copy(FIXTURES / name, d / f"{s}.txt")
    for s in r.sequences:
        shutil.copy(REPO / "tests/catalogs/demo_ko.csv", tmp_path / "data" / s / "catalog.csv")
    an = tmp_path / "out/analysis"
    an.mkdir()
    r.magnitude_and_latency(an)

    rows = list(csv.DictReader(open(an / "magnitude_latency.csv")))
    assert len(rows) == 2 * 2 * 8     # 2 detectors x 2 sets x 8 matched events
    m49 = [x for x in rows if x["detector"] == "learned" and x["set"] == "gr_nm_seq003" and x["magnitude"] == "4.9"]
    assert m49 == [{"detector": "learned", "set": "gr_nm_seq003", "origin": "18:20:51.00", "type": "MW",
                    "magnitude": "4.9", "depth_km": "9.64", "alarm_delay_s": "18.0", "pd_first": "4.7",
                    "pd_first_after_alarm_s": "3.5", "pd_final": "4.8", "pre_s_before": "2", "pre_s_used": "2",
                    "stations": "KULA, SIMA; later DKL, KCTX, TVSB, YAYO"}]
    # Depths from the catalogue by origin, type and magnitude: 18:07:17 ML 1.4
    # is at 4.39 km; the old build's rows carry no pre-S.
    old = {(x["origin"], x["depth_km"], x["pre_s_before"], x["pre_s_used"]) for x in rows if x["detector"] == "old"}
    assert ("18:07:17.00", "4.39", "", "") in old and len(old) == 8
    assert all(x["pd_first"] == "" and x["pd_first_after_alarm_s"] == "" for x in rows if x["magnitude"] != "4.9")

    leads = list(csv.DictReader(open(an / "warning_leads.csv")))
    assert len(leads) == 2 * 2 * 18
    assert {(x["station"], x["lead_s"], x["side"]) for x in leads if x["magnitude"] == "3.1"} == {
        ("SIMA", "-1.2", "after"), ("YAYO", "+6.4", "before"), ("DKL", "+8.9", "before"), ("KCTX", "+9.8", "before")}

    # The Pd summary: per detector, 7 ML events below M4 (none with a Pd
    # magnitude: the report gives Pd only for the M4.9) and the MW 4.9 in the
    # 4-5.5 band, 2 sets -> n = 14 and 2; its errors first 4.7 - 4.9 = -0.20,
    # final 4.8 - 4.9 = -0.10, both in the band: 2/2.
    lines = (an / "magnitude_summary.txt").read_text().splitlines()
    assert lines[0].split() == ["detector", "depth", "type", "band", "n", "with_Pd", "first:bias", "first:MAE",
                                "final:bias", "final:MAE", "class_first", "class_final"]
    assert [l.split() for l in lines[1:]] == [
        ["learned", "<=40km", "ML", "-10-4", "14", "0", "-", "-", "-", "-", "-", "-"],
        ["learned", "<=40km", "MW", "4-5.5", "2", "2", "-0.20", "0.20", "-0.10", "0.10", "2/2", "2/2"],
        ["old", "<=40km", "ML", "-10-4", "14", "0", "-", "-", "-", "-", "-", "-"],
        ["old", "<=40km", "MW", "4-5.5", "2", "2", "-0.20", "0.20", "-0.10", "0.10", "2/2", "2/2"]]


NEGATIVE_REPORT = """\
ayzek 0.0 replaying 2 files

Event #1: detected M-0.2 earthquake at 39.000N 28.000E, origin 10:01:00.00 UTC
  alarm      10:00:58.50 UTC by PA, PB
  Pd magn.   M-0.3 at the alarm, M-0.2 ± 0.1 (M<4 100%, 4-5.5 0%, M>=5.5 0%) final from 2 stations (2 above noise)
  AFAD       ML -0.3 at 10:01:00.00 UTC: alarm -1.5 s after origin, Pd magnitude +0.1, epicentre 1.0 km off, origin +0.5 s
  pre-S      2 of 2 alarm detections declared before the predicted S
  warning               distance   S wave arrives            warning
             PC            40 km   10:01:12.00 UTC (predicted)   +13.5 s  before S

SUMMARY
  AFAD events within 250 km                1
"""


def test_event_regex_negative_magnitudes():
    # Magnitudes printed with {:.1f} (src/pipeline/network.cpp:1090, 1134):
    # GeoNet's ML < 0 and a Pd magnitude below zero keep their sign.
    m = rh.EVENT.search(NEGATIVE_REPORT)
    assert m is not None
    assert m.groups() == ("-0.2", "10:00:58.50", "PA, PB", "-0.3", None, "-0.2", "ML", "-0.3", "10:01:00.00",
                          "-1.5", ", Pd magnitude +0.1, epicentre 1.0 km off, origin +0.5 s", "2", "2")


def test_magnitude_and_latency_negative_magnitude_and_depth(tmp_path):
    # The table and the warning leads carry the signed magnitude; the depth
    # comes from the catalogue, here -1.2 km (above sea level, as GeoNet
    # reports some).
    r = make_runner(tmp_path, settings={"learned": SETTINGS["learned"]})
    for s in r.sequences:
        d = tmp_path / "out/runs/learned"
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{s}.txt").write_text(NEGATIVE_REPORT)
        (tmp_path / "data" / s / "catalog.csv").write_text(
            "Date,Longitude,Latitude,Depth,Rms,Type,Magnitude,Location,EventID\n"
            "01/08/2026 10:01:00,28.0,39.0,-1.2,0.3,ML,-0.3,x,1\n")
    an = tmp_path / "out/analysis"
    an.mkdir()
    r.magnitude_and_latency(an)
    rows = list(csv.DictReader(open(an / "magnitude_latency.csv")))
    assert [(x["magnitude"], x["depth_km"], x["pd_first"], x["pd_final"], x["alarm_delay_s"]) for x in rows] == [
        ("-0.3", "-1.2", "-0.3", "-0.2", "-1.5")] * 2
    leads = list(csv.DictReader(open(an / "warning_leads.csv")))
    assert [(x["magnitude"], x["station"], x["lead_s"]) for x in leads] == [("-0.3", "PC", "+13.5")] * 2
    # Pd errors: first -0.3 - (-0.3) = 0.00, final -0.2 - (-0.3) = +0.10.
    line = (an / "magnitude_summary.txt").read_text().splitlines()[1].split()
    assert line == ["learned", "<=40km", "ML", "-10-4", "2", "2", "+0.00", "0.00", "+0.10", "0.10", "2/2", "2/2"]


# --- events no deeper than 40 km (protocol 7, Depth) -------------------------------------------

def test_neighbour_settings_from_the_grid(tmp_path):
    # SETTINGS' grid: threshold 0.75/0.85 x trigger-windows 6/10, in the
    # grid's key order, every combination once (by hand).
    r = make_runner(tmp_path)
    assert r.neighbour_settings() == [
        ("nb_threshold0.75_trigger-windows6", ["--threshold", "0.75", "--trigger-windows", "6"]),
        ("nb_threshold0.75_trigger-windows10", ["--threshold", "0.75", "--trigger-windows", "10"]),
        ("nb_threshold0.85_trigger-windows6", ["--threshold", "0.85", "--trigger-windows", "6"]),
        ("nb_threshold0.85_trigger-windows10", ["--threshold", "0.85", "--trigger-windows", "10"])]


def recorded_analysis(tmp_path, monkeypatch):
    """stage_analysis with the commands recorded instead of run."""
    split = tmp_path / "split.csv"
    split.write_text("station,split\nAAA,test\nBBB,training\nCCC,test\n")
    r = make_runner(tmp_path)
    r.a.station_split = str(split)
    for s in r.sequences:
        (r.out / "outcomes").mkdir(parents=True, exist_ok=True)
        (r.out / "outcomes" / f"learned_{s}.csv").write_text("x\n")
    (r.out / "neighbours" / "outcomes").mkdir(parents=True)
    cmds = []
    monkeypatch.setattr(r, "run", lambda cmd, log_to=None, check=True: cmds.append(([str(c) for c in cmd], str(log_to))))
    monkeypatch.setattr(r, "magnitude_and_latency", lambda an: None)
    r.stage_analysis()
    return r, cmds


def json_name(cmd):
    return cmd[cmd.index("--json") + 1].rsplit("/", 1)[1]


def test_every_detection_endpoint_also_for_events_no_deeper_than_40_km(tmp_path, monkeypatch):
    r, cmds = recorded_analysis(tmp_path, monkeypatch)
    by_name = {json_name(c): (c, log) for c, log in cmds if "significance.py" in c[1]}
    # Endpoints 1, 2 (significance), 5 (significance_pre_s), 4 (neighbours)
    # and 6 (station_split): each twice, the second with --max-depth 40 and
    # _depth40 names, otherwise the same command.
    for name in ("significance", "significance_pre_s", "neighbours", "station_split"):
        full, log = by_name[f"{name}.json"]
        shallow, log40 = by_name[f"{name}_depth40.json"]
        assert "--max-depth" not in full and log.endswith(f"/{name}.txt")
        assert shallow[shallow.index("--max-depth") + 1] == "40"
        assert log40.endswith(f"/{name}_depth40.txt")
        strip = lambda c: [x for x in c if x not in ("--max-depth", "40") and not x.endswith(".json")]
        assert strip(full) == strip(shallow)
    # The quiet windows have no events: no depth version.
    assert "quiet.json" in by_name and "quiet_depth40.json" not in by_name
    assert all("--max-depth" not in c for c, _ in cmds if "quiet_without_stations.py" in c[1])
    # The neighbours' comparison: the frozen rule's tables copied beside the
    # neighbours', every grid setting against it.
    nb, _ = by_name["neighbours.json"]
    assert nb[nb.index("--outcomes") + 1] == str(r.out / "neighbours" / "outcomes")
    assert nb[nb.index("--against") + 1] == ",".join(l for l, _ in r.neighbour_settings())
    assert sorted(p.name for p in (r.out / "neighbours" / "outcomes").glob("learned_*.csv")) == [
        "learned_gr_nm_seq003.csv", "learned_gr_nm_seq012.csv"]
    assert len(by_name) == 9


def test_manifest_lists_the_endpoint_files(tmp_path):
    r = make_runner(tmp_path)
    an = r.out / "analysis"
    an.mkdir(parents=True)
    for n in ("significance.json", "significance.txt", "significance_depth40.json",
              "significance_depth40.txt", "quiet.txt", "magnitude_summary.txt"):
        (an / n).write_text("")
    assert r.endpoint_files() == {
        "all_events": ["analysis/magnitude_summary.txt", "analysis/quiet.txt",
                       "analysis/significance.json", "analysis/significance.txt"],
        "events_no_deeper_than_40_km": ["analysis/significance_depth40.json",
                                        "analysis/significance_depth40.txt"]}
    (tmp_path / "ayzek").write_bytes(b"binary")          # hashed by the manifest
    r.manifest()
    m = json.loads((r.out / "manifest.json").read_text())
    assert m["endpoint_files"] == r.endpoint_files()
