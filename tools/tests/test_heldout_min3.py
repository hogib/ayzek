"""tools/heldout_min3.py (exploratory): the settings of its full stage (every
frozen detector with `--min-stations 3` appended) and the set names it takes
from a selection. No stage is run.

Expected values: the rule as the tool's docstring states it; a made-up settings
file; tools/heldout_settings.json read as it is (each entry compared with
itself plus the two appended arguments).
"""
import json

import heldout_min3 as hm
from conftest import REPO

FROZEN = {
    "_comment": ["ignored"],
    "learned": {"detector": "learned", "args": []},
    "learned_neighbours": {"detector": "learned", "grid": {"threshold": ["0.75", "0.85"]}},
    "stalta_station": {"detector": "stalta", "args": ["--stalta-band", "1,10"], "on": "8",
                       "ratios": {"gr": "tools/stalta_stations/stalta_on_gr.csv"}},
    "phasenet_causal_eqfa": {"detector": "learned", "model": "phasenet", "mode": "causal",
                             "args": ["--threshold", "0.8", "--min-stations", "2"]},
    "transformer": {"detector": "transformer", "args": []},
}


def test_every_detector_gets_min_stations_3_appended_and_keeps_everything_else():
    d = hm.min3_settings(FROZEN, "tools/heldout_settings.json", "abc123")
    assert list(d) == ["_comment", "learned", "stalta_station", "phasenet_causal_eqfa", "transformer"]
    assert d["learned"] == {"detector": "learned", "args": ["--min-stations", "3"]}
    assert d["stalta_station"] == {"detector": "stalta", "args": ["--stalta-band", "1,10", "--min-stations", "3"],
                                   "on": "8", "ratios": {"gr": "tools/stalta_stations/stalta_on_gr.csv"}}
    # appended at the end, also where the option already stands
    assert d["phasenet_causal_eqfa"]["args"] == ["--threshold", "0.8", "--min-stations", "2", "--min-stations", "3"]
    assert d["transformer"]["args"] == ["--min-stations", "3"]


def test_the_neighbour_grid_and_the_comments_are_left_out_and_the_source_is_named():
    d = hm.min3_settings(FROZEN, "tools/heldout_settings.json", "abc123")
    assert "learned_neighbours" not in d
    (c,) = d["_comment"]
    assert "tools/heldout_settings.json (sha256 abc123)" in c
    assert "Not a frozen setting" in c and "tools/heldout_min3.py" in c and "ignored" not in c


def test_only_keeps_the_named_labels():
    d = hm.min3_settings(FROZEN, "s.json", "x", only=["transformer", "learned", "learned_neighbours"])
    assert list(d) == ["_comment", "learned", "transformer"]


def test_the_frozen_settings_are_not_modified():
    before = json.dumps(FROZEN)
    d = hm.min3_settings(FROZEN, "s.json", "x")
    d["learned"]["args"].append("y")
    assert json.dumps(FROZEN) == before


def test_on_the_frozen_settings_file():
    frozen = json.loads((REPO / "tools/heldout_settings.json").read_text())
    d = hm.min3_settings(frozen, "tools/heldout_settings.json", "x")
    labels = [k for k, v in frozen.items() if not k.startswith("_") and "grid" not in v]
    assert list(d)[1:] == labels and len(labels) >= 5
    for k in labels:
        assert d[k]["args"] == frozen[k]["args"] + ["--min-stations", "3"]
        assert {x: y for x, y in d[k].items() if x != "args"} == {x: y for x, y in frozen[k].items() if x != "args"}


def test_the_changed_option_and_the_shifts():
    assert hm.MIN_STATIONS == "3"
    assert hm.SHIFTS == [0, 127, -127, 263, -263, 391, -391]


def test_sets_of_names_the_sets_as_run_heldout_does(tmp_path):
    (tmp_path / "selection.json").write_text(json.dumps({
        "group": "gr", "name": "now",
        "sequences": [{"index": 1, "selected": False}, {"index": 3, "selected": True}, {"index": 12, "selected": True}],
        "quiet": {"windows": [{}, {}, {}]}}))
    assert hm.sets_of(tmp_path) == ("gr", ["gr_now_seq003", "gr_now_seq012"],
                                    ["gr_now_quiet01", "gr_now_quiet02", "gr_now_quiet03"])
