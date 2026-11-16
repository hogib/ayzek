"""Grids on the three held-out test networks: one parameter per detector varied
around its frozen setting, everything else as frozen.

    python3 -u tools/heldout_grid.py --grid tradeoff [--jobs N] [--stages runs,outcomes]
    python3 -u tools/heldout_grid.py --grid instant  [--jobs N] [--stages runs,outcomes]

EXPLORATORY, OUTSIDE THE HELD-OUT PROTOCOL: nothing is retuned and no frozen setting
changes; the frozen point of each detector is on its curve and must reproduce the
release-time counts. The grid goes through tools/run_heldout.py itself (stages runs
and outcomes), with a settings file in which each grid point is one detector entry,
so every command is the protocol's own. Two grids (--grid; the former drivers
logs/grid_run.py and logs/instant_run.py, the second a copy of the first with
another grid and output directory):

 tradeoff, default output data/runs/heldout_grid: the trade-off curves
  learned (3-seed 6 s detector)   --threshold 0.4 ... 0.9    (8 windows, instant 0.9: frozen)
  STA/LTA rule 1 (stalta_eqfa)    --stalta-on 4 ... 20       (band, STA, LTA, off: frozen)
  causal PhaseNet (eqfa rule)     --threshold 0.3 ... 0.9    (1 window: frozen)
  causal EQTransformer (eqfa)     --threshold 0.3 ... 0.9    (8 windows: frozen)
  FilterPicker (eqfa)             --fp-thresholds S1,10, S1 = 10 ... 40  (windows, S2: frozen)

 instant, default output data/runs/heldout_instant: the learned detector's
 single-window branch
  learned                         --instant-threshold 0.85 ... 0.905 around its frozen
                                  0.9, and 2 (switched off)  (0.8 x 8 windows: frozen)

Inputs: tools/heldout_settings.json (the frozen settings), the selections and data of
the held-out tool (../ayzek_heldout/selections/<net>, ../ayzek_heldout/data), the
binaries build-linux/app/ayzek and build-linux/tools/network_subsets. The learned
detector and the published models replay from cached window scores (--scores-in): the
sequence sets' from the release run, the quiet windows' from
data/runs/heldout_quiet_scores (2026-11-14); OUT/<net>/scores/learned links to both.

The rule for a grid point: a copy of the frozen entry of its detector in which the
value after the varied option is replaced, or the option and value appended when the
frozen entry does not carry it; label g_<family>_<value>, the value up to its first
comma with `.` written as `p` (settings_for below).

Output: OUT/grid_settings.json (the settings, rewritten at every start),
OUT/<net>/{runs,quiet,outcomes}; a provenance record OUT/grid.json (written as
grid.json.tmp, renamed at the end); DONE written at the end. One log line per finished
unit with count, elapsed and ETA. Exits nonzero when a run fails or nothing is
produced. Standard library only (runs with the server's python3).
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# family: (frozen entry it copies, the option varied, its values). The order of the
# families and of the values is the order of the settings file, and with it the order
# of tools/significance.py's --against list in tools/heldout_grid_summary.py.
GRIDS = {
    "tradeoff": {
        "learned": ("learned", "--threshold", ["0.4", "0.5", "0.6", "0.7", "0.75", "0.8", "0.85", "0.9"]),
        "stalta": ("stalta_eqfa", "--stalta-on", ["4", "5", "6", "7", "8", "10", "12", "15", "20"]),
        "phasenet": ("phasenet_causal_eqfa", "--threshold", ["0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"]),
        "eqt": ("eqt_causal_eqfa", "--threshold", ["0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"]),
        # added 2026-11-15, appended last so that the earlier settings keep their place in the
        # list (and with it their bootstrap intervals): FilterPicker's first threshold, the
        # second as frozen (10); the label carries the first threshold only
        "fp": ("fp_eqfa", "--fp-thresholds", ["10,10", "12,10", "15,10", "20,10", "25,10", "30,10", "40,10"]),
    },
    "instant": {
        # the single-window ("instant") branch of the learned detector's rule: frozen 0.9;
        # 2 switches it off (ayzek: > 1 disables); everything else as frozen (0.8 x 8 windows)
        "instant": ("learned", "--instant-threshold", ["0.85", "0.88", "0.895", "0.9", "0.903", "0.905", "2"]),
    },
}
DEFAULT_OUT = {"tradeoff": "data/runs/heldout_grid", "instant": "data/runs/heldout_instant"}


def settings_for(frozen, grid):
    """The settings file of grid `grid` (a key of GRIDS) from the frozen settings: a
    comment, then one detector entry per grid point, in grid order."""
    settings = {"_comment": [f"Exploratory trade-off grid, outside the protocol; written by tools/heldout_grid.py --grid {grid}"]}
    for name, (base, opt, values) in GRIDS[grid].items():
        for v in values:
            e = json.loads(json.dumps(frozen[base]))
            args = list(e["args"])
            if opt in args:
                args[args.index(opt) + 1] = v
            else:
                args += [opt, v]
            e["args"] = args
            settings[f"g_{name}_{v.split(',')[0].replace('.', 'p')}"] = e
    return settings


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def stamp():
    return time.strftime("%FT%TZ", time.gmtime())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--grid", required=True, choices=list(GRIDS))
    ap.add_argument("--out", help="output directory, relative to ayzek_code/ (default: data/runs/heldout_grid for "
                                  "tradeoff, data/runs/heldout_instant for instant)")
    ap.add_argument("--jobs", type=int, default=20)
    ap.add_argument("--stages", default="runs,outcomes")
    ap.add_argument("--nets", default="ko_now,nz_now,greece_now")
    a = ap.parse_args()
    os.chdir(ROOT)
    OUT = Path(a.out or DEFAULT_OUT[a.grid])
    OUT.mkdir(parents=True, exist_ok=True)
    GRID = GRIDS[a.grid]

    frozen = json.loads(Path("tools/heldout_settings.json").read_text())
    settings = settings_for(frozen, a.grid)
    sfile = OUT / "grid_settings.json"
    sfile.write_text(json.dumps(settings, indent=1))
    n_det = len(settings) - 1

    prov = {"argv": sys.argv, "exploratory": True,
            "grid": {k: {"base": b, "option": o, "values": v} for k, (b, o, v) in GRID.items()},
            "code_commit": Path("COMMIT").read_text().strip() if Path("COMMIT").exists() else
            subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
            "python": sys.version.split()[0], "binary_sha256": sha("build-linux/app/ayzek"),
            "frozen_settings_sha256": sha("tools/heldout_settings.json"), "grid_settings_sha256": sha(sfile),
            "started_utc": time.strftime("%FT%TZ", time.gmtime()), "nets": a.nets.split(",")}
    (OUT / "grid.json.tmp").write_text(json.dumps(prov, indent=1))

    plan = []
    for net in a.nets.split(","):
        sel = json.loads(Path(f"../ayzek_heldout/selections/{net}/selection.json").read_text())
        tag = f"{sel['group']}_{sel['name']}"
        seqs = [f"{tag}_seq{q['index']:03d}" for q in sel["sequences"] if q["selected"]]
        quiet = [f"{tag}_quiet{j:02d}" for j in range(1, len(sel["quiet"]["windows"]) + 1)]
        link = OUT / net / "scores" / "learned"
        link.mkdir(parents=True, exist_ok=True)
        for s, src in [(s, f"data/runs/heldout/{net}/scores/learned/{s}") for s in seqs] + \
                      [(s, f"data/runs/heldout_quiet_scores/{net}/scores/learned/{s}") for s in quiet]:
            if not (Path(src) / ".done").exists():
                sys.exit(f"{src}: no complete score cache")
            if not (link / s).exists():
                (link / s).symlink_to(Path(src).resolve())
        units = n_det * (len(seqs) + len(quiet)) * ("runs" in a.stages) + n_det * len(seqs) * ("outcomes" in a.stages)
        plan.append((net, units))
    total = sum(u for _, u in plan)
    print(f"{stamp()} {n_det} settings, {total} units on {a.nets}, {a.jobs} jobs -> {OUT}", flush=True)

    start, done = time.time(), 0
    for net, units in plan:
        cmd = [sys.executable, "-u", "tools/run_heldout.py", "--selection", f"../ayzek_heldout/selections/{net}",
               "--data", "../ayzek_heldout/data", "--binary", "build-linux/app/ayzek",
               "--network-subsets", "build-linux/tools/network_subsets", "--settings", str(sfile),
               "--scores", "data/runs/heldout_scores", "--stages", a.stages, "--jobs", str(a.jobs),
               "--out", str(OUT / net)]
        before = done
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in p.stdout:
            line = line.rstrip()
            if line.startswith("[runs]") or line.startswith("[outcomes]"):
                done += 1
                el = time.time() - start
                print(f"{stamp()} {done}/{total} {line} | elapsed {el / 60:.0f} min, "
                      f"ETA {el / max(done, 1) * (total - done) / 60:.0f} min", flush=True)
            else:
                print(f"{stamp()} {line}", flush=True)
        if p.wait() != 0:
            sys.exit(f"{stamp()} {net}: run_heldout exited {p.returncode}; FAILED")
        done = before + units          # units skipped because their outputs existed
        print(f"{stamp()} {net} finished", flush=True)

    n_out = sum(1 for _ in OUT.glob("*/outcomes/*.csv")) + sum(1 for _ in OUT.glob("*/runs/*/*.txt"))
    if n_out == 0:
        sys.exit(f"{stamp()} nothing produced; FAILED")
    prov["finished_utc"] = stamp()
    prov["files"] = n_out
    (OUT / "grid.json.tmp").write_text(json.dumps(prov, indent=1))
    (OUT / "grid.json.tmp").rename(OUT / "grid.json")
    (OUT / "DONE").touch()
    print(f"{stamp()} DONE {n_out} files", flush=True)


if __name__ == "__main__":
    main()
