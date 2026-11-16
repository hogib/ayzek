"""Three-station sensitivity of the held-out tests: every frozen detector of
tools/heldout_settings.json rerun with `--min-stations 3` in place of the frozen 2.

    python3 -u tools/heldout_min3.py [--stages counts,full] [--jobs 20] [--out data/runs/heldout_min3]

Exploratory, outside the held-out protocol (docs/protocol-heldout.md): the one
changed option is --min-stations; no setting is retuned or re-selected (the
baselines keep the settings selected on the development data, protocol section
6), and nothing under data/runs/heldout is written. The transformer already
declares on 3 stations by default, so the option changes nothing for it.

Stages (each in its own directory under OUT, written as NAME.partial and
renamed with a .done mark when complete; a rerun resumes a .partial):

  counts  the cheap route, minutes: the release recordings (RELEASE/NET/runs and
          quiet, LABEL/SET.rec) replayed through the network stage alone with
          tools/network_subsets at the recorded --min-stations (control) and at
          3, against the set's catalogue and its copies shifted by +-127, 263,
          391 s (tools/sweep_trigger.py's write_catalogue), and --quiet for the
          quiet windows. Station processing does not depend on the other
          stations (recording.hpp), so this equals rerunning ayzek. It gives
          counts only: scored events, detected, unmatched alarms, quiet alarms;
          the control is compared with the release reports. The chance-corrected
          count is tools/sweep_trigger.py's, per set, summed.
          -> OUT/counts/counts_from_recordings.csv, summary.csv
  full    the per-event route, for the paired statistics: tools/run_heldout.py
          itself, stages runs,outcomes,analysis, with a copy of the settings in
          which every detector's `args` end with `--min-stations 3`
          (OUT/settings_min3.json; the learned detector's neighbour grid is
          left out, endpoint 4 is not rerun). The learned detector replays its
          cached window probabilities (--scores-in; sequence sets from
          RELEASE/NET/scores, quiet windows from data/runs/heldout_quiet_scores,
          linked into OUT/NET/scores/learned), the published models theirs
          (data/runs/heldout_scores); STA/LTA, FilterPicker and the transformer
          run from the waveforms. -> OUT/NET/{runs,quiet,outcomes,analysis,
          commands.log,manifest.json} as a release run, and OUT/summary_reports.csv
          (the reports' own counts, summed per detector).

One progress line per finished unit with elapsed and ETA; exit status nonzero
when a stage fails or writes nothing. OUT/provenance.json: argv, parameters,
commit, versions, input hashes. Standard library only (runs with the server's
python3; the tools it calls are standard library only too).
"""
import argparse
import csv
import json
import shutil
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from heldout_explore_common import (NETS, ROOT, Progress, guard_out, read_recording, say, sha256,  # noqa: E402
                                    tree_digest, write_provenance)
from sweep_trigger import parse_report, write_catalogue  # noqa: E402

MIN_STATIONS = "3"
SHIFTS = [0, 127, -127, 263, -263, 391, -391]     # tools/run_heldout.py SHIFTS


def min3_settings(frozen, source, source_sha256, only=None):
    """The settings of the full stage: the frozen detectors of `frozen` (the
    parsed settings file `source`), each with `--min-stations 3` appended to
    its args; comment keys and the learned detector's neighbour grid are left
    out, and with `only` (a list of labels) every other detector too."""
    derived = {"_comment": [f"Exploratory copy of {source} (sha256 {source_sha256}): every detector's "
                            f"args end with --min-stations {MIN_STATIONS}; learned_neighbours left out. "
                            "Not a frozen setting; written by tools/heldout_min3.py."]}
    for k, v in frozen.items():
        if k.startswith("_") or "grid" in v or (only and k not in only):
            continue
        derived[k] = {**v, "args": [*v["args"], "--min-stations", MIN_STATIONS]}
    return derived


def sets_of(selection):
    """The sequence sets and quiet windows of a selection, named as
    tools/run_heldout.py names them."""
    s = json.loads((selection / "selection.json").read_text())
    tag = f"{s['group']}_{s['name']}"
    seqs = [f"{tag}_seq{q['index']:03d}" for q in s["sequences"] if q["selected"]]
    quiet = [f"{tag}_quiet{j:02d}" for j, _ in enumerate(s["quiet"]["windows"], 1)]
    return s["group"], seqs, quiet


def subsets(binary, rec, n, extra):
    """(catalogue, declared, unmatched) of the full station set of a recording."""
    p = subprocess.run([str(binary), str(rec), "--min-k", str(n), "--max-k", str(n), *extra],
                       capture_output=True, text=True)
    if p.returncode:
        if "no catalogue events in the recording's span" in p.stderr:
            return 0, 0, ""
        sys.exit(f"network_subsets failed on {rec} {' '.join(extra)}: {p.stderr[-400:]}")
    rows = list(csv.DictReader(p.stdout.splitlines()))
    if len(rows) != 1:
        sys.exit(f"{rec}: {len(rows)} subset rows, expected the full set alone")
    return int(rows[0]["catalogue"]), int(rows[0]["declared"]), int(rows[0]["unmatched"])


def stage_counts(a, out, plan, labels):
    final = out / "counts"
    if (final / ".done").exists():
        say(f"counts: {final} is complete, kept")
        return
    work = out / "counts.partial"
    work.mkdir(parents=True, exist_ok=True)
    rows, mismatches = [], []
    prog = Progress(len(plan) * len(labels), "detector x network (counts)")
    with tempfile.TemporaryDirectory(prefix="min3_") as tmp:
        def one(job):
            net, label, s, kind = job
            sub = "runs" if kind == "seq" else "quiet"
            rec = next((d / net / sub / label / f"{s}.rec" for d in a.recordings
                        if (d / net / sub / label / f"{s}.rec").exists()), None)
            if rec is None:
                return [], [f"{net} {label} {s}: no recording"]
            r = read_recording(rec)
            n, got, bad = len(r["stations"]), [], []
            for m in ("recorded", MIN_STATIONS):
                ms = [] if m == "recorded" else ["--min-stations", m]
                if kind == "quiet":
                    _, _, um = subsets(a.network_subsets, rec, n, ["--quiet", *ms])
                    got.append([net, label, s, kind, m, 0, "", "", um])
                    continue
                for shift in SHIFTS:
                    cat = Path(tmp) / f"{s}_{shift}.csv"
                    if not cat.exists():
                        t = Path(tmp) / f"{s}_{shift}.{label}.tmp"
                        write_catalogue(a.data / s / "catalog.csv", t, (r["t_first"], r["t_last"]), shift)
                        t.replace(cat)
                    c, d, um = subsets(a.network_subsets, rec, n, ["--catalog", str(cat), *ms])
                    got.append([net, label, s, kind, m, shift, c, d, um if shift == 0 else ""])
            # Control: the recorded setting against the release report of the same run.
            rep = rec.with_suffix(".txt")
            if rep.exists():
                pr = parse_report(rep.read_text())
                ctl = next(g for g in got if g[4] == "recorded" and g[5] == 0)
                if kind == "quiet":
                    if int(pr["alarms"] or 0) != ctl[8]:
                        bad.append(f"{net} {label} {s}: report {pr['alarms']} alarms, replay {ctl[8]}")
                elif (int(pr["catalogue"] or 0), int(pr["detected"] or 0), int(pr["false"] or 0)) != tuple(ctl[6:9]):
                    bad.append(f"{net} {label} {s}: report {pr['catalogue']}/{pr['detected']}/{pr['false']} "
                               f"(scored/detected/unmatched), replay {ctl[6]}/{ctl[7]}/{ctl[8]}")
            return got, bad

        with ThreadPoolExecutor(a.jobs) as pool:
            for net, (_, seqs, quiet) in plan.items():
                for label in labels:
                    jobs = [(net, label, s, "seq") for s in seqs] + [(net, label, s, "quiet") for s in quiet]
                    for got, bad in pool.map(one, jobs):
                        rows += got
                        mismatches += bad
                    prog.step(f"{net} {label}")
    if not rows:
        sys.exit("counts: nothing produced (no recordings found)")
    with open(work / "counts_from_recordings.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["network", "detector", "set", "kind", "min_stations", "shift", "scored", "detected", "unmatched"])
        w.writerows(rows)
    # Per network, detector and min-stations: the sums; corrected as sweep_trigger.py.
    agg = defaultdict(lambda: {"scored": 0, "detected": 0, "corrected": 0.0, "unmatched": 0, "quiet_alarms": 0,
                               "sets": 0, "windows": 0})
    by = defaultdict(dict)
    for net, label, s, kind, m, shift, c, d, um in rows:
        if kind == "quiet":
            g = agg[(net, label, m)]
            g["quiet_alarms"] += um
            g["windows"] += 1
        else:
            by[(net, label, m, s)][shift] = (c, d, um)
    for (net, label, m, s), sh in by.items():
        c, d, um = sh[0]
        g = agg[(net, label, m)]
        g["sets"] += 1
        g["scored"] += c
        g["detected"] += d
        g["unmatched"] += um or 0
        if c:
            chance = statistics.mean(sh[x][1] / max(1, sh[x][0]) for x in SHIFTS[1:])
            g["corrected"] += (d / c - chance) / (1 - chance) * c if chance < 1 else 0.0
    with open(work / "summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["network", "detector", "min_stations", "sequence_sets", "scored", "detected",
                    "corrected_sum_of_sets", "unmatched", "quiet_windows", "quiet_alarms"])
        for (net, label, m), g in agg.items():
            w.writerow([net, label, m, g["sets"], g["scored"], g["detected"], f"{g['corrected']:.1f}",
                        g["unmatched"], g["windows"], g["quiet_alarms"]])
    (work / "control_mismatches.txt").write_text("\n".join(mismatches) + ("\n" if mismatches else ""))
    say(f"counts: {len(rows)} rows; control against the release reports: {len(mismatches)} differences"
        + (f" (see control_mismatches.txt; first: {mismatches[0]})" if mismatches else ""))
    if final.exists():
        shutil.rmtree(final)
    work.rename(final)
    (final / ".done").touch()
    say(f"counts -> {final}")


def stage_full(a, out, plan, labels, settings_path):
    total = sum(len(labels) * (len(s) + len(q)) + len(labels) * len(s) + 1 for _, s, q in plan.values())
    prog = Progress(total, "runs+outcomes+analysis")
    for net, (group, seqs, quiet) in plan.items():
        final = out / net
        if (final / ".done").exists():
            say(f"full: {final} is complete, kept")
            prog.step(f"{net} kept", len(labels) * (2 * len(seqs) + len(quiet)) + 1)
            continue
        work = out / f"{net}.partial"
        # The learned detector's cached window probabilities, linked where
        # run_heldout looks for them (OUT/scores/learned/SET with .done).
        links = work / "scores" / "learned"
        links.mkdir(parents=True, exist_ok=True)
        missing = []
        for s, src in ([(s, a.release / net / "scores/learned" / s) for s in seqs]
                       + [(s, a.quiet_scores / net / "scores/learned" / s) for s in quiet]):
            if not (src / ".done").exists():
                missing.append(str(src))
            elif not (links / s).exists():
                (links / s).symlink_to(src.resolve(), target_is_directory=True)
        if missing and not a.allow_networks:
            sys.exit(f"{net}: no complete learned score cache for {len(missing)} sets (first: {missing[0]}); "
                     "--allow-networks lets run_heldout run the detector's networks for them")
        cmd = [sys.executable, "-u", str(ROOT / "tools/run_heldout.py"),
               "--selection", str(a.selections / net), "--data", str(a.data),
               "--binary", str(a.binary), "--network-subsets", str(a.network_subsets),
               "--settings", str(settings_path), "--scores", str(a.model_scores),
               "--stages", "runs,outcomes,analysis", "--jobs", str(a.jobs), "--out", str(work)]
        if group == "ko":
            cmd += ["--station-split", str(ROOT / "tools/detector_station_split.csv")]
        say("full: " + " ".join(cmd))
        p = subprocess.Popen(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in p.stdout:
            line = line.rstrip("\n")
            if line.startswith(("[runs] ", "[outcomes] ", "[analysis] ")):
                prog.step(f"{net} {line}")
            elif line:
                say(f"{net} | {line}")
        if p.wait():
            sys.exit(f"full: run_heldout exited {p.returncode} on {net}; {work} kept for a resume")
        if not (work / "analysis" / "significance.json").exists():
            sys.exit(f"full: {net} wrote no analysis/significance.json")
        work.rename(final)
        (final / ".done").touch()
        say(f"full: {net} -> {final} (commands.log and manifest.json name it as {work.name})")
    # The reports' own counts, summed per network and detector.
    rows = []
    for net, (_, seqs, quiet) in plan.items():
        for label in labels:
            tot = defaultdict(int)
            for s in seqs:
                r = parse_report((out / net / "runs" / label / f"{s}.txt").read_text())
                for k in ("catalogue", "detected", "alarms", "false"):
                    tot[k] += int(r[k] or 0)
            for s in quiet:
                tot["quiet"] += int(parse_report((out / net / "quiet" / label / f"{s}.txt").read_text())["alarms"] or 0)
            rows.append([net, label, MIN_STATIONS, len(seqs), tot["catalogue"], tot["detected"], tot["alarms"],
                         tot["false"], len(quiet), tot["quiet"]])
    if not rows:
        sys.exit("full: nothing produced")
    tmp = out / "summary_reports.csv.tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["network", "detector", "min_stations", "sequence_sets", "scored", "detected", "alarms",
                    "unmatched", "quiet_windows", "quiet_alarms"])
        w.writerows(rows)
    tmp.rename(out / "summary_reports.csv")
    # Cross-check with the recordings' replay, where that stage ran.
    counts = out / "counts" / "summary.csv"
    if counts.exists():
        ref = {(r["network"], r["detector"]): (r["scored"], r["detected"], r["unmatched"], r["quiet_alarms"])
               for r in csv.DictReader(open(counts)) if r["min_stations"] == MIN_STATIONS}
        diff = [f"{r[0]} {r[1]}" for r in rows
                if ref.get((r[0], r[1])) != (str(r[4]), str(r[5]), str(r[7]), str(r[9]))]
        say(f"full: reports against the recordings' replay at {MIN_STATIONS} stations: "
            + (f"{len(diff)} detector x network differ: {', '.join(diff[:8])}" if diff else "all equal"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(ROOT / "data/runs/heldout_min3"))
    ap.add_argument("--stages", default="counts,full")
    ap.add_argument("--nets", default=",".join(NETS))
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--selections", type=Path, default=ROOT.parent / "ayzek_heldout/selections")
    ap.add_argument("--data", type=Path, default=ROOT.parent / "ayzek_heldout/data")
    ap.add_argument("--binary", type=Path, default=ROOT / "build-linux/app/ayzek")
    ap.add_argument("--network-subsets", type=Path, default=ROOT / "build-linux/tools/network_subsets")
    ap.add_argument("--settings", type=Path, default=ROOT / "tools/heldout_settings.json")
    ap.add_argument("--release", type=Path, default=ROOT / "data/runs/heldout",
                    help="the release runs: NET/scores/learned (full) and NET/runs (counts)")
    ap.add_argument("--recordings", type=Path, nargs="+",
                    default=[ROOT / "data/runs/heldout", ROOT / "data/runs/heldout_f05436a"],
                    help="where NET/{runs,quiet}/LABEL/SET.rec are looked for, in order (counts)")
    ap.add_argument("--quiet-scores", type=Path, default=ROOT / "data/runs/heldout_quiet_scores")
    ap.add_argument("--model-scores", type=Path, default=ROOT / "data/runs/heldout_scores")
    ap.add_argument("--labels", help="only these detector labels (comma-separated; default: every detector)")
    ap.add_argument("--allow-networks", action="store_true",
                    help="run the learned detector's networks where a score cache is missing (slow)")
    a = ap.parse_args()
    out = guard_out(a.out)
    out.mkdir(parents=True, exist_ok=True)
    stages = a.stages.split(",")
    if set(stages) - {"counts", "full"}:
        sys.exit("--stages: counts, full")
    nets = a.nets.split(",")

    frozen = json.loads(a.settings.read_text())
    derived = min3_settings(frozen, a.settings, sha256(a.settings), a.labels.split(",") if a.labels else None)
    labels = [k for k in derived if not k.startswith("_")]
    settings_path = out / "settings_min3.json"
    text = json.dumps(derived, indent=1) + "\n"
    if not settings_path.exists() or settings_path.read_text() != text:
        settings_path.with_suffix(".json.tmp").write_text(text)
        settings_path.with_suffix(".json.tmp").rename(settings_path)

    plan = {net: sets_of(a.selections / net) for net in nets}
    for net, (_, seqs, quiet) in plan.items():
        for s in seqs + quiet:
            if not (a.data / s / "catalog.csv").exists():
                sys.exit(f"{a.data / s}: not fetched")
    say(f"{len(labels)} detectors, " + ", ".join(f"{n}: {len(s)} sequences + {len(q)} quiet windows"
                                                 for n, (_, s, q) in plan.items())
        + f"; stages {','.join(stages)}; {a.jobs} jobs -> {out}")

    inputs = {"settings": {"path": str(a.settings), "sha256": sha256(a.settings)},
              "settings_min3": {"path": str(settings_path), "sha256": sha256(settings_path)},
              "binary": {"path": str(a.binary), "sha256": sha256(a.binary)},
              "network_subsets": {"path": str(a.network_subsets), "sha256": sha256(a.network_subsets)},
              "tools": {t: sha256(ROOT / "tools" / t) for t in
                        ("run_heldout.py", "event_outcomes.py", "sweep_trigger.py", "significance.py",
                         "quiet_without_stations.py", "detector_station_split.csv")},
              "selections": {n: sha256(a.selections / n / "selection.json") for n in nets},
              "datasets": {s: sha256(a.data / s / "dataset.json") for _, seqs, quiet in plan.values()
                           for s in seqs + quiet if (a.data / s / "dataset.json").exists()}}
    say("hashing the score caches and recordings (inputs) ...")
    if "full" in stages:
        inputs["model_scores"] = tree_digest(a.model_scores)
        inputs["learned_scores"] = {n: tree_digest(a.release / n / "scores/learned") for n in nets}
        inputs["quiet_scores"] = {n: tree_digest(a.quiet_scores / n / "scores/learned") for n in nets}
    if "counts" in stages:
        inputs["recordings"] = {f"{d.name}/{n}/{sub}": tree_digest(d / n / sub)
                                for d in a.recordings for n in nets for sub in ("runs", "quiet")
                                if (d / n / sub).is_dir()}
    params = {"min_stations": MIN_STATIONS, "shifts": SHIFTS, "stages": stages, "nets": nets, "jobs": a.jobs,
              "detectors": labels, "run_heldout_stages": "runs,outcomes,analysis",
              "paths": {k: str(v) for k, v in vars(a).items() if isinstance(v, Path)}}
    write_provenance(out / "provenance.json", __file__, params, inputs)

    if "counts" in stages:
        stage_counts(a, out, plan, labels)
    if "full" in stages:
        stage_full(a, out, plan, labels, settings_path)
    write_provenance(out / "provenance.json", __file__, {**params, "finished": True}, inputs)
    (out / "DONE").write_text(f"{','.join(stages)} complete\n")
    say("DONE")


if __name__ == "__main__":
    main()
