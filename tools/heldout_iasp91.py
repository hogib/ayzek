"""Catalogue matching of the held-out tests rescored with IASP91 travel times
for the catalogue events deeper than 40 km, from the release reports and
recordings: no detector is rerun.

    uv run --offline --project tools python -u tools/heldout_iasp91.py [--jobs 4] [--out data/runs/heldout_iasp91]

Exploratory, outside the held-out protocol (docs/protocol-heldout.md, which
keeps the half-space matcher: section 7, Depth). Nothing under data/runs/heldout
is written.

Why a wrapper: ayzek's matcher (Network::match, src/pipeline/network.cpp)
predicts P as hypot(D, depth) / 6 km/s with no option for another model, and
tools/event_outcomes.py only reads the detected/missed lists ayzek prints. So
the match is repeated here, in Python, on what the release runs left:

  alarms      RELEASE/NET/runs/LABEL/SET.txt: every alarm's time, the stations
              that raised it (`by` list), its final location;
  detections  RELEASE/NET/runs/LABEL/SET.rec: the stations' detections; the ones
              each alarm was raised by are found with the rule of
              heldout_explore_common.alarm_detections (the report does not
              print their window starts);
  catalogue   DATA/SET/catalog.csv, the events ayzek's reader takes, with origin
              in [max(first sample, --from) - 30 s, last sample - 30 s], and its
              copies shifted by +-127, 263, 391 s (protocol section 5).

The rule is Network::match and Network::score unchanged: a catalogue event fits
an alarm when the alarm's location agrees (epicentre within 30 km, origin
within 5 s; misfit |origin difference| / 100) or the predicted P at a station
that raised it is within 6 s of that detection's window start + 3.5 s (misfit
the smallest such difference); the alarm is credited to the event of smallest
misfit; an event scored (within 250 km of the stations' centroid) is detected
when an alarm is credited with it, the first such alarm counting; "detected
before S" as Network::pre_s. Two travel-time models, everything else equal:

  halfspace   P = hypot(D, depth) / 6, S = hypot(D, depth) / 3.5: ayzek's own.
              This is the control: every alarm's credit is compared with the
              one the report prints, and every event's outcome with the release
              outcome table (RELEASE/NET/outcomes, when present), per shift;
              the differences are counted in selfcheck.csv and listed.
  iasp91      the same, except that for catalogue events deeper than --deep-km
              (40) P and S are the first arrivals of obspy's TauP in iasp91
              (phases p, P, Pn, Pg, Pdiff / s, S, Sn, Sg, Sdiff; distance from
              ayzek's haversine on a 6371 km sphere). The location criterion is
              unchanged.

Output per network, OUT/NET (written as NET.partial, renamed and marked .done):
  outcomes_halfspace/, outcomes_iasp91/   LABEL_SET.csv in the format of
        tools/event_outcomes.py (all shifts), read by tools/significance.py;
  analysis_halfspace/, analysis_iasp91/   tools/significance.py as
        tools/run_heldout.py runs it: significance, significance_pre_s, each
        also with --max-depth 40 (_depth40);
  deep_events.csv   every scored event deeper than --deep-km, per detector: its
        outcome under either model;
  alarm_counts.csv  per detector and set: alarms, unmatched under either model;
  selfcheck.csv, selfcheck_differences.txt.
OUT/provenance.json: argv, parameters, commit, versions, input hashes. One
progress line per finished unit with elapsed and ETA; exit status nonzero when
nothing is produced, or when the control differs and --allow-differences is
not given (the outputs are then kept as NET.partial).
"""
import argparse
import csv
import json
import math
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from heldout_explore_common import (NETS, ROOT, Progress, alarm_detections, distance_km, guard_out, hms,  # noqa: E402
                                    km_to_deg, read_catalogue, read_recording, read_report, say, sha256,
                                    tree_digest, write_provenance)

SHIFTS = [0, 127, -127, 263, -263, 391, -391]       # tools/run_heldout.py SHIFTS
VP, VS, TOLERANCE_S, RADIUS_KM = 6.0, 3.5, 6.0, 250.0  # network.hpp / Network::match
EQUIVALENCE = {"phasenet_causal_eqfa": 0.10}        # tools/run_heldout.py
P_PHASES = ["p", "P", "Pn", "Pg", "Pdiff"]
S_PHASES = ["s", "S", "Sn", "Sg", "Sdiff"]


class TravelTimes:
    """(P, S) travel times in seconds from an event to a point `km` away."""

    def __init__(self, model, deep_km):
        self.model, self.deep_km, self.cache, self.fallbacks = model, deep_km, {}, 0
        self.taup = None
        if model == "iasp91":
            from obspy.taup import TauPyModel
            self.taup = TauPyModel("iasp91")

    def __call__(self, km, depth):
        if self.taup is None or not depth > self.deep_km:
            h = math.hypot(km, depth)
            return h / VP, h / VS
        key = (round(depth, 3), round(km, 3))
        if key not in self.cache:
            deg = km_to_deg(km)
            out = []
            for phases, v in ((P_PHASES, VP), (S_PHASES, VS)):
                arr = self.taup.get_travel_times(source_depth_in_km=depth, distance_in_degree=deg, phase_list=phases)
                if arr:
                    out.append(min(x.time for x in arr))
                else:                      # no such arrival: keep ayzek's own prediction, counted
                    self.fallbacks += 1
                    out.append(math.hypot(km, depth) / v)
            self.cache[key] = tuple(out)
        return self.cache[key]


def rematch(alarms, dets, events, stations, tt):
    """Network::match per alarm (index of the credited event or None) and
    Network::score per scored event: {event index: (detected, detected_pre_s,
    alarm number)}, and the alarms credited as no scored event's first."""
    lat0 = sum(v[0] for v in stations.values()) / len(stations)
    lon0 = sum(v[1] for v in stations.values()) / len(stations)
    pred = {}

    def times(i, st):
        if (i, st) not in pred:
            e = events[i]
            p, s = tt(distance_km(e["lat"], e["lon"], *stations[st]), e["depth"])
            pred[(i, st)] = (e["time"] + p, e["time"] + s)
        return pred[(i, st)]

    credit = []
    for a, d in zip(alarms, dets):
        best, best_miss = None, 1e18
        for i, e in enumerate(events):
            miss = 1e18
            loc = a["location"]
            if loc and distance_km(loc[0], loc[1], e["lat"], e["lon"]) <= 30.0 and abs(loc[2] - e["time"]) <= 5.0:
                miss = abs(loc[2] - e["time"]) / 100.0
            for st, (ws, _) in d.items():
                if st in stations:
                    miss = min(miss, abs(times(i, st)[0] - (ws + 3.5)))
            if miss <= TOLERANCE_S and miss < best_miss:
                best, best_miss = i, miss
        credit.append(best)
    scored, first = {}, {}
    for i, e in enumerate(events):
        if distance_km(lat0, lon0, e["lat"], e["lon"]) > RADIUS_KM:
            continue
        hit = next((k for k, c in enumerate(credit) if c == i), None)   # alarms are in number order
        if hit is None:
            scored[i] = (0, 0, "")
            continue
        first[hit] = i
        used = [(st, decl) for st, (_, decl) in dets[hit].items() if st in stations]
        before = sum(decl < times(i, st)[1] for st, decl in used)
        scored[i] = (1, int(len(used) > 0 and before == len(used)), alarms[hit]["number"])
    return first, scored, len(alarms) - len(first)


def event_key(e, shift):
    """The event's identity in the outcome tables: unshifted origin time of day,
    type, magnitude with one decimal (tools/event_outcomes.py)."""
    return hms(e["time"] - shift), e["type"], f"{e['mag']:.1f}"


def one(net, label, s, a, models):
    rec = read_recording(a.release / net / "runs" / label / f"{s}.rec")
    alarms = sorted(read_report(a.release / net / "runs" / label / f"{s}.txt", rec), key=lambda x: x["number"])
    dets, problems, used = [], [], set()
    for al in alarms:
        d, pr = alarm_detections(al, rec, used=used)
        dets.append(d)
        problems += [f"{net} {label} {s} alarm #{al['number']}: {x}" for x in pr]
    cat = read_catalogue(a.data / s / "catalog.csv")
    lo, hi = max(rec["t_first"], rec["from"]) - 30.0, rec["t_last"] - 30.0
    tables = {m: [] for m in models}
    deep, counts, diffs = [], {"alarms": len(alarms)}, list(problems)
    release = {}
    rel = a.release / net / "outcomes" / f"{label}_{s}.csv"
    if rel.exists():
        for r in csv.DictReader(open(rel)):
            release[(int(r["shift"]), r["origin"], r["type"], r["magnitude"])] = (r["detected"], r["detected_pre_s"])
    n_alarm_diff = n_event_diff = n_events = 0
    for shift in SHIFTS:
        events = [{**e, "time": e["time"] + shift} for e in cat if lo <= e["time"] + shift <= hi]
        res = {m: rematch(alarms, dets, events, rec["stations"], tt) for m, tt in models.items()}
        for m, (credit, scored, unmatched) in res.items():
            for i, (det, pre, _) in sorted(scored.items(), key=lambda x: event_key(events[x[0]], shift)):
                tables[m].append([s, label, shift, *event_key(events[i], shift), det, pre, events[i]["depth_text"]])
            if shift == 0:
                counts[f"unmatched_{m}"] = unmatched
        credit, scored, _ = res["halfspace"]
        if shift == 0:          # control 1: the report's AFAD line names the scored event an alarm is the
            for k, al in enumerate(alarms):     # first credited with (Network::report), nothing otherwise
                c = credit.get(k)
                mine = None if c is None else event_key(events[c], 0)
                theirs = None if al["credit"] is None else (al["credit"][0], al["credit"][1], al["credit"][2])
                if mine != theirs:
                    n_alarm_diff += 1
                    diffs.append(f"{net} {label} {s} alarm #{al['number']} {al['hms']}: report credits {theirs}, rematch {mine}")
            for i, e in enumerate(events):
                if i in scored and e["depth"] > a.deep_km:
                    deep.append([net, label, s, *event_key(e, 0), e["depth_text"], e["id"],
                                 f"{distance_km(*centroid(rec['stations']), e['lat'], e['lon']):.0f}",
                                 *scored[i][:2], *res["iasp91"][1][i][:2]] if "iasp91" in res else [])
        if release:             # control 2: each event's outcome against the release table
            for i, (det, pre, _) in scored.items():
                n_events += 1
                k = (shift, *event_key(events[i], shift))
                if release.get(k) != (str(det), str(pre)):
                    n_event_diff += 1
                    diffs.append(f"{net} {label} {s} shift {shift} event {k[1:]}: release {release.get(k)}, rematch {(det, pre)}")
            extra = len([k for k in release if k[0] == shift]) - len(scored)
            if extra:
                n_event_diff += abs(extra)
                diffs.append(f"{net} {label} {s} shift {shift}: release table has {extra:+d} events against the rematch")
    check = [net, label, s, len(alarms), n_alarm_diff, len(problems), n_events if release else "", n_event_diff if release else ""]
    return tables, [d for d in deep if d], counts, check, diffs


def centroid(stations):
    return (sum(v[0] for v in stations.values()) / len(stations), sum(v[1] for v in stations.values()) / len(stations))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(ROOT / "data/runs/heldout_iasp91"))
    ap.add_argument("--nets", default=",".join(NETS))
    ap.add_argument("--labels", help="only these detector labels (default: every detector of the settings)")
    ap.add_argument("--sets", help="only these sets (comma-separated; for a quick look)")
    ap.add_argument("--deep-km", type=float, default=40.0, help="IASP91 for events deeper than this (default 40)")
    ap.add_argument("--jobs", type=int, default=4, help="significance.py runs at a time")
    ap.add_argument("--stages", default="rematch,analysis")
    ap.add_argument("--control-only", action="store_true", help="half-space control alone: no obspy, no analysis")
    ap.add_argument("--allow-differences", action="store_true",
                    help="go on to the analysis although the half-space control differs from the release")
    ap.add_argument("--selections", type=Path, default=ROOT.parent / "ayzek_heldout/selections")
    ap.add_argument("--data", type=Path, default=ROOT.parent / "ayzek_heldout/data")
    ap.add_argument("--release", type=Path, default=ROOT / "data/runs/heldout")
    ap.add_argument("--settings", type=Path, default=ROOT / "tools/heldout_settings.json")
    a = ap.parse_args()
    out = guard_out(a.out)
    out.mkdir(parents=True, exist_ok=True)
    stages = ["rematch"] if a.control_only else a.stages.split(",")
    nets = a.nets.split(",")
    settings = json.loads(a.settings.read_text())
    labels = [k for k, v in settings.items() if not k.startswith("_") and "grid" not in v
              and (not a.labels or k in a.labels.split(","))]
    models = {"halfspace": TravelTimes("halfspace", a.deep_km)}
    versions = {}
    if not a.control_only:
        import obspy
        versions["obspy"] = obspy.__version__
        models["iasp91"] = TravelTimes("iasp91", a.deep_km)

    plan = {}
    for net in nets:
        sel = json.loads((a.selections / net / "selection.json").read_text())
        tag = f"{sel['group']}_{sel['name']}"
        seqs = [f"{tag}_seq{q['index']:03d}" for q in sel["sequences"] if q["selected"]]
        plan[net] = (sel["group"], [s for s in seqs if not a.sets or s in a.sets.split(",")])
    n_rematch = sum(len(labels) * len(s) for _, s in plan.values())
    n_analysis = 0 if "analysis" not in stages else 4 * len(models) * len(plan)
    params = {"deep_km": a.deep_km, "shifts": SHIFTS, "vp": VP, "vs": VS, "tolerance_s": TOLERANCE_S,
              "radius_km": RADIUS_KM, "p_phases": P_PHASES, "s_phases": S_PHASES, "taup_model": "iasp91",
              "labels": labels, "nets": nets, "sets": a.sets, "stages": stages, "control_only": a.control_only,
              "paths": {k: str(v) for k, v in vars(a).items() if isinstance(v, Path)}}
    say("hashing inputs ...")
    inputs = {"settings": sha256(a.settings),
              "tools": {t: sha256(ROOT / "tools" / t) for t in ("significance.py", "event_outcomes.py")},
              "selections": {n: sha256(a.selections / n / "selection.json") for n in nets},
              "catalogues": {s: sha256(a.data / s / "catalog.csv") for _, seqs in plan.values() for s in seqs},
              "release_runs": {n: tree_digest(a.release / n / "runs") for n in nets},
              "release_outcomes": {n: tree_digest(a.release / n / "outcomes") for n in nets
                                   if (a.release / n / "outcomes").is_dir()},
              "release_manifests": {n: sha256(a.release / n / "manifest.json") for n in nets
                                    if (a.release / n / "manifest.json").exists()}}
    write_provenance(out / "provenance.json", __file__, params, inputs, versions)
    say(f"{len(labels)} detectors, {sum(len(s) for _, s in plan.values())} sequence sets, models {list(models)}; "
        f"stages {','.join(stages)} -> {out}")
    prog = Progress(n_rematch + n_analysis, "units (detector x set rematches, then significance runs)")
    total_diff, produced = 0, 0

    for net, (group, seqs) in plan.items():
        final = out / net
        if (final / ".done").exists():
            say(f"{final} is complete, kept")
            prog.step(f"{net} kept", len(labels) * len(seqs) + (4 * len(models) if "analysis" in stages else 0))
            produced += 1
            continue
        work = out / f"{net}.partial"
        work.mkdir(exist_ok=True)
        for m in models:
            (work / f"outcomes_{m}").mkdir(exist_ok=True)
        deep_rows, count_rows, check_rows, diffs = [], [], [], []
        for label in labels:
            for s in seqs:
                tables, deep, counts, check, d = one(net, label, s, a, models)
                for m, rows in tables.items():
                    tmp = work / f"outcomes_{m}" / f"{label}_{s}.csv.tmp"
                    with open(tmp, "w", newline="") as fh:
                        w = csv.writer(fh)
                        w.writerow(["dataset", "detector", "shift", "origin", "type", "magnitude", "detected",
                                    "detected_pre_s", "depth_km"])
                        w.writerows(rows)
                    tmp.rename(tmp.with_suffix(""))
                    produced += 1
                deep_rows += deep
                count_rows.append([net, label, s, counts["alarms"], *(counts.get(f"unmatched_{m}", "") for m in models)])
                check_rows.append(check)
                diffs += d
                prog.step(f"{net} {label} {s}: {check[3]} alarms, control differences {check[4]} alarm credits"
                          f"{'' if check[7] == '' else f', {check[7]} of {check[6]} event outcomes'}")
        for name, head, rows in (
                ("deep_events.csv", ["network", "detector", "set", "origin", "type", "magnitude", "depth_km", "event_id",
                                     "km_from_centroid", "detected_halfspace", "pre_s_halfspace", "detected_iasp91",
                                     "pre_s_iasp91"], deep_rows),
                ("alarm_counts.csv", ["network", "detector", "set", "alarms", *(f"unmatched_{m}" for m in models)], count_rows),
                ("selfcheck.csv", ["network", "detector", "set", "alarms", "alarm_credit_differences",
                                   "alarms_with_association_problems", "events_compared", "event_outcome_differences"],
                 check_rows)):
            with open(work / name, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(head)
                w.writerows(rows)
        (work / "selfcheck_differences.txt").write_text("\n".join(diffs) + ("\n" if diffs else ""))
        nd = sum(int(c[4]) + int(c[7] or 0) for c in check_rows)
        total_diff += nd
        say(f"{net}: half-space control: {nd} differences from the release over {len(check_rows)} detector x set"
            + (f"; first: {diffs[0]}" if diffs else ""))
        if nd and not a.allow_differences and not a.control_only:
            sys.exit(f"{net}: the half-space control does not reproduce the release ({work}/selfcheck_differences.txt); "
                     "stopped before the analysis. --allow-differences goes on regardless.")
        if "analysis" in stages:
            base_groups = f"{group}={','.join(seqs)}"
            against = ",".join(l for l in labels if l != "learned")
            eq = [f"{k}={v}" for k, v in EQUIVALENCE.items() if k in labels]
            jobs = []
            for m in models:
                an = work / f"analysis_{m}"
                an.mkdir(exist_ok=True)
                base = [sys.executable, str(ROOT / "tools/significance.py"), "--outcomes", str(work / f"outcomes_{m}"),
                        "--groups", base_groups, "--against", against]
                for name, extra in (("significance", ["--equivalence", *eq] if eq else []),
                                    ("significance_pre_s", ["--outcome", "detected_pre_s"])):
                    for suffix, depth in (("", []), ("_depth40", ["--max-depth", "40"])):
                        jobs.append((m, f"{name}{suffix}", [*base, *extra, *depth], an))

            def run(job):
                m, name, cmd, an = job
                if (an / f"{name}.json").exists():
                    return m, name, 0, ""
                p = subprocess.run([*cmd, "--json", str(an / f"{name}.json.tmp")], capture_output=True, text=True, cwd=ROOT)
                if p.returncode == 0:
                    (an / f"{name}.txt").write_text(p.stdout)
                    (an / f"{name}.json.tmp").rename(an / f"{name}.json")
                return m, name, p.returncode, p.stderr[-600:]

            with open(work / "commands.log", "a") as log:
                for _, _, cmd, _ in jobs:
                    log.write(" ".join(cmd) + "\n")
            failed = []
            with ThreadPoolExecutor(a.jobs) as pool:
                for m, name, code, err in pool.map(run, jobs):
                    prog.step(f"{net} {m} {name}{'' if code == 0 else f' FAILED exit {code}: {err}'}")
                    failed += [name] if code else []
            if failed:
                sys.exit(f"{net}: significance.py failed for {failed}; {work} kept")
        if a.control_only or a.sets or a.labels:
            say(f"{net}: partial run (--control-only, --sets or --labels): left as {work}, not marked done")
            continue
        work.rename(final)
        (final / ".done").touch()
        say(f"{net} -> {final}")
    if not produced:
        sys.exit("nothing produced")
    write_provenance(out / "provenance.json", __file__,
                     {**params, "finished": True, "control_differences": total_diff,
                      "taup_fallbacks_to_halfspace": models["iasp91"].fallbacks if "iasp91" in models else None},
                     inputs, versions)
    say(f"DONE; half-space control differences in total: {total_diff}")


if __name__ == "__main__":
    main()
