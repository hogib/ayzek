"""Runs a held-out test end to end from a `heldout select` directory: every
frozen detector on every selected set, then every endpoint of the protocol
(docs/protocol-heldout.md, section 7), with a manifest of what ran.

    python3 tools/run_heldout.py --selection ../ayzek_heldout/selections/greece_now \\
        --data ../ayzek_heldout/data --binary build-linux/app/ayzek \\
        --settings tools/heldout_settings.json --out data/runs/heldout/greece_now

The sets are the selection's sequences and quiet windows, fetched with the
held-out tool (DATA/<set>/: STATION.mseed, stations.csv, pd_gain.csv,
catalog.csv, dataset.json); every run takes the set's own station table and
catalogue, and on Greece and New Zealand the set's own Pd gains (--pd-gains,
which zeroes the Pd station terms and drops the noise levels). KO's sets run
without --pd-gains: KO keeps the gains, station terms and noise levels of its
calibration in models/ (protocol, section 2), and a KO station without a gain
there is not measured for Pd. Stages, each skipped where its outputs exist
(--stages picks):

  scores     the published models' window probabilities on every set
             (tools/seisbench_scores.py, pinned weights), written to NAME.tmp
             and marked .done when complete, as the retuning did;
  runs       one full run per detector and set, magnitude on, recorded
             (--record), the quiet windows from 30 min after their start:
             the reports feed the magnitude, latency and quiet endpoints;
  outcomes   tools/event_outcomes.py per detector and sequence set, with the
             chance shifts: the per-event tables of the paired tests;
  neighbours tools/sweep_trigger.py over the learned detector's adjacent
             settings on the sequence sets (endpoint 4);
  analysis   tools/significance.py (primary and secondary comparisons,
             equivalence with causal PhaseNet, detection before S, quiet
             rates), tools/quiet_without_stations.py (the two noisiest
             stations left out), the magnitude and latency tables from the
             reports, and the station-split comparison when --station-split
             names the 6 s detector's split (KO): tools/detector_station_split.csv,
             station,split from the detector's training manifest
             (dataset_specdual_6s_catalog/manifest.csv, 183 station keys; KRBG
             and MDNY appear under two keys in two splits and count as
             training). A sequence set counts as "test-split stations only"
             when every station is a test-split station.

This script adds no statistic: every number comes from the tools named above,
run with the commands written to OUT/commands.log. OUT/manifest.json records
the binary's SHA-256, the settings, the selection, the tools' commit and the
versions. A setting marked TBD in the settings file stops the run before it
starts. Standard library only.
"""
import argparse
import csv
import datetime as dt
import hashlib
import itertools
import json
import re
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
QUIET_LEAD_S = 1800.0
SHIFTS = "127,263,391"
EQUIVALENCE = {"phasenet_causal_eqfa": 0.10}   # protocol 7.2
EXCLUDE_TOP, EXCLUDE_BY = 2, "stalta_eqfa"     # protocol 7.3
PRIMARY_AGAINST = "stalta_eqfa"                # protocol 7, primary
OWN_GAINS = {"ko"}       # groups that keep models/' Pd gains, terms and noise (protocol 2)
OUTCOME_DETECTOR = {"learned": "learned", "stalta": "stalta_eqfa", "filterpicker": "filterpicker",
                    "transformer": "transformer"}
EVENT = re.compile(
    r"^Event #\d+: detected (?:M([\d.]+) earthquake|earthquake \(no magnitude\)).*\n"
    r"  alarm\s+(\S+) UTC by (.*)\n"
    r"(?:  Pd magn\.\s+M([\d.]+) (?:at the alarm|([\d.]+) s after the alarm), M([\d.]+) .*\n)?"
    r"(?:  (?:location|distance|P loc\.) .*\n)*"
    r"  AFAD\s+(\w+) ([\d.]+) at (\S+) UTC: alarm ([-\d.]+) s after origin(.*)\n"
    r"(?:  pre-S\s+(\d+) of (\d+) alarm detections.*\n)?", re.M)
WARNING = re.compile(r"^\s+(\S+)\s+(\d+) km\s+\S+ UTC \(predicted\)\s+([+-][\d.]+) s\s+(before|after) S", re.M)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Runner:
    def __init__(self, a):
        self.a = a
        self.out = Path(a.out)
        self.out.mkdir(parents=True, exist_ok=True)
        self.log = open(self.out / "commands.log", "a")
        self.selection = json.loads((Path(a.selection) / "selection.json").read_text())
        self.group, self.name = self.selection["group"], self.selection["name"]
        tag = f"{self.group}_{self.name}"
        self.sequences = [f"{tag}_seq{q['index']:03d}" for q in self.selection["sequences"]
                          if q["selected"]]
        self.quiet = [(f"{tag}_quiet{j:02d}", w) for j, w in
                      enumerate(self.selection["quiet"]["windows"], 1)]
        self.settings = {k: v for k, v in json.loads(Path(a.settings).read_text()).items()
                         if not k.startswith("_")}
        self.detectors = {k: v for k, v in self.settings.items() if "grid" not in v}
        self.check_settings()
        self.data = Path(a.data)
        for s in self.sequences + [q for q, _ in self.quiet]:
            if not (self.data / s / "dataset.json").exists() or not (self.data / s / "catalog.csv").exists():
                sys.exit(f"{self.data / s}: not fetched (dataset.json and catalog.csv needed)")

    def check_settings(self):
        bad = []
        for k, v in self.settings.items():
            vals = [v.get("args"), v.get("weights"), *(v.get("ratios", {}).values())]
            if any(x == "TBD" for x in vals):
                bad.append(k)
            if "ratios" in v and self.group not in v["ratios"]:
                bad.append(f"{k} (no ratios for {self.group})")
        if bad:
            sys.exit("settings not final: " + ", ".join(bad))

    def run(self, cmd, log_to=None, check=True):
        self.log.write(" ".join(str(c) for c in cmd) + "\n")
        self.log.flush()
        p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
        if log_to:
            Path(log_to).write_text(p.stdout + ("\n--- stderr ---\n" + p.stderr if p.stderr else ""))
        if check and p.returncode:
            sys.exit(f"failed ({p.returncode}): {' '.join(str(c) for c in cmd[:8])} ...\n{p.stderr[-1500:]}")
        return p

    # ---- the sets ---------------------------------------------------------
    def files(self, s):
        return sorted(str(p) for p in (self.data / s).glob("*.mseed"))

    def set_args(self, s):
        d = self.data / s
        args = ["--stations", d / "stations.csv", "--catalog", d / "catalog.csv"]
        if (d / "pd_gain.csv").exists() and self.group not in OWN_GAINS:
            args += ["--pd-gains", d / "pd_gain.csv"]
        return args

    def quiet_from(self, w):
        return (dt.datetime.fromisoformat(w["start"]) + dt.timedelta(seconds=QUIET_LEAD_S)).isoformat()

    def write_quiet_index(self):
        q = self.out / "quiet"
        q.mkdir(exist_ok=True)
        with open(q / "quiet_index.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["name", "start_utc", "score_from_utc", "end_utc", "scored_hours", "stations"])
            for s, win in self.quiet:
                t0, t1 = dt.datetime.fromisoformat(win["start"]), dt.datetime.fromisoformat(win["end"])
                w.writerow([s, t0.isoformat() + "+00:00", self.quiet_from(win) + "+00:00",
                            t1.isoformat() + "+00:00",
                            f"{(t1 - t0).total_seconds() / 3600 - QUIET_LEAD_S / 3600:.4f}",
                            len(self.files(s))])
        return q / "quiet_index.csv"

    # ---- detector options -------------------------------------------------
    def cache(self, label, s):
        v = self.detectors[label]
        if "model" not in v:
            return None
        return Path(self.a.scores) / f"{v['model']}_{v['weights']}_{v['mode']}" / s

    def detector_args(self, label, s):
        v = self.detectors[label]
        args = []
        if v["detector"] != "learned":
            args += ["--detector", v["detector"]]
        args += list(v["args"])
        if "ratios" in v:
            args += ["--stalta-on", v["on"], "--stalta-on-stations", ROOT / v["ratios"][self.group]]
        c = self.cache(label, s)
        if c is not None:
            args += ["--scores-in", c]
        return args

    # ---- stages -----------------------------------------------------------
    def stage_scores(self):
        cands = sorted({(v["model"], v["weights"], v["mode"]) for v in self.detectors.values()
                        if "model" in v})
        for m, w, mode in cands:
            for s in self.sequences + [q for q, _ in self.quiet]:
                d = Path(self.a.scores) / f"{m}_{w}_{mode}" / s
                if (d / ".done").exists():
                    continue
                d.parent.mkdir(parents=True, exist_ok=True)
                tmp = d.with_name(d.name + ".tmp")
                subprocess.run(["rm", "-rf", str(tmp)])
                self.run(["uv", "run", "--quiet", "--project", TOOLS, "--with", "seisbench==0.12.6",
                          "--with", "scipy", "python", "-W", "ignore", TOOLS / "seisbench_scores.py",
                          "--model", m, "--weights", w, "--mode", mode, "--out", tmp, *self.files(s)],
                         log_to=d.parent / f"{s}.log")
                tmp.rename(d)
                (d / ".done").touch()
                print(f"[scores] {m} {w} {mode} {s}")

    def one_run(self, label, s, quiet_from=None):
        sub = "quiet" if quiet_from else "runs"
        d = self.out / sub / label
        d.mkdir(parents=True, exist_ok=True)
        rep, rec = d / f"{s}.txt", d / f"{s}.rec"
        if rep.exists() and rec.exists():
            return
        cmd = [self.a.binary, "--speed", "0", "--no-color", "--record", rec, *self.set_args(s),
               *self.detector_args(label, s)]
        if quiet_from:
            cmd += ["--from", quiet_from]
        p = self.run([*cmd, *self.files(s)], check=False)
        if p.returncode or "AFAD events within" not in p.stdout:
            sys.exit(f"{label} {s}: ayzek exited {p.returncode} without a report\n{p.stderr[-1500:]}")
        rep.write_text(p.stdout)
        print(f"[runs] {label} {s}")

    def stage_runs(self):
        jobs = [(l, s, None) for l in self.detectors for s in self.sequences]
        jobs += [(l, s, self.quiet_from(w)) for l in self.detectors for s, w in self.quiet]
        with ThreadPoolExecutor(self.a.jobs) as pool:
            list(pool.map(lambda j: self.one_run(*j), jobs))

    def stage_outcomes(self):
        d = self.out / "outcomes"
        d.mkdir(exist_ok=True)

        def one(label, s):
            out = d / f"{label}_{s}.csv"
            if out.exists():
                return
            v = self.detectors[label]
            c = self.cache(label, s)
            self.run([sys.executable, TOOLS / "event_outcomes.py", "--binary", self.a.binary,
                      "--catalog", self.data / s / "catalog.csv",
                      "--dataset", f"{s}:{c or ''}:{','.join(self.files(s))}",
                      "--detector", OUTCOME_DETECTOR[v["detector"]], "--label", label,
                      "--shifts", SHIFTS, "--jobs", "7", "--out", out, "--",
                      *[x for x in self.detector_args(label, s) if x != "--scores-in" and x != c],
                      *self.set_args(s)], log_to=d / f"{label}_{s}.log")
            print(f"[outcomes] {label} {s}")

        with ThreadPoolExecutor(max(1, self.a.jobs // 7)) as pool:
            list(pool.map(lambda j: one(*j), [(l, s) for l in self.detectors for s in self.sequences]))

    def stage_neighbours(self):
        g = self.settings.get("learned_neighbours")
        if not g:
            return
        d = self.out / "neighbours"
        d.mkdir(exist_ok=True)
        for s in self.sequences:
            out = d / f"sweep_{s}.csv"
            if out.exists():
                continue
            self.run([sys.executable, TOOLS / "sweep_trigger.py", "--binary", self.a.binary,
                      "--catalog", self.data / s / "catalog.csv",
                      "--dataset", f"{s}::{','.join(self.files(s))}",
                      "--args", " ".join(str(x) for x in self.set_args(s)),
                      *itertools.chain.from_iterable(("--set", f"{k}={' '.join(v)}") for k, v in g["grid"].items()),
                      "--shifts", SHIFTS, "--jobs", str(self.a.jobs), "--out", out],
                     log_to=d / f"sweep_{s}.log")
            print(f"[neighbours] {s}")

    # ---- analysis ----------------------------------------------------------
    def catalogue_depths(self, s):
        out = {}
        for r in csv.DictReader(open(self.data / s / "catalog.csv", encoding="utf-8-sig")):
            hms = r["Date"].split()[1][:8]
            out[(hms, r["Type"], r["Magnitude"])] = float(r["Depth"]) if r["Depth"] else None
        return out

    def magnitude_and_latency(self, an):
        """Per matched event and detector: Pd first and final against the catalogue
        magnitude (with its type and depth), alarm delay, pre-S; and per large
        event the stations' warning leads."""
        rows, leads = [], []
        for label in self.detectors:
            for s in self.sequences:
                rep = (self.out / "runs" / label / f"{s}.txt").read_text()
                depths = self.catalogue_depths(s)
                for m in EVENT.finditer(rep):
                    (final_head, alarm_hms, by, pd_first, pd_first_s, pd_final,
                     ctype, cmag, chms, delay, rest, pre_before, pre_used) = m.groups()
                    key = (chms[:8], ctype, cmag)
                    depth = depths.get(key)
                    if depth is None:  # the report prints magnitudes with one decimal
                        depth = next((v for (h, t, mg), v in depths.items()
                                      if h == chms[:8] and t.upper() == ctype.upper()
                                      and abs(float(mg) - float(cmag)) < 0.051), None)
                    rows.append({"detector": label, "set": s, "origin": chms, "type": ctype,
                                 "magnitude": cmag, "depth_km": depth, "alarm_delay_s": delay,
                                 "pd_first": pd_first or "", "pd_first_after_alarm_s": pd_first_s or ("0" if pd_first else ""),
                                 "pd_final": pd_final or "", "pre_s_before": pre_before or "",
                                 "pre_s_used": pre_used or "", "stations": by})
                for block in re.split(r"\nEvent #", rep)[1:]:
                    a = re.search(r"^  AFAD\s+(\w+) ([\d.]+) at (\S+) UTC", block, re.M)
                    if not a:
                        continue
                    for st, dist, lead, side in WARNING.findall(block):
                        leads.append({"detector": label, "set": s, "origin": a.group(3), "type": a.group(1),
                                      "magnitude": a.group(2), "station": st, "distance_km": dist,
                                      "lead_s": lead, "side": side})
        with open(an / "magnitude_latency.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["detector"])
            w.writeheader()
            w.writerows(rows)
        with open(an / "warning_leads.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(leads[0]) if leads else ["detector"])
            w.writeheader()
            w.writerows(leads)
        # Summary of the Pd endpoint (7.7): error by band and class, events no
        # deeper than 40 km apart from deeper ones, by catalogue magnitude type.
        bands = [(-10, 4), (4, 5.5), (5.5, 20)]
        lines = ["detector  depth  type   band        n  with_Pd  first:bias  first:MAE  final:bias  final:MAE  class_first  class_final"]
        for label in self.detectors:
            for deep in (False, True):
                for ctype in sorted({r["type"] for r in rows}):
                    for lo, hi in bands:
                        sel = [r for r in rows if r["detector"] == label and r["type"] == ctype
                               and lo <= float(r["magnitude"]) < hi
                               and (r["depth_km"] is not None and (r["depth_km"] > 40) == deep)]
                        if not sel:
                            continue
                        with_pd = [r for r in sel if r["pd_final"]]
                        def stats(key):
                            e = [float(r[key]) - float(r["magnitude"]) for r in with_pd]
                            return (f"{sum(e) / len(e):+.2f}", f"{sum(abs(x) for x in e) / len(e):.2f}") if e else ("-", "-")
                        def cls(key):
                            ok = sum((lo <= float(r[key]) < hi) for r in with_pd)
                            return f"{ok}/{len(with_pd)}" if with_pd else "-"
                        fb, fm = stats("pd_first")
                        lb, lm = stats("pd_final")
                        lines.append(f"{label:9s} {'>40km' if deep else '<=40km':6s} {ctype:5s} "
                                     f"{f'{lo:g}-{hi:g}' if hi < 20 else f'>={lo:g}':10s} {len(sel):3d} {len(with_pd):7d}  "
                                     f"{fb:>10} {fm:>9}  {lb:>10} {lm:>9}  {cls('pd_first'):>11} {cls('pd_final'):>11}")
        (an / "magnitude_summary.txt").write_text("\n".join(lines) + "\n")

    def stage_analysis(self):
        an = self.out / "analysis"
        an.mkdir(exist_ok=True)
        group_arg = f"{self.group}={','.join(self.sequences)}"
        against = [l for l in self.detectors if l != "learned"]
        base = [sys.executable, TOOLS / "significance.py", "--outcomes", self.out / "outcomes",
                "--groups", group_arg, "--against", ",".join(against)]
        eq = [f"{k}={v}" for k, v in EQUIVALENCE.items() if k in self.detectors]
        self.run([*base, *(["--equivalence", *eq] if eq else []), "--json", an / "significance.json"],
                 log_to=an / "significance.txt")
        self.run([*base, "--outcome", "detected_pre_s", "--json", an / "significance_pre_s.json"],
                 log_to=an / "significance_pre_s.txt")
        if self.quiet:
            index = self.write_quiet_index()
            runs = [f"{l}={self.out / 'quiet' / l}" for l in ["learned", *against]]
            self.run([*base, "--against", PRIMARY_AGAINST, "--quiet", index, "--quiet-runs", *runs,
                      "--json", an / "quiet.json"], log_to=an / "quiet.txt")
            self.run([sys.executable, TOOLS / "quiet_without_stations.py", "--index", index,
                      "--runs", *runs, "--exclude-top", str(EXCLUDE_TOP), "--by", EXCLUDE_BY,
                      "--binary", self.a.network_subsets], log_to=an / "quiet_without_stations.txt")
        self.magnitude_and_latency(an)
        if self.a.station_split:
            split = {r["station"]: r["split"] for r in csv.DictReader(open(self.a.station_split))}
            groups = defaultdict(list)
            for q in self.selection["sequences"]:
                if not q["selected"]:
                    continue
                s = f"{self.group}_{self.name}_seq{q['index']:03d}"
                codes = [c.split(".")[1] for c in q["stations"]]
                roles = {split.get(c, "unknown") for c in codes}
                groups["test-split stations only" if roles <= {"test"} else "with training-split stations"].append(s)
            garg = ";".join(f"{k}={','.join(v)}" for k, v in groups.items() if v)
            self.run([*base, "--groups", garg, "--against", PRIMARY_AGAINST,
                      "--json", an / "station_split.json"], log_to=an / "station_split.txt")
        print(f"[analysis] -> {an}")

    def manifest(self):
        tool_commit = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                                     capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "src", "app", "tools"],
                               capture_output=True, text=True).stdout.strip()
        m = {"selection": {"path": str(Path(self.a.selection).resolve()),
                           "sha256": sha256(Path(self.a.selection) / "selection.json"),
                           "tool": self.selection.get("tool")},
             "binary": {"path": str(Path(self.a.binary).resolve()), "sha256": sha256(self.a.binary)},
             "network_subsets": {"path": str(Path(self.a.network_subsets).resolve()),
                                 "sha256": sha256(self.a.network_subsets)} if Path(self.a.network_subsets).exists() else None,
             "settings": {"path": str(Path(self.a.settings).resolve()), "sha256": sha256(self.a.settings),
                          "content": self.settings},
             "ayzek_code": {"commit": tool_commit, "dirty": bool(dirty)},
             "sequences": self.sequences, "quiet_windows": [s for s, _ in self.quiet],
             "stages": self.a.stages, "run_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds")}
        (self.out / "manifest.json").write_text(json.dumps(m, indent=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--selection", required=True, help="a heldout select directory")
    ap.add_argument("--data", required=True, help="the fetched sets' parent directory")
    ap.add_argument("--binary", required=True, help="the frozen ayzek")
    ap.add_argument("--network-subsets", default=str(ROOT / "build-release/tools/network_subsets"))
    ap.add_argument("--settings", default=str(TOOLS / "heldout_settings.json"))
    ap.add_argument("--scores", default=str(ROOT / "data/runs/heldout_scores"))
    ap.add_argument("--station-split", help="CSV station,split (training/validation/test) of the 6 s detector; KO")
    ap.add_argument("--stages", default="scores,runs,outcomes,neighbours,analysis")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    a.stages = a.stages.split(",")
    r = Runner(a)
    r.manifest()
    for stage in a.stages:
        getattr(r, f"stage_{stage}")()
    r.manifest()
    print(f"-> {r.out}/manifest.json")


if __name__ == "__main__":
    main()
