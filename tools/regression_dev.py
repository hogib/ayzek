"""The frozen settings on the development replays: does a binary reproduce the
retuning's numbers? (held-out protocol, section 1: the release binary must.)

    python3 tools/regression_dev.py --binary build-linux/app/ayzek --write-reference
    python3 tools/regression_dev.py --binary build-linux/app/ayzek [--only events|quiet]

For every detector of tools/heldout_settings.json, the chance-corrected
sweep of tools/sweep_trigger.py on the five development event replays (demo,
marmara, sindirgi61, simav26, demo_ko; the published models from their score
caches under data/runs/sb) and tools/quiet_sweep.py on the 23 AFAD quiet
windows, with --no-onset-fp --no-p-locate as the retuning ran them. The rows
(alarms, unmatched, detected, catalogue, chance, corrected, alarm delay,
triggers; quiet alarms) are compared with tests/reference/heldout_dev.csv,
which --write-reference assembles from the retuning's own tables
(data/runs/tune_dev); the comparison exits 1 on any difference. The
per-station STA/LTA has no development reference (its ratios are the test
networks') and is left out.

Standard library only; the development data are not in the repository.
"""
import argparse
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
REF = ROOT / "tests/reference/heldout_dev.csv"
SETS = {"demo": ("data/demo/DEMI.mseed,data/demo/MANT.mseed,data/demo/BAND.mseed", "data/runs/cache_sindirgi"),
        "marmara": ("data/marmara/*.mseed", "data/runs/cache_marmara"),
        "sindirgi61": ("data/sindirgi61/*.mseed", "data/runs/cache_sindirgi61"),
        "simav26": ("data/simav26/*.mseed", "data/runs/cache_simav26"),
        "demo_ko": ("data/demo_ko/*.mseed", "data/runs/cache_demo_ko")}
FIELDS = ["alarms", "false", "detected", "catalogue", "chance", "corrected", "alarm_delay_s", "triggers"]
EXTRA = "--no-onset-fp --no-p-locate"


def settings():
    d = {k: v for k, v in json.loads((TOOLS / "heldout_settings.json").read_text()).items()
         if not k.startswith("_") and "grid" not in v and "ratios" not in v}
    return d


def sweep_dir(label, v):
    """The retuning's table directory of this detector."""
    if "model" in v:
        return ROOT / f"data/runs/tune_dev/sb_{v['model']}_{v['weights']}_{v['mode']}"
    return ROOT / "data/runs/tune_dev" / {"learned": "learned", "stalta": "stalta", "filterpicker": "fp",
                                          "transformer": "transformer"}[v["detector"]]


def quiet_table(label, v):
    if "model" in v:
        return ROOT / f"data/runs/tune_dev/quiet/sb_{v['model']}_{v['weights']}_{v['mode']}.csv"
    name = {"learned": "learned", "stalta": "stalta", "filterpicker": "fp", "transformer": "transformer"}
    return ROOT / "data/runs/tune_dev/quiet" / f"{name[v['detector']]}.csv"


def setting_keys(v):
    """The sweep columns and values that identify this setting."""
    args = list(v["args"])
    if v["detector"] == "learned" and "model" not in v and not args:
        return {"threshold": "0.8", "trigger-windows": "8", "instant-threshold": "0.9"}
    out = {}
    for i in range(0, len(args), 2):
        if args[i].startswith("--") and i + 1 < len(args) and not args[i + 1].startswith("--"):
            out[args[i][2:]] = args[i + 1]
        else:
            out[args[i][2:]] = "on"
    if v["detector"] == "stalta":
        out.setdefault("stalta-3c", "off")
        out.pop("stalta-off", None)
    return out


def match(rows, keys):
    hit = [r for r in rows if all(r.get(k) == val for k, val in keys.items())]
    if len(hit) != 1:
        sys.exit(f"{len(hit)} rows match {keys}")
    return hit[0]


def write_reference(labels):
    REF.parent.mkdir(parents=True, exist_ok=True)
    out = []
    for label, v in labels.items():
        keys = setting_keys(v)
        for name in SETS:
            r = match(list(csv.DictReader(open(sweep_dir(label, v) / f"sweep_{name}.csv"))), keys)
            out.append({"detector": label, "dataset": name, **{f: r[f] for f in FIELDS}, "quiet_alarms": ""})
        q = match(list(csv.DictReader(open(quiet_table(label, v)))),
                  {k: val for k, val in keys.items() if k not in ("stalta-off",)})
        out.append({"detector": label, "dataset": "quiet100", **{f: "" for f in FIELDS}, "quiet_alarms": q["alarms"]})
    with open(REF, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["detector", "dataset", *FIELDS, "quiet_alarms"], lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(f"reference: {len(out)} rows -> {REF}")


def run_events(binary, label, v, jobs, tmp):
    rows = {}
    base = []
    if v["detector"] != "learned":
        base += ["--detector", v["detector"]]
    for name, (files, cache) in SETS.items():
        c = (ROOT / f"data/runs/sb/{v['model']}_{v['weights']}_{v['mode']}/{name}") if "model" in v \
            else (ROOT / cache if v["detector"] == "learned" else "")
        out = Path(tmp) / f"{label}_{name}.csv"
        cmd = [sys.executable, TOOLS / "sweep_trigger.py", "--binary", binary,
               "--catalog", ROOT / f"tests/catalogs/{name}.csv",
               "--dataset", f"{name}:{c}:{','.join(str(ROOT / f) for f in files.split(','))}",
               "--args", " ".join([*base, EXTRA]), "--shifts", "127,263,391", "--jobs", str(jobs), "--out", out]
        for k, val in setting_keys(v).items():
            cmd += ["--set", f"{k}={val}"]
        p = subprocess.run([str(x) for x in cmd], capture_output=True, text=True, cwd=ROOT)
        if p.returncode:
            sys.exit(f"{label} {name}: sweep failed\n{p.stderr[-1200:]}")
        rows[name] = list(csv.DictReader(open(out)))[0]
    return rows


def run_quiet(binary, label, v, jobs, tmp):
    out = Path(tmp) / f"{label}_quiet.csv"
    args = ([] if v["detector"] == "learned" else ["--detector", v["detector"]]) + [EXTRA]
    cmd = [sys.executable, TOOLS / "quiet_sweep.py", "--binary", binary, "--args", " ".join(args),
           "--jobs", str(jobs), "--out", out]
    if "model" in v:
        cmd += ["--scores-root", ROOT / f"data/runs/sb_quiet/{v['model']}_{v['weights']}_{v['mode']}"]
    elif v["detector"] == "learned":
        cmd += ["--scores-root", ROOT / "data/runs/quiet100_scores"]
    for k, val in setting_keys(v).items():
        cmd += ["--set", f"{k}={val}"]
    p = subprocess.run([str(x) for x in cmd], capture_output=True, text=True, cwd=ROOT)
    if p.returncode:
        sys.exit(f"{label} quiet: sweep failed\n{p.stderr[-1200:]}")
    return list(csv.DictReader(open(out)))[0]["alarms"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--binary", required=True)
    ap.add_argument("--write-reference", action="store_true")
    ap.add_argument("--only", choices=("events", "quiet"))
    ap.add_argument("--detectors", help="comma-separated labels (default: all)")
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    labels = settings()
    if a.detectors:
        labels = {k: labels[k] for k in a.detectors.split(",")}
    if a.write_reference:
        write_reference(labels)
        return
    ref = {(r["detector"], r["dataset"]): r for r in csv.DictReader(open(REF))}
    bad = 0
    with tempfile.TemporaryDirectory(prefix="regression_") as tmp:
        for label, v in labels.items():
            if a.only != "quiet":
                for name, row in run_events(a.binary, label, v, a.jobs, tmp).items():
                    want = ref[(label, name)]
                    diff = {f: (want[f], row[f]) for f in FIELDS if want[f] != row[f]}
                    print(f"{label:28s} {name:11s} {'ok' if not diff else 'DIFFERS ' + str(diff)}")
                    bad += bool(diff)
            if a.only != "events":
                got = run_quiet(a.binary, label, v, a.jobs, tmp)
                want = ref[(label, "quiet100")]["quiet_alarms"]
                print(f"{label:28s} {'quiet100':11s} {'ok' if got == want else f'DIFFERS (reference {want}, got {got})'}")
                bad += got != want
    print("all rows reproduced" if not bad else f"{bad} rows differ")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
