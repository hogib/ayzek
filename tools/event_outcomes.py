"""Per-event outcomes of a detector on a replay, for paired statistics.

    python3 tools/event_outcomes.py --catalog tests/catalogs/sindirgi61.csv \\
        --dataset 'sindirgi61:data/runs/cache_sindirgi61:data/sindirgi61/*.mseed' \\
        --detector learned --shifts 127,263,391 --out outcomes.csv

Runs ayzek once against the catalogue and once against each copy shifted by
+S and -S seconds, exactly as tools/sweep_trigger.py does for its chance level,
and writes one row per catalogue event and run: detected or missed. The event
is identified by its unshifted origin time, so the rows of one event line up
across shifts and across detectors, which is what a paired test needs.

--detector is `learned` (the cached window probabilities, default trigger
rule), `stalta_eqfa` (--detector stalta), `filterpicker` or `transformer`;
a setting goes after `--` (anything there is passed to ayzek), with a --label.
`stalta_eqfa` alone is ayzek's default STA/LTA setting; the frozen settings of
the held-out tests are in tools/heldout_settings.json, which
tools/run_heldout.py passes after `--`. (The old `stalta_eqdet`, the equal
detections setting of docs/impl/09, was removed before the freeze: it was not
the frozen one; tools/reproduce_all.sh passes it after `--`.)
Another model's scores (tools/seisbench_scores.py) run as `learned` with their
own cache, trigger settings after `--`, and a --label. The dataset's CACHE may
be empty for any detector, the learned one then running its own networks.
Standard library only.
"""
import argparse
import concurrent.futures
import csv
import datetime as dt
import glob
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from sweep_trigger import data_span, write_catalogue  # noqa: E402

DETECTORS = {
    "learned": [],
    "stalta_eqfa": ["--detector", "stalta"],
    "filterpicker": ["--detector", "filterpicker"],
    "transformer": ["--detector", "transformer"],
}
# Detected: the AFAD line of an event block, with the pre-S line after it
# (absent in reports of older builds). Missed: the list after the counts.
DETECTED = re.compile(r"^  AFAD\s+(\w+) (-?[\d.]+) at (\S+) UTC: alarm[^\n]*\n"
                      r"(?:  pre-S\s+(\d+) of (\d+) alarm detections)?", re.M)
MISSED = re.compile(r"^    (\d\d:\d\d:\d\d\.\d\d) UTC  (\w+) (-?[\d.]+) ", re.M)
CATALOGUE = re.compile(r"AFAD events within \d+ km\s+(\d+)")


def unshift(hms, shift):
    t = dt.datetime.strptime(hms, "%H:%M:%S.%f") - dt.timedelta(seconds=shift)
    return t.strftime("%H:%M:%S.%f")[:-4]


def outcomes(out, shift):
    """{(origin, type, M): (detected, detected before S)} for one report, at
    unshifted times. Detected before S: every detection the match used was
    declared before the predicted S at its station (held-out protocol, section
    7); "" when the report predates the pre-S line."""
    ev = {}
    for kind, mag, hms, before, used in DETECTED.findall(out):
        pre_s = "" if not used else int(int(used) > 0 and before == used)
        ev[(unshift(hms, shift), kind, mag)] = (1, pre_s)
    tail = out[out.find("Missed AFAD events"):] if "Missed AFAD events" in out else ""
    tail = tail[:tail.find("\n\n", 20)] if tail else ""
    for hms, kind, mag in MISSED.findall(tail):
        ev[(unshift(hms, shift), kind, mag)] = (0, 0)
    m = CATALOGUE.search(out)
    if not m or int(m.group(1)) != len(ev):
        sys.exit(f"parsed {len(ev)} events, the report counts {m.group(1) if m else 'none'}")
    return ev


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--binary", default="build-release/app/ayzek")
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--dataset", required=True, help="NAME:CACHE:FILE[,FILE...]")
    ap.add_argument("--detector", choices=DETECTORS, required=True)
    ap.add_argument("--shifts", default="127,263,391")
    ap.add_argument("--jobs", type=int, default=7)
    ap.add_argument("--label", help="name written in the detector column (default: --detector); "
                    "for another model's scores, give `--detector learned` with its cache and a label")
    ap.add_argument("--out", required=True)
    ap.add_argument("extra", nargs="*", help="after --: more ayzek options")
    a = ap.parse_args()

    name, cache, files = a.dataset.split(":", 2)
    paths = sorted(p for f in files.split(",") for p in glob.glob(f))
    if not paths:
        sys.exit(f"{name}: no files")
    span = data_span(paths)
    shifts = [0, *(s * sign for s in (int(x) for x in a.shifts.split(",") if x) for sign in (1, -1))]

    with tempfile.TemporaryDirectory(prefix="outcomes_") as tmp:
        def run(shift):
            cat = Path(tmp) / f"{name}_{shift}.csv"
            write_catalogue(a.catalog, cat, span, shift)
            cmd = [a.binary, "--speed", "0", "--no-color", "--no-magnitude", "--catalog", str(cat),
                   *DETECTORS[a.detector], *a.extra]
            if cache:
                cmd += ["--scores-in", cache]
            p = subprocess.run([*cmd, *paths], capture_output=True, text=True)
            if p.returncode or not CATALOGUE.search(p.stdout):
                sys.exit(f"{name} shift {shift}: ayzek exited {p.returncode} without a report"
                         f"\n{p.stderr[-800:]}")
            return shift, outcomes(p.stdout, shift)

        with concurrent.futures.ThreadPoolExecutor(a.jobs) as pool:
            results = dict(pool.map(run, shifts))

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["dataset", "detector", "shift", "origin", "type", "magnitude", "detected",
                    "detected_pre_s"])
        for shift in shifts:
            for (origin, kind, mag), (d, pre_s) in sorted(results[shift].items()):
                w.writerow([name, a.label or a.detector, shift, origin, kind, mag, d, pre_s])
    obs = results[0]
    n = {s: sum(d for d, _ in results[s].values()) for s in shifts}
    pre = sum(p for _, p in obs.values() if p != "")
    print(f"{name} {a.label or a.detector}: {n[0]} of {len(obs)} detected, {pre} before S; "
          + ", ".join(f"{s:+d}: {n[s]}/{len(results[s])}" for s in shifts[1:]))


if __name__ == "__main__":
    main()
