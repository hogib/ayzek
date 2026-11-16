"""Where the learned detector's outputs level off, for the ensemble and for each of its
three separately trained members (seeds 42, 43, 44), on development replays.
EXPLORATORY: no test data, no setting changed, nothing retuned.

    python3 -u tools/plateau_seeds.py OUT [--sets sindirgi61 simav26]      (from ayzek_code/)

Inputs: build-linux/app/ayzek, models/, data/<set>/*.mseed.

ayzek loads detector_s42/43/44.ayzw from --models and averages them. A member is run
alone by giving it a models directory in which all three names link to that member's
file (the average of three equal outputs is that output); every other file links to
models/ unchanged. Each run writes its window probabilities with --scores (to
OUT/<run>/<set>.tmp, renamed and marked .done when the run ends well; a set marked
.done is not run again).

The rule, per run and set, over all windows of all stations (plateau_row below): the
count, the maximum, the 99th, 99.9th and 99.99th percentiles (the sorted value at
position round(q * n), counted from 1), and the windows at or above 0.895, 0.9 and
0.905. -> OUT/plateau.csv (through plateau.csv.tmp), OUT/provenance.json; scores and
reports under OUT/<run>/. One log line per finished run. Exits nonzero if nothing is
produced.
"""
import argparse
import csv
import glob
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

RUNS = {"ensemble": None, "s42": 42, "s43": 43, "s44": 44}


def log(m):
    print(f"{time.strftime('%FT%TZ', time.gmtime())} {m}", flush=True)


def plateau_row(run, st, probabilities):
    """The row of plateau.csv for one run and set from its window probabilities."""
    p = sorted(probabilities)
    n = len(p)
    q = lambda x: p[min(n - 1, max(0, int(round(x * n)) - 1))]  # noqa: E731
    return {"run": run, "set": st, "windows": n, "max": f"{p[-1]:.6f}", "p99": f"{q(0.99):.6f}",
            "p99.9": f"{q(0.999):.6f}", "p99.99": f"{q(0.9999):.6f}",
            "ge_0.895": sum(x >= 0.895 for x in p), "ge_0.9": sum(x >= 0.9 for x in p),
            "ge_0.905": sum(x >= 0.905 for x in p)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("out")
    ap.add_argument("--sets", nargs="+", default=["sindirgi61", "simav26"])
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    todo = [(r, s) for r in RUNS for s in a.sets]
    t0 = time.time()
    for i, (run, st) in enumerate(todo, 1):
        models = Path("models").resolve()
        if RUNS[run] is not None:
            md = out / f"models_{run}"
            md.mkdir(exist_ok=True)
            for f in models.iterdir():
                dst = md / f.name
                if dst.is_symlink() or dst.exists():
                    dst.unlink()
                src = models / f"detector_s{RUNS[run]}.ayzw" if f.name.startswith("detector_s") else f
                dst.symlink_to(src)
            models = md
        sc = out / run / st
        if not (sc / ".done").exists():
            sc.parent.mkdir(parents=True, exist_ok=True)
            tmp = sc.with_name(st + ".tmp")
            files = sorted(glob.glob(f"data/{st}/*.mseed"))
            if not files:
                sys.exit(f"data/{st}: no mseed files")
            with open(out / run / f"{st}.txt", "w") as fh:
                code = subprocess.run(["build-linux/app/ayzek", "--speed", "0", "--no-color", "--models", str(models),
                                       "--scores", str(tmp)] + files, stdout=fh, stderr=subprocess.STDOUT).returncode
            if code != 0:
                sys.exit(f"{run} {st}: ayzek exited {code}")
            if sc.exists():
                sys.exit(f"{sc} exists without .done")
            tmp.rename(sc)
            (sc / ".done").touch()
        p = sorted(float(r["probability"]) for f in sc.glob("*.csv") for r in csv.DictReader(open(f)))
        rows.append(plateau_row(run, st, p))
        el = time.time() - t0
        log(f"{i}/{len(todo)} {run} {st}: {len(p)} windows, max {p[-1]:.4f}, >=0.9: {rows[-1]['ge_0.9']} | "
            f"elapsed {el / 60:.1f} min, ETA {el / i * (len(todo) - i) / 60:.1f} min")
    if not rows:
        sys.exit("nothing produced")
    tmp = out / "plateau.csv.tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    tmp.rename(out / "plateau.csv")
    sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()  # noqa: E731
    (out / "provenance.json").write_text(json.dumps({
        "argv": sys.argv, "note": "exploratory; development replays only",
        "commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
        "binary_sha256": sha("build-linux/app/ayzek"),
        "models_sha256": {f"detector_s{s}.ayzw": sha(f"models/detector_s{s}.ayzw") for s in (42, 43, 44)},
        "sets": a.sets, "written_utc": time.strftime('%FT%TZ', time.gmtime())}, indent=1))
    log("DONE")


if __name__ == "__main__":
    main()
