"""Rerun the learned detector on the held-out quiet windows, keeping its window
scores, which the test runs did not cache (they cached the sequence sets only).

    python3 -u tools/heldout_quiet_scores.py JOBS OUT [--jobs N]      (from ayzek_code/)

Exploratory, outside the held-out protocol: no setting changes, nothing is
retuned.

JOBS holds the 62 original commands of the test runs (the lines of each run's
commands.log whose --record lies under quiet/learned/), run from ayzek_code/.
Each is repeated unchanged except for two things (plan below): the recording
goes under OUT, and --scores writes the window probabilities to
OUT/<net>/scores/learned/<set>.tmp, renamed to <set> and marked .done when the
run ends well (the convention of tools/run_heldout.py; a set marked .done is
kept). The report goes to OUT/<net>/reports/<set>.txt (its stderr to
<set>.err), to be compared with the release-time report of the same window.
A command without --record data/runs/heldout/<net>/quiet/learned/<set>.rec, or
one that already has --scores, stops the script before anything runs.

One log line per finished window: count, set, its time, elapsed, ETA. Exit
status 1 when a window failed. Standard library only.
"""
import argparse
import re
import shlex
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def plan(line, out):
    """One job from a command line of commands.log: (set, argv to run, the
    scores' temporary and final directories, the report file). Makes the
    directories under `out`."""
    argv = shlex.split(line)
    i = argv.index("--record")
    m = re.fullmatch(r"data/runs/heldout/(\w+)/quiet/learned/(\w+)\.rec", argv[i + 1])
    if not m or "--scores" in argv:
        sys.exit(f"unexpected command: {line[:120]}")
    net, name = m.groups()
    base = out / net
    for sub in ("rec", "reports", "scores/learned"):
        (base / sub).mkdir(parents=True, exist_ok=True)
    final = base / "scores/learned" / name
    tmp = final.with_name(name + ".tmp")
    argv[i + 1] = str(base / "rec" / f"{name}.rec")
    argv[i + 2:i + 2] = ["--scores", str(tmp)]
    return name, argv, tmp, final, base / "reports" / f"{name}.txt"


def run(job):
    name, argv, tmp, final, report = job
    if (final / ".done").exists():
        return name, 0.0, "kept"
    t0 = time.time()
    with open(report, "w") as fh, open(report.with_suffix(".err"), "w") as eh:
        code = subprocess.run(argv, stdout=fh, stderr=eh).returncode
    if code != 0:
        return name, time.time() - t0, f"FAILED exit {code}"
    tmp.rename(final)
    (final / ".done").touch()
    return name, time.time() - t0, "ok"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("jobs")
    ap.add_argument("out")
    ap.add_argument("--jobs-parallel", "--jobs", dest="n", type=int, default=8)
    a = ap.parse_args()
    out = Path(a.out)
    jobs = [plan(line, out) for line in Path(a.jobs).read_text().splitlines() if line.strip()]
    print(f"{time.strftime('%FT%TZ', time.gmtime())} {len(jobs)} windows, {a.n} at a time -> {out}", flush=True)
    start, done, failed = time.time(), 0, 0
    with ThreadPoolExecutor(a.n) as pool:
        for fut in as_completed([pool.submit(run, j) for j in jobs]):
            name, secs, status = fut.result()
            done += 1
            failed += status.startswith("FAILED")
            elapsed = time.time() - start
            eta = elapsed / done * (len(jobs) - done)
            print(f"{time.strftime('%FT%TZ', time.gmtime())} {done}/{len(jobs)} {name} {status} "
                  f"{secs:.0f} s | elapsed {elapsed / 60:.0f} min, ETA {eta / 60:.0f} min", flush=True)
    print(f"{time.strftime('%FT%TZ', time.gmtime())} {'DONE' if not failed else f'FINISHED WITH {failed} FAILURES'}",
          flush=True)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
