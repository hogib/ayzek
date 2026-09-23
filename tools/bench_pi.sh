#!/bin/sh
# Speed and accuracy on a Raspberry Pi, in one command.
#
#   tools/bench_pi.sh pi@raspberrypi.local     from here: copy, run there, bring
#                                              the results back and compare
#   tools/bench_pi.sh --local                  on the Pi itself, in this directory
#
# What it measures:
#   speed     ayzek_bench: per stage, per layer, and how many stations fit in
#             real time on that board (bench/*.json)
#   accuracy  the KO datasets through tools/scorecard.py, so the aarch64 build
#             is checked to give the same alarms as x86, not just to run
#
# It sends the aarch64 binaries, models/, the KO datasets and the catalogue
# excerpts: about 90 MB, or 25 MB with --no-quiet. Nothing is built on the Pi.
#
#   REMOTE=~/ayzek-bench   where to put it on the Pi
#   SETS="demo_ko"         datasets to score there (default: demo_ko quiet_ko)
set -eu
cd "$(dirname "$0")/.."
REMOTE=${REMOTE:-'~/ayzek-bench'}
SETS=${SETS:-"demo_ko quiet_ko"}
BENCH_ARGS=${BENCH_ARGS:-"--min-time 2"}

run_local() {
    # On the Pi: everything is beside this script's parent directory.
    mkdir -p bench
    echo "== machine"
    ./app/ayzek_bench --models models --group stages --min-time 0.1 2>&1 | head -1

    echo "\n== speed"
    for g in stages layers capacity; do
        ./app/ayzek_bench --models models --group "$g" $BENCH_ARGS \
            --json "bench/pi_$g.json" | tail -n +2
    done

    echo "\n== accuracy"
    only=$(echo "$SETS" | tr ' ' ',')
    python3 tools/scorecard.py --build . --only "$only" --jobs 2 \
        --json bench/pi_scorecard.json
    echo "\nresults: bench/pi_*.json"
}

if [ "${1:-}" = --local ]; then
    run_local
    exit 0
fi

HOST=${1:-}
[ -n "$HOST" ] || { echo "usage: tools/bench_pi.sh HOST | --local" >&2; exit 2; }
[ -x build-pi/app/ayzek ] || { echo "build build-pi first (README)" >&2; exit 1; }

echo "== sending to $HOST:$REMOTE"
ssh "$HOST" "mkdir -p $REMOTE/app $REMOTE/models $REMOTE/tools $REMOTE/tests/catalogs $REMOTE/data"
rsync -a build-pi/app/ayzek build-pi/bench/ayzek_bench "$HOST:$REMOTE/app/"
rsync -a models/ "$HOST:$REMOTE/models/"
rsync -a tools/scorecard.py tools/sweep_trigger.py tools/bench_pi.sh "$HOST:$REMOTE/tools/"
for s in $SETS; do
    [ -d "data/$s" ] || { echo "missing data/$s (tools/demo.sh fetch)" >&2; exit 1; }
    rsync -a "data/$s" "$HOST:$REMOTE/data/"
    rsync -a "tests/catalogs/$s.csv" "$HOST:$REMOTE/tests/catalogs/"
done

echo "== running on $HOST"
ssh "$HOST" "cd $REMOTE && SETS='$SETS' BENCH_ARGS='$BENCH_ARGS' sh tools/bench_pi.sh --local"

echo "== bringing the results back"
mkdir -p bench/results
rsync -a "$HOST:$REMOTE/bench/pi_*.json" bench/results/
echo "results in bench/results/pi_*.json"

if [ -f bench/results/x86.json ]; then
    echo "\n== x86 against the Pi"
    python3 tools/bench_compare.py bench/results/x86.json bench/results/pi_stages.json --threshold 0
fi
if [ -f bench/results/scorecard_baseline.json ]; then
    echo "\n== the same datasets, x86 against the Pi"
    python3 tools/scorecard.py --compare bench/results/scorecard_baseline.json \
        bench/results/pi_scorecard.json || true
fi
