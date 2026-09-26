#!/bin/sh
# Reproduces every measured claim this project makes, and says which ones this
# machine cannot check. Claims and their expected values are bench/claims.json;
# tools/reproduce_report.py compares them with what this run measured.
#
#   tools/reproduce_all.sh                 everything available here
#   tools/reproduce_all.sh --quick         skip the STA/LTA comparison runs
#   OUT=runs/xyz tools/reproduce_all.sh    results somewhere other than the default
#
# Stages, in order of what they depend on:
#
#   1 build and unit tests        agreement with PyTorch, scipy, torchaudio, ObsPy
#   2 scorecard, the detector     detection, unmatched alarms, false alarms,
#                                 magnitude and location on every dataset present
#   3 scorecard, STA/LTA          the same datasets with --detector stalta, at the
#                                 two tuned settings of docs/impl/09
#   4 benchmarks                  per stage, per layer, and real-time capacity
#   5 the aarch64 build           its numerical agreement under qemu, if built
#
# A dataset that is not on this machine is skipped and named in the report: the
# AFAD replays cannot be redistributed, the KO ones come with the release
# (tools/demo.sh fetch, or tools/fetch_ko_data.py).
set -eu
cd "$(dirname "$0")/.."
QUICK=${QUICK:-0}
[ "${1:-}" = --quick ] && QUICK=1
OUT=${OUT:-runs/reproduce-$(date -u +%Y%m%dT%H%M%SZ)}
BUILD=${BUILD:-build-release}
mkdir -p "$OUT"
echo "reproduce: results in $OUT"

ROOT=${ROOT:-$(cd "$(dirname "$0")/.." && pwd)}
say() { printf '\n=== %s\n' "$1"; }

say "0. published numbers against the files they came from"
# Cheap: reads artefacts earlier runs wrote and checks the published figure is
# still the one in the file. Catches drift in seconds; needs no replay.
python3 "$ROOT/tools/verify_claims.py" | tee "$OUT/verify_claims.txt" || VERIFY_FAILED=1

say "1. build and unit tests"
meson setup "$BUILD" --buildtype=release >/dev/null 2>&1 || true
ninja -C "$BUILD" > "$OUT/build.log" 2>&1
meson test -C "$BUILD" > "$OUT/tests.log" 2>&1 && echo "tests: all passed" || {
    echo "tests: FAILURES, see $OUT/tests.log"; tail -20 "$OUT/tests.log"; }

say "2. scorecard, the learned detector"
python3 tools/scorecard.py --build "$BUILD" --bench --jobs 3 \
    --json "$OUT/scorecard.json" > "$OUT/scorecard.txt" 2>&1 || true
tail -3 "$OUT/scorecard.txt"

if [ "$QUICK" = 0 ]; then
    say "3. scorecard, STA/LTA at the two tuned settings (docs/impl/09)"
    # Equal unmatched alarms: the defaults of --detector stalta.
    python3 tools/scorecard.py --build "$BUILD" --jobs 3 --args '--detector stalta' \
        --json "$OUT/scorecard_stalta_eqfa.json" > "$OUT/scorecard_stalta_eqfa.txt" 2>&1 || true
    # Equal detections: the setting that matches the model's detection rate.
    python3 tools/scorecard.py --build "$BUILD" --jobs 3 \
        --args '--detector stalta --stalta-band 2,10 --sta 0.5 --lta 10 --stalta-on 6' \
        --json "$OUT/scorecard_stalta_eqdet.json" > "$OUT/scorecard_stalta_eqdet.txt" 2>&1 || true
    tail -2 "$OUT/scorecard_stalta_eqfa.txt"
else
    echo "(3. STA/LTA runs skipped: --quick)"
fi

say "4. benchmarks"
for g in stages layers capacity; do
    "$BUILD/bench/ayzek_bench" --models models --group "$g" \
        --json "$OUT/bench_$g.json" > "$OUT/bench_$g.txt" 2>&1 || true
done
grep -E "per-window total|magnitude estimate" "$OUT/bench_stages.txt" || true

say "5. the aarch64 build"
if [ -d build-pi ]; then
    meson test -C build-pi --timeout-multiplier 10 > "$OUT/tests_pi.log" 2>&1 \
        && echo "aarch64 under qemu: all passed" \
        || echo "aarch64 under qemu: FAILURES, see $OUT/tests_pi.log"
else
    echo "no build-pi; cross-compile it to check the aarch64 agreement (README)"
fi

say "report"
if [ "${VERIFY_FAILED:-0}" = 1 ]; then
    echo "stage 0 found a published number that no longer matches its artefact;"
    echo "see $OUT/verify_claims.txt"
fi
python3 tools/reproduce_report.py "$OUT" | tee "$OUT/report.txt"
echo
echo "reproduce: done, $OUT/report.txt"
