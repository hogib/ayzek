#!/bin/sh
# Reproduces every measured claim this project makes, and says which ones this
# machine cannot check. Claims and their expected values are bench/claims.json;
# tools/reproduce_report.py compares them with what this run measured.
#
#   tools/reproduce_all.sh                 everything available here
#   tools/reproduce_all.sh --quick         skip the STA/LTA comparison runs and
#                                          stages 6-8
#   SEISBENCH=1 tools/reproduce_all.sh     also recompute the published models'
#                                          scores (GPU, uv), not only replay them
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
#   6 quiet windows               102.8 h of AFAD data without catalogued
#                                 seismicity, both detectors: alarm rates,
#                                 coincident pairs, signal duration behind each
#                                 alarm (the full model: a few hours)
#   7 paired statistics           per-event outcomes on the three AFAD replays
#                                 and the shifted catalogues; McNemar, paired
#                                 bootstrap, the quiet-window rate test
#   8 published models            EQTransformer and PhaseNet at the settings
#                                 selected on the tuning replays (paper §4.3),
#                                 replayed from their window scores
#   9 second held-out sequence    Simav Mw 4.8 (paper §4.5): every detector at
#                                 its frozen setting, and both held-out
#                                 sequences pooled
#  10 FilterPicker                its two selected settings (paper §4.4) on the
#                                 four AFAD replays and the quiet windows
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

# Every quiet window replayed on its own from the end of its lead-in, as
# data/runs/quiet100/run.sh: replay_quiet DIR [ayzek options]. With
# QUIET_SCORES set, each window's probabilities come from QUIET_SCORES/NAME.
QUIET=data/quiet100
QUIET_SCORES=""
replay_quiet() {
    d=$1; shift
    mkdir -p "$d"
    tail -n +2 "$QUIET/quiet_index.csv" | {
        n=0
        while IFS=, read -r name start from end hours stations; do
            "$BUILD/app/ayzek" --speed 0 --no-magnitude --no-color \
                --catalog tests/catalogs/quiet100.csv --from "${from%%+*}" \
                ${QUIET_SCORES:+--scores-in $QUIET_SCORES/$name} \
                --record "$d/$name.rec" "$@" "$QUIET/$name"/*.mseed > "$d/$name.txt" 2>&1 &
            n=$((n + 1))
            [ $((n % 4)) -eq 0 ] && wait
        done
        wait
    }
}

say "6. quiet windows"
if [ "$QUICK" = 1 ]; then
    echo "(skipped: --quick)"
elif [ ! -f "$QUIET/quiet_index.csv" ]; then
    echo "no $QUIET (AFAD, not redistributable; tools/extract_quiet.py)"
else
    replay_quiet "$OUT/quiet100"
    replay_quiet "$OUT/quiet100_stalta" --detector stalta
    python3 tools/misfire_review.py --records "$OUT/quiet100/*.rec" \
        --csv "$OUT/quiet100/misfires.csv" > "$OUT/quiet100/misfires.txt" 2>&1 || true
    if command -v uv > /dev/null; then
        for r in quiet100 quiet100_stalta; do
            uv run --quiet --with obspy --with scipy python tools/transient_duration.py \
                --runs "$OUT/$r" --data "$QUIET" --out "$OUT/$r/durations.csv" \
                --json "$OUT/durations_$r.json" > "$OUT/$r/durations.txt" 2>&1 || true
        done
    else
        echo "no uv: signal durations not measured (needs obspy and scipy)"
    fi
    for r in quiet100 quiet100_stalta; do
        cat "$OUT/$r"/*.txt | awk -v r=$r '/Alarms raised/ {n += $3} END {print r ":", n, "alarms"}'
    done
fi

# The three AFAD replays: name and files, as in tools/scorecard.py.
REPLAYS="demo:data/demo/DEMI.mseed,data/demo/MANT.mseed,data/demo/BAND.mseed
marmara:data/marmara/*.mseed
sindirgi61:data/sindirgi61/*.mseed"
QUIET_ARGS=""
[ -d "$OUT/quiet100_stalta" ] && QUIET_ARGS="--quiet $QUIET/quiet_index.csv --quiet-runs learned=$OUT/quiet100"

say "7. paired statistics"
if [ "$QUICK" = 1 ]; then
    echo "(skipped: --quick)"
elif [ ! -d data/sindirgi61 ]; then
    echo "no data/sindirgi61 (AFAD, not redistributable)"
else
    mkdir -p "$OUT/outcomes"
    for spec in $REPLAYS; do
        n=${spec%%:*}; files=${spec#*:}
        # The detector's window probabilities, once; the outcome runs replay them.
        "$BUILD/app/ayzek" --speed 0 --no-magnitude --no-color --scores "$OUT/cache_$n" \
            $(echo "$files" | tr , ' ') > /dev/null 2>&1 || echo "scores failed: $n"
        for d in learned stalta_eqfa stalta_eqdet; do
            python3 tools/event_outcomes.py --binary "$BUILD/app/ayzek" --catalog "tests/catalogs/$n.csv" \
                --dataset "$n:$OUT/cache_$n:$files" --detector "$d" \
                --out "$OUT/outcomes/${n}_$d.csv" >> "$OUT/outcomes.log" 2>&1 || echo "outcomes failed: $n $d"
        done
    done
    python3 tools/significance.py --outcomes "$OUT/outcomes" \
        ${QUIET_ARGS:+$QUIET_ARGS stalta=$OUT/quiet100_stalta} \
        --json "$OUT/significance.json" > "$OUT/significance.txt"
    sed -n '/== held out/,/== tuning/p' "$OUT/significance.txt" | grep -E "vs|corrected" || true
fi

# The settings selected on the tuning replays before any evaluation run
# (ayzek_paper revisions2.md R9.1): label:score set:threshold:windows.
SB_SETTINGS="eqt_causal_fa:eqtransformer_instance_causal:0.8:2
eqt_causal_det:eqtransformer_instance_causal:0.3:8
eqt_look:eqtransformer_instance_lookahead:0.5:8
pn_causal:phasenet_instance_causal:0.7:2
pn_look_fa:phasenet_instance_lookahead:0.4:8
pn_look_det:phasenet_instance_lookahead:0.6:8"
SB=${SB:-data/runs/sb}

say "8. published models (EQTransformer, PhaseNet)"
if [ "$QUICK" = 1 ]; then
    echo "(skipped: --quick)"
elif [ ! -d data/sindirgi61 ]; then
    echo "no data/sindirgi61 (AFAD, not redistributable)"
elif [ "${SEISBENCH:-0}" != 1 ] && [ ! -d "$SB" ]; then
    echo "no window scores in $SB; SEISBENCH=1 computes them (GPU, uv)"
else
    if [ "${SEISBENCH:-0}" = 1 ]; then
        SB=$OUT/sb
        for set in eqtransformer:causal eqtransformer:lookahead phasenet:causal phasenet:lookahead; do
            m=${set%%:*}; mode=${set#*:}; tag=${m}_instance_$mode
            score() {   # NAME FILES...
                n=$1; shift
                uv run --quiet --with seisbench --with scipy python -W ignore tools/seisbench_scores.py \
                    --model "$m" --weights instance --mode "$mode" --out "$SB/$tag/$n" "$@" \
                    > "$SB/$tag.$(echo "$n" | tr / _).log" 2>&1 || echo "scores failed: $tag $n"
            }
            mkdir -p "$SB/$tag"
            for spec in $REPLAYS; do score "${spec%%:*}" $(echo "${spec#*:}" | tr , ' '); done
            for n in demo_ko marmara_ko quiet_ko; do score "$n" data/$n/*.mseed; done
            if [ -f "$QUIET/quiet_index.csv" ]; then
                for w in $(tail -n +2 "$QUIET/quiet_index.csv" | cut -d, -f1); do
                    score "quiet/$w" "$QUIET/$w"/*.mseed
                done
            fi
        done
    fi
    labels=""
    sb_quiet=""
    for s in $SB_SETTINGS; do
        label=${s%%:*}; r=${s#*:}; tag=${r%%:*}; r=${r#*:}
        args="--threshold ${r%%:*} --trigger-windows ${r#*:} --instant-threshold 2"
        labels="$labels${labels:+,}$label"
        for spec in $REPLAYS; do
            n=${spec%%:*}; files=${spec#*:}
            python3 tools/event_outcomes.py --binary "$BUILD/app/ayzek" --catalog "tests/catalogs/$n.csv" \
                --dataset "$n:$SB/$tag/$n:$files" --detector learned --label "$label" \
                --out "$OUT/outcomes/${n}_$label.csv" -- $args >> "$OUT/outcomes.log" 2>&1 \
                || echo "outcomes failed: $n $label"
        done
        python3 tools/scorecard.py --build "$BUILD" --only sindirgi61 --jobs 1 \
            --args="--scores-in $SB/$tag/sindirgi61 $args" \
            --json "$OUT/sb_${label}_sindirgi61.json" > "$OUT/sb_${label}_sindirgi61.txt" 2>&1 || true
        for n in demo_ko marmara_ko quiet_ko; do
            [ -d "$SB/$tag/$n" ] || continue
            python3 tools/scorecard.py --build "$BUILD" --only "$n" --jobs 1 \
                --args="--scores-in $SB/$tag/$n $args" \
                --json "$OUT/sb_${label}_$n.json" > "$OUT/sb_${label}_$n.txt" 2>&1 || true
        done
        if [ -n "$QUIET_ARGS" ] && [ -d "$SB/$tag/quiet" ]; then
            QUIET_SCORES=$SB/$tag/quiet
            replay_quiet "$OUT/sb_$label/quiet" $args
            QUIET_SCORES=""
            sb_quiet="$sb_quiet $label=$OUT/sb_$label/quiet"
        fi
    done
    python3 tools/significance.py --outcomes "$OUT/outcomes" --against "$labels" \
        ${QUIET_ARGS:+$QUIET_ARGS$sb_quiet} \
        --json "$OUT/significance_sb.json" > "$OUT/significance_sb.txt"
    sed -n '/== held out/,/== tuning/p' "$OUT/significance_sb.txt" | grep -E "vs|corrected" || true
fi

say "9. second held-out sequence (Simav)"
if [ "$QUICK" = 1 ]; then
    echo "(skipped: --quick)"
elif [ ! -d data/simav26 ] || [ ! -d data/sindirgi61 ]; then
    echo "no data/simav26 or data/sindirgi61 (AFAD, not redistributable)"
else
    files="data/simav26/*.mseed"
    "$BUILD/app/ayzek" --speed 0 --no-magnitude --no-color --scores "$OUT/cache_simav26" \
        $files > /dev/null 2>&1 || echo "scores failed: simav26"
    mkdir -p "$OUT/outcomes"
    for d in learned stalta_eqfa stalta_eqdet; do
        python3 tools/event_outcomes.py --binary "$BUILD/app/ayzek" --catalog tests/catalogs/simav26.csv \
            --dataset "simav26:$OUT/cache_simav26:$files" --detector "$d" \
            --out "$OUT/outcomes/simav26_$d.csv" >> "$OUT/outcomes.log" 2>&1 || echo "outcomes failed: simav26 $d"
    done
    against=stalta_eqfa,stalta_eqdet
    for s in $SB_SETTINGS; do
        label=${s%%:*}; r=${s#*:}; tag=${r%%:*}; r=${r#*:}
        [ -d "$SB/$tag/simav26" ] || continue
        python3 tools/event_outcomes.py --binary "$BUILD/app/ayzek" --catalog tests/catalogs/simav26.csv \
            --dataset "simav26:$SB/$tag/simav26:$files" --detector learned --label "$label" \
            --out "$OUT/outcomes/simav26_$label.csv" \
            -- --threshold "${r%%:*}" --trigger-windows "${r#*:}" --instant-threshold 2 \
            >> "$OUT/outcomes.log" 2>&1 || echo "outcomes failed: simav26 $label"
        against="$against,$label"
    done
    python3 tools/significance.py --outcomes "$OUT/outcomes" --against "$against" \
        --groups "simav26=simav26;both held out=sindirgi61,simav26" \
        --json "$OUT/significance_simav.json" > "$OUT/significance_simav.txt"
    python3 tools/significance.py --outcomes "$OUT/outcomes" --against "$against" \
        --groups "simav26=simav26;both held out=sindirgi61,simav26" --block-minutes 10 \
        --json "$OUT/significance_simav_block10.json" > "$OUT/significance_simav_block10.txt"
    grep -E "^==|vs|corrected" "$OUT/significance_simav.txt" | head -12 || true
fi

# FilterPicker's settings selected on the tuning replays (ayzek_paper
# revisions2.md R15.1): label:windows:thresholds.
FP_SETTINGS="fp_fa:6,10,0.4:15,10
fp_det:3,5,0.2:4,10"

say "10. FilterPicker"
if [ "$QUICK" = 1 ]; then
    echo "(skipped: --quick)"
elif [ ! -d data/sindirgi61 ]; then
    echo "no data/sindirgi61 (AFAD, not redistributable)"
else
    fp_quiet=""
    for s in $FP_SETTINGS; do
        label=${s%%:*}; r=${s#*:}
        args="--fp-windows ${r%%:*} --fp-thresholds ${r#*:}"
        for spec in $REPLAYS simav26:data/simav26/*.mseed; do
            n=${spec%%:*}; files=${spec#*:}
            [ -d "data/$n" ] || [ "$n" = demo ] || continue
            python3 tools/event_outcomes.py --binary "$BUILD/app/ayzek" --catalog "tests/catalogs/$n.csv" \
                --dataset "$n::$files" --detector filterpicker --label "$label" \
                --out "$OUT/outcomes/${n}_$label.csv" -- $args >> "$OUT/outcomes.log" 2>&1 \
                || echo "outcomes failed: $n $label"
        done
        if [ -n "$QUIET_ARGS" ]; then
            replay_quiet "$OUT/$label/quiet" --detector filterpicker $args
            fp_quiet="$fp_quiet $label=$OUT/$label/quiet"
        fi
    done
    python3 tools/significance.py --outcomes "$OUT/outcomes" --against fp_fa,fp_det \
        --groups "held out=sindirgi61;simav26=simav26;tuning=demo,marmara;both held out=sindirgi61,simav26" \
        ${QUIET_ARGS:+$QUIET_ARGS$fp_quiet} \
        --json "$OUT/significance_fp.json" > "$OUT/significance_fp.txt"
    grep -E "^==|vs|corrected|per day" "$OUT/significance_fp.txt" | head -16 || true
fi

say "report"
if [ "${VERIFY_FAILED:-0}" = 1 ]; then
    echo "stage 0 found a published number that no longer matches its artefact;"
    echo "see $OUT/verify_claims.txt"
fi
python3 tools/reproduce_report.py "$OUT" | tee "$OUT/report.txt"
echo
echo "reproduce: done, $OUT/report.txt"
