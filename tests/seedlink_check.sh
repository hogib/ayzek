#!/bin/sh
# End-to-end tests of the live SeedLink input: ayzek against
# tools/seedlink_mock.py serving the public KO demo set (data/demo_ko,
# catalogue tests/catalogs/demo_ko.csv) at AYZEK_LIVE_SPEED times real time
# (default 20).
#   seedlink_check.sh AYZEK ROOT MODE
# MODE live: the mock serves the six stations without faults. The live run's
#   report (alarms, magnitudes, locations, catalogue comparison) and its
#   --record file must equal those of the file replay of the same files, and so
#   must those of a replay of its --archive.
# MODE reconnect: the mock drops the connection three times halfway through a
#   packet, stalls twice (answering ayzek's keepalive with a two-packet INFO
#   reply the first time, ignoring it the second, past ayzek's read timeout),
#   sends packets again after every resume, and numbers packets from FFFFF0 so
#   that they wrap. The archive must hold every record of the files exactly
#   once, and the report and record file must equal the file replay's.
# Exits 77 (skipped) if the data, python3 or a binary that runs here is
# missing. Set AYZEK_KEEP=1 to keep the run directory.
set -eu
AYZEK=$1
ROOT=$2
MODE=$3
SPEED=${AYZEK_LIVE_SPEED:-20}
DATA=$ROOT/data/demo_ko
CAT=$ROOT/tests/catalogs/demo_ko.csv
MOCK=$ROOT/tools/seedlink_mock.py

FILES=$(ls "$DATA"/*.mseed 2>/dev/null || true)
[ -n "$FILES" ] || { echo "skipped: no $DATA/*.mseed"; exit 77; }
for f in "$CAT" "$MOCK" "$ROOT/models/detector_s42.ayzw"; do
    [ -e "$f" ] || { echo "skipped: $f not found"; exit 77; }
done
command -v python3 >/dev/null || { echo "skipped: no python3 for the mock server"; exit 77; }
# Live input runs natively only: a cross-compiled binary under qemu would not
# keep up with the mock.
"$AYZEK" --help >/dev/null 2>&1 || { echo "skipped: cannot run $AYZEK here"; exit 77; }

T=$(mktemp -d "${TMPDIR:-/tmp}/ayzek_seedlink.XXXXXX")
MOCK_PID=
LIVE_PID=
cleanup() {
    [ -n "$LIVE_PID" ] && kill "$LIVE_PID" 2>/dev/null || true
    [ -n "$MOCK_PID" ] && kill "$MOCK_PID" 2>/dev/null || true
    if [ "${AYZEK_KEEP:-0}" = 1 ]; then echo "run directory kept: $T"; else rm -rf "$T"; fi
}
trap cleanup EXIT

STREAMS=
for f in $FILES; do
    s=$(basename "$f" .mseed)
    STREAMS=${STREAMS:+$STREAMS,}KO.$s
done

# The report section (alarms through the summary) and the recorded station
# outputs, which must be identical between runs apart from compute times (the
# compute_ms of D and P lines). The span line of a recording is written last
# by a live run and first by a replay, so it is compared on its own.
report() { sed -n '/^REPORT/,/^Station processing/p' "$1" | grep -v -e '^Station processing' -e '^Recording written'; }
recorded() {
    awk -F '\t' 'BEGIN { OFS = "\t" } $1 == "D" { $6 = "-" } $1 == "P" { $9 = "-" } $1 != "T"' "$1"
}
span() { grep '^T' "$1"; }

echo "file replay of $(echo "$FILES" | wc -w) stations"
# shellcheck disable=SC2086
"$AYZEK" --speed 0 --no-color --models "$ROOT/models" --catalog "$CAT" \
    --record "$T/replay.rec" $FILES > "$T/replay.txt"

case "$MODE" in
live)
    FAULTS=""
    LIVE_OPTS=""
    ;;
reconnect)
    FAULTS="--drop-after 2000 --drop-after 1500 --drop-after 1000 --seq-start FFFFF0 \
            --resume-overlap 3 --info-packets 2 --stall-after 6000:6 --stall-after 8000:9:mute"
    LIVE_OPTS="--seedlink-keepalive 2 --seedlink-timeout 5"
    ;;
*)
    echo "unknown mode $MODE"; exit 2 ;;
esac

# shellcheck disable=SC2086
python3 "$MOCK" --port 0 --port-file "$T/port" --speed "$SPEED" --exit-when-done \
    $FAULTS $FILES 2> "$T/mock.log" &
MOCK_PID=$!
for _ in $(seq 150); do [ -s "$T/port" ] && break; sleep 0.2; done
[ -s "$T/port" ] || { cat "$T/mock.log"; echo "FAIL: the mock server did not start"; exit 1; }
PORT=$(cat "$T/port")

echo "live run against the mock at ${SPEED}x on port $PORT ($MODE)"
# --live-wait well above any wait here, so that the network stage keeps the
# replay's order even on a loaded machine; the report counts what it did not.
# shellcheck disable=SC2086
"$AYZEK" --no-color --models "$ROOT/models" --catalog "$CAT" --record "$T/live.rec" \
    --archive "$T/archive" --seedlink "127.0.0.1:$PORT" --streams "$STREAMS" \
    --live-wait 300 --duration 3600 $LIVE_OPTS > "$T/live.txt" 2>&1 &
LIVE_PID=$!

# The mock exits once it has sent every record and closed the connection
# behind them; ayzek has then received them all, and is stopped.
wait "$MOCK_PID" || true
MOCK_PID=
for _ in $(seq 300); do
    grep -q 'connection closed by the server' "$T/live.txt" && break
    sleep 0.2
done
kill -TERM "$LIVE_PID"
status=0
wait "$LIVE_PID" || status=$?
LIVE_PID=
if [ "$status" != 0 ]; then
    tail -30 "$T/live.txt"; echo "FAIL: live run exited with status $status"; exit 1
fi
grep -E '^(seedlink|warn|ALARM|ISSUED)' "$T/live.txt" | grep -v 'crosses a gap' | head -60
sed -n '/^SeedLink input/,$p' "$T/live.txt"

fail=0
check() { if eval "$1"; then echo "ok:   $2"; else echo "FAIL: $2"; fail=1; fi; }

report "$T/replay.txt" > "$T/replay.report"
report "$T/live.txt" > "$T/live.report"
recorded "$T/replay.rec" > "$T/replay.d"
recorded "$T/live.rec" > "$T/live.d"
check '[ -s "$T/replay.report" ] && cmp -s "$T/replay.report" "$T/live.report"' \
    "the live report equals the file replay's ($(grep -c '^Event #' "$T/replay.report" || true) alarms)"
check 'cmp -s "$T/replay.d" "$T/live.d" && [ "$(span "$T/replay.rec")" = "$(span "$T/live.rec")" ]' \
    "the live --record file equals the file replay's ($(wc -l < "$T/replay.d") lines, compute times aside)"
check 'grep -q "network stage: 0 messages processed out of order" "$T/live.txt"' \
    "no message processed out of order"
# The input's station rows: code, received, accepted, rejected, dropped, ...
check 'awk "/^SeedLink input/ {on=1}
            on && \$1 ~ /^[A-Z][A-Z0-9]*\$/ && \$2 ~ /^[0-9]+\$/ && \$5 ~ /^[0-9]+\$/ {rows++; if (\$5 != 0) bad=1}
            END {exit bad || rows == 0}" "$T/live.txt"' \
    "no record dropped"
check 'grep -q "^  AFAD .* at 18:20:51" "$T/live.txt"' "the 2025-11-10 M4.9 alarmed live"

# Every record of the files received exactly once: the archive holds the
# same records as the files, none twice and none missing.
check 'python3 - "$T/archive" $FILES <<EOF
import sys, pathlib, collections
def records(paths):
    out = collections.Counter()
    for p in paths:
        b = pathlib.Path(p).read_bytes()
        out.update(b[i:i + 512] for i in range(0, len(b), 512))
    return out
archive = records(sorted(pathlib.Path(sys.argv[1]).glob("*.mseed")))
files = records(sys.argv[2:])
print(f"      archive {sum(archive.values())} records, files {sum(files.values())}, "
      f"{sum((archive - files).values())} extra, {sum((files - archive).values())} missing, "
      f"{sum(1 for n in archive.values() if n > 1)} received twice")
sys.exit(archive != files)
EOF' "the archive holds every record of the files exactly once"

if [ "$MODE" = live ]; then
    echo "replay of the live run's archive"
    "$AYZEK" --speed 0 --no-color --models "$ROOT/models" --catalog "$CAT" \
        --record "$T/archive.rec" "$T"/archive/*.mseed > "$T/archive.txt"
    report "$T/archive.txt" > "$T/archive.report"
    recorded "$T/archive.rec" > "$T/archive.d"
    check 'cmp -s "$T/live.report" "$T/archive.report"' \
        "the replay of the archive gives the live run's report"
    check 'cmp -s "$T/live.d" "$T/archive.d" && [ "$(span "$T/archive.rec")" = "$(span "$T/live.rec")" ]' \
        "the replay of the archive gives the live run's --record file"
else
    check 'grep -Eq "^  [0-9]+ connections" "$T/live.txt" &&
           [ "$(sed -n "s/^  \([0-9]*\) connections.*/\1/p" "$T/live.txt")" -ge 5 ]' \
        "reconnected after every drop and after the read timeout"
    check 'grep -q "nothing received for 5 s" "$T/live.txt"' \
        "the read timeout ended the connection whose keepalive went unanswered"
    check 'sed -n "s/.* \([0-9]*\) duplicates discarded.*/\1/p" "$T/live.txt" | grep -qv "^0$"' \
        "packets sent again after a resume were discarded"
    check 'sed -n "s/^  \([0-9]*\) keepalives sent, \([0-9]*\) INFO packets.*/\1 \2/p" "$T/live.txt" |
           awk "{exit !(\$1 >= 1 && \$2 >= 2)}"' \
        "keepalives were sent and the two-packet INFO replies read"
fi

if [ "$fail" != 0 ]; then
    echo "--- diff of the reports (file replay, live)"
    diff "$T/replay.report" "$T/live.report" | head -40 || true
    AYZEK_KEEP=1
    exit 1
fi
echo "PASS: live input ($MODE) gives the file replay's alarms"
