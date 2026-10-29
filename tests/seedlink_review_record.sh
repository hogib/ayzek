#!/bin/sh
# Review test (branch seedlink-review): an error raised in main() after the
# threads have started must end a live run. Here --record names a file that
# cannot be written; the Recorder is constructed after the live input and
# processor threads have started, its exception unwinds through the
# std::jthread vector, whose destructors join threads that only stop on
# SIGINT/SIGTERM or --duration: the run hangs, retrying the server, and never
# reports the error. The server is a closed local port; nothing is contacted.
#   seedlink_review_record.sh AYZEK ROOT
set -u
AYZEK=$1
ROOT=$2
[ -e "$ROOT/models/detector_s42.ayzw" ] || { echo "skipped: no models"; exit 77; }
"$AYZEK" --help >/dev/null 2>&1 || { echo "skipped: cannot run $AYZEK here"; exit 77; }
T=$(mktemp -d "${TMPDIR:-/tmp}/ayzek_review.XXXXXX")
trap 'rm -rf "$T"' EXIT

"$AYZEK" --no-color --models "$ROOT/models" --seedlink 127.0.0.1:1 \
    --streams KO.DKL --record "$T/missing/dir/x.rec" > "$T/out.txt" 2>&1 &
PID=$!
for _ in $(seq 120); do
    kill -0 "$PID" 2>/dev/null || break
    sleep 0.5
done
if kill -0 "$PID" 2>/dev/null; then
    kill -INT "$PID"; sleep 5; kill -KILL "$PID" 2>/dev/null
    tail -3 "$T/out.txt"
    echo "FAIL: ayzek still running 60 s after an unwritable --record"
    exit 1
fi
wait "$PID"
status=$?
if [ "$status" = 0 ] || ! grep -q "cannot write" "$T/out.txt"; then
    tail -3 "$T/out.txt"; echo "FAIL: exit $status without the error"; exit 1
fi
echo "PASS: the run ended with the error"
