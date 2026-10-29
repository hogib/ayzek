#!/bin/sh
# End-to-end test: replays a demo set through ayzek and requires the
# 2025-11-10 M4.9 to be declared, located within 10 km of the catalogue
# epicentre, and, with MAG=1, assigned a magnitude between 3.9 and 5.9.
#   demo_check.sh AYZEK ROOT SET MAG
# SET is demo (AFAD DEMI, MANT and BAND; not redistributable) or demo_ko (six
# KO stations, attached to the release). MAG=0 for demo until the TU stations
# have gains for the Pd magnitude (models/pd_gain.csv), which needs them.
# Exits 77 (skipped) if the data or models are missing. The catalogue is
# tests/catalogs/SET.csv; AYZEK_CATALOG overrides it.
set -eu
AYZEK=$1
ROOT=$2
SET=$3
MAG=$4
CATALOG=${AYZEK_CATALOG:-$ROOT/tests/catalogs/$SET.csv}
if [ "$SET" = demo ]; then
    FILES="$ROOT/data/demo/DEMI.mseed $ROOT/data/demo/MANT.mseed $ROOT/data/demo/BAND.mseed"
else
    FILES=$(ls "$ROOT/data/$SET"/*.mseed 2>/dev/null || true)
    [ -n "$FILES" ] || { echo "skipped: no $ROOT/data/$SET/*.mseed"; exit 77; }
fi
for f in $FILES "$ROOT/models/detector_s42.ayzw" "$CATALOG"; do
    [ -e "$f" ] || { echo "skipped: $f not found"; exit 77; }
done
# A cross-compiled binary cannot run natively. Meson's exe_wrapper applies only to
# the test executable itself, so this script uses qemu-aarch64 directly.
RUN="$AYZEK"
if ! "$AYZEK" --help >/dev/null 2>&1; then
    command -v qemu-aarch64 >/dev/null || { echo "skipped: cannot run $AYZEK"; exit 77; }
    RUN="qemu-aarch64 $AYZEK"
fi
out=$($RUN --speed 0 --no-color --models "$ROOT/models" --from 2025-11-10T18:19:00 --catalog "$CATALOG" \
      $FILES)
echo "$out" | grep -E '^(ALARM|PDMAG|LOCATE|PLOC)|detected|AFAD|SUMMARY|correct|false|missed'
# The report block for the M4.9: the headline gives ayzek's magnitude, the AFAD
# line gives the epicentre error.
block=$(echo "$out" | awk '/^Event #/ {head=$0} /AFAD .* at 18:20:51/ {print head; print $0}')
echo "$block" | grep -q 'AFAD' || { echo "FAIL: M4.9 not detected"; exit 1; }
km=$(echo "$block" | sed -n 's/.*epicentre \([0-9.]*\) km off.*/\1/p')
[ -n "$km" ] || { echo "FAIL: M4.9 not located"; exit 1; }
awk -v k="$km" 'BEGIN { exit !(k < 10) }' || { echo "FAIL: M4.9 located $km km off"; exit 1; }
if [ "$MAG" = 0 ]; then
    echo "PASS: M4.9 detected, located $km km from the catalogue epicentre (magnitude not checked)"
    exit 0
fi
mag=$(echo "$block" | sed -n 's/.*detected M\([0-9.]*\) earthquake.*/\1/p')
[ -n "$mag" ] || { echo "FAIL: M4.9 has no magnitude"; exit 1; }
awk -v m="$mag" 'BEGIN { exit !(m > 3.9 && m < 5.9) }' || { echo "FAIL: M4.9 estimated M$mag"; exit 1; }
echo "PASS: M4.9 detected, located $km km from the catalogue epicentre, magnitude M$mag"
