#!/bin/sh
# End-to-end test: replays DEMI, MANT and BAND through ayzek and requires the
# 2025-11-10 M4.9 to be declared, located within 10 km of the AFAD epicentre,
# and assigned a magnitude between 3.9 and 5.9.
# Exits 77 (skipped) if the demo data, models or catalogue are missing.
set -eu
AYZEK=$1
ROOT=$2
CATALOG=${AYZEK_CATALOG:-$HOME/Projects/sismokaos/data_downloader/catalogs/catalog_afad_full_2026-08-30.csv}
for f in "$ROOT/data/demo/DEMI.mseed" "$ROOT/data/demo/MANT.mseed" "$ROOT/data/demo/BAND.mseed" \
         "$ROOT/models/detector_s42.ayzw" "$CATALOG"; do
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
      "$ROOT/data/demo/DEMI.mseed" "$ROOT/data/demo/MANT.mseed" "$ROOT/data/demo/BAND.mseed")
echo "$out" | grep -E 'EVENT|LOCATE|MAG|^catalog'
line=$(echo "$out" | grep -E '^catalog +18:20:51' || true)
echo "$line" | grep -q 'alert' || { echo "FAIL: M4.9 not declared"; exit 1; }
km=$(echo "$line" | sed -n 's/.*located \([0-9.]*\) km off.*/\1/p')
[ -n "$km" ] || { echo "FAIL: M4.9 not located"; exit 1; }
awk -v k="$km" 'BEGIN { exit !(k < 10) }' || { echo "FAIL: M4.9 located $km km off"; exit 1; }
mag=$(echo "$line" | sed -n 's/.*, M\([0-9.]*\) (.*/\1/p')
[ -n "$mag" ] || { echo "FAIL: M4.9 has no magnitude"; exit 1; }
awk -v m="$mag" 'BEGIN { exit !(m > 3.9 && m < 5.9) }' || { echo "FAIL: M4.9 estimated M$mag"; exit 1; }
echo "PASS: M4.9 declared, located $km km from the catalogue epicentre, magnitude M$mag"
