#!/bin/sh
# Demo: replays the 2025-11-10 Sındırgı Mw 4.9 and the aftershocks around it,
# and compares the alarms with the AFAD catalogue.
#
#   tools/demo.sh                 the KO stations, 10x, triggering from 18:20
#   tools/demo.sh all             the same at full speed over the whole span
#   tools/demo.sh fetch           download the KO dataset (4 MB) and replay it
#   BUILD=build-pi tools/demo.sh  the aarch64 binary (under qemu on x86)
#   AYZEK_DATA=data/demo tools/demo.sh    a particular dataset directory
#
# It uses `data/demo_ko` (six KOERI stations, attached to the release) when it
# is there, and the AFAD `data/demo` otherwise: the AFAD recordings cannot be
# redistributed, so a fresh clone has only the KO one. The numbers differ
# between them, because the networks differ; both are in the README.
set -eu
cd "$(dirname "$0")/.."
BUILD=${BUILD:-build-release}
RELEASE=https://github.com/hogib/ayzek/releases/download/v0.1.0-beta.1/ayzek-demo_ko.tar.gz

if [ "${1:-}" = fetch ]; then
    shift
    mkdir -p data/demo_ko
    echo "fetching the KO demo dataset into data/demo_ko"
    curl -fsSL "$RELEASE" | tar xz -C data/demo_ko
fi

# Dataset: what the caller asked for, else KO, else AFAD.
if [ -n "${AYZEK_DATA:-}" ]; then
    DATA=$AYZEK_DATA
elif ls data/demo_ko/*.mseed >/dev/null 2>&1; then
    DATA=data/demo_ko
elif ls data/demo/*.mseed >/dev/null 2>&1; then
    DATA=data/demo
else
    echo "no demo data: run 'tools/demo.sh fetch' for the KO dataset (4 MB)," >&2
    echo "or build data/demo from an AFAD archive (README, 'Data and models')" >&2
    exit 1
fi
[ -f models/detector_s42.ayzw ] || { echo "missing models/" >&2; exit 1; }

# Catalogue: the excerpt that matches the dataset, unless one is given.
NAME=$(basename "$DATA")
CATALOG=${AYZEK_CATALOG:-tests/catalogs/$NAME.csv}
[ -f "$CATALOG" ] || CATALOG=$DATA/$NAME.csv          # as shipped in the tarball
[ -f "$CATALOG" ] || { echo "no catalogue for $DATA" >&2; exit 1; }

# The three nearest stations of the AFAD set are the published three-station
# replay; the KO set is small enough to use whole.
if [ "$NAME" = demo ] && [ "${1:-}" != all ]; then
    FILES="$DATA/DEMI.mseed $DATA/MANT.mseed $DATA/BAND.mseed"
else
    FILES=$(ls "$DATA"/*.mseed)
fi

echo "demo: $DATA against $CATALOG"
if [ "${1:-}" = all ]; then
    exec "$BUILD/app/ayzek" --speed 0 --catalog "$CATALOG" $FILES
fi
exec "$BUILD/app/ayzek" --speed 10 --from 2025-11-10T18:20:00 --catalog "$CATALOG" $FILES
