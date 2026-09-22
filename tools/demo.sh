#!/bin/sh
# Demo: replays the 2025-11-10 Sındırgı M4.9 and aftershocks from AFAD stations
# and compares the results with the AFAD catalogue (tests/catalogs/demo.csv, or
# $AYZEK_CATALOG).
#
#   tools/demo.sh                 3 stations, 10x, from 18:20
#   tools/demo.sh all             all 7 stations, full speed, whole 30 minutes
#   BUILD=build-pi tools/demo.sh  the aarch64 binary (under qemu on x86)
#
# Requires data/demo and models/; see docs/IMPLEMENTATION.md.
set -eu
cd "$(dirname "$0")/.."
BUILD=${BUILD:-build-release}
CATALOG=${AYZEK_CATALOG:-tests/catalogs/demo.csv}

if [ ! -f data/demo/DEMI.mseed ] || [ ! -f models/detector_s42.ayzw ]; then
    echo "missing data/demo or models: see docs/IMPLEMENTATION.md, 'First run'" >&2
    exit 1
fi

if [ "${1:-}" = all ]; then
    exec "$BUILD/app/ayzek" --speed 0 --catalog "$CATALOG" data/demo/*.mseed
fi
exec "$BUILD/app/ayzek" --speed 10 --from 2025-11-10T18:20:00 --catalog "$CATALOG" \
    data/demo/DEMI.mseed data/demo/MANT.mseed data/demo/BAND.mseed
