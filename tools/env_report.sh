#!/bin/sh
# Regenerates the table in docs/impl/12-environment.md.
#
#   tools/env_report.sh            print it
#   tools/env_report.sh > /tmp/env.md
#
# Reports what a number here was produced with: the build toolchain, the two
# Python environments (they differ, deliberately reported apart), the machine,
# and the repository commits. Anything missing prints "-" rather than failing,
# so this runs on a fresh checkout too.
set -u
cd "$(dirname "$0")/.."

v() { "$@" 2>/dev/null | head -1 || echo "-"; }

echo "## Build toolchain"
echo
printf '| component | version |\n|---|---|\n'
printf '| zig | %s |\n'        "$(v zig version)"
printf '| clang inside zig | %s |\n' "$(zig cc --version 2>/dev/null | head -1 | sed 's/.*version //')"
printf '| musl (bundled) | %s |\n'  "$(sed -n 's/.*VERSION "\(.*\)".*/\1/p' /usr/lib/zig/libc/musl/src/internal/version.h 2>/dev/null || echo -)"
printf '| mingw-w64 (bundled) | %s |\n' "$(grep -hoE '__MINGW64_VERSION_MAJOR[[:space:]]+[0-9]+' /usr/lib/zig/libc/include/any-windows-any/_mingw_mac.h 2>/dev/null | grep -oE '[0-9]+$' || echo -)"
printf '| gcc | %s |\n'        "$(v gcc -dumpfullversion)"
printf '| clang (system) | %s |\n' "$(clang --version 2>/dev/null | head -1 | sed 's/.*version //')"
printf '| glibc | %s |\n'      "$(ldd --version 2>/dev/null | head -1 | sed 's/.*) //')"
printf '| meson | %s |\n'      "$(v meson --version)"
printf '| ninja | %s |\n'      "$(v ninja --version)"
printf '| qemu-aarch64 | %s |\n' "$(qemu-aarch64 --version 2>/dev/null | head -1 | sed 's/.*version //')"

echo
echo "## Configured build directories"
echo
printf '| dir | host compiler |\n|---|---|\n'
for b in build-*; do
    [ -f "$b/meson-info/intro-compilers.json" ] || continue
    printf '| %s | %s |\n' "$b" "$(python3 - "$b" <<'PY'
import json, sys
d = json.load(open(f"{sys.argv[1]}/meson-info/intro-compilers.json"))
c = d.get("host", {}).get("cpp") or d.get("host", {}).get("c") or {}
print(f'{c.get("id","?")} {c.get("version","?")}')
PY
)"
done

echo
echo "## Python environments"
echo
for venv in "$@"; do
    [ -x "$venv/.venv/bin/python" ] || continue
    echo
    echo "### $(basename "$venv")"
    "$venv/.venv/bin/python" - <<'PY'
import importlib.metadata as m, sys
print(f"python {sys.version.split()[0]}")
for p in ("torch","numpy","scipy","obspy","matplotlib","pandas","h5py","scikit-learn"):
    try:
        print(f"{p} {m.version(p)}")
    except Exception:
        pass
PY
done

echo
echo "## Machine"
echo
printf '| CPU | %s |\n' "$(lscpu 2>/dev/null | sed -n 's/^Model name: *//p' | head -1)"
printf '| RAM | %s GB |\n' "$(free -g 2>/dev/null | awk 'NR==2{print $2}')"
printf '| GPU | %s |\n' "$(nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader 2>/dev/null | head -1)"
printf '| kernel | %s |\n' "$(uname -sr)"

echo
echo "## Repositories"
echo
for r in . ../ayzek_paper ../../../sismokaos/cascade_impl ../../../sismokaos/cnn_earthquake; do
    [ -d "$r/.git" ] || continue
    printf '| %-14s | %s | %s uncommitted |\n' "$(basename "$(cd "$r" && pwd)")" \
        "$(git -C "$r" rev-parse --short HEAD 2>/dev/null)" \
        "$(git -C "$r" status --porcelain 2>/dev/null | wc -l)"
done
