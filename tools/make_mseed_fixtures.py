"""Writes the small miniSEED files that tests/test_mseed.cpp decodes, one per
encoding the decoder supports, with ObsPy's reading of each as the reference.

Each file is one 100 Hz HHZ channel of synthetic counts in 512-byte records:
a random walk in stretches whose differences need 4 to 30 bits (Steim2), 8,
16 or 32 bits (Steim1), or stay within 16-bit range (INT16), so that every
packing the encodings have occurs. ObsPy (libmseed) writes the file and reads
it back; the samples it reads are stored as NAME.i32 (little-endian int32),
the reference the C++ decoder must equal sample for sample.

The little-endian integer files have a big-endian header and little-endian
data (blockette 1000 word order 0), the layout of CATL's uncompressed records
in the TDVMS archive, which ObsPy does not write: they are ObsPy's big-endian
records with the data swapped and the word order changed here, and ObsPy's
reading of them is the reference all the same (it warns "Inconsistent word
order" and reads them as written). A record whose header is little-endian too
is not supported by the decoder.

    uv run --with obspy python tools/make_mseed_fixtures.py --out tests/reference/mseed

Writes fixtures.json beside the files: the command, every parameter, the git
commit, the library versions and, per file, its records, samples and SHA-256.
Exits nonzero if a file comes out empty or does not read back as written.
"""
import argparse
import hashlib
import json
import platform
import struct
import subprocess
import sys
from pathlib import Path

import numpy as np
import obspy
from obspy import Trace, UTCDateTime

START = UTCDateTime(2025, 11, 10, 18, 0, 0)

# name: ObsPy encoding, data byte order, difference bounds of the stretches
FIXTURES = {
    "steim1": ("STEIM1", ">", [100, 30000, 1_000_000_000, 7, 32767, 2**31 - 1]),
    "steim2": ("STEIM2", ">", [7, 15, 31, 127, 511, 16383, 2**29 - 1]),
    "int16_be": ("INT16", ">", [100, 3000, 30000]),
    "int16_le": ("INT16", "<", [100, 3000, 30000]),
    "int32_be": ("INT32", ">", [100, 100_000, 1_000_000_000]),
    "int32_le": ("INT32", "<", [100, 100_000, 1_000_000_000]),
}


def walk(rng, bounds, stretch, limit):
    """Random walk: `stretch` samples per bound, differences within +-bound,
    values kept within +-limit by reflecting the step."""
    x, out = 0, []
    for b in bounds:
        for _ in range(stretch):
            d = int(rng.integers(-b, b, endpoint=True))
            if abs(x + d) > limit:
                d = -d
            if abs(x + d) > limit:
                d = 0
            x += d
            out.append(x)
    return np.array(out, dtype=np.int64)


def swap_data(path, width):
    """Big-endian integer records -> the same with little-endian data and
    blockette 1000 word order 0; the header stays big-endian."""
    raw = bytearray(path.read_bytes())
    fmt = "h" if width == 2 else "i"
    for off in range(0, len(raw), 512):
        n, = struct.unpack(">H", raw[off + 30:off + 32])
        data, = struct.unpack(">H", raw[off + 44:off + 46])
        b, = struct.unpack(">H", raw[off + 46:off + 48])
        while b:
            typ, nxt = struct.unpack(">HH", raw[off + b:off + b + 4])
            if typ == 1000:
                raw[off + b + 5] = 0
            b = nxt
        lo, hi = off + data, off + data + n * width
        raw[lo:hi] = struct.pack(f"<{n}{fmt}", *struct.unpack(f">{n}{fmt}", raw[lo:hi]))
    path.write_bytes(raw)


def git_state():
    root = Path(__file__).resolve().parent.parent
    try:
        commit = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", str(root), "status", "--porcelain"],
                                    capture_output=True, text=True, check=True).stdout.strip())
    except (OSError, subprocess.CalledProcessError):
        commit, dirty = None, None
    return commit, dirty


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20251110)
    ap.add_argument("--stretch", type=int, default=200, help="samples per stretch")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    files, ok = {}, True
    for name, (encoding, order, bounds) in FIXTURES.items():
        limit = 32767 if encoding == "INT16" else 2**30
        x = walk(rng, bounds, a.stretch, limit)
        dtype = np.int16 if encoding == "INT16" else np.int32
        tr = Trace(data=x.astype(dtype), header={
            "network": "XX", "station": "FIXT", "location": "", "channel": "HHZ",
            "sampling_rate": 100.0, "starttime": START})
        path = out / f"{name}.mseed"
        tr.write(str(path), format="MSEED", encoding=encoding, reclen=512, byteorder=">")
        if order == "<":
            swap_data(path, 2 if encoding == "INT16" else 4)
        st = obspy.read(str(path))
        ref = np.concatenate([t.data.astype(np.int32) for t in st])
        if len(st) != 1 or not np.array_equal(ref, x.astype(np.int32)):
            print(f"{name}: does not read back as written", file=sys.stderr)
            ok = False
        ref.astype("<i4").tofile(out / f"{name}.i32")
        size = path.stat().st_size
        if size == 0 or size % 512 or len(ref) == 0:
            print(f"{name}: empty or partial file", file=sys.stderr)
            ok = False
        files[name] = {
            "encoding": encoding, "header": "big-endian",
            "data": "big-endian" if order == ">" else "little-endian",
            "difference_bounds": bounds, "records": size // 512, "samples": int(len(ref)),
            "start": str(st[0].stats.starttime), "sha256_mseed": hashlib.sha256(path.read_bytes()).hexdigest(),
            "sha256_i32": hashlib.sha256((out / f"{name}.i32").read_bytes()).hexdigest(),
        }
        print(f"{name:9} {encoding:7} {files[name]['data']:13} "
              f"{files[name]['records']:3} records {len(ref):5} samples")

    commit, dirty = git_state()
    record = {
        "command": "tools/make_mseed_fixtures.py",
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "argv": sys.argv[1:],
        "parameters": {"seed": a.seed, "stretch": a.stretch, "start": str(START),
                       "id": "XX.FIXT..HHZ", "sampling_rate": 100.0, "reclen": 512},
        "git_commit": commit, "git_dirty": dirty,
        "versions": {"python": platform.python_version(), "obspy": obspy.__version__,
                     "numpy": np.__version__},
        "files": files,
    }
    (out / "fixtures.json").write_text(json.dumps(record, indent=1) + "\n")
    sys.exit(0 if ok and files else 1)


if __name__ == "__main__":
    main()
