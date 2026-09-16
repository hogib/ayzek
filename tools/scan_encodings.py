"""Counts the data encodings of all miniSEED records in archive chunk files.

For each zip, reports the number of records per blockette-1000 encoding and
byte order, and the times of the first and last non-Steim2 records. Used to
find how often dataloggers write uncompressed integer records. Standard
library only.

    python3 tools/scan_encodings.py --workers 8 --out encodings.csv ~/Projects/sismokaos/tdvms/afad_raw/*/*.zip
"""
import argparse
import collections
import csv
import datetime as dt
import multiprocessing
import struct
import zipfile
from pathlib import Path

REC = 512
NAMES = {1: "int16", 3: "int32", 10: "steim1", 11: "steim2"}


def record_time(rec):
    year, doy, hour, minute, sec, _, frac = struct.unpack(">HHBBBBH", rec[20:30])
    t = dt.datetime(year, 1, 1, tzinfo=dt.UTC) + dt.timedelta(
        days=doy - 1, hours=hour, minutes=minute, seconds=sec, microseconds=frac * 100)
    return t


def scan(path):
    counts = collections.Counter()
    first = last = None
    try:
        with zipfile.ZipFile(path) as zf:
            member = next(n for n in zf.namelist() if n.lower().endswith(".mseed"))
            with zf.open(member) as src:
                while (rec := src.read(REC)) and len(rec) == REC:
                    off = struct.unpack(">H", rec[46:48])[0]
                    key = "no b1000"
                    for _ in range(16):
                        if off == 0 or off + 8 > REC:
                            break
                        typ, nxt = struct.unpack(">HH", rec[off:off + 4])
                        if typ == 1000:
                            enc, order = rec[off + 4], rec[off + 5]
                            key = f"{NAMES.get(enc, f'enc{enc}')}-{'be' if order else 'le'}"
                            break
                        off = nxt
                    counts[key] += 1
                    if not key.startswith("steim2"):
                        t = record_time(rec)
                        first = first or t
                        last = t
    except (zipfile.BadZipFile, StopIteration, OSError) as e:
        return {"file": str(path), "error": type(e).__name__}
    return {"file": str(path), "counts": dict(counts), "first": first, "last": last}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", required=True)
    ap.add_argument("zips", nargs="+")
    a = ap.parse_args()
    zips = [z for z in a.zips if ".dup" not in z and Path(z).name.count("_") == 1]
    with multiprocessing.Pool(a.workers) as pool, open(a.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["station", "chunk", "records", "non_steim2", "encodings", "first_non_steim2", "last_non_steim2"])
        for i, r in enumerate(pool.imap_unordered(scan, zips), 1):
            name = Path(r["file"]).stem
            if "error" in r:
                w.writerow([name.split("_")[0], name, "", "", r["error"], "", ""])
                continue
            total = sum(r["counts"].values())
            other = total - r["counts"].get("steim2-be", 0)
            w.writerow([name.split("_")[0], name, total, other,
                        " ".join(f"{k}:{v}" for k, v in sorted(r["counts"].items())),
                        r["first"].isoformat() if r["first"] else "", r["last"].isoformat() if r["last"] else ""])
            f.flush()
            if other:
                print(f"[{i}/{len(zips)}] {name}: {other} of {total} records not Steim2 "
                      f"({r['first']:%Y-%m-%d %H:%M} .. {r['last']:%Y-%m-%d %H:%M})", flush=True)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
