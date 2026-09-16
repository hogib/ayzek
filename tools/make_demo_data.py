"""Extracts a time range from archive chunks into one miniSEED file per station.

Records whose start time lies in [start - 60 s, end) are copied unchanged, in
file order. Only HH? channels are kept. Uses the standard library only.

    python3 tools/make_demo_data.py --out data/demo \
        --start 2025-11-10T18:05:00 --end 2025-11-10T18:35:00 \
        ~/Projects/sismokaos/tdvms/afad_raw/{DEMI,MANT,BAND}/*_2025-10-29.zip

Chunk files are named by their start date.
"""
import argparse
import datetime as dt
import struct
import sys
import zipfile
from pathlib import Path

REC = 512


def record_time(rec):
    year, doy, hour, minute, sec, _, frac = struct.unpack(">HHBBBBH", rec[20:30])
    t = dt.datetime(year, 1, 1, tzinfo=dt.UTC) + dt.timedelta(
        days=doy - 1, hours=hour, minutes=minute, seconds=sec, microseconds=frac * 100)
    return t.timestamp()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--start", required=True, help="UTC, ISO format")
    ap.add_argument("--end", required=True, help="UTC, ISO format")
    ap.add_argument("zips", nargs="+")
    a = ap.parse_args()
    t0 = dt.datetime.fromisoformat(a.start).replace(tzinfo=dt.UTC).timestamp()
    t1 = dt.datetime.fromisoformat(a.end).replace(tzinfo=dt.UTC).timestamp()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    for z in a.zips:
        with zipfile.ZipFile(z) as zf:
            member = next(n for n in zf.namelist() if n.lower().endswith(".mseed"))
            kept = seen = 0
            station = None
            with zf.open(member) as src, open(out / "tmp.mseed", "wb") as dst:
                while rec := src.read(REC):
                    if len(rec) < REC:
                        break
                    seen += 1
                    if rec[15:17] != b"HH":
                        continue
                    # Records are ~4 s long; the 60 s margin keeps complete records at the start.
                    if t0 - 60 <= record_time(rec) < t1:
                        dst.write(rec)
                        kept += 1
                        station = station or rec[8:13].decode().strip()
        if not kept:
            print(f"{z}: no records in range", file=sys.stderr)
            (out / "tmp.mseed").unlink()
            continue
        path = out / f"{station}.mseed"
        (out / "tmp.mseed").rename(path)
        print(f"{path}: {kept} of {seen} records")


if __name__ == "__main__":
    main()
