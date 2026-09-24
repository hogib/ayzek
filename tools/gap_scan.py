#!/usr/bin/env python3
"""Report the data gaps in miniSEED files, per channel.

    tools/gap_scan.py data/noise_a1/*.mseed
    tools/gap_scan.py --min-seconds 60 --csv gaps.csv data/marmara/*.mseed

For each file it prints the span of every channel, the number of gaps, and the
largest one in samples and seconds, with the ring capacity for comparison: a
gap longer than the ring is the case that used to deadlock ingest
(`src/pipeline/writer.hpp`). Only 512-byte miniSEED-2 record headers are read,
so a multi-gigabyte file costs a pass over 1/128th of its bytes and nothing is
decoded.

No dependencies beyond the standard library, on purpose: this is the tool that
checks the data the pipeline is about to be blamed for.
"""

import argparse
import csv
import datetime as dt
import glob
import struct
import sys
from pathlib import Path

RECORD = 512
RING_SAMPLES = 1 << 17  # kRingCapacity, src/pipeline/station.hpp


def btime(raw):
    """SEED BTIME -> POSIX seconds."""
    year, day, hour, minute, sec, _, frac = struct.unpack(">HHBBBBH", raw)
    t = dt.datetime(year, 1, 1, tzinfo=dt.timezone.utc) + dt.timedelta(
        days=day - 1, hours=hour, minutes=minute, seconds=sec, microseconds=frac * 100
    )
    return t.timestamp()


def sample_rate(factor, multiplier):
    """SEED sample-rate factor/multiplier -> Hz (the two usual encodings)."""
    if factor > 0 and multiplier > 0:
        return float(factor * multiplier)
    if factor > 0 and multiplier < 0:
        return factor / -multiplier
    if factor < 0 and multiplier > 0:
        return -multiplier / factor
    if factor < 0 and multiplier < 0:
        return 1.0 / (factor * multiplier)
    return 0.0


def records(path):
    """Yields (station, channel, start, samples, rate) per record header."""
    with open(path, "rb") as f:
        while True:
            raw = f.read(RECORD)
            if len(raw) < RECORD:
                return
            sta = raw[8:13].decode("ascii", "replace").strip()
            cha = raw[15:18].decode("ascii", "replace").strip()
            nsamp = struct.unpack(">H", raw[30:32])[0]
            rate = sample_rate(*struct.unpack(">hh", raw[32:36]))
            if rate <= 0 or nsamp == 0:
                continue
            yield sta, cha, btime(raw[20:30]), nsamp, rate


def utc(t):
    return dt.datetime.fromtimestamp(t, dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--min-seconds", type=float, default=1.0,
                    help="report gaps at least this long (default 1)")
    ap.add_argument("--csv", help="also write every reported gap here")
    args = ap.parse_args()

    paths = sorted({p for pattern in args.files for p in glob.glob(pattern)})
    if not paths:
        sys.exit("no files matched")

    rows = []
    worst = 0.0
    for path in paths:
        chans = {}
        station = ""
        for sta, cha, start, nsamp, rate in records(path):
            station = station or sta
            chans.setdefault(cha, []).append((start, nsamp, rate))
        if not chans:
            print(f"{Path(path).name}: no readable record headers")
            continue
        print(f"{Path(path).name} ({station})")
        for cha, recs in sorted(chans.items()):
            recs.sort()
            gaps = []
            for (t0, n, rate), (t1, _, _) in zip(recs, recs[1:]):
                seconds = t1 - (t0 + n / rate)
                if seconds >= args.min_seconds:
                    gaps.append((seconds, t0 + n / rate, rate))
            end = recs[-1][0] + recs[-1][1] / recs[-1][2]
            longest = max(gaps, default=(0.0, None, recs[0][2]))
            worst = max(worst, longest[0])
            samples = round(longest[0] * longest[2])
            flag = "  <-- longer than the ring" if samples > RING_SAMPLES else ""
            print(f"  {cha}  {utc(recs[0][0])} .. {utc(end)}  "
                  f"{len(recs):>7} records  {len(gaps):>3} gaps  "
                  f"longest {longest[0]:8.2f} s ({samples:>7} samples)"
                  f"{'' if longest[1] is None else ' at ' + utc(longest[1])}{flag}")
            for seconds, at, rate in gaps:
                rows.append({"file": Path(path).name, "station": station,
                             "channel": cha, "start": utc(at),
                             "seconds": round(seconds, 3),
                             "samples": round(seconds * rate)})

    print(f"\nring capacity {RING_SAMPLES} samples "
          f"({RING_SAMPLES / 100:.1f} s at 100 Hz); "
          f"longest gap seen {worst:.2f} s ({round(worst * 100)} samples at 100 Hz)")
    if args.csv:
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["file", "station", "channel",
                                              "start", "seconds", "samples"])
            w.writeheader()
            w.writerows(rows)
        print(f"{len(rows)} gaps written to {args.csv}")


if __name__ == "__main__":
    main()
