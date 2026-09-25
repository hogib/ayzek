#!/usr/bin/env python3
"""Cut catalogue-quiet windows out of the TDVMS archive, for false-alarm runs.

    tools/extract_quiet.py --archive ~/…/afad_raw --catalog CAT.csv \
        --stations-csv tests/stations.csv --hours 100 --out data/quiet

A window is quiet when no catalogued event lies within `--radius` km of any of
the chosen stations for `--guard` minutes either side of it. The first
`--lead-in` minutes of each window are extracted but not scored: they warm the
STA/LTA and the noise baseline, and because they are taken from inside the
quiet window itself, every scored second is still `--guard` minutes clear of
any catalogued event.

Why it works this way:

- **Records are copied, not decoded.** The archive holds one 512-byte-record
  miniSEED file per station-chunk; a record is kept when it overlaps the
  window. Steim2 stays Steim2 and nothing is re-encoded.
- **One pass per station-chunk.** Chunks are 21 days and ~900 MB, so the read
  dominates. Windows are grouped by chunk and the chunks holding the most
  quiet time are taken first, which is why `--hours` is a target rather than
  an exact figure.
- **Each window is its own replay.** They are minutes-to-hours apart, so
  joining them would invent gaps of days; `run.sh` replays each separately and
  the scorecard sums their hours.

Writes `dataset.json` per window and `quiet_index.csv` over all of them, and
exits non-zero if nothing was produced.
"""

import argparse
import bisect
import csv
import datetime as dt
import json
import math
import re
import struct
import sys
import zipfile
from pathlib import Path

RECORD = 512


def hav(a, b, c, d):
    r = 6371.0
    p = math.radians
    x = (math.sin(p(c - a) / 2) ** 2
         + math.cos(p(a)) * math.cos(p(c)) * math.sin(p(d - b) / 2) ** 2)
    return 2 * r * math.asin(math.sqrt(x))


def btime(raw):
    year, day, hour, minute, sec, _, frac = struct.unpack(">HHBBBBH", raw)
    return (dt.datetime(year, 1, 1, tzinfo=dt.timezone.utc)
            + dt.timedelta(days=day - 1, hours=hour, minutes=minute,
                           seconds=sec, microseconds=frac * 100)).timestamp()


def sample_rate(factor, mult):
    if factor > 0 and mult > 0:
        return float(factor * mult)
    if factor > 0 and mult < 0:
        return factor / -mult
    if factor < 0 and mult > 0:
        return -mult / factor
    if factor < 0 and mult < 0:
        return 1.0 / (factor * mult)
    return 0.0


def load_stations(path, wanted):
    out = {}
    with open(path, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            code = (r.get("station") or r.get("Code") or r.get("code") or "").strip()
            if code and (not wanted or code in wanted):
                lat = r.get("latitude") or r.get("Latitude") or r.get("lat")
                lon = r.get("longitude") or r.get("Longitude") or r.get("lon")
                out.setdefault(code, (float(lat), float(lon)))
    return out


def load_events(path, stations, radius):
    """Origin times of catalogued events within `radius` km of any station."""
    out = []
    with open(path, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            raw = (r.get("Date") or r.get("date") or "")[:19]
            for fmt in ("%d/%m/%Y %H:%M:%S", "%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S"):
                try:
                    t = dt.datetime.strptime(raw, fmt).replace(tzinfo=dt.timezone.utc)
                    break
                except ValueError:
                    t = None
            if t is None:
                continue
            try:
                la, lo = float(r["Latitude"]), float(r["Longitude"])
            except (KeyError, ValueError):
                continue
            if any(hav(la, lo, a, b) < radius for a, b in stations.values()):
                out.append(t.timestamp())
    out.sort()
    return out


def quiet_windows(events, lo, hi, guard_s, min_len_s):
    free, cur = [], lo
    for t in events:
        if t - guard_s > cur:
            free.append((cur, t - guard_s))
        cur = max(cur, t + guard_s)
    if hi > cur:
        free.append((cur, hi))
    return [(a, b) for a, b in free if b - a >= min_len_s]


def chunk_spans(archive, stations):
    """{station: [(start, end, zip path)]} from the chunk file names."""
    out = {}
    for st in stations:
        for z in sorted((Path(archive) / st).glob(f"{st}_*.zip")):
            m = re.search(r"_(\d{4})-(\d{2})-(\d{2})\.zip$", z.name)
            if not m:
                continue
            start = dt.datetime(*map(int, m.groups()), tzinfo=dt.timezone.utc)
            out.setdefault(st, []).append(
                (start.timestamp(), (start + dt.timedelta(days=21)).timestamp(), z))
    return out


def extract(zpath, windows, writers):
    """Copies every record of `zpath` overlapping any window to its writer."""
    kept = 0
    starts = [w[0] for w in windows]
    with zipfile.ZipFile(zpath) as z:
        name = z.namelist()[0]
        with z.open(name) as f:
            while True:
                raw = f.read(RECORD)
                if len(raw) < RECORD:
                    break
                n = struct.unpack(">H", raw[30:32])[0]
                fs = sample_rate(*struct.unpack(">hh", raw[32:36]))
                if not n or fs <= 0:
                    continue
                t0 = btime(raw[20:30])
                t1 = t0 + n / fs
                # Windows are disjoint and sorted, so one lookup finds the
                # only one a record can fall in.
                i = bisect.bisect_right(starts, t1) - 1
                if i >= 0 and t0 < windows[i][1] and t1 > windows[i][0]:
                    writers[i].write(raw)
                    kept += 1
    return kept


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--archive", required=True, help="afad_raw root, one directory per station")
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--stations-csv", required=True)
    ap.add_argument("--stations", required=True, help="comma-separated station codes")
    ap.add_argument("--out", required=True)
    ap.add_argument("--hours", type=float, default=100.0, help="target scored hours")
    ap.add_argument("--guard", type=float, default=30.0, help="minutes clear of any event")
    ap.add_argument("--lead-in", type=float, default=30.0, help="unscored warm-up minutes")
    ap.add_argument("--min-window", type=float, default=4.0, help="hours, before the lead-in")
    ap.add_argument("--radius", type=float, default=250.0, help="km from any station")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    want = [s.strip() for s in a.stations.split(",") if s.strip()]
    stations = load_stations(a.stations_csv, set(want))
    missing = [s for s in want if s not in stations]
    if missing:
        sys.exit(f"no coordinates for {', '.join(missing)}")
    chunks = chunk_spans(a.archive, want)
    absent = [s for s in want if s not in chunks]
    if absent:
        sys.exit(f"no archive chunks for {', '.join(absent)}")

    # Only where every station has a chunk.
    lo = max(min(s for s, _, _ in v) for v in chunks.values())
    hi = min(max(e for _, e, _ in v) for v in chunks.values())
    print(f"{len(want)} stations, common archive span "
          f"{dt.datetime.fromtimestamp(lo, dt.timezone.utc):%Y-%m-%d} to "
          f"{dt.datetime.fromtimestamp(hi, dt.timezone.utc):%Y-%m-%d}")

    # A window is only usable when EVERY station has an archive chunk over it:
    # a false-alarm rate measured on two stations is not the same quantity as
    # one measured on eight, and the archive is missing whole chunk periods at
    # several stations.
    def covered(w):
        return all(any(s <= w[0] and w[1] <= e for s, e, _ in chunks[st]) for st in want)

    events = load_events(a.catalog, stations, a.radius)
    print(f"{len(events):,} catalogued events within {a.radius:.0f} km of any station")
    lead = a.lead_in * 60
    wins = quiet_windows(events, lo, hi, a.guard * 60, a.min_window * 3600 + lead)
    before = len(wins)
    wins = [w for w in wins if covered(w)]
    print(f"{before - len(wins)} window(s) dropped: not every station has archive over them")
    print(f"{len(wins)} quiet windows of at least {a.min_window} h "
          f"({sum(b - x for x, b in wins) / 3600:,.0f} h total, "
          f"{sum(b - x - lead for x, b in wins) / 3600:,.0f} h scored)")

    # Spread the selection across the archive rather than packing it into the
    # fewest chunks: a false-alarm rate is a property of the noise, which
    # varies with season, weather and cultural activity, so sixteen windows
    # from one fortnight would measure that fortnight. Reading a chunk costs
    # about five seconds here (headers only, no decoding), so there is nothing
    # to save by clustering. The span is cut into as many bins as windows are
    # needed and the longest window in each bin is taken.
    lead_h = lead / 3600
    need = max(1, int(a.hours / max(1e-9, (sum(b - x for x, b in wins) / len(wins)) / 3600 - lead_h)))
    picked, scored = [], 0.0
    edges = [lo + (hi - lo) * i / need for i in range(need + 1)]
    bins = [[w for w in wins if edges[i] <= w[0] < edges[i + 1]] for i in range(need)]
    for b in bins:
        if not b:
            continue
        w = max(b, key=lambda w: w[1] - w[0])
        picked.append(w)
        scored += w[1] - w[0] - lead
    # Top up from the longest windows still unused if the bins fell short.
    for w in sorted(wins, key=lambda w: -(w[1] - w[0])):
        if scored >= a.hours * 3600:
            break
        if w not in picked:
            picked.append(w)
            scored += w[1] - w[0] - lead
    picked.sort()
    print(f"selected {len(picked)} windows -> {scored / 3600:.1f} scored hours")
    for i, (x, b) in enumerate(picked):
        print(f"  {i:>3} {dt.datetime.fromtimestamp(x, dt.timezone.utc):%Y-%m-%d %H:%M} "
              f"+{(b - x) / 3600:.1f} h  (scored {(b - x - lead) / 3600:.1f} h)")
    if a.dry_run:
        return
    if not picked:
        sys.exit("no windows selected")

    out = Path(a.out)
    names = [f"quiet_{dt.datetime.fromtimestamp(x, dt.timezone.utc):%Y%m%d_%H%M}"
             for x, _ in picked]
    for nm in names:
        (out / nm).mkdir(parents=True, exist_ok=True)

    total = 0
    for st in want:
        for s, e, z in chunks[st]:
            idx = [i for i, w in enumerate(picked) if s <= w[0] and w[1] <= e]
            if not idx:
                continue
            idx.sort(key=lambda i: picked[i][0])
            paths = [out / names[i] / f"{st}.mseed" for i in idx]
            handles = [open(p, "wb") for p in paths]
            try:
                n = extract(z, [picked[i] for i in idx], handles)
            finally:
                for h in handles:
                    h.close()
            total += n
            print(f"  {z.name}: {n:,} records -> {len(idx)} window(s)")
    if total == 0:
        sys.exit("no records extracted; refusing to leave an empty dataset")

    rows = []
    for nm, (x, b) in zip(names, picked):
        d = out / nm
        sizes = {p.stem: p.stat().st_size for p in sorted(d.glob("*.mseed"))}
        empty = [k for k, v in sizes.items() if v == 0]
        for k in empty:
            (d / f"{k}.mseed").unlink()
        meta = {
            "generator": "tools/extract_quiet.py",
            "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "argv": sys.argv[1:],
            "window_start_utc": dt.datetime.fromtimestamp(x, dt.timezone.utc).isoformat(),
            "window_end_utc": dt.datetime.fromtimestamp(b, dt.timezone.utc).isoformat(),
            "score_from_utc": dt.datetime.fromtimestamp(x + lead, dt.timezone.utc).isoformat(),
            "scored_hours": round((b - x - lead) / 3600, 4),
            "lead_in_minutes": a.lead_in,
            "guard_minutes": a.guard,
            "radius_km": a.radius,
            "stations": sorted(k for k in sizes if k not in empty),
            "bytes": {k: v for k, v in sizes.items() if v},
        }
        (d / "dataset.json").write_text(json.dumps(meta, indent=2) + "\n")
        rows.append({"name": nm, "start_utc": meta["window_start_utc"],
                     "score_from_utc": meta["score_from_utc"],
                     "end_utc": meta["window_end_utc"],
                     "scored_hours": meta["scored_hours"],
                     "stations": len(meta["stations"])})
    with open(out / "quiet_index.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\n{len(rows)} windows, {sum(r['scored_hours'] for r in rows):.1f} scored hours "
          f"-> {out}/quiet_index.csv")


if __name__ == "__main__":
    main()
