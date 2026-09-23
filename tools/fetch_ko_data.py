"""Pulls continuous KO (KOERI) waveforms into the replay layout ayzek reads.

The demo and evaluation datasets this repository ships are AFAD recordings,
which cannot be redistributed. KOERI's network KO is open data under citation
(doi:10.7914/SN/KO), so the published datasets are built from it instead. The
detector and the magnitude regressor were both trained on KO waveforms, so this
is also the domain they were fitted on.

    # which stations actually have the span, before downloading it
    uv run --project tools python tools/fetch_ko_data.py probe --set demo_ko

    # pull what the probe chose
    uv run --project tools python tools/fetch_ko_data.py pull --set demo_ko

`probe` asks the archive for one minute per candidate station, ranked by
distance from the target, and keeps those that return three components. It
writes `data/<set>/stations.json`, which `pull` then reads: the two steps are
separate because the archive is patchy and a station list is worth inspecting
before spending an hour of requests on it.

`pull` writes one file per station, `data/<set>/<STATION>.mseed`, hour by hour
so that a failure costs one hour and not the whole span, plus a `dataset.json`
recording the request, the service, the station list and what arrived.

Run in the tools environment (obspy).
"""
import argparse
import datetime as dt
import json
import sys
from pathlib import Path

from obspy import UTCDateTime, Stream, read
from obspy.clients.fdsn import Client
from obspy.geodetics import locations2degrees

ROOT = Path(__file__).resolve().parents[1]
SERVICE = "KOERI"
NETWORK = "KO"
CHANNEL = "HH?"

# set: (start, end, centre lat, centre lon, max distance km, how many stations)
SETS = {
    # The 2025-11-10 Sındırgı Mw 4.9, the event the README and the end-to-end
    # test use. Same span as the AFAD demo, so the numbers stay comparable.
    "demo_ko": ("2025-11-10T18:04:00", "2025-11-10T18:35:00",
                39.22333, 28.16556, 150.0, 6),
    # The 2025-04-23 Marmara Mw 6.2, for the large-event and magnitude story.
    # Starts 49 min before the Mw 6.2: the magnitude regressor needs the
    # station noise baselines warm (30 windows over 5 min, §3.3.5 of the paper),
    # and a 4-minute lead-in cost 0.7 magnitude units on the mainshock.
    "marmara_ko": ("2025-04-23T09:00:00", "2025-04-23T10:45:00",
                   40.87, 28.20, 200.0, 8),
    # A catalogue-quiet span on the demo stations; `quiet-scan` chooses it.
    # 2024-01-03 11:14 to 2024-01-04 02:33 holds no catalogue event within
    # 250 km of these stations (`quiet-scan`); this window sits inside it with
    # margin at both ends. Triggering starts 10 min in, so the first windows
    # build each station's noise baseline.
    "quiet_ko": ("2024-01-03T11:25:00", "2024-01-03T17:25:00",
                 39.22333, 28.16556, 150.0, 6),
}


def client():
    return Client(SERVICE, timeout=300)


def candidates(name):
    """KO stations with HH channels active over the span, nearest first."""
    t0, t1, lat, lon, max_km, _ = SETS[name]
    inv = client().get_stations(network=NETWORK, channel=CHANNEL, level="station",
                                starttime=UTCDateTime(t0), endtime=UTCDateTime(t1))
    rows = []
    for net in inv:
        for s in net:
            km = locations2degrees(lat, lon, s.latitude, s.longitude) * 111.195
            if km <= max_km:
                rows.append({"station": s.code, "lat": float(s.latitude),
                             "lon": float(s.longitude),
                             "elevation_m": float(s.elevation), "distance_km": round(km, 1)})
    rows.sort(key=lambda r: r["distance_km"])
    return rows


def cmd_probe(args):
    """One minute per candidate; keeps the stations that return three components."""
    name = args.set
    t0, t1, *_rest, want = SETS[name]
    rows = candidates(name)
    print(f"[probe] {len(rows)} KO stations with {CHANNEL} within "
          f"{SETS[name][4]:.0f} km, active {t0} to {t1}")
    c = client()
    keep, probed = [], []
    for r in rows:
        if len(keep) >= want and not args.all:
            break
        try:
            st = c.get_waveforms(NETWORK, r["station"], "*", CHANNEL,
                                 UTCDateTime(t0), UTCDateTime(t0) + 60)
        except Exception as e:
            r["probe"] = type(e).__name__
            probed.append(r)
            print(f"  {r['station']:6s} {r['distance_km']:6.1f} km  -- {r['probe']}")
            continue
        comps = sorted({tr.stats.channel[-1] for tr in st})
        rates = sorted({tr.stats.sampling_rate for tr in st})
        r["probe"] = f"{len(st)} traces {''.join(comps)} at {'/'.join(f'{x:g}' for x in rates)} Hz"
        r["ok"] = len(comps) == 3 and rates == [100.0]
        probed.append(r)
        print(f"  {r['station']:6s} {r['distance_km']:6.1f} km  {r['probe']}"
              f"{'' if r['ok'] else '   SKIPPED'}")
        if r["ok"]:
            keep.append(r)
    out = ROOT / "data" / name
    out.mkdir(parents=True, exist_ok=True)
    (out / "stations.json").write_text(json.dumps(
        {"set": name, "span": [t0, t1], "service": SERVICE, "network": NETWORK,
         "probed": probed, "chosen": keep}, indent=1) + "\n")
    print(f"[probe] {len(keep)} usable -> {out / 'stations.json'}")
    if keep:
        hours = (UTCDateTime(t1) - UTCDateTime(t0)) / 3600
        print(f"[probe] pull estimate: {len(keep)} stations x {hours:.1f} h "
              f"~ {1.7 * len(keep) * hours:.0f} MB")
    return 0


def cmd_pull(args):
    """Downloads the chosen stations hour by hour and writes one file each."""
    name = args.set
    out = ROOT / "data" / name
    spec = json.loads((out / "stations.json").read_text())
    t0, t1 = (UTCDateTime(x) for x in spec["span"])
    c = client()
    got = []
    for r in spec["chosen"]:
        dest = out / f"{r['station']}.mseed"
        if dest.exists() and not args.force:
            print(f"  {r['station']}: present")
            continue
        st = Stream()
        t = t0
        while t < t1:
            end = min(t + 3600, t1)
            try:
                st += c.get_waveforms(NETWORK, r["station"], "*", CHANNEL, t, end)
            except Exception as e:
                print(f"  {r['station']} {t.strftime('%H:%M')}-{end.strftime('%H:%M')}: "
                      f"{type(e).__name__}")
            t = end
        if not len(st):
            print(f"  {r['station']}: nothing")
            continue
        # Merge touching segments, then split again: `merge` alone returns a
        # masked array across gaps, which miniSEED cannot store. Splitting
        # leaves one trace per continuous stretch, which is what the archive
        # files look like and what the replay expects -- gaps stay gaps.
        st.merge(method=1, fill_value=None)
        st = st.split()
        st.write(str(dest), format="MSEED", encoding="STEIM2")
        comps = sorted({tr.stats.channel[-1] for tr in st})
        secs = sum(tr.stats.npts / tr.stats.sampling_rate for tr in st) / max(1, len(comps))
        got.append({**r, "file": dest.name, "seconds_per_component": round(secs),
                    "megabytes": round(dest.stat().st_size / 1e6, 1)})
        print(f"  {r['station']}: {''.join(comps)} {secs / 3600:.2f} h "
              f"-> {dest.name} ({got[-1]['megabytes']} MB)")
    (out / "dataset.json").write_text(json.dumps(
        {"set": name, "span": spec["span"], "service": SERVICE, "network": NETWORK,
         "channel": CHANNEL, "retrieved_utc": dt.datetime.now(dt.UTC).isoformat(),
         "citation": "Kandilli Observatory and Earthquake Research Institute (KOERI), "
                     "Bogazici University. International Federation of Digital "
                     "Seismograph Networks. doi:10.7914/SN/KO",
         "stations": got}, indent=1) + "\n")
    print(f"[pull] {len(got)} stations -> {out}")
    return 0


def cmd_quiet_scan(args):
    """Spans with no catalogue event within 250 km of a set's stations."""
    import csv
    import math
    spec = json.loads((ROOT / "data" / args.stations_from / "stations.json").read_text())
    stations = [(r["lat"], r["lon"]) for r in spec["chosen"]]
    if not stations:
        sys.exit(f"no chosen stations in {args.stations_from}")

    def km(a, b, c, d):
        p, q = math.radians(a), math.radians(c)
        return 6371 * 2 * math.asin(math.sqrt(
            math.sin((q - p) / 2) ** 2 + math.cos(p) * math.cos(q)
            * math.sin((math.radians(d) - math.radians(b)) / 2) ** 2))

    lo, hi = dt.datetime.fromisoformat(args.start), dt.datetime.fromisoformat(args.end)
    near = []
    with open(args.catalog, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            try:
                t = dt.datetime.strptime(row["Date"], "%d/%m/%Y %H:%M:%S")
            except (ValueError, KeyError):
                continue
            if not (lo <= t < hi):
                continue
            la, ln = float(row["Latitude"]), float(row["Longitude"])
            if min(km(la, ln, s[0], s[1]) for s in stations) <= args.radius_km:
                near.append(t)
    near.sort()
    edges = [lo] + near + [hi]
    gaps = sorted(((edges[i + 1] - edges[i], edges[i]) for i in range(len(edges) - 1)),
                  reverse=True)[:args.top]
    print(f"[quiet] {len(near)} events within {args.radius_km:g} km of "
          f"{len(stations)} stations, {lo:%F} to {hi:%F}")
    for length, start in gaps:
        print(f"  {start:%F %T} + {length.total_seconds() / 3600:5.1f} h")
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("probe", help="which stations have the span")
    q.add_argument("--set", required=True, choices=sorted(SETS))
    q.add_argument("--all", action="store_true", help="probe every candidate, not just enough")
    q.set_defaults(fn=cmd_probe)
    f = sub.add_parser("pull", help="download the chosen stations")
    f.add_argument("--set", required=True, choices=sorted(SETS))
    f.add_argument("--force", action="store_true", help="re-download files already present")
    f.set_defaults(fn=cmd_pull)
    s = sub.add_parser("quiet-scan", help="find catalogue-quiet spans for a set's stations")
    s.add_argument("--stations-from", required=True, help="set whose stations.json to use")
    s.add_argument("--catalog", required=True)
    s.add_argument("--start", required=True)
    s.add_argument("--end", required=True)
    s.add_argument("--radius-km", type=float, default=250.0)
    s.add_argument("--top", type=int, default=10)
    s.set_defaults(fn=cmd_quiet_scan)
    a = p.parse_args()
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
