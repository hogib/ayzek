#!/usr/bin/env python3
"""Look up network alarms in external earthquake catalogues (EMSC, ISC).

For each alarm of a classes file (the paper's figures/supp/quiet_alarms.csv:
alarm, window, alarm_utc, declared_by), the FDSN event services are asked for
events with origin from --before seconds before the alarm to the alarm, within
--radius km of the midpoint of the first two declaring stations, any magnitude
(ayzek_paper revisions2.md R17.0). Every event found is listed with its
distances to the two stations; it counts as a candidate source when the P
travel time from it (--vp) puts P within --tol seconds of each station's
reported P (the latest detection in the replay's .rec file in the 120 s before
the alarm).

    uv run --with obspy python tools/alarm_catalogues.py \\
        --alarms ../ayzek_paper/figures/supp/quiet_alarms.csv \\
        --runs data/runs/quiet100 --out data/runs/r5/alarm_catalogues.csv
"""
import argparse
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from misfire_review import hav, read_record  # noqa: E402

# KOERI's FDSN node serves waveforms but no events; the ISC bulletin carries
# KOERI's own solutions (agency ISK) with every other agency's, so it stands in.
SERVICES = {"EMSC": "https://www.seismicportal.eu",
            "ISC": "http://www.isc.ac.uk"}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--alarms", required=True)
    ap.add_argument("--runs", required=True, help="directory of the replays' .rec files")
    ap.add_argument("--before", type=float, default=180.0)
    ap.add_argument("--radius", type=float, default=300.0)
    ap.add_argument("--vp", type=float, default=6.0)
    ap.add_argument("--tol", type=float, default=10.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    from obspy import UTCDateTime
    from obspy.clients.fdsn import Client
    from obspy.clients.fdsn.header import FDSNNoDataException

    clients = {}
    for name, url in SERVICES.items():
        try:
            clients[name] = Client(url, timeout=60)
        except Exception as e:  # a service that is down is reported, not fatal
            print(f"{name}: {e}", file=sys.stderr)

    rows = []
    for al in csv.DictReader(open(a.alarms)):
        t = UTCDateTime(al["alarm_utc"].replace(" ", "T"))
        stations, det = read_record(f"{a.runs}/{al['window']}.rec")[:2]
        A, B = al["declared_by"].split()[:2]
        refp = {}
        for s in (A, B):
            before = [x for x in det.get(s, []) if t.timestamp - 120 <= x <= t.timestamp]
            refp[s] = before[-1] if before else math.nan
        lat = (stations[A][0] + stations[B][0]) / 2
        lon = (stations[A][1] + stations[B][1]) / 2
        found = 0
        for name, cl in clients.items():
            try:
                kw = {"includeallorigins": True} if name == "ISC" else {}
                cat = cl.get_events(starttime=t - a.before, endtime=t, latitude=lat, longitude=lon,
                                    maxradius=a.radius / 111.19, **kw)
            except FDSNNoDataException:
                continue
            except Exception as e:
                print(f"alarm {al['alarm']} {name}: {e}", file=sys.stderr)
                rows.append({"alarm": al["alarm"], "catalogue": name, "origin_utc": "", "magnitude": "",
                             "dist_a_km": "", "dist_b_km": "", "p_residual_a_s": "", "p_residual_b_s": "",
                             "candidate": f"query failed: {e}"[:80], "agencies": ""})
                continue
            for ev in cat:
                o = ev.preferred_origin() or ev.origins[0]
                agencies = sorted({x.creation_info.agency_id for x in ev.origins
                                   if x.creation_info and x.creation_info.agency_id})
                m = ev.preferred_magnitude() or (ev.magnitudes[0] if ev.magnitudes else None)
                da = hav(o.latitude, o.longitude, *stations[A])
                db = hav(o.latitude, o.longitude, *stations[B])
                ra = refp[A] - (o.time.timestamp + da / a.vp)
                rb = refp[B] - (o.time.timestamp + db / a.vp)
                cand = abs(ra) <= a.tol and abs(rb) <= a.tol
                found += 1
                rows.append({"alarm": al["alarm"], "catalogue": name,
                             "origin_utc": o.time.strftime("%Y-%m-%d %H:%M:%S.%f")[:-4],
                             "magnitude": f"{m.mag:.1f} {m.magnitude_type or ''}".strip() if m else "",
                             "dist_a_km": f"{da:.0f}", "dist_b_km": f"{db:.0f}",
                             "p_residual_a_s": f"{ra:.1f}", "p_residual_b_s": f"{rb:.1f}",
                             "candidate": "yes" if cand else "no", "agencies": " ".join(agencies)})
                print(f"alarm {al['alarm']:>2} {A}-{B} {name}: {rows[-1]['origin_utc']} "
                      f"M {rows[-1]['magnitude']}, {da:.0f}/{db:.0f} km, P residuals {ra:+.1f}/{rb:+.1f} s"
                      f"{'  <- candidate' if cand else ''}  [{' '.join(agencies)}]")
        if not found:
            print(f"alarm {al['alarm']:>2} {A}-{B}: no event in {', '.join(clients)}")
            rows.append({"alarm": al["alarm"], "catalogue": "+".join(clients), "origin_utc": "",
                         "magnitude": "", "dist_a_km": "", "dist_b_km": "", "p_residual_a_s": "",
                         "p_residual_b_s": "", "candidate": "none found", "agencies": ""})
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
