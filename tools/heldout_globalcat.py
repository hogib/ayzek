"""Global-catalogue check of the held-out quiet windows' alarms: could a
teleseism or regional earthquake of a global catalogue explain an alarm of the
learned detector or of STA/LTA (stalta_eqfa) on the 62 quiet windows?

    uv run --project tools python -u tools/heldout_globalcat.py [--out data/runs/heldout_globalcat]

Exploratory, outside the held-out protocol (docs/protocol-heldout.md): the
quiet windows were drawn with each network's own catalogue; this looks at other
catalogues after the fact. No detector is rerun and nothing under
data/runs/heldout is written. Needs the internet (FDSN event services), so it
runs on this machine, not on the server; the requests carry a time span, a
magnitude, and for the regional one the stations' centroid rounded to 0.1
degree and a radius, with a generic User-Agent and nothing else.

Inputs: RELEASE/NET/quiet/quiet_index.csv (the windows), and per detector
RELEASE/NET/quiet/LABEL/WINDOW.txt (the alarms: time and the stations that
raised each) and WINDOW.rec (the stations' coordinates and detections; an
alarm's detections are found by heldout_explore_common.alarm_detections).

Catalogues, per window and service (--services; default usgs,isc; emsc known):
  global    M >= 4.0 anywhere, origin from window start - 7200 s to window end + 1800 s;
  regional  M >= 2.5 within 10 degrees of the stations' centroid, same span.
Every response is kept as received (OUT/catalogues/SERVICE/WINDOW_TIER.txt,
with its URL, status and SHA-256 in OUT/queries.json) and is not asked for
again on a rerun. Catalogues are not merged: every figure is per service.

The rule, fixed before any catalogue was looked at (--tiers, --tolerance,
--min-explained change it, and the provenance records them):
  1. Candidates at a station: events of M >= 4.5 at any distance, M >= 4.0
     within 30 degrees, M >= 2.5 within 10 degrees of that station.
  2. A detection's time is its window start + 3.5 s, the P date ayzek's matcher
     uses (Network::match), for either detector.
  3. Predicted arrivals: obspy TauP, iasp91, the first P-type arrival (p, P,
     Pn, Pg, Pdiff, PKP, PKIKP) and the first S-type arrival (s, S, Sn, Sg,
     Sdiff, SKS), from the catalogue's origin, epicentre and depth (10 km where
     the catalogue gives none, 700 km at most).
  4. A station is explained by an event and phase type when its detection time
     is within +-6 s of that predicted arrival (the matcher's tolerance).
  5. An alarm is explained by an event when at least 2 of the stations that
     raised it are explained by that event, with the same phase type; where
     several events do, the one with the smallest mean absolute residual is
     listed.
  6. Chance: the same rule with every detection time moved by +-600, +-1200 and
     +-1800 s; the explained alarms under each shift are counted beside the
     observed count (later phases, surface waves and events absent from the
     catalogue are not covered by the rule; the nearest miss of each alarm is
     listed so that they can be looked for).

Output (OUT, each file through a temporary name; OUT/DONE at the end):
  alarms.csv   one row per alarm and service: explained or not, the event, the
               residuals, and the nearest miss;
  summary.csv  per network, detector and service: alarms, explained, explained
               under each time shift;
  events.csv   the catalogue events that explain an alarm;
  queries.json, provenance.json (argv, parameters, commit, versions, hashes).
One progress line per finished unit (a catalogue request, then a window of a
detector) with elapsed and ETA; exit status nonzero when no alarm is checked
or a service returned nothing usable for every window (--allow-missing lets
a window's failed request pass, recorded).
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from heldout_explore_common import (NETS, ROOT, Progress, alarm_detections, distance_km, guard_out, hms,  # noqa: E402
                                    km_to_deg, read_recording, read_report, say, sha256, tree_digest,
                                    write_provenance)

SERVICES = {"usgs": "https://earthquake.usgs.gov/fdsnws/event/1/query",
            "isc": "https://www.isc.ac.uk/fdsnws/event/1/query",
            "emsc": "https://www.seismicportal.eu/fdsnws/event/1/query"}
USER_AGENT = "fdsn-event-client/1.0"
BEFORE_S, AFTER_S = 7200, 1800
GLOBAL_MIN_MAG, REGIONAL_MIN_MAG, REGIONAL_RADIUS_DEG = 4.0, 2.5, 10.0
P_PHASES = ["p", "P", "Pn", "Pg", "Pdiff", "PKP", "PKIKP"]
S_PHASES = ["s", "S", "Sn", "Sg", "Sdiff", "SKS"]
SHIFTS = [0, 600, -600, 1200, -1200, 1800, -1800]
MAX_TT_S = 2400.0          # no first P or S arrival is later than this after the origin


def iso(t):
    return dt.datetime.fromtimestamp(t, dt.UTC).strftime("%Y-%m-%dT%H:%M:%S")


def fetch(url, tries=3):
    """(status, body) of a GET; 204 and FDSN's 404-for-no-data are an empty catalogue."""
    last = None
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": USER_AGENT}), timeout=90) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            if e.code in (204, 404):
                return e.code, b""
            last = f"HTTP {e.code}"
        except Exception as e:  # noqa: BLE001  (timeouts, DNS, resets: retried, then recorded)
            last = repr(e)
        time.sleep(5 * (k + 1))
    return None, (last or "failed").encode()


def parse_fdsn_text(body):
    """Events of an fdsnws-event text listing: id, time (epoch), lat, lon, depth
    (km or None), magnitude, magnitude type, place. Rows without a time,
    position or magnitude are left out."""
    out, head = [], None
    for line in body.decode("utf-8", "replace").splitlines():
        if line.startswith("#"):
            if head is None and "|" in line:
                head = [c.strip().lower() for c in line.lstrip("#").split("|")]
            continue
        if not line.strip() or head is None:
            continue
        f = dict(zip(head, (c.strip() for c in line.split("|"))))
        try:
            t = dt.datetime.fromisoformat(f["time"].replace("Z", "")).replace(tzinfo=dt.UTC).timestamp()
            ev = {"id": f.get("eventid", ""), "time": t, "lat": float(f["latitude"]), "lon": float(f["longitude"]),
                  "mag": float(f["magnitude"]), "magtype": f.get("magtype", ""),
                  "place": f.get("eventlocationname", "")}
        except (KeyError, ValueError):
            continue
        try:
            ev["depth"] = float(f.get("depth/km", ""))
        except ValueError:
            ev["depth"] = None
        out.append(ev)
    return out


class Arrivals:
    def __init__(self):
        from obspy.taup import TauPyModel
        self.model, self.cache = TauPyModel("iasp91"), {}

    def __call__(self, ev, lat, lon):
        """(first P-type, first S-type) arrival epochs at a point, None where no such phase arrives."""
        key = (ev["id"], ev["time"], round(lat, 4), round(lon, 4))
        if key not in self.cache:
            deg = km_to_deg(distance_km(ev["lat"], ev["lon"], lat, lon))
            depth = 10.0 if ev["depth"] is None else min(max(ev["depth"], 0.0), 700.0)
            out = []
            for phases in (P_PHASES, S_PHASES):
                arr = self.model.get_travel_times(source_depth_in_km=depth, distance_in_degree=deg, phase_list=phases)
                out.append(ev["time"] + min(x.time for x in arr) if arr else None)
            self.cache[key] = tuple(out)
        return self.cache[key]


def explain(dets, stations, events, arrivals, tiers, tol, need, shift):
    """Rule 1-5 for one alarm, its detection times moved by `shift`. Returns
    (best explaining event or None, nearest miss): each (event, phase type,
    {station: residual}, mean |residual|)."""
    best, near = None, None
    obs = {st: ws + 3.5 + shift for st, (ws, _) in dets.items() if st in stations}
    if not obs:
        return None, None
    lo, hi = min(obs.values()) - MAX_TT_S, max(obs.values()) + tol
    for ev in events:
        if not lo <= ev["time"] <= hi:
            continue
        res = {"P": {}, "S": {}}
        for st, t in obs.items():
            deg = km_to_deg(distance_km(ev["lat"], ev["lon"], *stations[st]))
            if not any(ev["mag"] >= m and deg <= r for m, r in tiers):
                continue
            for kind, pred in zip("PS", arrivals(ev, *stations[st])):
                if pred is not None:
                    res[kind][st] = t - pred
        for kind, r in res.items():
            if not r:
                continue
            ok = {st: x for st, x in r.items() if abs(x) <= tol}
            if len(ok) >= need:
                cand = (ev, kind, ok, sum(abs(x) for x in ok.values()) / len(ok))
                if best is None or cand[3] < best[3]:
                    best = cand
            if len(r) >= need:      # nearest miss: the need-th smallest |residual| over the stations
                worst = sorted(abs(x) for x in r.values())[need - 1]
                if near is None or worst < near[3]:
                    near = (ev, kind, r, worst)
    return best, near


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(ROOT / "data/runs/heldout_globalcat"))
    ap.add_argument("--nets", default=",".join(NETS))
    ap.add_argument("--detectors", default="learned,stalta_eqfa")
    ap.add_argument("--services", default="usgs,isc")
    ap.add_argument("--tiers", default="4.5:180,4.0:30,2.5:10", help="MINMAG:MAXDEG,... (rule 1)")
    ap.add_argument("--tolerance", type=float, default=6.0, help="seconds (rule 4)")
    ap.add_argument("--min-explained", type=int, default=2, help="stations (rule 5)")
    ap.add_argument("--pause", type=float, default=1.0, help="seconds between requests")
    ap.add_argument("--windows", help="only these windows (comma-separated; for a quick look, not marked done)")
    ap.add_argument("--allow-missing", action="store_true", help="go on when a request fails after its retries")
    ap.add_argument("--release", type=Path, default=ROOT / "data/runs/heldout")
    a = ap.parse_args()
    out = guard_out(a.out)
    (out / "catalogues").mkdir(parents=True, exist_ok=True)
    nets, detectors, services = a.nets.split(","), a.detectors.split(","), a.services.split(",")
    tiers = [tuple(float(x) for x in t.split(":")) for t in a.tiers.split(",")]
    import obspy

    # The windows, their stations (from the first detector's recording) and alarms.
    windows = []
    for net in nets:
        for w in csv.DictReader(open(a.release / net / "quiet" / "quiet_index.csv")):
            if a.windows and w["name"] not in a.windows.split(","):
                continue
            t0 = dt.datetime.fromisoformat(w["start_utc"]).timestamp()
            t1 = dt.datetime.fromisoformat(w["end_utc"]).timestamp()
            rec = read_recording(a.release / net / "quiet" / detectors[0] / f"{w['name']}.rec")
            lat = sum(v[0] for v in rec["stations"].values()) / len(rec["stations"])
            lon = sum(v[1] for v in rec["stations"].values()) / len(rec["stations"])
            windows.append({"net": net, "name": w["name"], "t0": t0, "t1": t1, "lat": lat, "lon": lon})
    if not windows:
        sys.exit("no quiet windows")
    params = {"services": {s: SERVICES[s] for s in services}, "user_agent": USER_AGENT, "tiers": tiers,
              "tolerance_s": a.tolerance, "min_explained": a.min_explained, "shifts": SHIFTS,
              "span": {"before_s": BEFORE_S, "after_s": AFTER_S}, "global_min_mag": GLOBAL_MIN_MAG,
              "regional": {"min_mag": REGIONAL_MIN_MAG, "radius_deg": REGIONAL_RADIUS_DEG},
              "p_phases": P_PHASES, "s_phases": S_PHASES, "taup_model": "iasp91", "detectors": detectors,
              "nets": nets, "windows": a.windows, "release": str(a.release)}
    inputs = {"quiet_index": {n: sha256(a.release / n / "quiet" / "quiet_index.csv") for n in nets},
              "quiet_runs": {f"{n}/{d}": tree_digest(a.release / n / "quiet" / d) for n in nets for d in detectors},
              "release_manifests": {n: sha256(a.release / n / "manifest.json") for n in nets
                                    if (a.release / n / "manifest.json").exists()}}
    versions = {"obspy": obspy.__version__}
    write_provenance(out / "provenance.json", __file__, params, inputs, versions)
    prog = Progress(len(windows) * len(services) * 2 + len(windows) * len(detectors),
                    "units (catalogue requests, then windows x detectors)")
    say(f"{len(windows)} windows, detectors {detectors}, services {services} -> {out}")

    # ---- catalogues ----------------------------------------------------------
    qfile = out / "queries.json"
    queries = json.loads(qfile.read_text()) if qfile.exists() else {}
    catalogue = defaultdict(dict)      # (service, window) -> {event id: event}
    failed = defaultdict(int)
    for w in windows:
        span = {"starttime": iso(w["t0"] - BEFORE_S), "endtime": iso(w["t1"] + AFTER_S), "format": "text"}
        tiersq = {"global": {**span, "minmagnitude": GLOBAL_MIN_MAG},
                  "regional": {**span, "minmagnitude": REGIONAL_MIN_MAG, "latitude": f"{w['lat']:.1f}",
                               "longitude": f"{w['lon']:.1f}", "maxradius": REGIONAL_RADIUS_DEG}}
        for svc in services:
            for tier, q in tiersq.items():
                key = f"{svc}/{w['name']}_{tier}"
                f = out / "catalogues" / svc / f"{w['name']}_{tier}.txt"
                f.parent.mkdir(exist_ok=True)
                url = SERVICES[svc] + "?" + urllib.parse.urlencode(q)
                if not (f.exists() and queries.get(key, {}).get("url") == url and queries[key].get("ok")):
                    status, body = fetch(url)
                    ok = status is not None
                    if ok:
                        tmp = f.with_suffix(".txt.tmp")
                        tmp.write_bytes(body)
                        tmp.rename(f)
                    queries[key] = {"url": url, "status": status, "ok": ok, "asked_utc": iso(time.time()),
                                    "bytes": len(body) if ok else 0,
                                    "sha256": hashlib.sha256(body).hexdigest() if ok else None,
                                    "error": None if ok else body.decode("utf-8", "replace")}
                    tmpq = qfile.with_suffix(".json.tmp")
                    tmpq.write_text(json.dumps(queries, indent=1) + "\n")
                    tmpq.rename(qfile)
                    time.sleep(a.pause)
                if queries[key]["ok"]:
                    evs = parse_fdsn_text(f.read_bytes())
                    queries[key]["events"] = len(evs)
                    for e in evs:
                        catalogue[(svc, w["name"])].setdefault(e["id"] or f"{e['time']}", e)
                else:
                    failed[svc] += 1
                    if not a.allow_missing:
                        sys.exit(f"{key}: request failed ({queries[key]['error']}); rerun to retry (answers so far are "
                                 "kept), or --allow-missing")
                prog.step(f"{key}: {queries[key].get('events', 'FAILED')} events")
    qfile.with_suffix(".json.tmp").write_text(json.dumps(queries, indent=1) + "\n")
    qfile.with_suffix(".json.tmp").rename(qfile)
    for svc in services:
        if failed[svc] == 2 * len(windows):
            sys.exit(f"{svc}: every request failed")

    # ---- the rule --------------------------------------------------------------
    arrivals = Arrivals()
    rows, used_events = [], {}
    summary = defaultdict(lambda: {"alarms": 0, "no_detections": 0, **{f"explained_shift_{s:+d}": 0 for s in SHIFTS}})
    for w in windows:
        for det in detectors:
            rec = read_recording(a.release / w["net"] / "quiet" / det / f"{w['name']}.rec")
            alarms = sorted(read_report(a.release / w["net"] / "quiet" / det / f"{w['name']}.txt", rec),
                            key=lambda x: x["number"])
            used = set()    # detections are reserved across a report's alarms, as in heldout_iasp91.py;
            # checked 2026-11-16: with and without it, alarms.csv of the 62 quiet windows is byte-identical
            for al in alarms:
                dets, problems = alarm_detections(al, rec, used=used)
                for svc in services:
                    events = sorted(catalogue[(svc, w["name"])].values(), key=lambda e: e["time"])
                    g = summary[(w["net"], det, svc)]
                    g["alarms"] += 1
                    g["no_detections"] += len(dets) < a.min_explained
                    hit = near = None
                    for shift in SHIFTS:
                        b, n = explain(dets, rec["stations"], events, arrivals, tiers, a.tolerance, a.min_explained, shift)
                        g[f"explained_shift_{shift:+d}"] += b is not None
                        if shift == 0:
                            hit, near = b, n
                    def cols(x):
                        if x is None:
                            return ["", "", "", "", "", "", "", "", ""]
                        ev, kind, res, score = x
                        st0 = next(iter(res))
                        return [ev["id"], iso(ev["time"]), f"{ev['mag']:.1f}", ev["magtype"],
                                "" if ev["depth"] is None else f"{ev['depth']:.0f}",
                                f"{km_to_deg(distance_km(ev['lat'], ev['lon'], *rec['stations'][st0])):.1f}", kind,
                                " ".join(f"{s}:{r:+.1f}" for s, r in res.items()), ev["place"]]
                    rows.append([w["net"], det, w["name"], al["number"], iso(al["time"]) + "." + hms(al["time"])[-2:],
                                 "+".join(al["by"]), len(dets), "; ".join(problems), svc,
                                 int(hit is not None), *cols(hit), *cols(near)])
                    if hit is not None:
                        used_events[(svc, hit[0]["id"])] = hit[0]
            prog.step(f"{w['name']} {det}: {len(alarms)} alarms")
    if not rows:
        sys.exit("no alarm checked")

    def write(name, head, body):
        tmp = out / (name + ".tmp")
        with open(tmp, "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(head)
            wr.writerows(body)
        tmp.rename(out / name)

    ev_cols = ["event_id", "origin_utc", "magnitude", "magnitude_type", "depth_km", "distance_deg", "phase",
               "residuals_s", "place"]
    write("alarms.csv", ["network", "detector", "window", "alarm", "alarm_utc", "raised_by", "detections_found",
                         "association_problems", "service", "explained", *ev_cols,
                         *(f"nearest_{c}" for c in ev_cols)], rows)
    write("summary.csv", ["network", "detector", "service", "alarms", "alarms_with_too_few_detections",
                          *(f"explained_shift_{s:+d}s" for s in SHIFTS)],
          [[n, d, s, g["alarms"], g["no_detections"], *(g[f"explained_shift_{x:+d}"] for x in SHIFTS)]
           for (n, d, s), g in sorted(summary.items())])
    write("events.csv", ["service", "event_id", "origin_utc", "latitude", "longitude", "depth_km", "magnitude",
                         "magnitude_type", "place"],
          [[s, e["id"], iso(e["time"]), e["lat"], e["lon"], e["depth"], e["mag"], e["magtype"], e["place"]]
           for (s, _), e in sorted(used_events.items(), key=lambda x: (x[0][0], x[1]["time"]))])
    write_provenance(out / "provenance.json", __file__,
                     {**params, "finished": True, "failed_requests": dict(failed)},
                     {**inputs, "queries": {k: {"url": v["url"], "sha256": v["sha256"], "status": v["status"]}
                                            for k, v in queries.items()}}, versions)
    for (n, d, s), g in sorted(summary.items()):
        say(f"{n} {d} {s}: {g['explained_shift_+0']} of {g['alarms']} alarms explained; under shifts "
            + ", ".join(str(g[f"explained_shift_{x:+d}"]) for x in SHIFTS[1:]))
    if a.windows:
        say("partial run (--windows): not marked done")
    else:
        (out / "DONE").write_text("complete\n")
    say("DONE")


if __name__ == "__main__":
    main()
