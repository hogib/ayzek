"""Shared helpers of the three exploratory held-out reruns of 2026-11-14
(tools/heldout_min3.py, tools/heldout_globalcat.py, tools/heldout_iasp91.py).

Exploratory, outside the held-out protocol (docs/protocol-heldout.md): nothing
here changes a frozen setting or a release output. Standard library only.
"""
import datetime as dt
import hashlib
import json
import math
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NETS = ("greece_now", "nz_now", "ko_now")
RELEASE = ROOT / "data/runs/heldout"


def stamp():
    return time.strftime("%FT%TZ", time.gmtime())


def say(msg):
    print(f"{stamp()} {msg}", flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tree_digest(d):
    """One SHA-256 over a directory: every file's relative path and content
    hash, in path order (symlinks followed). {"files": n, "sha256": ...}."""
    d = Path(d)
    h, n = hashlib.sha256(), 0
    for p in sorted(x for x in d.rglob("*") if x.is_file()):
        h.update(f"{p.relative_to(d).as_posix()} {sha256(p)}\n".encode())
        n += 1
    return {"files": n, "sha256": h.hexdigest()}


def code_commit():
    """The tree's commit: git's HEAD, or the COMMIT file of a `git archive` copy
    (the server mirror), with whether src, app or tools differ from it."""
    if (ROOT / ".git").exists():
        head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "src", "app", "tools"],
                               capture_output=True, text=True).stdout.strip()
        return {"commit": head, "dirty": bool(dirty), "from": "git"}
    f = ROOT / "COMMIT"
    return {"commit": f.read_text().strip() if f.exists() else "", "dirty": None, "from": "COMMIT file"}


def guard_out(out):
    """The output directory, refused when it is, or lies inside, a release
    output (data/runs/heldout and the caches the release runs read)."""
    out = Path(out).resolve()
    for keep in ("heldout", "heldout_scores", "heldout_quiet_scores", "heldout_r38", "heldout_f05436a"):
        k = (ROOT / "data/runs" / keep).resolve()
        if out == k or k in out.parents:
            sys.exit(f"{out}: inside {k}, which this script must not write to")
    return out


def write_provenance(path, script, params, inputs, extra=None):
    """The provenance record: argv, parameters, code commit, versions, input
    hashes. Written through a temporary name."""
    rec = {"script": str(Path(script).resolve()), "script_sha256": sha256(script),
           "common_sha256": sha256(__file__), "argv": sys.argv, "parameters": params,
           "ayzek_code": code_commit(),
           "versions": {"python": sys.version, "platform": platform.platform(), **(extra or {})},
           "inputs": inputs, "written_utc": stamp(),
           "note": "exploratory, outside the held-out protocol; no frozen setting or release output changed"}
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(rec, indent=1, default=str) + "\n")
    tmp.rename(path)


class Progress:
    """One line per finished unit: count, label, elapsed, ETA (elapsed per unit
    so far times the units left: rough when units differ in cost)."""

    def __init__(self, total, what="units"):
        self.total, self.done, self.t0, self.what = total, 0, time.time(), what

    def step(self, label, n=1):
        self.done += n
        el = time.time() - self.t0
        eta = el / self.done * max(0, self.total - self.done) if self.done else float("nan")
        say(f"{self.done}/{self.total} {self.what} | {label} | elapsed {el / 60:.1f} min, ETA {eta / 60:.1f} min")


# ---- geometry, recordings and reports (items 2 and 3) -------------------------

def distance_km(lat1, lon1, lat2, lon2):
    """As ayzek's distance_km (src/pipeline/network.cpp): haversine, R 6371 km."""
    r = math.pi / 180.0
    a = (math.sin((lat2 - lat1) * r / 2) ** 2
         + math.cos(lat1 * r) * math.cos(lat2 * r) * math.sin((lon2 - lon1) * r / 2) ** 2)
    return 2 * 6371.0 * math.asin(math.sqrt(a))


def km_to_deg(km):
    return km / (6371.0 * math.pi / 180.0)


def read_recording(path):
    """An ayzek recording (src/pipeline/recording.hpp): stations {code: (lat,
    lon)}, t_first, t_last, the run's --from (0 without), and the station
    detections per station in file order: (window_start, declared_at)."""
    stations, det, t_first, t_last, frm = {}, {}, None, None, 0.0
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] == "S":
                stations[p[1]] = (float(p[2]), float(p[3]))
            elif p[0] == "T":
                t_first, t_last = float(p[1]), float(p[2])
            elif p[0] == "D":
                det.setdefault(p[1], []).append((float(p[2]), float(p[3])))
            elif line.startswith("# from "):
                frm = float(line.split()[2])
    if t_first is None:
        sys.exit(f"{path}: no T line")
    return {"stations": stations, "t_first": t_first, "t_last": t_last, "from": frm, "detections": det}


HEAD = re.compile(r"^Event #(\d+): detected (?:M(-?[\d.]+) earthquake|earthquake \(no magnitude\))"
                  r"(?: at (-?[\d.]+)N (-?[\d.]+)E, origin (\S+) UTC)?", re.M)
ALARM_LINE = re.compile(r"^  alarm\s+(\S+) UTC by ([^;\n]+)(?:; later (.*))?$", re.M)
CREDIT = re.compile(r"^  AFAD\s+(?:(no matching catalogue event)|(\w+) (-?[\d.]+) at (\S+) UTC:)", re.M)


def epoch_near(hms, ref, lo=None, hi=None):
    """The epoch of a time of day HH:MM:SS.ss nearest `ref` (a replay is a few
    hours long, so the day is the one that puts it nearest)."""
    h, m, s = hms.split(":")
    sec = int(h) * 3600 + int(m) * 60 + float(s)
    day = math.floor(ref / 86400.0) * 86400.0
    return min((day + k * 86400.0 + sec for k in (-1, 0, 1)), key=lambda t: abs(t - ref))


def read_report(path, rec):
    """The alarms of an ayzek report, in alarm-number order: number, alarm epoch,
    the stations that raised it (the `by` list: those whose detection was
    declared by the alarm time, which is what Network::match uses), the later
    ones, the final location (lat, lon, origin epoch) or None, the final Pd
    magnitude of the headline or None, and the report's own credit: None or
    (origin HH:MM:SS.ss, type, magnitude as printed)."""
    text = Path(path).read_text()
    body = text[text.find("\nREPORT "):] if "\nREPORT " in text else text
    mid = 0.5 * (rec["t_first"] + rec["t_last"])
    alarms = []
    for block in re.split(r"(?m)^(?=Event #\d+: )", body)[1:]:
        block = block.split("\n\n")[0]
        h, a, c = HEAD.match(block), ALARM_LINE.search(block), CREDIT.search(block)
        if not (h and a and c):
            sys.exit(f"{path}: cannot read event block: {block[:120]!r}")
        t = epoch_near(a.group(1), mid)
        loc = None
        if h.group(3):
            loc = (float(h.group(3)), float(h.group(4)), epoch_near(h.group(5), t))
        alarms.append({"number": int(h.group(1)), "time": t, "hms": a.group(1),
                       "by": [s.strip() for s in a.group(2).split(",")],
                       "later": [s.strip() for s in (a.group(3) or "").split(",") if s.strip()],
                       "location": loc, "magnitude": float(h.group(2)) if h.group(2) else None,
                       "credit": None if c.group(1) else (c.group(4), c.group(2), c.group(3))})
    m = re.search(r"Alarms raised\s+(\d+)", text)
    if not m or int(m.group(1)) != len(alarms):
        sys.exit(f"{path}: read {len(alarms)} alarms, the report counts {m.group(1) if m else 'none'}")
    return alarms


def alarm_detections(alarm, rec, vp=6.0, slack=3.0, used=None):
    """The station detections an alarm was raised by: {station: (window_start,
    declared_at)}, found in the recording (the report does not print them).
    Rule: the alarm time is the declaration time of the last of them, so the
    station(s) of the `by` list with a detection declared at that time (within
    the 0.01 s the report prints) anchor it; every other `by` station gets its earliest
    detection declared by the alarm time that is compatible with every one
    already taken (window starts within the inter-station distance / 6 km/s
    + 3 s, Network::compatible). A station without one is left out and named
    in the returned list of problems. This reconstructs the association and
    is checked by the caller where it matters (heldout_iasp91.py reproduces the
    report's own credits with it before changing anything)."""
    T, out, problems = alarm["time"], {}, []
    D, S = rec["detections"], rec["stations"]
    for st in alarm["by"]:
        hit = [d for d in D.get(st, []) if -0.001 <= d[1] - T <= 0.011]   # the report cuts times to 0.01 s
        if hit:
            out[st] = hit[0]
    if not out:
        problems.append("no detection declared at the alarm time")

    def limit(a, b):
        return distance_km(*S[a], *S[b]) / vp + slack if a in S and b in S else 30.0

    def fits(st, d, chosen):
        return all(abs(d[0] - od[0]) <= limit(st, o) + 1e-9 for o, od in chosen.items())

    # The other stations: one detection each, pairwise compatible with every one
    # taken (depth-first over the candidates). Candidates not already taken by
    # an earlier alarm come first (`used`, kept by the caller over the alarms in
    # number order); then the earliest, or, when no detection anchors the alarm
    # (the transformer's merge replaced it), the latest.
    anchored = bool(out)
    rest = [st for st in alarm["by"] if st not in out]
    cands = {}
    for st in rest:
        c = [d for d in D.get(st, []) if d[1] <= T + 0.011 and d[0] > T - 300 and fits(st, d, out)]
        c.sort(key=lambda d: ((st, d[0]) in used if used is not None else False, d[0] if anchored else -d[0]))
        cands[st] = c

    def search(k, chosen):
        if k == len(rest):
            return chosen
        st = rest[k]
        for d in cands[st]:
            if fits(st, d, chosen):
                got = search(k + 1, {**chosen, st: d})
                if got is not None:
                    return got
        return None

    full = search(0, dict(out))
    if full is not None:
        out = full
    else:                       # no complete compatible choice: station by station, the gaps named
        for st in rest:
            c = [d for d in cands[st] if fits(st, d, out)]
            if c:
                out[st] = c[0]
            else:
                problems.append(f"{st}: no compatible detection")
    if used is not None:
        used.update((st, d[0]) for st, d in out.items())
        # The alarm's `later` stations hold a detection too (it joined after the
        # alarm), which a later alarm must not be given: one per station is
        # reserved in `used`, in declaration order, declared after the alarm
        # time, compatible with everything taken for this alarm. Not returned:
        # Network::match does not use them. (Without this, greece_now_seq020
        # transformer alarm #11 was given HORT's detection of alarm #10.)
        held = dict(out)
        later = sorted((d[1], st, d) for st in alarm["later"] for d in D.get(st, [])
                       if d[1] > T + 0.011 and d[0] < T + 300 and (st, d[0]) not in used)
        for _, st, d in later:
            if st not in held and fits(st, d, held):
                held[st] = d
                used.add((st, d[0]))
    return out, problems


def read_catalogue(path):
    """The set's catalogue as ayzek's reader takes it (load_afad_catalog): rows
    with a readable time, latitude, longitude, depth and magnitude, in time
    order: dicts with time (epoch), lat, lon, depth, mag, type, depth_text."""
    out = []
    with open(path, encoding="utf-8-sig") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split(",")
            if len(p) < 7:
                continue
            try:
                t = dt.datetime.strptime(p[0], "%d/%m/%Y %H:%M:%S").replace(tzinfo=dt.UTC).timestamp()
                ev = {"time": t, "lon": float(p[1]), "lat": float(p[2]), "depth": float(p[3]),
                      "mag": float(p[6]), "type": p[5], "depth_text": p[3].strip(),
                      "id": p[8].strip() if len(p) > 8 else ""}
            except ValueError:
                continue
            out.append(ev)
    out.sort(key=lambda e: e["time"])
    return out


def hms(t):
    """HH:MM:SS.ss of an epoch, as ayzek prints times."""
    c = int(t * 1000.0 + 0.5) // 10     # milliseconds, cut to 0.01 s (common.hpp hms)
    return f"{c // 360000 % 24:02d}:{c // 6000 % 60:02d}:{c % 6000 // 100:02d}.{c % 100:02d}"
