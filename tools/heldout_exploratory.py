"""Exploratory figures from a finished held-out test, outside the frozen protocol.

    python3 tools/heldout_exploratory.py --run data/runs/heldout/greece_now \\
        --selection docs/heldout-selections/greece_now --quiet-sets 13,17,4 \\
        --json exploratory_greece.json

Nothing here is an endpoint of docs/protocol-heldout.md: every figure was
thought of after the tests were run and read, and is exploratory. The tool
reads the outputs of tools/run_heldout.py (RUN/outcomes, RUN/runs, RUN/quiet,
RUN/analysis/magnitude_latency.csv) and the selection (SELECTION/selection.json:
station coordinates, each selected sequence's `first` event, and each quiet
window's `sequence`, the sequence whose station set the window was drawn on);
it runs nothing and writes only the report to stdout and, with --json, the
same figures to a file (with --mainshocks, item 9's table to a CSV file).

--detectors names the labels (default learned, stalta_eqfa,
phasenet_causal_eqfa); the first two are the pair of the paired figures.

1. Paired latency: over the events both detectors of the pair detected (joined
   on set, origin, type and magnitude in magnitude_latency.csv), the median and
   interquartile range (linear interpolation between order statistics) of the
   first's alarm delay minus the second's, and how often each was first.
2. Per sequence and detector: scored events, detected, the chance level c (the
   mean over the shifted catalogues of the share of shifted events detected),
   the chance-corrected count n (o - c) / (1 - c) (tools/significance.py's
   corrected(), called here), and the unmatched alarms (the entries of the
   report's "False alarms" list); then the number of sequences on which the
   first detector detected more, the same and fewer events than the second.
3. Per detector over the network: detected, the sum of the sequences'
   corrected counts, and the effective chance level (raw - corrected) /
   (n - corrected), the single c that would give that sum.
4. Quiet windows by station set: the windows grouped by their `sequence`, with
   the scored hours (quiet/quiet_index.csv) and each detector's alarms (the
   ALARM lines of quiet/LABEL/WINDOW.txt) and rate per day; and the alarms of
   every window. --quiet-sets N,N,... adds one line for those sets together
   and one for the remaining windows.
5. Closest station pair: per set (sequence or quiet window) the two stations
   nearest each other (haversine on the selection's coordinates of the
   stations that ran: those of the report's "Station processing" table), and
   per detector how many matched, unmatched and quiet alarms were declared by
   exactly that pair (the stations on the ALARM line). An alarm of a sequence
   set is unmatched when the report's "False alarms" list has it and matched
   otherwise: the "(AFAD ...)" note on an ALARM line is the match at the time
   of the alarm, which the report's final matching can withdraw.
6. Detection by catalogue-magnitude band (<2, 2-3, 3-4, >=4; lower edge
   included) for the pair, on the magnitudes of the outcome tables. Those are
   the report's, printed with one decimal, so an event whose catalogue
   magnitude is within 0.05 of an edge can fall in the neighbouring band.
7. Pd estimates of the unmatched alarms (sequence sets) and of the quiet
   windows' alarms: alarms, those with an estimate (the event block's
   "Pd magn." line), those whose final estimate is at least --large (4.0),
   each of the latter listed.
8. With --other-events FILE: events of another agency's bulletin that the
   network's catalogue lacks, and whether an unmatched alarm follows them.
   FILE is an FDSN event service text listing (fdsnws-event, format=text;
   the columns `Time` and `Author` are read); --other-author keeps the
   entries of one author, --other-set N names the sequence set whose report
   is read. Per detector: the events followed within --follow-seconds (60) by
   an unmatched alarm of that set, and the distinct alarms that do so (one
   alarm can follow two events that lie close together). Whether a listed
   event is missing from the network's catalogue is for the caller to
   establish; the tool counts what it is given.
9. The sequence-initiating events: per selected sequence the selection's
   `first` event (origin, magnitude and type, depth, epicentre, all the
   selection's) with the epicentral distance of the nearest and the second
   nearest station (haversine on the selection's coordinates, over the
   stations that ran: those of the first detector's "Station processing"
   table), and per detector:
   - the scored event: the event of the outcome table (shift 0) of the same
     type whose one-decimal magnitude is within 0.051 of the selection's and
     whose origin is within 1 s of it (the set's catalogue gives whole
     seconds); the nearest in time if several;
   - status, the outcome table's: `detected` or `missed` (its `detected`),
     `not scored` when the table has no such event; and its `detected_pre_s`;
   - from the report's event block whose "AFAD" line names the scored event
     (the earliest alarm if several): alarm time, the delay after origin as
     printed (against the catalogue's whole-second origin) and against the
     selection's origin, the declaring stations (the "alarm" line up to
     "; later"), the "pre-S" line's "n of m", the Pd estimate (first, seconds
     after the alarm, 0 when "at the alarm", and final), the rows of the
     warning table, those "before S", and the warning at the nearest station
     before S (at equal printed distance the smaller warning, as
     tools/run_heldout.py's warning_summary) with its S arrival picked or
     predicted;
   - note, by rule: the scored events earlier in the set (how many, and the
     last one); the further events the selection joined to the sequence;
     every other alarm of the report from --alarm-window B,A (30,60) seconds
     before to after the selection's origin, with the event its block's
     "AFAD" line credits it to (or none) and its final Pd estimate; any
     disagreement between the outcome table and the report.
   Then per detector over the network: scored, detected, detected before S
   (outcome table), with a Pd estimate at the alarm (and at any time), with a
   station before S at the alarm, and the medians of the printed alarm delay
   and of the warning at the nearest station before S over the detected
   events that have one. A missed event has no block, so no delay, estimate
   or warning: an alarm credited to another event is in the note only.

Standard library only.
"""

import argparse
import collections
import csv
import json
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from significance import corrected  # noqa: E402

BANDS = [
    ("<2", -math.inf, 2.0),
    ("2-3", 2.0, 3.0),
    ("3-4", 3.0, 4.0),
    (">=4", 4.0, math.inf),
]
ALARM = re.compile(
    r"^ALARM\s+#(\d+)\s+(\S+)\s+earthquake detected by ([^(\n]+?)\s*(?:\(AFAD.*\))?$",
    re.M,
)
# The final Pd estimate of an event block, as tools/run_heldout.py's EVENT reads it.
PD_FINAL = re.compile(
    r"^  Pd magn\.\s+M-?[\d.]+ (?:at the alarm|[\d.]+ s after the alarm), M(-?[\d.]+) ",
    re.M,
)
# The lines of an event block, as tools/run_heldout.py's EVENT and WARNING read them.
BLOCK_ALARM = re.compile(r"^  alarm\s+(\S+) UTC by ([^;\n]+)", re.M)
BLOCK_PD = re.compile(
    r"^  Pd magn\.\s+M(-?[\d.]+) (?:at the alarm|([\d.]+) s after the alarm), M(-?[\d.]+) ",
    re.M,
)
BLOCK_AFAD = re.compile(
    r"^  AFAD\s+(\w+) (-?[\d.]+) at (\S+) UTC: alarm ([-\d.]+) s after origin", re.M
)
BLOCK_PRE_S = re.compile(r"^  pre-S\s+(\d+) of (\d+) alarm detections", re.M)
WARNING = re.compile(
    r"^\s+(\S+)\s+(\d+) km\s+\S+ UTC \((picked|predicted)\)\s+([+-][\d.]+) s\s+(before|after) S",
    re.M,
)
MAINSHOCK_FIELDS = [
    "set",
    "origin",
    "magnitude",
    "type",
    "depth_km",
    "lat",
    "lon",
    "nearest_station",
    "nearest_km",
    "second_station",
    "second_km",
    "scored_as",
    "detector",
    "status",
    "detected_pre_s",
    "alarm",
    "alarm_time",
    "alarm_delay_s",
    "alarm_after_selection_origin_s",
    "declaring",
    "pre_s",
    "pd_first",
    "pd_first_after_alarm_s",
    "pd_final",
    "warning_stations",
    "stations_ahead",
    "nearest_ahead_station",
    "nearest_ahead_km",
    "nearest_ahead_warning_s",
    "nearest_ahead_s_source",
    "note",
]
UNMATCHED = re.compile(r"\s+\S+ UTC\s+#(\d+)\s")
STATION_ROW = re.compile(r"  (\S+)\s+\d+\s+\d+\s+\d+")


def quantile(values, q):
    """The q-quantile by linear interpolation between order statistics (the
    position q (n - 1) in the sorted values)."""
    v = sorted(values)
    pos = q * (len(v) - 1)
    lo = math.floor(pos)
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (pos - lo)


def paired_latency(rows, a, b):
    """From the rows of magnitude_latency.csv: the events both a and b detected,
    the median and quartiles of a's alarm delay minus b's, and who was first."""
    delay = {
        d: {
            (r["set"], r["origin"], r["type"], r["magnitude"]): float(
                r["alarm_delay_s"]
            )
            for r in rows
            if r["detector"] == d
        }
        for d in (a, b)
    }
    diffs = [delay[a][e] - delay[b][e] for e in sorted(delay[a]) if e in delay[b]]
    out = {
        "events": len(diffs),
        f"{a}_first": sum(d < 0 for d in diffs),
        f"{b}_first": sum(d > 0 for d in diffs),
        "ties": sum(d == 0 for d in diffs),
    }
    if diffs:
        out.update(
            median=quantile(diffs, 0.5),
            q1=quantile(diffs, 0.25),
            q3=quantile(diffs, 0.75),
        )
    return out


def load_outcomes(d, labels, sets, column="detected"):
    """{detector: {shift: {event: detected}}} from OUTCOMES/LABEL_SET.csv, the
    event (dataset, origin, type, magnitude) as tools/significance.py keys it.
    `column` names another 0/1 column (detected_pre_s); an empty cell is None."""
    out = {label: collections.defaultdict(dict) for label in labels}
    for label in labels:
        for s in sets:
            f = Path(d) / f"{label}_{s}.csv"
            if not f.exists():
                sys.exit(f"{f}: no outcome table")
            for r in csv.DictReader(open(f)):
                ev = (r["dataset"], r["origin"], r["type"], r["magnitude"])
                out[label][int(r["shift"])][ev] = (
                    int(r[column]) if r[column] != "" else None
                )
    return out


def sequence_row(det, s):
    """One detector's figures on sequence s: scored, detected, c, corrected."""
    obs = {e: v for e, v in det[0].items() if e[0] == s}
    shifted = {
        k: {e: v for e, v in ev.items() if e[0] == s} for k, ev in det.items() if k != 0
    }
    cs = [sum(ev.values()) / len(ev) for ev in shifted.values() if ev]
    corr, raw = corrected(
        obs, shifted, sorted(obs), {k: sorted(ev) for k, ev in shifted.items()}
    )
    return {
        "scored": len(obs),
        "detected": raw,
        "chance": sum(cs) / len(cs) if cs else 0.0,
        "corrected": corr,
    }


def more_same_fewer(table, a, b):
    """Over the sequences of `table` ({set: {detector: row}}): on how many a
    detected more events than b, the same number, fewer."""
    d = [row[a]["detected"] - row[b]["detected"] for row in table.values()]
    return {
        "more": sum(x > 0 for x in d),
        "same": sum(x == 0 for x in d),
        "fewer": sum(x < 0 for x in d),
    }


def network_chance(table, label):
    """Detected, corrected (the sum over sequences) and the effective chance
    level (raw - corrected) / (n - corrected) of one detector."""
    n = sum(row[label]["scored"] for row in table.values())
    raw = sum(row[label]["detected"] for row in table.values())
    corr = sum(row[label]["corrected"] for row in table.values())
    return {
        "scored": n,
        "detected": raw,
        "corrected": corr,
        "effective_chance": (raw - corr) / (n - corr) if n > corr else float("nan"),
    }


def haversine_km(lat1, lon1, lat2, lon2):
    p = math.radians
    h = (
        math.sin(p(lat2 - lat1) / 2) ** 2
        + math.cos(p(lat1)) * math.cos(p(lat2)) * math.sin(p(lon2 - lon1) / 2) ** 2
    )
    return 6371.0 * 2 * math.asin(math.sqrt(h))


def closest_pair(coords):
    """((station, station), km) of the two stations of {station: (lat, lon)}
    nearest each other, the pair in alphabetical order; None with fewer than two."""
    names = sorted(coords)
    pairs = [
        (haversine_km(*coords[x], *coords[y]), (x, y))
        for i, x in enumerate(names)
        for y in names[i + 1 :]
    ]
    if not pairs:
        return None
    km, pair = min(pairs)
    return pair, km


def parse_report(text):
    """An ayzek report's alarms {number: {time, stations, pd_final}},
    the numbers of its "False alarms" list, and the stations that ran."""
    alarms = {
        int(m.group(1)): {
            "time": m.group(2),
            "pd_final": None,
            "stations": [x.strip() for x in m.group(3).split(",")],
        }
        for m in ALARM.finditer(text)
    }
    for block in re.split(r"^Event #", text, flags=re.M)[1:]:
        number = int(re.match(r"\d+", block).group())
        m = PD_FINAL.search(block.split("\n\n")[0])
        if m and number in alarms:
            alarms[number]["pd_final"] = float(m.group(1))
    unmatched, stations = [], []
    lines = text.split("\n")
    for i, line in enumerate(lines):
        if line.startswith("  False alarms"):
            for row in lines[i + 1 :]:
                m = UNMATCHED.match(row)
                if not m:
                    break
                unmatched.append(int(m.group(1)))
        if line.startswith("Station processing"):
            for row in lines[i + 2 :]:
                m = STATION_ROW.match(row)
                if not m:
                    break
                stations.append(m.group(1))
    return {"alarms": alarms, "unmatched": unmatched, "stations": stations}


def magnitude_bands(det, labels):
    """Per band: scored events and, per detector, those detected (shift 0)."""
    out = []
    for name, lo, hi in BANDS:
        ev = [e for e in det[labels[0]][0] if lo <= float(e[3]) < hi]
        out.append(
            {
                "band": name,
                "scored": len(ev),
                **{label: sum(det[label][0][e] for e in ev) for label in labels},
            }
        )
    return out


def alarm_magnitudes(alarms, large):
    """Of (set, alarm) pairs: how many, with a Pd estimate, with a final
    estimate of at least `large`, and the latter listed."""
    with_pd = [(s, a) for s, a in alarms if a["pd_final"] is not None]
    big = [
        {"set": s, "time": a["time"], "pd_final": a["pd_final"]}
        for s, a in with_pd
        if a["pd_final"] >= large
    ]
    return {
        "alarms": len(alarms),
        "with_pd": len(with_pd),
        "large": len(big),
        "large_cases": big,
    }


def quiet_groups(windows, hours, counts, labels):
    """Windows grouped by station set (the selection's `sequence`), largest
    group first: windows, scored hours, alarms and rate per day per detector.
    `windows` is {name: sequence}."""
    groups = collections.defaultdict(list)
    for name, seq in windows.items():
        groups[seq].append(name)
    return [
        quiet_total(f"{seq:03d}", names, hours, counts, labels)
        for seq, names in sorted(groups.items(), key=lambda g: (-len(g[1]), g[0]))
    ]


def quiet_total(name, names, hours, counts, labels):
    h = sum(hours[n] for n in names)
    alarms = {label: sum(counts[label][n] for n in names) for label in labels}
    return {
        "set": name,
        "windows": len(names),
        "hours": h,
        "alarms": alarms,
        "per_day": {
            label: 24 * v / h if h else float("nan") for label, v in alarms.items()
        },
    }


def pair_shares(reports, pairs):
    """Per kind of alarm (quiet; on a sequence set unmatched, those of the
    report's "False alarms" list, or else matched): how many there were and
    how many were declared by exactly the set's closest pair. `reports` is
    [(set, kind of set, parsed report)], `pairs` {set: (pair, km) or None}."""
    out = {
        k: {"alarms": 0, "by_closest_pair": 0}
        for k in ("matched", "unmatched", "quiet")
    }
    for s, kind, rep in reports:
        pair = set(pairs[s][0]) if pairs[s] else None
        listed = set(rep["unmatched"])
        for number, a in rep["alarms"].items():
            k = (
                "quiet"
                if kind == "quiet"
                else "unmatched" if number in listed else "matched"
            )
            out[k]["alarms"] += 1
            out[k]["by_closest_pair"] += set(a["stations"]) == pair
    return out


def day_seconds(hms):
    """Seconds of the day of an HH:MM:SS[.ss] string."""
    h, m, s = hms.split(":")
    return int(h) * 3600 + int(m) * 60 + float(s)


def read_other_events(path, author=None):
    """The times (HH:MM:SS.sss) of an FDSN text listing's events, optionally
    of one author only."""
    lines = [l.rstrip("\n") for l in open(path, encoding="utf-8") if l.strip()]
    head = [c.strip() for c in lines[0].lstrip("#").split("|")]
    it, ia = head.index("Time"), head.index("Author")
    out = []
    for line in lines[1:]:
        if line.startswith("#"):
            continue
        f = line.split("|")
        if author is None or f[ia].strip() == author:
            out.append(f[it].strip().split("T")[1])
    return out


def followed_events(events, alarms, within):
    """Of `events` (times of day), those followed within `within` seconds by one
    of `alarms` ({number: time of day}), and the distinct alarms that follow one.
    A set is a replay of three hours, so times are compared within the day,
    across midnight if need be."""
    followed, used = [], set()
    for e in events:
        hits = [
            n
            for n, t in alarms.items()
            if 0 <= (day_seconds(t) - day_seconds(e)) % 86400 <= within
        ]
        if hits:
            followed.append(e)
            used.update(hits)
    return {
        "events": len(events),
        "followed": len(followed),
        "distinct_alarms": len(used),
        "followed_events": followed,
    }


def seconds_after(t, origin):
    """t minus origin in seconds, both times of day, taken within half a day
    either way (a set is a replay of three hours)."""
    return (day_seconds(t) - day_seconds(origin) + 43200) % 86400 - 43200


def parse_events(text):
    """An ayzek report's event blocks {number: {time, declaring, pd_first,
    pd_first_after_alarm_s, pd_final, afad, delay, pre_s, warning}}: `afad` is
    the (type, magnitude, origin) of the block's "AFAD" line or None,
    `pre_s` (before, used) or None, `warning` the warning table's rows
    (station, km, picked or predicted, warning in seconds, before or after)."""
    out = {}
    for block in re.split(r"^Event #", text, flags=re.M)[1:]:
        block = block.split("\n\n")[0]
        alarm = BLOCK_ALARM.search(block)
        if not alarm:
            continue
        pd, afad, pre = (x.search(block) for x in (BLOCK_PD, BLOCK_AFAD, BLOCK_PRE_S))
        out[int(re.match(r"\d+", block).group())] = {
            "time": alarm.group(1),
            "declaring": [x.strip() for x in alarm.group(2).split(",")],
            "pd_first": float(pd.group(1)) if pd else None,
            "pd_first_after_alarm_s": float(pd.group(2) or 0) if pd else None,
            "pd_final": float(pd.group(3)) if pd else None,
            "afad": afad.groups()[:3] if afad else None,
            "delay": float(afad.group(4)) if afad else None,
            "pre_s": (int(pre.group(1)), int(pre.group(2))) if pre else None,
            "warning": [
                (st, int(km), source, float(lead), side)
                for st, km, source, lead, side in WARNING.findall(block)
            ],
        }
    return out


def scored_mainshock(first, events):
    """The event (dataset, origin, type, magnitude) of an outcome table that is
    the selection's `first`: the same type, the one-decimal magnitude within
    0.051, the origin within 1 s; the nearest in time, None without one."""
    hms = first["time"].split("T")[1]
    near = [
        (abs(seconds_after(e[1], hms)), e)
        for e in events
        if e[2].upper() == first["type"].upper()
        and abs(float(e[3]) - first["magnitude"]) < 0.051
        and abs(seconds_after(e[1], hms)) <= 1.0
    ]
    return min(near)[1] if near else None


def station_distances(first, coords, stations):
    """[(km, station)] of the stations from the event's epicentre, nearest first."""
    return sorted(
        (haversine_km(first["lat"], first["lon"], *coords[s]), s) for s in stations
    )


def earlier_events_note(first, key, events, joined):
    """The scored events of the set before the mainshock, and the events the
    selection joined to its sequence, as a note ("" when there are none)."""
    hms = first["time"].split("T")[1]
    notes = []
    before = sorted(
        (-seconds_after(e[1], hms), e)
        for e in events
        if e != key and seconds_after(e[1], hms) < 0
    )
    if before:
        gap, e = before[0]
        notes.append(
            f"{len(before)} scored event{'s' * (len(before) > 1)} earlier in the set, "
            f"the last {e[2]} {e[3]} at {e[1]} ({gap:.1f} s before)"
        )
    if joined:
        notes.append(
            f"the selection joined {len(joined)} later event{'s' * (len(joined) > 1)}: "
            + ", ".join(
                f"{j['type']} {j['magnitude']:.1f} at {j['time'].split('T')[1][:11]}"
                for j in joined
            )
        )
    return notes


def mainshock_detector(first, key, detected, pre_s, events, window):
    """One detector's row for a sequence-initiating event: `key` its scored
    event or None, `detected` and `pre_s` the outcome table's, `events` the
    report's blocks (parse_events), `window` the seconds before and after the
    selection's origin within which other alarms are noted."""
    hms = first["time"].split("T")[1]
    row = {
        "status": (
            "not scored" if key is None else "detected" if detected else "missed"
        ),
        "detected_pre_s": None if key is None else pre_s,
    }
    own = sorted(
        (seconds_after(b["time"], hms), n)
        for n, b in events.items()
        if key is not None and b["afad"] == (key[2], key[3], key[1])
    )
    notes = []
    taken = None
    if own:
        taken = own[0][1]
        b = events[taken]
        ahead = [(km, lead, st, src) for st, km, src, lead, side in b["warning"] if side == "before"]
        row.update(
            alarm=taken,
            alarm_time=b["time"],
            alarm_delay_s=b["delay"],
            alarm_after_selection_origin_s=round(own[0][0], 2),
            declaring=b["declaring"],
            pre_s=f"{b['pre_s'][0]} of {b['pre_s'][1]}" if b["pre_s"] else None,
            pd_first=b["pd_first"],
            pd_first_after_alarm_s=b["pd_first_after_alarm_s"],
            pd_final=b["pd_final"],
            warning_stations=len(b["warning"]),
            stations_ahead=len(ahead),
        )
        if ahead:
            km, lead, st, src = min(ahead)  # nearest; at equal distance the smaller
            row.update(
                nearest_ahead_station=st,
                nearest_ahead_km=km,
                nearest_ahead_warning_s=lead,
                nearest_ahead_s_source=src,
            )
        if len(own) > 1:
            notes.append(
                f"{len(own)} alarms credited to it ("
                + ", ".join(f"#{n}" for _, n in own)
                + "), the earliest taken"
            )
    if key is not None and bool(detected) != bool(own):
        notes.append(
            f"outcome table {'detected' if detected else 'missed'}, the report has "
            f"{'an' if own else 'no'} event block credited to it"
        )
    for n, b in sorted(events.items()):
        d = seconds_after(b["time"], hms)
        if n == taken or not -window[0] <= d <= window[1]:
            continue
        notes.append(
            f"alarm #{n} at {b['time']} ({d:+.1f} s) by {', '.join(b['declaring'])} "
            + (
                "credited to it too"
                if key is not None and (d, n) in own
                else (
                    f"credited to {b['afad'][0]} {b['afad'][1]} at {b['afad'][2]}"
                    if b["afad"]
                    else "credited to no catalogue event"
                )
            )
            + (f", Pd final M{b['pd_final']:.1f}" if b["pd_final"] is not None else "")
        )
    row["note"] = notes
    return row


def mainshock_summary(rows):
    """Over one detector's rows (mainshock_detector): the counts and medians
    of the network summary."""
    det = [r for r in rows if r["status"] == "detected"]
    delays = [r["alarm_delay_s"] for r in det if r.get("alarm_delay_s") is not None]
    leads = [
        r["nearest_ahead_warning_s"]
        for r in det
        if r.get("nearest_ahead_warning_s") is not None
    ]
    return {
        "mainshocks": len(rows),
        "scored": sum(r["status"] != "not scored" for r in rows),
        "detected": len(det),
        "detected_pre_s": sum(bool(r["detected_pre_s"]) for r in det),
        "pd_at_alarm": sum(r.get("pd_first_after_alarm_s") == 0 for r in det),
        "pd_any": sum(r.get("pd_final") is not None for r in det),
        "station_ahead": len(leads),
        "median_alarm_delay_s": quantile(delays, 0.5) if delays else None,
        "median_nearest_ahead_warning_s": quantile(leads, 0.5) if leads else None,
    }


def mainshocks(sel, tag, run, labels, coords, det, pre, stations, window):
    """Item 9: the selection's sequence-initiating events with each detector's
    row and the network summary. `det` and `pre` are load_outcomes' tables of
    `detected` and `detected_pre_s`, `stations` {set: the stations that ran}."""
    events = []
    for q in sel["sequences"]:
        if not q["selected"]:
            continue
        s, first = f"{tag}_seq{q['index']:03d}", q["first"]
        scored = [e for e in det[labels[0]][0] if e[0] == s]
        key = scored_mainshock(first, scored)
        dist = station_distances(first, coords, stations[s])
        common = earlier_events_note(first, key, scored, q.get("joined", []))
        if key is None:
            common.insert(0, "no such event in the outcome table")
        ev = {
            "set": s,
            "origin": first["time"],
            "magnitude": first["magnitude"],
            "type": first["type"],
            "depth_km": first["depth_km"],
            "lat": first["lat"],
            "lon": first["lon"],
            "nearest": [{"station": x, "km": km} for km, x in dist[:2]],
            "scored_as": f"{key[1]} {key[2]} {key[3]}" if key else None,
            "note": common,
            "detectors": {},
        }
        for label in labels:
            text = (run / "runs" / label / f"{s}.txt").read_text()
            ev["detectors"][label] = mainshock_detector(
                first,
                key,
                det[label][0].get(key),
                pre[label][0].get(key),
                parse_events(text),
                window,
            )
        events.append(ev)
    return {
        "alarm_window_s": list(window),
        "events": events,
        "summary": {
            label: mainshock_summary([ev["detectors"][label] for ev in events])
            for label in labels
        },
    }


def mainshock_rows(m):
    """Item 9 as flat rows, one per event and detector (MAINSHOCK_FIELDS)."""
    rows = []
    for ev in m["events"]:
        near = ev["nearest"] + [{"station": None, "km": None}] * 2
        for label, d in ev["detectors"].items():
            row = {
                **{k: ev[k] for k in ("set", "origin", "type", "lat", "lon", "scored_as")},
                "magnitude": f"{ev['magnitude']:.2f}",
                "depth_km": f"{ev['depth_km']:.1f}",
                "nearest_station": near[0]["station"],
                "nearest_km": near[0]["km"],
                "second_station": near[1]["station"],
                "second_km": near[1]["km"],
                "detector": label,
                **d,
                "declaring": " ".join(d.get("declaring", [])),
                "note": "; ".join(ev["note"] + d["note"]),
            }
            for k in ("nearest_km", "second_km"):
                if row[k] is not None:
                    row[k] = f"{row[k]:.1f}"
            rows.append({k: "" if row.get(k) is None else row[k] for k in MAINSHOCK_FIELDS})
    return rows


def analyse(
    run, selection, labels, quiet_sets=(), large=4.0, other=None, window=(30.0, 60.0)
):
    run = Path(run)
    sel = json.loads((Path(selection) / "selection.json").read_text())
    tag = f"{sel['group']}_{sel['name']}"
    sequences = [
        f"{tag}_seq{q['index']:03d}" for q in sel["sequences"] if q["selected"]
    ]
    windows = {
        f"{tag}_quiet{j:02d}": w["sequence"]
        for j, w in enumerate(sel["quiet"]["windows"], 1)
    }
    coords = {s["station"]: (s["lat"], s["lon"]) for s in sel["stations"]["usable"]}
    a, b = labels[:2]
    res = {"network": tag, "detectors": labels}

    with open(run / "analysis" / "magnitude_latency.csv") as fh:
        res["paired_latency"] = paired_latency(list(csv.DictReader(fh)), a, b)

    reports = {
        label: [
            (s, kind, parse_report((run / sub / label / f"{s}.txt").read_text()))
            for sub, kind, sets in (
                ("runs", "sequence", sequences),
                ("quiet", "quiet", windows),
            )
            for s in sets
        ]
        for label in labels
    }
    unmatched = {
        label: {
            s: len(rep["unmatched"])
            for s, kind, rep in reports[label]
            if kind == "sequence"
        }
        for label in labels
    }

    det = load_outcomes(run / "outcomes", labels, sequences)
    for label in labels[1:]:
        if set(det[label][0]) != set(det[a][0]):
            sys.exit(f"{label}: not the events of {a} in the outcome tables")
    table = {
        s: {
            label: {**sequence_row(det[label], s), "unmatched": unmatched[label][s]}
            for label in labels
        }
        for s in sequences
    }
    res["sequences"] = table
    res["sequence_comparison"] = more_same_fewer(table, a, b)
    res["network_chance"] = {label: network_chance(table, label) for label in labels}

    with open(run / "quiet" / "quiet_index.csv") as fh:
        hours = {r["name"]: float(r["scored_hours"]) for r in csv.DictReader(fh)}
    counts = {
        label: {
            s: len(rep["alarms"]) for s, kind, rep in reports[label] if kind == "quiet"
        }
        for label in labels
    }
    res["quiet_windows"] = {
        s: {"set": f"{windows[s]:03d}", **{label: counts[label][s] for label in labels}}
        for s in windows
    }
    res["quiet_sets"] = quiet_groups(windows, hours, counts, labels)
    if quiet_sets:
        named = [s for s in windows if windows[s] in quiet_sets]
        rest = [s for s in windows if windows[s] not in quiet_sets]
        res["quiet_named_sets"] = [
            quiet_total(
                "+".join(f"{q:03d}" for q in quiet_sets), named, hours, counts, labels
            ),
            quiet_total("other", rest, hours, counts, labels),
        ]

    pairs = {}
    for s, kind, rep in reports[a]:
        missing = [x for x in rep["stations"] if x not in coords]
        if missing:
            sys.exit(f"{s}: no coordinates in the selection for {', '.join(missing)}")
        pairs[s] = closest_pair({x: coords[x] for x in rep["stations"]})
    res["closest_pairs"] = {
        s: {"pair": list(p[0]), "km": p[1]} if p else None for s, p in pairs.items()
    }
    res["closest_pair_alarms"] = {
        label: pair_shares(reports[label], pairs) for label in labels
    }

    res["magnitude_bands"] = magnitude_bands(det, [a, b])

    res["mainshocks"] = mainshocks(
        sel,
        tag,
        run,
        labels,
        coords,
        det,
        load_outcomes(run / "outcomes", labels, sequences, "detected_pre_s"),
        {s: rep["stations"] for s, kind, rep in reports[a] if kind == "sequence"},
        window,
    )

    res["alarm_magnitudes"] = {}
    for label in labels:
        um = [
            (s, rep["alarms"][n])
            for s, kind, rep in reports[label]
            if kind == "sequence"
            for n in rep["unmatched"]
        ]
        qu = [
            (s, al)
            for s, kind, rep in reports[label]
            if kind == "quiet"
            for al in rep["alarms"].values()
        ]
        res["alarm_magnitudes"][label] = {
            "unmatched": alarm_magnitudes(um, large),
            "quiet": alarm_magnitudes(qu, large),
        }

    if other:
        name = f"{tag}_seq{other['set']:03d}"
        events = read_other_events(other["file"], other.get("author"))
        res["other_events"] = {
            "set": name,
            "author": other.get("author"),
            "within_s": other["within"],
            "detectors": {},
        }
        for label in labels:
            rep = next(r for s, kind, r in reports[label] if s == name)
            alarms = {n: rep["alarms"][n]["time"] for n in rep["unmatched"]}
            res["other_events"]["detectors"][label] = followed_events(
                events, alarms, other["within"]
            )
    return res


def report(res, large):
    labels = res["detectors"]
    a, b = labels[:2]
    short = lambda s: s.rsplit("_", 1)[-1]
    print(f"== {res['network']}: exploratory, outside the frozen protocol")

    p = res["paired_latency"]
    print(f"\n1. Paired latency, {a} minus {b}: {p['events']} events both detected")
    if p["events"]:
        print(
            f"   median {p['median']:+.2f} s (interquartile range {p['q1']:+.2f} to {p['q3']:+.2f}); "
            f"{b} first {p[f'{b}_first']}, {a} first {p[f'{a}_first']}, ties {p['ties']}"
        )

    print("\n2. Per sequence: detected / chance level c / corrected / unmatched alarms")
    print(
        f"   {'set':8s} {'scored':>6s}  "
        + "  ".join(f"{label:>28s}" for label in labels)
    )
    for s, row in res["sequences"].items():
        print(
            f"   {short(s):8s} {row[a]['scored']:6d}  "
            + "  ".join(
                f"{r['detected']:5d} {r['chance']:6.3f} {r['corrected']:7.1f} {r['unmatched']:6d}"
                for r in (row[label] for label in labels)
            )
        )
    c = res["sequence_comparison"]
    print(
        f"   sequences where {a} detected more / the same / fewer events than {b}: "
        f"{c['more']} / {c['same']} / {c['fewer']}"
    )

    print(
        "\n3. Network: detected, corrected (sum over sequences), effective c = (raw - corr) / (n - corr)"
    )
    for label in labels:
        r = res["network_chance"][label]
        print(
            f"   {label:22s} {r['detected']:4d} of {r['scored']}  corrected {r['corrected']:6.1f}  "
            f"effective c {r['effective_chance']:.3f}"
        )

    print(
        "\n4. Quiet windows by station set: windows, scored hours; alarms (per day) of "
        + ", ".join(labels)
    )
    for g in res["quiet_sets"] + res.get("quiet_named_sets", []):
        print(
            f"   set {g['set']:12s} {g['windows']:3d} {g['hours']:6.1f} h  "
            + "  ".join(
                f"{g['alarms'][label]:4d} ({g['per_day'][label]:5.1f})"
                for label in labels
            )
        )
    print("   per window (station set): alarms of " + ", ".join(labels))
    for s, w in res["quiet_windows"].items():
        print(
            f"   {short(s):8s} ({w['set']})  "
            + "  ".join(f"{w[label]:4d}" for label in labels)
        )

    print(
        "\n5. Closest station pair per set, and the alarms declared by exactly that pair"
    )
    seq_pairs = collections.Counter()
    for s, p in res["closest_pairs"].items():
        print(
            f"   {short(s):8s} "
            + (
                f"{p['pair'][0]}-{p['pair'][1]} {p['km']:.1f} km"
                if p
                else "fewer than two stations"
            )
        )
        if p and "_seq" in s:
            seq_pairs[f"{p['pair'][0]}-{p['pair'][1]}"] += 1
    print(
        "   sequence sets per pair: "
        + ", ".join(f"{k} {v}" for k, v in seq_pairs.most_common())
    )
    for label in labels:
        print(
            f"   {label:22s} "
            + ", ".join(
                f"{k} {v['by_closest_pair']} of {v['alarms']}"
                for k, v in res["closest_pair_alarms"][label].items()
            )
        )

    print(
        f"\n6. Detection by catalogue-magnitude band (the outcome tables' one-decimal magnitudes): "
        f"scored, {a} vs {b}"
    )
    for r in res["magnitude_bands"]:
        print(f"   {r['band']:4s} {r['scored']:4d}  {r[a]:4d} vs {r[b]}")

    print(
        f"\n7. Pd estimates of unmatched and quiet alarms: alarms / with an estimate / final >= {large:g}"
    )
    for label in labels:
        for kind, m in res["alarm_magnitudes"][label].items():
            print(
                f"   {label:22s} {kind:9s} {m['alarms']:4d} / {m['with_pd']} / {m['large']}"
            )
            for x in m["large_cases"]:
                print(f"      {x['set']} {x['time']} M{x['pd_final']:.1f}")

    if "other_events" in res:
        o = res["other_events"]
        print(
            f"\n8. Other bulletin's events in {o['set']}"
            f"{' (author ' + o['author'] + ')' if o['author'] else ''}: followed within "
            f"{o['within_s']:g} s by an unmatched alarm / distinct alarms"
        )
        for label, m in o["detectors"].items():
            print(
                f"   {label:22s} {m['followed']} of {m['events']} / {m['distinct_alarms']}"
            )

    m = res["mainshocks"]
    fmt = lambda v, f: "none" if v is None else format(v, f)
    print(
        "\n9. Sequence-initiating events (the selection's first event of each sequence): "
        "origin, magnitude, depth, epicentre, the two nearest stations (epicentral km, "
        "selection coordinates); per detector the outcome table's status and the "
        "report's block; notes on other alarms from "
        f"{m['alarm_window_s'][0]:g} s before to {m['alarm_window_s'][1]:g} s after the origin"
    )
    for ev in m["events"]:
        print(
            f"   {short(ev['set']):7s} {ev['origin'].replace('T', ' ')[:23]}  "
            f"{ev['type']} {ev['magnitude']:.2f}  {ev['depth_km']:.1f} km  "
            f"{ev['lat']:.3f} {ev['lon']:.3f}  "
            + ", ".join(f"{x['station']} {x['km']:.1f}" for x in ev["nearest"])
            + f"  scored as {ev['scored_as'] or 'nothing'}"
        )
        for label, d in ev["detectors"].items():
            line = f"      {label:22s} {d['status']:10s}"
            if "alarm" in d:
                pd = (
                    "no Pd"
                    if d["pd_final"] is None
                    else f"Pd M{d['pd_first']:.1f} "
                    + (
                        "at the alarm"
                        if d["pd_first_after_alarm_s"] == 0
                        else f"{d['pd_first_after_alarm_s']:g} s after"
                    )
                    + f", final M{d['pd_final']:.1f}"
                )
                ahead = (
                    f"nearest {d['nearest_ahead_station']} {d['nearest_ahead_km']} km "
                    f"{d['nearest_ahead_warning_s']:+.1f} s ({d['nearest_ahead_s_source']})"
                    if "nearest_ahead_station" in d
                    else "none ahead"
                )
                line += (
                    f" alarm {d['alarm_delay_s']:.1f} s after origin by {', '.join(d['declaring'])}; "
                    f"pre-S {d['pre_s'] or 'not printed'} (outcome table {fmt(d['detected_pre_s'], 'd')}); {pd}; "
                    f"before S {d['stations_ahead']} of {d['warning_stations']} stations, {ahead}"
                )
            print(line)
            for note in d["note"]:
                print(f"         {note}")
        for note in ev["note"]:
            print(f"      note: {note}")
    print(
        "   network: mainshocks / scored / detected / detected before S / Pd at the alarm "
        "(any Pd) / a station before S; median alarm delay; median warning at the nearest "
        "station before S"
    )
    for label in labels:
        s = m["summary"][label]
        print(
            f"   {label:22s} {s['mainshocks']} / {s['scored']} / {s['detected']} / "
            f"{s['detected_pre_s']} / {s['pd_at_alarm']} ({s['pd_any']}) / {s['station_ahead']}; "
            f"{fmt(s['median_alarm_delay_s'], '.2f')} s; "
            f"{fmt(s['median_nearest_ahead_warning_s'], '.2f')} s"
        )


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--run", required=True, help="a tools/run_heldout.py output directory"
    )
    ap.add_argument(
        "--selection", required=True, help="the heldout select directory of that run"
    )
    ap.add_argument(
        "--detectors",
        default="learned,stalta_eqfa,phasenet_causal_eqfa",
        help="labels, comma-separated; the first two are the pair compared",
    )
    ap.add_argument(
        "--quiet-sets",
        default="",
        metavar="N,N,...",
        help="station sets (sequence numbers) also totalled together, with the rest",
    )
    ap.add_argument(
        "--large",
        type=float,
        default=4.0,
        help="final Pd estimate from which an unmatched or quiet alarm is listed",
    )
    ap.add_argument(
        "--other-events",
        metavar="FILE",
        help="FDSN text listing of another bulletin's events (item 8)",
    )
    ap.add_argument(
        "--other-set",
        type=int,
        metavar="N",
        help="the sequence set those events lie in",
    )
    ap.add_argument(
        "--other-author", help="keep the listing's entries of this author only"
    )
    ap.add_argument(
        "--follow-seconds",
        type=float,
        default=60.0,
        help="an unmatched alarm this long after an event follows it",
    )
    ap.add_argument(
        "--alarm-window",
        default="30,60",
        metavar="B,A",
        help="item 9 notes the other alarms from B seconds before to A after a mainshock's origin",
    )
    ap.add_argument(
        "--mainshocks",
        metavar="FILE",
        help="also write item 9's table here (CSV, one row per event and detector)",
    )
    ap.add_argument("--json", help="also write the figures here")
    a = ap.parse_args()
    if bool(a.other_events) != (a.other_set is not None):
        sys.exit("--other-events and --other-set go together")
    other = (
        {
            "file": a.other_events,
            "set": a.other_set,
            "author": a.other_author,
            "within": a.follow_seconds,
        }
        if a.other_events
        else None
    )
    labels = a.detectors.split(",")
    if len(labels) < 2:
        sys.exit("--detectors: at least two labels")
    res = analyse(
        a.run,
        a.selection,
        labels,
        [int(x) for x in a.quiet_sets.split(",") if x],
        a.large,
        other,
        tuple(float(x) for x in a.alarm_window.split(",")),
    )
    report(res, a.large)
    if a.mainshocks:
        with open(a.mainshocks, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=MAINSHOCK_FIELDS)
            w.writeheader()
            w.writerows(mainshock_rows(res["mainshocks"]))
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=1, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
