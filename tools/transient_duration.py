"""How long the signal behind each alarm lasts: spikes against earthquakes.

    uv run --with obspy --with scipy python tools/transient_duration.py \\
        --runs data/runs/quiet100_stalta --data data/quiet100 --out durations.csv

For every alarm in the replay reports (RUN/*.txt) and each station that
declared it, this takes the station's detection nearest before the alarm from
the recording (RUN/*.rec, D lines: station, P time, trigger time, ...) and measures, on
the three components band-passed 1-20 Hz, how long the envelope stays above
four times the noise RMS of the 25 s before P, within [P - 2 s, P + 20 s].

An instrument spike or a glitch lasts well under a second, a cultural
transient a second or two, and an earthquake recorded at a few tens of
kilometres several seconds or more. The classes written are a triage, not a
verdict: spike (< 0.5 s), transient (0.5-3 s), extended (>= 3 s).
"""
import argparse
import csv
import glob
import json
import re
from pathlib import Path

import numpy as np
from obspy import UTCDateTime, read
from scipy.signal import butter, hilbert, sosfiltfilt

ALARM = re.compile(r"^ALARM\s+#\d+\s+(\S+)\s+earthquake detected by ([A-Z0-9, ]+?)\s*(?:\(|$)", re.M)


def detections(rec):
    """{station: [P time, ...]} from the D lines of a recording."""
    out = {}
    for line in open(rec):
        f = line.split("\t")
        if f[0] == "D":
            out.setdefault(f[1], []).append(float(f[2]))  # the reported P time
    return out


def duration(path, p, fs=100.0):
    st = read(str(path), starttime=UTCDateTime(p - 40), endtime=UTCDateTime(p + 25))
    st.merge(fill_value=0)
    sos = butter(4, [1, 20], "bandpass", fs=fs, output="sos")
    env = None
    for tr in st:
        x = tr.data.astype(np.float64)
        x -= x.mean()
        e = np.abs(hilbert(sosfiltfilt(sos, x)))
        env = e if env is None else np.maximum(env[: len(e)], e[: len(env)])
    if env is None:
        return None, None
    t0 = st[0].stats.starttime.timestamp
    i = lambda t: int(round((t - t0) * fs))
    noise = env[max(0, i(p - 27)): i(p - 2)]
    sig = env[i(p - 2): i(p + 20)]
    if len(noise) < 100 or len(sig) < 100:
        return None, None
    level = 4 * np.sqrt(np.mean(noise ** 2))
    return float(np.sum(sig > level) / fs), float(sig.max() / max(level / 4, 1e-9))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", required=True, help="directory of NAME.txt and NAME.rec")
    ap.add_argument("--data", required=True, help="directory of NAME/STATION.mseed")
    ap.add_argument("--out", required=True)
    ap.add_argument("--json", help="also write the counts here, for tools/reproduce_report.py")
    a = ap.parse_args()

    rows = []
    for txt in sorted(glob.glob(f"{a.runs}/*.txt")):
        name = Path(txt).stem
        det = detections(f"{a.runs}/{name}.rec")
        _, day, hhmm = name.split("_")
        start = UTCDateTime(f"{day[:4]}-{day[4:6]}-{day[6:]}T{hhmm[:2]}:{hhmm[2:]}:00").timestamp
        prev = start
        for hms, stations in ALARM.findall(open(txt).read()):
            # Alarm times carry no date; a window can run past midnight.
            t = UTCDateTime(f"{day[:4]}-{day[4:6]}-{day[6:]}T{hms}").timestamp
            while t < prev:
                t += 86400
            prev = t
            for s in [x.strip() for x in stations.split(",")][:2]:
                ps = [p for p in det.get(s, []) if t - 60 <= p <= t]
                if not ps:
                    continue
                d, snr = duration(Path(a.data) / name / f"{s}.mseed", max(ps))
                if d is None:
                    continue
                cls = "spike" if d < 0.5 else "transient" if d < 3 else "extended"
                rows.append({"record": name, "alarm_utc": f"{t:.2f}", "station": s,
                             "p_utc": f"{max(ps):.2f}", "seconds_above": f"{d:.2f}",
                             "peak_over_noise": f"{snr:.1f}", "class": cls})
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print(f"{len(rows)} station detections behind the alarms")
    for s, c in sorted(Counter((r["station"], r["class"]) for r in rows).items()):
        print(f"  {s[0]:5s} {s[1]:9s} {c}")
    print("  all  ", dict(Counter(r["class"] for r in rows)))
    alarms = {}
    for r in rows:
        alarms.setdefault((r["record"], r["alarm_utc"]), set()).add(r["class"])
    if a.json:
        summary = {"station detections": len(rows), "alarms": len(alarms),
                   "alarms with an extended signal": sum("extended" in c for c in alarms.values()),
                   **{k: v for k, v in Counter(r["class"] for r in rows).items()}}
        Path(a.json).write_text(json.dumps({"summary": summary}, indent=1) + "\n")


if __name__ == "__main__":
    main()
