#!/usr/bin/env python3
"""Classify network alarms with an independent picker and the S-P test.

    uv run --with seisbench --with obspy python tools/classify_alarms.py \\
        --runs data/runs/quiet100 --data data/quiet100 --weights ethz --out learned_ethz.csv

For every alarm in RUN/*.txt (the replay reports), PhaseNet or EQTransformer
(SeisBench, --model) picks P
and S on the first two declaring stations, from 150 s before to 120 s after the
alarm, at probability >= 0.5. Per station the P/S pair with the highest summed
probability is kept, among P picks from 90 s before to 5 s after the alarm and S
picks within 90 s after the P. The pair of stations then goes through
tools/misfire_review.py's locate_check unchanged: the two S-P distances
(VP 6.0, VS 3.5 km/s) must fit the station separation and predict the observed
difference of the two P picks to within 1.5 s.

Classes: confirmed earthquake; inconsistent (picks at both stations, test
fails); no picks (no P/S pair at one or both stations). The picker is the
pipeline's competitor, not its own, so the classification does not rest on the
pipeline's picker; the same rule applies to every detector's alarms.
"""
import argparse
import csv
import glob
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from misfire_review import locate_check, read_record  # noqa: E402

ALARM = re.compile(r"^ALARM\s+#\d+\s+(\S+)\s+earthquake detected by ([A-Z0-9, ]+?)\s*(?:\(|$)", re.M)
PROB = 0.5


def alarms(txt):
    """[(epoch, [stations])] of one report; alarm times carry no date."""
    from obspy import UTCDateTime
    name = Path(txt).stem
    day, hhmm = name.split("_")[1:3]
    start = UTCDateTime(f"{day[:4]}-{day[4:6]}-{day[6:]}T{hhmm[:2]}:{hhmm[2:]}:00").timestamp
    out, prev = [], start
    for hms, st in ALARM.findall(Path(txt).read_text()):
        t = UTCDateTime(f"{day[:4]}-{day[4:6]}-{day[6:]}T{hms}").timestamp
        while t < prev:
            t += 86400
        prev = t
        out.append((t, [s.strip() for s in st.split(",")]))
    return out


def best_pair(picks, t):
    """The P/S pair with the highest summed probability, or None."""
    ps = [(p.peak_time.timestamp, p.peak_value) for p in picks if p.phase == "P"
          and t - 90 <= p.peak_time.timestamp <= t + 5]
    ss = [(p.peak_time.timestamp, p.peak_value) for p in picks if p.phase == "S"]
    best = None
    for pt, pp in ps:
        for st, sp in ss:
            if 0 < st - pt <= 90 and (best is None or pp + sp > best[2] + best[3]):
                best = (pt, st, pp, sp)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", required=True, help="directory of NAME.txt and NAME.rec")
    ap.add_argument("--data", required=True, help="directory of NAME/STATION.mseed")
    ap.add_argument("--model", choices=["phasenet", "eqtransformer"], default="phasenet")
    ap.add_argument("--weights", default="ethz")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    import seisbench.models as sbm
    import torch
    from obspy import UTCDateTime, read
    model = {"phasenet": sbm.PhaseNet, "eqtransformer": sbm.EQTransformer}[a.model].from_pretrained(a.weights)
    if torch.cuda.is_available():
        model.cuda()

    rows = []
    for txt in sorted(glob.glob(f"{a.runs}/*.txt")):
        name = Path(txt).stem
        stations = read_record(f"{a.runs}/{name}.rec")[0]
        for t, declared in alarms(txt):
            A, B = declared[:2]
            picks = {}
            for s in (A, B):
                st = read(f"{a.data}/{name}/{s}.mseed", starttime=UTCDateTime(t - 150),
                          endtime=UTCDateTime(t + 120))
                st.merge(fill_value=0)
                out = model.classify(st, P_threshold=PROB, S_threshold=PROB)
                pair = best_pair(out.picks, t)
                if pair:
                    picks[s] = [pair]
            if len(picks) < 2:
                verdict, detail = "no picks", f"P/S pair at {', '.join(picks) or 'neither station'} only"
            else:
                dt_obs = picks[B][0][0] - picks[A][0][0]
                verdict, detail = locate_check(picks, stations, A, B, t, dt_obs, min_prob=PROB - 1e-9)
                verdict = {"earthquake": "confirmed earthquake", "misfire": "inconsistent",
                           "unclassified": "no picks"}[verdict]
            rows.append({"record": name, "alarm_utc": UTCDateTime(t).strftime("%Y-%m-%d %H:%M:%S.%f")[:-4],
                         "stations": f"{A} {B}", "class": verdict, "detail": detail})
            print(f"{rows[-1]['alarm_utc']} {A}-{B}: {verdict}  {detail}", flush=True)

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    c = Counter(r["class"] for r in rows)
    print(f"\n{len(rows)} alarms: " + ", ".join(f"{k} {v}" for k, v in c.most_common()))


if __name__ == "__main__":
    main()
