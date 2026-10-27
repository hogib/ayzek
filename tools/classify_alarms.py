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
fails); no picks (no P/S pair at one or both stations).

--pairs any (referee round 5, M7; protocol in ayzek_paper revisions2.md R17.0)
widens the test to every pair of candidate stations: the declaring stations and
every station with a detection (reported P) compatible with a declaring
station's detection under the pipeline's association rule, |dt| <= d / vp +
slack. The declaring station's detection is its latest with reported P in the
120 s before the alarm. An alarm is confirmed if any pair passes; otherwise
inconsistent if any pair has picks at both stations; otherwise no picks. The
picker is the
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
from misfire_review import hav, locate_check, read_record  # noqa: E402

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


def associated(declared, det, stations, t, a):
    """The declaring stations, then every other station with a detection that the
    pipeline's association rule makes compatible with a declaring station's."""
    ref = {}
    for s in declared:
        before = [x for x in det.get(s, []) if t - 120 <= x <= t]
        if before:
            ref[s] = before[-1]
    out = list(declared)
    for s, times in det.items():
        if s in out:
            continue
        for d, td in ref.items():
            w = hav(*stations[s], *stations[d]) / a.vp + a.slack
            if any(abs(x - td) <= w for x in times):
                out.append(s)
                break
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", required=True, help="directory of NAME.txt and NAME.rec")
    ap.add_argument("--data", required=True, help="directory of NAME/STATION.mseed")
    ap.add_argument("--model", choices=["phasenet", "eqtransformer"], default="phasenet")
    ap.add_argument("--weights", default="ethz")
    ap.add_argument("--pairs", choices=["declared", "any"], default="declared",
                    help="the first two declaring stations (the pre-registered rule), or every "
                         "pair of stations associated with the alarm")
    ap.add_argument("--vp", type=float, default=6.0)
    ap.add_argument("--slack", type=float, default=3.0)
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
        stations, det = read_record(f"{a.runs}/{name}.rec")[:2]
        for t, declared in alarms(txt):
            cand = list(declared[:2]) if a.pairs == "declared" else associated(declared, det, stations, t, a)
            picks = {}
            for s in cand:
                st = read(f"{a.data}/{name}/{s}.mseed", starttime=UTCDateTime(t - 150),
                          endtime=UTCDateTime(t + 120))
                st.merge(fill_value=0)
                out = model.classify(st, P_threshold=PROB, S_threshold=PROB)
                pair = best_pair(out.picks, t)
                if pair:
                    picks[s] = [pair]
            tried = []
            for i, A in enumerate(cand):
                for B in cand[i + 1:]:
                    if A not in picks or B not in picks:
                        continue
                    dt_obs = picks[B][0][0] - picks[A][0][0]
                    v, d = locate_check(picks, stations, A, B, t, dt_obs, min_prob=PROB - 1e-9)
                    tried.append((v, f"{A}-{B}: {d}"))
            if any(v == "earthquake" for v, _ in tried):
                verdict, detail = "confirmed earthquake", next(d for v, d in tried if v == "earthquake")
            elif any(v == "misfire" for v, _ in tried):
                verdict, detail = "inconsistent", next(d for v, d in tried if v == "misfire")
            elif tried:
                verdict, detail = "no picks", tried[0][1]
            else:
                verdict, detail = "no picks", f"P/S pair at {', '.join(picks) or 'no station'} only"
            rows.append({"record": name, "alarm_utc": UTCDateTime(t).strftime("%Y-%m-%d %H:%M:%S.%f")[:-4],
                         "stations": " ".join(cand), "class": verdict, "detail": detail})
            print(f"{rows[-1]['alarm_utc']} {'-'.join(cand)}: {verdict}  {detail}", flush=True)

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    c = Counter(r["class"] for r in rows)
    print(f"\n{len(rows)} alarms: " + ", ".join(f"{k} {v}" for k, v in c.most_common()))


if __name__ == "__main__":
    main()
