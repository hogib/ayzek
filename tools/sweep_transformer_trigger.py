"""Re-scores the transformer's trigger rule from saved token streams, without
replaying: rising edges only, and the dt restart over a grid of settings.

    build-release/app/ayzek --speed 0 --detector transformer --no-dt-reset \\
        --no-pick --no-magnitude --scores cmp/scores_marmara_ko data/marmara_ko/*.mseed
    uv run --project tools python tools/sweep_transformer_trigger.py \\
        cmp/scores_demo_ko:tests/catalogs/demo_ko.csv \\
        cmp/scores_marmara_ko:tests/catalogs/marmara_ko.csv cmp/scores_quiet_ko

Each argument is a --scores directory (one CSV per station: token_end,
probability, dt), optionally with the catalogue to score against after a
colon. The rule is a copy of src/pipeline/trigger.hpp, with the processor's
15 s minimum between triggers. An arrival is P = origin + hypot(distance,
depth) / 6 km/s, visible by compare_detectors.py's rule, and counts as
detected by a trigger decided in [P - 3 s, P + 20 s]; a trigger with no
arrival in [P - 3 s, P + 60 s] is unmatched. That P is cruder than
compare_detectors.py's TauP + AIC, so use this to choose settings and
compare_detectors.py for the verdict.

For each missed arrival it also reports why: p never reached the threshold,
or it did while the trigger was not re-armed and dt did not restart.
"""
import argparse
import datetime as dtm
import glob
import itertools
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
VISIBILITY = ((50.0, 0.0), (150.0, 2.0), (400.0, 3.0), (1500.0, 4.5))
GATE_S = 15.0


def distance_km(a, b, c, d):
    r = np.pi / 180
    x = np.sin((c - a) * r / 2) ** 2 + np.cos(a * r) * np.cos(c * r) * np.sin((d - b) * r / 2) ** 2
    return 2 * 6371 * np.arcsin(np.sqrt(x))


def trigger(t, p, dt, thr, rel, below, frm, ntok, reset=True, max_dt=10.0):
    """(decision time, P, kind) per trigger, as TokenTrigger plus the gate."""
    out, armed, peak, low, last = [], True, 0.0, 0, -1e18
    for i in range(len(p)):
        if p[i] < rel:
            armed, peak, low = True, 0.0, 0
            continue
        if p[i] < thr:
            peak, low = max(peak, dt[i]), 0
            continue
        kind = None
        if armed:
            kind = "rising"
        elif reset and peak >= frm and dt[i] <= below:
            low += 1
            if low >= ntok:
                kind = "dt restart"
        else:
            low = 0
        armed = False
        if kind is None:
            peak = max(peak, dt[i])
            continue
        peak, low = dt[i], 0
        P = t[i] - min(dt[i], max_dt)
        if P - last >= GATE_S:
            last = P
            out.append((t[i], P, kind))
    return out


def load(scores, catalog, stations):
    streams = {}
    files = sorted(glob.glob(f"{scores}/*.csv"))
    if not files:
        raise SystemExit(f"{scores}: no score CSVs (ayzek --scores {scores} ...)")
    for f in files:
        s = pd.read_csv(f)
        if "dt" not in s:
            raise SystemExit(f"{f}: no dt column; re-run ayzek --scores with this build")
        if len(s):
            streams[Path(f).stem] = (s.token_end.values, s.probability.values, s.dt.values)
    arrivals = {code: np.array([]) for code in streams}
    if catalog:
        c = pd.read_csv(catalog, encoding="utf-8-sig")
        c["t"] = [dtm.datetime.strptime(d, "%d/%m/%Y %H:%M:%S")
                  .replace(tzinfo=dtm.timezone.utc).timestamp() for d in c.Date]
        for code, (t, _, _) in streams.items():
            if code not in stations:
                continue
            a = []
            for _, e in c.iterrows():
                d = distance_km(e.Latitude, e.Longitude, *stations[code])
                if not any(d <= r and e.Magnitude >= m for r, m in VISIBILITY):
                    continue
                P = e.t + np.hypot(d, e.Depth) / 6.0
                if t[0] + 30 < P < t[-1] - 20:
                    a.append(P)
            arrivals[code] = np.array(a)
    return streams, arrivals


def score(streams, arrivals, thr, rel, **rule):
    r = dict(arrivals=0, detected=0, never_reached=0, not_rearmed=0, triggers=0,
             restarts=0, unmatched=0, hours=0.0)
    for code, (t, p, dt) in streams.items():
        trig = trigger(t, p, dt, thr, rel, **rule)
        dec = np.array([x[0] for x in trig])
        a = arrivals[code]
        r["triggers"] += len(trig)
        r["restarts"] += sum(k == "dt restart" for *_, k in trig)
        r["hours"] += len(t) * (t[1] - t[0]) / 3600 if len(t) > 1 else 0
        for P in a:
            r["arrivals"] += 1
            if ((dec >= P - 3) & (dec <= P + 20)).any():
                r["detected"] += 1
            elif p[(t >= P - 3) & (t <= P + 20)].max() < thr:
                r["never_reached"] += 1
            else:
                r["not_rearmed"] += 1
        r["unmatched"] += sum(1 for x in dec if not ((a - 3 <= x) & (x <= a + 60)).any())
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sets", nargs="+", help="SCORES_DIR[:CATALOG.csv]")
    ap.add_argument("--threshold", type=float, required=True,
                    help="the run's threshold (the ayzek header line prints it)")
    ap.add_argument("--release", type=float, default=None, help="default: threshold / 2")
    ap.add_argument("--below", default="1,2,3,4")
    ap.add_argument("--from", dest="frm", default="5,7,9")
    ap.add_argument("--tokens", default="1,2,3")
    a = ap.parse_args()
    rel = a.release if a.release is not None else a.threshold / 2
    stations = {l.split(",")[0]: (float(l.split(",")[1]), float(l.split(",")[2]))
                for l in (ROOT / "models/stations.csv").read_text().splitlines()[1:]}
    data = {}
    for s in a.sets:
        scores, _, cat = s.partition(":")
        if not scores or (cat and not os.path.isfile(cat)):
            raise SystemExit(f"{s}: expected SCORES_DIR[:CATALOG.csv], no space after the colon")
        data[os.path.basename(scores.rstrip("/"))] = load(scores, cat or None, stations)
    floats = lambda s: [float(x) for x in s.split(",")]
    grid = [dict(below=0.0, frm=0.0, ntok=1, reset=False)] + [
        dict(below=b, frm=f, ntok=int(k)) for b, f, k in
        itertools.product(floats(a.below), floats(a.frm), floats(a.tokens))]
    rows = []
    for g in grid:
        row = {"rule": "rising edges only" if not g.get("reset", True)
               else f"dt <= {g['below']:g} after {g['frm']:g}, {g['ntok']} tok"}
        for name, (streams, arrivals) in data.items():
            r = score(streams, arrivals, a.threshold, rel, **g)
            if r["arrivals"]:
                row[f"{name} detected"] = f"{r['detected']}/{r['arrivals']}"
                row[f"{name} missed: p<thr / not re-armed"] = f"{r['never_reached']} / {r['not_rearmed']}"
            row[f"{name} unmatched/h"] = round(r["unmatched"] / r["hours"], 2) if r["hours"] else None
            row[f"{name} restarts"] = r["restarts"]
        rows.append(row)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
