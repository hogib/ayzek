#!/usr/bin/env python3
"""Separate probable uncatalogued earthquakes from genuine misfires.

    tools/misfire_review.py --records data/runs/quiet100/*.rec \
        --csv data/runs/quiet100/misfires.csv

An alarm on catalogue-quiet data is not proof of a false positive: the
catalogue is incomplete for small events, so some alarms are real earthquakes
nobody listed. This tool asks, for each alarm, how surprising its coincidence
is — because that is the one thing that distinguishes the two.

ayzek declares an event when `--min-stations` stations detect within the
inter-station P travel time plus a slack (`d / vp + slack`, network.hpp), so
every alarm already passes a moveout test. What that test cannot say is how
often independent station noise would pass it by chance. Two measures of that
are reported:

**Analytic.** With station i detecting at rate λ_i, the expected number of
accidental pairs in T seconds is Σ λ_i λ_j w_ij T, where w_ij = 2(d_ij/vp +
slack) is the coincidence window. Per alarm, `p_chance` is the probability
that the partner station fired inside that window by chance.

**Time-slide.** Each station's detection series is shifted by a multiple of
`--slide` seconds, wrapped inside the run, and the coincidence count is
recomputed. Real coincidences are destroyed, each station's own rate and
clustering are preserved, and the spread over slides is an empirical null.
This is the same construction as the network coincidence analysis.

Detections at one station closer together than `--coda` seconds are counted
once, as the pipeline does.

The verdict column is a triage aid, not a conclusion: `review` means the
coincidence is too tight and too improbable to dismiss as noise, and the
waveform deserves eyes.
"""

import argparse
import csv
import glob
import itertools
import math
import statistics
import sys
from pathlib import Path


def hav(a, b, c, d):
    r = 6371.0
    p = math.radians
    x = (math.sin(p(c - a) / 2) ** 2
         + math.cos(p(a)) * math.cos(p(c)) * math.sin(p(d - b) / 2) ** 2)
    return 2 * r * math.asin(math.sqrt(x))


def read_record(path):
    """-> (stations, detections, span, picks). `picks` is {code: [(p, s, p_prob,
    s_prob), …]} — the S arrival is what makes a coincidence testable."""
    stations, det, span, picks = {}, {}, None, {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] == "S" and len(p) >= 4:
                stations[p[1]] = (float(p[2]), float(p[3]))
            elif p[0] == "T" and len(p) >= 3:
                span = (float(p[1]), float(p[2]))
            elif p[0] == "D" and len(p) >= 3:
                det.setdefault(p[1], []).append(float(p[2]))
            elif p[0] == "P" and len(p) >= 7:
                pt, st_, pp, sp = (float(p[3]), float(p[4]), float(p[5]), float(p[6]))
                if st_ > pt:
                    picks.setdefault(p[1], []).append((pt, st_, pp, sp))
    for v in det.values():
        v.sort()
    return stations, det, span, picks


def cluster(times, coda):
    """One entry per burst: the pipeline treats later detections within
    `coda` seconds at the same station as the same event."""
    out = []
    for t in times:
        if not out or t - out[-1] > coda:
            out.append(t)
    return out


def coincidences(det, stations, vp, slack, shifts=None, span=None):
    """Pairs (i, j) detecting within d_ij/vp + slack. `shifts` moves each
    station's series by that many seconds, wrapped inside `span`."""
    series = {}
    for s, ts in det.items():
        off = (shifts or {}).get(s, 0.0)
        if off and span:
            lo, hi = span
            length = hi - lo
            ts = sorted(lo + ((t - lo + off) % length) for t in ts)
        series[s] = ts
    found = []
    for a, b in itertools.combinations(sorted(series), 2):
        if a not in stations or b not in stations:
            continue
        w = hav(*stations[a], *stations[b]) / vp + slack
        ta, tb = series[a], series[b]
        j = 0
        for t in ta:
            while j < len(tb) and tb[j] < t - w:
                j += 1
            k = j
            while k < len(tb) and tb[k] <= t + w:
                found.append((min(t, tb[k]), a, b, tb[k] - t, w))
                k += 1
    return found


VP, VS = 6.0, 3.5


def locate_check(picks, stations, a, b, t, dt_obs, min_prob=0.5, tol=1.5):
    """Can one source explain both stations' picks?

    S−P at a station gives its distance, d = (t_S − t_P) / (1/vs − 1/vp). Two
    such distances over-determine the problem: they must satisfy the triangle
    inequality against the station separation, and they must predict the
    observed difference in P arrival times. Correlated noise does not produce a
    confident S phase at a consistent S−P at two stations *and* an
    inter-station arrival difference that matches the implied distances.

    Returns (verdict, detail).
    """
    k = 1 / VS - 1 / VP
    best = {}
    for code in (a, b):
        for pt, st_, pp, sp in picks.get(code, ()):
            if abs(pt - t) > 90:
                continue
            if code not in best or pp + sp > best[code][2] + best[code][3]:
                best[code] = (pt, st_, pp, sp)
    if len(best) < 2:
        return "unclassified", "no P+S pick at both stations"
    da = (best[a][1] - best[a][0]) / k
    db = (best[b][1] - best[b][0]) / k
    sep = hav(*stations[a], *stations[b])
    geom = abs(da - db) <= sep + 5 and sep <= da + db + 5
    pred = (da - db) / VP
    err = min(abs(pred - dt_obs), abs(-pred - dt_obs))
    conf = (min(best[a][2], best[b][2]) > min_prob
            and min(best[a][3], best[b][3]) > min_prob)
    detail = (f"da={da:.0f}km db={db:.0f}km sep={sep:.0f}km "
              f"dt_pred={pred:+.1f} err={err:.1f}s "
              f"P={best[a][2]:.2f}/{best[b][2]:.2f} S={best[a][3]:.2f}/{best[b][3]:.2f}")
    if not geom:
        return "misfire", "no single source fits: " + detail
    if not conf:
        return "unclassified", "picks not confident: " + detail
    if err > tol:
        return "misfire", "moveout inconsistent: " + detail
    return "earthquake", detail


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--records", nargs="+", required=True)
    ap.add_argument("--vp", type=float, default=6.0)
    ap.add_argument("--slack", type=float, default=3.0)
    ap.add_argument("--coda", type=float, default=40.0)
    ap.add_argument("--slide", type=float, default=600.0, help="time-slide step, seconds")
    ap.add_argument("--slides", type=int, default=50)
    ap.add_argument("--csv", help="write the per-alarm table here")
    a = ap.parse_args()

    paths = sorted({p for g in a.records for p in glob.glob(g)})
    if not paths:
        sys.exit("no record files matched")

    rows, tot_hours, obs_total, null_means, null_sds = [], 0.0, 0, [], []
    per_station_rate = {}
    for path in paths:
        stations, det, span, picks = read_record(path)
        if not span:
            print(f"{Path(path).name}: no T line, skipped")
            continue
        hours = (span[1] - span[0]) / 3600
        tot_hours += hours
        det = {s: cluster(ts, a.coda) for s, ts in det.items()}
        # Every station the record declares (S lines) contributes its hours to
        # the denominator, whether or not it detected anything: crediting only
        # the stations that fired would inflate the rate of the quiet ones.
        for s in stations:
            r = per_station_rate.setdefault(s, [0, 0.0])
            r[0] += len(det.get(s, ()))
            r[1] += hours

        obs = coincidences(det, stations, a.vp, a.slack)
        obs_total += len(obs)
        nulls = []
        codes = sorted(det)
        # Fewer than two stations with data cannot produce a coincidence, and
        # the slide has nothing to shift against.
        for k in range(1, a.slides + 1) if len(codes) > 1 else ():
            sh = {c: (i + 1) * k * a.slide for i, c in enumerate(codes)}
            sh[codes[0]] = 0.0
            nulls.append(len(coincidences(det, stations, a.vp, a.slack, sh, span)))
        mu = statistics.fmean(nulls) if nulls else 0.0
        sd = statistics.pstdev(nulls) if len(nulls) > 1 else 0.0
        null_means.append(mu)
        null_sds.append(sd)
        print(f"{Path(path).name}: {hours:5.2f} h, "
              f"{sum(len(v) for v in det.values()):>4} station detections, "
              f"{len(obs):>3} coincident pairs, null {mu:5.2f} +- {sd:4.2f}")

        for t, s1, s2, dt_, w in obs:
            # Chance that the partner fell inside this pair's window.
            n2 = len(det[s2])
            lam2 = n2 / (hours * 3600) if hours else 0.0
            p_chance = 1.0 - math.exp(-lam2 * 2 * w)
            verdict, detail = locate_check(picks, stations, s1, s2, t, dt_)
            rows.append({
                "record": Path(path).name, "time_utc": f"{t:.2f}",
                "station_a": s1, "station_b": s2,
                "dt_s": round(dt_, 2), "window_s": round(w, 2),
                "tightness": round(abs(dt_) / w, 3) if w else "",
                "km_apart": round(hav(*stations[s1], *stations[s2]), 1),
                "p_chance": round(p_chance, 4),
                "verdict": verdict, "detail": detail,
            })

    print(f"\n{len(paths)} record(s), {tot_hours:.1f} h")
    print(f"coincident pairs observed : {obs_total}")
    print(f"expected by chance (slides): {sum(null_means):.2f} "
          f"+- {math.sqrt(sum(s * s for s in null_sds)):.2f}")
    excess = obs_total - sum(null_means)
    print(f"chance-corrected excess    : {excess:+.2f}"
          f"  ({'consistent with noise' if excess <= math.sqrt(sum(s*s for s in null_sds) + 1e-9) * 2 else 'more than noise'})")
    if rows:
        import collections as _c
        tally = _c.Counter(r["verdict"] for r in rows)
        print(f"by the S-P location test   : " +
              ", ".join(f"{k} {v}" for k, v in sorted(tally.items())))
        print("  'earthquake' means a single source fits both stations' S-P distances")
        print("  and predicts their arrival difference; it is strong evidence, not a")
        print("  substitute for looking at the waveform.")
    print("\nper-station detection rate (clustered):")
    for s, (n, h) in sorted(per_station_rate.items()):
        print(f"  {s:<6} {n:>5} in {h:6.1f} h = {n / h * 24:6.2f}/day" if h else f"  {s}: -")
    if a.csv and rows:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"\n{len(rows)} rows -> {a.csv}")


if __name__ == "__main__":
    main()
