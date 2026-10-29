"""Uncertainty of the detector comparison: paired tests on events, and on the quiet windows.

    python3 tools/significance.py --outcomes data/runs/stats \\
        --quiet data/quiet100/quiet_index.csv \\
        --quiet-runs learned=data/runs/quiet100 stalta=data/runs/quiet100_stalta \\
        --json significance.json

Events (tools/event_outcomes.py writes the inputs): for each pair of detectors
on the same replays,

- an exact McNemar test on the events one detected and the other missed;
- a paired bootstrap over events of the difference in chance-corrected
  detections. Each resample draws the observed events, and separately the
  events of each shifted catalogue, with the same draws for both detectors, and
  recomputes o, c and (o - c) / (1 - c) as tools/sweep_trigger.py does;
- Clopper-Pearson intervals of each detector's observed detection rate.

Quiet windows: the exact conditional test of two Poisson counts over the same
hours, with the exact interval of their ratio, and a bootstrap over windows,
which keeps the clustering of alarms within a window (one station raises most).

--equivalence OTHER=FRACTION adds, for that comparison, an equivalence test
by two one-sided tests at 5%: the learned detector and OTHER are shown
equivalent when the 90% paired bootstrap interval of the difference in
chance-corrected detections lies within +-FRACTION of the group's events.

--block-minutes resamples blocks of origin time rather than single events,
which keeps together events whose detectability is shared (noise state, the
stations up, overlapping codas). --quiet-runs compares every run after the first
with the first. --json writes
the figures under "summary" for tools/reproduce_report.py: "GROUP: OTHER
corrected" and so on for the events, "quiet: LABEL per day" and "quiet: OTHER
ratio" for the quiet windows.

Standard library only; the seed is fixed so every number reproduces.
"""
import argparse
import collections
import csv
import glob
import json
import math
import random
import re
from pathlib import Path

B = 10000


def binom_cdf(k, n, p=0.5):
    return sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(k + 1))


def mcnemar(b, c):
    n = b + c
    return 1.0 if n == 0 else min(1.0, 2 * binom_cdf(min(b, c), n))


def beta_ppf(q, a, b):
    """Quantile of Beta(a, b) by bisection on the regularized incomplete beta."""
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if betainc(a, b, mid) < q:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def betainc(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    # Continued fraction (Numerical Recipes betacf), with the symmetry switch.
    if x > (a + 1) / (a + b + 2):
        return 1 - betainc(b, a, 1 - x)
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x)
    f, c, d = 1.0, 1.0, 0.0
    for i in range(400):
        m = i // 2
        if i == 0:
            num = 1.0
        elif i % 2 == 0:
            num = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
        else:
            num = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1 + num * d
        d = 1 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1 + num / c if abs(c) > 1e-300 else 1e-300
        f *= c * d
        if abs(1 - c * d) < 1e-12:
            break
    return math.exp(lbeta) * (f - 1) / a


def clopper_pearson(k, n, alpha=0.05):
    lo = 0.0 if k == 0 else beta_ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta_ppf(1 - alpha / 2, k + 1, n - k)
    return lo, hi


def pct(v, q):
    v = sorted(v)
    return v[min(len(v) - 1, max(0, int(round(q * (len(v) - 1)))))]


def load_outcomes(d):
    """{(dataset, detector): {shift: {event: detected}}}"""
    out = collections.defaultdict(lambda: collections.defaultdict(dict))
    for f in glob.glob(f"{d}/*.csv"):
        for r in csv.DictReader(open(f)):
            ev = (r["dataset"], r["origin"], r["type"], r["magnitude"])
            out[r["detector"]][int(r["shift"])][ev] = int(r["detected"])
    return out


def corrected(obs, shifted, take_obs, take_shift):
    """Chance-corrected detections for one detector and one resample: per
    replay, then summed, as tools/sweep_trigger.py reports them."""
    total = 0.0
    for ds in {e[0] for e in take_obs}:
        to = [e for e in take_obs if e[0] == ds]
        o = sum(obs[e] for e in to) / len(to)
        cs = []
        for s in shifted:
            ts = [e for e in take_shift[s] if e[0] == ds]
            if ts:
                cs.append(sum(shifted[s][e] for e in ts) / len(ts))
        c = sum(cs) / len(cs)
        total += (o - c) / (1 - c) * len(to)
    return total, sum(obs[e] for e in take_obs)


def by_ds(events):
    g = collections.defaultdict(list)
    for e in events:
        g[e[0]].append(e)
    return g


def seconds(origin):
    h, m, sec = origin.split(":")
    return int(h) * 3600 + int(m) * 60 + float(sec)


def resample(events, rng, block):
    """One bootstrap draw within each replay: events one by one, or with
    block > 0, whole blocks of `block` seconds of origin time, as many blocks as
    the replay has, so that events close in time are drawn together."""
    out = []
    for g in by_ds(events).values():
        if not block:
            out += [rng.choice(g) for _ in g]
            continue
        blocks = collections.defaultdict(list)
        for e in g:
            blocks[int(seconds(e[1]) // block)].append(e)
        bs = list(blocks.values())
        for _ in bs:
            out += rng.choice(bs)
    return out


def compare(det, a, b, datasets, rng, block=0):
    def subset(d):
        return {s: {e: v for e, v in ev.items() if e[0] in datasets} for s, ev in det[d].items()}
    A, Bd = subset(a), subset(b)
    obs_ev = sorted(A[0])
    shifts = sorted(s for s in A if s != 0)
    shift_ev = {s: sorted(A[s]) for s in shifts}
    shA = {s: A[s] for s in shifts}
    shB = {s: Bd[s] for s in shifts}
    ca, _ = corrected(A[0], shA, obs_ev, shift_ev)
    cb, _ = corrected(Bd[0], shB, obs_ev, shift_ev)
    nb = sum(A[0][e] and not Bd[0][e] for e in obs_ev)
    nc = sum(Bd[0][e] and not A[0][e] for e in obs_ev)
    diffs = []
    for _ in range(B):
        # Resample within each replay, so every replay keeps its size.
        to = resample(obs_ev, rng, block)
        ts = {s: resample(shift_ev[s], rng, block) for s in shifts}
        diffs.append(corrected(A[0], shA, to, ts)[0] - corrected(Bd[0], shB, to, ts)[0])
    ka, kb, n = sum(A[0].values()), sum(Bd[0].values()), len(obs_ev)
    return {
        "events": n, "a_detected": ka, "b_detected": kb,
        "a_corrected": ca, "b_corrected": cb, "difference": ca - cb,
        "diff_ci": (pct(diffs, 0.025), pct(diffs, 0.975)),
        "diff_ci90": (pct(diffs, 0.05), pct(diffs, 0.95)),
        "p_boot_le0": sum(d <= 0 for d in diffs) / B,
        "a_only": nb, "b_only": nc, "mcnemar_p": mcnemar(nb, nc),
        "a_cp": clopper_pearson(ka, n), "b_cp": clopper_pearson(kb, n),
    }


def quiet_counts(index, runs):
    hours = {r["name"]: float(r["scored_hours"]) for r in csv.DictReader(open(index))}
    counts = {}
    for label, d in runs.items():
        counts[label] = {}
        for name in hours:
            m = re.search(r"Alarms raised\s+(\d+)", Path(f"{d}/{name}.txt").read_text())
            counts[label][name] = int(m.group(1)) if m else 0
    return hours, counts


def quiet(hours, counts, rng):
    (la, a), (lb, b) = counts.items()
    T = sum(hours.values())
    xa, xb = sum(a.values()), sum(b.values())
    # Conditional on the total, xa ~ Binomial(xa + xb, 1/2) under equal rates.
    p = mcnemar(xa, xb)
    lo, hi = clopper_pearson(xa, xa + xb)
    ratio_ci = (lo / (1 - lo), hi / (1 - hi))
    names = sorted(hours)
    boot = []
    for _ in range(B):
        take = [rng.choice(names) for _ in names]
        h = sum(hours[n] for n in take)
        ra, rb = 24 * sum(a[n] for n in take) / h, 24 * sum(b[n] for n in take) / h
        boot.append((ra, rb, ra / rb if rb else float("inf")))
    return {
        "hours": T, la: xa, lb: xb,
        f"{la}_per_day": 24 * xa / T, f"{lb}_per_day": 24 * xb / T,
        "ratio": xa / xb, "ratio_exact_ci": ratio_ci, "p_equal_rates": p,
        f"{la}_per_day_boot_ci": (pct([x[0] for x in boot], 0.025), pct([x[0] for x in boot], 0.975)),
        f"{lb}_per_day_boot_ci": (pct([x[1] for x in boot], 0.025), pct([x[1] for x in boot], 0.975)),
        "ratio_boot_ci": (pct([x[2] for x in boot], 0.025), pct([x[2] for x in boot], 0.975)),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--outcomes", required=True)
    ap.add_argument("--quiet")
    ap.add_argument("--quiet-runs", nargs="+", metavar="LABEL=DIR",
                    help="the reference first, then each run compared with it")
    ap.add_argument("--against", default="stalta_eqfa,stalta_eqdet",
                    help="detectors compared with `learned`, comma-separated labels")
    ap.add_argument("--groups", default="held out=sindirgi61;tuning=demo,marmara;all three=sindirgi61,demo,marmara",
                    help="NAME=DATASET,...;... replays pooled per comparison")
    ap.add_argument("--seed", type=int, default=20261022)
    ap.add_argument("--json", help="also write the figures here")
    ap.add_argument("--equivalence", nargs="*", default=[], metavar="OTHER=FRACTION",
                    help="equivalence margin, as a fraction of the group's events, "
                         "for the comparison with OTHER")
    ap.add_argument("--block-minutes", type=float, default=0,
                    help="bootstrap whole blocks of this many minutes of origin time "
                         "instead of single events (0: events)")
    a = ap.parse_args()
    rng = random.Random(a.seed)

    det = load_outcomes(a.outcomes)
    margins = {k: float(v) for k, v in (x.split("=", 1) for x in a.equivalence)}
    summary = {}
    groups = {g.split("=")[0]: set(g.split("=")[1].split(",")) for g in a.groups.split(";")}
    for label, ds in groups.items():
        print(f"\n== {label}")
        for other in a.against.split(","):
            if not any(e[0] in ds for e in det[other][0]):
                print(f"learned vs {other}: no outcomes for {sorted(ds)}")
                continue
            r = compare(det, "learned", other, ds, rng, 60 * a.block_minutes)
            k = f"{label}: {other}"
            summary.update({f"{label}: learned corrected": round(r["a_corrected"], 1),
                            f"{k} detected": r["b_detected"],
                            f"{k} corrected": round(r["b_corrected"], 1),
                            f"{k} difference": round(r["difference"], 1),
                            f"{k} difference low": round(r["diff_ci"][0], 1),
                            f"{k} difference high": round(r["diff_ci"][1], 1),
                            f"{k} mcnemar p": float(f"{r['mcnemar_p']:.2g}")})
            print(f"learned vs {other}: {r['events']} events")
            print(f"  detected            {r['a_detected']} vs {r['b_detected']}  "
                  f"(95% CP {r['a_cp'][0]:.2f}-{r['a_cp'][1]:.2f} vs {r['b_cp'][0]:.2f}-{r['b_cp'][1]:.2f})")
            print(f"  chance-corrected    {r['a_corrected']:.1f} vs {r['b_corrected']:.1f}, "
                  f"difference {r['difference']:+.1f} (95% bootstrap {r['diff_ci'][0]:+.1f} to {r['diff_ci'][1]:+.1f}, "
                  f"P(diff <= 0) = {r['p_boot_le0']:.4f})")
            print(f"  discordant          learned only {r['a_only']}, {other} only {r['b_only']}, "
                  f"exact McNemar p = {r['mcnemar_p']:.2g}")
            if other in margins:
                m = margins[other] * r["events"]
                lo, hi = r["diff_ci90"]
                eq = -m < lo and hi < m
                summary.update({f"{k} equivalence margin": round(m, 1),
                                f"{k} difference 90% low": round(lo, 1),
                                f"{k} difference 90% high": round(hi, 1),
                                f"{k} equivalent": eq})
                print(f"  equivalence         margin +-{m:.1f} events ({margins[other]:.0%} of "
                      f"{r['events']}); 90% interval {lo:+.1f} to {hi:+.1f}: "
                      + ("equivalent (both one-sided tests p < 0.05)" if eq else "not shown"))

    if a.quiet and a.quiet_runs:
        runs = dict(x.split("=", 1) for x in a.quiet_runs)
        hours, counts = quiet_counts(a.quiet, runs)
        ref, *others = counts
        for other in others:
            q = quiet(hours, {ref: counts[ref], other: counts[other]}, rng)
            print("\n== quiet windows")
            for k, v in q.items():
                print(f"  {k:28s} " + (f"{v[0]:.3g} to {v[1]:.3g}" if isinstance(v, tuple) else f"{v:.4g}"))
            for lab in (ref, other):
                summary[f"quiet: {lab} alarms"] = q[lab]
                summary[f"quiet: {lab} per day"] = round(q[f"{lab}_per_day"], 2)
            summary["quiet: hours"] = round(q["hours"], 1)
            summary[f"quiet: {other} ratio"] = round(q["ratio"], 3)
            summary[f"quiet: {other} p"] = float(f"{q['p_equal_rates']:.2g}")

    if a.json:
        Path(a.json).write_text(json.dumps({"summary": summary}, indent=1, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
