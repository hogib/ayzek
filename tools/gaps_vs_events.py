"""Tests whether real-time data gaps at DEMI coincide with earthquakes.

Real-time gaps are recomputed from mseed_dump's index.csv with the reorderer's
max_lateness = 0 rule (the counts match tools/replay_check on the DEMI chunk).

  1. catalogue  fraction of events whose P-to-S window at DEMI overlaps a gap,
                compared with a null distribution from random circular time
                shifts of the catalogue (shifts preserve the clustering of
                aftershocks and of gaps)
  2. station    DEMI vertical-component amplitude at the midpoint of each
                horizontal gap, as a percentile within its UTC hour, compared
                with circularly shifted gap times
  alignment     amplitude ratio after the predicted P time for M >= 3.5 events,
                to confirm that catalogue times are UTC

The gap rate changes during the chunk. Tests 1 and 2 are therefore run
separately before and after a change point, estimated by maximising a two-rate
Poisson likelihood over hourly gap counts, and over the whole chunk.

    uv run --with numpy --with pandas python tools/gaps_vs_events.py DUMP_DIR CATALOG.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DEMI_LON, DEMI_LAT = 28.7162, 39.0428
FS = 100
BIN = 10 * FS                       # amplitude bins, in samples
BIN_S = BIN / FS
VP, VS = 6.0, 3.5                   # km/s, uniform; windows are padded for the error
PAD_BEFORE, PAD_AFTER = 5.0, 20.0   # seconds before P, after S
N_SHIFTS = 10_000
RNG = np.random.default_rng(20251028)

# Event selections, defined before running the tests. The first is the primary one.
TIERS = [
    ("M>=3.0 within 150 km *", 3.0, 150.0),
    ("M>=2.0 within 100 km", 2.0, 100.0),
    ("M>=4.0 within 300 km", 4.0, 300.0),
]


def realtime_gaps(pos, n):
    """The gaps Reorderer emits at max_lateness = 0, as [from, to) in samples."""
    gaps, nxt = [], None
    for p, k in zip(pos.tolist(), n.tolist()):
        end = p + k
        if nxt is None:
            nxt = end
        elif end > nxt:                 # not stale
            if p > nxt:
                gaps.append((nxt, p))
            nxt = end
    return gaps


def merge(intervals):
    out = []
    for a, b in sorted(intervals):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return np.array(out, dtype=np.float64).reshape(-1, 2)


def touches(gaps, ws, we):
    """Whether each window [ws, we) overlaps any gap."""
    starts, ends = gaps[:, 0], gaps[:, 1]
    i = np.searchsorted(ends, ws, side="right")                 # first gap ending after ws
    return (i < len(starts)) & (starts[np.minimum(i, len(starts) - 1)] < we)


def covered(gaps, lo, hi):
    """Fraction of [lo, hi) inside a gap."""
    clipped = np.clip(gaps, lo, hi)
    return (clipped[:, 1] - clipped[:, 0]).sum() / (hi - lo)


def pvalue(obs, null):
    return (1 + np.sum(null >= obs)) / (1 + len(null))


def stamp(t):
    return pd.Timestamp(t, unit="s").strftime("%Y-%m-%d %H:%M")


# --- real-time gaps ------------------------------------------------------------

dump, catalog = Path(sys.argv[1]), Path(sys.argv[2])
idx = pd.read_csv(dump / "index.csv")
idx["pos"] = idx.start_ns // 10_000_000

gaps = {}
for cha, g in idx.groupby("channel", sort=True):                # groupby keeps file order
    gaps[cha] = realtime_gaps(g.pos.to_numpy(), g.num_samples.to_numpy())

t0 = idx.pos.min() / FS
t1 = (idx.pos + idx.num_samples).max() / FS
span = t1 - t0
lost = merge([(a / FS, b / FS) for g in gaps.values() for a, b in g])
horiz = merge([(a / FS, b / FS) for c in ("HHE", "HHN") for a, b in gaps[c]])
# DEMI chunk: 485, 513 and 25 gaps, the same as tools/replay_check.
print(f"chunk {stamp(t0)} .. {stamp(t1)} UTC, {span / 86400:.1f} days; real-time gaps "
      + ", ".join(f"{c} {len(g)}" for c, g in gaps.items()))

# One change in the rate of horizontal gaps, by maximum Poisson likelihood over hours.
hourly = np.bincount(((horiz[:, 0] - t0) // 3600).astype(int), minlength=int(span // 3600))


def loglik(c):
    lam = c.mean()
    return 0.0 if lam == 0 else float((c * np.log(lam)).sum() - lam * len(c))


split = max(range(24, len(hourly) - 24), key=lambda k: loglik(hourly[:k]) + loglik(hourly[k:]))
tc = t0 + split * 3600
PERIODS = [("before", t0, tc), ("after", tc, t1), ("whole chunk", t0, t1)]
print(f"gap rate changes at {stamp(tc)} UTC: "
      f"{24 * hourly[:split].mean():.1f}/day before, {24 * hourly[split:].mean():.1f}/day after")
for name, lo, hi in PERIODS:
    print(f"   time inside a gap on some channel, {name:11s} {100 * covered(lost, lo, hi):.3f}%")

# --- catalogue -------------------------------------------------------------------

cat = pd.read_csv(catalog, encoding="utf-8-sig")
ts = pd.to_datetime(cat.Date, format="%d/%m/%Y %H:%M:%S", utc=True)
cat["t"] = (ts - pd.Timestamp(0, tz="UTC")).dt.total_seconds()
cat = cat[(cat.t >= t0) & (cat.t < t1 - 120)].copy()
la1, lo1 = np.radians(cat.Latitude), np.radians(cat.Longitude)
la2, lo2 = np.radians(DEMI_LAT), np.radians(DEMI_LON)
h = np.sin((la1 - la2) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo1 - lo2) / 2) ** 2
cat["R"] = np.hypot(2 * 6371 * np.arcsin(np.sqrt(h)), cat.Depth)
near = cat[cat.R <= 150]
top = near.loc[near.Magnitude.idxmax()]
print(f"\ncatalogue: {len(cat)} events in the chunk, {len(near)} within 150 km; "
      f"largest within 150 km M{top.Magnitude} at {top.R:.0f} km")


def windows(t, R):
    return t + R / VP - PAD_BEFORE, t + R / VS + PAD_AFTER


# --- vertical amplitude, 10 s bins ------------------------------------------------

z = np.memmap(dump / "HHZ.i32", dtype=np.int32, mode="r")
zi = idx[idx.channel == "HHZ"]
zpos, zn = zi.pos.to_numpy(), zi.num_samples.to_numpy()
zoff = np.concatenate([[0], np.cumsum(zn)[:-1]])
brk = np.flatnonzero(zpos[1:] != zpos[:-1] + zn[:-1]) + 1
b0 = int(zpos.min() // BIN)
nb = int((zpos[-1] + zn[-1]) // BIN) - b0 + 1
sumsq, cnt = np.zeros(nb), np.zeros(nb)
BLOCK = 86400 * FS
for s, e in zip(np.concatenate([[0], brk]), np.concatenate([brk, [len(zpos)]])):
    f_begin, f_end, p_begin = zoff[s], zoff[e - 1] + zn[e - 1], zpos[s]
    for f in range(int(f_begin), int(f_end), BLOCK):
        x = np.diff(z[f:min(f + BLOCK + 1, f_end)].astype(np.float64))   # first difference as a high-pass
        bins = (p_begin + (f - f_begin) + 1 + np.arange(len(x))) // BIN - b0
        sumsq += np.bincount(bins, weights=x * x, minlength=nb)
        cnt += np.bincount(bins, minlength=nb)

amp = pd.DataFrame({"amp": np.where(cnt >= 0.8 * BIN, np.sqrt(sumsq / np.maximum(cnt, 1)), np.nan)})
amp["t"] = (b0 + np.arange(nb)) * BIN_S
amp["hour"] = (amp.t // 3600).astype(int) % 24
amp_t, amp_v = amp.t.to_numpy(), amp.amp.to_numpy()

big = cat[(cat.Magnitude >= 3.5) & (cat.R <= 150)]
print(f"\nalignment: amplitude just after predicted P over the minute before, "
      f"median of {len(big)} events M>=3.5 within 150 km")
for label, shift in [("catalogue as given (UTC)", 0), ("shifted -3 h", -10800), ("shifted +3 h", 10800)]:
    kb = ((big.t.to_numpy() + big.R.to_numpy() / VP + shift - amp_t[0]) // BIN_S).astype(int)
    kb = kb[(kb >= 8) & (kb + 1 < nb)]
    ratio = amp_v[kb + 1] / np.array([np.nanmedian(amp_v[i - 8:i - 2]) for i in kb])
    print(f"   {label:26s} {np.nanmedian(ratio):6.1f}x")

# --- the tests --------------------------------------------------------------------


def catalogue_test(lo, hi):
    L = hi - lo
    for name, mmin, rmax in TIERS:
        ev = cat[(cat.Magnitude >= mmin) & (cat.R <= rmax) & (cat.t >= lo) & (cat.t < hi - 120)]
        t, R = ev.t.to_numpy(), ev.R.to_numpy()
        if len(t) < 5:
            print(f"     {name:22s} {len(t):5d} events -- too few")
            continue
        obs_hits = int(touches(lost, *windows(t, R)).sum())
        obs = obs_hits / len(t)
        null = np.array([touches(lost, *windows(lo + (t - lo + s) % L, R)).mean()
                         for s in RNG.uniform(0, L, N_SHIFTS)])
        # Shifts by whole days preserve time of day, controlling for a daily cycle.
        days = np.array([touches(lost, *windows(lo + (t - lo + d * 86400) % L, R)).mean()
                         for d in range(1, int(L // 86400))])
        ratio = f"{obs / null.mean():4.2f}" if null.mean() > 0 else " n/a"
        q95 = np.quantile(null, 0.95)
        detect = f"{q95 / null.mean():4.1f}x" if q95 > 0 else "  n/a"      # undefined with too few events
        print(f"     {name:22s} {len(t):5d} events {obs_hits:3d} in a gap   "
              f"{100 * obs:5.2f}% vs chance {100 * null.mean():5.2f}%   ratio {ratio}   "
              f"p = {pvalue(obs, null):.3f}   day shifts beaten {int((obs > days).sum()):2d}/{len(days):<2d}   "
              f"could detect >= {detect}")


def amplitude_test(lo, hi):
    sel = np.flatnonzero((amp_t >= lo) & (amp_t < hi))
    pct = amp.iloc[sel].groupby("hour").amp.rank(pct=True).to_numpy()
    n = len(sel)
    # One bin per gap, at its midpoint. Using all bins inside gaps would weight
    # gaps by length, and gap length depends on amplitude (Steim2 records hold
    # more samples when the signal is quiet; see the correlation printed below).
    mids = (horiz[:, 0] + horiz[:, 1]) / 2
    mb = ((mids[(mids >= lo) & (mids < hi)] - amp_t[sel[0]]) // BIN_S).astype(int)
    mb = mb[(mb >= 0) & (mb < n)]
    mb = mb[~np.isnan(pct[mb])]
    shifts = RNG.integers(1, n, N_SHIFTS)
    for label, stat in [("mean amplitude rank", np.nanmean),
                        ("share in top 5% for hour", lambda v: np.nanmean(v > 0.95)),
                        ("share in top 1% for hour", lambda v: np.nanmean(v > 0.99))]:
        obs = stat(pct[mb])
        null = np.array([stat(pct[(mb + s) % n]) for s in shifts])
        print(f"     {label:25s} {obs:5.3f} vs chance {null.mean():5.3f}   "
              f"louder p = {pvalue(obs, null):.3f}   quieter p = {pvalue(-obs, -null):.3f}")
    return len(mb)


print("\n1. do event windows at DEMI (P - 5 s to S + 20 s) overlap real-time gaps more than chance?")
print(f"   null: {N_SHIFTS:,} circular time shifts within the period.  * = primary\n")
for name, lo, hi in PERIODS:
    print(f"   {name} ({stamp(lo)} .. {stamp(hi)})")
    catalogue_test(lo, hi)

print("\n2. is DEMI's vertical louder than usual for its hour when the horizontals drop out?")
print(f"   null: {N_SHIFTS:,} circular shifts of the gap midpoints within the period\n")
for name, lo, hi in PERIODS:
    print(f"   {name}")
    n_gaps = amplitude_test(lo, hi)
    print(f"     ({n_gaps} gaps)")

# --- context --------------------------------------------------------------------------

zrank = amp.groupby("hour").amp.rank(pct=True).to_numpy()
rec_rank = pd.Series(zrank[np.clip(zpos // BIN - b0, 0, nb - 1)])
print(f"\nsamples per HHZ record vs amplitude rank: Spearman "
      f"{pd.Series(zn).rank().corr(rec_rank.rank()):.2f}  (quiet signal, longer records)")

print("\nhorizontal gap starts by UTC hour, after the change (local = UTC+3):")
after = horiz[horiz[:, 0] >= tc, 0]
print("   " + " ".join(f"{h:>3d}" for h in range(24)))
print("   " + " ".join(f"{c:>3d}" for c in np.bincount((after // 3600).astype(int) % 24, minlength=24)))

n_days = int(round(span / 86400))
day = lambda t: ((np.asarray(t) - t0) // 86400).astype(int)
daily = pd.DataFrame({
    "gaps": np.bincount(day(horiz[:, 0]), minlength=n_days)[:n_days],
    "events": np.bincount(day(cat[(cat.Magnitude >= 2.0) & (cat.R <= 100)].t), minlength=n_days)[:n_days],
})
print(f"\nper day from {stamp(t0)[:10]}: horizontal gaps, and M>=2 events within 100 km")
print("   gaps   " + " ".join(f"{v:>3d}" for v in daily.gaps))
print("   events " + " ".join(f"{v:>3d}" for v in daily.events))
