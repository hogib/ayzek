"""6 s detector vs the streaming transformer, on the same continuous replays.

    uv run --project tools python tools/compare_detectors.py --only demo_ko,marmara_ko,quiet_ko
    uv run --project tools python tools/compare_detectors.py --json compare.json      # every dataset

Runs ayzek twice per dataset, `--detector 6s` and `--detector transformer`, with
everything else identical (picker, magnitude, association, the datasets and
catalogue excerpts of tools/scorecard.py), and reports two levels:

**Station level: the detector itself.** For every visible catalogued arrival at
every station (the visibility rule of onset's catalog.VISIBILITY), P is predicted
with TauP (iasp91) and refined with an AIC pick on that station's vertical
(accepted when the second after it has twice the RMS of the second before; else
the TauP time is kept and flagged). Against that P:

  latency     decision time (a Detection's declared_at) minus P, for the first
              detection at the station in [P - 0.5 s, P + 20 s]
  recall@d    arrivals detected with latency <= d
  P dating    the detection's own P estimate (window start + 3.5 s: the STA/LTA
              anchor or 3.5 s convention for 6 s, t - dt for the transformer)
              minus P; this is what places the picker window
  unmatched   station triggers with no visible arrival in [P - 0.5, P + 60 s],
              per station-hour; not necessarily false, the catalogue is incomplete

**Network level: the whole pipeline** (tools/scorecard.py's parser): catalogue
events detected, unmatched alarms, median alarm delay after origin, false
alarms per day on the quiet sets, magnitude error.

Both detectors must see identical input, so this compares what each detector
does with a real stream, gaps and all, rather than their test-set numbers:
the transformer's test set is FDSN KO windows within 55 km, and several of these
replays are AFAD stations or events farther out. The per-distance table shows
where each one's recall comes from.
"""
import argparse
import concurrent.futures
import json
import re
import shlex
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from obspy import UTCDateTime, read
from obspy.taup import TauPyModel

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from scorecard import DATASETS, files_of, score_dataset  # noqa: E402

DETECTORS = ("6s", "transformer")
VISIBILITY = ((50.0, 0.0), (150.0, 2.0), (400.0, 3.0), (1500.0, 4.5))
DELAYS = (0.5, 1.0, 2.0, 4.0, 8.0)
DIST_BINS = ((0, 30), (30, 60), (60, 100), (100, 200), (200, 1500))
FS = 100.0


# --- running ayzek ---------------------------------------------------------------

def run(name, detector, build, models, extra, workdir, reuse=False):
    files, cat, start, _end, _role = DATASETS[name]
    paths = files_of(files)
    if not paths:
        return None
    rec = Path(workdir) / f"{name}_{detector}.rec"
    out = Path(workdir) / f"{name}_{detector}.txt"
    if reuse and rec.exists() and out.exists():
        return {"stdout": out.read_text(), "record": rec.read_text(), "files": paths}
    cmd = [str(ROOT / build / "app/ayzek"), "--speed", "0", "--no-color",
           "--models", str(models), "--detector", detector, "--record", str(rec),
           "--catalog", str(ROOT / "tests/catalogs" / f"{cat}.csv")]
    if start:
        cmd += ["--from", start]
    cmd += extra + paths
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(f"{name} {detector}: ayzek exited {p.returncode}: {p.stderr[-400:]}")
    out.write_text(p.stdout)
    return {"stdout": p.stdout, "record": rec.read_text(), "files": paths}


def parse_record(text):
    """Station coordinates, the replay span and every detection."""
    stations, dets, span = {}, defaultdict(list), None
    for line in text.splitlines():
        f = line.split("\t")
        if f[0] == "S":
            stations[f[1]] = (float(f[2]), float(f[3]))
        elif f[0] == "T":
            span = (float(f[1]), float(f[2]))
        elif f[0] == "D":
            dets[f[1]].append({"window_start": float(f[2]), "declared_at": float(f[3]),
                               "p": float(f[4]), "ms": float(f[5])})
    return stations, dets, span


def cpu_per_hour(stdout):
    """Detector milliseconds per station-hour, from the station table: windows
    (or tokens) times mean ms per window, over the replay's stream time."""
    rows = re.findall(r"^  (\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\S+\s+([\d.]+) \(", stdout, re.M)
    total_ms = sum(int(r[1]) * float(r[5]) for r in rows)
    return total_ms, len(rows)


# --- ground truth ---------------------------------------------------------------

def load_catalog(path):
    out = []
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines()[1:]:
        f = line.split(",")
        try:
            t = UTCDateTime.strptime(f[0], "%d/%m/%Y %H:%M:%S").timestamp
            out.append({"t": t, "lon": float(f[1]), "lat": float(f[2]),
                        "depth": float(f[3]), "mag": float(f[6]), "id": f[-1]})
        except (ValueError, IndexError):
            continue
    return out


def km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return float(2 * 6371.0 * np.arcsin(np.sqrt(a)))


def visible(mag, dist):
    return any(dist <= d and mag >= m for d, m in VISIBILITY)


def aic_pick(z, i_guess, i_s):
    """onset.labels.refine_p: AIC on [P - 2 s, P + 3 s], stopping before S."""
    lo, hi = max(int(1.5 * FS), int(i_guess - 2 * FS)), min(len(z), int(i_guess + 3 * FS))
    if i_s is not None:
        hi = min(hi, int(i_s - 0.5 * FS))
    if hi - lo < int(1.5 * FS) or not np.all(np.isfinite(z[lo:hi])):
        return None
    x = z[lo:hi] - z[lo:hi].mean()
    n = len(x)
    c1, c2 = np.cumsum(x), np.cumsum(x * x)
    k = np.arange(2, n - 2)
    va = c2[k - 1] / k - (c1[k - 1] / k) ** 2
    vb = (c2[-1] - c2[k - 1]) / (n - k) - ((c1[-1] - c1[k - 1]) / (n - k)) ** 2
    eps = 1e-12 * max(c2[-1] / n, 1e-30)
    pick = lo + int(k[np.argmin(k * np.log(va + eps) + (n - k - 1) * np.log(vb + eps))])
    pre, post = z[max(0, pick - 100):pick], z[pick:pick + 100]
    if len(pre) < 50 or len(post) < 50:
        return None
    snr = np.sqrt(np.mean(post ** 2) / max(np.mean(pre ** 2), 1e-30))
    return pick if snr >= 2.0 else None


class Truth:
    """Visible catalogued P arrivals per station, AIC-refined where possible."""

    def __init__(self):
        self.taup = TauPyModel("iasp91")
        self.cache = {}

    def travel(self, dist, depth, phases):
        key = (round(dist), round(max(depth, 0)), phases)
        if key not in self.cache:
            a = self.taup.get_travel_times(source_depth_in_km=key[1],
                                           distance_in_degree=key[0] / 111.195,
                                           phase_list=list(phases))
            self.cache[key] = a[0].time if a else None
        return self.cache[key]

    def arrivals(self, files, stations, catalog, span):
        out = defaultdict(list)
        z_by_station = {}
        for f in files:
            st = read(f)
            z = st.select(component="Z")
            if not len(z):
                continue
            code = z[0].stats.station
            for tr in z:
                tr.data = tr.data.astype(np.float64)
            z.merge(method=0, fill_value=np.nan)
            tr = z[0]
            ok = np.isfinite(tr.data)
            tr.data = np.where(ok, tr.data - np.nanmean(tr.data), 0.0)
            tr.filter("bandpass", freqmin=1.0, freqmax=45.0, corners=4, zerophase=True)
            d = tr.data.copy()
            d[~ok] = np.nan
            z_by_station[code] = (tr.stats.starttime.timestamp, d)
        for code, (lat, lon) in stations.items():
            if code not in z_by_station:
                continue
            t0, z = z_by_station[code]
            for ev in catalog:
                dist = km(ev["lat"], ev["lon"], lat, lon)
                if not visible(ev["mag"], dist):
                    continue
                tp = self.travel(dist, ev["depth"], ("p", "P", "Pn", "Pg"))
                if tp is None:
                    continue
                t_p = ev["t"] + tp
                if not (span[0] + 30 <= t_p <= span[1] - 20):
                    continue
                ts = self.travel(dist, ev["depth"], ("s", "S", "Sn", "Sg"))
                i = (t_p - t0) * FS
                i_s = (ev["t"] + ts - t0) * FS if ts else None
                pick = aic_pick(z, i, i_s) if 0 <= i < len(z) else None
                out[code].append({"p": t0 + pick / FS if pick is not None else t_p,
                                  "refined": pick is not None, "dist": dist,
                                  "mag": ev["mag"], "event": ev["id"]})
        return out


# --- scoring --------------------------------------------------------------------

def station_level(arrivals, dets, hours):
    lat, dating, dist, recs = [], [], [], []
    unmatched = 0
    for code, arrs in arrivals.items():
        ds = sorted(dets.get(code, []), key=lambda d: d["declared_at"])
        for a in arrs:
            hit = next((d for d in ds if a["p"] - 0.5 <= d["declared_at"] <= a["p"] + 20), None)
            l = hit["declared_at"] - a["p"] if hit else np.nan
            lat.append(l)
            dist.append(a["dist"])
            dating.append(hit["window_start"] + 3.5 - a["p"] if hit else np.nan)
            recs.append({**a, "station": code, "latency": l})
        for d in ds:
            if not any(a["p"] - 0.5 <= d["declared_at"] <= a["p"] + 60 for a in arrs):
                unmatched += 1
    lat, dating, dist = np.asarray(lat), np.asarray(dating), np.asarray(dist)
    out = {"arrivals": len(lat), "detected": int(np.isfinite(lat).sum()),
           "latency_p50": float(np.nanmedian(lat)) if np.isfinite(lat).any() else None,
           "dating_abs_p50": float(np.nanmedian(np.abs(dating))) if np.isfinite(dating).any() else None,
           "unmatched": unmatched, "station_hours": hours,
           "unmatched_per_station_hour": unmatched / hours if hours else None}
    for d in DELAYS:
        out[f"recall@{d:g}s"] = float(np.mean(lat <= d)) if len(lat) else None
    by = {}
    for lo, hi in DIST_BINS:
        m = (dist >= lo) & (dist < hi)
        if m.any():
            by[f"{lo}-{hi} km"] = {"n": int(m.sum()), "recall@2s": float(np.mean(lat[m] <= 2)),
                                   "detected": float(np.mean(np.isfinite(lat[m]))),
                                   "latency_p50": float(np.nanmedian(lat[m])) if np.isfinite(lat[m]).any() else None}
    out["by_distance"] = by
    return out, recs


def evaluate(name, runs, truth):
    files, cat, start, _end, role = DATASETS[name]
    catalog = load_catalog(ROOT / "tests/catalogs" / f"{cat}.csv")
    res = {"role": role}
    arrivals = None
    for det, r in runs.items():
        stations, dets, span = parse_record(r["record"])
        if start:
            span = (max(span[0], UTCDateTime(start).timestamp), span[1])
        if arrivals is None:
            arrivals = truth.arrivals(r["files"], stations, catalog, span)
        hours = len(stations) * (span[1] - span[0]) / 3600
        st, recs = station_level(arrivals, dets, hours)
        ms, _ = cpu_per_hour(r["stdout"])
        res[det] = {"station": st, "network": score_dataset(name, r["stdout"]),
                    "detector_ms_per_station_hour": ms / hours if hours else None,
                    "arrivals": recs}
    res["refined_fraction"] = float(np.mean([a["refined"] for v in arrivals.values() for a in v])) \
        if arrivals and any(arrivals.values()) else None
    return res


# --- report ---------------------------------------------------------------------

def fmt(v, spec="{:.2f}", none="-"):
    return none if v is None or (isinstance(v, float) and not np.isfinite(v)) else spec.format(v)


def print_report(results):
    print("\nSTATION LEVEL (each detector on each station's stream)")
    head = f"  {'dataset':11s} {'detector':11s} {'arrivals':>8} {'det.':>5} " + \
           " ".join(f"{'≤' + format(d, 'g') + 's':>6}" for d in DELAYS) + \
           f" {'lat p50':>7} {'P-date':>7} {'unmatched/st-h':>14}"
    print(head)
    for name, r in results.items():
        for det in DETECTORS:
            s = r[det]["station"]
            print(f"  {name:11s} {det:11s} {s['arrivals']:>8} {s['detected']:>5} " +
                  " ".join(f"{fmt(s[f'recall@{d:g}s'], '{:.0%}'):>6}" for d in DELAYS) +
                  f" {fmt(s['latency_p50']):>7} {fmt(s['dating_abs_p50']):>7} "
                  f"{fmt(s['unmatched_per_station_hour']):>14}")
    print("\n  recall within 2 s of P, by epicentral distance (all datasets)")
    pooled = {det: defaultdict(list) for det in DETECTORS}
    for r in results.values():
        for det in DETECTORS:
            for a in r[det]["arrivals"]:
                for lo, hi in DIST_BINS:
                    if lo <= a["dist"] < hi:
                        pooled[det][f"{lo}-{hi} km"].append(a["latency"])
    for lo, hi in DIST_BINS:
        k = f"{lo}-{hi} km"
        cells = []
        for det in DETECTORS:
            v = np.asarray(pooled[det][k], dtype=float)
            cells.append(f"{det} {fmt(float(np.mean(v <= 2)) if len(v) else None, '{:.0%}'):>4}"
                         f" (ever {fmt(float(np.mean(np.isfinite(v))) if len(v) else None, '{:.0%}'):>4})")
        n = len(pooled[DETECTORS[0]][k])
        print(f"    {k:>12s}  n={n:<5d} " + "   ".join(cells))

    print("\nNETWORK LEVEL (the whole pipeline: association, picker, magnitude)")
    print(f"  {'dataset':11s} {'detector':11s} {'detected':>9} {'alarms':>7} {'unmatched':>9} "
          f"{'delay':>6} {'mag MAE':>8} {'false/day':>9} {'det. ms/st-h':>12}")
    for name, r in results.items():
        for det in DETECTORS:
            n = r[det]["network"]
            print(f"  {name:11s} {det:11s} {str(n['detected']) + '/' + str(n['catalogue']):>9} "
                  f"{n['alarms']:>7} {n['unmatched']:>9} {fmt(n['alarm_delay_median_s'], '{:.1f}'):>6} "
                  f"{fmt(n['magnitude_mae']):>8} {fmt(n.get('false_per_day'), '{:.1f}'):>9} "
                  f"{fmt(r[det]['detector_ms_per_station_hour'], '{:.0f}'):>12}")
    ref = [r["refined_fraction"] for r in results.values() if r["refined_fraction"] is not None]
    if ref:
        print(f"\n  P reference: AIC-refined for {np.mean(ref):.0%} of arrivals (TauP otherwise)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", default="build-release")
    ap.add_argument("--models", default="models")
    ap.add_argument("--args", default="", help="extra ayzek options for both runs")
    ap.add_argument("--only", help="comma-separated dataset names (tools/scorecard.py)")
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--json", help="write everything here")
    ap.add_argument("--workdir", help="keep each run's output and recording here "
                                      "(default: a temporary directory)")
    ap.add_argument("--reuse", action="store_true",
                    help="use the runs already in --workdir instead of replaying")
    a = ap.parse_args()

    names = a.only.split(",") if a.only else list(DATASETS)
    unknown = [n for n in names if n not in DATASETS]
    if unknown:
        sys.exit(f"unknown dataset(s): {', '.join(unknown)}; known: {', '.join(DATASETS)}")
    if not (ROOT / a.build / "app/ayzek").exists():
        sys.exit(f"no {a.build}/app/ayzek; build first")
    if not (Path(a.models) / "transformer.ayzw").exists():
        sys.exit(f"no {a.models}/transformer.ayzw; run tools/export_transformer.py")
    models = Path(a.models).resolve()
    extra = shlex.split(a.args)
    truth = Truth()
    results = {}
    tmpdir = None if a.workdir else tempfile.TemporaryDirectory()
    work = Path(a.workdir) if a.workdir else Path(tmpdir.name)
    work.mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:
        jobs = {(n, d): ex.submit(run, n, d, a.build, models, extra, work, a.reuse)
                for n in names for d in DETECTORS}
        for n in names:
            runs = {d: jobs[(n, d)].result() for d in DETECTORS}
            if any(r is None for r in runs.values()):
                print(f"  {n}: skipped (no data files)", file=sys.stderr)
                continue
            results[n] = evaluate(n, runs, truth)
            print(f"  {n}: done", file=sys.stderr)
    print_report(results)
    if a.json:
        Path(a.json).write_text(json.dumps(results, indent=1, default=float))
        print(f"\n  -> {a.json}")


if __name__ == "__main__":
    main()
