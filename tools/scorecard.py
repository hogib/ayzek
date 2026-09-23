"""The scorecard: one set of numbers for judging a change to ayzek or its models.

Runs ayzek over fixed datasets with fixed catalogues and reports detection,
unmatched alarms, false alarms on quiet data, magnitude and location on the
large events, and (with --bench) the speed of each stage. Two scorecards can be
compared side by side.

    python3 tools/scorecard.py --json score.json                 # current models
    python3 tools/scorecard.py --models /tmp/new --json new.json # candidate models
    python3 tools/scorecard.py --args '--detector stalta' --json stalta.json
    python3 tools/scorecard.py --compare score.json new.json

Datasets (data/, not tracked; a missing one is skipped and says so):

  demo        Sindirgi 2025-11-10, 30 min, DEMI MANT BAND     tuning (AFAD)
  marmara     Marmara Sea 2025-04-23, 3 h, 8 stations         tuning (AFAD)
  sindirgi61  Sindirgi 2025-08-10, 3 h, 6 stations            held out (AFAD)
  noise_*     four catalogue-quiet windows, 33 h in total     false alarms (AFAD)
  demo_ko     the same Mw 4.9 on six KO stations              public
  marmara_ko  the Mw 6.2 and the hour after it, eight KO      public
  quiet_ko    6 h on six KO stations, no event within 250 km  false alarms (public)

The AFAD sets are the ones the published numbers were measured on and cannot be
redistributed; the KO sets ship with the release, so a score on them is
reproducible by anyone (`tools/fetch_ko_data.py`).

Catalogues are the excerpts in tests/catalogs/ (AFAD events within a day of each
dataset), so a score does not move when AFAD revises its catalogue.

Standard library only.
"""
import argparse
import concurrent.futures
import datetime as dt
import glob
import json
import re
import shlex
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from sweep_trigger import BANDS, PATTERNS, parse_report  # noqa: E402

# name: (files, catalogue, --from or None, end of the quiet period or None, role)
DATASETS = {
    "demo": (["data/demo/DEMI.mseed", "data/demo/MANT.mseed", "data/demo/BAND.mseed"],
             "demo", None, None, "tuning"),
    "marmara": (["data/marmara/*.mseed"], "marmara", None, None, "tuning"),
    "sindirgi61": (["data/sindirgi61/*.mseed"], "sindirgi61", None, None, "held out"),
    "noise_a1": (["data/noise_a1/*.mseed"], "noise_a1",
                 "2024-06-05T01:20:00", "2024-06-05T07:16:00", "noise"),
    "noise_a2": (["data/noise_a2/*.mseed"], "noise_a2",
                 "2024-06-05T07:56:00", "2024-06-05T12:15:00", "noise"),
    "noise_b": (["data/noise_b/*.mseed"], "noise_b",
                "2024-10-02T13:00:00", "2024-10-03T01:08:00", "noise"),
    "noise_c": (["data/noise_c/*.mseed"], "noise_c",
                "2025-03-31T04:17:00", "2025-03-31T15:10:00", "noise"),
    # KOERI (network KO) sets, the ones the release ships: open data under
    # citation, so these are the datasets anyone can reproduce a score on.
    # `tools/fetch_ko_data.py` rebuilds them from the KOERI FDSN service.
    "demo_ko": (["data/demo_ko/*.mseed"], "demo_ko", None, None, "public"),
    "marmara_ko": (["data/marmara_ko/*.mseed"], "marmara_ko", None, None, "public"),
    "quiet_ko": (["data/quiet_ko/*.mseed"], "quiet_ko",
                 "2024-01-03T11:35:00", "2024-01-03T17:25:00", "public noise"),
}

# One matched alarm's report block: ayzek's magnitudes, and the AFAD comparison.
EVENT = re.compile(
    r"^Event #\d+: detected M([\d.]+) .*\n"
    r"  alarm\s+(\S+) UTC by (.*)\n"
    r"  magnitude\s+M([\d.]+) (?:at the alarm|[\d.]+ s after the alarm), M([\d.]+) final.*\n"
    r"(?:  location .*\n)?"
    r"  AFAD\s+(\w+) ([\d.]+) at (\S+) UTC: alarm ([\d.]+) s after origin, "
    r"magnitude ([+-][\d.]+)(?:, epicentre ([\d.]+) km off)?", re.M)


def files_of(patterns):
    out = []
    for p in patterns:
        out += sorted(glob.glob(str(ROOT / p)))
    return out


def run_ayzek(name, build, models, extra):
    files, cat, start, end, role = DATASETS[name]
    paths = files_of(files)
    if not paths:
        return name, None, f"no files for {', '.join(files)}"
    cmd = [str(ROOT / build / "app/ayzek"), "--speed", "0", "--no-color",
           "--models", str(models), "--catalog", str(ROOT / "tests/catalogs" / f"{cat}.csv")]
    if start:
        cmd += ["--from", start]
    cmd += extra + paths
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        return name, None, f"ayzek exited {p.returncode}: {p.stderr.strip()[-300:]}"
    return name, p.stdout, None


def score_dataset(name, out):
    files, _, start, end, role = DATASETS[name]
    row = parse_report(out)
    num = lambda k: int(row[k]) if row.get(k, "") != "" else 0
    s = {"role": role, "stations": len(files_of(files)),
         "alarms": num("alarms"), "matched": num("correct"), "unmatched": num("false"),
         "catalogue": num("catalogue"), "detected": num("detected"),
         "alarm_delay_median_s": float(row["alarm_delay_s"]) if row["alarm_delay_s"] else None,
         "bands": {b: row[b] for b in BANDS if row.get(b)}}
    m = re.search(r"magnitude error\s+mean ([\d.]+) units", out)
    s["magnitude_mae"] = float(m.group(1)) if m else None
    m = re.search(r"alarm before the S wave\s+at (\d+) of (\d+) station arrivals, median warning ([-\d.]+) s", out)
    s["before_s"] = f"{m.group(1)}/{m.group(2)}" if m else None
    s["warning_median_s"] = float(m.group(3)) if m else None
    events = []
    for e in EVENT.finditer(out):
        events.append({"afad_type": e.group(6), "afad_M": float(e.group(7)),
                       "origin": e.group(8), "alarm_delay_s": float(e.group(9)),
                       "M_alarm": float(e.group(4)), "M_final": float(e.group(5)),
                       "M_error": float(e.group(10)),
                       "epicentre_km": float(e.group(11)) if e.group(11) else None,
                       "stations": e.group(3).replace("; later", " +")})
    s["events"] = events
    if "noise" in role:
        hours = (dt.datetime.fromisoformat(end) - dt.datetime.fromisoformat(start)).total_seconds() / 3600
        s["hours"] = hours
        s["false_per_day"] = 24 * s["unmatched"] / hours
    s["dropped_warnings"] = row.get("dropped_warnings", 0)
    return s


def bench(build):
    """p50 of the stages that set the real-time budget, from ayzek_bench."""
    exe = ROOT / build / "bench/ayzek_bench"
    if not exe.exists():
        return None
    out = ROOT / build / "bench/scorecard_stages.json"
    p = subprocess.run([str(exe), "--models", str(ROOT / "models"), "--group", "stages",
                        "--min-time", "0.5", "--json", str(out)], capture_output=True, text=True)
    if p.returncode != 0:
        return {"error": p.stderr.strip()[-300:]}
    r = json.loads(out.read_text())
    keep = {"per-window total", "picker", "magnitude estimate (3 models)"}
    return {"cpu": r["cpu"], "version": r["version"],
            **{x["name"]: x["p50_ms"] for x in r["results"] if x["name"] in keep}}


def summarise(sets):
    """The headline numbers over all datasets that ran."""
    def total(role, key):
        return sum(v[key] for v in sets.values() if v["role"] == role)
    h = {}
    roles = list(dict.fromkeys(v["role"] for v in sets.values()))
    for role in [r for r in roles if "noise" not in r]:
        h[f"detected ({role})"] = f"{total(role, 'detected')}/{total(role, 'catalogue')}"
        h[f"unmatched alarms ({role})"] = total(role, "unmatched")
    # False alarms are pooled per role, so the AFAD windows and the KO one are
    # not mixed into a rate that belongs to neither.
    for role in [r for r in roles if "noise" in r]:
        noise = [v for v in sets.values() if v["role"] == role]
        hrs = sum(v["hours"] for v in noise)
        h[f"false alarms per day ({role})"] = round(24 * sum(v["unmatched"] for v in noise) / hrs, 2)
        h[f"hours ({role})"] = round(hrs, 1)
    evs = [e for v in sets.values() if "noise" not in v["role"] for e in v["events"]]
    if evs:
        h["magnitude MAE, all matched"] = round(statistics.mean(abs(e["M_error"]) for e in evs), 3)
        big = [e for e in evs if e["afad_M"] >= 4.0]
        if big:
            h["magnitude MAE, AFAD M>=4"] = round(statistics.mean(abs(e["M_error"]) for e in big), 3)
            h["magnitude bias, AFAD M>=4"] = round(statistics.mean(e["M_error"] for e in big), 3)
        loc = [e["epicentre_km"] for e in evs if e["epicentre_km"] is not None]
        if loc:
            h["epicentre error median km"] = round(statistics.median(loc), 1)
    return h


def print_card(card):
    print(f"\nscorecard  models={card['models']}  args={card['args'] or '-'}  "
          f"({card['date']}, {card['version']})")
    print("\n  " + "\n  ".join(f"{k:34s} {v}" for k, v in card["summary"].items()))
    print(f"\n  {'dataset':11s} {'role':13s} {'st':>3} {'detected':>9} {'alarms':>7} "
          f"{'unmatched':>9} {'delay':>6} {'mag MAE':>8} {'false/day':>9}")
    for n, v in card["datasets"].items():
        print(f"  {n:11s} {v['role']:13s} {v['stations']:>3} "
              f"{str(v['detected']) + '/' + str(v['catalogue']):>9} {v['alarms']:>7} "
              f"{v['unmatched']:>9} "
              f"{(str(v['alarm_delay_median_s']) if v['alarm_delay_median_s'] else '-'):>6} "
              f"{(str(v['magnitude_mae']) if v['magnitude_mae'] else '-'):>8} "
              f"{(f'{v['false_per_day']:.1f}' if 'false_per_day' in v else '-'):>9}")
    for n, why in card["skipped"].items():
        print(f"  {n:11s} skipped: {why}")
    big = [(n, e) for n, v in card["datasets"].items() for e in v["events"] if e["afad_M"] >= 4.0]
    if big:
        print("\n  large events (AFAD M>=4)")
        for n, e in big:
            print(f"    {n:11s} {e['afad_type']} {e['afad_M']:.1f} at {e['origin']}: alarm +{e['alarm_delay_s']} s, "
                  f"M{e['M_alarm']} at the alarm, M{e['M_final']} final ({e['M_error']:+.1f}), "
                  f"epicentre {e['epicentre_km'] if e['epicentre_km'] is not None else '-'} km")
    if card.get("bench"):
        b = card["bench"]
        print(f"\n  speed ({b.get('cpu', '?')}): " + ", ".join(
            f"{k} {v:.1f} ms" for k, v in b.items() if isinstance(v, float)))


def compare(a, b):
    A, B = json.loads(Path(a).read_text()), json.loads(Path(b).read_text())
    print(f"{'':36s} {'A':>14} {'B':>14}")
    print(f"{'':36s} {Path(a).name[:14]:>14} {Path(b).name[:14]:>14}")
    keys = list(dict.fromkeys(list(A["summary"]) + list(B["summary"])))
    for k in keys:
        print(f"  {k:34s} {str(A['summary'].get(k, '-')):>14} {str(B['summary'].get(k, '-')):>14}")
    ea = {(n, e["origin"]): e for n, v in A["datasets"].items() for e in v["events"] if e["afad_M"] >= 4}
    eb = {(n, e["origin"]): e for n, v in B["datasets"].items() for e in v["events"] if e["afad_M"] >= 4}
    if ea or eb:
        print("\n  large events: final magnitude error, A -> B")
        for key in sorted(set(ea) | set(eb)):
            e = ea.get(key) or eb[key]
            fa = f"{ea[key]['M_error']:+.1f}" if key in ea else "missed"
            fb = f"{eb[key]['M_error']:+.1f}" if key in eb else "missed"
            print(f"    {key[0]:11s} {e['afad_type']} {e['afad_M']:.1f} {key[1]}: {fa:>7} -> {fb}")
    for side, card in (("A", A), ("B", B)):
        if card.get("bench"):
            print(f"\n  speed {side}: " + ", ".join(f"{k} {v:.1f} ms" for k, v in card["bench"].items()
                                                if isinstance(v, float)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", default="build-release")
    ap.add_argument("--models", default="models")
    ap.add_argument("--args", default="", help="extra ayzek options, e.g. '--detector stalta'")
    ap.add_argument("--only", help="comma-separated dataset names")
    ap.add_argument("--jobs", type=int, default=2, help="datasets run at once")
    ap.add_argument("--bench", action="store_true", help="also time the stages with ayzek_bench")
    ap.add_argument("--json", help="write the scorecard here")
    ap.add_argument("--save-output", metavar="DIR",
                    help="also keep each dataset's ayzek output, so a parser fix "
                         "does not need the replays run again "
                         "(default: <json>.out/ when --json is given)")
    ap.add_argument("--compare", nargs=2, metavar=("A", "B"), help="compare two scorecards and exit")
    a = ap.parse_args()
    if a.compare:
        return compare(*a.compare)

    names = a.only.split(",") if a.only else list(DATASETS)
    unknown = [n for n in names if n not in DATASETS]
    if unknown:
        sys.exit(f"unknown dataset(s): {', '.join(unknown)}; known: {', '.join(DATASETS)}")
    if not (ROOT / a.build / "app/ayzek").exists():
        sys.exit(f"no {a.build}/app/ayzek; build first")
    models = Path(a.models).resolve()
    extra = shlex.split(a.args)
    version = subprocess.run(["git", "-C", str(ROOT), "describe", "--always", "--dirty"],
                             capture_output=True, text=True).stdout.strip()

    outdir = Path(a.save_output) if a.save_output else (Path(a.json + ".out") if a.json else None)
    if outdir:
        outdir.mkdir(parents=True, exist_ok=True)
    sets, skipped = {}, {}
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:
        for name, out, err in ex.map(lambda n: run_ayzek(n, a.build, models, extra), names):
            if err:
                skipped[name] = err
                print(f"  {name}: skipped ({err})", file=sys.stderr)
            else:
                sets[name] = score_dataset(name, out)
                if outdir:
                    (outdir / f"{name}.txt").write_text(out)
                print(f"  {name}: done", file=sys.stderr)
    card = {"version": version, "date": dt.datetime.now(dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "models": str(models), "args": a.args, "datasets": sets, "skipped": skipped,
            "summary": summarise(sets)}
    if a.bench:
        card["bench"] = bench(a.build)
    print_card(card)
    if a.json:
        Path(a.json).write_text(json.dumps(card, indent=1) + "\n")
        print(f"\nscorecard -> {a.json}")


if __name__ == "__main__":
    sys.exit(main())
