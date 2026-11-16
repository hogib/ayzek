"""Matched pairs on the exploratory trade-off grid: EXPLORATORY, outside the held-out protocol.

    uv run --offline --project tools python tools/heldout_grid_matched.py [--scratch DIR] [--jobs 4]
                                                                  (from ayzek_code/)

Reads data/runs/heldout_grid/summary.csv (tools/heldout_grid_summary.py --grid tradeoff)
and, per network, for each learned grid setting L and each baseline family B (stalta,
phasenet, eqt), picks baseline settings by a fixed rule (RULE below; applied mechanically,
nothing hand-picked), then compares L with them with tools/significance.py, unchanged.

significance.py always compares the label `learned`. So per network and L a scratch
outcomes directory is built: rewritten copies of L's outcome tables under the label
`learned` (file learned_<set>.csv, detector column `learned`), next to links to every
grid setting's own tables. It is then called once per network and L with --against all
grid settings in the order of grid_settings.json, exactly as tools/heldout_grid_summary.py
calls it; the random stream is shared along that list, so the same list and order keep
every interval at the position it has in summary.csv.

Self-check, before anything else is written: with L = g_learned_0p8 the detected and
corrected counts, the difference and its interval against every setting must equal
summary.csv's (frozen learned minus setting) exactly; otherwise the script stops nonzero.

Writes under data/runs/heldout_grid/matched/: <net>/<L>.txt and .json (raw significance
outputs), matched_pairs.csv, provenance.json, DONE. Files are written to temp names and
renamed; a finished <net>/<L>.txt is reused on a rerun. Exits nonzero if nothing is produced.
"""
import argparse, csv, hashlib, json, math, platform, re, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

G = Path("data/runs/heldout_grid")
OUT = G / "matched"
NETS = ("greece_now", "nz_now", "ko_now")
FAMILIES = ("stalta", "phasenet", "eqt")
LO, HI = 0.8, 1.25
RULE = {
    "scope": "per network, per learned grid setting L (g_learned_*), per baseline family B in stalta, phasenet, eqt; "
             "counts from summary.csv; ratio = count(L) / count(B); grid order = order of grid_settings.json",
    "matched_unmatched": "the B setting whose unmatched_alarms is closest to L's on a log scale (|ln ratio| smallest) among "
                         f"the B settings with a count > 0; flag 'matched' if {LO} <= ratio <= {HI}, else 'not matched (ratio r)'. "
                         "If L's count is 0 or no B setting has a count > 0: closest by absolute difference among all B "
                         "settings, flag 'zero count'.",
    "matched_quiet": "the same with quiet_alarms.",
    "bracket_below / bracket_above": "the B settings with the largest unmatched_alarms strictly below L's and the smallest "
                                     "strictly above L's; both rows are written only if both exist (flag 'bracket').",
    "ties": "equal distance, or equal counts: the first in grid order (lowest threshold); for bracket_above the last in "
            "grid order with that count. A tie is noted in the column `tie`.",
    "statistics": "tools/significance.py, default seed, `learned` = L, --against every grid setting in grid order "
                  "(the call of tools/heldout_grid_summary.py); difference = L minus B, chance-corrected detections; "
                  "95% paired bootstrap over events.",
}
SIG = re.compile(r"^learned vs (\S+): (\d+) events\n  detected\s+(\d+) vs (\d+).*\n"
                 r"  chance-corrected\s+([\d.]+) vs ([\d.]+), difference ([+-][\d.]+) "
                 r"\(95% bootstrap ([+-][\d.]+) to ([+-][\d.]+), P\(diff <= 0\) = ([\d.]+)\)\n"
                 r"  discordant\s+learned only (\d+), \S+ only (\d+), exact McNemar p = (\S+)", re.M)
KEYS = ("events", "detected_L", "detected_B", "corrected_L", "corrected_B", "difference", "ci_low", "ci_high",
        "p_diff_le0", "L_only", "B_only", "mcnemar_p")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def pick(L, cands, col):
    """(setting, flag, tie) for the matched rule on column `col`; cands in grid order."""
    x = int(L[col])
    pos = [c for c in cands if int(c[col]) > 0]
    if x == 0 or not pos:
        d = [abs(x - int(c[col])) for c in cands]
        i = d.index(min(d))
        return cands[i], "zero count", d.count(min(d)) > 1
    d = [abs(math.log(x / int(c[col]))) for c in pos]
    i = d.index(min(d))
    r = x / int(pos[i][col])
    flag = "matched" if LO <= r <= HI else f"not matched (ratio {r:.2f})"
    return pos[i], flag, sum(abs(v - min(d)) < 1e-12 for v in d) > 1


def bracket(L, cands):
    x = int(L["unmatched_alarms"])
    below = [c for c in cands if int(c["unmatched_alarms"]) < x]
    above = [c for c in cands if int(c["unmatched_alarms"]) > x]
    if not below or not above:
        return []
    vb = max(int(c["unmatched_alarms"]) for c in below)
    va = min(int(c["unmatched_alarms"]) for c in above)
    b = [c for c in below if int(c["unmatched_alarms"]) == vb]
    a = [c for c in above if int(c["unmatched_alarms"]) == va]
    return [("bracket_below", b[0], "bracket", len(b) > 1), ("bracket_above", a[-1], "bracket", len(a) > 1)]


def run(net, L, settings, scratch):
    """One significance.py call: L as `learned` against every grid setting. Returns {setting: figures}."""
    txt = OUT / net / f"{L}.txt"
    if not txt.exists():
        d = scratch / net / L
        d.mkdir(parents=True, exist_ok=True)
        src = G / net / "outcomes"
        for f in src.glob("*.csv"):
            if not (d / f.name).exists():
                (d / f.name).symlink_to(f.resolve())
        seqs = []
        for f in sorted(src.glob(f"{L}_*.csv")):
            seq = f.stem[len(L) + 1:]
            seqs.append(seq)
            rows = list(csv.DictReader(open(f)))
            if not rows or {r["detector"] for r in rows} != {L}:
                sys.exit(f"{f}: empty, or a detector label other than {L}")
            with open(d / f"learned_{seq}.csv", "w", newline="") as fh:
                w = csv.DictWriter(fh, list(rows[0]))
                w.writeheader()
                w.writerows({**r, "detector": "learned"} for r in rows)
        if not seqs:
            sys.exit(f"{net} {L}: no outcome tables")
        (OUT / net).mkdir(parents=True, exist_ok=True)
        js = OUT / net / f"{L}.json"
        p = subprocess.run([sys.executable, "tools/significance.py", "--outcomes", str(d),
                            "--groups", f"{net}={','.join(seqs)}", "--against", ",".join(settings),
                            "--json", str(js) + ".tmp"], capture_output=True, text=True)
        if p.returncode:
            sys.exit(f"{net} {L}: significance.py failed\n{p.stderr[-2000:]}")
        Path(str(js) + ".tmp").rename(js)
        Path(str(txt) + ".tmp").write_text(p.stdout)
        Path(str(txt) + ".tmp").rename(txt)
    sig = {m.group(1): dict(zip(KEYS, m.groups()[1:])) for m in SIG.finditer(txt.read_text())}
    if set(sig) != set(settings):
        sys.exit(f"{txt}: {len(sig)} comparisons parsed, {len(settings)} expected")
    print(net, L, "done", flush=True)
    return sig


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--scratch", help="directory for the relabelled outcome tables (default: a new temp directory)")
    ap.add_argument("--jobs", type=int, default=4, help="significance.py processes at a time")
    a = ap.parse_args()
    scratch = Path(a.scratch or tempfile.mkdtemp(prefix="grid_matched_"))
    settings = [k for k in json.loads((G / "grid_settings.json").read_text()) if not k.startswith("_")]
    summ = {(r["network"], r["setting"]): r for r in csv.DictReader(open(G / "summary.csv"))}
    learned = [s for s in settings if s.startswith("g_learned_")]
    OUT.mkdir(exist_ok=True)

    # Self-check first: g_learned_0p8 as `learned` must reproduce summary.csv.
    with ThreadPoolExecutor(a.jobs) as pool:
        chk = dict(zip(NETS, pool.map(lambda n: run(n, "g_learned_0p8", settings, scratch), NETS)))
    bad = []
    for net in NETS:
        for s in settings:
            g, r = chk[net][s], summ[net, s]
            got = (g["events"], g["detected_B"], g["corrected_B"], g["corrected_L"], g["difference"], g["ci_low"], g["ci_high"])
            want = (r["events"], r["detected"], r["corrected"], r["frozen_learned_corrected"], r["difference"], r["ci_low"], r["ci_high"])
            if got != want:
                bad.append(f"{net} {s}: got {got}, summary.csv {want}")
    if bad:
        sys.exit("self-check FAILED (g_learned_0p8 as learned vs summary.csv):\n" + "\n".join(bad))
    check = f"g_learned_0p8 as learned reproduces summary.csv exactly: {len(NETS) * len(settings)} comparisons x 7 fields"
    print("self-check passed:", check, flush=True)

    jobs = [(n, L) for n in NETS for L in learned]
    with ThreadPoolExecutor(a.jobs) as pool:
        sig = dict(zip(jobs, pool.map(lambda j: run(j[0], j[1], settings, scratch), jobs)))

    rows = []
    for net in NETS:
        for L in learned:
            rl = summ[net, L]
            for fam in FAMILIES:
                cands = [summ[net, s] for s in settings if s.startswith(f"g_{fam}_")]
                picks = [("matched_unmatched", *pick(rl, cands, "unmatched_alarms")),
                         ("matched_quiet", *pick(rl, cands, "quiet_alarms"))] + bracket(rl, cands)
                for rule, rb, flag, tie in picks:
                    col = "quiet_alarms" if rule == "matched_quiet" else "unmatched_alarms"
                    x, y = int(rl[col]), int(rb[col])
                    g = sig[net, L][rb["setting"]]
                    if (g["detected_L"], g["detected_B"]) != (rl["detected"], rb["detected"]):
                        sys.exit(f"{net} {L} vs {rb['setting']}: detected counts differ from summary.csv")
                    rows.append({"network": net, "learned_setting": L, "family": fam, "rule": rule,
                                 "baseline_setting": rb["setting"], "matched_on": col,
                                 "ratio_L_over_B": f"{x / y:.3f}" if y else "", "flag": flag, "tie": "tie" if tie else "",
                                 "unmatched_alarms_L": rl["unmatched_alarms"], "unmatched_alarms_B": rb["unmatched_alarms"],
                                 "quiet_alarms_L": rl["quiet_alarms"], "quiet_alarms_B": rb["quiet_alarms"],
                                 "quiet_hours": rl["quiet_hours"], **g})
    if not rows:
        sys.exit("nothing produced")
    tmp = OUT / "matched_pairs.csv.tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    git = lambda *c: subprocess.run(["git", *c], capture_output=True, text=True).stdout.strip()
    prov = {"note": "exploratory, outside the held-out protocol", "argv": sys.argv, "rule": RULE,
            "commit": git("rev-parse", "HEAD"), "uncommitted_changes_tools_src": git("status", "--porcelain", "tools", "src"),
            "python": sys.version, "platform": platform.platform(),
            "significance_py": {"seed": "default (20261022)", "against_order": settings,
                                "sha256": sha("tools/significance.py")},
            "inputs_sha256": {"summary.csv": sha(G / "summary.csv"), "grid_settings.json": sha(G / "grid_settings.json"),
                              "tools/heldout_grid_matched.py": sha(__file__)},
            "self_check": check, "rows": len(rows), "significance_calls": len(jobs)}
    ptmp = OUT / "provenance.json.tmp"
    ptmp.write_text(json.dumps(prov, indent=1, ensure_ascii=False) + "\n")
    tmp.rename(OUT / "matched_pairs.csv")
    ptmp.rename(OUT / "provenance.json")
    (OUT / "DONE").write_text("")
    print("wrote", OUT / "matched_pairs.csv", len(rows), "rows")


if __name__ == "__main__":
    main()
