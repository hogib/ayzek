"""Compares two ayzek_bench results: per benchmark, the p50 of each and the ratio.

    python3 tools/bench_compare.py before.json after.json
    python3 tools/bench_compare.py x86.json pi.json --threshold 0     # just the ratios

A benchmark more than --threshold percent slower in the second file is flagged,
and the exit status is 1 if any is, so this can gate a change. Machines and
builds are printed first: a comparison across different CPUs is a port, not a
regression.

Standard library only.
"""
import argparse
import json
import sys


def load(path):
    with open(path) as f:
        return json.load(f)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--threshold", type=float, default=10.0,
                    help="flag benchmarks this many percent slower in B (0: never flag)")
    args = ap.parse_args()
    A, B = load(args.a), load(args.b)
    for tag, r in (("A", A), ("B", B)):
        print(f"{tag}: {r['cpu']}, {r['cores']} cores, {r['backend']}, "
              f"{'optimised' if r['optimized'] else 'UNOPTIMISED'}, {r['version']}, {r['date']}")
    if A["cpu"] != B["cpu"]:
        print("   (different CPUs: ratios compare machines, not code)")

    key = lambda r: (r["group"], r["name"])
    a = {key(r): r for r in A["results"]}
    b = {key(r): r for r in B["results"]}
    flagged = 0
    print(f"\n{'':42s} {'A p50 ms':>11} {'B p50 ms':>11} {'B/A':>7}")
    group = None
    for k in list(dict.fromkeys(list(a) + list(b))):
        if k[0] != group:
            group = k[0]
            print(group)
        ra, rb = a.get(k), b.get(k)
        if not (ra and rb):
            print(f"  {k[1][:40]:40s} {'-' if not ra else f'{ra['p50_ms']:.4f}':>11} "
                  f"{'-' if not rb else f'{rb['p50_ms']:.4f}':>11}")
            continue
        ratio = rb["p50_ms"] / ra["p50_ms"] if ra["p50_ms"] else float("nan")
        slow = args.threshold > 0 and ratio > 1 + args.threshold / 100
        flagged += slow
        print(f"  {k[1][:40]:40s} {ra['p50_ms']:>11.4f} {rb['p50_ms']:>11.4f} {ratio:>6.2f}x"
              f"{'  SLOWER' if slow else ''}")

    ca = {c["threads"]: c for c in A.get("capacity", [])}
    cb = {c["threads"]: c for c in B.get("capacity", [])}
    if ca or cb:
        print(f"\ncapacity (stations at one window per 0.5 s)\n  {'threads':>7} {'A':>8} {'B':>8}")
        for t in sorted(set(ca) | set(cb)):
            fa = f"{ca[t]['stations']:.0f}" if t in ca else "-"
            fb = f"{cb[t]['stations']:.0f}" if t in cb else "-"
            print(f"  {t:>7} {fa:>8} {fb:>8}")
    for tag, r in (("A", A), ("B", B)):
        m = r.get("magnitude_latency_ms") or {}
        if m.get("loaded_p50"):
            print(f"magnitude under load {tag}: p50 {m['loaded_p50']:.1f} ms, p99 {m['loaded_p99']:.1f} ms")
    if flagged:
        print(f"\n{flagged} benchmark(s) more than {args.threshold:g}% slower in B")
    return 1 if flagged else 0


if __name__ == "__main__":
    sys.exit(main())
