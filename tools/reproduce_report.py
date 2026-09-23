"""Compares the claims in bench/claims.json with what a reproduction run measured.

    python3 tools/reproduce_report.py runs/reproduce-20260923T101500Z

Reads the JSON files `tools/reproduce_all.sh` wrote in that directory and
prints one line per claim: what is published, what this machine measured, and
whether they agree within the claim's tolerance. Exits 1 if any claim that
could be checked disagrees, so the script can gate a change.

A claim whose data is missing is reported as "not checked" rather than passed:
the AFAD replays cannot be redistributed, so most machines will check the KO
claims only.

Key syntax, against the JSON of the named source:
    summary.<name>                     a scorecard headline
    event.<dataset>.<type M>.<field>   one matched event's field
    result.<benchmark>.<field>         one ayzek_bench result
    capacity.<threads>.<field>         one capacity row

Standard library only.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(out, name):
    p = out / f"{name}.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except ValueError:
        return None


def lookup(doc, key):
    """The value a claim's key names, or None when it is not in this run."""
    kind, _, rest = key.partition(".")
    if doc is None:
        return None
    if kind == "summary":
        return doc.get("summary", {}).get(rest)
    if kind == "event":
        dataset, _, rest2 = rest.partition(".")
        label, _, field = rest2.rpartition(".")
        for e in doc.get("datasets", {}).get(dataset, {}).get("events", []):
            if f"{e['afad_type']} {e['afad_M']:.1f}" == label:
                return e.get(field)
        return None
    if kind == "result":
        name, _, field = rest.rpartition(".")
        for r in doc.get("results", []):
            if r["name"] == name:
                return r.get(field)
        return None
    if kind == "capacity":
        threads, _, field = rest.partition(".")
        for c in doc.get("capacity", []):
            if str(c["threads"]) == threads:
                return c.get(field)
        return None
    return None


def agrees(expect, got, tol):
    if got is None:
        return None
    if isinstance(expect, str):
        return str(got) == expect
    try:
        return abs(float(got) - float(expect)) <= float(tol)
    except (TypeError, ValueError):
        return False


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: reproduce_report.py RUN_DIR")
    out = Path(sys.argv[1])
    spec = json.loads((ROOT / "bench" / "claims.json").read_text())
    docs = {name: load(out, name) for name in
            ("scorecard", "scorecard_stalta_eqfa", "scorecard_stalta_eqdet",
             "bench_stages", "bench_layers", "bench_capacity")}

    ok = bad = skipped = 0
    print(f"{'claim':52s} {'published':>12} {'measured':>12}  where")
    for c in spec["claims"]:
        got = lookup(docs.get(c["source"]), c["key"])
        verdict = agrees(c["expect"], got, c.get("tolerance", 0))
        shown = ("-" if got is None else
                 f"{got:.2f}" if isinstance(got, float) else str(got))
        mark = {True: "ok", False: "DIFFERS", None: "not checked"}[verdict]
        ok += verdict is True
        bad += verdict is False
        skipped += verdict is None
        print(f"{c['claim'][:52]:52s} {str(c['expect']):>12} {shown:>12}  {mark}"
              f"{'' if verdict is not False else '  <- ' + c['where']}")

    print(f"\n{ok} claims reproduced, {bad} differ, {skipped} not checked here")
    if skipped:
        print("\nnot checked, most likely because the dataset is not on this machine:")
        for c in spec["claims"]:
            if lookup(docs.get(c["source"]), c["key"]) is None:
                print(f"  {c['id']:24s} needs {c['source']}")
    print("\nclaims no run on this machine can check:")
    for n in spec["not_reproducible_here"]:
        print(f"  {n['claim']}\n      {n['how']}")

    sc = docs.get("scorecard")
    if sc and sc.get("skipped"):
        print("\ndatasets skipped by the scorecard:")
        for name, why in sc["skipped"].items():
            print(f"  {name:12s} {why}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
