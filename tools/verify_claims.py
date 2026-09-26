#!/usr/bin/env python3
"""Checks every published claim against the artefact it came from.

    tools/verify_claims.py                     check everything present here
    tools/verify_claims.py --list              what each claim needs, and how to
                                               regenerate it
    tools/verify_claims.py --id ring_gap_len   one claim
    tools/verify_claims.py --missing           only the ones this machine cannot check

`tools/reproduce_all.sh` re-measures claims by running the pipeline again, which
takes hours and needs data that cannot be redistributed. This does the cheaper
half: it reads the files those runs already wrote and checks that the numbers
published in the paper and the docs are the numbers in the files. A drifted
figure is caught in seconds rather than after a five-hour rescan.

Each claim in `bench/claims.json` may carry:

    artifact   path to the file holding the measurement, relative to a repo
    probe      how to read the number out of it (below)
    repro      the command that regenerates the artefact from scratch
    expect     the published value
    tolerance  how far the measured value may sit from it

Probe syntax:

    json:a.b.c            nested key in a JSON document
    csv:COL@ROW=VAL       column COL of the row whose ROW column equals VAL
    csv:COL#N             column COL of row N (0-based)
    csvcount:             number of data rows
    csvcount:COL>=V       rows whose column COL passes the comparison
    grep:PATTERN          number of lines matching
    regex:PATTERN         first capture group of the first match
    lines:                number of lines

Claims with no `artifact` are left to `tools/reproduce_report.py`, which checks
them against a fresh run.
"""

import argparse
import csv as csvmod
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPOS = {
    "": ROOT,
    "ayzek_code": ROOT,
    "ayzek_paper": ROOT.parent / "ayzek_paper",
    "cascade_impl": ROOT.parent.parent.parent / "sismokaos" / "cascade_impl",
    "cnn_earthquake": ROOT.parent.parent.parent / "sismokaos" / "cnn_earthquake",
}


def resolve(artifact):
    """'repo:path' or a bare path relative to ayzek_code."""
    repo, _, rel = artifact.partition(":")
    if not rel:
        repo, rel = "", repo
    base = REPOS.get(repo)
    if base is None:
        return None
    return base / rel


def _rows(path):
    with open(path, newline="") as f:
        return list(csvmod.DictReader(f))


def _cmp(value, op, target):
    try:
        v = float(value)
        t = float(target)
    except (TypeError, ValueError):
        return str(value) == str(target)
    return {">=": v >= t, ">": v > t, "<=": v <= t, "<": v < t, "==": v == t}[op]


def probe(path, spec):
    """The number a probe names, or None if the file cannot answer."""
    kind, _, rest = spec.partition(":")
    if kind == "json":
        doc = json.loads(Path(path).read_text())
        for part in rest.split("."):
            if isinstance(doc, list):
                doc = doc[int(part)]
            else:
                doc = doc[part]
        return doc
    if kind == "lines":
        return sum(1 for _ in open(path))
    if kind == "grep":
        pat = re.compile(rest)
        return sum(1 for line in open(path, errors="replace") if pat.search(line))
    if kind == "regex":
        m = re.search(rest, Path(path).read_text(errors="replace"))
        return m.group(1) if m else None
    if kind == "csvcount":
        rows = _rows(path)
        if not rest:
            return len(rows)
        m = re.match(r"(\w+)(>=|<=|==|>|<)(.+)", rest)
        col, op, val = m.groups()
        return sum(1 for r in rows if _cmp(r.get(col), op, val))
    if kind == "csv":
        rows = _rows(path)
        if "#" in rest:
            col, _, idx = rest.partition("#")
            return rows[int(idx)][col]
        if "@" in rest:
            col, _, cond = rest.partition("@")
            key, _, val = cond.partition("=")
            for r in rows:
                if str(r.get(key)) == val:
                    return r[col]
            return None
        return rows[0][rest]
    raise ValueError(f"unknown probe {spec!r}")


def agrees(expect, got, tol):
    if got is None:
        return None
    if isinstance(expect, str):
        return str(got).strip() == expect
    try:
        return abs(float(got) - float(expect)) <= float(tol)
    except (TypeError, ValueError):
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--claims", default=str(ROOT / "bench" / "claims.json"))
    ap.add_argument("--id", action="append", help="check only these claim ids")
    ap.add_argument("--list", action="store_true", help="print what each claim needs")
    ap.add_argument("--missing", action="store_true", help="only unavailable claims")
    a = ap.parse_args()

    doc = json.loads(Path(a.claims).read_text())
    claims = doc if isinstance(doc, list) else doc.get("claims", doc)
    claims = [c for c in claims if c.get("artifact")]
    if a.id:
        claims = [c for c in claims if c["id"] in set(a.id)]

    if a.list:
        for c in claims:
            p = resolve(c["artifact"])
            print(f"{c['id']}\n    claim    {c['claim']}\n    artifact {c['artifact']}"
                  f"  [{'present' if p and p.exists() else 'MISSING'}]"
                  f"\n    repro    {c.get('repro', '-')}\n")
        return 0

    ok = bad = skipped = 0
    print(f"{'claim':<54} {'published':>13} {'in the file':>13}")
    print("-" * 84)
    for c in claims:
        p = resolve(c["artifact"])
        got = None
        if p and p.exists():
            try:
                got = probe(p, c["probe"])
            except Exception as e:
                got = f"probe failed: {type(e).__name__}"
        verdict = agrees(c["expect"], got, c.get("tolerance", 0))
        if a.missing and verdict is not None:
            continue
        if verdict is None:
            mark, skipped = "not checked", skipped + 1
        elif verdict:
            mark, ok = "ok", ok + 1
        else:
            mark, bad = "DISAGREES", bad + 1
        shown = got if got is not None else "-"
        if isinstance(shown, float):
            shown = f"{shown:.4g}"
        print(f"{c['claim'][:54]:<54} {str(c['expect']):>13} {str(shown)[:13]:>13}  {mark}")
        if verdict is None and not a.missing:
            print(f"{'':<54} needs {c['artifact']}")

    print("-" * 84)
    print(f"{ok} agree, {bad} disagree, {skipped} not checkable here")
    if skipped:
        print("\nregenerate a missing artefact with its `repro` command "
              "(tools/verify_claims.py --list)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
