# Every claim, and how to reproduce it

A number that appears in the paper, the README or these docs is registered in
`bench/claims.json` with the artefact it came from and the command that
regenerates that artefact. Nothing is quoted from memory.

There are two ways to check a claim, cheap and expensive, and both are wired up.

## Cheap: is the published number still the number in the file?

```
tools/verify_claims.py             check every artefact present on this machine
tools/verify_claims.py --list      what each claim needs, and its repro command
tools/verify_claims.py --missing   only the claims this machine cannot check
tools/verify_claims.py --id X      one claim
```

Seconds, no replay, no GPU. It reads the CSVs and reports that earlier runs
wrote and compares them with the published figures. This catches the common
failure — a figure edited in one place and not another, or an artefact
regenerated while the text stayed put — which a five-hour rescan would also
catch but much later. It runs as stage 0 of `tools/reproduce_all.sh`.

A claim whose artefact is absent is reported "not checked", never "passed".

## Expensive: measure it again from the waveforms

```
tools/reproduce_all.sh             stages 0-5: verify, build, scorecard,
                                   STA/LTA, benchmarks, aarch64
tools/reproduce_all.sh --quick     skip the STA/LTA comparison
```

`tools/reproduce_report.py` then compares the fresh run with `claims.json`
and exits non-zero on a disagreement.

## What each family costs, and what it needs

| family | claims | cost | needs |
|---|---|---|---|
| detector scorecard, alarms, location | 13 | ~20 min | the replay datasets |
| STA/LTA comparison | 4 | ~15 min | same |
| benchmarks, capacity | 4 | ~5 min | a release build |
| ring-buffer fix (R1) | 8 | ~40 min | `data/noise_a1`, `noise_a2`, `data/marmara`, and a build of `b489a48` for the before-numbers |
| 100 h quiet set, misfires (R4) | 6 | ~3 h | the TDVMS archive (not redistributable) |
| magnitude corpus (R3) | 4 | ~10 min to rebuild the manifest | `cascade_impl`, the FDSN pull |
| amplitude audit (R6.M1) | 2 | ~15 min | `cascade_impl`, KOERI station metadata |
| Table 6, MANT scan (R2) | 5 | **~5 h GPU** | the MANT archive, an 8 GB card, **three** scan processes (four exhaust it) |

The archives that cannot be redistributed are named in the report rather than
silently skipped. The KOERI datasets ship with the release, so a fresh machine
can check the KO claims without any private data.

## Adding a claim

Put the number in `bench/claims.json` with:

```json
{"id": "...", "where": "where it is published", "claim": "one sentence",
 "artifact": "repo:path/to/file", "probe": "csvcount:magnitude>=6.0",
 "expect": 413, "tolerance": 0,
 "repro": "the command that regenerates the artefact"}
```

`artifact` may name any of the four repositories (`ayzek_code`,
`ayzek_paper`, `cascade_impl`, `cnn_earthquake`); a bare path means this one.
Probe syntax is in `tools/verify_claims.py`. A claim with no `artifact` is
checked only by a full reproduction run.

**A number that moves must move in `claims.json` in the same commit.** Stage 0
will fail otherwise, which is the intended behaviour and not a nuisance.
