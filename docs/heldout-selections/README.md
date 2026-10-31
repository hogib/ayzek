# The held-out selections (after the freeze)

`heldout select` on the three pools of the protocol (section 3), run after
the freeze tag v1.0.0 (commit 3e43064, pushed 2026-10-30 12:55:36 UTC, the pool
end). Each directory holds `selection.json`, the set files of the selected
sequences and quiet windows, and the tool's log. No waveform was fetched
before these selections existed; the fetch of each followed its selection.

| pool | period (UTC) | usable stations | catalogue events | qualifying / now / selected sequences | scorable events (target) | quiet windows |
|---|---|---|---|---|---|---|
| `ko_now` | 2026-08-30 to 2026-10-30T12:55:36 | 123 of 155 | 5,693 (EMSC) | 3 / 2 / 2 | **36 (240), under the target** | **14 of 24**, after the relaxed pass |
| `nz_now` | 2023-01-01 to 2026-10-30T12:55:36 | 36 of 134 | 81,990 (GeoNet) | 29 / 16 / 13 | 154 (150) | 24 after 551 draws |
| `greece_now` | 2023-01-01 to 2026-10-30T12:55:36 | 81 of 200 | 90,382 (NOA) | 84 / 42 / 11 | 178 (150) | 24 after 2,059 draws |

KO's two-month pool holds too few M >= 4.5 events near five usable stations
to reach its target; the protocol uses all of its "now" sequences and reports
the result as under the target. Its quiet windows stop at 14 because 10,000
draws in each pass found no more windows free of catalogue events (KO's
catalogue is EMSC's, dense in small events; `selection.json`, `quiet`).

## The fetch tool's commits

The release names the fetch tool at commit `440a2d9`. The selections ran as
follows (the tool stamps its commit when it writes, from the repository's
HEAD at that moment, so a stamp can postdate the code that ran):

- **KO**: run and fetched entirely at `440a2d9`.
- **New Zealand**: the process started at `440a2d9` (2026-10-30 13:00 UTC)
  and ran that code to the end; its `selection.json` is stamped `6fda3e2`
  because HEAD had moved by the time it wrote (23:03 UTC). None of the
  changes below is on its path: GeoNet has no availability service, so the
  listing by years and the one-at-a-time requests never ran for it, and no
  request of its exhausted the earlier retry waits. Its sets were fetched at
  `6fda3e2`.
- **Greece**: run and fetched at `6fda3e2`, after three changes to how the tool
  talks to NOA's node, each recorded in the tool's README ("Changes after the
  freeze") and none to what is selected or written:
  1. `658ae03`: fdsnws/availability asked one year at a time. Four runs died
     in the availability stage because a station's four-year listing never
     came back within the 300 s read timeout (`failed-attempts/`; a survey of
     all 200 stations found five such stations).
  2. `5d088c5`: retries wait 33 minutes in all per request (was 8), and the
     probe's availability requests go one at a time (was two): with yearly
     pieces, NOA still answered 502 or nothing to both threads for stretches
     longer than 8 minutes, about an hour into each run.
  3. `6fda3e2`: the retry messages name the station and piece, and a year's
     listing the service never answers is asked for in halves down to a week.
     The successful run needed no halving: the stubborn station was HP.KSTR,
     whose 2024 and 2025 listings came back on the sixth and seventh retries.
  During the successful run's catalogue stage, NOA's event service returned
  503 to everything for a few minutes (2026-10-31 02:40 UTC); the process was
  paused with SIGSTOP at 02:43, so that the retries did not run out, and
  resumed at 02:47 when the service answered again. The pause changes no
  request and no answer.

The rules and seeds are the protocol's; the pool end is the UTC time of the
protocol's public commit (section 10). The stopping rule's counts are in
`selection.json` (`scorable`, `cumulative`).
