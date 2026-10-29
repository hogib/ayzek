# Dry runs of the held-out selection

`heldout select` (ayzek_heldout 440a2d9) on periods before the test pools,
run on 2026-10-29 before the freeze, as the protocol (section 3) requires.
Each directory holds the `selection.json`, the set files of the selected
sequences and quiet windows, and the tool's log.

| run | period | outcome |
|---|---|---|
| `greece_dry2022` | 2022 | 17 qualifying sequences, 10 selected, 161 scorable events (target 150 reached); 24 quiet windows after 452 draws |
| `nz_dry2022` | 2022 | 23 qualifying sequences, 12 selected, 63 scorable events (under the target 150); 24 quiet windows after 453 draws |
| `ko_dry2026` | 2026-07-01 to 2026-08-30 | 1 qualifying sequence, none selected by its coin, 0 scorable events; no sets |

They test the tool, not the method: no waveform was fetched and nothing was
replayed. New Zealand's 2022 count shows that a year of its "now" sequences
may not reach the target; the protocol (section 3) then uses all of them and
reports the result as under the target.
