# Release notes: the method freeze (draft)

Draft written before the tag; the commit, tag and date are filled in at the
freeze. The held-out protocol (`docs/protocol-heldout.md`, section 1) fixes the
method by this release: every default of ayzek as in the tagged commit, the
models below by their SHA-256, the baselines' settings in
`tools/heldout_settings.json` and the published models' weight versions below.

- **Tag**: _at the freeze_
- **Commit**: _at the freeze_
- **Binary used by the tests**: `build-linux/app/ayzek`, the static x86-64
  build (zig/musl, `cross/x86_64-linux.ini`); its SHA-256 _at the freeze_.

## Models (`models/`)

| file | SHA-256 |
|---|---|
| `bandpass.ayzw` | `16e95eff0a22a2bea85d93e42a18cdfd037ae8382b70dadeb0c2288ecb4e279e` |
| `detector_s42.ayzw` | `9da4692c77d1fd250517a6ca070173e24346cd50193256429e76033bfb2242cb` |
| `detector_s43.ayzw` | `dfa6e7bedef17c0f3e20184209fb9911ad9301fb19a46a8126eaa222f34dd23e` |
| `detector_s44.ayzw` | `0887ce81b7a4126f8ecf694de61324ad1d054a5c604f93d57bb2aed1d904d18c` |
| `pd_gain.csv` | `75b59961836e0cb191be9eb270972ba201ca6e536065c67647104bb4f69b6f78` |
| `pd_noise.csv` | `5515181475301d73a8c21f27140243edd55b9d560e888a288c69bdf447b1d822` |
| `pd_relation.csv` | `52fab92cdbd49d41bc98dd0de4f4379e9a12869bdcd47706c7b01c3064b836e0` |
| `pd_terms.csv` | `cead7ab2406e86812f63a4072dbf61ca95003702ce7ff5de5d92dffe65ea8b1a` |
| `spicker.ayzw` | `f0ddbc33fb12c4fdfda2fcb04b8158396e6e431bb860404a2a49fc8f42282d5a` |
| `stations.csv` | `fd47434ade1ac1d9b2f7b1d2c12ee57dff0cf00c3cab2203baa55a0aa3fa7ae2` |
| `transformer.ayzw` | `d04cf8230a34b4381df1c0e21e54cdc503eb6c7750598acc146a87c0ffd940c7` |

Taken on 2026-10-30 at commit 5ed97f8; to be checked again at the tag
(`sha256sum models/*`).

## Published models (SeisBench 0.12.6)

`tools/seisbench_scores.py` pins each weight version (`VERSIONS`), since
`from_pretrained` would otherwise take the latest upload. The files, as
SeisBench caches them (`~/.seisbench/models/v3/`), are byte-identical on the
two machines that computed score caches (checked 2026-10-30). The weights
the frozen settings use are marked.

| model | weights | version | SHA-256 of the `.pt` file | frozen settings |
|---|---|---|---|---|
| EQTransformer | original | 3 | `c363c59dfb4562ba3b1356fb03dc9c776bbcfe3384b5db22e5ee3689de0a4df9` | look-ahead, equal detections and equal quiet alarms |
| EQTransformer | STEAD | 2 | `a7c0dead4bfc07df2e316e9230c28250f356841135ef13ccd4db33ed1a15c62e` | |
| EQTransformer | INSTANCE | 2 | `78f7b98b030019e352bcbdb75fa49755c704252747176a0e1f02df3bb39d6261` | causal (all three), look-ahead equal unmatched |
| EQTransformer | ETHZ | 2 | `3c5b7d687fd3572027719655ac3f10fa03e14f4250da51ada499a2d3b2b329b4` | |
| PhaseNet | STEAD | 2 | `761c46496139e9917a7480ff56f660914124c806511ea794f0005518358b09ec` | |
| PhaseNet | INSTANCE | 2 | `dddcc3d872106a8cd2013ce7a035c972a7f435749b5c02a015efd5ce3767da6d` | all six |
| PhaseNet | ETHZ | 2 | `89cb5ff3272ea91df7d6a3c347fdd067f5a6c2aab91aa884c448a93a482cd848` | |

PhaseNet's original weights were not a candidate (protocol, section 6).

## Changes since the development replays were scored

Merged into the release branch before the freeze, each checked on the
development data (`tools/regression_dev.py`):

- Crediting an event goes to the first alarm declared (lowest id), not the
  last. Three alarm delays among the 100 event rows changed by 0.2 to 0.3 s;
  no count changed, so no selection changed
  (`tests/reference/heldout_dev.csv`).
- `tools/network_subsets` scores as ayzek does: the full catalogue is
  matched, then the radius is taken from the full network's centroid.
- The analysis tools: fixes for empty quiet totals, the exact binomial at
  p = 1/2, reports without an "Alarms raised" line, and negative
  magnitudes and delays.
- The miniSEED decoder refuses dates past 2262 rather than overflowing.
- `tools/run_heldout.py`: the chance-corrected counts used the set's own
  catalogue instead of the shifted one (every count was 0); fixed before any
  test data exists.
- Live SeedLink input (`docs/impl/19-seedlink.md`), with Steim1 decoding.
  Replays are unchanged by it.
