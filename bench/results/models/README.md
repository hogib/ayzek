# Benchmarks of every model, 21 October 2026

`ayzek_bench` (stages, layers, capacity) and `tools/scorecard.py` (demo_ko and
quiet_ko) for each detector configuration, at 7748168 plus the `ayzek_bench`
change of 68e0c47, the code of v0.1.0-beta.3. Summary and tables:
`docs/impl/10-benchmarks.md`; the paper's log: `ayzek_paper/revisions2.md` R7.

| path | what |
|---|---|
| `x86_{stages,layers,capacity}.*` | AMD Ryzen 5 5600, `build-release` (release, scalar backend), idle machine; the load average at the start of each group is in `x86_load.txt` |
| `x86_scorecard_*.json` | the same machine: `6s` (the default detector), `transformer` (Pd magnitude), `transformer_cnnmag` (`--no-pd-magnitude`), `stalta_eqfa` and `stalta_eqdet` (the two STA/LTA settings of docs/impl/09) |
| `pi_cool/` | Raspberry Pi 4 Model B (Cortex-A72, 1.8 GHz, NEON), static `build-pi`. Each speed group started below 55 °C; `thermal.csv` is time, temperature (°C), ARM clock (Hz) and throttle flags every second, `marks.txt` the start and end of each group. The clock held 1.8 GHz throughout. These are the Pi numbers the docs and the paper quote. |
| `pi_first/` | the first Pi run (`tools/bench_pi.sh`), without cooling or a thermal log; the board reached 84 °C during its accuracy pass. Its speed agrees with `pi_cool` (36.8 against 36.3 ms per window); kept for the scorecard, which agrees with x86 field for field. |

The Pi was borrowed and has been returned, so these cannot be remeasured on it.
