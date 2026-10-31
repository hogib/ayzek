# ayzek: command-line options

The complete option table of the `ayzek` binary, moved here from the README.
Defaults are those of release v1.0.0. The documents named in the right-hand
column are under `docs/impl/`.

| option | default | meaning |
|---|---|---|
| `--models DIR` | `models` | weights directory |
| `--stations FILE` | `MODELS/stations.csv` | station coordinates (`code,lat,lon[,...]`), for a replay of another network: the default table is keyed by station code alone and holds the TU and KO stations |
| `--pd-gains FILE` | `MODELS/pd_gain.csv` | the replay's own Pd gains (`station,channel,start,end,sensitivity`); KO's station terms and median noise levels are then not used, so terms are zero and the noise-band exclusions do not apply (`docs/impl/17-pd-magnitude.md`) |
| `--speed X` | 1 | replay speed; 1 = real time, 0 = as fast as possible |
| `--from TIME` | start of data | no triggers before this UTC time; earlier data still runs through the detectors' and the Pd chain's filters |
| `--threshold P` | 0.8 | trigger when `--trigger-windows` consecutive windows reach P |
| `--trigger-windows N` | 8 | see `--threshold`; 8 windows add 3.5 s to a trigger |
| `--instant-threshold P` | 0.9 | also trigger on a single window at or above P; a value > 1 disables this |
| `--release P` | 0.3 | a trigger re-arms after two windows below P |
| `--detector KIND` | `6s` | `6s` (also `model`): the 3-seed 6 s window detector. `transformer`: the streaming onset transformer, one output per 0.1 s, P dated by its own `dt` (`docs/impl/14-transformer.md`). `stalta`: a recursive STA/LTA. `filterpicker`: Lomax's FilterPicker5 on the vertical (`docs/impl/18-filterpicker.md`). Association and magnitude are shared; the transformer also places the picker on its own P and, with the geometry head, locates events itself (`--locate`) |
| `--transformer FILE` | `MODELS/transformer.ayzw` | transformer weights; its `--threshold` and `--release` default to the model's validation operating point |
| `--dt-reset BELOW,FROM` | the model's (`1,5` for the current one; `2,5` for models exported without one) | transformer: also trigger while p stays at or above the threshold when dt falls back to BELOW s or less (2 tokens in a row) after reaching FROM s: a new onset in the coda of the last event (`docs/impl/14-transformer.md`) |
| `--no-dt-reset` | | transformer: rising edges only, the rule before `--dt-reset` |
| `--geo-sd-scale X`, `--geo-sigma-p S`, `--geo-max-z Z`, `--geo-max-err-km KM` | 2.1, 2.75, off, off | geometry locator: widen every station's distance sd, the P-time uncertainty of its misfit, leave out or refuse stations whose distance disagrees with the solution, refuse a location with a 68% radius over KM (`docs/impl/15-geometry-location.md`) |
| `--sta S`, `--lta S` | 1, 30 | STA/LTA averaging lengths in seconds |
| `--stalta-on R`, `--stalta-off R` | 8, 1.5 | STA/LTA ratio that triggers, and below which a trigger re-arms |
| `--stalta-on-stations FILE` | | per-station STA/LTA trigger ratios (`station,ratio`), the per-station baseline; stations not listed keep `--stalta-on` |
| `--stalta-band LO,HI` | 2,20 | STA/LTA pass band in Hz (4th-order Butterworth high- and low-pass) |
| `--stalta-3c` | | STA/LTA on the energy of all three components instead of the vertical |
| `--fp-windows F,L,U` | 3,5,0.2 | FilterPicker filter window, long-term window and tUpEvent in seconds (Lomax et al. 2012 at 100 Hz; Earthworm's `pick_FP` uses 6,10,0.4) |
| `--fp-thresholds T1,T2` | 10,10 | FilterPicker trigger threshold and the mean characteristic function over tUpEvent needed to accept a pick |
| `--no-anchor` | | do not take the P time of a model trigger from the STA/LTA onset |
| `--anchor-on R` | 3 | STA/LTA ratio taken as the onset when anchoring |
| `--anchor-picker` | | also place the picker window from the anchored P; it makes the picker worse (`docs/impl/09-sta-lta.md`) |
| `--no-onset-fp` | | 6 s detector: do not date the P-time location and Pd by the FilterPicker pick within 3 s of the trigger's P (`docs/impl/17-pd-magnitude.md`, section 2.4) |
| `--no-p-locate`, `--p-sigma S` | , 1.35 | window detectors: no location from the detections' P times alone; its P-time uncertainty, which sets its 68% radius |
| `--require-onset` | | a detector trigger also needs an STA/LTA onset; measured, and worse than not requiring one (`docs/impl/09-sta-lta.md`) |
| `--step N` | 50 | samples between detector windows (0.5 s at 100 Hz) |
| `--min-stations N` | 2; transformer 3 | station detections needed to declare an event; 3 halves the unmatched alarms in a dense network but costs warning time (`docs/impl/05-pipeline.md`); the transformer, which also fires in codas, needs 3 (`docs/impl/14-transformer.md`) |
| `--slack S` | 3 | tolerance on the inter-station P travel time, seconds |
| `--locate KIND` | | `geometry`: locate from the transformer's geometry head, from the trigger on; the picker is not run. `picks`: from P and S picks. Default: `geometry` when `--detector transformer` and the model has the head, else `picks` |
| `--assess` | | diagnostic: judge every alarm earthquake, possible, misfire or unclassified from its S-P picks (runs the picker even with the geometry locator); `ASSESS` lines while running, a report section (`docs/impl/16-alarm-assessment.md`) |
| `--assess-csv FILE` | | with `--assess`: one row per alarm, readable by `tools/plot_candidates.py` |
| `--no-pick` | | no P/S picker, and so no location unless `--locate geometry` |
| `--no-magnitude` | | no magnitude estimate, as `--no-pd-magnitude`; with `--from`, a detector that needs no earlier data (the transformer) also starts reading there |
| `--no-pd-magnitude` | | no magnitude from peak P displacement. On by default, for every detector, when `models/pd_relation.csv` exists and the station has a gain in `models/pd_gain.csv` (`tools/export_pd.py`, `docs/impl/17-pd-magnitude.md`) |
| `--catalog CSV` | | AFAD catalogue export; each alarm is compared with it |
| `--site NAME,LAT,LON` | | also report the S-wave warning time at this place (repeatable) |
| `--scores DIR` | | write every window's detector probability to `DIR/STATION.csv` |
| `--scores-in DIR` | | use the probabilities of an earlier `--scores` run instead of running the detector |
| `--record FILE` | | write all station outputs, for `tools/network_subsets` |
| `--verbose` | | print every station detection, pick and Pd value |
| `--no-color` | | |
| `--seedlink HOST:PORT` | | live input from this SeedLink 3 server (port 18000 if omitted), in place of the files; not in the Windows build (`docs/impl/19-seedlink.md`) |
| `--streams LIST` | | with `--seedlink`: the stations, `NET.STA` or `NET.STA.LOC` (`--` for the empty location), comma-separated; each needs coordinates in the station table, and its HHZ, HHN and HHE at 100 Hz are used |
| `--archive DIR` | | with `--seedlink`: append every record received, as received, to `DIR/STATION.mseed`, which replays like any other file |
| `--duration S` | until SIGINT or SIGTERM | with `--seedlink`: stop the input after S seconds; the first SIGINT or SIGTERM stops the input and the run ends with the report once what arrived is processed, a second one stops at once |
| `--live-wait S` | 10 | with `--seedlink`: about how long a station's message waits for the other stations' data before the network stage takes it; a station that holds one back that long is not waited for until it catches up. With nothing out of order (the report counts it), a replay of the archive gives the live run's alarms |
| `--seedlink-timeout S` | 120 | with `--seedlink`: reconnect after S seconds with nothing received; reconnection is retried after 1, 2, 4, ... up to 30 s, and every station resumes after its last packet |
| `--seedlink-keepalive S` | 60 | with `--seedlink`: send `INFO ID` after S seconds with nothing received (0: never); IRIS asks for no more than one every 4 minutes on an idle connection |


## Running on Windows

`ayzek.exe` is built on Linux with zig and needs no runtime DLLs beyond what
Windows 10 and 11 already have (`docs/impl/11-windows.md`):

```bash
meson setup build-win --cross-file cross/x86_64-windows.ini --buildtype=release
ninja -C build-win                       # -> build-win/app/ayzek.exe, 1.6 MB
```

To run it, put `ayzek.exe` in a folder with the `models` folder beside it, and
the waveforms wherever you like:

```
C:\ayzek\
    ayzek.exe
    models\          the whole folder from this repository
    data\demo\       one .mseed per station
    demo.csv         optional AFAD catalogue export (tests\catalogs\demo.csv)
```

In PowerShell or Command Prompt, from `C:\ayzek`:

```
ayzek.exe --speed 10 --from 2025-11-10T18:20:00 data\demo\DEMI.mseed data\demo\MANT.mseed data\demo\BAND.mseed
```

- `--speed 10` replays at ten times real time; `--speed 0` runs as fast as the
  machine allows, `--speed 1` in real time.
- `--catalog demo.csv` adds the comparison with the AFAD catalogue. Without it,
  alarms are still reported, just not scored.
- `--models DIR` points at the weights if they are not in `.\models`.
- Colour is disabled automatically when the output is redirected to a file.
  The colours are ANSI escapes, which Windows Terminal renders and the legacy
  console may print as `[32m` instead; `--no-color` turns them off.
- Everything else is the same as on Linux, including the option table below,
  except live input: `--seedlink` needs POSIX sockets, which the Windows build
  does not have, and says so.

**Not yet run on Windows.** The binaries are built and linked, and the
numerical tests have not been executed on a Windows machine or under wine. If
you run them, `meson test -C build-win` needs wine on Linux; on Windows the
test executables in `build-win\tests\` can be run directly, and `test_models`,
`test_magnitude`, `test_dsp` and `test_stalta` are the ones that would catch a
floating-point difference.

