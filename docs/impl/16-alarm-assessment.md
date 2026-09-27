# 16 · Alarm assessment (`--assess`)

An alarm the catalogue does not list is not necessarily false: the catalogue
misses small events, and aftershock sequences most of all. `--assess` judges
every declared event from its own waveforms' picks, so the unmatched alarms
can be split into probable uncatalogued earthquakes and real misfires. It is a
diagnostic, for replays and benchmarks. It runs the P/S picker at every
trigger, also with the geometry locator, where the picks are then used for the
assessment only, and after the replay it picks at the stations that did not
trigger (below).

```bash
build-release/app/ayzek --speed 0 --detector transformer --assess \
    --assess-csv alarms.csv --record run.rec \
    --catalog tests/catalogs/marmara_ko.csv data/marmara_ko/*.mseed
```

## The test

**S-P.** A station's S-P interval gives its distance to the source,
d = (tS − tP)·vp·vs/(vp − vs), about 8.4 km per second at vp 6 and vs 3.5 km/s:
a circle around the station. With confident P and S picks (the pick QC of the
network stage) at two or more stations, the event is consistent when

- every pair's circles can meet: |d_a − d_b| ≤ separation + 5 km and
  separation ≤ d_a + d_b + 5 km;
- one source fits every P and S time (the picks locator's grid search) with an
  rms of at most 1.5 s (`assess_max_rms`);
- no station's P or S is more than 2 s off that source (`assess_station_tol`).

With three or more stations, one station may fail and be left out, since one
bad pick should not condemn an event the others agree on. Noise that happens
to trigger two stations together does not produce a confident S at each, at
S-P intervals whose circles meet, with P times that match the distances.

**Support and objections**, each listed as a reason:

| evidence | for | against |
|---|---|---|
| the transformer's distances (geometry head) against S-P | median within 2 sd, epicentres within 30 km | median beyond 4 sd (not counted when S-P puts a station beyond 60 km, outside the head's training range) |
| station magnitudes | sd ≤ 0.5 | sd > 1 |
| chance coincidence: the other declaring stations each triggering inside their association window, at their own trigger rate in the run | p ≤ 0.01 | p > 0.1 |

**Verdicts:**

| verdict | when |
|---|---|
| `earthquake` | S-P consistent at 3+ stations, or at 2 with support and nothing against |
| `possible` | S-P consistent at 2 stations, without support or with an objection |
| `misfire` | S-P picks at 3+ stations that no single source explains, even leaving one out; or at 2 with an objection besides (two stations alone cannot say which pick is wrong) |
| `unclassified` | fewer than 2 stations with confident P and S |

**Context**, with `--catalog` (loaded 30 days back for this): an unmatched
alarm whose source lies within 100 km and 60 s of a catalogued event is
probably that event, missed by the matching; one within 50 km of an M ≥ 3
catalogued in the 30 days before is in its aftershock zone.

## Picks at stations that did not trigger

Most alarms are declared by exactly two stations, so one failed pick leaves
nothing to test. After the replay, for every alarm with a source estimate
(the S-P source if its picks agree, else its location), the network lists the
stations within 200 km that have no pick for it (`Network::pick_requests`),
with their predicted P. The main thread cuts that 60 s window from the
station's replay file (`read_window`), and the picker looks for P within 4 s
of the prediction and for S only before the station's next detected onset.
The picks pass the usual QC around the predicted P, join the assessment
(marked "not triggered" in the reasons) and go into the recording with the
event they were made for. They never change an event's location.

This is done after the replay and on one thread, so it does not depend on
thread timing. A live client would need to keep the data for it instead.

## The picker window

For the transformer the picker looks for P within 3 s of the transformer's
own P and for S only before the station's next trigger (`Picker::Search`,
`14-transformer.md`; `--pick-anywhere` turns this off). Unconstrained, in an
aftershock sequence it often picked a later event's P or S: 39 of marmara's
201 picks had P outside the trigger window, and 26 an S not after P.

## Output

- **While running**, an `ASSESS` line whenever an event's verdict or its
  number of picked stations changes, i.e. as the picks arrive, about a minute
  after the alarm. Trigger rates are those so far.
- **Report**: an `assessed` line in every event block, with its reasons and
  context, and a summary of verdicts for catalogued alarms (how often a real
  event passes) and for unmatched ones.
- `--assess-csv FILE`: one row per alarm. It carries `record`, `time_utc`,
  `station_a`, `station_b` and `dt_s` as `tools/misfire_review.py` does, so
  `tools/plot_candidates.py --csv FILE --records DIR --data DIR` draws the
  candidates (the record's name without `.rec` is the waveform folder under
  `--data`).
- `tools/network_subsets --assess` re-assesses a recording without replaying
  (the recording holds the picks), and `tools/scorecard.py` and
  `tools/compare_detectors.py` report the counts when ayzek ran with
  `--args '--assess'`.

## Relation to tools/misfire_review.py

That script tests pairs of stations offline, from recordings, with the same
circle and moveout test and a time-slide null. `--assess` does it inside the
network stage, for any number of stations, and adds the geometry head and the
station magnitudes as independent evidence. The time slide is not ported: the
analytic chance probability is.

## Results

The v2-geo transformer (geo_v1) on the public KO sets, `--assess`:

| set | catalogued alarms: real (earthquake or possible) / misfire / unclassified | unmatched alarms: real / misfire / unclassified |
|---|---|---|
| demo_ko | 4 of 6 / 0 / 2 | 0 / 0 / 7 |
| marmara_ko | 14 of 19 / 1 / 4 | 3 / 4 / 30 |
| quiet_ko | – | 0 / 0 / 1 |

How it got there, on marmara's 19 catalogued alarms:

| version | real | misfire |
|---|---:|---:|
| first: rms 1.5 s, 2 s per station, no S-P cap, two disagreeing stations a misfire | 8 | 6 |
| rms 2.5 s, 3.5 s per station, S-P beyond 200 km set aside, two stations a misfire only with an objection | 11 | 3 |
| + the picker placed on the transformer's P, and picks at stations that did not trigger | 14 | 1 |

The misfires of the first version came from one implausible S-P (38 s, "321
km": the picker's window held a later event) or from fits just over 1.5 s in a
uniform half-space at 50–100 km.

**Most unmatched alarms are still unclassified.** The test needs confident P
and S at two stations, and the picker is rarely confident on small events near
their true P. Of marmara's 201 picks after triggers, 82 pass the pick QC; most
of the rest have a P probability below 0.5 (113). Of the 243 picks made at
stations that did not trigger, 42 pass. Placing the picker on the right P
removed the confident picks on later events, which had passed QC for the wrong
event. What would classify more: a picker trained on small events, or an S
pick from the transformer itself.
