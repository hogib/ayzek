# ayzek

A real-time earthquake early-warning pipeline in C++23: detect, pick S, and
estimate magnitude on live continuous streams from three stations, with
hand-written inference and no ML runtime dependency.

Named after a bird.

```bash
meson setup build
meson compile -C build
meson test -C build
```

See [docs/DESIGN.md](docs/DESIGN.md) for the architecture and the staged plan.

**Status:** stage 0. Scaffold only.

## Shape

```
SeedLink ─→ decode ─→ [RING 1: raw @ 100 Hz]
                            │  continuous IIR, state carried across packets
                            ↓
                      [RING 2: filtered @ target rate]
                            │  mdspan views, zero copy
                            ↓
                      detector ─→ trigger ─→ picker + magnitude ─→ location
```

One ring pair per station, three stations -- the minimum that locates an event,
since differential arrival times cancel origin time and leave latitude and
longitude with depth fixed.

Models come from the sibling projects (`cnn_earthquake`, `sphase`) as exported
weights. Every kernel is checked against the Python it came from before it is
composed into anything.
