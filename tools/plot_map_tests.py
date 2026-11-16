"""Maps of the three held-out test regions, one per network group: Greece (HL,
HT, HA, HP, HC, coloured by network), New Zealand (NZ) and Türkiye (KO). A
station is usable when the held-out repository's pre-pool tuning sets found
three HH channels at 100 Hz with at least 90% coverage: it is drawn filled
when any of a group's four sets chose it, and open when it was probed but
never chosen, so the map shows the network's footprint. No event or sequence
is drawn: the test sequences are not selected yet.

    uv run --project tools --with pygmt --with xarray --with netcdf4 \\
        python tools/plot_map_tests.py --heldout ../ayzek_heldout \\
        --out ../ayzek_paper/figures/fig_map

Stations: HELDOUT/data/tune_<group>_{1..4}/stations.json ("probed", with
"ok", and "chosen"). Where HELDOUT/selections/<name>/selection.json exists for
a group (--selection, default greece_dry2022, nz_dry2022 and ko_dry2026), its
stations.usable and stations.rejected lists are drawn instead and the CSV's
provenance column says so; nothing but the stations is read from a selection.
Relief, coastlines and borders as in plot_map_development.py: ETOPO 2022 cut
over OPeNDAP once into data/maps/, GSHHG and DCW from --gshhg and --dcw.
Writes OUT_<group>.pdf, OUT_<group>.png (300 dpi) and OUT_<group>.csv (every
station drawn, with its role, the sets that chose or probed it, and its
provenance). GeoNet's Chatham Islands station CTZ (176.6 W, east of the
antimeridian) is drawn in an inset of the New Zealand map, and the Tongariro
and Ruapehu volcano stations, too close to label at the map's scale, are
labelled in a second, zoomed inset.

With --tests the maps show the tests as they were run, from the final
selections (HELDOUT/selections/<group>_now/selection.json, --tests-selection):
the stations of the selected sequences filled and labelled, the other usable
stations as small open grey symbols, a ring around every station that takes
part in one of the protocol's quiet windows (its pen in three classes of the
number of windows), and a star with the sequence
number at the first event of each selected sequence. The CSV then lists the
stations (role, sequences, number of quiet windows) and the first events.

    uv run --project tools --with pygmt --with xarray --with netcdf4 \\
        python tools/plot_map_tests.py --heldout ../ayzek_heldout --tests \\
        --out ../ayzek_paper/figures/fig_map
"""
import argparse
import csv
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ETOPO = ("https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/30s/"
         "30s_surface_elev_netcdf/ETOPO_2022_v1_30s_N90W180_surface.nc")

# Symbol (GMT -S code), fill and legend name of each network. The five Greek
# networks take Okabe-Ito colours, which stay apart under colour-vision
# deficiency, and a shape each, so a network is not told by colour alone; KO
# keeps its symbol from the development map, and NZ shares it.
NETWORKS = {
    "HL": ("i", "0/114/178", "HL, NOA"),
    "HT": ("s", "213/94/0", "HT, AUTH"),
    "HA": ("c", "0/158/115", "HA, NKUA"),
    "HP": ("d", "204/121/167", "HP, Univ. of Patras"),
    "HC": ("t", "230/159/0", "HC, Crete"),
    "NZ": ("i", "40/110/190", "NZ, GeoNet"),
    "KO": ("i", "40/110/190", "KO, KOERI"),
}

# The regions are a little wider than the networks' own areas, so that every
# usable station lies inside: HL.OTHN (19.4 E) and HL.KSL (29.6 E) for Greece,
# KO's four Cyprus stations (35.1-35.3 N) for Türkiye.
GROUPS = {
    "greece": dict(networks=["HL", "HT", "HA", "HP", "HC"], selection="greece_dry2022",
                   region=[19.0, 30.0, 34.0, 42.0], projection="M12c", frame=["xa2f1", "ya1f0.5"],
                   symbol=0.2, labels=False),
    "nz": dict(networks=["NZ"], selection="nz_dry2022",
               region=[165.0, 180.0, -48.0, -34.0], projection="M12c", frame=["xa5f1", "ya2f1"],
               symbol=0.22, labels=True,
               insets=[
                   # Off the map: GeoNet's Chatham Islands station, east of the antimeridian.
                   dict(name="chatham", kind="off-map", title="Chatham Islands",
                        region=[-177.5, -175.7, -44.5, -43.5], projection="M2.6c",
                        position="jBR+w2.6c/2.2c+o0.3c/0.5c"),
                   # A zoom: the nine Tongariro and Ruapehu volcano stations lie within
                   # 20 km of each other and can be labelled only enlarged.
                   dict(name="tongariro", kind="zoom", title="Tongariro and Ruapehu",
                        region=[175.42, 175.80, -39.36, -39.05], projection="M5.5c",
                        position="jTL+w5.5c/5.8c+o0.3c/1.2c"),
               ],
               # --tests: the North Island stations and first events lie too close to
               # label at the map's scale; the inset sits over the Tasman Sea.
               # The Chatham Islands inset is left out: its one station, CTZ, is in no
               # selected sequence (it stays in the CSV, placed "outside the frame").
               tests_insets=[
                   dict(name="north_island", kind="zoom", title="North Island",
                        region=[174.3, 178.6, -39.5, -36.3], projection="M5.2c",
                        position="jTL+w5.2c/4.76c+o0.3c/0.3c"),
               ]),
    "ko": dict(networks=["KO"], selection="ko_dry2026",
               region=[25.0, 45.0, 34.6, 42.5], projection="M17c", frame=["xa5f1", "ya2f1"],
               symbol=0.2, labels=False),
}

# Station labels moved off their neighbours (dx/dy); the rest sit below their
# symbol. WIZ and WSRZ share White Island; the volcano stations are labelled
# in the Tongariro inset.
LABEL_OFFSET = {
    "WIZ": "-0.3c/0.05c", "WSRZ": "0.35c/0.05c", "WVZ": "0.32c/0.05c", "WHSZ": "0.12c/-0.24c",
    "URZ": "0.12c/0.24c", "MWZ": "-0.05c/-0.25c", "RTZ": "0c/-0.24c", "KHEZ": "-0.4c/0c",
    "RATZ": "0.42c/0c",
    "COVZ": "-0.42c/0c", "FWVZ": "-0.42c/0c", "MAVZ": "0.42c/0c", "WHVZ": "0.44c/0c",
    "TRVZ": "0c/-0.26c", "OTVZ": "0c/-0.26c", "ETVZ": "0.42c/0c", "TMVZ": "-0.42c/0.05c",
    "NTVZ": "0c/0.26c",
}


# --tests: the stations of the selected sequences that are labelled (None: all
# of them), and labels and sequence numbers moved off their neighbours.
USED, OTHER = "in a selected sequence", "usable, in no selected sequence"
TESTS_LABELS = {"greece": set(), "nz": None, "ko": None}
TESTS_LABEL_OFFSET = {
    ("nz", "HIZ"): "-0.45c/0c", ("nz", "VRZ"): "0c/-0.3c", ("nz", "KPRZ"): "0c/0.33c",
    ("nz", "RATZ"): "0c/-0.36c", ("nz", "TLZ"): "-0.45c/0c", ("nz", "URZ"): "-0.45c/0c",
    ("nz", "MWZ"): "0.42c/0c", ("nz", "WVZ"): "0c/0.33c", ("nz", "WHSZ"): "-0.5c/0.15c",
    ("ko", "AKDG"): "0c/0.33c",
}
TESTS_EVENT_OFFSET = {
    ("nz", "008"): "-0.05c/0.3c", ("nz", "017"): "0.32c/-0.15c",
    ("nz", "009"): "-0.36c/-0.1c", ("nz", "013"): "0.34c/0c",
}
STAR = "a0.28c"
# Ring pen by the number of the protocol's quiet windows a station takes part
# in: (upper limit, pen width, legend text). The classes follow the counts of
# the three final selections (1 and 3; 7, 10 and 14; 17 to 23).
RING_CLASSES = [(3, "0.4p", "1 to 3 quiet windows"), (14, "1p", "4 to 14 quiet windows"),
                (10 ** 9, "1.9p", "15 or more quiet windows")]


def ring_class(n):
    return next(i for i, (hi, _, _) in enumerate(RING_CLASSES) if n <= hi)


def from_tests(path, name, networks):
    """The usable stations of a final selection, each with the selected
    sequences that use it and the number of quiet windows it takes part in,
    and the first event of every selected sequence."""
    sel = json.loads(path.read_text())
    seqs = [q for q in sel["sequences"] if q.get("selected")]
    src = f"selections/{name}/selection.json"
    rows = []
    for s in sel["stations"]["usable"]:
        code = f"{s['network']}.{s['station']}"
        of = [q["index"] for q in seqs if code in q["stations"]]
        nq = sum(1 for w in sel["quiet"]["windows"] if code in w["stations"])
        rows.append(dict(role=USED if of else OTHER, network=s["network"], station=s["station"],
                         lat=s["lat"], lon=s["lon"], sequences=" ".join(f"{i:03d}" for i in of),
                         quiet_windows=nq, magnitude="", depth_km="", time="",
                         provenance=f"{src}: stations.usable, sequences, quiet.windows"))
    known = {f"{r['network']}.{r['station']}" for r in rows}
    for q in seqs:
        missing = set(q["stations"]) - known
        if missing:
            raise SystemExit(f"{src}: sequence {q['index']} uses {sorted(missing)}, not in stations.usable")
    if any(r["network"] not in networks for r in rows):
        raise SystemExit(f"{src}: a network outside {networks}")
    events = [dict(role="first event of a selected sequence", network="", station=f"{q['index']:03d}",
                   lat=q["first"]["lat"], lon=q["first"]["lon"], sequences=f"{q['index']:03d}",
                   quiet_windows="", magnitude=round(q["first"]["magnitude"], 1),
                   depth_km=round(q["first"]["depth_km"]), time=q["first"]["time"][:19],
                   provenance=f"{src}: sequences[].first") for q in seqs]
    return rows, events


def from_tune_sets(heldout, group):
    """Every station probed by the group's four tuning sets: usable when any
    set chose it, else probed but not usable."""
    seen = {}
    for i in range(1, 5):
        d = json.loads((heldout / "data" / f"tune_{group}_{i}" / "stations.json").read_text())
        chosen = {(s["network"], s["station"]) for s in d["chosen"]}
        for s in d["probed"] + [s for s in d["chosen"]
                                if (s["network"], s["station"]) not in
                                {(p["network"], p["station"]) for p in d["probed"]}]:
            k = (s["network"], s["station"])
            r = seen.setdefault(k, dict(network=k[0], station=k[1], lat=s["lat"], lon=s["lon"],
                                        elevation_m=s.get("elevation_m"), chosen=[], probed=[],
                                        coverage_min=None, reason=""))
            r["probed"].append(i)
            if k in chosen:
                r["chosen"].append(i)
            c = s.get("coverage_min")
            if c is not None and (r["coverage_min"] is None or c > r["coverage_min"]):
                r["coverage_min"], r["reason"] = c, s.get("reason", "")
    src = f"data/tune_{group}_{{1..4}}/stations.json"
    rows = []
    for r in seen.values():
        usable = bool(r["chosen"])
        rows.append(dict(role="usable" if usable else "probed, not usable",
                         network=r["network"], station=r["station"], lat=r["lat"], lon=r["lon"],
                         elevation_m=r["elevation_m"],
                         sets_chosen=" ".join(map(str, r["chosen"])),
                         sets_probed=" ".join(map(str, r["probed"])),
                         coverage_min=r["coverage_min"], reason="" if usable else r["reason"],
                         provenance=f"{src}: {'chosen' if usable else 'probed, not chosen'}"))
    return rows


def from_selection(path, name, probed):
    """The usable and rejected stations of a selection; its sequences and
    catalogue are left unread. A rejected station carries no coordinates:
    they come from the tuning sets' probed list, and one probed by no set is
    left out."""
    st = json.loads(path.read_text())["stations"]
    known = {(r["network"], r["station"]): r for r in probed}
    rows = []
    for role, key in (("usable", "usable"), ("probed, not usable", "rejected")):
        for s in st.get(key, []):
            k = (s["network"], s["station"])
            if "lat" not in s and k not in known:
                print(f"   {k[0]}.{k[1]}: {key} in the selection, coordinates unknown, left out")
                continue
            p = known.get(k, {})
            rows.append(dict(role=role, network=k[0], station=k[1],
                             lat=s.get("lat", p.get("lat")), lon=s.get("lon", p.get("lon")),
                             elevation_m=s.get("elevation_m", p.get("elevation_m")),
                             sets_chosen="", sets_probed="", coverage_min=s.get("coverage_min"),
                             reason=s.get("reason", ""),
                             provenance=f"selections/{name}/selection.json: stations.{key}"))
    return rows


def inside(lat, lon, region):
    lo, hi, la, lb = region
    return la <= lat <= lb and any(lo <= x <= hi for x in (lon, lon - 360, lon + 360))


def relief(path, region):
    """ETOPO 2022 over `region`, cut over OPeNDAP the first time."""
    import xarray as xr
    if path.exists():
        with xr.open_dataset(path) as ds:
            lon, lat = ds["lon"].values, ds["lat"].values
        if (lon.min() <= region[0] + 0.01 and lon.max() >= region[1] - 0.01
                and lat.min() <= region[2] + 0.01 and lat.max() >= region[3] - 0.01):
            return path
    path.parent.mkdir(parents=True, exist_ok=True)
    z = xr.open_dataset(ETOPO)["z"].sel(lat=slice(region[2], region[3]),
                                        lon=slice(region[0], region[1])).load()
    z = z.astype("float32")
    z.attrs = {"units": "m", "long_name": "elevation (ETOPO 2022, 30 arc-seconds)"}
    z.to_dataset(name="z").to_netcdf(path, format="NETCDF3_64BIT")
    return path


def legend_spec(path, rows, g, events=None):
    """A legend file: one filled entry per network with its usable count, and
    one open entry for the stations probed but not usable. With events
    (--tests): the stations of the selected sequences per network, the other
    usable stations, the stations of the quiet windows and the first events."""
    size = g["symbol"]
    lines = []
    if events is not None:
        for net in g["networks"]:
            sym, fill, name = NETWORKS[net]
            n = sum(1 for r in rows if r["network"] == net and r["role"] == USED)
            lines.append(f"S 0.3c {sym} {size}c {fill} 0.3p,black 0.65c {name} ({n})")
        sym = NETWORKS[g["networks"][0]][0]
        n = sum(1 for r in rows if r["role"] == OTHER and r["drawn"])   # the symbols drawn
        lines.append(f"S 0.3c {sym} {size * 0.6:.2f}c white 0.4p,gray35 0.65c usable, in no sequence ({n})")
        for i, (_, pen, text) in enumerate(RING_CLASSES):
            n = sum(1 for r in rows if r["quiet_windows"] and r["drawn"] and ring_class(r["quiet_windows"]) == i)
            if n:
                lines.append(f"S 0.3c c {size * 1.75:.2f}c - {pen},black 0.65c in {text} ({n})")
        lines.append(f"S 0.3c a 0.3c 235/200/40 0.4p,black 0.65c first event, with its sequence ({len(events)})")
        path.write_text("\n".join(lines) + "\n")
        return path
    for net in g["networks"]:
        sym, fill, name = NETWORKS[net]
        n = sum(1 for r in rows if r["network"] == net and r["role"] == "usable")
        lines.append(f"S 0.3c {sym} {size}c {fill} 0.3p,black 0.65c {name} ({n})")
    n = sum(1 for r in rows if r["role"] != "usable")
    if n:
        sym, fill, _ = NETWORKS[g["networks"][0]]
        # An open symbol keeps its network's shape; the entry shows the first
        # network's, and stays short so the box covers no station.
        lines.append(f"S 0.3c {sym} {size * 0.85:.2f}c white 0.5p,{fill} 0.65c probed, not usable ({n})")
    path.write_text("\n".join(lines) + "\n")
    return path


def draw_tests(fig, rows, events, g, group, placed, size):
    """--tests: the other usable stations small, open and grey; the rings of
    the quiet-window stations; the stations of the selected sequences filled;
    the first events as stars; then the labels and the sequence numbers."""
    for net in g["networks"]:
        sym, fill, _ = NETWORKS[net]
        sel = [r for r in rows if r["network"] == net and r["role"] == OTHER and placed in r["drawn"]]
        if sel:
            fig.plot(x=[r["lon"] for r in sel], y=[r["lat"] for r in sel],
                     style=f"{sym}{size * 0.6:.2f}c", fill="white", pen="0.4p,gray35")
    for i in reversed(range(len(RING_CLASSES))):        # the thick rings first, the thin over them
        ring = [r for r in rows if r["quiet_windows"] and placed in r["drawn"]
                and ring_class(r["quiet_windows"]) == i]
        if ring:
            fig.plot(x=[r["lon"] for r in ring], y=[r["lat"] for r in ring],
                     style=f"c{size * 1.75:.2f}c", pen=f"{RING_CLASSES[i][1]},black")
    for net in g["networks"]:
        sym, fill, _ = NETWORKS[net]
        sel = [r for r in rows if r["network"] == net and r["role"] == USED and placed in r["drawn"]]
        if sel:
            fig.plot(x=[r["lon"] for r in sel], y=[r["lat"] for r in sel],
                     style=f"{sym}{size}c", fill=fill, pen="0.3p,black")
    ev = [e for e in events if placed in e["drawn"]]
    if ev:
        fig.plot(x=[e["lon"] for e in ev], y=[e["lat"] for e in ev], style=STAR,
                 fill="235/200/40", pen="0.4p,black")
    only = TESTS_LABELS[group]
    for r in rows:
        if r["labelled"] == placed and r["role"] == USED and (only is None or r["station"] in only):
            fig.text(x=r["lon"], y=r["lat"], text=r["station"], font="5p,Helvetica,40/80/150",
                     offset=TESTS_LABEL_OFFSET.get((group, r["station"]),
                                                   "0c/-0.33c"))
    for e in ev:
        if e["labelled"] == placed:
            fig.text(x=e["lon"], y=e["lat"], text=e["station"], font="6p,Helvetica-Bold,black",
                     offset=TESTS_EVENT_OFFSET.get((group, e["station"]), "0.28c/0.12c"),
                     fill="white@30", clearance="0.03c/0.03c")


def draw_stations(fig, pygmt, rows, g, placed, size, labels):
    """Open symbols first, then the filled ones over them (HA.MRTH is usable
    where HL.MRTH is not), then the labels."""
    for role in ("probed, not usable", "usable"):
        for net in g["networks"]:
            sel = [r for r in rows if r["network"] == net and r["role"] == role and placed in r["drawn"]]
            if not sel:
                continue
            sym, fill, _ = NETWORKS[net]
            if role == "usable":
                fig.plot(x=[r["lon"] for r in sel], y=[r["lat"] for r in sel],
                         style=f"{sym}{size}c", fill=fill, pen="0.3p,black")
            else:
                fig.plot(x=[r["lon"] for r in sel], y=[r["lat"] for r in sel],
                         style=f"{sym}{size * 0.85:.2f}c", fill="white", pen=f"0.5p,{fill}")
    if labels:
        for r in rows:
            if r["labelled"] == placed and r["role"] == "usable":
                fig.text(x=r["lon"], y=r["lat"], text=r["station"], font="5p,Helvetica,40/80/150",
                         offset=LABEL_OFFSET.get(r["station"], "0c/-0.24c"))


def plot(group, rows, out, a, pygmt, events=None):
    """Each station is drawn on the map, in an off-map inset, or nowhere
    (outside the frame); one inside a zoom inset is drawn there too, and
    labelled only there. With events (--tests) the first events are placed
    the same way and the tests' roles are drawn."""
    g = GROUPS[group]
    insets = g.get("insets", []) if events is None else g.get("tests_insets", g.get("insets", []))
    for r in rows + (events or []):
        r["drawn"], r["labelled"], r["placed"] = [], None, "outside the frame"
        for ins in insets:
            if inside(r["lat"], r["lon"], ins["region"]):
                r["drawn"].append(ins["name"])
                r["labelled"] = ins["name"]
                r["placed"] = f"{ins['title']} inset"
                if ins["kind"] == "zoom":
                    r["drawn"].append("map")
                    r["placed"] = f"map, labelled in the {ins['title']} inset"
                break
        else:
            if inside(r["lat"], r["lon"], g["region"]):
                r["drawn"], r["labelled"], r["placed"] = ["map"], "map", "map"
    rows.sort(key=lambda r: (g["networks"].index(r["network"]), r["station"]))
    with open(f"{out}.csv", "w", newline="") as f:
        if events is None:
            w = csv.DictWriter(f, fieldnames=["role", "network", "station", "lat", "lon", "elevation_m",
                                              "placed", "sets_chosen", "sets_probed", "coverage_min",
                                              "reason", "provenance"], extrasaction="ignore")
        else:
            w = csv.DictWriter(f, fieldnames=["role", "network", "station", "lat", "lon", "placed",
                                              "sequences", "quiet_windows", "magnitude", "depth_km",
                                              "time", "provenance"], extrasaction="ignore")
        w.writeheader()
        w.writerows(rows + (events or []))

    def stations(placed):
        if events is None:
            draw_stations(fig, pygmt, rows, g, placed, g["symbol"], g["labels"])
        else:
            draw_tests(fig, rows, events, g, group, placed, g["symbol"])

    grid = relief(ROOT / "data/maps" / f"etopo2022_30s_{group}.nc", g["region"])
    with pygmt.config(DIR_GSHHG=a.gshhg, DIR_DCW=a.dcw, GMT_DATA_UPDATE_INTERVAL="off",
                      FONT_ANNOT_PRIMARY="7p,Helvetica", FONT_LABEL="7p,Helvetica",
                      MAP_FRAME_TYPE="plain", MAP_FRAME_PEN="0.5p", MAP_TICK_LENGTH_PRIMARY="2p"):
        fig = pygmt.Figure()
        shade = pygmt.grdgradient(grid=str(grid), azimuth=315, normalize="e0.5")
        pygmt.makecpt(cmap="gray", series=[-3500, 3500], continuous=True)
        fig.grdimage(grid=str(grid), shading=shade, region=g["region"], projection=g["projection"],
                     cmap=True, transparency=35, frame=["WSne"] + g["frame"])
        fig.coast(shorelines="0.3p,gray20", borders="1/0.4p,gray25,-", resolution="h")
        stations("map")
        fig.legend(spec=str(legend_spec(Path(f"{out}.legend"), rows, g, events)),
                   position="jTR+o0.15c", box="+gwhite@10+p0.3p")
        # basemap's map_scale places the bar inside the frame; Figure.scalebar,
        # which replaces it in newer PyGMT, put its annotations outside.
        fig.basemap(map_scale="jBL+o0.4c/0.5c+w150k+f+l", frame=False)
        for ins in insets:
            if not any(ins["name"] in r["drawn"] for r in rows + (events or [])):
                continue
            if ins["kind"] == "zoom":
                lo, hi, la, lb = ins["region"]
                fig.plot(x=[lo, hi, hi, lo, lo], y=[la, la, lb, lb, la], pen="0.5p,black")
            ipath = ROOT / "data/maps" / f"etopo2022_30s_{ins['name']}.nc"
            if ins["kind"] == "zoom" and not ipath.exists():
                igrid = pygmt.grdcut(grid=str(grid), region=ins["region"])  # from the map's own cut
            else:
                igrid = relief(ipath, ins["region"])
            if ins["kind"] == "zoom":
                # At a zoom's scale the 30" cells would show as blocks: the
                # cut is resampled to 3" with bicubic interpolation first.
                igrid = pygmt.grdsample(grid=igrid if not isinstance(igrid, Path) else str(igrid),
                                        spacing="3s" if ipath.exists() else "10s", interpolation="c",
                                        region=ins["region"])
            with fig.inset(position=ins["position"], box="+gwhite+p0.5p,black"):
                ishade = pygmt.grdgradient(grid=igrid, azimuth=315, normalize="e0.5")
                fig.grdimage(grid=igrid, shading=ishade, region=ins["region"],
                             projection=ins["projection"], cmap=True, transparency=35)
                fig.coast(shorelines="0.3p,gray20", resolution="f")
                stations(ins["name"])
                title = ins["title"] + (" (boxed on the map)" if ins["kind"] == "zoom" else "")
                fig.text(position="TL", text=title, font="5p,Helvetica,black",
                         offset="0.08c/-0.08c", justify="TL", fill="white@25")
        fig.savefig(f"{out}.pdf")
        fig.savefig(f"{out}.png", dpi=300)
    Path(f"{out}.legend").unlink()
    print(f"-> {out}.pdf, {out}.png, {out}.csv")
    if events is not None:
        print(f"   {sum(1 for r in rows if r['role'] == USED)} stations in the selected sequences, "
              f"{sum(1 for r in rows if r['role'] == OTHER)} other usable "
              f"({sum(1 for r in rows if r['role'] == OTHER and r['drawn'])} drawn), {len(events)} sequences, "
              f"{sum(1 for r in rows if r['quiet_windows'])} stations in a quiet window: "
              + ", ".join(f"{t} {sum(1 for r in rows if r['quiet_windows'] and ring_class(r['quiet_windows']) == i)}"
                          for i, (_, _, t) in enumerate(RING_CLASSES)))
    for net in [] if events is not None else g["networks"]:
        n_use = sum(1 for r in rows if r["network"] == net and r["role"] == "usable")
        n_all = sum(1 for r in rows if r["network"] == net)
        print(f"   {net}: {n_use} usable of {n_all} probed")
    for r in rows + (events or []):
        if r["placed"] != "map":
            print(f"   {r['network']}.{r['station']} ({r['lat']:.2f}, {r['lon']:.2f}): {r['placed']}")
    for r in rows + (events or []):
        del r["drawn"], r["labelled"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--heldout", required=True, type=Path, help="the ayzek_heldout checkout")
    ap.add_argument("--out", required=True, help="path prefix: writes PREFIX_<group>.{pdf,png,csv}")
    ap.add_argument("--group", nargs="*", default=list(GROUPS), choices=list(GROUPS))
    ap.add_argument("--selection", action="append", default=[], metavar="GROUP=NAME",
                    help="draw HELDOUT/selections/NAME/selection.json's stations for GROUP "
                         "when it exists (an empty NAME turns the default off)")
    ap.add_argument("--tests", action="store_true",
                    help="the tests as run: stations of the selected sequences, quiet-window "
                         "stations and first events, from the final selections")
    ap.add_argument("--tests-selection", action="append", default=[], metavar="GROUP=NAME",
                    help="with --tests: the final selection of GROUP (default GROUP_now)")
    ap.add_argument("--gshhg", default=os.path.expanduser("~/.gmt/gshhg/gshhg-gmt-2.3.7"))
    ap.add_argument("--dcw", default=os.path.expanduser("~/.gmt/dcw/dcw-gmt-2.2.0"))
    a = ap.parse_args()

    import pygmt

    selection = {k: v["selection"] for k, v in GROUPS.items()}
    for s in a.selection:
        k, _, v = s.partition("=")
        selection[k] = v or None
    final = {k: f"{k}_now" for k in GROUPS}
    for s in a.tests_selection:
        k, _, v = s.partition("=")
        final[k] = v
    for group in a.group if a.tests else []:
        name = final[group]
        print(f"{group}: the tests of selections/{name}/selection.json")
        rows, events = from_tests(a.heldout / "selections" / name / "selection.json", name,
                                  GROUPS[group]["networks"])
        plot(group, rows, f"{a.out}_{group}", a, pygmt, events)
    for group in [] if a.tests else a.group:
        name = selection[group]
        path = a.heldout / "selections" / (name or "") / "selection.json"
        rows = from_tune_sets(a.heldout, group)
        if name and path.exists():
            print(f"{group}: stations from selections/{name}/selection.json")
            rows = from_selection(path, name, rows)
        plot(group, rows, f"{a.out}_{group}", a, pygmt)


if __name__ == "__main__":
    main()
