"""Map of the development data: the TU stations of the AFAD archive, the KO
stations of the public KOERI sets, the epicentres of the replayed sequences,
and the area within 250 km of the quiet-window stations, where no catalogued
event may lie during a quiet window.

    uv run --project tools --with pygmt --with xarray --with netcdf4 --with shapely \\
        python tools/plot_map_development.py --out ../ayzek_paper/figures/fig_map_development

Relief: ETOPO 2022 (30 arc-seconds, NOAA NCEI), cut over OPeNDAP once and kept
in data/maps/. Coastlines and borders: GSHHG and DCW, from --gshhg and --dcw
(default ~/.gmt/gshhg/gshhg-gmt-2.3.7 and ~/.gmt/dcw/dcw-gmt-2.2.0), since GMT's
own data server may be unreachable. Writes OUT.pdf, OUT.png and OUT.csv (every
point drawn, with its role).
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGION = [24.3, 33.2, 36.2, 44.2]
ETOPO = ("https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/30s/"
         "30s_surface_elev_netcdf/ETOPO_2022_v1_30s_N90W180_surface.nc")

TU = ["ARNA", "BAND", "CATL", "CMH", "DEMI", "ELBA", "GCAM", "KAND", "KURT", "KIRK",
      "MANT", "SEMS", "VIZE"]
# The largest catalogued event of each replayed sequence (tests/catalogs),
# with the label's position: the two Sındırgı events lie 7 km apart and share
# one label.
EVENTS = [
    ("Marmara Mw 6.2", "2025-04-23", "marmara", "0.35c/0c", "LM"),
    ("Sındırgı Mw 6.1", "2025-08-10", "sindirgi61", None, None),
    ("Sındırgı Mw 6.1 and 4.9", "2025-11-10", "demo", "-0.35c/0.1c", "RM"),
    ("Simav Mw 4.8", "2026-04-11", "simav26", "0.3c/-0.35c", "LM"),
]
# Station labels moved off their neighbours (dx/dy); the rest sit above (TU)
# or below (KO) their symbol.
LABEL_OFFSET = {
    "CATL": "-0.3c/0.22c", "ARNA": "0.4c/0.12c", "ELBA": "-0.12c/-0.25c",
    "SLVT": "-0.35c/0.12c", "CTKS": "0c/0.3c", "BGKT": "0.38c/-0.15c",
    "BOTS": "0c/-0.25c", "TKR": "0c/-0.25c", "MRMT": "-0.38c/0c", "EDC": "-0.3c/-0.2c",
    "KCTX": "0.38c/-0.05c", "CMH": "-0.33c/0c", "BAND": "0.36c/0.05c", "ARMT": "-0.42c/0c",
    "KURT": "0.36c/0.08c", "DEMI": "-0.38c/0c", "SIMA": "0.05c/-0.25c", "TVSB": "0.38c/0c",
    "GEDZ": "0.38c/0c", "MANT": "-0.4c/0c", "KULA": "0.38c/-0.05c",
}


def stations_table():
    out = {}
    for r in csv.DictReader(open(ROOT / "models/stations.csv")):
        out.setdefault(r["code"], (float(r["lat"]), float(r["lon"])))
    return out


def largest(cat):
    rows = list(csv.DictReader(open(ROOT / "tests/catalogs" / f"{cat}.csv", encoding="utf-8-sig")))
    m = max(rows, key=lambda r: float(r["Magnitude"]))
    return float(m["Latitude"]), float(m["Longitude"])


def circle(lat, lon, km, n=180):
    """Points of the small circle of radius `km` about (lat, lon), on a sphere."""
    R = 6371.0
    d = km / R
    la, lo = math.radians(lat), math.radians(lon)
    pts = []
    for i in range(n + 1):
        b = 2 * math.pi * i / n
        la2 = math.asin(math.sin(la) * math.cos(d) + math.cos(la) * math.sin(d) * math.cos(b))
        lo2 = lo + math.atan2(math.sin(b) * math.sin(d) * math.cos(la),
                              math.cos(d) - math.sin(la) * math.sin(la2))
        pts.append((math.degrees(lo2), math.degrees(la2)))
    return pts


def relief(path):
    if path.exists():
        return path
    import xarray as xr
    path.parent.mkdir(parents=True, exist_ok=True)
    z = xr.open_dataset(ETOPO)["z"].sel(lat=slice(REGION[2], REGION[3]),
                                        lon=slice(REGION[0], REGION[1])).load()
    z = z.astype("float32")
    z.attrs = {"units": "m", "long_name": "elevation (ETOPO 2022, 30 arc-seconds)"}
    z.to_dataset(name="z").to_netcdf(path, format="NETCDF3_64BIT")
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True, help="path without extension")
    ap.add_argument("--gshhg", default=os.path.expanduser("~/.gmt/gshhg/gshhg-gmt-2.3.7"))
    ap.add_argument("--dcw", default=os.path.expanduser("~/.gmt/dcw/dcw-gmt-2.2.0"))
    a = ap.parse_args()

    import pygmt
    from shapely.geometry import Polygon
    from shapely.ops import unary_union

    st = stations_table()
    quiet = []
    for line in (ROOT / "data/runs/quiet100").glob("*.rec"):
        quiet += [l.split("\t")[1] for l in line.read_text().splitlines() if l.startswith("S\t")]
    quiet = sorted(set(quiet))
    ko = {}
    for d in ("demo_ko", "marmara_ko", "quiet_ko"):
        for s in json.loads((ROOT / "data" / d / "dataset.json").read_text())["stations"]:
            ko.setdefault(s["station"], (s["lat"], s["lon"]))
    zone = unary_union([Polygon(circle(*st[c], 250.0)) for c in quiet])

    rows = [("TU", c, *st[c]) for c in TU] + [("KO", c, *ll) for c, ll in sorted(ko.items())]
    rows += [("epicentre", name, *largest(cat)) for name, _, cat, _, _ in EVENTS]
    with open(f"{a.out}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["role", "name", "lat", "lon"])
        w.writerows(rows)

    grid = relief(ROOT / "data/maps/etopo2022_30s_development.nc")
    with pygmt.config(DIR_GSHHG=a.gshhg, DIR_DCW=a.dcw, GMT_DATA_UPDATE_INTERVAL="off",
                      FONT_ANNOT_PRIMARY="7p,Helvetica", FONT_LABEL="7p,Helvetica",
                      MAP_FRAME_TYPE="plain", MAP_FRAME_PEN="0.5p", MAP_TICK_LENGTH_PRIMARY="2p"):
        fig = pygmt.Figure()
        shade = pygmt.grdgradient(grid=str(grid), azimuth=315, normalize="e0.5")
        pygmt.makecpt(cmap="gray", series=[-3500, 3500], continuous=True)
        fig.grdimage(grid=str(grid), shading=shade, region=REGION, projection="M12c", cmap=True,
                     transparency=35, frame=["WSne", "xa2f1", "ya1f0.5"])
        fig.coast(shorelines="0.3p,gray20", borders="1/0.4p,gray25,-", resolution="h")
        x, y = zone.exterior.xy
        fig.plot(x=list(x), y=list(y), pen="0.8p,40/100/170,--",
                 label="within 250 km of a quiet-window station")
        tu = [r for r in rows if r[0] == "TU"]
        fig.plot(x=[r[3] for r in tu], y=[r[2] for r in tu], style="t0.26c", fill="200/40/40",
                 pen="0.3p,black", label="TU (AFAD archive)")
        kk = [r for r in rows if r[0] == "KO"]
        fig.plot(x=[r[3] for r in kk], y=[r[2] for r in kk], style="i0.24c", fill="40/110/190",
                 pen="0.3p,black", label="KO (public KOERI sets)")
        ev = [r for r in rows if r[0] == "epicentre"]
        fig.plot(x=[r[3] for r in ev], y=[r[2] for r in ev], style="a0.38c", fill="250/200/20",
                 pen="0.4p,black", label="largest event of a replay")
        for r in tu:
            fig.text(x=r[3], y=r[2], text=r[1], font="5p,Helvetica,black",
                     offset=LABEL_OFFSET.get(r[1], "0c/0.24c"))
        for r in kk:
            fig.text(x=r[3], y=r[2], text=r[1], font="5p,Helvetica,40/80/150",
                     offset=LABEL_OFFSET.get(r[1], "0c/-0.24c"))
        for name, _, _, off, just in EVENTS:
            if off is None:
                continue
            r = next(x for x in ev if x[1] == name)
            fig.text(x=r[3], y=r[2], text=name, font="6p,Helvetica-Bold,black", justify=just,
                     offset=off, fill="white@25")
        fig.legend(position="jTR+o0.15c", box="+gwhite@10+p0.3p")
        # basemap's map_scale places the bar inside the frame; Figure.scalebar,
        # which replaces it in newer PyGMT, put its annotations outside.
        fig.basemap(map_scale="jBL+o0.4c/0.5c+w100k+f+l", frame=False)
        fig.savefig(f"{a.out}.pdf")
        fig.savefig(f"{a.out}.png", dpi=300)
    print(f"-> {a.out}.pdf, {a.out}.png, {a.out}.csv")


if __name__ == "__main__":
    main()
