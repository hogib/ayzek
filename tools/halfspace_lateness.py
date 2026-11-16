#!/usr/bin/env python3
"""How late the matcher's half-space P prediction is against IASP91 (EXPLORATORY).

    uv run --offline --project tools python tools/halfspace_lateness.py

ayzek's matcher (Network::match, src/pipeline/network.cpp) predicts P as
hypot(D, depth) / 6 km/s. This prints, per source depth and epicentral distance,
that prediction minus the first P arrival of obspy's TauP in iasp91 (phases p, P,
Pn, Pg, Pdiff; distance in degrees = km / 111.19, a 6371 km sphere), in seconds.
Nothing is read or written; the table goes to stdout.
"""
import math

from obspy.taup import TauPyModel

VP = 6.0                                    # network.hpp
P_PHASES = ["p", "P", "Pn", "Pg", "Pdiff"]  # as tools/heldout_iasp91.py
DEPTHS_KM = (40, 100, 150, 175, 200, 211)
DISTANCES_KM = (0, 50, 100, 150, 200, 250)


def main():
    model = TauPyModel("iasp91")
    print("depth_km," + ",".join(f"{d}km" for d in DISTANCES_KM))
    for z in DEPTHS_KM:
        row = []
        for d in DISTANCES_KM:
            arrivals = model.get_travel_times(source_depth_in_km=z, distance_in_degree=d / 111.19,
                                              phase_list=P_PHASES)
            row.append(math.hypot(d, z) / VP - min(a.time for a in arrivals))
        print(f"{z}," + ",".join(f"{x:.1f}" for x in row))


if __name__ == "__main__":
    main()
