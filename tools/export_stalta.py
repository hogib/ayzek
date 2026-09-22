"""Exports reference outputs for the STA/LTA trigger (src/stalta.hpp).

From 4 minutes of DEMI around the 2025-11-10 18:20:51 M4.9:
  raw         (n, 3) Z/N/E counts
  filtered    (n, 3) raw minus its first sample, then scipy sosfilt with
              butter(4, f_lo, 'highpass') and butter(4, f_hi, 'lowpass')
  ratio_z     obspy.signal.trigger.recursive_sta_lta of the filtered Z component
  ratio_zne   the same recursion on the sum of the squares of all three components
              (ObsPy has no three-component version; written out here)
  sos_hp, sos_lp  scipy's second-order sections

Run from the ayzek root:

    uv run --project tools python tools/export_stalta.py

Writes data/fixtures/stalta.ayzw.
"""
import sys
from pathlib import Path

import numpy as np
from obspy import UTCDateTime, read
from obspy.signal.trigger import recursive_sta_lta
from scipy import signal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ayzw import write_ayzw  # noqa: E402

FS = 100.0
F_LO, F_HI = 1.0, 10.0
NSTA, NLTA = 100, 3000
START, N = "2025-11-10T18:18:00", 24000


def main():
    st = read("data/demo/DEMI.mseed")
    st.merge(method=1, fill_value=None)
    st = st.split()
    t = UTCDateTime(START)
    raw = np.empty((N, 3))
    for k, comp in enumerate("ZNE"):
        tr = next(tr for tr in st if tr.stats.channel.endswith(comp) and tr.stats.starttime <= t < tr.stats.endtime)
        i = int(round((t - tr.stats.starttime) * FS))
        if i + N > tr.stats.npts:
            raise ValueError(f"{comp}: gap within the fixture span")
        raw[:, k] = tr.data[i:i + N]

    sos_hp = signal.butter(4, F_LO, "highpass", fs=FS, output="sos")
    sos_lp = signal.butter(4, F_HI, "lowpass", fs=FS, output="sos")
    filtered = np.stack([signal.sosfilt(sos_lp, signal.sosfilt(sos_hp, raw[:, c] - raw[0, c])) for c in range(3)],
                        axis=1)

    ratio_z = recursive_sta_lta(filtered[:, 0], NSTA, NLTA)

    e = (filtered ** 2).sum(axis=1)
    ratio_zne = np.zeros(N)
    sta, lta = 0.0, np.finfo(0.0).tiny
    for i in range(1, N):
        sta = e[i] / NSTA + (1 - 1 / NSTA) * sta
        lta = e[i] / NLTA + (1 - 1 / NLTA) * lta
        ratio_zne[i] = sta / lta
    ratio_zne[:NLTA] = 0

    write_ayzw(Path("data/fixtures/stalta.ayzw"),
               {"raw": raw, "filtered": filtered, "ratio_z": ratio_z, "ratio_zne": ratio_zne,
                "sos_hp": sos_hp, "sos_lp": sos_lp},
               {"start": START, "fs": FS, "f_lo": F_LO, "f_hi": F_HI, "nsta": NSTA, "nlta": NLTA,
                "max_ratio_z": float(ratio_z.max()), "max_ratio_zne": float(ratio_zne.max())})


if __name__ == "__main__":
    main()
