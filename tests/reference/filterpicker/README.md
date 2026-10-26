# FilterPicker5, Anthony Lomax's C reference

Used only by `tests/test_filterpicker.cpp`, to check the port in
`src/filterpicker.cpp` against the original. Nothing here is linked into
`ayzek`.

- Source: `picker_20130212.tgz` from <http://alomax.free.fr/FilterPicker/archive/>,
  sha256 `3e0943e3046e9263ff863a6bc601596206d810456af7dd46bae95750f1cc1222`,
  downloaded 2026-10-26. FilterPicker5.c v5.0.3 (2012-10-19).
- Files copied unchanged: `FilterPicker5.c`, `FilterPicker5.h`,
  `FilterPicker5_Memory.c`, `FilterPicker5_Memory.h`, `PickData.c`,
  `PickData.h`, `ew_bridge.h`, `COPYING`.
- License: GPL-2.0-or-later (`COPYING`), Copyright (C) 2006-2012 Anthony Lomax.
  ayzek is GPL-3.0-or-later.
- Added here: `fp_reference.c` and `fp_reference.h`, which run `Pick()` over a
  trace in packets (the first packet one sample, so that the C function starts
  its first difference from the first sample, as the port does) and return the
  picks and the characteristic function.

Reference: Lomax, A., C. Satriano and M. Vassallo (2012), Automatic picker
developments and optimization: FilterPicker, a robust, broadband picker for
real-time seismic monitoring and earthquake early warning, Seismol. Res. Lett.
83(3), 531-540, doi:10.1785/gssrl.83.3.531.
