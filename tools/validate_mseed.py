"""Check mseed_dump against ObsPy, which wraps libmseed, the reference decoder.

Steim2 is lossless integer compression, so there is no tolerance: every sample
must be identical. Also checks that sample counts per channel agree and that
every trace ObsPy assembles begins at a record start we also saw.

    uv run --with obspy --with numpy python tools/validate_mseed.py FILE.mseed DUMP_DIR
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from obspy import read

path, dump = Path(sys.argv[1]), Path(sys.argv[2])
idx = pd.read_csv(dump / "index.csv")
st = read(str(path))

ok = True
for cha in sorted({t.stats.channel for t in st}):
    traces = [t for t in st if t.stats.channel == cha]
    ref = np.concatenate([t.data.astype(np.int32) for t in traces])
    ours = np.fromfile(dump / f"{cha}.i32", dtype=np.int32)

    same_len = len(ref) == len(ours)
    identical = same_len and np.array_equal(ref, ours)
    starts = set(idx[idx.channel == cha].start_ns)
    trace_starts_known = all(t.stats.starttime.ns in starts for t in traces)

    first_diff = ""
    if same_len and not identical:
        i = int(np.argmax(ref != ours))
        first_diff = f"  first mismatch at sample {i}: obspy {ref[i]} vs ours {ours[i]}"

    print(f"{cha}  obspy {len(ref):>11,}  ours {len(ours):>11,}  "
          f"{'IDENTICAL' if identical else 'MISMATCH'}  "
          f"({len(traces)} traces, starts {'match' if trace_starts_known else 'DO NOT match'})"
          f"{first_diff}")
    ok &= identical and trace_starts_known

print("PASS: bit-exact against libmseed" if ok else "FAIL")
sys.exit(0 if ok else 1)
