"""The AYZW tensor file writer, shared by the export scripts. See docs/impl/02-weights.md."""
import json
import struct

import numpy as np

DTYPES = {np.dtype(np.float32): 0, np.dtype(np.float64): 1, np.dtype(np.int64): 2, np.dtype(np.uint8): 3}


def write_ayzw(path, tensors, meta):
    """magic, version, count, then (name, dtype, dims, bytes) per tensor."""
    path.parent.mkdir(parents=True, exist_ok=True)
    items = dict(tensors)
    items["__meta__"] = np.frombuffer(json.dumps(meta, indent=1).encode(), dtype=np.uint8)
    with open(path, "wb") as f:
        f.write(b"AYZW" + struct.pack("<II", 1, len(items)))
        for name, arr in items.items():
            arr = np.ascontiguousarray(arr)
            if arr.dtype not in DTYPES:
                arr = arr.astype(np.float32)
            nb = name.encode()
            f.write(struct.pack("<H", len(nb)) + nb)
            f.write(struct.pack("<BB", DTYPES[arr.dtype], arr.ndim))
            f.write(struct.pack(f"<{arr.ndim}I", *arr.shape))
            f.write(struct.pack("<Q", arr.nbytes))
            f.write(arr.astype(arr.dtype.newbyteorder("<")).tobytes())
    print(f"  {path}  ({len(items)} tensors, {path.stat().st_size / 1e3:.0f} kB)")


