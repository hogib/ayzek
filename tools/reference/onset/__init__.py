"""The onset transformer's model code, for `tools/export_transformer.py`.

Copied from the onset repository (github.com/hogib/onset, `src/onset/`,
branch `onset-sp` on commit 21db821, the distance-only geometry head):
`config.py`, `model.py`, `conditioning.py` and `dsp.py`,
unchanged except that `from onset.config import` became
`from .config import`. The export must load a checkpoint with the code it
was trained with, so copy these again when onset's model, conditioning or
filter changes; `test_transformer` then fails if they no longer agree with
the C++ implementation.
"""
