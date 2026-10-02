"""The onset code that ayzek's exports and fixtures are made with.

Copied from the onset repository (github.com/hogib/onset, `src/onset/`):

- for `tools/export_transformer.py`: `conditioning.py` and `dsp.py` from
  branch `onset-sp`, commit d06de22 (the distance-only geometry head), and
  `config.py` and `model.py` from commit 5a31ce8 (which add the optional
  binned dt head; the export leaves that head out);
- for `tools/export_pd.py`: `pd.py` and `pd_fit.py`, from commit badbfd7 (the
  Pd magnitude relation, with the rules for channels that do not record and
  for values inside a coda).

They are unchanged except that `from onset.config import` became
`from .config import` and `from onset.pd import` became `from .pd import`.
An export must use the code its model was trained or fitted with, so copy
these again when onset changes them; `test_transformer` and `test_pd` then
fail if they no longer agree with the C++ implementation.
"""
