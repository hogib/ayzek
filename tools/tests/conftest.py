"""Shared helpers for the tests of the analysis tools in tools/.

    uv run --project tools --with pytest pytest tools/tests -q

The tools are standard-library scripts, not a package: tests import them from
tools/ (put on sys.path here), and end-to-end tests run a tool as a subprocess
with a stand-in for the binary it drives (write_fake_binary), which prints what
the test prescribes. No expected value comes from running the code under test:
each test says where its numbers come from (arithmetic by hand, scipy, the
report's own summary, or an invariant).

Fixtures (tools/tests/fixtures/): stdout of real ayzek runs on data/demo_ko
(the six KO stations, open KOERI data, 2025-11-10 18:04-18:35 UTC), captured on
2026-10-29 with the main checkout's binaries and data, read-only, and the
catalogue excerpt tests/catalogs/demo_ko.csv (identical to
tests/catalogs/demo.csv):

  report_new_learned_demo_ko.txt      the main checkout's build-dev/app/ayzek
                                      (built 2026-10-29 11:45; prints the pre-S
                                      line documented in f54e958),
                                      --scores-in data/runs/cache_demo_ko
  report_new_stalta3_demo_ko.txt      the same build, --detector stalta
                                      --stalta-band 1,10 --sta 0.5 --lta 10
                                      --stalta-on 3 --stalta-off 1.5
  report_frozen_learned_demo_ko.txt   data/runs/tune_dev/ayzek (commit
                                      9c78e4e, before the pre-S line),
                                      --scores-in data/runs/cache_demo_ko
  report_new_learned_demo_ko_p127.txt build-dev, learned, against the
  report_new_learned_demo_ko_m263.txt catalogue shifted by +127 s and -263 s
                                      (shifted by an independent script with
                                      calendar.timegm, not by the tools)

Every run: `ayzek --speed 0 --no-color --catalog CAT [options]
data/demo_ko/{DKL,KCTX,KULA,SIMA,TVSB,YAYO}.mseed` from a checkout of this
branch (models/ identical to the main checkout's), magnitude on; stdout only,
stderr was empty.
"""
import os
import stat
import struct
import subprocess
import sys
import textwrap
import time
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[1]
REPO = TOOLS.parent
FIXTURES = Path(__file__).resolve().parent / "fixtures"

if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))


def fixture_text(name):
    return (FIXTURES / name).read_text()


def run_tool(name, *args, env=None, cwd=None):
    """Runs tools/NAME.py with this interpreter; returns the CompletedProcess."""
    e = dict(os.environ)
    e.update(env or {})
    return subprocess.run([sys.executable, str(TOOLS / f"{name}.py"), *map(str, args)],
                          capture_output=True, text=True, env=e, cwd=cwd)


def write_fake_binary(path, body):
    """An executable Python script standing in for ayzek or network_subsets."""
    path.write_text(f"#!{sys.executable}\n" + textwrap.dedent(body))
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def write_mseed_stub(path, first, last):
    """Two 512-byte records whose fixed headers carry the start times `first`
    and `last` (epoch seconds). Only the BTime field (bytes 20-30: year, day of
    year, hour, minute, second, unused, 0.0001 s) is filled; that is all
    tools/sweep_trigger.data_span reads. Fields from time.gmtime, independently
    of the tool's datetime arithmetic."""
    recs = b""
    for t in (first, last):
        g = time.gmtime(int(t))
        frac = round((t - int(t)) * 10000)
        hdr = bytearray(512)
        hdr[20:30] = struct.pack(">HHBBBBH", g.tm_year, g.tm_yday, g.tm_hour, g.tm_min,
                                 g.tm_sec, 0, frac)
        recs += bytes(hdr)
    path.write_bytes(recs)
    return path
