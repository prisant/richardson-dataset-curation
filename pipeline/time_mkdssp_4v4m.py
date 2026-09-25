"""Scratch diagnostic: time run_dssp._clean_and_run_dssp on 4v4m with no timeout.

Wraps subprocess.run to drop the 120 s timeout and log each mkdssp call's
wall time, return code and last stderr line. Output .dssp goes to scratchpad.
"""
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, "/home/prisant/Desktop/Projects/richardson-dataset-curation/pipeline")
import run_dssp  # noqa: E402

_real_run = subprocess.run


def timed_run(cmd, *args, **kwargs):
    kwargs.pop("timeout", None)
    t0 = time.monotonic()
    r = _real_run(cmd, *args, **kwargs)
    last = (r.stderr or "").strip().split("\n")[-1] if kwargs.get("capture_output") else ""
    print(f"  mkdssp call: {time.monotonic() - t0:7.1f} s  rc={r.returncode}  stderr: {last[:120]}", flush=True)
    return r


subprocess.run = timed_run

pdb = Path.home() / "Desktop/Data/top2018_pdbs_full_filtered_hom70/4v/4v4m/4v4m.pdb"
out = Path.cwd() / "4v4m_timed.dssp"
t0 = time.monotonic()
ok, msg = run_dssp._clean_and_run_dssp(pdb, out, "/usr/bin/mkdssp")
print(f"total {time.monotonic() - t0:.1f} s  ok={ok}  {msg}")
