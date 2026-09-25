"""Scratch driver: run pipeline/run_dssp.py cleaning with an alternate mkdssp.

Calls run_dssp._clean_and_run_dssp (production cleaning + chain-batch
fallback) but with a chosen mkdssp binary and a separate flat output
directory, so the dataset's own {pdbid}.dssp files are untouched.
Checks each output is classic DSSP format (line 1 starts with '====').

Usage: python run_dssp_alt.py SRC_DIR OUT_DIR MKDSSP [-j N]
"""
import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, "/home/prisant/Desktop/Projects/richardson-dataset-curation/pipeline")
import run_dssp  # noqa: E402


def one(args):
    entry, out_dir, mkdssp = args
    pdb = entry / f"{entry.name}.pdb"
    out = out_dir / f"{entry.name}.dssp"
    if not pdb.exists():
        return entry.name, False, "missing PDB file"
    try:
        ok, msg = run_dssp._clean_and_run_dssp(pdb, out, mkdssp)
    except Exception as exc:  # noqa: BLE001
        return entry.name, False, str(exc)
    if ok and not out.read_text(errors="replace").startswith("===="):
        return entry.name, False, "output is not classic DSSP format"
    return entry.name, ok, msg


def main():
    p = argparse.ArgumentParser()
    p.add_argument("src_dir", type=Path)
    p.add_argument("out_dir", type=Path)
    p.add_argument("mkdssp")
    p.add_argument("-j", "--jobs", type=int, default=8)
    a = p.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    entries = run_dssp.collect_entries(a.src_dir)
    fails = []
    with ProcessPoolExecutor(a.jobs) as ex:
        for n, (pid, ok, msg) in enumerate(ex.map(one, [(e, a.out_dir, a.mkdssp) for e in entries]), 1):
            if not ok:
                fails.append((pid, msg))
                print(f"  ERROR {pid}: {msg[:200]}", flush=True)
            elif msg != "ok":
                print(f"  NOTE {pid}: {msg[:200]}", flush=True)
            if n % 500 == 0:
                print(f"  {n}/{len(entries)}", flush=True)
    print(f"Done: {len(entries) - len(fails)} generated, {len(fails)} errors")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
