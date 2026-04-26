#!/usr/bin/env python3
"""Generate DSSP files for all PDB structures in a dataset directory.

Runs mkdssp on each original PDB file after cleaning records that
cause mkdssp 4.x to fail:

  - ANISOU records (cause incomplete residue assignments)
  - SIGATM/SIGUIJ records (cause zero-residue output)
  - Blank chain IDs in column 22 (cause parse failures)

Usage:
    python pipeline/run_dssp.py SRC_DIR [-j JOBS]

Output per entry:
    {pdbid}.dssp — classic DSSP output from the original PDB
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Record prefixes that cause mkdssp 4.x failures
_STRIP_PREFIXES = ("ANISOU", "SIGATM", "SIGUIJ")


def _find_mkdssp() -> str | None:
    """Return the path to the mkdssp executable, or None."""
    return shutil.which("mkdssp")


def _clean_and_run_dssp(
    pdb_path: Path,
    dssp_path: Path,
    mkdssp: str,
) -> tuple[bool, str]:
    """Clean a PDB file and run mkdssp on it.

    Builds a minimal PDB containing only HEADER, CRYST1, and
    coordinate records (ATOM/HETATM/TER) with ANISOU/SIGATM/SIGUIJ
    stripped and blank chain IDs filled.  This avoids all known
    mkdssp 4.x parsing issues with old PDB files.

    Returns (success, message).
    """
    # Find default chain ID from ATOM records
    default_chain = "A"
    with open(pdb_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("ATOM  ") and len(line) > 21:
                ch = line[21]
                if ch != " ":
                    default_chain = ch
                    break

    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".pdb")
    try:
        with os.fdopen(tmp_fd, "w") as tmp:
            with open(pdb_path, encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith(("HEADER", "CRYST1")):
                        tmp.write(line)
                    elif line.startswith("SEQRES"):
                        # Fill blank chain ID at column 12 (0-indexed 11)
                        out = line
                        if len(line) > 11 and line[11] == " ":
                            out = line[:11] + default_chain + line[12:]
                        tmp.write(out)
                    elif line.startswith(("ATOM  ", "HETATM", "TER   ")):
                        if line.startswith(_STRIP_PREFIXES):
                            continue
                        out = line
                        if len(line) > 21 and line[21] == " ":
                            out = line[:21] + default_chain + line[22:]
                        tmp.write(out)
            tmp.write("END\n")

        result = subprocess.run(
            [mkdssp, tmp_path, str(dssp_path)],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
    finally:
        os.unlink(tmp_path)

    if result.returncode != 0:
        stderr = result.stderr.strip().split("\n")[-1] if result.stderr else ""
        return False, f"mkdssp failed: {stderr}"

    return True, "ok"


def process_entry(
    entry_dir: Path,
    mkdssp: str,
    force: bool,
) -> tuple[str, bool, str]:
    """Generate DSSP file for one PDB entry.

    Returns (pdb_id, success, message).
    """
    pdb_id = entry_dir.name
    pdb_path = entry_dir / f"{pdb_id}.pdb"
    dssp_path = entry_dir / f"{pdb_id}.dssp"

    if not pdb_path.exists():
        return (pdb_id, False, "missing PDB file")

    if dssp_path.exists() and not force:
        return (pdb_id, True, "exists (skipped)")

    try:
        success, msg = _clean_and_run_dssp(pdb_path, dssp_path, mkdssp)
        return (pdb_id, success, msg)
    except Exception as exc:
        return (pdb_id, False, str(exc))


def collect_entries(src_dir: Path) -> list[Path]:
    """Find all entry directories (two-level: prefix/pdbid/)."""
    entries: list[Path] = []
    for prefix_dir in sorted(src_dir.iterdir()):
        if not prefix_dir.is_dir():
            continue
        for entry_dir in sorted(prefix_dir.iterdir()):
            if not entry_dir.is_dir():
                continue
            entries.append(entry_dir)
    return entries


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate DSSP files for all PDB structures in a dataset. "
            "Cleans ANISOU/SIGATM/SIGUIJ records and blank chain IDs "
            "before running mkdssp."
        )
    )
    parser.add_argument(
        "src_dir", type=Path,
        help="Dataset directory (two-level: prefix/pdbid/)",
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=8,
        help="Number of parallel workers (default: 8)",
    )
    parser.add_argument(
        "-f", "--force", action="store_true",
        help="Regenerate existing .dssp files",
    )
    args = parser.parse_args()

    mkdssp = _find_mkdssp()
    if mkdssp is None:
        print(
            "Error: mkdssp not found. Install with: apt install dssp",
            file=sys.stderr,
        )
        sys.exit(1)

    if not args.src_dir.is_dir():
        print(f"Error: {args.src_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    print(f"Scanning {args.src_dir} ...")
    entries = collect_entries(args.src_dir)
    print(f"Found {len(entries)} entries")
    print(f"mkdssp: {mkdssp}")

    done = 0
    errors = 0
    skipped = 0

    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(process_entry, entry, mkdssp, args.force): entry
            for entry in entries
        }
        for future in as_completed(futures):
            pdb_id, success, msg = future.result()
            done += 1
            if not success:
                errors += 1
                print(f"  ERROR {pdb_id}: {msg}", file=sys.stderr)
            elif msg == "exists (skipped)":
                skipped += 1
            if done % 500 == 0 or done == len(entries):
                print(
                    f"  {done}/{len(entries)} entries processed",
                    flush=True,
                )

    generated = done - errors - skipped
    print(f"\nDone: {generated} generated, {skipped} skipped, {errors} errors")


if __name__ == "__main__":
    main()
