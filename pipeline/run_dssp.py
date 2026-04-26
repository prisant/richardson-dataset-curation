#!/usr/bin/env python3
"""Generate DSSP files for all PDB structures in a dataset directory.

Runs mkdssp on each original PDB file after building a minimal
protein-only PDB that avoids the records mkdssp 4.x mishandles:

  - ANISOU records (cause incomplete residue assignments)
  - SIGATM/SIGUIJ records (cause zero-residue output)
  - HETATM records (collide with SEQRES protein alignment, silently
    nulling the structure's DSSP output for ligand-heavy entries
    like 1a1y, 1gdo, 1tax, 2myr in top500)
  - Blank chain IDs in column 12 (SEQRES) and column 22 (ATOM)

The cleaned input contains HEADER + CRYST1 + SEQRES + ATOM + TER only.

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


def _find_mkdssp() -> str | None:
    """Return the path to the mkdssp executable, or None."""
    return shutil.which("mkdssp")


def _clean_and_run_dssp(
    pdb_path: Path,
    dssp_path: Path,
    mkdssp: str,
) -> tuple[bool, str]:
    """Clean a PDB file and run mkdssp on it.

    Builds a minimal PDB containing HEADER + CRYST1 + SEQRES + ATOM +
    TER + protein-only HETATM records, with blank chain IDs filled.

    HETATM records are filtered: only those whose residue name appears
    in SEQRES or ATOM (i.e., modified amino-acid residues like MSE,
    PCA, MEN, PTR, LOV that some old PDBs deposit as HETATM) are
    written.  Non-protein HETATM (HOH, SO4, GOL, NAG, ZN, ...) are
    omitted because mkdssp 4.x silently produces zero-residue output
    when those heteroatoms collide with the SEQRES protein alignment
    (observed on top500: 1a1y, 1gdo, 1tax, 2myr).

    Returns (success, message).
    """
    # First pass: collect default chain + protein residue names from
    # SEQRES and ATOM records.
    default_chain = "A"
    chain_seen = False
    protein_residues: set[str] = set()
    with open(pdb_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("ATOM  ") and len(line) > 21:
                if not chain_seen:
                    ch = line[21]
                    if ch != " ":
                        default_chain = ch
                        chain_seen = True
                if len(line) >= 20:
                    protein_residues.add(line[17:20].strip())
            elif line.startswith("SEQRES") and len(line) > 19:
                idx = 19
                while idx + 3 <= len(line):
                    name = line[idx:idx + 3].strip()
                    if name:
                        protein_residues.add(name)
                    idx += 4

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
                    elif line.startswith(("ATOM  ", "TER   ")):
                        out = line
                        if len(line) > 21 and line[21] == " ":
                            out = line[:21] + default_chain + line[22:]
                        tmp.write(out)
                    elif line.startswith("HETATM"):
                        resname = line[17:20].strip() if len(line) >= 20 else ""
                        if resname not in protein_residues:
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
