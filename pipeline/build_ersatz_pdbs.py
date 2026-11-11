#!/usr/bin/env python3
"""Build ersatz full-structure PDB files for Richardson Lab datasets.

For each entry in the dataset, creates an "ersatz" PDB by starting with
the Reduce-processed FH file (which has correct NQH flips on all chains)
and stripping hydrogens and problematic headers.  Also generates per-chain
residue mask files from the pruned files.

The ersatz files preserve multi-chain context needed for correct DSSP
secondary structure assignments, with NQH corrections applied to all
chains.

Prerequisites:
    Run scripts/run_reduce.py first to generate {pdbid}FH.pdb files.

Usage:
    python scripts/build_ersatz_pdbs.py SRC_DIR [-j JOBS] [-s SUFFIX]
                                        [--masks-only] [--rcsb-remap]

Arguments:
    SRC_DIR     Dataset directory with two-level structure (prefix/pdbid/)
    -j, --jobs  Number of parallel workers (default: 8)
    -s, --suffix
                Pruned file suffix to match (default: pruned_all).
                Use "pruned_mc" for mainchain-filtered datasets.

                The script searches for files matching the pattern
                {pdbid}_{chain}_{suffix}.pdb in each entry directory.
    --masks-only
                Rebuild the .mask files only; leave the ersatz PDBs as
                they are (they have been modified since they were built).
    --rcsb-remap
                Apply the RCSB remaps in _RCSB_REMAPS (top2018_full only;
                see that table).

Examples:
    # Full-filtered dataset (default suffix "pruned_all")
    python scripts/build_ersatz_pdbs.py top2018_pdbs_full_filtered_hom70/

    # Mainchain-filtered dataset
    python scripts/build_ersatz_pdbs.py top2018_pdbs_mc_filtered_hom70/ -s pruned_mc

Input per entry (e.g., 1a/1a2z/):
    1a2z.pdb              - original full multi-chain PDB
    1a2zFH.pdb            - Reduce output with H + NQH flips
    1a2z_C_{suffix}.pdb   - pruned chain (for mask and fragment records)

Output per entry (in-place):
    1a2z_ersatz.pdb       - full PDB with NQH corrections, no H, clean headers
    1a2z_C.mask           - pruned residue list (chain resnum icode)

Note: Some structures contribute multiple quality-filtered chains
(e.g., 5t5i has chains A, C, G, J, L).  The script processes ALL
pruned files in each entry directory and generates a separate mask
file for each chain.

Masks hold ATOM residues and polymer HETATM residues (modified residues
such as MSE or SMC that the lab kept); see polymer_residues.py.  Before
2026-11-08 they held ATOM residues only.
"""

from __future__ import annotations

import argparse
import importlib.util
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "polymer_residues", Path(__file__).with_name("polymer_residues.py"))
_polymer = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_polymer)

# PDB record types to keep in ersatz files (same as clean_pdb_headers.py)
_KEEP_PREFIXES = (
    "HEADER", "CRYST1",
    "SCALE1", "SCALE2", "SCALE3",
    "ORIGX1", "ORIGX2", "ORIGX3",
    "MTRIX1", "MTRIX2", "MTRIX3",
    "MODEL ", "ENDMDL", "END   ",
    "ATOM  ", "HETATM", "ANISOU", "TER   ",
)


# Entries whose RCSB file (used for the ersatz PDB and hence for pydangle
# output) has different chain IDs or numbering from the lab's pruned file.
# Found and remapped by hand on 2026-03-15 (commit 9e288f2,
# top2018_full/top2018full_issues.md section 8); recorded here so that a
# mask rebuild reproduces them.  Keyed by (pdb_id, chain in the pruned file
# name); the function maps a mask key in pruned-file coordinates to RCSB
# coordinates.  The pruned files' USER INC records were edited in place at
# the same time and are not touched by this script.
_RCSB_REMAPS = {
    # RCSB renamed chain D to E; the pruned file was copied to
    # 4ct3_E_pruned_all.pdb (atoms still chain D).
    ("4ct3", "E"): lambda k: ("E", k[1], k[2]) if k[0] == "D" else k,
    # RCSB renumbered chain A: lab 87-353 is RCSB 402-668 (+315).
    ("5tvo", "A"): lambda k: (k[0], k[1] + 315, k[2]) if k[0] == "A" else k,
    # RCSB split chain B: residues 66 and above are now chain C.
    ("1xmz", "B"): lambda k: ("C", k[1], k[2]) if k[0] == "B" and k[1] >= 66 else k,
}
# The original 4ct3 chain-D pruned file is kept on disk but superseded by
# 4ct3_E_pruned_all.pdb; no mask is built from it.
_RCSB_SKIP = {("4ct3", "D")}


def _build_mask(pruned_path: Path) -> list[tuple[str, int, str]]:
    """Extract residue identifiers from pruned file for mask output.

    ATOM residues and polymer HETATM residues of the kept chain, in file
    order (polymer_residues.mask_lines).

    Returns list of (chain_id, resnum, icode) tuples, one per residue.
    """
    seen: set[tuple[str, int, str]] = set()
    residues: list[tuple[str, int, str]] = []
    with open(pruned_path) as fh:
        lines = fh.read().splitlines()
    for line in _polymer.mask_lines(lines):
        key = (line[21], int(line[22:26]), line[26])
        if key not in seen:
            seen.add(key)
            residues.append(key)
    return residues


def _is_hydrogen(line: str) -> bool:
    """Check if an ATOM/HETATM line is a hydrogen atom."""
    # Element symbol is in columns 77-78 (0-indexed 76-77)
    if len(line) >= 78:
        element = line[76:78].strip()
        return element in ("H", "D")
    # Fallback: check atom name (col 13-16) for hydrogen naming conventions
    atom_name = line[12:16]
    first_nonspace = atom_name.lstrip()[0] if atom_name.strip() else ""
    return first_nonspace in ("H", "D")


def process_entry(
    entry_dir: Path, suffix: str,
    masks_only: bool = False, rcsb_remap: bool = False,
) -> tuple[str, bool, str, int]:
    """Process a single entry directory.

    Creates {pdbid}_ersatz.pdb and {pdbid}_{chain}.mask in the entry dir.

    Parameters
    ----------
    entry_dir : Path
        Entry directory (e.g., .../1a/1a2z/)
    suffix : str
        Pruned file suffix (e.g., "pruned_all" or "pruned_mc")
    masks_only : bool
        Write the mask files only; do not rebuild the ersatz PDB.
    rcsb_remap : bool
        Apply _RCSB_REMAPS and _RCSB_SKIP (top2018_full).

    Returns (pdb_id, success, message, n_mask_residues).
    """
    pdb_id = entry_dir.name

    # Find input files
    fh_pdb = entry_dir / f"{pdb_id}FH.pdb"
    pruned_files = list(
        entry_dir.glob(f"{pdb_id}_*_{suffix}.pdb")
    )

    if not masks_only and not fh_pdb.exists():
        return (pdb_id, False, "missing FH file", 0)
    if not pruned_files:
        return (pdb_id, False, f"missing pruned PDB (*_{suffix}.pdb)", 0)

    try:
        # Build ersatz PDB: FH file with H stripped and headers cleaned
        # (one ersatz per entry, shared across all chains)
        ersatz_path = entry_dir / f"{pdb_id}_ersatz.pdb"
        atom_count = 0

        if not masks_only:
            _write_ersatz(fh_pdb, ersatz_path)
            atom_count = sum(1 for ln in open(ersatz_path) if ln.startswith("ATOM  "))

        # Build mask files from ALL pruned PDBs in this entry
        # (some entries have multiple chains, e.g. 5t5i has A,C,G,J,L)
        total_mask = 0
        chains: list[str] = []
        pattern = re.compile(
            rf"{re.escape(pdb_id)}_(\w)_{re.escape(suffix)}\.pdb"
        )
        for pruned_pdb in sorted(pruned_files):
            chain_match = pattern.match(pruned_pdb.name)
            if not chain_match:
                continue
            target_chain = chain_match.group(1)
            if rcsb_remap and (pdb_id, target_chain) in _RCSB_SKIP:
                continue
            chains.append(target_chain)

            mask_residues = _build_mask(pruned_pdb)
            remap = _RCSB_REMAPS.get((pdb_id, target_chain)) if rcsb_remap else None
            if remap is not None:
                mask_residues = [remap(k) for k in mask_residues]
            mask_path = entry_dir / f"{pdb_id}_{target_chain}.mask"
            with open(mask_path, "w") as fout:
                for chain_id, resnum, icode in mask_residues:
                    fout.write(f"{chain_id} {resnum} {icode}\n")
            total_mask += len(mask_residues)

        msg = (
            f"ok: chains={','.join(chains)}, "
            + ("masks only, " if masks_only else f"{atom_count} heavy atoms, ")
            + f"{total_mask} mask residues"
        )
        return (pdb_id, True, msg, total_mask)

    except Exception as exc:
        return (pdb_id, False, str(exc), 0)


def _write_ersatz(fh_pdb: Path, ersatz_path: Path) -> None:
    """Write the ersatz PDB: FH file with hydrogens and unwanted headers removed."""
    with open(fh_pdb) as fin, open(ersatz_path, "w") as fout:
        for line in fin:
            # Header filtering
            if not line.startswith(_KEEP_PREFIXES):
                continue

            # Skip hydrogen atoms
            if line.startswith(("ATOM  ", "HETATM")) and _is_hydrogen(line):
                continue

            fout.write(line)


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
            "Build ersatz PDB files from Reduce FH output. "
            "Strips hydrogens and problematic headers, generates "
            "per-chain mask files from pruned PDB files."
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
        "-s", "--suffix", type=str, default="pruned_all",
        help=(
            "Pruned file suffix to match (default: pruned_all). "
            "Use 'pruned_mc' for mainchain-filtered datasets. "
            "Matches files named {pdbid}_{chain}_{suffix}.pdb"
        ),
    )
    parser.add_argument(
        "--masks-only", action="store_true",
        help="Rebuild the .mask files only; do not rewrite the ersatz PDBs",
    )
    parser.add_argument(
        "--rcsb-remap", action="store_true",
        help="Apply the RCSB chain/numbering remaps (top2018_full only)",
    )
    args = parser.parse_args()

    if not args.src_dir.is_dir():
        print(f"Error: {args.src_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    print(f"Scanning {args.src_dir} ...")
    print(f"Pruned file suffix: {args.suffix}")
    print(f"Masks only: {args.masks_only}  RCSB remap: {args.rcsb_remap}")
    entries = collect_entries(args.src_dir)
    print(f"Found {len(entries)} entries")

    if not entries:
        print("Nothing to do.")
        return

    done = 0
    errors = 0
    total_mask = 0

    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(process_entry, entry, args.suffix,
                        args.masks_only, args.rcsb_remap): entry
            for entry in entries
        }
        for future in as_completed(futures):
            pdb_id, success, msg, n_mask = future.result()
            done += 1
            if not success:
                errors += 1
                print(f"  ERROR {pdb_id}: {msg}", file=sys.stderr)
            else:
                total_mask += n_mask
            if done % 500 == 0 or done == len(entries):
                print(
                    f"  {done}/{len(entries)} entries processed",
                    flush=True,
                )

    print(f"\nDone: {done - errors} succeeded, {errors} errors")
    print(f"Total mask residues: {total_mask}")
    print(f"Output: in-place alongside originals in {args.src_dir}")


if __name__ == "__main__":
    main()
