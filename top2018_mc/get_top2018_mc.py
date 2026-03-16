#!/usr/bin/env python3
"""Download full PDB files from RCSB for the top2018 MC-filtered dataset.

Sets up the two-level directory structure and downloads full PDB files
for each entry in the MC-filtered chain list. The pruned MC files
({pdbid}_{chain}_pruned_mc.pdb) must be obtained separately from the
Richardson Lab top2018 Zenodo deposit (doi:10.5281/zenodo.5115233) and
placed in the corresponding entry directories.

Usage:
    python get_top2018_mc.py [-j JOBS]

The chain list (top2018mc_chain_list.txt) must be in the same directory.
Each line has the format: {pdbid}_{chain} (e.g., 19hc_A).
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import urlretrieve


def download_pdb(pdb_id: str, dst_dir: Path) -> tuple[str, bool, str]:
    """Download a single PDB file from RCSB."""
    prefix = pdb_id[:2]
    entry_dir = dst_dir / prefix / pdb_id
    entry_dir.mkdir(parents=True, exist_ok=True)
    dst_path = entry_dir / f"{pdb_id}.pdb"

    if dst_path.exists() and dst_path.stat().st_size > 0:
        return (pdb_id, True, "already exists")

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        urlretrieve(url, dst_path)
        size_kb = dst_path.stat().st_size / 1024
        return (pdb_id, True, f"ok ({size_kb:.0f} KB)")
    except HTTPError as e:
        return (pdb_id, False, f"HTTP {e.code}")
    except Exception as exc:
        return (pdb_id, False, str(exc))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download full PDB files for top2018 MC-filtered dataset"
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=8,
        help="Number of parallel download threads (default: 8)",
    )
    args = parser.parse_args()

    # Read chain list
    script_dir = Path(__file__).parent
    chain_list = script_dir / "top2018mc_chain_list.txt"
    if not chain_list.exists():
        print(f"Error: {chain_list} not found", file=sys.stderr)
        sys.exit(1)

    pdb_ids = set()
    with open(chain_list) as f:
        for line in f:
            entry = line.strip()
            if entry:
                pdb_ids.add(entry.split("_")[0].lower())

    pdb_ids_sorted = sorted(pdb_ids)
    print(f"Unique PDB IDs: {len(pdb_ids_sorted)}")

    dst_dir = Path("top2018_pdbs_mc_filtered_hom70")
    dst_dir.mkdir(exist_ok=True)

    done = 0
    errors = 0

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(download_pdb, pid, dst_dir): pid
            for pid in pdb_ids_sorted
        }
        for future in as_completed(futures):
            pdb_id, success, msg = future.result()
            done += 1
            if not success:
                errors += 1
                print(f"  FAIL {pdb_id}: {msg}", file=sys.stderr)
            if done % 500 == 0 or done == len(pdb_ids_sorted):
                print(f"  {done}/{len(pdb_ids_sorted)} processed", flush=True)

    print(f"\nDone: {done - errors} downloaded, {errors} errors")
    print(
        "\nNote: You must also place the pruned MC files "
        "({pdbid}_{chain}_pruned_mc.pdb) from the Richardson Lab "
        "top2018 Zenodo deposit (doi:10.5281/zenodo.5115233) in "
        "the corresponding entry directories."
    )


if __name__ == "__main__":
    main()
