#!/usr/bin/env python3
"""Download full PDB files from RCSB for the top500 dataset.

Reads the top500 chain list CSV and downloads each unique PDB structure
into a two-level directory structure (prefix/pdbid/pdbid.pdb).

Usage:
    python top500/download_pdbs.py top500/top500_chain_list.csv DST_DIR [-j JOBS]

Example:
    python top500/download_pdbs.py \
        top500/top500_chain_list.csv \
        top500_pdbs/ -j 8
"""

from __future__ import annotations

import argparse
import csv
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
        description="Download full PDB files from RCSB for the top500 dataset"
    )
    parser.add_argument(
        "chain_list", type=Path,
        help="CSV file with columns pdb_id,chain",
    )
    parser.add_argument(
        "dst_dir", type=Path,
        help="Destination directory (two-level structure)",
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=8,
        help="Number of parallel download threads (default: 8)",
    )
    args = parser.parse_args()

    # Extract unique PDB IDs from chain list
    pdb_ids = set()
    with open(args.chain_list) as f:
        reader = csv.DictReader(f)
        for row in reader:
            pdb_ids.add(row["pdb_id"].lower().strip())

    pdb_ids = sorted(pdb_ids)
    print(f"PDB IDs to download: {len(pdb_ids)}")

    done = 0
    errors = 0

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(download_pdb, pid, args.dst_dir): pid
            for pid in pdb_ids
        }
        for future in as_completed(futures):
            pdb_id, success, msg = future.result()
            done += 1
            if not success:
                errors += 1
                print(f"  FAIL {pdb_id}: {msg}", file=sys.stderr)
            if done % 100 == 0 or done == len(pdb_ids):
                print(f"  {done}/{len(pdb_ids)} processed", flush=True)

    print(f"\nDone: {done - errors} downloaded, {errors} errors")


if __name__ == "__main__":
    main()
