#!/usr/bin/env python3
"""Download full PDB files from RCSB for the top8000 dataset.

Reads a list of PDB IDs and downloads each from RCSB into the
two-level directory structure (prefix/pdbid/pdbid.pdb).

Usage:
    python top8000/download_pdbs.py DOWNLOAD_LIST DST_DIR [-j JOBS]

Example:
    python top8000/download_pdbs.py \
        ~/Desktop/top8000_pdbs_hom70/need_download.txt \
        ~/Desktop/top8000_pdbs_hom70/ -j 8
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import urlretrieve


def download_pdb(pdb_id: str, dst_dir: Path) -> tuple[str, bool, str]:
    """Download a single PDB file from RCSB.

    Returns (pdb_id, success, message).
    """
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
        # Try .pdb.gz fallback
        if e.code == 404:
            return (pdb_id, False, f"404 not found on RCSB")
        return (pdb_id, False, f"HTTP {e.code}")
    except Exception as exc:
        return (pdb_id, False, str(exc))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download full PDB files from RCSB"
    )
    parser.add_argument(
        "download_list", type=Path,
        help="File with one PDB ID per line",
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

    if not args.download_list.exists():
        print(f"Error: {args.download_list} not found", file=sys.stderr)
        sys.exit(1)

    pdb_ids = [
        line.strip() for line in open(args.download_list)
        if line.strip()
    ]
    print(f"PDB IDs to download: {len(pdb_ids)}")

    done = 0
    errors = 0
    skipped = 0
    failed_ids: list[str] = []

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(download_pdb, pdb_id, args.dst_dir): pdb_id
            for pdb_id in pdb_ids
        }
        for future in as_completed(futures):
            pdb_id, success, msg = future.result()
            done += 1
            if not success:
                errors += 1
                failed_ids.append(pdb_id)
                print(f"  FAIL {pdb_id}: {msg}", file=sys.stderr)
            elif "already exists" in msg:
                skipped += 1
            if done % 200 == 0 or done == len(pdb_ids):
                print(
                    f"  {done}/{len(pdb_ids)} processed "
                    f"({errors} errors)",
                    flush=True,
                )

    print(f"\nDone: {done - errors - skipped} downloaded, "
          f"{skipped} already existed, {errors} errors")

    if failed_ids:
        fail_path = args.dst_dir / "download_failures.txt"
        with open(fail_path, "w") as f:
            for pid in sorted(failed_ids):
                f.write(f"{pid}\n")
        print(f"Failed IDs saved to: {fail_path}")


if __name__ == "__main__":
    main()
