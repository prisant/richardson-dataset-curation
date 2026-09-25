#!/usr/bin/env python3
"""Fetch legacy-format DSSP files from the PDB-REDO DSSP databank.

Provides an independent DSSP reference for ``validate_dssp.py``: the
databank's legacy files are computed by the DSSP authors (NKI, DSSP 4.x)
directly from the deposited wwPDB coordinates, with none of this
pipeline's input cleaning.  Verified 2026-09-24 on top100 entries 1aac,
1ads, 1aky: CA coordinates in the databank files are identical to those
in ``run_dssp.py`` output from the deposited PDB files.

For every entry directory in SRC_DIR (two-level: prefix/pdbid/), fetches

    https://pdb-redo.eu/dssp/db/{pdbid}/legacy

and writes ``OUT_DIR/{pdbid}.dssp`` only when the response is HTTP 200
and contains a DSSP residue table.  The databank serves legacy files
only for structures that fit the classic format (e.g. not > 26 chains),
so absent entries are expected; they are recorded, never written as
empty files, and ``validate_dssp.py`` reports them as missing-dssp.

A manifest ``OUT_DIR/manifest.tsv`` records, per pdbid: status, HTTP
code, bytes, SHA-256, and the DSSP version and DATE from line 1.

Usage:
    python pipeline/fetch_pdbredo_dssp.py SRC_DIR OUT_DIR [--delay S] [--force]

Then, e.g.:
    python pipeline/validate_dssp.py top100/top100_measures.jsonl OUT_DIR
"""

from __future__ import annotations

import argparse
import hashlib
import http.client
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

URL = "https://pdb-redo.eu/dssp/db/{pdbid}/legacy"
HEADER_RE = re.compile(r"DSSP, NKI version (\S+).*DATE=(\S+)")


def collect_pdbids(src_dir: Path) -> list[str]:
    """Entry directory names (two-level: prefix/pdbid/), sorted."""
    return sorted(
        entry.name
        for prefix in src_dir.iterdir() if prefix.is_dir()
        for entry in prefix.iterdir() if entry.is_dir()
    )


RETRY_CODES = {0, 429, 500, 502, 503, 504}


def fetch_once(pdbid: str) -> tuple[int, bytes]:
    """One request. Returns (http_code, body); code 0 = transport failure."""
    req = urllib.request.Request(
        URL.format(pdbid=pdbid.lower()),
        headers={"User-Agent": "richardson-dataset-curation/fetch_pdbredo_dssp"},
    )
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            return resp.status, resp.read()
    except urllib.error.HTTPError as exc:
        return exc.code, b""
    except (urllib.error.URLError, http.client.HTTPException, TimeoutError, OSError) as exc:
        print(f"  ERROR {pdbid}: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 0, b""


def fetch(pdbid: str, retries: int = 3, backoff: float = 10.0) -> tuple[int, bytes]:
    """Fetch with retries on transient failures (5xx, 429, transport errors)."""
    code, data = fetch_once(pdbid)
    for attempt in range(1, retries + 1):
        if code not in RETRY_CODES:
            break
        time.sleep(backoff * attempt)
        code, data = fetch_once(pdbid)
    return code, data


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("src_dir", type=Path, help="Dataset directory (prefix/pdbid/)")
    parser.add_argument("out_dir", type=Path, help="Directory for {pdbid}.dssp + manifest.tsv")
    parser.add_argument("--delay", type=float, default=0.5,
                        help="Seconds between requests (default: 0.5)")
    parser.add_argument("-f", "--force", action="store_true",
                        help="Re-fetch files that already exist")
    args = parser.parse_args()

    if not args.src_dir.is_dir():
        print(f"Error: {args.src_dir} is not a directory", file=sys.stderr)
        return 2
    args.out_dir.mkdir(parents=True, exist_ok=True)

    pdbids = collect_pdbids(args.src_dir)
    print(f"{len(pdbids)} entries in {args.src_dir}")

    rows = []
    counts = {"ok": 0, "cached": 0, "absent": 0, "no-table": 0, "error": 0}
    try:
        _fetch_all(pdbids, args, rows, counts)
    finally:
        # Always write the manifest, even if the loop dies unexpectedly.
        with open(args.out_dir / "manifest.tsv", "w", encoding="utf-8") as fh:
            fh.write("pdbid\tstatus\thttp\tbytes\tsha256\tdssp_version\tdssp_date\n")
            for row in rows:
                fh.write("\t".join(row) + "\n")
        print(f"manifest: {len(rows)}/{len(pdbids)} entries recorded")

    print("Done: " + ", ".join(f"{k}={v}" for k, v in counts.items()))
    return 1 if counts["error"] or counts["no-table"] else 0


def _fetch_all(pdbids, args, rows, counts) -> None:
    for n, pdbid in enumerate(pdbids, 1):
        out = args.out_dir / f"{pdbid}.dssp"
        if out.exists() and not args.force:
            data, code, status = out.read_bytes(), 200, "cached"
        else:
            code, data = fetch(pdbid)
            if code == 200 and b"  #  RESIDUE" in data:
                out.write_bytes(data)
                status = "ok"
            elif code == 404:
                status = "absent"
            elif code == 200:
                status = "no-table"
            else:
                status = "error"
            time.sleep(args.delay)
        counts[status] += 1

        version = date = ""
        if data:
            m = HEADER_RE.search(data.split(b"\n", 1)[0].decode("latin-1"))
            if m:
                version, date = m.groups()
        sha = hashlib.sha256(data).hexdigest() if status in ("ok", "cached") else ""
        rows.append((pdbid, status, str(code), str(len(data)), sha, version, date))
        if n % 100 == 0:
            print(f"  {n}/{len(pdbids)}")


if __name__ == "__main__":
    sys.exit(main())
