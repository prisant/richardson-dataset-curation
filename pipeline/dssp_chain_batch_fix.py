"""DSSP chain-batching fallback for structures with >26 case-distinguished chain IDs.

Background — known mkdssp 4.x edge case
---------------------------------------

The classic DSSP output format reserves a single column (PDB column 22 of
each ATOM record) for the chain identifier.  mkdssp 4.x will silently cap
its legacy DSSP output at the first 26 distinct chain IDs it encounters
(reading uppercase before lowercase).  Any further chains are emitted
with a stderr warning ::

    This file contains data that won't fit in the original DSSP format

and a returncode of 1 — but the partial .dssp file is written to disk
**with the surplus chains entirely absent**.  No record at all of the
dropped chains' residues is preserved.  This is an mkdssp 4.x format
limitation, not a chain-naming issue.

In the top2018_full / top2018_mc datasets, 15 unique structures hit this
limit (chain counts 28 to 60).  Without the fallback in this module they
silently lose between 2 and 34 chains' worth of residues from the
ground-truth .dssp the validation gate consumes.

What this module does
---------------------

When ``_clean_and_run_dssp`` (or its caller) reports the
"won't fit in the original DSSP format" stderr, this module:

1. Splits the structure's chains into ordered batches of ≤26 chains.
2. Builds a per-batch cleaned PDB (HEADER + CRYST1 + SEQRES/ATOM/TER/HETATM
   for the batch's chains only, matching ``_clean_and_run_dssp``'s own
   filter rules so REMARK 800 / SITE references don't trip mkdssp).
3. Runs mkdssp on each batch — every batch fits the 26-chain cap so
   mkdssp returns RC=0 with no warnings.
4. Concatenates the per-batch residue-table records into a single
   ``.dssp`` file with renumbered global sequential indices and BP1/BP2
   references shifted to match.

Information loss at batch boundaries
------------------------------------

Cross-batch hydrogen bonds are NOT detected by per-batch mkdssp runs.
The per-residue 4-column H-bond partner offsets+energies will lose any
H-bond whose partner residue lies in a different batch.  Empirically on
4fmh (20 chains, 13/7 split, an inter-chain β-sheet structure), this
causes:

- 0.4% of residues to lose their secondary-structure-code at the strand
  boundary that crosses the split (e.g. an "E" reverts to coil because
  its sheet partner is in the other batch).
- ~4% of residues to have at least one altered H-bond partner column.

For homo-oligomers where the inter-batch interface is small or absent,
the loss is correspondingly smaller (1nfv 8/8 split: 0 sec-struct
mismatches, 0.3% H-bond column changes).

This loss is documented per dataset in ``*_issues.md`` and per
structure in the JSON output of this module.

What is fully preserved
-----------------------

- All chains, all residues (no chain dropping).
- Per-residue secondary-structure code for any residue whose H-bond
  partners are all in the same batch (vast majority).
- Per-residue ACC, TCO, KAPPA, ALPHA, PHI, PSI, X-CA/Y-CA/Z-CA — these
  are local geometry, computed entirely within the residue's own batch.

What is NOT fully preserved
---------------------------

- H-bond partner offset+energy columns for residues with cross-batch
  partners.
- Bridge/sheet labeling (cols 23-25) when a sheet spans the batch
  boundary — labels will diverge between full and batched runs.

Authority
---------

The 26-chain cap is empirically derived (mkdssp 4.2.2 NKI build,
2026-04-27 testing on 5b66, 5v2c and other top2018 entries).  All 15
affected entries in top2018_full + top2018_mc have >26 distinct
case-distinguished chain IDs; all entries with >26 chains in either
dataset hit this cap; the biconditional was confirmed by direct
``_clean_and_run_dssp`` invocation on every candidate.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

CHAIN_BATCH_SIZE = 26
WONT_FIT_MARKER = "won't fit in the original DSSP format"


def find_mkdssp() -> str | None:
    return shutil.which("mkdssp")


def chains_in_pdb(pdb_path: Path) -> list[str]:
    """Return distinct case-sensitive chain IDs from ATOM records in file order."""
    seen, ordered = set(), []
    with open(pdb_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("ATOM") and len(line) > 21:
                c = line[21]
                if c not in seen:
                    seen.add(c)
                    ordered.append(c)
    return ordered


def write_clean_batch_pdb(
    src_pdb: Path,
    out_pdb: Path,
    keep_chains: set[str],
) -> None:
    """Write a cleaned PDB with the same record filter as ``_clean_and_run_dssp``,
    restricted to ``keep_chains``.

    Filters: keep HEADER/CRYST1; keep SEQRES if column-12 chain in keep_chains;
    keep ATOM/TER/HETATM if column-22 chain in keep_chains.
    """
    with open(src_pdb, encoding="utf-8") as fin, open(out_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("HEADER", "CRYST1")):
                fout.write(line)
            elif line.startswith("SEQRES"):
                if len(line) > 11 and line[11] in keep_chains:
                    fout.write(line)
            elif line.startswith(("ATOM  ", "TER   ", "HETATM")):
                if len(line) > 21 and line[21] in keep_chains:
                    fout.write(line)
        fout.write("END\n")


def split_residue_table(dssp_path: Path) -> tuple[list[str], list[str], list[str]]:
    """Split a .dssp file into (header_lines, residue_table_lines, trailing_lines).

    The residue table starts at the line that begins with '  #  RESIDUE'.
    Everything before is header, everything after is residue records.
    """
    header, body = [], []
    in_body = False
    with open(dssp_path, encoding="utf-8") as fh:
        for line in fh:
            if not in_body:
                if line.startswith("  #  RESIDUE"):
                    in_body = True
                    header.append(line)
                else:
                    header.append(line)
            else:
                body.append(line)
    return header, body, []


def renumber_record(line: str, new_idx: int, bp_shift: int) -> str:
    """Rewrite the global sequential index (cols 1-5) and shift BP1/BP2 (cols 27-33).

    The DSSP record format:
      cols 1-5:   global sequential index (right-justified, 5 wide)
      cols 6-10:  PDB residue number
      col 11:     PDB insertion code
      col 12:     chain ID
      col 14:     amino acid 1-letter
      cols 17-25: STRUCTURE field
      cols 27-29: BP1 (3 wide, right-justified)
      cols 31-33: BP2 (3 wide, right-justified)
      ...

    Chain-break records (col 14 == '!') are not residues; pass them through unchanged.
    """
    if len(line) < 33 or (len(line) > 13 and line[13] == "!"):
        return line
    new_idx_str = f"{new_idx:>5}"
    rest_after_idx = line[5:]

    # BP1 in cols 26-29 (0-indexed 26-29), BP2 in cols 30-33 (0-indexed 30-33).
    # Use 0-indexed slices: cols 27-29 = indices 26..29 (4 chars including leading)
    # Actually classic DSSP: col 27-30 BP1, col 31-34 BP2. We'll be conservative:
    # parse cols 26-29 and 30-33 as BP1/BP2 (4 chars each).
    new_line = new_idx_str + rest_after_idx
    if bp_shift:
        try:
            bp1_str = new_line[26:30]
            bp2_str = new_line[30:34]
            bp1 = int(bp1_str.strip()) if bp1_str.strip() else 0
            bp2 = int(bp2_str.strip()) if bp2_str.strip() else 0
            if bp1 > 0:
                bp1 += bp_shift
            if bp2 > 0:
                bp2 += bp_shift
            bp1_new = f"{bp1:>4}" if bp1 else "   0"
            bp2_new = f"{bp2:>4}" if bp2 else "   0"
            new_line = new_line[:26] + bp1_new + bp2_new + new_line[34:]
        except ValueError:
            pass
    return new_line


def run_dssp_chain_batched(
    pdb_path: Path,
    dssp_path: Path,
    mkdssp: str,
    batch_size: int = CHAIN_BATCH_SIZE,
) -> dict:
    """Run mkdssp in chain-batches and concatenate the residue tables.

    Returns a stats dict with per-batch counts and any warnings.
    Writes the merged .dssp to ``dssp_path``.
    """
    chains = chains_in_pdb(pdb_path)
    n_chains = len(chains)
    batches = [chains[i:i + batch_size] for i in range(0, n_chains, batch_size)]
    stats = {
        "pdb_id": pdb_path.stem,
        "n_chains_input": n_chains,
        "chains_input": "".join(chains),
        "n_batches": len(batches),
        "batches": [],
        "warnings": [],
    }

    with tempfile.TemporaryDirectory(prefix=f"dssp_batch_{pdb_path.stem}_") as tmp:
        tmp = Path(tmp)
        per_batch_records: list[list[str]] = []
        first_header: list[str] | None = None

        for i, batch in enumerate(batches):
            batch_pdb = tmp / f"batch{i}.pdb"
            batch_dssp = tmp / f"batch{i}.dssp"
            write_clean_batch_pdb(pdb_path, batch_pdb, set(batch))
            res = subprocess.run(
                [mkdssp, str(batch_pdb), str(batch_dssp)],
                capture_output=True, text=True, timeout=300, check=False,
            )
            if res.returncode != 0:
                stats["warnings"].append(
                    f"batch {i} (chains={''.join(batch)}) mkdssp failed: "
                    f"{res.stderr.strip().split(chr(10))[-1] if res.stderr else 'unknown'}"
                )
                stats["batches"].append({
                    "i": i, "chains": "".join(batch),
                    "rc": res.returncode, "n_residues": 0, "succeeded": False,
                })
                continue
            header, body, _ = split_residue_table(batch_dssp)
            if first_header is None:
                first_header = header
            per_batch_records.append(body)
            n_res = sum(
                1 for L in body
                if len(L) > 13 and L[13] != "!" and L.strip()
            )
            stats["batches"].append({
                "i": i, "chains": "".join(batch),
                "rc": 0, "n_residues": n_res, "succeeded": True,
            })

        if first_header is None:
            return {**stats, "ok": False, "reason": "all batches failed"}

        # Concatenate, renumbering global index and shifting BP1/BP2.
        merged_body: list[str] = []
        global_idx = 0
        running_offset = 0
        for body in per_batch_records:
            batch_residue_count = 0
            for L in body:
                if len(L) > 13 and L[13] == "!":
                    merged_body.append(L)
                    continue
                if not L.strip() or len(L) < 33:
                    continue
                global_idx += 1
                merged_body.append(renumber_record(L, global_idx, running_offset))
                batch_residue_count += 1
            running_offset += batch_residue_count

        with open(dssp_path, "w") as fout:
            fout.writelines(first_header)
            fout.writelines(merged_body)

        stats["n_residues_output"] = global_idx
        stats["ok"] = True
        return stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--src-dir", required=True,
                    help="Dataset PDB directory (e.g. top2018_pdbs_mc_filtered_hom70)")
    ap.add_argument("--ids", required=True,
                    help="Comma-separated PDB IDs to reprocess (the >26-chain set)")
    ap.add_argument("--report", required=True,
                    help="JSON output file with per-structure stats")
    args = ap.parse_args()

    mkdssp = find_mkdssp()
    if mkdssp is None:
        sys.exit("mkdssp not found")
    src = Path(args.src_dir)
    ids = [s.strip() for s in args.ids.split(",") if s.strip()]
    report = []
    for pdb_id in ids:
        prefix = pdb_id[:2]
        entry_dir = src / prefix / pdb_id
        pdb_path = entry_dir / f"{pdb_id}.pdb"
        dssp_path = entry_dir / f"{pdb_id}.dssp"
        if not pdb_path.exists():
            print(f"[skip] {pdb_id}: missing {pdb_path}", file=sys.stderr)
            continue
        stats = run_dssp_chain_batched(pdb_path, dssp_path, mkdssp)
        report.append(stats)
        print(
            f"{pdb_id}: {stats.get('n_chains_input')} chains in {stats.get('n_batches')} batches → "
            f"{stats.get('n_residues_output', 0)} residues, ok={stats.get('ok')}",
            file=sys.stderr,
        )

    with open(args.report, "w") as f:
        json.dump(report, f, indent=2)
    print(f"# wrote {args.report}", file=sys.stderr)


if __name__ == "__main__":
    main()
