#!/usr/bin/env python3
"""Count residues with backbone alternate conformations and their DSSP fate.

Finding (2026-09-24, top500 vs PDB-REDO DSSP databank): mkdssp (4.2.2 and
4.6.1) omits residues whose backbone atoms carry altloc codes when it reads
PDB-format input (e.g. 2hmz A-D 21, 1ejg A22/A25, 1din A123); the same
structures read from mmCIF keep them.  This audit measures the scale.

For each entry directory (two-level prefix/pdbid/) it reads the original
``{pdbid}.pdb`` and flags every residue (chain, resnum, icode) that has an
ATOM/HETATM record for N, CA, C or O with a non-blank altloc (column 17).
For each flagged residue it reports whether it is present in the entry's
``{pdbid}.dssp`` (run_dssp.py output) and whether it is present in the
dataset JSONL (i.e. whether the omission reaches the published data).

Usage:
    python pipeline/audit_backbone_altlocs.py SRC_DIR MEASURES.jsonl[.gz] [-o REPORT.md]
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from run_dssp import collect_entries  # noqa: E402

BACKBONE = {"N", "CA", "C", "O"}


def backbone_altloc_residues(pdb: Path) -> set[tuple[str, int, str]]:
    found: set[tuple[str, int, str]] = set()
    with open(pdb, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break  # model 1 only, as validate_dssp.py
            if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 27:
                continue
            if line[16] == " " or line[12:16].strip() not in BACKBONE:
                continue
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            found.add((line[21], resnum, line[26].strip()))
    return found


def dssp_keys(dssp: Path) -> set[tuple[str, int, str]] | None:
    if not dssp.exists():
        return None
    keys: set[tuple[str, int, str]] = set()
    in_table = False
    for line in dssp.read_text(errors="replace").splitlines():
        if line.startswith("  #  RESIDUE"):
            in_table = True
            continue
        if not in_table or len(line) < 17 or line[13] == "!":
            continue
        try:
            keys.add((line[11], int(line[5:10]), line[10].strip()))
        except ValueError:
            continue
    return keys


def jsonl_keys(path: Path) -> dict[str, set[tuple[str, int, str]]]:
    op = gzip.open if path.suffix == ".gz" else open
    out: dict[str, set[tuple[str, int, str]]] = {}
    with op(path, "rt", encoding="utf-8") as fh:
        for raw in fh:
            rec = json.loads(raw)
            if rec.get("_meta") or rec.get("model", 1) != 1:
                continue
            pid = Path(rec["file"]).name.split(".", 1)[0].split("_", 1)[0]
            out.setdefault(pid, set()).add(
                (rec["chain"], int(rec["resnum"]), (rec.get("ins") or "").strip()))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("src_dir", type=Path)
    ap.add_argument("jsonl", type=Path)
    ap.add_argument("-o", "--output", type=Path)
    a = ap.parse_args()

    jk = jsonl_keys(a.jsonl)
    c: Counter[str] = Counter()
    examples: list[str] = []
    for entry in collect_entries(a.src_dir):
        pid = entry.name
        pdb = entry / f"{pid}.pdb"
        if not pdb.exists():
            c["entries without .pdb"] += 1
            continue
        c["entries scanned"] += 1
        alt = backbone_altloc_residues(pdb)
        if not alt:
            continue
        c["entries with backbone altlocs"] += 1
        dk = dssp_keys(entry / f"{pid}.dssp")
        in_json = jk.get(pid, set())
        for k in alt:
            c["backbone-altloc residues"] += 1
            inj = k in in_json
            c["  ... of which in JSONL"] += inj
            if dk is None:
                c["  ... no .dssp file for entry"] += 1
            elif k in dk:
                c["  ... present in .dssp"] += 1
            else:
                c["  ... MISSING from .dssp"] += 1
                c["  ... MISSING from .dssp AND in JSONL"] += inj
                if inj and len(examples) < 20:
                    examples.append(f"{pid} {k[0]}{k[1]}{k[2]}")

    lines = [f"# Backbone-altloc audit: `{a.jsonl.name}`", "",
             f"- Source dir: `{a.src_dir}`", "", "| Count | Value |", "|---|---:|"]
    lines += [f"| {k} | {v:,} |" for k, v in c.items()]
    lines += ["", "Examples (missing from .dssp, present in JSONL): " + ", ".join(examples)]
    text = "\n".join(lines) + "\n"
    print(text)
    if a.output:
        a.output.write_text(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
