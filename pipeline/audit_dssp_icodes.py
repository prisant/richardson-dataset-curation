#!/usr/bin/env python3
"""Find DSSP residue keys that do not exist in the deposited PDB file.

Finding (2026-09-24): mkdssp 4.2.2 labels 2hft A165 as ``165A`` although
the deposited file has no insertion code there; 4.6.1 and the PDB-REDO
databank write ``165``.  If pydangle's own mkdssp call does the same,
the residue's lookup by (chain, resnum, icode) fails and the JSONL gets
``dssp = null`` -- invisible to validate_dssp.py when the reference
.dssp files come from the same mkdssp.

For every entry: residue keys (chain, resnum, icode) from ATOM/HETATM
records of model 1 in ``{pdbid}.pdb``; keys from ``{pdbid}.dssp`` (the
reference made by run_dssp.py).  A DSSP key absent from the PDB is
"phantom"; it is "icode-shifted" when the PDB has the same (chain,
resnum) with a different insertion code.  For each phantom key the
report gives the JSONL record(s) for that (chain, resnum): insertion
code, ``dssp`` value, and the code in ALT_DIR (e.g. DSSP 4.6.1 output)
for the PDB's real key.

Usage:
    python pipeline/audit_dssp_icodes.py SRC_DIR MEASURES.jsonl[.gz] [--alt ALT_DIR] [-o REPORT.md]
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from pydangle_biopython.dssp import parse_dssp_output  # noqa: E402
from run_dssp import collect_entries  # noqa: E402


def pdb_keys(pdb: Path) -> set[tuple[str, int, str]]:
    """Residue keys of model 1, with blank chain IDs filled as run_dssp.py does
    (first non-blank ATOM chain ID, else 'A')."""
    raw: list[tuple[str, int, str]] = []
    default_chain = None
    with open(pdb, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 27:
                if default_chain is None and line.startswith("ATOM  ") and line[21] != " ":
                    default_chain = line[21]
                try:
                    raw.append((line[21], int(line[22:26]), line[26].strip()))
                except ValueError:
                    pass
    fill = default_chain or "A"
    return {((fill if ch == " " else ch), n, ic) for ch, n, ic in raw}


def jsonl_index(path: Path) -> dict[str, dict[tuple[str, int], list[tuple[str, object]]]]:
    op = gzip.open if path.suffix == ".gz" else open
    out: dict[str, dict[tuple[str, int], list[tuple[str, object]]]] = {}
    with op(path, "rt", encoding="utf-8") as fh:
        for raw in fh:
            rec = json.loads(raw)
            if rec.get("_meta") or rec.get("model", 1) != 1:
                continue
            pid = Path(rec["file"]).name.split(".", 1)[0].split("_", 1)[0]
            key = (rec["chain"], int(rec["resnum"]))
            out.setdefault(pid, {}).setdefault(key, []).append(
                ((rec.get("ins") or "").strip(), rec.get("dssp")))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("src_dir", type=Path)
    ap.add_argument("jsonl", type=Path)
    ap.add_argument("--alt", type=Path, help="Directory of alternate {pdbid}.dssp (e.g. DSSP 4.6.1)")
    ap.add_argument("-o", "--output", type=Path)
    a = ap.parse_args()

    jx = jsonl_index(a.jsonl)
    c: Counter[str] = Counter()
    rows = []
    for entry in collect_entries(a.src_dir):
        pid = entry.name
        pdb, dssp = entry / f"{pid}.pdb", entry / f"{pid}.dssp"
        if not pdb.exists() or not dssp.exists():
            continue
        c["entries checked"] += 1
        pk = pdb_keys(pdb)
        dk = parse_dssp_output(dssp.read_text(errors="replace"))
        phantom = sorted(k for k in dk if k not in pk)
        if not phantom:
            continue
        c["entries with phantom keys"] += 1
        alt = None
        if a.alt and (a.alt / f"{pid}.dssp").exists():
            alt = parse_dssp_output((a.alt / f"{pid}.dssp").read_text(errors="replace"))
        for k in phantom:
            c["phantom DSSP keys"] += 1
            real = sorted(r for r in pk if r[0] == k[0] and r[1] == k[1])
            shifted = bool(real)
            c["  icode-shifted (PDB has same chain+resnum)"] += shifted
            js = jx.get(pid, {}).get((k[0], k[1]), [])
            in_json = bool(js)
            c["  ... with JSONL record at chain+resnum"] += in_json
            nulls = [j for j in js if j[1] is None]
            c["  ... JSONL dssp null there"] += bool(nulls)
            alt_codes = [alt.get(r, "ABSENT") for r in real] if alt is not None else []
            c["  ... alt has a real code for PDB key"] += any(x not in (None, "ABSENT") for x in alt_codes)
            c["  ... JSONL null AND alt has real code (lost assignment)"] += bool(nulls) and any(
                x not in (None, "ABSENT") for x in alt_codes)
            rows.append((pid, f"{k[0]}{k[1]}{k[2]}", dk[k],
                         ",".join(f"{r[0]}{r[1]}{r[2]}" for r in real) or "-",
                         ";".join(f"ins={j[0]!r} dssp={j[1]!r}" for j in js) or "-",
                         ",".join(repr(x) for x in alt_codes) or "-"))

    lines = [f"# Phantom DSSP residue keys: `{a.jsonl.name}`", "", f"- Source dir: `{a.src_dir}`",
             f"- Alt dir: `{a.alt}`", "", "| Count | Value |", "|---|--:|"]
    lines += [f"| {k} | {v:,} |" for k, v in c.items()]
    lines += ["", "| Entry | DSSP key | DSSP code | PDB key(s) same chain+resnum | JSONL record(s) | Alt code(s) for PDB key |",
              "|---|---|---|---|---|---|"]
    lines += ["| " + " | ".join(str(v) for v in r) + " |" for r in rows]
    text = "\n".join(lines) + "\n"
    print("\n".join(lines[:16 + len(c)]))
    if a.output:
        a.output.write_text(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
