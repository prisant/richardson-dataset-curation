#!/usr/bin/env python3
"""Classify per-entry DSSP differences by proximity to dropped backbone-altloc residues.

Question (2026-09-25): are the residual differences between our pipeline's
DSSP 4.6.1 output (PDB-format input) and the PDB-REDO DSSP databank
(4.6.1, mmCIF input) all explained by residues with backbone alternate
conformations that mkdssp omits from PDB-format input?

For every entry present in both OURS and REF, a "difference" is a residue
whose code differs or that is present in only one of the two.  A
difference is "near" if it lies in the same chain within +/-WINDOW
residue numbers of a backbone-altloc residue (from the original PDB)
that is absent from OURS.  Entries are classified:

  explained     all differences near a dropped altloc residue
  partial       some near, some not
  unexplained   differences but no dropped altloc residue nearby

Usage:
    python pipeline/diagnose_dssp_diffs_altloc.py OURS_DIR REF_DIR SRC_DIR [--window N] [-o REPORT.md]
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from audit_backbone_altlocs import backbone_altloc_residues  # noqa: E402
from pydangle_biopython.dssp import parse_dssp_output  # noqa: E402
from run_dssp import collect_entries  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("ours", type=Path)
    ap.add_argument("ref", type=Path)
    ap.add_argument("src_dir", type=Path)
    ap.add_argument("--window", type=int, default=4)
    ap.add_argument("-o", "--output", type=Path)
    a = ap.parse_args()

    entries = {e.name: e for e in collect_entries(a.src_dir)}
    ours = {p.stem: p for p in a.ours.glob("*.dssp")}
    ref = {p.stem: p for p in a.ref.glob("*.dssp")}
    classes: Counter[str] = Counter()
    diff_total: Counter[str] = Counter()
    rows = []
    for pid in sorted(set(ours) & set(ref)):
        o = parse_dssp_output(ours[pid].read_text(errors="replace"))
        r = parse_dssp_output(ref[pid].read_text(errors="replace"))
        diffs = [k for k in o.keys() | r.keys() if o.get(k, "ABSENT") != r.get(k, "ABSENT")]
        if not diffs:
            continue
        pdb = entries[pid] / f"{pid}.pdb" if pid in entries else None
        alt = backbone_altloc_residues(pdb) if pdb and pdb.exists() else set()
        dropped = {k for k in alt if k not in o}
        near = [k for k in diffs
                if any(d[0] == k[0] and abs(d[1] - k[1]) <= a.window for d in dropped)]
        if len(near) == len(diffs):
            cls = "explained"
        elif near:
            cls = "partial"
        else:
            cls = "unexplained"
        classes[cls] += 1
        diff_total[cls] += len(diffs)
        far = sorted(set(diffs) - set(near))
        rows.append((cls, pid, len(diffs), len(near), len(dropped),
                     ", ".join(f"{k[0]}{k[1]}{k[2]}:{o.get(k, 'ABSENT')}/{r.get(k, 'ABSENT')}"
                               for k in far[:8])))

    lines = [f"# DSSP differences vs dropped backbone altlocs (window ±{a.window})", "",
             f"- ours: `{a.ours}`", f"- ref: `{a.ref}`", "",
             "| Class | Entries | Differing residues |", "|---|--:|--:|"]
    lines += [f"| {c} | {classes[c]} | {diff_total[c]} |" for c in ("explained", "partial", "unexplained")]
    lines += ["", "| Class | Entry | Diffs | Near dropped altloc | Dropped altloc residues | Far diffs (ours/ref, first 8) |",
              "|---|---|--:|--:|--:|---|"]
    order = {"unexplained": 0, "partial": 1, "explained": 2}
    for row in sorted(rows, key=lambda x: (order[x[0]], -x[2])):
        lines.append("| " + " | ".join(str(v) for v in row) + " |")
    text = "\n".join(lines) + "\n"
    print(text)
    if a.output:
        a.output.write_text(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
