"""Scratch diagnostic: compare per-residue DSSP codes between two .dssp directories.

Matches files by {pdbid}.dssp and residues by (chain, resnum, icode) using
pydangle's parse_dssp_output (same parser as validate_dssp.py).

Usage: python compare_dssp_dirs.py DIR_A DIR_B [LABEL_A LABEL_B]
"""
import sys
from collections import Counter
from pathlib import Path

from pydangle_biopython.dssp import parse_dssp_output


def index(d):
    return {p.stem: p for p in Path(d).rglob("*.dssp")}


def load(p):
    return parse_dssp_output(p.read_text(errors="replace"))


a_dir, b_dir = sys.argv[1], sys.argv[2]
la, lb = (sys.argv[3], sys.argv[4]) if len(sys.argv) > 4 else ("A", "B")
ia, ib = index(a_dir), index(b_dir)
both = sorted(set(ia) & set(ib))
print(f"files: {la}={len(ia)} {lb}={len(ib)} both={len(both)}")
print(f"  only in {la}: {sorted(set(ia) - set(ib))[:20]}")
print(f"  only in {lb}: {sorted(set(ib) - set(ia))[:20]}")
same = diff = only_a = only_b = 0
pairs = Counter()
per_entry = Counter()
for pid in both:
    a, b = load(ia[pid]), load(ib[pid])
    for k in a.keys() | b.keys():
        if k not in b:
            only_a += 1
        elif k not in a:
            only_b += 1
        elif a[k] == b[k]:
            same += 1
        else:
            diff += 1
            pairs[(a[k], b[k])] += 1
            per_entry[pid] += 1
print(f"residues: same={same:,} different={diff:,} only-{la}={only_a:,} only-{lb}={only_b:,}")
print(f"entries with any difference: {len(per_entry)}  {per_entry.most_common(10)}")
for (x, y), n in pairs.most_common(12):
    print(f"  {la} {x!r:>5} vs {lb} {y!r:<5} {n:,}")
