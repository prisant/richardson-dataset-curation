#!/usr/bin/env python3
"""List added mask residues that pydangle measures but that are absent from NEW.

Helper for the HETATM mask fix regression.  "Measured" = resname in pydangle's
_PROTEIN_NAMES.  For each such residue absent from the new filtered JSONL,
report whether pydangle produced a record at all (in the unfiltered RAW).

Usage: b1_missing_measured.py ADDED_TSV NEW_FILTERED.jsonl RAW.jsonl
"""
import csv
import json
import sys
from collections import Counter

from pydangle_biopython.measure import _PROTEIN_NAMES


def keys(path):
    out = set()
    for line in open(path):
        r = json.loads(line)
        if not r.get("_meta"):
            out.add((r["file"].split("_")[0], r["chain"], r["resnum"], r["ins"] or " "))
    return out


added = [r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")
         if r["resname"] in _PROTEIN_NAMES]
new, raw = keys(sys.argv[2]), keys(sys.argv[3])
miss = [r for r in added if (r["pdb"], r["chain"], int(r["resnum"]), r["icode"] or " ") not in new]
print(f"added & measured-type: {len(added)}  absent from new filtered: {len(miss)}")
print("  by resname:", dict(Counter(r["resname"] for r in miss)))
for r in miss:
    k = (r["pdb"], r["chain"], int(r["resnum"]), r["icode"] or " ")
    print("  ", r["pdb"], r["chain"], r["resnum"], r["resname"], "in raw" if k in raw else "NOT in raw")
