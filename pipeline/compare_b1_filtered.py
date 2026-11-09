#!/usr/bin/env python3
"""Compare filtered JSONL before and after the HETATM mask fix (register B1).

OLD and NEW are the same pydangle output filtered by filter_pruned_residues.py
with the old and the new masks.  Expected:

  * no record is lost;
  * every new record is a residue added to the masks (ADDED_TSV, written by
    verify_hetatm_masks.py --added-tsv);
  * a record present in both changes only in neighbour-dependent fields
    (phi, psi, omega, is_cis, is_trans, peptide_bond, rama*), only between a
    value and null, and only when the residue is next to an added residue
    (adjacent in residue order within the chain, over the old records plus
    the added residues, so numbering jumps such as a fluorescent-protein
    chromophore standing for three residues, 2dd7 CR2 56 -> PHE 58, count
    as adjacent; for pruned_bfilt an added
    residue that fails the B-factor rule now ends a fragment, one that passes
    can join two fragments).  For lab-pruned datasets (top2018) the fragment
    records do not change, so no field should change at all.

Usage:
    compare_b1_filtered.py OLD.jsonl NEW.jsonl ADDED_TSV [--list 15]
"""

import argparse
import csv
import json
from collections import Counter

NEIGHBOUR_FIELDS = {"phi", "psi", "omega", "is_cis", "is_trans", "peptide_bond",
                    "rama_category", "rama5", "rama4", "rama3"}
KEYF = ("file", "model", "chain", "resnum", "ins")


def key(r):
    return (r["file"].split("_")[0], r["model"], r["chain"], r["resnum"], r["ins"] or " ")


def load(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            r = json.loads(line)
            if not r.get("_meta"):
                out[key(r)] = r
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("added_tsv")
    ap.add_argument("--list", type=int, default=15)
    a = ap.parse_args()

    added = {}
    for r in csv.DictReader(open(a.added_tsv), delimiter="\t"):
        added[(r["pdb"], r["chain"], int(r["resnum"]), r["icode"] or " ")] = (r["resname"], r["bfilt"])
    old, new = load(a.old), load(a.new)
    only_old = sorted(set(old) - set(new))
    only_new = sorted(set(new) - set(old))
    print(f"old {len(old)}  new {len(new)}  only old {len(only_old)}  only new {len(only_new)}")
    for k in only_old[: a.list]:
        print("   LOST", k, old[k]["resname"])

    not_added = [k for k in only_new if (k[0], k[2], k[3], k[4]) not in added]
    by_name = Counter(new[k]["resname"] for k in only_new)
    print(f"new records by resname: {dict(by_name.most_common(12))}")
    print(f"new records not in the added list: {len(not_added)}")
    for k in not_added[: a.list]:
        print("   ", k, new[k]["resname"])

    # Residue order per (pdb, chain): old records plus added residues.
    order = {}
    for k in set(old) | {(p, 1, c, n, i) for (p, c, n, i) in added}:
        order.setdefault((k[0], k[2]), set()).add((k[3], k[4]))
    neighbours = {}
    for (p, c), res in order.items():
        res = sorted(res)
        for i, r in enumerate(res):
            neighbours[(p, c) + r] = [x for x in (res[i - 1] if i else None,
                                                  res[i + 1] if i + 1 < len(res) else None) if x]

    def next_to_added(k):
        return any((k[0], k[2]) + nb in added for nb in neighbours.get((k[0], k[2], k[3], k[4]), []))

    change = Counter()
    unexplained = []
    for k in set(old) & set(new):
        o, n = old[k], new[k]
        for f in set(o) | set(n):
            if f in KEYF or o.get(f) == n.get(f):
                continue
            kind = ("unnulled" if o.get(f) is None else
                    "nulled" if n.get(f) is None else "value_changed")
            change[(f, kind)] += 1
            if f not in NEIGHBOUR_FIELDS or kind == "value_changed" or not next_to_added(k):
                unexplained.append((k, o["resname"], f, o.get(f), n.get(f)))
    print("field changes in common records:",
          dict(sorted(change.items())) if change else "none")
    rec_changed = len({k for k in set(old) & set(new) if old[k] != new[k]})
    print(f"common records changed: {rec_changed}")
    print(f"unexplained changes: {len(unexplained)}")
    for row in unexplained[: a.list]:
        print("   ", row)


if __name__ == "__main__":
    main()
