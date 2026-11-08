#!/usr/bin/env python3
"""Check the Ramachandran-category precedence fix (register A2).

Compares pydangle JSONL before (OLD) and after (NEW) the fix, which moves
the Ile/Val test below the pre-Pro test (Read et al. 2011: Gly, Pro >
pre-Pro > Ile/Val > general).

1. Every field other than rama_category / rama5 / rama4 must be identical.
2. Category changes are tallied by (resname, old -> new) for rama_category,
   rama5 and rama4; the only expected change is Ile/Val IleVal -> PrePro
   (rama5 IleVal -> PrePro, rama4 General -> PrePro).
3. The rule is checked directly on NEW, using record order within a
   (file, model, chain) as fragment order: a residue that is not Gly or Pro
   and is followed by a Pro whose omega is non-null and whose residue number
   is the same or one higher (so the two are bonded, with no unreported
   modified residue between them, e.g. 1qfm A56 GLN - CSX 57 - PRO 58)
   must be PrePro; a PrePro residue must be followed by a Pro.  Residues
   whose rama_category is null (the post-filter nulls it when a neighbour
   is not in the filtered set) are skipped.  The same adjacency test also
   counts, in OLD, Ile/Val labelled IleVal before a bonded Pro -- the
   residues the fix should move.

Usage:
    check_rama_precedence.py OLD.jsonl NEW.jsonl [--list N]
"""

import argparse
import json
from collections import Counter

CATS = ("rama_category", "rama5", "rama4")


def load(path):
    recs = {}
    order = []
    with open(path) as fh:
        for line in fh:
            r = json.loads(line)
            if r.get("_meta"):
                continue
            k = (r["file"], r["model"], r["chain"], r["resnum"], r["ins"] or " ")
            recs[k] = r
            order.append(k)
    return recs, order


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--list", type=int, default=20)
    args = ap.parse_args()

    old, _ = load(args.old)
    new, order = load(args.new)
    print(f"old {len(old)}  new {len(new)}  only old {len(set(old) - set(new))}"
          f"  only new {len(set(new) - set(old))}")

    other = Counter()
    trans = {c: Counter() for c in CATS}
    for k in set(old) & set(new):
        o, n = old[k], new[k]
        for f in set(o) | set(n):
            if f in CATS or f in ("file", "model", "chain", "resnum", "ins"):
                continue
            if o.get(f) != n.get(f):
                other[f] += 1
        for c in CATS:
            if o.get(c) != n.get(c):
                trans[c][(o["resname"], o.get(c), n.get(c))] += 1
    print("non-category field differences:", dict(other) or "none")
    for c in CATS:
        print(f"{c} changes:")
        for (rn, a, b), cnt in sorted(trans[c].items(), key=lambda x: -x[1]):
            print(f"  {rn} {a} -> {b}: {cnt}")
        print(f"  total {sum(trans[c].values())}")

    _, old_order = load(args.old)
    old_ileval_before_pro = 0
    for i, k in enumerate(old_order):
        r = old[k]
        nxt = old[old_order[i + 1]] if i + 1 < len(old_order) else None
        if (r["resname"] in ("ILE", "VAL") and r.get("rama_category") == "IleVal"
                and nxt is not None and old_order[i + 1][:3] == k[:3]
                and nxt["resname"] == "PRO" and nxt.get("omega") is not None
                and nxt["resnum"] - r["resnum"] in (0, 1)):
            old_ileval_before_pro += 1
    print(f"OLD: Ile/Val labelled IleVal before a bonded Pro: {old_ileval_before_pro}")

    missing_prepro = []
    bad_prepro = []
    for i, k in enumerate(order):
        r = new[k]
        nxt = new[order[i + 1]] if i + 1 < len(order) else None
        same_chain = nxt is not None and order[i + 1][:3] == k[:3]
        next_is_bonded_pro = (same_chain and nxt["resname"] == "PRO"
                              and nxt.get("omega") is not None
                              and nxt["resnum"] - r["resnum"] in (0, 1))
        if (r["resname"] not in ("GLY", "PRO") and next_is_bonded_pro
                and r.get("rama_category") not in ("PrePro", None)):
            missing_prepro.append((k, r["resname"], r.get("rama_category")))
        if r.get("rama_category") == "PrePro" and not (same_chain and nxt["resname"] == "PRO"):
            bad_prepro.append((k, r["resname"]))
    print(f"rule check on NEW: non-Gly/Pro before bonded Pro but not PrePro: {len(missing_prepro)};"
          f" PrePro not followed by Pro: {len(bad_prepro)}")
    for row in (missing_prepro + bad_prepro)[: args.list]:
        print("  ", row)


if __name__ == "__main__":
    main()
