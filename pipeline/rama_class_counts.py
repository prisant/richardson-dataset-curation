#!/usr/bin/env python3
"""Residues per Ramachandran class: our JSONL vs the lab's Top8000 lists.

Our counts come from rama_category in a pydangle JSONL (a filtered dataset
file).  "Lab omega rule" applies the rule of the lab's get_srcdata.py
(lines 116-125, 150-154) to those residues: omega defined; |omega| > 150
for every class except CisPro, |omega| < 30 for CisPro (so twisted Pro,
pydangle TransPro, and non-Pro after a cis or twisted peptide are left
out).  The lab's counts are the rows of residues/rama8000-<class>.csv in
Top8000_Rama_2013 (gly, not gly-sym).

Usage:
    rama_class_counts.py --jsonl OURS.jsonl > ours.tsv
    rama_class_counts.py --lab RESIDUES_DIR --ours ours.tsv
"""

import argparse
import json
import os
from collections import Counter

CLASSES = (("General", "general-noGPIVpreP"), ("Gly", "gly"), ("IleVal", "ileval-nopreP"),
           ("PrePro", "prepro-noGP"), ("TransPro", "transpro"), ("CisPro", "cispro"))


def ours(path):
    all_, kept = Counter(), Counter()
    with open(path) as fh:
        for line in fh:
            r = json.loads(line)
            cat = r.get("rama_category")
            if r.get("_meta") or cat is None:
                continue
            all_[cat] += 1
            om = r.get("omega")
            if om is None:
                continue
            ok = abs(om) < 30 if cat == "CisPro" else abs(om) > 150
            kept[cat] += ok
    for cat, _ in CLASSES:
        print(f"{cat}\t{all_[cat]}\t{kept[cat]}")


def lab(directory, ours_tsv):
    mine = {}
    with open(ours_tsv) as fh:
        for line in fh:
            cat, a, k = line.split("\t")
            mine[cat] = (int(a), int(k))
    print(f"{'class':10s}{'lab lists':>12s}{'ours':>12s}{'ours-lab':>10s}"
          f"{'ours, lab omega rule':>22s}{'-lab':>8s}")
    tot = [0, 0, 0]
    for cat, fname in CLASSES:
        with open(os.path.join(directory, f"rama8000-{fname}.csv")) as fh:
            n = sum(1 for row in fh if row.strip() and not row.startswith("#"))
        a, k = mine[cat]
        tot = [tot[0] + n, tot[1] + a, tot[2] + k]
        print(f"{cat:10s}{n:12,d}{a:12,d}{a - n:+10,d}{k:22,d}{k - n:+8,d}")
    n, a, k = tot
    print(f"{'total':10s}{n:12,d}{a:12,d}{a - n:+10,d}{k:22,d}{k - n:+8,d}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jsonl")
    ap.add_argument("--lab")
    ap.add_argument("--ours")
    a = ap.parse_args()
    if a.jsonl:
        ours(a.jsonl)
    if a.lab:
        lab(a.lab, a.ours)


if __name__ == "__main__":
    main()
