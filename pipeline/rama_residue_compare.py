#!/usr/bin/env python3
"""Residue-by-residue comparison: our top8000 Ramachandran classes vs the lab's lists.

The lab's lists are residues/rama8000-<class>.csv of Top8000_Rama_2013
(gly, not gly-sym); each row is "<pdb> <chain> <resnum><icode> <resname>
...,phi,psi,omega,tau".  Ours are rama_category in a filtered pydangle
JSONL (model 1); the unfiltered JSONL of the same run tells whether a
residue we lack was measured at all.  Residues are matched on (entry,
chain, residue number, insertion code).

Every lab residue falls in one group:
  same class            in ours with the same class
  other class           in ours with another class (crosstab)
  ours: no phi / no psi / neither
                        in our filtered file without a category because the
                        angle is null; the lab's value is given in the TSV
  filtered out by us    measured by us (unfiltered file) but not in our
                        filtered file
  not measured: chain   chain absent from our unfiltered output
  not measured: residue chain present, residue absent
Every residue of ours with a class that the lab does not list falls in:
  lab: chain absent     the lab lists no residue of that chain
  lab: omega rule       fails the lab's omega test (|omega| > 150, CisPro
                        |omega| < 30, omega defined)
  lab: other            chain listed, omega passes
Writes the summary to stdout and every residue not in "same class" to
--tsv.

Usage:
    rama_residue_compare.py LAB_RESIDUES_DIR FILTERED.jsonl ALL.jsonl [--tsv OUT]
"""

import argparse
import json
import os
import re
from collections import Counter

CLASSES = (("General", "general-noGPIVpreP"), ("Gly", "gly"), ("IleVal", "ileval-nopreP"),
           ("PrePro", "prepro-noGP"), ("TransPro", "transpro"), ("CisPro", "cispro"))
NUM = re.compile(r"^(-?\d+)([A-Za-z]?)$")


def read_lab(directory):
    lab = {}
    for cls, fname in CLASSES:
        with open(os.path.join(directory, f"rama8000-{fname}.csv")) as fh:
            for row in fh:
                if not row.strip() or row.startswith("#"):
                    continue
                f = row.split(",")
                pdb, chain, num, resname = f[0].split()[:4]
                m = NUM.match(num)
                key = (pdb.lower(), chain, int(m.group(1)), m.group(2) or " ")
                if key in lab:
                    raise SystemExit(f"duplicate lab residue {key}")
                lab[key] = (cls, resname, f[1], f[2], f[3])
    return lab


def key_of(r):
    return (r["file"][:4].lower(), r["chain"], r["resnum"], r["ins"] or " ")


def omega_ok(cls, om):
    if om is None:
        return False
    return abs(om) < 30 if cls == "CisPro" else abs(om) > 150


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("lab")
    ap.add_argument("filtered")
    ap.add_argument("all")
    ap.add_argument("--tsv")
    a = ap.parse_args()
    lab = read_lab(a.lab)
    lab_chains = {k[:2] for k in lab}

    ours = {}
    with open(a.filtered) as fh:
        for line in fh:
            r = json.loads(line)
            if r.get("_meta") or r["model"] != 1:
                continue
            ours[key_of(r)] = r
    measured, chains = set(), set()
    with open(a.all) as fh:
        for line in fh:
            r = json.loads(line)
            if r.get("_meta") or r["model"] != 1:
                continue
            k = key_of(r)
            measured.add(k)
            chains.add(k[:2])

    groups = Counter()
    cross = Counter()
    by_class = Counter()
    rows = []
    for k, (cls, resname, phi, psi, om) in lab.items():
        r = ours.get(k)
        if r is not None and r.get("rama_category") == cls:
            g = "same class"
        elif r is not None and r.get("rama_category"):
            g = "other class"
            cross[(cls, r["rama_category"])] += 1
        elif r is not None:
            g = ("ours: neither" if r.get("phi") is None and r.get("psi") is None
                 else "ours: no phi" if r.get("phi") is None else "ours: no psi")
        elif k in measured:
            g = "filtered out by us"
        elif k[:2] in chains:
            g = "not measured: residue"
        else:
            g = "not measured: chain"
        groups[g] += 1
        by_class[(cls, g)] += 1
        if g != "same class":
            rows.append(("lab", g, *k[:3], k[3].strip(), resname, cls,
                          r.get("rama_category") if r else "", phi, psi, om,
                          "" if r is None else r.get("phi"), "" if r is None else r.get("psi"),
                          "" if r is None else r.get("omega")))

    extra = Counter()
    for k, r in ours.items():
        cls = r.get("rama_category")
        if not cls or k in lab:
            continue
        g = ("lab: chain absent" if k[:2] not in lab_chains
             else "lab: omega rule" if not omega_ok(cls, r.get("omega")) else "lab: other")
        extra[g] += 1
        by_class[(cls, g)] += 1
        rows.append(("ours", g, *k[:3], k[3].strip(), r["resname"], "", cls, "", "", "",
                     r.get("phi"), r.get("psi"), r.get("omega")))

    lab_groups = ("same class", "other class", "ours: no phi", "ours: no psi", "ours: neither",
                  "filtered out by us", "not measured: residue", "not measured: chain")
    our_groups = ("lab: chain absent", "lab: omega rule", "lab: other")
    print(f"lab residues {len(lab):,d}; our filtered residues with a class "
          f"{sum(1 for r in ours.values() if r.get('rama_category')):,d}")
    head = "".join(f"{c:>10s}" for c, _ in CLASSES)
    print(f"{'lab residues':24s}{'total':>10s}{head}")
    for g in lab_groups:
        print(f"  {g:22s}{groups[g]:10,d}" + "".join(f"{by_class[(c, g)]:10,d}" for c, _ in CLASSES))
    print(f"{'our residues not listed':24s}")
    for g in our_groups:
        print(f"  {g:22s}{extra[g]:10,d}" + "".join(f"{by_class[(c, g)]:10,d}" for c, _ in CLASSES))
    print("other class, lab -> ours:")
    for (c1, c2), n in cross.most_common():
        print(f"  {c1:10s} -> {c2:10s} {n:8,d}")
    if a.tsv:
        with open(a.tsv, "w") as out:
            out.write("side\tgroup\tpdb\tchain\tresnum\tins\tresname\tlab_class\tour_class"
                      "\tlab_phi\tlab_psi\tlab_omega\tour_phi\tour_psi\tour_omega\n")
            for row in sorted(rows, key=lambda x: (x[0], x[1], x[2], x[3], x[4])):
                out.write("\t".join("" if v is None else str(v) for v in row) + "\n")


if __name__ == "__main__":
    main()
