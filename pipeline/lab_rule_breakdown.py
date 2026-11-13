#!/usr/bin/env python3
"""Which rule of the lab's Top8000 lists excludes each of our extra residues?

Input: the TSV of rama_residue_compare.py; the residues of ours that the
lab's lists do not include (side "ours"; by default group "lab: other",
i.e. chain listed by the lab and omega test passed).  Each residue is read
from its structure file (model 1, heavy atoms only: the lab's filters use
deposited heavy-atom data, not Reduce hydrogens) and tested against the
rules of get_srcdata.py (lines 110-125) and those found in the rama8000
reconciliation, in this order; the first rule it fails is reported:

  phi/psi null  our phi or psi is null (the lab requires phi, psi, omega and
                tau defined; pydangle assigns rama_category from residue
                identity even when phi or psi is null, e.g. at chain ends)
  altloc        any heavy atom of the residue has an altloc label
  occupancy<1   any heavy atom has occupancy < 1 (min_res_q = 1)
  B>=30         mainchain max B >= 30 (lab: max_mc_b < 30; ours keeps 30.00)
                (for MSE the max also includes SE, as in the lab's lists)
  B<=1          mainchain max B <= 1 (lab: max_mc_b > 1)
  prevC alt     not Pro, and the preceding residue's C has two or more
                positions more than 0.001 A apart
  none          passes all of these; with --identical (identical_inputs.tsv
                of the rama8000 reconciliation) split by whether the
                residue's filtering inputs are identical in the lab's 2011
                chain files and our current files

Usage:
    lab_rule_breakdown.py COMPARE.tsv DATA_DIR [--group "lab: other"]
        [--identical identical_inputs.tsv] [--tsv OUT]
"""

import argparse
import csv
import os
from collections import Counter, defaultdict

MC = ("N", "CA", "C", "O")
RULES = ("phi/psi null", "altloc", "occupancy<1", "B>=30", "B<=1", "prevC alt", "none")
NONE = ("none, inputs identical", "none, inputs differ", "none, not in 2011 files")


def heavy(line):
    el = line[76:78].strip()
    if el:
        return el != "H" and el != "D"
    return not line[12:16].strip().startswith(("H", "D"))


def read(path):
    """{(chain, resnum, ins): {'name', 'atoms': {atom: [(alt, occ, b, xyz)]}}}, order."""
    res, order = {}, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if line[:6] not in ("ATOM  ", "HETATM") or not heavy(line):
                continue
            key = (line[21], int(line[22:26]), line[26])
            if key not in res:
                res[key] = {"name": line[17:20], "atoms": defaultdict(list)}
                order.append(key)
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            res[key]["atoms"][line[12:16].strip()].append(
                (line[16], float(line[54:60]), float(line[60:66]), xyz))
    return res, order


def rule(res, order, key):
    r = res[key]
    recs = [a for v in r["atoms"].values() for a in v]
    if any(a[0] != " " for a in recs):
        return "altloc"
    if min(a[1] for a in recs) < 1:
        return "occupancy<1"
    names = MC + (("SE",) if r["name"] == "MSE" else ())
    bs = [a[2] for n in names for a in r["atoms"].get(n, [])]
    if bs and max(bs) >= 30:
        return "B>=30"
    if bs and max(bs) <= 1:
        return "B<=1"
    i = order.index(key)
    if r["name"] != "PRO" and i > 0 and order[i - 1][0] == key[0]:
        cs = [a[3] for a in res[order[i - 1]]["atoms"].get("C", [])]
        if any(sum((p - q) ** 2 for p, q in zip(a, b)) ** 0.5 > 0.001 for a in cs for b in cs):
            return "prevC alt"
    return "none"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("compare_tsv")
    ap.add_argument("data_dir")
    ap.add_argument("--group", default="lab: other")
    ap.add_argument("--tsv")
    ap.add_argument("--identical")
    a = ap.parse_args()
    ident = {}
    if a.identical:
        with open(a.identical) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                ident[(row["pdb"], row["chain"], int(row["resnum"]), row["ins"].strip())] = row["identical"] == "1"
    todo = defaultdict(list)
    with open(a.compare_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["side"] == "ours" and row["group"] == a.group:
                todo[row["pdb"]].append(row)
    tally, by_class = Counter(), Counter()
    out = open(a.tsv, "w") if a.tsv else None
    for pdb, rows in sorted(todo.items()):
        res, order = read(os.path.join(a.data_dir, pdb[:2], pdb, f"{pdb}_ersatz.pdb"))
        for row in rows:
            key = (row["chain"], int(row["resnum"]), row["ins"] or " ")
            if not row["our_phi"] or not row["our_psi"]:
                why = "phi/psi null"
            else:
                why = rule(res, order, key) if key in res else "not in file"
            if why == "none" and a.identical:
                same = ident.get((pdb, row["chain"], int(row["resnum"]), row["ins"].strip()))
                why = NONE[2] if same is None else NONE[0] if same else NONE[1]
            tally[why] += 1
            by_class[(row["our_class"], why)] += 1
            if out:
                out.write(f"{pdb}\t{row['chain']}\t{row['resnum']}\t{row['ins']}"
                          f"\t{row['resname']}\t{row['our_class']}\t{why}\n")
    classes = ("General", "Gly", "IleVal", "PrePro", "TransPro", "CisPro")
    print(f"{a.group}: {sum(tally.values()):,d} residues in {len(todo):,d} entries")
    print(f"{'first rule failed':18s}{'total':>9s}" + "".join(f"{c:>10s}" for c in classes))
    for why in RULES + NONE + ("not in file",):
        if tally[why] or why in RULES:
            print(f"{why:18s}{tally[why]:9,d}" + "".join(f"{by_class[(c, why)]:10,d}" for c in classes))


if __name__ == "__main__":
    main()
