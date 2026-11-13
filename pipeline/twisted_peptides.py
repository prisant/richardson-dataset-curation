#!/usr/bin/env python3
"""Count cis and twisted peptides among residues with a Ramachandran category.

The lab's Top8000 residue lists (Top8000_Rama_2013, written by
get_srcdata.py) require a trans peptide, |omega| > 150, for every class
except CisPro, and |omega| < 30 for CisPro (get_srcdata.py lines 123-125,
150-154).  A twisted peptide, 30 <= |omega| <= 150, is in no class, and
neither is a cis peptide before a residue other than Pro.

For each JSONL file (pydangle output, filtered) this counts chains (file,
model, chain), residues (records), the residues that have a Ramachandran
category (rama_category not null) and, among those, the ones the lab's rule
leaves out (and the residues without a category, split by phi/psi null):
  twisted Pro       Pro, 30 <= |omega| <= 150 (pydangle: TransPro)
  twisted other     not Pro, 30 <= |omega| <= 150
  cis other         not Pro, |omega| < 30
Optionally LABEL=FILTERED.jsonl,ALL.jsonl also counts the residues (records)
in ALL, the unfiltered output (total residues).
With --lab DIR (the Top8000_Rama_2013 residues/ directory) it also checks
the lab's lists: rows per class whose omega lies outside the class's window
(expected 0).

Usage:
    twisted_peptides.py [--lab DIR] LABEL=FILTERED.jsonl[,ALL.jsonl] ...
"""

import argparse
import csv
import glob
import json
import os
from collections import Counter


def window(omega):
    a = abs(omega)
    return "cis" if a < 30 else "trans" if a > 150 else "twisted"


def count_jsonl(path):
    c = Counter()
    chains = set()
    with open(path) as fh:
        for line in fh:
            r = json.loads(line)
            if r.get("_meta"):
                continue
            chains.add((r["file"], r["model"], r["chain"]))
            c["residues"] += 1
            if r.get("rama_category") is None:
                phi, psi = r.get("phi") is None, r.get("psi") is None
                c["nc both null" if phi and psi else "nc phi null" if phi
                  else "nc psi null" if psi else "nc phi+psi"] += 1
                continue
            c["categorised"] += 1
            if r.get("omega") is None:
                c["omega null"] += 1
                continue
            w = window(r["omega"])
            pro = r["resname"] == "PRO"
            if w == "twisted":
                c["twisted Pro" if pro else "twisted other"] += 1
            elif w == "cis" and not pro:
                c["cis other"] += 1
    c["chains"] = len(chains)
    return c


def check_lab(directory):
    print("lab lists: rows whose omega lies outside the class's window")
    for path in sorted(glob.glob(os.path.join(directory, "rama8000-*.csv"))):
        cls = os.path.basename(path)[len("rama8000-"):-len(".csv")]
        want = "cis" if cls == "cispro" else "trans"
        n = bad = 0
        with open(path) as fh:
            for row in csv.reader(fh):
                if not row or row[0].startswith("#"):
                    continue
                omega = float(row[3])          # (-90, 270), get_srcdata.py line 88
                omega = omega - 360 if omega > 180 else omega
                n += 1
                bad += window(omega) != want
        print(f"  {cls:22s} rows {n:8d}  outside {bad}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--lab")
    ap.add_argument("files", nargs="*")
    a = ap.parse_args()
    if a.lab:
        check_lab(a.lab)
    if a.files:
        cols = ("total", "chains", "residues", "no category", "twisted Pro", "twisted other",
                "cis other", "omega null")
        print(f"{'dataset':14s}" + "".join(f"{c:>14s}" for c in cols))
        for spec in a.files:
            label, paths = spec.split("=", 1)
            path, _, total = paths.partition(",")
            c = count_jsonl(path)
            c["no category"] = c["residues"] - c["categorised"]
            if total:
                with open(total) as fh:
                    c["total"] = sum(1 for line in fh if '"_meta"' not in line)
            print(f"{label:14s}" + "".join(f"{c[k]:14,d}" for k in cols))
            print(f"{'':14s}  no category: " + ", ".join(
                f"{k[3:]} {c[k]:,d}" for k in ("nc phi null", "nc psi null", "nc both null",
                                                "nc phi+psi")))


if __name__ == "__main__":
    main()
