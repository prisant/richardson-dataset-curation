#!/usr/bin/env python3
"""Measure polymer HETATM residues missing from residue masks (register B1).

The mask builders (`build_ersatz_pdbs._build_mask` for top2018,
`build_bfactor_masks._parse_residues` for top100/500/8000) read ATOM lines
only, so a modified residue written as HETATM (MSE, SMC, KCX, ...) never
enters a mask even when the residue is part of the kept polymer.

A HETATM residue counts as polymer here when its residue name appears in the
SEQRES records of its chain (the same self-evidencing rule the DSSP cleanup
uses) and it is not water.  Residues that satisfy that rule but sit after the
chain's first TER are counted separately ("after TER") as a check that the
rule is not picking up ligands.

Source file per mask:
  pruned_all / pruned_mc (top2018): the lab's pruned file; deleted residues
      are physically absent, so presence means "kept by the lab".
  pruned_bfilt (top100/500/8000): the ersatz file, chain from the mask name;
      the residue is classified by the mask builder's own B-factor rule
      (>= 3 of N, CA, C, O present, all with B <= BMAX).

Usage:
    measure_hetatm_mask_gap.py DATASET_DIR SUFFIX [--bmax 30] [--top 15]
"""

import argparse
from collections import Counter
from pathlib import Path

from pydangle_biopython.measure import _PROTEIN_NAMES

MC = {" N  ", " CA ", " C  ", " O  "}
WATER = {"HOH", "DOD", "WAT"}


def seqres_names(lines):
    names = {}
    for ln in lines:
        if ln.startswith("SEQRES"):
            names.setdefault(ln[11], set()).update(ln[19:].split())
    return names


def scan(pdb_path, chain, mask, bmax, bfilt):
    lines = open(pdb_path).read().splitlines()
    seq = seqres_names(lines).get(chain, set())
    after_ter = False
    res = {}  # key -> [resname, after_ter, mc_b list]
    for ln in lines:
        if ln.startswith("TER") and len(ln) > 21 and ln[21] == chain:
            after_ter = True
            continue
        if not ln.startswith("HETATM") or ln[21] != chain:
            continue
        rn = ln[17:20].strip()
        if rn in WATER:
            continue
        key = (chain, int(ln[22:26]), ln[26])
        r = res.setdefault(key, [rn, after_ter, [], rn in seq])
        if ln[12:16] in MC:
            try:
                r[2].append(float(ln[60:66]))
            except ValueError:
                pass
    out = []
    for key, (rn, at, bs, in_seq) in res.items():
        if key in mask:
            continue
        passes = len(bs) >= 3 and all(b <= bmax for b in bs)
        out.append((key, rn, at, passes, in_seq))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dataset")
    ap.add_argument("suffix")
    ap.add_argument("--bmax", type=float, default=30.0)
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--dump", help="write every missing polymer HETATM residue to this TSV")
    a = ap.parse_args()
    bfilt = a.suffix == "pruned_bfilt"

    by_name = Counter()
    n_masks = n_files_hit = 0
    after_ter = Counter()
    bf = Counter()
    not_seq = Counter()
    examples = []
    dump = []
    for mask_file in sorted(Path(a.dataset).rglob("*.mask")):
        pdbid, chain = mask_file.stem.rsplit("_", 1)
        mask = set()
        for ln in open(mask_file):
            # "C 254 A" / "C 254  ": chain is column 1 (may be blank)
            toks = ln[1:].split()
            if toks and toks[0].lstrip("-").isdigit():
                mask.add((ln[0], int(toks[0]), toks[1] if len(toks) > 1 else " "))
        n_masks += 1
        src = (mask_file.parent / f"{pdbid}_ersatz.pdb" if bfilt
               else mask_file.parent / f"{pdbid}_{chain}_{a.suffix}.pdb")
        if not src.exists():
            continue
        miss = scan(src, chain, mask, a.bmax, bfilt)
        if any(in_seq and not at for _, _, at, _, in_seq in miss):
            n_files_hit += 1
        for key, rn, at, passes, in_seq in miss:
            if not in_seq:
                # Rejected by the SEQRES test.  Count the ones before TER
                # (non-water) to see what the rule loses.
                if not at:
                    not_seq[rn] += 1
                continue
            if at:
                after_ter[rn] += 1
                continue
            if bfilt and not passes:
                bf["fail"] += 1
                continue
            by_name[rn] += 1
            dump.append((pdbid, *key, rn))
            if len(examples) < 5:
                examples.append((pdbid, key, rn))
    total = sum(by_name.values())
    measured = sum(c for n, c in by_name.items() if n in _PROTEIN_NAMES)
    print(f"{a.dataset}  suffix={a.suffix}  masks={n_masks}  chains affected={n_files_hit}")
    print(f"  missing polymer HETATM residues{' (pass B rule)' if bfilt else ''}: {total}"
          f"  (pydangle measures {measured})")
    if bfilt:
        print(f"  missing but fail B rule: {bf['fail']}")
    print(f"  after TER (excluded, ligand check): {sum(after_ter.values())} {dict(after_ter.most_common(8))}"
          f"  MSE after TER: {after_ter['MSE']}")
    print(f"  MSE missing, all positions: {by_name['MSE'] + after_ter['MSE']}")
    print("  by resname:", ", ".join(f"{n}{'*' if n in _PROTEIN_NAMES else ''} {c}"
                                     for n, c in by_name.most_common(a.top)))
    print(f"  not in chain SEQRES, before TER (excluded): {sum(not_seq.values())} "
          f"{dict(not_seq.most_common(10))}")
    if a.dump:
        with open(a.dump, "w") as fh:
            fh.write("pdb\tchain\tresnum\ticode\tresname\n")
            for row in dump:
                fh.write("\t".join(map(str, row)) + "\n")
    print("  (* = in pydangle _PROTEIN_NAMES)  examples:", examples)


if __name__ == "__main__":
    main()
