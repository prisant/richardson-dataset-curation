#!/usr/bin/env python3
"""Which altloc does Williams' Top2018 SQLite (cctbx/MolProbity) use for phi/psi?

For residues whose ``altlocs`` field is non-empty, recompute phi and psi from
each entry's whole-entry ``{pdbid}FH.pdb`` file (the lab's pruned chain files
omit residues with alternates, so they cannot be used) under the three atom-selection rules of
``altloc_java_diffs.py`` (maxocc, first, altA) and count which rule
reproduces the database value.  Only residues where the rules disagree with
each other are informative; they are counted separately.

Usage:
    altloc_vs_williams.py DB DATASET_DIR [--entries N] [--seed S]
        [--suffix FH] [--filter none|all|mc]

--filter all / mc restricts to residues passing all_pass_filters /
mc_pass_filters, and lists every residue where the rules differ.
"""

import argparse
import importlib.util
import random
import sqlite3
from collections import Counter
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "altloc_java_diffs", Path(__file__).with_name("altloc_java_diffs.py"))
_m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_m)
RULES = ("maxocc", "first", "altA")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("db")
    ap.add_argument("dataset_dir")
    ap.add_argument("--entries", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--suffix", default="FH")
    ap.add_argument("--filter", choices=("none", "all", "mc"), default="none")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    by_entry = {}
    for row in con.execute(
            "select r.pdbid, r.chain, r.resseq, r.icode, r.altlocs, a.phi, a.psi "
            "from residues r join ramachandran a using (chain_i, residue_i) "
            + {"none": "", "all": "join allatom_filters f using (chain_i, residue_i) ",
               "mc": "join mc_filters f using (chain_i, residue_i) "}[args.filter]
            + "where r.altlocs != ''"
            + {"none": "", "all": " and f.all_pass_filters", "mc": " and f.mc_pass_filters"}[args.filter]):
        by_entry.setdefault(row[0], []).append(row[1:])
    pdbids = sorted(by_entry)
    random.Random(args.seed).shuffle(pdbids)
    tally = Counter()
    informative = Counter()
    files_used = 0
    for pid in pdbids[: args.entries]:
        files = list(Path(args.dataset_dir).glob(f"*/{pid}/{pid}{args.suffix}.pdb"))
        if not files:
            tally["no_file"] += len(by_entry[pid])
            continue
        files_used += 1
        res, order = _m.read_pdb(files[0])
        index = {k: i for i, k in enumerate(order)}
        for chain, resseq, icode, _alts, phi, psi in by_entry[pid]:
            key = (chain, int(resseq), icode or " ")
            if key not in index:
                tally["residue_not_in_file"] += 1
                continue
            for field, dbv in (("phi", phi), ("psi", psi)):
                if dbv is None:
                    continue
                vals = {ru: _m.angle(res, order, key, field, ru, index) for ru in RULES}
                hit = "+".join(ru for ru in RULES if _m.close(vals[ru], float(dbv), 0.01)) or "none"
                tally[hit] += 1
                distinct = {round(v, 2) for v in vals.values() if v is not None}
                if len(distinct) > 1:
                    informative[hit] += 1
                    if args.filter != "none":
                        print(f"  {pid} {chain}{resseq.strip()}{(icode or '').strip()} {_alts:4s} {field}"
                              f" db {float(dbv):8.2f} "
                              + " ".join(f"{ru} {vals[ru]:8.2f}" for ru in RULES))
    print(f"entries sampled {min(args.entries, len(pdbids))} of {len(pdbids)} with altlocs; "
          f"files used {files_used}")
    print("all altloc residue angles, DB value reproduced by:")
    for k, n in tally.most_common():
        print(f"  {k:22s} {n}")
    print("angles where the rules give different values, DB value reproduced by:")
    for k, n in informative.most_common():
        print(f"  {k:22s} {n}")


if __name__ == "__main__":
    main()
