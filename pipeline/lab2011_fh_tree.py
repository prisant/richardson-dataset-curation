#!/usr/bin/env python3
"""Lay out the lab's 2011 Top8000 chain files as FH files of our pipeline.

Input: the extracted top8000_chains_70.tgz (Richardson lab, 7,957
single-chain files <pdb>FH_<chain>.pdb, Reduce 3.15.091106 with flips;
sha256 4a6bb1b9bb340cd42c14cd9054688732ac09317dfe0af49f3fd453f77ba486e6).
Each chain file carries the entry's full header (SEQRES of every chain,
CRYST1, ...) and the coordinates of its one chain.

Output under DST_DIR, in the layout of our top8000 tree:
    <pdb[:2]>/<pdb>/<pdb>FH.pdb   header records of the entry's first chain
                                  file (everything before its first
                                  coordinate record), then the coordinate
                                  records (ATOM, HETATM, ANISOU, SIGATM,
                                  SIGUIJ, TER) of each chain file in chain
                                  order, up to ENDMDL or END; then END
    Top8000-best_hom70_pdb_chain_list.csv   pdb_id,chain
The FH files then go through the steps of top8000/reproduce.sh from step 3
on (build_ersatz_only.py, build_bfactor_masks.py, pydangle, the filter),
so the lab's input files are measured and filtered exactly as ours.

Unlike rama8000_reconciliation/scripts/lab2011_tree.py, which keeps
coordinates only, this keeps the headers: build_bfactor_masks.py admits
HETATM polymer residues (MSE, ...) by the chain's SEQRES and TER records,
and DSSP needs CRYST1.

Usage:
    lab2011_fh_tree.py SRC_DIR DST_DIR
"""

import os
import sys
from collections import defaultdict

COORD = ("ATOM  ", "HETATM", "ANISOU", "SIGATM", "SIGUIJ", "TER   ")


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]
    files = defaultdict(list)
    for f in sorted(os.listdir(src)):
        if not f.endswith(".pdb") or "FH_" not in f:
            continue
        p, c = f[:-4].split("FH_")
        files[p.lower()].append((c, f))
    os.makedirs(dst, exist_ok=True)
    n_chains = 0
    with open(os.path.join(dst, "Top8000-best_hom70_pdb_chain_list.csv"), "w") as cl:
        cl.write("pdb_id,chain\n")
        for p, chains in sorted(files.items()):
            d = os.path.join(dst, p[:2], p)
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, f"{p}FH.pdb"), "w") as out:
                for i, (c, f) in enumerate(chains):
                    cl.write(f"{p},{c}\n")
                    n_chains += 1
                    in_coords = False
                    with open(os.path.join(src, f)) as fh:
                        for line in fh:
                            if line.startswith(("ENDMDL", "END   ")) or line.rstrip() == "END":
                                break
                            if line.startswith(COORD):
                                in_coords = True
                                out.write(line)
                            elif i == 0 and not in_coords and not line.startswith("MODEL "):
                                out.write(line)
                out.write("END\n")
    print(f"{n_chains} chain files, {len(files)} entries", file=sys.stderr)


if __name__ == "__main__":
    main()
