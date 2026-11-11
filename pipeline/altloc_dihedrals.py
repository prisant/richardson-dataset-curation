#!/usr/bin/env python3
"""phi/psi of one residue for every combination of alternate conformers.

Diagnostic for values that disagree with the Richardson Lab's Top2018
database: prints phi = C(i-1),N,CA,C and psi = N,CA,C,N(i+1) for each
choice of altloc per atom of the window, so the database value can be
matched to the conformer(s) it was measured from.  Residues i-1, i, i+1 are
the residues before and after in the file's chain order (BioPython parse).

With --resname, residue i is the residue of that name at the position
(for split microheterogeneity, e.g. 1xpm D111 SCY/CYS), and residues
i-1 and i+1 are found by number, skipping other residues at position i.

Usage:
    altloc_dihedrals.py PDB CHAIN RESSEQ [ICODE] [--resname CYS]
"""

import argparse
import itertools

from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral
import math


def _alts(atom):
    if atom.is_disordered():
        return [(k, atom.disordered_get(k)) for k in sorted(atom.disordered_get_id_list())]
    return [(" ", atom)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("pdb")
    ap.add_argument("chain")
    ap.add_argument("resseq", type=int)
    ap.add_argument("icode", nargs="?", default=" ")
    ap.add_argument("--resname")
    a = ap.parse_args()
    chain = PDBParser(QUIET=True).get_structure("x", a.pdb)[0][a.chain]
    residues = [r for r in chain if r.has_id("CA")]
    here = (a.resseq, a.icode)
    idx = next(i for i, r in enumerate(residues) if r.id[1:] == here
               and (a.resname is None or r.get_resname() == a.resname))
    res = residues[idx]
    before = [r for r in residues[:idx] if r.id[1:] != here]
    after = [r for r in residues[idx + 1:] if r.id[1:] != here]
    prev = before[-1] if before else None
    nxt = after[0] if after else None
    for name, atoms in (("phi", [prev["C"], res["N"], res["CA"], res["C"]] if prev else None),
                        ("psi", [res["N"], res["CA"], res["C"], nxt["N"]] if nxt else None)):
        if atoms is None:
            continue
        for combo in itertools.product(*(_alts(x) for x in atoms)):
            labels = "".join(k for k, _ in combo)
            occ = " ".join(f"{x.get_occupancy():.2f}" for _, x in combo)
            ang = math.degrees(calc_dihedral(*(x.get_vector() for _, x in combo)))
            print(f"{name} {labels!r:8s} occ {occ}  {ang:9.3f}")


if __name__ == "__main__":
    main()
