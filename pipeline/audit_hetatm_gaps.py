#!/usr/bin/env python3
"""Audit false residue adjacencies and HETATM omissions around modified residues.

Two defects reported 2026-11-07 (Christopher Williams, 1gk8 C 257):

1. pydangle ``_build_fragments()`` (measure.py, since 9b25c27) builds
   fragments with ``CaPPBuilder(aa_only=0)`` and then drops residues not in
   the ``_PROTEIN_NAMES`` whitelist.  The residues on either side of a
   dropped *interior* residue become list-adjacent, so every i-1 / i+1 term
   (phi, omega of the residue after the gap; psi of the residue before it;
   neighbour-dependent labels) uses atoms from the wrong residue.
2. ``_build_mask()`` (build_ersatz_pdbs.py) reads ``ATOM`` records only,
   so HETATM residues the lab kept in its pruned files (e.g. MSE, SMC) are
   absent from our masks and hence from the filtered dataset.

For each dropped interior residue this script reports the nearest kept
residue on each side, whether each is in our mask (i.e. in the published
filtered dataset), and the C(before)-N(after) distance across the false
link.  A distance >= 2.0 A means the link is not a peptide bond and the
neighbour-dependent values are wrong.  Links < 2.0 A are counted separately
and are not yet explained.  It also lists non-water HETATM residues present
in the lab's pruned files but absent from our masks.

Requires the pydangle-biopython venv (imports its whitelist so the audit
uses exactly the code that produced the data).

Usage (from the Data directory):
    PY=~/Desktop/Projects/pydangle-biopython/.venv/bin/python
    $PY -I audit_hetatm_gaps.py top2018_pdbs_full_filtered_hom70 OUT.tsv [-j 8]

Outputs:
    OUT.tsv               one row per false link
    OUT_hetdropped.tsv    one row per lab-kept HETATM residue not in our mask
    summary on stdout
"""

import argparse
import csv
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import CaPPBuilder
from pydangle_biopython.measure import _is_polymer_residue

PEPTIDE_CN_MAX = 2.0  # Angstrom; above this C(i)-N(i+1) is not a bond


def _key(residue):
    return (residue.get_parent().id, residue.id[1], residue.id[2])


def _label(residue):
    return f"{residue.get_resname()}{residue.id[1]}{residue.id[2].strip()}"


def _load_mask(entry_dir):
    """Same parse as filter_pruned_residues.load_masks_and_fragments."""
    mask = set()
    for mask_file in entry_dir.glob("*.mask"):
        for line in mask_file.read_text().splitlines():
            parts = line.strip().split()
            if len(parts) >= 3:
                mask.add((parts[0], int(parts[1]), parts[2]))
            elif len(parts) == 2:
                mask.add((parts[0], int(parts[1]), " "))
    return mask


def scan_entry(entry_dir_str):
    entry_dir = Path(entry_dir_str)
    pdb_id = entry_dir.name
    ersatz = entry_dir / f"{pdb_id}_ersatz.pdb"
    mask = _load_mask(entry_dir)
    structure = PDBParser(QUIET=True).get_structure(pdb_id, str(ersatz))
    ppb = CaPPBuilder()
    links = []
    for chain in structure[0]:
        for pp in ppb.build_peptides(chain, aa_only=0):
            res = list(pp)
            i = 0
            while i < len(res):
                if _is_polymer_residue(res[i]):
                    i += 1
                    continue
                j = i
                while j < len(res) and not _is_polymer_residue(res[j]):
                    j += 1
                if i > 0 and j < len(res):  # interior run -> false link
                    a, b = res[i - 1], res[j]
                    dist = (a["C"] - b["N"]) if ("C" in a and "N" in b) else math.nan
                    links.append({
                        "pdb": pdb_id,
                        "chain": chain.id,
                        "dropped": "+".join(r.get_resname() for r in res[i:j]),
                        "before": _label(a),
                        "before_in_mask": _key(a) in mask,
                        "after": _label(b),
                        "after_in_mask": _key(b) in mask,
                        "C_N_dist": round(float(dist), 2),
                    })
                i = j
    het_dropped = []
    for pruned in entry_dir.glob(f"{pdb_id}_*_pruned_*.pdb"):
        seen = set()
        with open(pruned) as fh:
            for line in fh:
                if line.startswith("HETATM") and line[17:20] != "HOH":
                    k = (line[21], int(line[22:26]), line[26])
                    if k not in seen and k not in mask:
                        seen.add(k)
                        het_dropped.append({
                            "pdb": pdb_id, "chain": k[0], "resnum": k[1],
                            "icode": k[2], "resname": line[17:20],
                        })
    return links, het_dropped


def _write_tsv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def summarize(links, het_dropped, n_entries):
    def hit(rows):
        b = sum(r["before_in_mask"] for r in rows)
        a = sum(r["after_in_mask"] for r in rows)
        e = len({r["pdb"] for r in rows if r["before_in_mask"] or r["after_in_mask"]})
        return b, a, e

    broken = [r for r in links if not math.isnan(r["C_N_dist"]) and r["C_N_dist"] >= PEPTIDE_CN_MAX]
    short = [r for r in links if not math.isnan(r["C_N_dist"]) and r["C_N_dist"] < PEPTIDE_CN_MAX]
    undef = [r for r in links if math.isnan(r["C_N_dist"])]
    print(f"entries scanned: {n_entries}")
    print(f"false links (interior dropped runs): {len(links)}")
    for name, rows in (("broken, C-N >= 2.0 A", broken),
                       ("C-N < 2.0 A (unexplained)", short),
                       ("C or N missing", undef)):
        b, a, e = hit(rows)
        print(f"  {name}: {len(rows)} links; filtered residues before gap "
              f"(psi wrong) {b}, after gap (phi/omega wrong) {a}; entries {e}")
    mse = sum(h["resname"] == "MSE" for h in het_dropped)
    print(f"lab-kept non-water HETATM residues absent from our mask: "
          f"{len(het_dropped)} (MSE: {mse}; includes ligands/ions)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dataset_dir", type=Path)
    ap.add_argument("out_tsv")
    ap.add_argument("-j", "--jobs", type=int, default=8)
    args = ap.parse_args()

    entries = sorted(str(p.parent) for p in args.dataset_dir.glob("*/*/*_ersatz.pdb"))
    links, het_dropped = [], []
    with ProcessPoolExecutor(args.jobs) as ex:
        for lk, hd in ex.map(scan_entry, entries, chunksize=8):
            links += lk
            het_dropped += hd

    _write_tsv(args.out_tsv, links,
               ["pdb", "chain", "dropped", "before", "before_in_mask",
                "after", "after_in_mask", "C_N_dist"])
    _write_tsv(args.out_tsv.replace(".tsv", "_hetdropped.tsv"), het_dropped,
               ["pdb", "chain", "resnum", "icode", "resname"])
    summarize(links, het_dropped, len(entries))


if __name__ == "__main__":
    main()
