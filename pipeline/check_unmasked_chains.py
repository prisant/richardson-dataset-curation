#!/usr/bin/env python3
"""List chain-list entries that have no mask, and say why (register B7).

For top100/500/8000 (build_bfactor_masks.py).  For every chain in the
chain list with no {pdbid}_{chain}.mask, reports:

  ersatz_res   polymer residues of that chain in {pdbid}_ersatz.pdb
               (polymer_residues.mask_lines, the rule the mask builder uses)
  best_mc_b    lowest, over residues with at least N, CA, C, of the
               residue's highest mainchain B (N, CA, C, O); the chain is
               genuinely fully B-filtered when this exceeds --bmax
  orig_chain   ATOM lines with that chain ID in the original {pdbid}.pdb
  orig_blank   ATOM lines with a blank chain ID in the original {pdbid}.pdb

Usage:
    check_unmasked_chains.py DATASET_DIR CHAIN_LIST [--bmax 30]
"""

import argparse
import csv
import importlib.util
from collections import defaultdict
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "polymer_residues", Path(__file__).with_name("polymer_residues.py"))
_polymer = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_polymer)

_MC = {" N  ", " CA ", " C  ", " O  "}


def _lines(path):
    return path.read_text().splitlines() if path.exists() else []


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dataset_dir", type=Path)
    ap.add_argument("chain_list", type=Path)
    ap.add_argument("--bmax", type=float, default=30.0)
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.chain_list)))
    missing = []
    for r in rows:
        pdb, chain = r["pdb_id"].lower().strip(), r["chain"].strip()
        d = a.dataset_dir / pdb[:2] / pdb
        if not (d / f"{pdb}_{chain}.mask").exists():
            missing.append((pdb, chain, d))
    print(f"{a.dataset_dir.name}: {len(rows)} chains in list, "
          f"{len(missing)} without a mask (bmax {a.bmax})")
    print("pdb   chain  ersatz_res  best_mc_b  orig_chain  orig_blank  verdict")
    counts = defaultdict(int)
    for pdb, chain, d in missing:
        ers = _lines(d / f"{pdb}_ersatz.pdb")
        mc = defaultdict(list)
        res = set()
        for ln in _polymer.mask_lines(ers, [chain]):
            key = (int(ln[22:26]), ln[26])
            res.add(key)
            if ln[12:16] in _MC:
                mc[key].append(float(ln[60:66]))
        best = min((max(v) for v in mc.values() if len(v) >= 3), default=None)
        orig = [ln for ln in _lines(d / f"{pdb}.pdb") if ln.startswith("ATOM  ")]
        n_chain = sum(ln[21] == chain for ln in orig)
        n_blank = sum(ln[21] == " " for ln in orig)
        if not res:
            verdict = "chain absent from ersatz"
        elif best is None:
            verdict = "no residue with N, CA, C (not protein)"
        elif best <= a.bmax:
            verdict = "HAS PASSING RESIDUES"
        else:
            verdict = "fully B-filtered"
        if n_chain == 0 and n_blank:
            verdict += "; original chain ID blank"
        counts[verdict] += 1
        print(f"{pdb}  {chain:5}  {len(res):10}  "
              f"{best if best is not None else '-':>9}  {n_chain:10}  {n_blank:10}  {verdict}")
    for v, n in sorted(counts.items()):
        print(f"  {n:5}  {v}")


if __name__ == "__main__":
    main()
