"""Mainchain atom records (N, CA, C, O) of the 7,957 Top8000-best_hom70 chains: the lab's 2011 chain files
(top8000_chains_70.tgz, <pdb>FH_<chain>.pdb) vs our current RCSB downloads (<pdb>.pdb, as deposited).

Per chain, first model, ATOM and HETATM records of residues that have both N and CA (so waters and ligands drop out).
Records are matched on (residue number, insertion code, atom name, altloc) and compared field by field as written
in the files: residue name, x/y/z, occupancy, B. Hydrogens play no part.
Usage (from rama8000_reconciliation/):
  python scripts/mc_atom_compare.py > results/mc_atom_compare.txt
"""
import csv, os
from collections import Counter
from multiprocessing import Pool
import numpy as np

OLD = os.path.expanduser('~/Desktop/Data/top8000_chains_70_2011/top8000_chains_70')
NEW = os.path.expanduser('~/Desktop/Data/top8000_pdbs_hom70')
MC = {'N', 'CA', 'C', 'O'}


def records(fn, chain):
    res, out = {}, {}
    for ln in open(fn):
        if ln.startswith('ENDMDL'): break
        if not ln.startswith(('ATOM', 'HETATM')) or ln[21] != chain: continue
        rk = (ln[22:26].strip(), ln[26].strip())
        name = ln[12:16].strip()
        res.setdefault(rk, set()).add(name)
        if name in MC:
            out[rk + (name, ln[16].strip())] = (ln[17:20].strip(), ln[30:38].strip(), ln[38:46].strip(), ln[46:54].strip(),
                                                ln[54:60].strip(), ln[60:66].strip())
    return {k: v for k, v in out.items() if {'N', 'CA'} <= res[k[:2]]}


def work(pc):
    pdb, ch = pc
    o = records(f'{OLD}/{pdb}FH_{ch}.pdb', ch)
    n = records(f'{NEW}/{pdb[:2]}/{pdb}/{pdb}.pdb', ch)
    c = Counter(old=len(o), new=len(n))
    dmax = 0.0
    ro, rn = {k[:2] for k in o}, {k[:2] for k in n}
    for k in o.keys() | n.keys():
        if k not in n: c['only 2011'] += 1; c['only 2011, residue in both'] += k[:2] in rn; continue
        if k not in o: c['only current'] += 1; c['only current, residue in both'] += k[:2] in ro; continue
        a, b = o[k], n[k]
        d = []
        if a[0] != b[0]: d.append('residue name')
        if a[1:4] != b[1:4]:
            d.append('coordinates')
            dist = float(np.linalg.norm(np.array(a[1:4], float) - np.array(b[1:4], float)))
            dmax = max(dmax, dist)
            c['coord > 0.01 Å'] += dist > 0.01; c['coord > 0.1 Å'] += dist > 0.1; c['coord > 1 Å'] += dist > 1
        if float(a[4]) != float(b[4]): d.append('occupancy')
        if float(a[5]) != float(b[5]): d.append('B')
        c['identical' if not d else 'differs'] += 1
        for x in d: c['field: ' + x] += 1
    return pdb, ch, c, dmax


if __name__ == '__main__':
    pcs = [(r['pdb_id'], r['chain']) for r in csv.DictReader(open(f'{NEW}/Top8000-best_hom70_pdb_chain_list.csv'))]
    with Pool(8) as p: R = p.map(work, pcs, chunksize=20)
    T = Counter()
    for _, _, c, _ in R: T.update(c)
    ent = Counter(); ent_any = set()
    for pdb, ch, c, _ in R:
        nd = c['differs'] + c['only 2011'] + c['only current']
        if nd: ent[pdb] += nd
    print(f'chains {len(R):,}; entries {len({p for p, *_ in R}):,}')
    print(f'mainchain atom records: 2011 files {T["old"]:,}; current files {T["new"]:,}')
    print(f'  matched (residue number, insertion code, atom, altloc): {T["identical"] + T["differs"]:,}')
    print(f'    identical in every field: {T["identical"]:,}')
    print(f'    differing in at least one field: {T["differs"]:,}')
    for f in ('residue name', 'coordinates', 'occupancy', 'B'):
        print(f'      {f}: {T["field: " + f]:,}')
    print(f'      coordinates moved > 0.01 Å: {T["coord > 0.01 Å"]:,}; > 0.1 Å: {T["coord > 0.1 Å"]:,}; > 1 Å: {T["coord > 1 Å"]:,};'
          f' max {max(d for *_, d in R):.2f} Å')
    print(f'  only in the 2011 files: {T["only 2011"]:,}')
    print(f'  only in the current files: {T["only current"]:,}')
    print(f'  unmatched records whose residue (number, insertion code) is in both files (atom or altloc label differs):'
          f' 2011 {T["only 2011, residue in both"]:,}; current {T["only current, residue in both"]:,}')
    tot = sum(ent.values())
    print(f'\nentries with any differing, missing or extra mainchain record: {len(ent)} of 7,737; records {tot:,}')
    for p, n in ent.most_common():
        if n >= 20: print(f'  {p} {n:,}')
    small = [n for n in ent.values() if n < 20]
    print(f'  {len(small)} entries with < 20 records each, {sum(small):,} records')
    fld = Counter()
    for pdb, ch, c, _ in R:
        if 0 < c['differs'] + c['only 2011'] + c['only current'] < 20:
            for f in ('residue name', 'coordinates', 'occupancy', 'B'): fld[f] += c['field: ' + f]
            fld['only 2011'] += c['only 2011']; fld['only current'] += c['only current']
    print(f'  in those small entries, by field: {dict(fld)}')
