"""Compact per-residue records (all residue types) for the listed top8000 chains, as TSV.

Columns: pdb chain resnum ins resname next_resname rama6 phi psi omega
         own_mc_max own_mc_min own_minocc own_alt d_prevC_N d_C_nextN
rama6 uses Read 2011 precedence (Gly, Pro > pre-Pro > Ile/Val > General); Pro cis if |omega| < 30.
Only residues with phi and psi computable (prev/next present with C/N atoms) are written.
"""
import sys, os
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from recount3 import chains, dih
from cis_records_ext import parse

MC = ('N', 'CA', 'C', 'O')

def rama6(rn, nxt_rn, om):
    if rn == 'GLY': return 'Gly'
    if rn == 'PRO': return 'CisPro' if om is not None and abs(om) < 30 else 'TransPro'
    if nxt_rn == 'PRO': return 'PrePro'
    if rn in ('ILE', 'VAL'): return 'IleVal'
    return 'General'

def work(pdb):
    out = []
    res, order = parse(pdb)
    for ch, keys in order.items():
        for i in range(1, len(keys) - 1):
            prev, cur, nxt = res[keys[i-1]], res[keys[i]], res[keys[i+1]]
            if not all(a in cur['xyz'] for a in ('N', 'CA', 'C')) or 'C' not in prev['xyz'] or 'N' not in nxt['xyz']: continue
            X = lambda r, a: r['xyz'][a]
            om = dih(X(prev,'CA'), X(prev,'C'), X(cur,'N'), X(cur,'CA')) if 'CA' in prev['xyz'] else None
            mcb = [b for a in MC if a in cur['b'] for b in cur['b'][a]]
            k = keys[i]
            out.append('\t'.join(map(str, (
                pdb, ch, k[1], k[2], cur['rn'], nxt['rn'], rama6(cur['rn'], nxt['rn'], om),
                round(dih(X(prev,'C'), X(cur,'N'), X(cur,'CA'), X(cur,'C')), 2),
                round(dih(X(cur,'N'), X(cur,'CA'), X(cur,'C'), X(nxt,'N')), 2),
                'nan' if om is None else round(om, 2), max(mcb), min(mcb),
                min(min(v) for v in cur['occ'].values()), int(any(cur['alt'].values())),
                round(float(np.linalg.norm(X(prev,'C') - X(cur,'N'))), 3),
                round(float(np.linalg.norm(X(cur,'C') - X(nxt,'N'))), 3)))))
    return out

if __name__ == '__main__':
    with Pool(8) as p, open(sys.argv[1], 'w') as f:
        f.write('pdb\tchain\tresnum\tins\tresname\tnext_resname\trama6\tphi\tpsi\tomega\town_mc_max\town_mc_min\town_minocc\town_alt\td_prevC_N\td_C_nextN\n')
        for recs in p.imap_unordered(work, sorted(chains), chunksize=20):
            for line in recs: f.write(line + '\n')
