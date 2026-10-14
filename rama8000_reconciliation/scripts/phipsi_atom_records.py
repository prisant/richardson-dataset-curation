"""Occupancy and altloc of the neighbour mainchain atoms that define phi, psi and omega, per residue (ersatz PDBs, model 1).

TSV columns: pdb chain resnum ins, then for each of prevCA prevC prevO prevN nextN nextCA:
  <atom>_occ (min occupancy over altlocs; nan if absent) and <atom>_alt (1 if the atom has an altloc label).
phi = prevC-N-CA-C, psi = N-CA-C-nextN, omega = prevCA-prevC-N-CA. Same residues as residue_records.py.
"""
import sys, os
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from recount3 import chains
from cis_records_ext import parse

ATOMS = (('prev', 'CA'), ('prev', 'C'), ('prev', 'O'), ('prev', 'N'), ('next', 'N'), ('next', 'CA'))

def work(pdb):
    out = []
    res, order = parse(pdb)
    for ch, keys in order.items():
        for i in range(1, len(keys) - 1):
            prev, cur, nxt = res[keys[i-1]], res[keys[i]], res[keys[i+1]]
            if not all(a in cur['xyz'] for a in ('N', 'CA', 'C')) or 'C' not in prev['xyz'] or 'N' not in nxt['xyz']: continue
            r = dict(prev=prev, next=nxt); k = keys[i]; f = [pdb, ch, k[1], k[2]]
            for w, a in ATOMS:
                f += [min(r[w]['occ'][a]) if a in r[w]['occ'] else 'nan', int(bool(r[w]['alt'].get(a)))]
            out.append('\t'.join(map(str, f)))
    return out

if __name__ == '__main__':
    with Pool(8) as p, open(sys.argv[1], 'w') as f:
        f.write('\t'.join(['pdb', 'chain', 'resnum', 'ins'] + [f'{w}{a}_{x}' for w, a in ATOMS for x in ('occ', 'alt')]) + '\n')
        for recs in p.imap_unordered(work, sorted(chains), chunksize=20):
            for line in recs: f.write(line + '\n')
