"""Spread of the preceding residue's mainchain alternates, per residue, from the Reduce output (*FH.pdb), listed chains, model 1.

The lab's 2013 residue lists (Top8000_Rama_2013) drop a residue when the preceding C has alternates at different
positions (phi ambiguous), and keep it when the alternates coincide. Columns: largest distance (Å) between any two
alternates of the preceding residue's CA, C, O (0 when none). The preceding residue is the previous residue record
in file order within the chain.
TSV columns: pdb chain resnum ins prevCA_spread prevC_spread prevO_spread
usage: prevC_spread.py data/prevC_spread.tsv
"""
import sys, os
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from recount3 import chains, ROOT

def work(pdb):
    R, order = {}, []
    for line in open(f'{ROOT}/{pdb[:2]}/{pdb}/{pdb}FH.pdb'):
        if line.startswith('ENDMDL'): break
        if not line.startswith(('ATOM', 'HETATM')) or line[21] not in chains[pdb]: continue
        k = (line[21], int(line[22:26]), line[26])
        if k not in R: R[k] = {}; order.append(k)
        R[k].setdefault(line[12:16].strip(), []).append(np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    out = []
    for prev, k in zip(order, order[1:]):
        if prev[0] != k[0]: continue
        sp = [max((float(np.linalg.norm(a - b)) for a in R[prev].get(n, []) for b in R[prev].get(n, [])), default=0.0)
              for n in ('CA', 'C', 'O')]
        out.append(f'{pdb}\t{k[0]}\t{k[1]}\t{k[2]}\t' + '\t'.join(f'{s:.3f}' for s in sp))
    return out

if __name__ == '__main__':
    with Pool(8) as p, open(sys.argv[1], 'w') as f:
        f.write('pdb\tchain\tresnum\tins\tprevCA_spread\tprevC_spread\tprevO_spread\n')
        for recs in p.imap_unordered(work, sorted(chains), chunksize=20):
            for line in recs: f.write(line + '\n')
