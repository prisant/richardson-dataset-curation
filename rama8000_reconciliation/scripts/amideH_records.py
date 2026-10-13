"""Amide-H attributes per residue from the Reduce output (*FH.pdb), listed chains, model 1.

TSV columns: pdb chain resnum ins has_H H_maxB H_minocc
"""
import sys, os
from collections import defaultdict
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from recount3 import chains, ROOT

def work(pdb):
    fn = f'{ROOT}/{pdb[:2]}/{pdb}/{pdb}FH.pdb'
    seen, H = {}, defaultdict(list)
    for line in open(fn):
        if line.startswith('ENDMDL'): break
        if not line.startswith(('ATOM', 'HETATM')): continue
        ch = line[21]
        if ch not in chains[pdb]: continue
        key = (ch, int(line[22:26]), line[26])
        seen[key] = 1
        if line[12:16].strip() == 'H':
            H[key].append((float(line[60:66]), float(line[54:60])))
    out = []
    for k in seen:
        h = H.get(k)
        out.append(f"{pdb}\t{k[0]}\t{k[1]}\t{k[2]}\t{int(bool(h))}\t{max(b for b, o in h) if h else 'nan'}\t{min(o for b, o in h) if h else 'nan'}")
    return out

if __name__ == '__main__':
    with Pool(8) as p, open(sys.argv[1], 'w') as f:
        f.write('pdb\tchain\tresnum\tins\thas_H\tH_maxB\tH_minocc\n')
        for recs in p.imap_unordered(work, sorted(chains), chunksize=20):
            for line in recs: f.write(line + '\n')
