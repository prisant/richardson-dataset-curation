"""All-atom occupancy and side-chain B per residue from the Reduce output (*FH.pdb), listed chains, model 1.

R-lab (2026-10-14): B cutoff 30 on mainchain (N, CA, C, O), 40 on side-chain heavy (and other) atoms; the occupancy
filter may have run on the hydrogenated files. Columns:
  all_minocc    min occupancy over every atom, hydrogens included
  Hsc_minocc    min occupancy over hydrogens other than the amide H
  sc_maxB_heavy max B over heavy atoms other than N, CA, C, O, OXT
  sc_maxB_all   the same over all atoms, hydrogens included (amide H and HA count as side chain: not N, CA, C, O, OXT)
TSV columns: pdb chain resnum ins all_minocc Hsc_minocc sc_maxB_heavy sc_maxB_all
"""
import sys, os
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from recount3 import chains, ROOT

BB = {'N', 'CA', 'C', 'O', 'OXT'}
NAN = float('nan')

def is_h(line):
    el = line[76:78].strip()
    if el: return el in ('H', 'D')
    return line[12:16].strip().lstrip('0123456789')[:1] in ('H', 'D')

def work(pdb):
    fn = f'{ROOT}/{pdb[:2]}/{pdb}/{pdb}FH.pdb'
    R = {}
    for line in open(fn):
        if line.startswith('ENDMDL'): break
        if not line.startswith(('ATOM', 'HETATM')): continue
        ch = line[21]
        if ch not in chains[pdb]: continue
        key = (ch, int(line[22:26]), line[26])
        r = R.setdefault(key, [NAN] * 4)
        name, occ, b, h = line[12:16].strip(), float(line[54:60]), float(line[60:66]), is_h(line)
        r[0] = occ if r[0] != r[0] else min(r[0], occ)
        if h and name != 'H': r[1] = occ if r[1] != r[1] else min(r[1], occ)
        if name not in BB:
            if not h: r[2] = b if r[2] != r[2] else max(r[2], b)
            r[3] = b if r[3] != r[3] else max(r[3], b)
    return [f'{pdb}\t{k[0]}\t{k[1]}\t{k[2]}\t' + '\t'.join(str(x) for x in v) for k, v in R.items()]

if __name__ == '__main__':
    with Pool(8) as p, open(sys.argv[1], 'w') as f:
        f.write('pdb\tchain\tresnum\tins\tall_minocc\tHsc_minocc\tsc_maxB_heavy\tsc_maxB_all\n')
        for recs in p.imap_unordered(work, sorted(chains), chunksize=20):
            for line in recs: f.write(line + '\n')
