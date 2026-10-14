"""Side-chain completeness and side-chain B/occupancy per residue, listed top8000 chains, model 1 (ersatz PDBs).

TSV columns: pdb chain resnum ins resname n_heavy n_missing missing sc_maxB sc_minocc
n_missing / missing: expected heavy atoms (standard residues + MSE) absent from the file; -1 / '' for other residue types.
sc_maxB / sc_minocc: over heavy atoms other than N, CA, C, O, OXT ('nan' if none).
"""
import sys, os
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from recount3 import chains
from cis_records_ext import parse

SC = dict(ALA='CB', ARG='CB CG CD NE CZ NH1 NH2', ASN='CB CG OD1 ND2', ASP='CB CG OD1 OD2', CYS='CB SG',
          GLN='CB CG CD OE1 NE2', GLU='CB CG CD OE1 OE2', GLY='', HIS='CB CG ND1 CD2 CE1 NE2', ILE='CB CG1 CG2 CD1',
          LEU='CB CG CD1 CD2', LYS='CB CG CD CE NZ', MET='CB CG SD CE', MSE='CB CG SE CE', PHE='CB CG CD1 CD2 CE1 CE2 CZ',
          PRO='CB CG CD', SER='CB OG', THR='CB OG1 CG2', TRP='CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
          TYR='CB CG CD1 CD2 CE1 CE2 CZ OH', VAL='CB CG1 CG2')
EXPECT = {k: ['N', 'CA', 'C', 'O'] + v.split() for k, v in SC.items()}
BB = {'N', 'CA', 'C', 'O', 'OXT'}

def work(pdb):
    out = []
    res, order = parse(pdb)
    for ch, keys in order.items():
        for k in keys:
            r = res[k]
            heavy = [a for a in r['xyz'] if not a.startswith(('H', 'D')) or r['rn'] == 'HG']
            exp = EXPECT.get(r['rn'])
            miss = [a for a in exp if a not in r['xyz']] if exp else []
            sc = [a for a in heavy if a not in BB]
            out.append('\t'.join(map(str, (pdb, ch, k[1], k[2], r['rn'], len(heavy), len(miss) if exp else -1, ','.join(miss),
                max((max(r['b'][a]) for a in sc), default='nan'), min((min(r['occ'][a]) for a in sc), default='nan')))))
    return out

if __name__ == '__main__':
    with Pool(8) as p, open(sys.argv[1], 'w') as f:
        f.write('pdb\tchain\tresnum\tins\tresname\tn_heavy\tn_missing\tmissing\tsc_maxB\tsc_minocc\n')
        for recs in p.imap_unordered(work, sorted(chains), chunksize=20):
            for line in recs: f.write(line + '\n')
