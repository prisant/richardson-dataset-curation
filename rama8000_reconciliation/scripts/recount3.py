"""Recount rama6 classes in top8000 under alternative residue filters.
Parses ersatz PDBs directly (first altloc / highest occupancy per atom).
"""
import csv, math, sys, os
from collections import Counter, defaultdict
from multiprocessing import Pool
import numpy as np

ROOT = os.path.expanduser('~/Desktop/Data/top8000_pdbs_hom70')
chains = defaultdict(set)
for r in csv.DictReader(open(f'{ROOT}/Top8000-best_hom70_pdb_chain_list.csv')):
    k = list(r.keys())
    chains[r[k[0]].lower()].add(r[k[1]])

MC = ('N', 'CA', 'C', 'O')
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())

def dih(p0, p1, p2, p3):
    b0 = p0 - p1; b1 = p2 - p1; b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - b0.dot(b1) * b1; w = b2 - b2.dot(b1) * b1
    return math.degrees(math.atan2(np.cross(b1, v).dot(w), v.dot(w)))

def parse(pdb):
    fn = f'{ROOT}/{pdb[:2]}/{pdb}/{pdb}_ersatz.pdb'
    if not os.path.exists(fn):
        fn = f'{ROOT}/{pdb[1:3]}/{pdb}/{pdb}_ersatz.pdb'
    res = {}  # (chain,resseq,icode) -> dict
    order = defaultdict(list)
    for line in open(fn):
        if line.startswith('ENDMDL'): break
        if not line.startswith(('ATOM', 'HETATM')): continue
        ch = line[21]
        if ch not in chains[pdb]: continue
        key = (ch, int(line[22:26]), line[26])
        name = line[12:16].strip(); alt = line[16]; rn = line[17:20]
        b = float(line[60:66]); occ = float(line[54:60])
        r = res.get(key)
        if r is None:
            r = res[key] = dict(rn=rn, atoms={}, maxb_all=0.0, maxb_mc=0.0, alts=set(), b={})
            order[ch].append(key)
        if alt != ' ':
            r['alts'].add(alt)
            if name in MC: r['mcalt'] = True
        # keep first altloc encountered (A) for coordinates, but take max B over all alts
        if name not in r['atoms']:
            r['atoms'][name] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        r['maxb_all'] = max(r['maxb_all'], b)
        if name in MC: r['maxb_mc'] = max(r['maxb_mc'], b)
        r['b'][name] = max(r['b'].get(name, 0.0), b)
    return res, order

def linked(a, b):
    return a and b and 'C' in a['atoms'] and 'N' in b['atoms'] and \
        np.linalg.norm(a['atoms']['C'] - b['atoms']['N']) < 2.0

def classify(prev, cur, nxt):
    rn = cur['rn']
    if rn == 'GLY': return 'Gly'
    if rn == 'PRO':
        om = None
        if linked(prev, cur) and 'CA' in prev['atoms']:
            om = dih(prev['atoms']['CA'], prev['atoms']['C'], cur['atoms']['N'], cur['atoms']['CA'])
        return 'CisPro' if om is not None and abs(om) < 30 else 'TransPro'
    if nxt and nxt['rn'] == 'PRO': return 'PrePro'
    if rn in ('ILE', 'VAL'): return 'IleVal'
    return 'General'

def work(pdb):
    out = Counter()
    try:
        res, order = parse(pdb)
    except Exception as e:
        return Counter({('ERR',): 1})
    for ch, keys in order.items():
        for i, k in enumerate(keys):
            cur = res[k]
            if not all(a in cur['atoms'] for a in ('N', 'CA', 'C')): continue
            prev = res[keys[i-1]] if i > 0 else None
            nxt = res[keys[i+1]] if i + 1 < len(keys) else None
            if not (linked(prev, cur) and linked(cur, nxt)): continue  # phi & psi defined
            std = cur['rn'] in STD
            if not std and cur['rn'] != 'MSE': continue
            cls = classify(prev, cur, nxt)
            if cur['rn'] == 'MSE': out[('MSE_present',)] += 1
            own_mc = cur['maxb_mc']
            tri_mc = max(prev['maxb_mc'], own_mc, nxt['maxb_mc'])
            alt = bool(cur['alts']); mcalt = cur.get('mcalt', False)
            nbalt = prev.get('mcalt', False) or nxt.get('mcalt', False)
            dihB = max(prev['b'].get('C',99), cur['b']['N'], cur['b']['CA'], cur['b']['C'], nxt['b'].get('N',99))
            dihBO = max(dihB, cur['b'].get('O',99))
            for tag, ok in (('dih<30', dihB < 30), ('dih+O<30', dihBO < 30), ('dih+O<=30', dihBO <= 30),
                            ('own_mc+prevC+nextN<30', max(own_mc, prev['b'].get('C',99), nxt['b'].get('N',99)) < 30),('own_mc<30', own_mc < 30), ('own_mc<=30', own_mc <= 30),
                            ('tri_mc<=30', tri_mc <= 30), ('tri_mc<30', tri_mc < 30),
                            ('own_all<30', cur['maxb_all'] < 30)):
                if ok:
                    out[(tag, 'std+MSE', cls)] += 1
                    if std: out[(tag, 'std', cls)] += 1
                    if std and not alt: out[(tag, 'std_noalt', cls)] += 1
                    mse = std or cur['rn'] == 'MSE'
                    if mse and not alt: out[(tag, 'std+MSE_noalt', cls)] += 1
                    if std and not mcalt: out[(tag, 'std_nomcalt', cls)] += 1
                    if mse and not mcalt: out[(tag, 'std+MSE_nomcalt', cls)] += 1
                    if mse and not (mcalt or nbalt): out[(tag, 'std+MSE_nomcalt3', cls)] += 1
    return out

if __name__ == '__main__':
    tot = Counter()
    with Pool(8) as p:
        for c in p.imap_unordered(work, sorted(chains), chunksize=20):
            tot.update(c)
    import json
    json.dump({'|'.join(k): v for k, v in tot.items()}, open(sys.argv[1], 'w'), indent=0)
    print('done', tot[('ERR',)], 'errors;', tot[('MSE_present',)], 'MSE residues w/ phi,psi')
