"""Extended per-residue attributes for every Pro with |omega| < OMEGA_MAX (default 45) in the listed top8000 chains.

Parses the ersatz PDBs directly (model 1). Keeps every altloc's B and occupancy so that
edge rules (partial occupancy, neighbour altlocs, near-threshold B, weak peptide links)
can be tested against the archived cispro.ndft grid.
"""
import json, sys, os, math
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from recount3 import chains, ROOT, dih

MC = ('N', 'CA', 'C', 'O')
OMEGA_MAX = 45.0  # override with argv[2]; 999 keeps every Pro

def parse(pdb):
    fn = f'{ROOT}/{pdb[:2]}/{pdb}/{pdb}_ersatz.pdb'
    res, order = {}, defaultdict(list)
    for line in open(fn):
        if line.startswith('ENDMDL'): break
        if not line.startswith(('ATOM', 'HETATM')): continue
        ch = line[21]
        if ch not in chains[pdb]: continue
        key = (ch, int(line[22:26]), line[26])
        name, alt, rn = line[12:16].strip(), line[16], line[17:20]
        r = res.get(key)
        if r is None:
            r = res[key] = dict(rn=rn, xyz={}, b=defaultdict(list), occ=defaultdict(list), alt=defaultdict(set))
            order[ch].append(key)
        if name not in r['xyz']:
            r['xyz'][name] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        r['b'][name].append(float(line[60:66])); r['occ'][name].append(float(line[54:60]))
        if alt != ' ': r['alt'][name].add(alt)
    return res, order

def work(pdb):
    out = []
    res, order = parse(pdb)
    for ch, keys in order.items():
        for i, k in enumerate(keys):
            cur = res[k]
            if cur['rn'] != 'PRO' or i == 0 or i + 1 == len(keys): continue
            prev, nxt = res[keys[i-1]], res[keys[i+1]]
            need = [(prev, 'CA'), (prev, 'C'), (cur, 'N'), (cur, 'CA'), (cur, 'C'), (nxt, 'N')]
            if not all(a in r['xyz'] for r, a in need): continue
            X = lambda r, a: r['xyz'][a]
            om = dih(X(prev,'CA'), X(prev,'C'), X(cur,'N'), X(cur,'CA'))
            if abs(om) >= OMEGA_MAX: continue
            mxb = lambda r, atoms: max((max(r['b'][a]) for a in atoms if a in r['b']), default=99.0)
            out.append(dict(
                pdb=pdb, chain=ch, resnum=k[1], ins=k[2], resname_prev=prev['rn'], resname_next=nxt['rn'],
                phi=round(dih(X(prev,'C'), X(cur,'N'), X(cur,'CA'), X(cur,'C')), 2),
                psi=round(dih(X(cur,'N'), X(cur,'CA'), X(cur,'C'), X(nxt,'N')), 2),
                omega=round(om, 2),
                d_prevC_N=round(float(np.linalg.norm(X(prev,'C') - X(cur,'N'))), 3),
                d_C_nextN=round(float(np.linalg.norm(X(cur,'C') - X(nxt,'N'))), 3),
                own_mc=mxb(cur, MC), own_ncac=mxb(cur, ('N','CA','C')), own_O=mxb(cur, ('O',)),
                own_all=mxb(cur, list(cur['b'])), own_cb=mxb(cur, ('CB',)),
                prevC=mxb(prev, ('C',)), nextN=mxb(nxt, ('N',)),
                prev_mc=mxb(prev, MC), next_mc=mxb(nxt, MC),
                has_O='O' in cur['xyz'],
                own_alt=sorted({x for s in cur['alt'].values() for x in s}),
                own_mcalt=any(cur['alt'].get(a) for a in MC),
                prev_alt=any(prev['alt'].values()), next_alt=any(nxt['alt'].values()),
                prev_mcalt=any(prev['alt'].get(a) for a in MC), next_mcalt=any(nxt['alt'].get(a) for a in MC),
                own_minocc=min(min(v) for v in cur['occ'].values()),
                nbr_minocc=min(min(min(v) for v in prev['occ'].values()), min(min(v) for v in nxt['occ'].values())),
                n_atoms=len(cur['xyz'])))
    return out

if __name__ == '__main__':
    if len(sys.argv) > 2: OMEGA_MAX = float(sys.argv[2])
    with Pool(8) as p, open(sys.argv[1], 'w') as f:
        for recs in p.imap_unordered(work, sorted(chains), chunksize=20):
            for r in recs: f.write(json.dumps(r) + '\n')
