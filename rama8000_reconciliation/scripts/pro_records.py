"""Emit one record per Pro (with phi & psi defined) in the listed top8000 chains."""
import json, sys
from multiprocessing import Pool
import recount3 as R

def work(pdb):
    out = []
    res, order = R.parse(pdb)
    for ch, keys in order.items():
        for i, k in enumerate(keys):
            cur = res[k]
            if cur['rn'] != 'PRO' or not all(a in cur['atoms'] for a in ('N', 'CA', 'C')): continue
            prev = res[keys[i-1]] if i > 0 else None
            nxt = res[keys[i+1]] if i + 1 < len(keys) else None
            if not (R.linked(prev, cur) and R.linked(cur, nxt)) or 'CA' not in prev['atoms']: continue
            a = lambda r, n: r['atoms'][n]
            om = R.dih(a(prev,'CA'), a(prev,'C'), a(cur,'N'), a(cur,'CA'))
            phi = R.dih(a(prev,'C'), a(cur,'N'), a(cur,'CA'), a(cur,'C'))
            psi = R.dih(a(cur,'N'), a(cur,'CA'), a(cur,'C'), a(nxt,'N'))
            b = cur['b']
            out.append(dict(pdb=pdb, chain=ch, resnum=k[1], ins=k[2], phi=round(phi,2), psi=round(psi,2),
                omega=round(om,2), own_mc=cur['maxb_mc'], prevC=prev['b'].get('C',99), nextN=nxt['b'].get('N',99),
                tri_mc=max(prev['maxb_mc'], cur['maxb_mc'], nxt['maxb_mc']),
                alt=sorted(cur['alts']), mcalt=bool(cur.get('mcalt')), maxb_all=cur['maxb_all']))
    return out

if __name__ == '__main__':
    with Pool(8) as p, open(sys.argv[1], 'w') as f:
        for recs in p.imap_unordered(work, sorted(R.chains), chunksize=20):
            for r in recs: f.write(json.dumps(r) + '\n')
