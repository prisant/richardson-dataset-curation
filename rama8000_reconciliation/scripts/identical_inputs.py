"""Residues whose filtering inputs are identical in the lab's 2011 chain files and our current RCSB downloads.

Inputs of a residue = every heavy-atom record of the residue itself, plus the mainchain records (N, CA, C, O) and the
residue name of the residues immediately before and after it in the chain (file order, first model). Records are
compared as written: atom name, altloc, residue name, x/y/z, occupancy, B -- except that side-chain x/y/z are not
compared: no filtering rule uses them, and the 2011 files carry Reduce's Asn/Gln/His flips, which move only side-chain
atoms. Hydrogens (element H or D) are ignored.
2011: <pdb>FH_<chain>.pdb from top8000_chains_70.tgz; current: <pdb>.pdb as deposited.
Writes data/identical_inputs.tsv (pdb chain resnum ins identical) for every residue of the 7,957 chains that has N and CA
in either file set.
Usage (from rama8000_reconciliation/):  python scripts/identical_inputs.py
"""
import csv, os
from multiprocessing import Pool

OLD = os.path.expanduser('~/Desktop/Data/top8000_chains_70_2011/top8000_chains_70')
NEW = os.path.expanduser('~/Desktop/Data/top8000_pdbs_hom70')
MC = ('N', 'CA', 'C', 'O')


def residues(fn, chain):
    order, rec = [], {}
    for ln in open(fn):
        if ln.startswith('ENDMDL'): break
        if not ln.startswith(('ATOM', 'HETATM')) or ln[21] != chain: continue
        name = ln[12:16].strip()
        el = ln[76:78].strip() or name[0]
        if el in ('H', 'D'): continue
        k = (ln[22:26].strip(), ln[26].strip())
        if k not in rec: rec[k] = []; order.append(k)
        xyz = (ln[30:38].strip(), ln[38:46].strip(), ln[46:54].strip()) if name in MC else ()
        rec[k].append((name, ln[16].strip(), ln[17:20].strip()) + xyz + (float(ln[54:60]), float(ln[60:66])))
    return order, {k: tuple(sorted(v)) for k, v in rec.items()}


def signatures(order, rec):
    mc = lambda k: tuple(r for r in rec[k] if r[0] in MC)
    out = {}
    for i, k in enumerate(order):
        if not {'N', 'CA'} <= {r[0] for r in rec[k]}: continue
        p = order[i - 1] if i > 0 else None; n = order[i + 1] if i + 1 < len(order) else None
        out[k] = (rec[k], p, p and mc(p), n, n and mc(n), n and tuple(sorted({r[2] for r in rec[n]})))
    return out


def work(pc):
    pdb, ch = pc
    a = signatures(*residues(f'{OLD}/{pdb}FH_{ch}.pdb', ch))
    b = signatures(*residues(f'{NEW}/{pdb[:2]}/{pdb}/{pdb}.pdb', ch))
    return [(pdb, ch, k[0], k[1], int(k in a and k in b and a[k] == b[k])) for k in sorted(a.keys() | b.keys())]


if __name__ == '__main__':
    pcs = [(r['pdb_id'], r['chain']) for r in csv.DictReader(open(f'{NEW}/Top8000-best_hom70_pdb_chain_list.csv'))]
    with Pool(8) as p: R = p.map(work, pcs, chunksize=20)
    with open('data/identical_inputs.tsv', 'w') as f:
        f.write('pdb\tchain\tresnum\tins\tidentical\n')
        for rows in R:
            for r in rows: f.write('\t'.join(map(str, r)) + '\n')
    n = sum(len(r) for r in R); s = sum(x[-1] for r in R for x in r)
    print(f'residues with N and CA: {n:,}; inputs identical: {s:,}; not identical: {n - s:,}')
