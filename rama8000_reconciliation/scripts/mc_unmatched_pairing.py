"""Can the N, CA, C, O records that match on no (residue number, insertion code, atom, altloc) key
(mc_atom_compare.py: 900 only in the 2011 files, 917 only in the current files) be paired another way?

Per chain, unmatched records are paired in order:
  1. same atom name, residue name and x/y/z (as written): renumbered residue or relabelled altloc, same position
  2. same residue number, insertion code and atom name; altloc label differs (coordinates may differ)
  3. same atom name and x/y/z, residue name differs
Whatever is left has no counterpart in the other file (residue added or removed, or both moved and renumbered).
Completeness: for each residue holding paired records (2011 and current residue number and insertion code, any altloc),
check that it has N, CA, C and O in each file and that none of its records is left unpaired.
Then all matched records (key-matched as in mc_atom_compare.py, plus these pairs) are compared field by field:
residue name, x/y/z, occupancy, B.
Optional argument 'noalt': leave out every record with an altloc label, in both file sets, before comparing.
Usage (from rama8000_reconciliation/):
  python scripts/mc_unmatched_pairing.py > results/mc_unmatched_pairing.txt
  python scripts/mc_unmatched_pairing.py noalt > results/mc_unmatched_pairing_noalt.txt
"""
import csv, os, sys
from collections import Counter
from multiprocessing import Pool
sys.path.insert(0, 'scripts')
from mc_atom_compare import OLD, NEW, MC, records   # key (resnum, ins, atom, altloc) -> (resname, x, y, z, occ, B)
NOALT = len(sys.argv) > 1 and sys.argv[1] == 'noalt'


def pair(o, n, key):
    """Pair records of o and n with equal key(k, v); remove them from o and n; return the (old, new) value pairs."""
    idx = {}
    for k, v in n.items(): idx.setdefault(key(k, v), []).append(k)
    out = []
    for k, v in list(o.items()):
        cand = idx.get(key(k, v))
        if cand:
            kn = cand.pop(); out.append((o.pop(k), n.pop(kn), k, kn))
    return out


def fields(c, pairs, tag):
    """Count identical / differing records and differing fields for (old, new) value pairs."""
    for a, b, *_ in pairs:
        d = [f for f, x in (('residue name', a[0] != b[0]), ('coordinates', a[1:4] != b[1:4]),
                            ('occupancy', float(a[4]) != float(b[4])), ('B', float(a[5]) != float(b[5]))) if x]
        c[f'{tag} identical' if not d else f'{tag} differs'] += 1
        for f in d: c[f'{tag} field: {f}'] += 1
        c[f'{tag} field: coordinates or B'] += ('coordinates' in d) or ('B' in d)


def work(pc):
    pdb, ch = pc
    o = records(f'{OLD}/{pdb}FH_{ch}.pdb', ch)
    n = records(f'{NEW}/{pdb[:2]}/{pdb}/{pdb}.pdb', ch)
    if NOALT:
        c0 = Counter(**{'altloc records left out, 2011': sum(k[3] != '' for k in o),
                        'altloc records left out, current': sum(k[3] != '' for k in n)})
        o = {k: v for k, v in o.items() if k[3] == ''}; n = {k: v for k, v in n.items() if k[3] == ''}
    O, N = dict(o), dict(n)
    shared = o.keys() & n.keys()
    c = c0 if NOALT else Counter()
    c['records 2011'], c['records current'] = len(o), len(n)
    fields(c, [(o[k], n[k]) for k in shared], 'key')
    o = {k: v for k, v in o.items() if k not in shared}
    n = {k: v for k, v in n.items() if k not in shared}
    c.update(**{'only 2011': len(o), 'only current': len(n)})
    allP = []
    for nm, key in (('1 same position, renumbered/relabelled', lambda k, v: (k[2], v[0]) + v[1:4]),
                    ('2 same residue and atom, altloc label differs', lambda k, v: k[:3]),
                    ('3 same position, residue name differs', lambda k, v: (k[2],) + v[1:4])):
        P = pair(o, n, key); c[nm] = len(P); fields(c, P, 'paired'); allP += P
    c['left 2011'], c['left current'] = len(o), len(n)
    # completeness per residue (residue number, insertion code; any altloc) for residues holding paired records
    G = {}
    for *_, ko, kn in allP: G.setdefault((ko[:2], kn[:2]), None)
    RO, RN = {}, {}
    for k in O: RO.setdefault(k[:2], set()).add(k[2])
    for k in N: RN.setdefault(k[:2], set()).add(k[2])
    names = lambda D, r: (RO if D is O else RN).get(r, set())
    # completeness of every matched residue (key-matched or paired): residue pair -> atom names matched
    M = {}
    for k in shared: M.setdefault((k[:2], k[:2]), set()).add(k[2])
    for *_, ko, kn in allP: M.setdefault((ko[:2], kn[:2]), set()).add(ko[2])
    for (ro, rn), at in M.items():
        c['matched residues'] += 1
        c['matched: 2011 residue has N, CA, C, O'] += RO[ro] == MC
        c['matched: current residue has N, CA, C, O'] += RN[rn] == MC
        c['matched: all of N, CA, C, O matched'] += at == MC
        if at != MC: c['matched, incomplete: ' + '+'.join(a for a in ('N', 'CA', 'C', 'O') if a in at)] += 1
    for ro, rn in G:
        c['residues'] += 1
        c['2011 residue has N, CA, C, O'] += names(O, ro) == MC
        c['current residue has N, CA, C, O'] += names(N, rn) == MC
        lo, ln = sorted(k[2:] for k in o if k[:2] == ro), sorted(k[2:] for k in n if k[:2] == rn)
        c['all its records matched (none left unpaired)'] += not (lo or ln)
        if names(O, ro) != MC or names(N, rn) != MC or lo or ln:
            c[f'exception: {pdb} {ch} 2011 {"".join(ro)} {sorted(names(O, ro))} -> current {"".join(rn)} {sorted(names(N, rn))};'
              f' unpaired 2011 {lo}, current {ln}'] += 1
    return pdb, c


if __name__ == '__main__':
    pcs = [(r['pdb_id'], r['chain']) for r in csv.DictReader(open(f'{NEW}/Top8000-best_hom70_pdb_chain_list.csv'))]
    with Pool(8) as p: R = p.map(work, pcs, chunksize=20)
    T = Counter()
    for _, c in R: T.update(c)
    if NOALT:
        print(f'altloc records left out: 2011 {T["altloc records left out, 2011"]:,}; current {T["altloc records left out, current"]:,}')
    print(f'records compared: 2011 {T["records 2011"]:,}; current {T["records current"]:,}')
    print(f'unmatched records: only 2011 {T["only 2011"]:,}; only current {T["only current"]:,}; total {T["only 2011"] + T["only current"]:,}')
    print('paired (each pair = one 2011 + one current record):')
    for k in sorted(k for k in T if k[0].isdigit() and k[1] == ' '):
        print(f'  {k[2:]:48s} {T[k]:5,d} pairs = {2 * T[k]:5,d} records')
    print(f'left without a counterpart: 2011 {T["left 2011"]:,}; current {T["left current"]:,}')
    E = Counter(); E2 = Counter()
    for pdb, c in R:
        if c['left 2011']: E[pdb] += c['left 2011']
        if c['left current']: E2[pdb] += c['left current']
    print(f'  2011 leftovers by entry ({len(E)}):', dict(E.most_common(12)))
    print(f'  current leftovers by entry ({len(E2)}):', dict(E2.most_common(12)))
    print('\ncompleteness of all matched residues (residues with at least one matched record; key-matched or paired):')
    for k in ('matched residues', 'matched: 2011 residue has N, CA, C, O', 'matched: current residue has N, CA, C, O',
              'matched: all of N, CA, C, O matched'):
        print(f'  {k[9:] if k.startswith("matched:") else k:46s} {T[k]:9,d}')
    print('  incomplete, by the atoms matched:', {k[21:]: T[k] for k in sorted(T, key=lambda k: -T[k]) if k.startswith('matched, incomplete: ')})
    print('\ncompleteness of the residues holding paired records (N, CA, C, O; any altloc unless noalt):')
    for k in ('residues', '2011 residue has N, CA, C, O', 'current residue has N, CA, C, O',
              'all its records matched (none left unpaired)'):
        print(f'  {k:46s} {T[k]:5,d}')
    for k in sorted(k for k in T if k.startswith('exception:')): print('   ', k[11:])
    print('\nmatched records compared field by field (a record may differ in several fields):')
    print(f'{"":28s}{"key-matched":>13s}{"paired":>10s}{"all":>13s}')
    for r in ('identical', 'differs'):
        a, b = T[f'key {r}'], T[f'paired {r}']
        print(f'{r:28s}{a:13,d}{b:10,d}{a + b:13,d}')
    a, b = sum(T[f'key {r}'] for r in ('identical', 'differs')), sum(T[f'paired {r}'] for r in ('identical', 'differs'))
    print(f'{"total matched":28s}{a:13,d}{b:10,d}{a + b:13,d}')
    for f in ('residue name', 'coordinates', 'occupancy', 'B', 'coordinates or B'):
        a, b = T[f'key field: {f}'], T[f'paired field: {f}']
        print(f'{"  differs in " + f:28s}{a:13,d}{b:10,d}{a + b:13,d}')
