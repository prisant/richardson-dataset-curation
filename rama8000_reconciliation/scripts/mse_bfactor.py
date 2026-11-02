"""Which atoms give the B value the lab printed for MSE residues in the 2013 lists (label "... MSE B<x ...")?

For every MSE in the lab's 2013 lists, read the residue's heavy-atom records from the lab's own 2011 chain file
(top8000_chains_70.tgz, <pdb>FH_<chain>.pdb, first model, all altlocs) and compare the printed B with candidate maxima.
A candidate matches when it is within 0.05 of the printed value (one decimal; ties such as 21.55 print either way). Hydrogens play no part.
Usage (from rama8000_reconciliation/):  python scripts/mse_bfactor.py > results/mse_bfactor.txt
"""
import os
from collections import Counter
exec(open('scripts/current_files_reconcile.py').read().split('\ndef load')[0])   # K, L (lab lists), LBL, FILES, LAB

OLD = os.path.expanduser('~/Desktop/Data/top8000_chains_70_2011/top8000_chains_70')
B = {}
for cls, f in FILES.items():
    for ln in open(f'{LAB}/residues/rama8000-{f}.csv'):
        m = LBL.match(ln.rstrip('\n').rsplit(',', 4)[0])
        if m.group(5) == 'MSE':
            B[(m.group(1), m.group(2), m.group(3), m.group(4).strip())] = float(m.group(6))


def atoms(pdb, ch, rn, ins):
    out = []
    for ln in open(f'{OLD}/{pdb}FH_{ch}.pdb'):
        if ln.startswith('ENDMDL'): break
        if (ln.startswith(('ATOM', 'HETATM')) and ln[21] == ch and ln[22:26].strip() == rn and ln[26].strip() == ins):
            name = ln[12:16].strip(); el = ln[76:78].strip() or name[0]
            if el not in ('H', 'D'): out.append((name, ln[16].strip(), float(ln[60:66])))
    return out


MC = {'N', 'CA', 'C', 'O'}
CAND = {'mainchain N CA C O': lambda a: [b for n, _, b in a if n in MC],
        'all heavy atoms': lambda a: [b for *_, b in a],
        'all heavy atoms except SE': lambda a: [b for n, _, b in a if n != 'SE'],
        'N CA C O CB': lambda a: [b for n, _, b in a if n in MC | {'CB'}],
        'SE only': lambda a: [b for n, _, b in a if n == 'SE'],
        'N CA C O SE': lambda a: [b for n, _, b in a if n in MC | {'SE'}],
        'N CA C O CB SE': lambda a: [b for n, _, b in a if n in MC | {'CB', 'SE'}]}
rows = []
for k, pb in B.items():
    a = atoms(*k)
    hit = {c: bool(f(a)) and abs(max(f(a)) - pb) <= 0.05 + 1e-9 for c, f in CAND.items()}
    allmax = max(b for *_, b in a)
    top = [n for n, _, b in a if b == allmax]
    rows.append((k, pb, hit, allmax, top, Counter(n for n, *_ in a)))

def report(sel, title):
    print(f'\n{title}: {len(sel):,} residues')
    for c in CAND: print(f'  printed B = max over {c:28s} {sum(r[2][c] for r in sel):6,d}')
    pat = Counter(tuple(c for c in CAND if r[2][c]) or ('none',) for r in sel)
    print('  combinations matched:')
    for p, n in pat.most_common(): print(f'    {n:6,d}  {" | ".join(p)}')
    print('  atom(s) holding the all-heavy-atom maximum:', dict(Counter('+'.join(sorted(set(r[4]))) for r in sel).most_common(8)))

report(rows, 'all MSE in the lab lists (2011 files)')
hi = [r for r in rows if r[3] >= 30]
report(hi, 'MSE the lab kept although the all-heavy-atom max B >= 30 (2011 files)')
print('\n  examples (printed B; heavy-atom B values):')
for k, pb, hit, m, top, _ in hi[:12]:
    print(f'    {" ".join(map(str, k)).strip():16s} B<{pb:<5}', {n + (alt and f".{alt}"): b for n, alt, b in atoms(*k)})
print('\n  atom-name sets of these residues:', dict(Counter(tuple(sorted(r[5])) for r in hi).most_common(4)))

none = [r for r in rows if not any(r[2].values())]
print(f'\nno candidate matches: {len(none)}; examples:')
for k, pb, *_ in none[:8]:
    print(f'    {" ".join(map(str, k)).strip():16s} B<{pb:<5}', {n + (alt and f".{alt}"): b for n, alt, b in atoms(*k)})
mcse = lambda a: max(b for n, _, b in a if n in MC | {'SE'})
low = [(k, pb) for k, pb in B.items() if max(b for n, _, b in atoms(*k) if n in MC) <= 1]
print(f'\nlab-listed MSE with mainchain max B <= 1 (fails 1 < B on N CA C O alone): {len(low)}')
for k, pb in low:
    print(f'    {" ".join(map(str, k)).strip():16s} B<{pb:<5} max(N CA C O SE) {mcse(atoms(*k)):.2f}')
