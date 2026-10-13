"""Leave-one-out over set S (own mc<30, no alt) and add-one over remaining candidates, scored vs cispro.ndft."""
import sys, json, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, 'scripts')
import cis_group_test as G   # reuses score(), data (prints its table once)
S = [r for r in G.cis if r['own_mc'] < 30 and not r['alt']]
rest = [r for r in G.cis if r not in S]
base = G.score(S)
print(f'\nset S: n={len(S)} score={base}')
loo = lambda i: G.score(S[:i] + S[i+1:])[0]
add = lambda r: G.score(S + [r])[0]
with ThreadPoolExecutor(8) as ex:
    d_loo = list(ex.map(loo, range(len(S))))
    d_add = list(ex.map(add, rest))
out = []
for r, s in zip(S, d_loo): out.append(dict(kind='remove', delta=s - base[0], **r))
for r, s in zip(rest, d_add): out.append(dict(kind='add', delta=s - base[0], **r))
json.dump(out, open('results/cispro_single_tests.json', 'w'))
for kind in ('remove', 'add'):
    xs = sorted([o for o in out if o['kind'] == kind], key=lambda o: o['delta'])
    print(f'\n{kind}: {sum(o["delta"] < 0 for o in xs)} of {len(xs)} improve the fit; best 15:')
    for o in xs[:15]:
        print(f"  {o['delta']:+.4f}  {o['pdb']} {o['chain']}{o['resnum']}{o['ins'].strip():1s} phi={o['phi']:7.1f} psi={o['psi']:7.1f} "
              f"om={o['omega']:6.1f} ownB={o['own_mc']:.1f} prevC={o['prevC']:.1f} nextN={o['nextN']:.1f} alt={''.join(o['alt']) or '-'}")
