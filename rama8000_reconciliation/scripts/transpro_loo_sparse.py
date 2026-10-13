"""TransPro leave-one-out over residues in sparse regions (our grid fraction < 0.02), rule |ω|>=150."""
import sys, json, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.argv = ['x', 'data/pro_records_ext.jsonl', 'transpro', 'data/transpro.ndft', '8']
exec(open('scripts/class_rule_check.py').read().split("rules = {")[0])
F = [r for r in P if r['own_mc'] < 30 and r['own_minocc'] >= 1 and abs(r['omega']) >= 150]
g0 = silk(F); base = np.abs(g0 - T).mean() * 1e4
sparse = [i for i, r in enumerate(F) if g0[ndft.binof(r['phi']), ndft.binof(r['psi'])] < 0.02]
print(f'n={len(F)} base={base:.4f}; sparse residues tested: {len(sparse)}')
def loo(i):
    g = silk(F[:i] + F[i+1:]); return i, np.abs(g - T).mean() * 1e4, int(((g > 0) != (T > 0)).sum())
with ThreadPoolExecutor(8) as ex: res = list(ex.map(loo, sparse))
res.sort(key=lambda x: x[1])
json.dump([dict(F[i], delta=s - base, supp=m) for i, s, m in res], open('results/transpro_loo_sparse.json', 'w'))
print(f'improve: {sum(s < base for _, s, _ in res)}; best 25:')
for i, s, m in res[:25]:
    r = F[i]
    print(f"  {s-base:+.4f} supp={m:3d}  {r['pdb']} {r['chain']}{r['resnum']}{r['ins'].strip():1s} φ={r['phi']:7.1f} ψ={r['psi']:7.1f} ω={r['omega']:7.1f} "
          f"ownB={r['own_mc']:.1f} prev={r['resname_prev']} next={r['resname_next']}")
