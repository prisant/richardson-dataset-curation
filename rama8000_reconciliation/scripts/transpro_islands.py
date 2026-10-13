"""TransPro: clusters of ours-only / ref-only support at Makefile params; residues of our set inside each cluster."""
import sys, json, numpy as np
from scipy import ndimage
sys.argv = ['x', 'data/pro_records_ext.jsonl', 'transpro', 'data/transpro.ndft', '8']
exec(open('scripts/class_rule_check.py').read().split("rules = {")[0])
F = [r for r in P if r['own_mc'] < 30 and r['own_minocc'] >= 1]
g = silk(F)
np.save('results/transpro_rule_c5_t8.npy', g)
for name, m in [('OURS-ONLY', (g > 0) & (T == 0)), ('REF-ONLY', (g == 0) & (T > 0))]:
    lab, n = ndimage.label(m, structure=np.ones((3, 3)))
    sizes = ndimage.sum(m, lab, range(1, n + 1))
    print(f'\n{name}: {int(m.sum())} bins in {n} clusters (sizes ≥5 listed)')
    for k in np.argsort(-sizes):
        if sizes[k] < 5: break
        ii, jj = np.nonzero(lab == k + 1)
        pc, qc = ndft.center(int(ii.mean().round())), ndft.center(int(jj.mean().round()))
        inside = [r for r in F if lab[ndft.binof(r['phi']), ndft.binof(r['psi'])] == k + 1]
        near = [r for r in F if abs(((r['phi'] - pc + 180) % 360) - 180) < 12 and abs(((r['psi'] - qc + 180) % 360) - 180) < 12]
        print(f'  {int(sizes[k]):4d} bins centred ~({pc},{qc}); our residues inside: {len(inside)}; within ±12°: {len(near)}')
        for r in (inside or near)[:4]:
            print(f"      {r['pdb']} {r['chain']}{r['resnum']}{r['ins'].strip()} phi={r['phi']:7.1f} psi={r['psi']:7.1f} om={r['omega']:7.1f} "
                  f"ownB={r['own_mc']:.1f} prev={r['resname_prev']} next={r['resname_next']} Cn={r['d_prevC_N']:.2f}/{r['d_C_nextN']:.2f}")
