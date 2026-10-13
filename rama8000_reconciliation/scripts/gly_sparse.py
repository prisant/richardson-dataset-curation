"""Sparse-region Gly: which of our Gly residues the lab's set lacks, and what they share.

Our Gly set (reconciled rule, B < 30) has +76 residues over glycine.ndft and 2,084 support bins the NDFT lacks
(support_diag.py). Here:
  1. our grid (mirrored, Silk cosine 5 / twopass 8) vs the NDFT at each residue's own bin;
  2. property profile of the sparsest residues vs all;
  3. count + grid for candidate exclusions.
Neighbour values come from the adjacent row of residue_records.tsv (same chain, peptide-linked); residues whose
neighbour has no phi/psi row get NaN and are never excluded by a neighbour rule.
Run from rama8000_reconciliation/.
"""
import sys
sys.path.insert(0, 'scripts')
import numpy as np, pandas as pd
sys.argv = ['x', 'Gly', 'data/glycine.ndft', '5', '8']
src = open('scripts/silk_param_scan.py').read()
exec(src[:src.index("print(f'{cls}: realcount")])   # d, common, rc, T, silk_ct
import ndft

# neighbour columns from adjacent rows (file is in chain order)
same_next = (d.pdb.shift(-1) == d.pdb) & (d.chain.shift(-1) == d.chain) & (d.d_C_nextN < 2.0)
same_prev = (d.pdb.shift(1) == d.pdb) & (d.chain.shift(1) == d.chain) & (d.d_prevC_N < 2.0)
d['next_omega'] = d.omega.shift(-1).where(same_next)
d['next_mc_max'] = d.own_mc_max.shift(-1).where(same_next)
d['prev_mc_max'] = d.own_mc_max.shift(1).where(same_prev)
d['next_minocc'] = d.own_minocc.shift(-1).where(same_next)
d['prev_minocc'] = d.own_minocc.shift(1).where(same_prev)

base = common & (d.own_mc_max < 30) & (d.rama6 == 'Gly')
G = d[base].copy()
n, g = silk_ct(G, '5', '8', True)
bi = np.array([ndft.binof(x) for x in G.phi]); bj = np.array([ndft.binof(x) for x in G.psi])
G['ours'] = g[bi, bj] * rc; G['ndft'] = T[bi, bj] * rc
G['ratio'] = G.ndft / G.ours
print(f'Gly set {len(G)} (x2 mirrored = {n}); NDFT realcount {rc} (= {rc // 2} residues); Δ {len(G) - rc // 2:+d}')

def profile(df):
    aw = df.omega.abs(); an = df.next_omega.abs()
    return {'n': len(df),
            'prev ω non-trans %': 100 * (aw < 150).mean(),
            'next ω non-trans %': 100 * (an < 150).mean(),
            'next ω missing %': 100 * an.isna().mean(),
            'own B median': df.own_mc_max.median(),
            'prev B≥30 %': 100 * (df.prev_mc_max >= 30).mean(),
            'next B≥30 %': 100 * (df.next_mc_max >= 30).mean(),
            'prev occ<1 %': 100 * (df.prev_minocc < 1).mean(),
            'next occ<1 %': 100 * (df.next_minocc < 1).mean(),
            'next is PRO %': 100 * (df.next_resname == 'PRO').mean()}

sparse = G.nsmallest(300, 'ours')
print('\nProfile: all Gly vs the 300 sparsest (lowest our-grid density at own bin)')
pa, ps = profile(G), profile(sparse)
for k in pa: print(f'  {k:22s} {pa[k]:>10.2f} {ps[k]:>10.2f}')
print('\nLowest NDFT/ours ratio among sparse residues (NDFT density much below ours at the point):')
print(sparse.nsmallest(25, 'ratio')[['pdb', 'chain', 'resnum', 'ins', 'phi', 'psi', 'omega', 'next_omega',
                                      'own_mc_max', 'prev_mc_max', 'next_mc_max', 'ours', 'ndft']].round(2).to_string(index=False))

aw, an = G.omega.abs(), G.next_omega.abs()
variants = {
    'reconciled rule':                      pd.Series(True, index=G.index),
    'drop prev ω non-trans (|ω|<150)':      ~(aw < 150),
    'drop next ω non-trans':                ~(an < 150),
    'drop prev or next non-trans':          ~((aw < 150) | (an < 150)),
    'drop neighbour mc B ≥ 30':             ~((G.prev_mc_max >= 30) | (G.next_mc_max >= 30)),
    'drop neighbour occ < 1':               ~((G.prev_minocc < 1) | (G.next_minocc < 1)),
}
print(f'\n{"variant":34s}{"n":>8s}{"Δn":>6s}{"meanΔ×1e4":>11s}{"supp":>6s}{"ours-only":>10s}{"rank-identical":>16s}')
for name, m in variants.items():
    S = G[m]; n2, g2 = silk_ct(S, '5', '8', True); dr = np.rint(g2 * rc) - np.rint(T * rc)
    print(f'{name:34s}{len(S):8d}{len(S) - rc // 2:+6d}{np.abs(g2 - T).mean() * 1e4:11.4f}'
          f'{int(((g2 > 0) != (T > 0)).sum()):6d}{int(((g2 > 0) & (T == 0)).sum()):10d}{int((dr == 0).sum()):>10d}/{dr.size}', flush=True)
