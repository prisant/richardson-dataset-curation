"""Preceding-peptide (omega) exclusion for the non-Pro classes: cis (|ω|<30), twisted (30≤|ω|<150), or both.

gly_sparse.py found that dropping Gly whose preceding peptide is non-trans takes Gly from +76 to +7 and the
support mismatch from 2,084 bins to 2. This tests the split and the other non-Pro classes, at their Makefile
smoothing (General twopass 6.5; Gly, PrePro, IleVal 8; cosine 5). Reconciled rule, B < 30.
usage: nontrans_nonpro.py [CLASS ...]   Run from rama8000_reconciliation/.
"""
import sys
sys.path.insert(0, 'scripts')
import numpy as np
classes = sys.argv[1:] or ['Gly', 'PrePro', 'IleVal', 'General']
sys.argv = ['x', 'Gly', 'data/glycine.ndft', '5', '8']
src = open('scripts/silk_param_scan.py').read()
exec(src[:src.index("print(f'{cls}: realcount")])   # d, common, silk_ct
import ndft
setup = {'Gly': ('glycine', '8'), 'PrePro': ('prepro', '8'), 'IleVal': ('ileval', '8'), 'General': ('general', '6.5')}
aw = d.omega.abs()
variants = {'reconciled rule': aw.isna() | (aw >= 0), 'drop cis (|ω|<30)': ~(aw < 30), 'drop twisted (30–150)': ~((aw >= 30) & (aw < 150)),
            'drop non-trans (|ω|<150)': ~(aw < 150)}
print(f'{"class":8s}{"variant":28s}{"n":>9s}{"Δn":>7s}{"meanΔ×1e4":>11s}{"supp":>6s}{"ours-only":>10s}{"NDFT-only":>10s}{"rank-identical":>16s}{"max|Δrank|":>11s}')
for cls in classes:
    f, tp = setup[cls]; rc, T = ndft.load(f'data/{f}.ndft'); sym = cls == 'Gly'
    base = common & (d.own_mc_max < 30) & (d.rama6 == cls)
    for name, m in variants.items():
        S = d[base & m]; n, g = silk_ct(S, '5', tp, sym); dr = np.rint(g * rc) - np.rint(T * rc)
        print(f'{cls:8s}{name:28s}{len(S):9d}{len(S) - (rc // 2 if sym else rc):+7d}{np.abs(g - T).mean() * 1e4:11.4f}'
              f'{int(((g > 0) != (T > 0)).sum()):6d}{int(((g > 0) & (T == 0)).sum()):10d}{int(((T > 0) & (g == 0)).sum()):10d}'
              f'{int((dr == 0).sum()):>10d}/{dr.size}{int(np.abs(dr).max()):11d}', flush=True)
