"""Non-standard residue types that pass the rest of the best rule, by type and class, and the General / PrePro grid
when each frequent type is added on top of standard + MSE.
Run from rama8000_reconciliation/.
"""
import sys
sys.path.insert(0, 'scripts')
import numpy as np, pandas as pd
sys.argv = ['x', 'Gly', 'data/glycine.ndft', '5', '8']
src = open('scripts/silk_param_scan.py').read()
exec(src[:src.index("print(f'{cls}: realcount")])   # d (with amide-H cols), silk_ct
import ndft
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
aw = d.omega.abs()
rule = ((d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0) & (d.own_mc_max < 30) & (d.own_mc_max > 0) & (d.own_minocc >= 1)
        & ~(d.H_minocc < 1) & pd.Series(np.where(d.rama6 == 'CisPro', aw < 30, aw.isna() | (aw >= 150)), index=d.index))
ns = d[rule & ~d.resname.isin(STD)]
print(f'non-standard residues passing the rest of the rule: {len(ns)}; types: {ns.resname.nunique()}')
print(ns.groupby(['resname', 'rama6']).size().unstack(fill_value=0).assign(total=lambda x: x.sum(1)).sort_values('total', ascending=False).head(25).to_string())
top = [t for t in ns.resname.value_counts().index if t != 'MSE'][:8]
print(f'\n{"class":8s}{"added to std+MSE":22s}{"n":>9s}{"Δn":>7s}{"meanΔ×1e4":>11s}{"rank-identical":>16s}')
for cls, f, tp in (('General', 'general', '6.5'), ('PrePro', 'prepro', '8')):
    rc, T = ndft.load(f'data/{f}.ndft')
    for name, extra in [('(none)', set()), ('all non-standard', set(ns.resname))] + [(t, {t}) for t in top]:
        S = d[rule & d.resname.isin(STD | {'MSE'} | extra) & (d.rama6 == cls)]
        n, g = silk_ct(S, '5', tp, False); dr = np.rint(g * rc) - np.rint(T * rc)
        print(f'{cls:8s}{name:22s}{len(S):9d}{len(S) - rc:+7d}{np.abs(g - T).mean() * 1e4:11.4f}{int((dr == 0).sum()):>10d}/{dr.size}', flush=True)
