"""Which of our exclusions is too broad? Relax each one in turn under the full rule (incl. preceding-peptide trans).

Full rule: peptide-linked both sides (C-N < 2.0 Å), own mainchain 0 < B < 30, all heavy atoms occ 1, amide H occ 1,
standard residue, preceding peptide trans (|ω| ≥ 150) except CisPro (|ω| < 30).
After the trans rule: Gly +7, IleVal +109, PrePro -68, General -1,283, TransPro -1, CisPro 0.
Part 1: count Δ per class for each relaxation. Part 2: General / PrePro grids for the relaxations that add residues.
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
R = {  # each component of the full rule
    'link':     (d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0),
    'B<30':     d.own_mc_max < 30,
    'B>0':      d.own_mc_max > 0,
    'occ1':     d.own_minocc >= 1,
    'Hocc1':    ~(d.H_minocc < 1),
    'std':      d.resname.isin(STD),
    'trans':    np.where(d.rama6 == 'CisPro', aw < 30, aw.isna() | (aw >= 150)),
}
R['trans'] = pd.Series(R['trans'], index=d.index)
full = np.logical_and.reduce(list(R.values()))
nd = dict(General=1061639, Gly=116589, IleVal=195083, PrePro=61568, TransPro=66364, CisPro=3654)
relax = {'full rule': full}
for k in R:
    relax[f'relax {k}'] = np.logical_and.reduce([v for kk, v in R.items() if kk != k])
# narrower variants of the amide-H rule
relax['Hocc: drop only H occ == 0'] = np.logical_and.reduce([v for kk, v in R.items() if kk != 'Hocc1'] + [~(d.H_minocc == 0)])
relax['Hocc: drop only 0 < H occ < 1'] = np.logical_and.reduce([v for kk, v in R.items() if kk != 'Hocc1'] + [~((d.H_minocc < 1) & (d.H_minocc > 0))])
relax['link: 2.5 Å'] = np.logical_and.reduce([v for kk, v in R.items() if kk != 'link'] + [(d.d_prevC_N < 2.5) & (d.d_C_nextN < 2.5)])
relax['std + MSE'] = np.logical_and.reduce([v for kk, v in R.items() if kk != 'std'] + [d.resname.isin(STD | {'MSE'})])

print(f'{"variant":42s}' + ''.join(f'{c:>10s}' for c in nd) + f'{"total":>10s}')
for name, m in relax.items():
    dl = [int((m & (d.rama6 == c)).sum()) - n for c, n in nd.items()]
    print(f'{name:42s}' + ''.join(f'{x:>+10,d}' for x in dl) + f'{sum(dl):>+10,d}', flush=True)

print(f'\n{"class":8s}{"variant":42s}{"n":>9s}{"Δn":>7s}{"meanΔ×1e4":>11s}{"supp":>6s}{"ours-only":>10s}{"NDFT-only":>10s}{"rank-identical":>16s}')
for cls, f, tp in (('PrePro', 'prepro', '8'), ('General', 'general', '6.5')):
    rc, T = ndft.load(f'data/{f}.ndft')
    for name in ('full rule', 'relax Hocc1', 'Hocc: drop only H occ == 0', 'relax link', 'relax B>0', 'relax std', 'std + MSE'):
        S = d[relax[name] & (d.rama6 == cls)]
        n, g = silk_ct(S, '5', tp, False); dr = np.rint(g * rc) - np.rint(T * rc)
        print(f'{cls:8s}{name:42s}{len(S):9d}{len(S) - rc:+7d}{np.abs(g - T).mean() * 1e4:11.4f}{int(((g > 0) != (T > 0)).sum()):6d}'
              f'{int(((g > 0) & (T == 0)).sum()):10d}{int(((T > 0) & (g == 0)).sum()):10d}{int((dr == 0).sum()):>10d}/{dr.size}', flush=True)
