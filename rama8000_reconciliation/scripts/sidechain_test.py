"""MSE inclusion plus side-chain exclusions (incomplete side chains, side-chain B).

relax_exclusions.py: adding MSE to the standard residues gives the best General/PrePro grids but overshoots
(General +1,552, PrePro +82). Candidate extra exclusions, from data/sidechain_records.tsv
(sidechain_records.py): any expected heavy atom missing; CB missing; side-chain B ≥ 30 (i.e. all-atom B < 30).
Part 1: count Δ per class. Part 2: grids for General (twopass 6.5), PrePro, IleVal (8).
Run from rama8000_reconciliation/.
"""
import sys
sys.path.insert(0, 'scripts')
import numpy as np, pandas as pd
sys.argv = ['x', 'Gly', 'data/glycine.ndft', '5', '8']
src = open('scripts/silk_param_scan.py').read()
exec(src[:src.index("print(f'{cls}: realcount")])   # d (with amide-H cols), silk_ct
import ndft
s = pd.read_csv('data/sidechain_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str, 'missing': str})
d = d.merge(s.drop(columns='resname'), on=['pdb', 'chain', 'resnum', 'ins'], how='left')
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
aw = d.omega.abs()
rule = ((d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0) & (d.own_mc_max < 30) & (d.own_mc_max > 0) & (d.own_minocc >= 1)
        & ~(d.H_minocc < 1) & pd.Series(np.where(d.rama6 == 'CisPro', aw < 30, aw.isna() | (aw >= 150)), index=d.index))
std, mse = d.resname.isin(STD), d.resname.isin(STD | {'MSE'})
miss_any = d.n_missing > 0
miss_O = d.missing.fillna('').str.split(',').apply(lambda x: 'O' in x)
miss_cb = d.missing.fillna('').str.split(',').apply(lambda x: 'CB' in x)
scB = d.sc_maxB >= 30
V = {'full rule (std only)': rule & std,
     'std + MSE': rule & mse,
     '+MSE, drop any atom missing': rule & mse & ~miss_any,
     '+MSE, drop O missing': rule & mse & ~miss_O,
     '+MSE, drop CB missing': rule & mse & ~miss_cb,
     '+MSE, drop side-chain B ≥ 30': rule & mse & ~scB,
     '+MSE, drop missing or sc B ≥ 30': rule & mse & ~miss_any & ~scB}
nd = dict(General=1061639, Gly=116589, IleVal=195083, PrePro=61568, TransPro=66364, CisPro=3654)
print('incomplete residues under std + MSE, by class:')
print(d[rule & mse & miss_any].groupby('rama6').size().to_string())
print('most common missing-atom sets (std + MSE):')
print(d[rule & mse & miss_any].missing.value_counts().head(10).to_string())
print(f'\n{"variant":38s}' + ''.join(f'{c:>10s}' for c in nd) + f'{"total":>10s}')
for name, m in V.items():
    dl = [int((m & (d.rama6 == c)).sum()) - n for c, n in nd.items()]
    print(f'{name:38s}' + ''.join(f'{x:>+10,d}' for x in dl) + f'{sum(dl):>+10,d}', flush=True)

print(f'\n{"class":8s}{"variant":38s}{"n":>9s}{"Δn":>7s}{"meanΔ×1e4":>11s}{"supp":>6s}{"ours-only":>10s}{"NDFT-only":>10s}{"rank-identical":>16s}')
for cls, f, tp in (('General', 'general', '6.5'), ('PrePro', 'prepro', '8'), ('IleVal', 'ileval', '8')):
    rc, T = ndft.load(f'data/{f}.ndft')
    for name, m in V.items():
        S = d[m & (d.rama6 == cls)]
        n, g = silk_ct(S, '5', tp, False); dr = np.rint(g * rc) - np.rint(T * rc)
        print(f'{cls:8s}{name:38s}{len(S):9d}{len(S) - rc:+7d}{np.abs(g - T).mean() * 1e4:11.4f}{int(((g > 0) != (T > 0)).sum()):6d}'
              f'{int(((g > 0) & (T == 0)).sum()):10d}{int(((T > 0) & (g == 0)).sum()):10d}{int((dr == 0).sum()):>10d}/{dr.size}', flush=True)
