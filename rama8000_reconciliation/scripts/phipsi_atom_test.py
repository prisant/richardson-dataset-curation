"""Replace the Reduce amide-H rule with an occupancy / altloc rule on the neighbour mainchain atoms of phi, psi, omega.

Hypothesis: the amide-H rule (H occ 1) is a proxy for 'every atom defining phi/psi is fully occupied, no altloc',
since Reduce splits the amide H only when the preceding C has alternates. Base: best rule (std + MSE) without the H rule.
Atom data from data/phipsi_atom_records.tsv (phipsi_atom_records.py); own atoms are already covered by own_minocc.
Run from rama8000_reconciliation/.
"""
import sys
sys.path.insert(0, 'scripts')
import numpy as np, pandas as pd
sys.argv = ['x', 'Gly', 'data/glycine.ndft', '5', '8']
src = open('scripts/silk_param_scan.py').read()
exec(src[:src.index("print(f'{cls}: realcount")])   # d (with amide-H cols), silk_ct
import ndft
a = pd.read_csv('data/phipsi_atom_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
d = d.merge(a, on=['pdb', 'chain', 'resnum', 'ins'], how='left')
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
aw = d.omega.abs()
base = ((d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0) & (d.own_mc_max < 30) & (d.own_mc_max > 0) & (d.own_minocc >= 1)
        & d.resname.isin(STD | {'MSE'}) & pd.Series(np.where(d.rama6 == 'CisPro', aw < 30, aw.isna() | (aw >= 150)), index=d.index))
ok = lambda w: (d[f'{w}_occ'] >= 1) & (d[f'{w}_alt'] == 0)   # missing atom (nan occ) fails
okn = lambda w: ~(d[f'{w}_occ'] < 1) & (d[f'{w}_alt'] == 0)  # missing atom passes
H = ~(d.H_minocc < 1)
V = {'H rule (current best)': H,
     'no H rule': True,
     'prev C': ok('prevC'),
     'prev C + next N (phi/psi)': ok('prevC') & ok('nextN'),
     'prev CA, C + next N (+omega)': ok('prevCA') & ok('prevC') & ok('nextN'),
     'prev N, CA, C, O + next N, CA': ok('prevN') & ok('prevCA') & ok('prevC') & okn('prevO') & ok('nextN') & ok('nextCA'),
     'H rule + next N': H & ok('nextN'),
     'H rule + prev CA + next N': H & ok('prevCA') & ok('nextN')}
print('H rule vs prev-C rule (base, all classes): '
      f'H fails {int((base & ~H).sum())}; prev C fails {int((base & ~ok("prevC")).sum())}; '
      f'H fails & prev C ok {int((base & ~H & ok("prevC")).sum())}; H ok & prev C fails {int((base & H & ~ok("prevC")).sum())}')
nd = dict(General=1061639, Gly=116589, IleVal=195083, PrePro=61568, TransPro=66364, CisPro=3654)
print(f'\n{"variant":32s}' + ''.join(f'{c:>10s}' for c in nd) + f'{"total":>10s}')
for name, m in V.items():
    dl = [int((base & m & (d.rama6 == c)).sum()) - n for c, n in nd.items()]
    print(f'{name:32s}' + ''.join(f'{x:>+10,d}' for x in dl) + f'{sum(dl):>+10,d}', flush=True)
print(f'\n{"class":9s}{"variant":32s}{"n":>9s}{"Δn":>7s}{"meanΔ×1e4":>11s}{"supp":>6s}{"ours-only":>10s}{"NDFT-only":>10s}{"rank-identical":>16s}{"max|Δrank|":>11s}')
for cls, f, tp in (('CisPro', 'cispro', '8'), ('TransPro', 'transpro', '8'), ('Gly', 'glycine', '8'),
                   ('IleVal', 'ileval', '8'), ('PrePro', 'prepro', '8'), ('General', 'general', '6.5')):
    rc, T = ndft.load(f'data/{f}.ndft'); sym = cls == 'Gly'; nref = rc // 2 if sym else rc
    for name, m in V.items():
        if name == 'no H rule': continue
        S = d[base & m & (d.rama6 == cls)]
        n, g = silk_ct(S, '5', tp, sym); dr = np.rint(g * rc) - np.rint(T * rc)
        print(f'{cls:9s}{name:32s}{len(S):9d}{len(S) - nref:+7d}{np.abs(g - T).mean() * 1e4:11.4f}{int(((g > 0) != (T > 0)).sum()):6d}'
              f'{int(((g > 0) & (T == 0)).sum()):10d}{int(((T > 0) & (g == 0)).sum()):10d}{int((dr == 0).sum()):>10d}/{dr.size}{int(np.abs(dr).max()):11d}', flush=True)
