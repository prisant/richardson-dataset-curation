"""Per-class summary under the best rule so far, every class with every metric.

Best rule: peptide-linked both sides (C-N < 2.0 Å); own mainchain 0 < B < 30; all heavy atoms occ 1; amide H occ 1;
standard residues + MSE; preceding peptide trans (|ω| ≥ 150), except CisPro (|ω| < 30). Read 2011 precedence.
Silk cosine 5, twopass 8 (General 6.5); Gly mirrored. Also: the standard-only rule; the best rule without the amide-H
rule; and the best rule with the amide-H rule replaced by an explicit preceding-residue rule (previous residue has no
altloc and all atoms occ 1). Part 1 checks why a Reduce amide H has occ < 1 when the residue's own heavy atoms are occ 1.
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
same_prev = (d.pdb.shift(1) == d.pdb) & (d.chain.shift(1) == d.chain) & (d.d_prevC_N < 2.0)
prev_alt = d.own_alt.shift(1).where(same_prev); prev_occ = d.own_minocc.shift(1).where(same_prev)
prev_bad = (prev_alt > 0) | (prev_occ < 1)
hlow = (d.own_minocc >= 1) & (d.H_minocc < 1)
print(f'own heavy atoms occ 1 and amide H occ < 1: {int(hlow.sum())} residues')
print(f'  H occ == 0: {int((hlow & (d.H_minocc == 0)).sum())}; 0 < H occ < 1: {int((hlow & (d.H_minocc > 0)).sum())}')
print(f'  of 0 < H occ < 1, previous residue has altloc or occ < 1: {int((hlow & (d.H_minocc > 0) & prev_bad).sum())}; '
      f'previous residue not in records (no phi/psi row): {int((hlow & (d.H_minocc > 0) & prev_alt.isna()).sum())}')
print(f'previous residue has altloc or occ < 1, own heavy occ 1, H occ 1: {int(((d.own_minocc >= 1) & prev_bad & ~(d.H_minocc < 1)).sum())}\n')
rule = ((d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0) & (d.own_mc_max < 30) & (d.own_mc_max > 0) & (d.own_minocc >= 1)
        & pd.Series(np.where(d.rama6 == 'CisPro', aw < 30, aw.isna() | (aw >= 150)), index=d.index))
setup = [('CisPro', 'cispro', '8'), ('TransPro', 'transpro', '8'), ('Gly', 'glycine', '8'),
         ('IleVal', 'ileval', '8'), ('PrePro', 'prepro', '8'), ('General', 'general', '6.5')]
print(f'{"rule":10s}{"class":9s}{"NDFT n":>10s}{"ours":>10s}{"Δn":>8s}{"Δ%":>8s}{"meanΔ×1e4":>11s}{"supp":>6s}{"ours-only":>10s}{"NDFT-only":>10s}'
      f'{"rank-identical":>16s}{"max|Δrank|":>11s}')
hrule = ~(d.H_minocc < 1)
for rname, res, extra in (('std+MSE', STD | {'MSE'}, hrule), ('std only', STD, hrule), ('no H rule', STD | {'MSE'}, True),
                          ('prev rule', STD | {'MSE'}, ~prev_bad)):
    for cls, f, tp in setup:
        rc, T = ndft.load(f'data/{f}.ndft'); sym = cls == 'Gly'; nref = rc // 2 if sym else rc
        S = d[rule & extra & d.resname.isin(res) & (d.rama6 == cls)]
        n, g = silk_ct(S, '5', tp, sym); dr = np.rint(g * rc) - np.rint(T * rc)
        print(f'{rname:10s}{cls:9s}{nref:10d}{len(S):10d}{len(S) - nref:+8d}{100 * (len(S) - nref) / nref:+8.3f}{np.abs(g - T).mean() * 1e4:11.4f}'
              f'{int(((g > 0) != (T > 0)).sum()):6d}{int(((g > 0) & (T == 0)).sum()):10d}{int(((T > 0) & (g == 0)).sum()):10d}'
              f'{int((dr == 0).sum()):>10d}/{dr.size}{int(np.abs(dr).max()):11d}', flush=True)
