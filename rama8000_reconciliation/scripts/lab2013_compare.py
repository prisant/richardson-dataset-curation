"""Residue-level comparison: the lab's 2013 rama8000 residue lists vs our best rule on the top8000 files.

Lab lists: Top8000_Rama_2013.zip (recovered from the lab archive, October 2026), extracted to data/Top8000_Rama_2013/;
residues/rama8000-<class>.csv written 2013-02-16 by get_srcdata.py from the lab's top8000 MySQL database
(label "pdb chain resnum[ins] restype B<.. wc.. (res MP=..)", phi, psi, omega(-90..270), tau).
Best rule: see contour_compare.py. Per class: lab-only, ours-only, shared; shared residues in a different class;
phi/psi agreement on shared residues; for ours-only and lab-only, the filter values that could explain them.
Run from rama8000_reconciliation/.
"""
import sys, re
sys.path.insert(0, 'scripts')
import numpy as np, pandas as pd
sys.argv = ['x', 'Gly', 'data/glycine.ndft', '5', '8']
src = open('scripts/silk_param_scan.py').read()
exec(src[:src.index("print(f'{cls}: realcount")])   # d (with amide-H cols)
d['ins'] = d.ins.fillna('').astype(str).str.strip(); d['resnum'] = d.resnum.astype(int)

LAB = 'data/Top8000_Rama_2013'
FILES = dict(General='general-noGPIVpreP', IleVal='ileval-nopreP', Gly='gly', PrePro='prepro-noGP',
             TransPro='transpro', CisPro='cispro')
LBL = re.compile(r'^(\w{4}) (.) (-?\d+)(\D?) (\w{3}) B<([\d.]+) wc(\S+) \(([\d.]+)A MP=([\d.]+)\)$')
rows = []
for cls, f in FILES.items():
    for ln in open(f'{LAB}/residues/rama8000-{f}.csv'):
        lab, phi, psi, om, tau = ln.rstrip('\n').rsplit(',', 4)
        m = LBL.match(lab)
        if not m: sys.exit(f'unparsed label: {lab!r}')
        p, ch, rn, ins, rt, b, wc, res, mp = m.groups()
        rows.append((p, ch, int(rn), ins.strip(), rt, float(b), float(res), cls, float(phi), float(psi), float(om), float(tau)))
L = pd.DataFrame(rows, columns=['pdb', 'chain', 'resnum', 'ins', 'lab_resname', 'lab_mcB', 'lab_res', 'lab_cls',
                                'lab_phi', 'lab_psi', 'lab_omega', 'lab_tau'])
lab_chains = pd.read_csv(f'{LAB}/chains/rama8000_chains.tsv', sep='\t', dtype={'chain': str})
print(f'lab residues {len(L):,} in {L.groupby(["pdb", "chain"]).ngroups:,} chains (chains.tsv {len(lab_chains):,});'
      f' duplicate keys {int(L.duplicated(["pdb", "chain", "resnum", "ins"]).sum())}')

STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
aw = d.omega.abs()
rule = ((d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0) & (d.own_mc_max < 30) & (d.own_mc_max > 0) & (d.own_minocc >= 1)
        & ~(d.H_minocc < 1) & d.resname.isin(STD | {'MSE'})
        & pd.Series(np.where(d.rama6 == 'CisPro', aw < 30, aw.isna() | (aw >= 150)), index=d.index))
our_chains = set(zip(d.pdb, d.chain)); labc = set(zip(L.pdb, L.chain))
print(f'our listed chains {len(our_chains):,}; lab chains not in ours {len(labc - our_chains)};'
      f' our chains with a best-rule residue not in lab chains {len(set(zip(d[rule].pdb, d[rule].chain)) - labc)}')

K = ['pdb', 'chain', 'resnum', 'ins']
M = d.assign(ours=rule).merge(L, on=K, how='outer', indicator=True)
inL, inO = M._merge != 'left_only', M.ours.fillna(False).astype(bool)
print(f'lab residues not in our residue table at all: {int((M._merge == "right_only").sum()):,}')
sh = inL & inO
print(f'\nshared {int(sh.sum()):,}; lab-only {int((inL & ~inO).sum()):,}; ours-only {int((inO & ~inL).sum()):,}')
S = M[sh]
print(f'shared, resname differs {int((S.resname != S.lab_resname).sum())}; class differs {int((S.rama6 != S.lab_cls).sum())}')
dphi = np.abs((S.phi - S.lab_phi + 180) % 360 - 180); dpsi = np.abs((S.psi - S.lab_psi + 180) % 360 - 180)
print(f'shared |Δφ| max {dphi.max():.3f}, |Δψ| max {dpsi.max():.3f}; >0.1° {int(((dphi > 0.1) | (dpsi > 0.1)).sum()):,};'
      f' >1° {int(((dphi > 1) | (dpsi > 1)).sum()):,}')
print(f'shared, lab B<x vs our own_mc_max: max |Δ| {np.abs(S.lab_mcB - S.own_mc_max.round(1)).max():.2f}')

nd = dict(General=1061639, Gly=116589, IleVal=195083, PrePro=61568, TransPro=66364, CisPro=3654)
print(f'\n{"class":9s}{"lab":>10s}{"ours":>10s}{"shared":>10s}{"lab-only":>10s}{"ours-only":>10s}{"cls≠":>7s}')
for c in nd:
    lc, oc = inL & (M.lab_cls == c), inO & (M.rama6 == c)
    print(f'{c:9s}{int(lc.sum()):10,d}{int(oc.sum()):10,d}{int((lc & oc).sum()):10,d}{int((lc & ~inO).sum()):10,d}'
          f'{int((oc & ~inL).sum()):10,d}{int((lc & inO & ~oc).sum()):7d}')
print('\nshared residues in a different class (lab → ours):')
print(S[S.rama6 != S.lab_cls].groupby(['lab_cls', 'rama6']).size().to_string())

O = M[inO & ~inL].copy()
O['why'] = np.select([~O.set_index(['pdb', 'chain']).index.isin(list(labc)),
                      O.own_mc_max <= 1, ~O.resname.isin(STD)],
                     ['chain not in lab set', 'own mc B ≤ 1', 'MSE'], 'unexplained')
print('\nours-only by reason:'); print(O.groupby(['why', 'rama6']).size().unstack(fill_value=0).to_string())
Lo = M[inL & ~inO].copy()
inD = Lo._merge == 'both'
Lo['why'] = np.select([~inD, Lo.d_prevC_N.isna() | Lo.d_C_nextN.isna() | (Lo.d_prevC_N >= 2) | (Lo.d_C_nextN >= 2),
                       Lo.own_minocc < 1, Lo.H_minocc < 1, ~Lo.resname.isin(STD | {'MSE'}), Lo.own_mc_max >= 30,
                       Lo.own_mc_max <= 0],
                      ['not in our table', 'not peptide-linked', 'heavy-atom occ < 1', 'amide H occ < 1',
                       'non-standard', 'own mc B ≥ 30', 'own mc B ≤ 0'], 'ω rule')
print('\nlab-only by first failing rule of ours:'); print(Lo.groupby(['why', 'lab_cls']).size().unstack(fill_value=0).to_string())
print('\nlab-only non-standard residue types:', Lo[Lo.why == 'non-standard'].lab_resname.value_counts().to_dict())
print('ours-only unexplained, first 15:'); print(O[O.why == 'unexplained'][K + ['resname', 'rama6', 'own_mc_max', 'omega', 'H_minocc']].head(15).to_string())

# ---- Refined rules from the lab lists -------------------------------------------------------------------------------
# (1) phi ambiguity: drop a non-Pro residue when the preceding C has alternates at different positions
#     (data/prevC_spread.tsv, prevC_spread.py); the amide-H rule was a proxy for this.
# (2) MSE: the lab's B value is the max over all heavy atoms, not N, CA, C, O (data/allatom_records.tsv).
sp = pd.read_csv('data/prevC_spread.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
a = pd.read_csv('data/allatom_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
for t in (sp, a): t['ins'] = t.ins.str.strip()
E = M.merge(sp, on=K, how='left').merge(a[K + ['sc_maxB_heavy']], on=K, how='left')
bad = set(S[S.resname != S.lab_resname].pdb)
lab = E._merge != 'left_only'
ew = E.omega.abs()
base = ((E.d_prevC_N < 2.0) & (E.d_C_nextN < 2.0) & (E.own_mc_max < 30) & (E.own_mc_max > 0) & (E.own_minocc >= 1)
        & E.resname.isin(STD | {'MSE'})
        & pd.Series(np.where(E.rama6 == 'CisPro', ew < 30, ew.isna() | (ew >= 150)), index=E.index))
H = ~(E.H_minocc < 1)
phiamb = (E.prevC_spread > 0) & (E.resname != 'PRO')
mseB = ~((E.resname == 'MSE') & (np.fmax(E.own_mc_max, E.sc_maxB_heavy) >= 30))
V = {'best rule (amide-H rule)': base & H,
     'best + MSE all-heavy B < 30': base & H & mseB,
     'prev-C rule instead of H rule': base & ~phiamb,
     'prev-C rule + MSE B': base & ~phiamb & mseB,
     'prev-C rule + H rule + MSE B': base & ~phiamb & H & mseB}
print(f'\nrefined rules (Δ vs NDFT per class; ours-only / lab-only residues, excluding the {len(bad)} renumbered entries {sorted(bad)})')
print(f'{"variant":32s}' + ''.join(f'{c:>9s}' for c in nd) + f'{"total":>8s}{"ours-only":>11s}{"lab-only":>10s}')
ok = ~E.pdb.isin(bad)
for k, v in V.items():
    dl = [int((v & (E.rama6 == c)).sum()) - n for c, n in nd.items()]
    print(f'{k:32s}' + ''.join(f'{x:>+9,d}' for x in dl) + f'{sum(dl):>+8,d}{int((v & ~lab & ok).sum()):11,d}{int((lab & ~v & ok).sum()):10,d}')
best = V['prev-C rule + H rule + MSE B']
for nm, sel, cls in (('ours-only', best & ~lab & ok, 'rama6'), ('lab-only', lab & ~best & ok, 'lab_cls')):
    x = E[sel]
    print(f'\n{nm} under prev-C + H + MSE B, by class:', x[cls].value_counts().to_dict())
    print(f'{nm} by entry (top 10):', x.pdb.value_counts().head(10).to_dict())
S2 = S[~S.pdb.isin(bad)]
dphi = np.abs((S2.phi - S2.lab_phi + 180) % 360 - 180); dpsi = np.abs((S2.psi - S2.lab_psi + 180) % 360 - 180)
print(f'\nshared outside the renumbered entries: {len(S2):,}; resname differs {int((S2.resname != S2.lab_resname).sum())};'
      f' class differs {int((S2.rama6 != S2.lab_cls).sum())}; |Δφ| max {dphi.max():.3f}, |Δψ| max {dpsi.max():.3f};'
      f' >0.01° {int(((dphi > 0.01) | (dpsi > 0.01)).sum())}')
print('class differences there (lab → ours):', S2[S2.rama6 != S2.lab_cls].groupby(['lab_cls', 'rama6']).size().to_dict())
big = S2[(dphi > 0.01) | (dpsi > 0.01)]
print(f'φ/ψ differences > 0.01°: {len(big)} residues in {big.pdb.nunique()} entries; by entry:', big.pdb.value_counts().head(12).to_dict())
print('  |Δ| quantiles (max of φ, ψ):', np.percentile(np.maximum(dphi, dpsi)[big.index], [50, 90, 100]).round(2).tolist())

print('\nabsolute class counts: lab lists (= NDFT realcounts) vs ours under each rule')
print(f'{"rule":32s}' + ''.join(f'{c:>11s}' for c in nd) + f'{"total":>11s}')
print(f'{"lab 2013 lists":32s}' + ''.join(f'{n:11,d}' for n in nd.values()) + f'{sum(nd.values()):11,d}')
for k, v in V.items():
    n = [int((v & (E.rama6 == c)).sum()) for c in nd]
    print(f'{k:32s}' + ''.join(f'{x:11,d}' for x in n) + f'{sum(n):11,d}')

# ---- Evidence for the refined rules (outside the renumbered entries) ------------------------------------------------
pa = pd.read_csv('data/phipsi_atom_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
pa['ins'] = pa.ins.str.strip()
F = E.merge(pa, on=K, how='left').merge(a[K + ['all_minocc', 'sc_maxB_all']], on=K, how='left')
flab = F._merge != 'left_only'; fok = ~F.pdb.isin(bad)
fbase = pd.Series(base.values, index=F.index)
prevmc = sum((F[f'prev{x}_alt'] > 0) for x in ('N', 'CA', 'C', 'O')) > 0
T = F[fbase & fok & prevmc & (F.resname != 'PRO')]
tl = flab[T.index]
print(f'\nbase rule (best rule without the H rule), non-Pro, preceding residue with a mainchain altloc: {len(T):,}; lab kept {int(tl.sum()):,}')
for nm, m in (('preceding C alternates at different positions', T.prevC_spread > 0),
              ('preceding C alternates coincide', T.prevC_spread == 0),
              ('  of which CA or O alternates differ', (T.prevC_spread == 0) & ((T.prevCA_spread > 0) | (T.prevO_spread > 0)))):
    print(f'  {nm}: {int(m.sum()):,}; lab kept {int(tl[m].sum()):,}, dropped {int((~tl[m]).sum()):,}')
x = T[T.pdb == '1agy']
print(f'  1agy: {len(x)} such residues, lab kept {int(flab[x.index].sum())}')
y = F[best & flab & fok]
print(f'shared under prev-C + H + MSE B with a non-amide H (or any atom) below occupancy 1: {int((y.all_minocc < 1).sum()):,}')
z = F[fbase & fok & (F.resname == 'MSE')].merge(L[K + ['lab_mcB']], on=K, how='left', suffixes=('', '_l'))
zl = z.lab_mcB_l.notna()
print(f'MSE under the base rule: {len(z):,}; lab kept {int(zl.sum()):,}; dropped {int((~zl).sum()):,},'
      f' of which side-chain heavy B >= 30: {int(((~zl) & (z.sc_maxB_heavy >= 30)).sum()):,}')
zk = z[zl]
for nm, v in (('N, CA, C, O', zk.own_mc_max), ('all heavy atoms', np.fmax(zk.own_mc_max, zk.sc_maxB_heavy))):
    print(f'  lab printed B matches max over {nm} (±0.05): {int((np.abs(v - zk.lab_mcB_l) <= 0.051).sum()):,} / {len(zk):,}')
mt = F[fbase & fok & flab & (F.resname == 'MET')]
print(f'  MET: lab printed B matches max over N, CA, C, O: {int((np.abs(mt.own_mc_max - mt.lab_mcB) <= 0.051).sum()):,} / {len(mt):,}')
