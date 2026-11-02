"""Residue-level differences between the lab's 2011 chain files (top8000_chains_70.tgz, Reduce 3.15)
and our current RCSB-derived files (~/Desktop/Data/top8000_pdbs_hom70, Reduce 4.16) for the same chains.

Inputs: data/{residue_records,amideH_records,allatom_records}.tsv (current files) and the same under data/lab2011/.
Usage: python scripts/files2011_vs_current.py > results/files2011_vs_current.txt
"""
import numpy as np, pandas as pd

K = ['pdb', 'chain', 'resnum', 'ins']


def load(d):
    rd = lambda f: pd.read_csv(f'{d}/{f}', sep='\t', keep_default_na=False, na_values=['nan'],
                               dtype={'ins': str, 'chain': str, 'pdb': str})
    r = rd('residue_records.tsv').merge(rd('amideH_records.tsv'), on=K, how='left').merge(rd('allatom_records.tsv'), on=K, how='left')
    r['ins'] = r.ins.str.strip()
    return r


N, O = load('data'), load('data/lab2011')
print(f'residues: 2011 files {len(O):,}; current files {len(N):,}')
print(f'chains:   2011 files {O.groupby(["pdb", "chain"]).ngroups:,}; current files {N.groupby(["pdb", "chain"]).ngroups:,}')
M = O.merge(N, on=K, how='outer', suffixes=('_o', '_n'), indicator=True)
print(f'same residue key in both: {int((M._merge == "both").sum()):,}; only 2011: {int((M._merge == "left_only").sum()):,}; '
      f'only current: {int((M._merge == "right_only").sum()):,}')
ent = lambda m: M[m].groupby('pdb').size().sort_values(ascending=False)
for nm, m in (('only 2011', M._merge == 'left_only'), ('only current', M._merge == 'right_only')):
    e = ent(m)
    print(f'  {nm}: {e.size} entries; top: {dict(e.head(10))}')

S = M[M._merge == 'both'].copy()
dang = lambda a, b: np.abs((a - b + 180) % 360 - 180)
S['dphi'], S['dpsi'], S['dome'] = dang(S.phi_o, S.phi_n), dang(S.psi_o, S.psi_n), dang(S.omega_o, S.omega_n)
S['dB'] = (S.own_mc_max_o - S.own_mc_max_n).abs()
print(f'\nshared residues {len(S):,}')
print(f'  residue name differs: {int((S.resname_o != S.resname_n).sum()):,}; '
      f'examples: {S.loc[S.resname_o != S.resname_n, ["pdb", "resname_o", "resname_n"]].drop_duplicates(["resname_o", "resname_n"]).head(8).values.tolist()}')
print(f'  Ramachandran class differs: {int((S.rama6_o != S.rama6_n).sum()):,}')
for c, nm in (('dphi', 'φ'), ('dpsi', 'ψ'), ('dome', 'ω')):
    d = S[c].dropna()
    print(f'  |Δ{nm}| > 0.01°: {int((d > 0.01).sum()):,} in {S.loc[d[d > 0.01].index, "pdb"].nunique()} entries; > 1°: {int((d > 1).sum()):,}; max {d.max():.2f}')
nd = lambda a, b: (a.isna() != b.isna())
print(f'  φ/ψ/ω defined in one file only: {int((nd(S.phi_o, S.phi_n) | nd(S.psi_o, S.psi_n) | nd(S.omega_o, S.omega_n)).sum()):,}')
print(f'  mainchain max B differs > 0.01: {int((S.dB > 0.01).sum()):,} in {S.loc[S.dB > 0.01, "pdb"].nunique()} entries; max {S.dB.max():.2f}')
print(f'  own heavy-atom min occupancy differs: {int((S.own_minocc_o != S.own_minocc_n).sum()):,}')
print(f'  mainchain altloc flag differs: {int((S.own_alt_o != S.own_alt_n).sum()):,}')
print(f'  peptide link (C–N < 2.0 Å) differs, prev or next: '
      f'{int((((S.d_prevC_N_o < 2) != (S.d_prevC_N_n < 2)) | ((S.d_C_nextN_o < 2) != (S.d_C_nextN_n < 2))).sum()):,}')
print(f'  amide H present differs: {int((S.has_H_o.fillna(0) != S.has_H_n.fillna(0)).sum()):,}')
print(f'  amide H min occupancy differs: {int((S.H_minocc_o.fillna(-1) != S.H_minocc_n.fillna(-1)).sum()):,}')
print(f'  all-atom min occupancy differs: {int((S.all_minocc_o.fillna(-1) != S.all_minocc_n.fillna(-1)).sum()):,}')
x = S[S.all_minocc_o.fillna(-1) != S.all_minocc_n.fillna(-1)]
print('    (2011, current) all-atom min occupancy, top pairs:',
      dict(x.groupby(['all_minocc_o', 'all_minocc_n']).size().sort_values(ascending=False).head(6)))
print(f'    of these, heavy-atom min occupancy also differs: {int((x.own_minocc_o != x.own_minocc_n).sum()):,}')

# entries with any heavy-atom-level difference (hydrogens ignored)
hv = ((S.resname_o != S.resname_n) | (S.dphi > 0.01) | (S.dpsi > 0.01) | (S.dome > 0.01) | (S.dB > 0.01)
      | (S.own_minocc_o != S.own_minocc_n) | (S.own_alt_o != S.own_alt_n)
      | ((S.d_prevC_N_o < 2) != (S.d_prevC_N_n < 2)) | ((S.d_C_nextN_o < 2) != (S.d_C_nextN_n < 2)))
e1 = set(S[hv].pdb) | set(M[M._merge != 'both'].pdb)
big = pd.concat([S[hv].pdb, M[M._merge != 'both'].pdb]).value_counts()
print(f'\nentries with any heavy-atom difference (number, name, φ/ψ/ω, B, occupancy, altloc, link): {len(e1)} of '
      f'{M.pdb.nunique():,}; residues {int(hv.sum()) + int((M._merge != "both").sum()):,}')
print(f'  entries with ≥ 20 differing residues: {dict(big[big >= 20].sort_values(ascending=False))}')
print(f'  entries with < 20: {int((big < 20).sum())}, residues {int(big[big < 20].sum())}')
