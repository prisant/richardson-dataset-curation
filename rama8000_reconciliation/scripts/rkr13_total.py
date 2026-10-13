"""Which residue rule reproduces Richardson, Keedy & Richardson (2013) "1,573,349 total residues for analysis"?

RKR13 sect. 2.1: "we omit any residue that has a backbone atom with B-factor >=30 or an alternate conformation."
Counts over residue_records.tsv (residues with phi and psi in the 7,957 Top8000-best_hom70 chains).
Run from rama8000_reconciliation/.
"""
import pandas as pd
d = pd.read_csv('data/residue_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
h = pd.read_csv('data/amideH_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
d = d.merge(h, on=['pdb', 'chain', 'resnum', 'ins'], how='left')
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
TARGET = 1573349
link = (d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0)
b = d.own_mc_max < 30
noalt = d.own_alt == 0
occ = d.own_minocc >= 1
std = d.resname.isin(STD)
nz = d.own_mc_max > 0
noH = ~(d.H_minocc < 1)
notw = ~((d.rama6 == 'TransPro') & (d.omega.abs() < 150))
V = {
  'all rows':                                 pd.Series(True, index=d.index),
  'B<30':                                     b,
  'B<30, no altloc':                          b & noalt,
  'B<30, no altloc, linked':                  b & noalt & link,
  'B<=30, no altloc':                         (d.own_mc_max <= 30) & noalt,
  'B<30, occ 1':                              b & occ,
  'B<30, occ 1, linked':                      b & occ & link,
  'B<30, occ 1, linked, std':                 b & occ & link & std,
  'B<30, occ 1, linked, std, B>0':            b & occ & link & std & nz,
  'B<30, occ 1, linked, std, B>0, H occ 1':   b & occ & link & std & nz & noH,
  '  + no twisted Pro (= reconciled rule)':   b & occ & link & std & nz & noH & notw,
}
for k, m in V.items():
    n = int(m.sum()); print(f'{k:42s} n={n:>10,d}  Δ vs 1,573,349 = {n - TARGET:>+9,d}')
