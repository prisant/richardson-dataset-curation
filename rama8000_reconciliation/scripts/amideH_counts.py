"""Class counts under the reconciled rule plus amide-H exclusion variants."""
import pandas as pd
d = pd.read_csv('data/residue_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
h = pd.read_csv('data/amideH_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
d = d.merge(h, on=['pdb', 'chain', 'resnum', 'ins'], how='left')
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
rule = (d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0) & (d.own_mc_max < 30) & (d.own_minocc >= 1) & (d.own_mc_max > 0) & d.resname.isin(STD)
rule &= ~((d.rama6 == 'TransPro') & (d.omega.abs() < 150))
nd = dict(General=1061639, Gly=116589, IleVal=195083, PrePro=61568, TransPro=66364, CisPro=3654)
V = {'rule': rule, 'drop no amide H': rule & ~((d.rama6.isin(['TransPro','CisPro']) == False) & (d.has_H != 1)),
     'drop H B>=30': rule & ~(d.H_maxB >= 30), 'drop H occ<1': rule & ~(d.H_minocc < 1)}
V['drop no-H, H B>=30, H occ<1'] = V['drop no amide H'] & V['drop H B>=30'] & V['drop H occ<1']
print(f'{"variant":30s}' + ''.join(f'{c:>11s}' for c in nd))
for k, m in V.items():
    print(f'{k:30s}' + ''.join(f'{int((m & (d.rama6 == c)).sum()) - n:>+11,d}' for c, n in nd.items()))
