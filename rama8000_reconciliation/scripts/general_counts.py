"""Class counts from residue_records.tsv under the reconciled rule and General-specific variants."""
import pandas as pd, numpy as np
d = pd.read_csv('data/residue_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
link = (d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0)
rule = link & (d.own_mc_max < 30) & (d.own_minocc >= 1) & (d.own_mc_max > 0)
aom = d.omega.abs()
print('sanity CisPro  :', int((rule & (d.rama6 == 'CisPro')).sum()), '(NDFT 3654)')
print('sanity TransPro:', int((rule & (d.rama6 == 'TransPro') & (aom >= 150)).sum()), '(NDFT 66364; expect 66363)')
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
std = d.resname.isin(STD)
gen = rule & (d.rama6 == 'General')
nd = 1061639
v = {
 'rule, standard residues only':                 gen & std,
 'rule, + MSE':                                  gen & (std | (d.resname == 'MSE')),
 'rule, + all non-standard':                     gen,
 'std, drop preceding cis/twisted (|ω|<150)':    gen & std & (aom >= 150),
 'std+MSE, drop preceding |ω|<150':              gen & (std | (d.resname == 'MSE')) & (aom >= 150),
}
for k, m in v.items(): print(f'General {k:44s} n={int(m.sum()):>9,d}  Δ={int(m.sum())-nd:+,d}')
nonstd = d[gen & ~std].resname.value_counts()
print('\nnon-standard General-class residues passing rule:', dict(nonstd.head(12)))
