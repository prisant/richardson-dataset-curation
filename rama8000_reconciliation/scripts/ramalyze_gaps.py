"""Why does ramalyze omit residues that our records have phi/psi for (and that pass the peptide-link test)?

Profiles the residues in our records but not in ramalyze's output (ersatz files), restricted to C-N < 2.0 Å on both
sides, against neighbour properties: previous / next residue non-standard (not in the 20 + MSE), neighbour or own
altloc, missing atoms (sidechain_records.tsv), whether the residue passes the best rule, by class.
usage: ramalyze_gaps.py [FH|ersatz]   Run from rama8000_reconciliation/.
"""
import sys
sys.path.insert(0, 'scripts')
import numpy as np, pandas as pd
from recount3 import chains
kind = sys.argv[1] if len(sys.argv) > 1 else 'ersatz'
rd = lambda f, **k: pd.read_csv(f, sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str}, **k)
d = rd('data/residue_records.tsv').merge(rd('data/amideH_records.tsv'), on=['pdb', 'chain', 'resnum', 'ins'], how='left')
s = rd('data/sidechain_records.tsv', usecols=['pdb', 'chain', 'resnum', 'ins', 'n_missing', 'missing'])
d = d.merge(s, on=['pdb', 'chain', 'resnum', 'ins'], how='left'); d['ins'] = d.ins.str.strip()
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
OK = STD | {'MSE'}
key = ['pdb', 'chain', 'resnum', 'ins']
# previous residue name and altloc from adjacent rows of the full sidechain table (all residues, chain order)
allr = s.copy(); allr['ins'] = allr.ins.str.strip()
full = rd('data/sidechain_records.tsv', usecols=['pdb', 'chain', 'resnum', 'ins', 'resname']); full['ins'] = full.ins.str.strip()
same = lambda k: (full.pdb.shift(k) == full.pdb) & (full.chain.shift(k) == full.chain)
full['prev_rn'] = full.resname.shift(1).where(same(1)); full['next_rn2'] = full.resname.shift(-1).where(same(-1))
d = d.merge(full[key + ['prev_rn', 'next_rn2']], on=key, how='left')
r = pd.read_csv(f'data/ramalyze_{kind}.tsv', sep='\t', keep_default_na=False, dtype=str)
r = r[(r.pdb != '#ERR') & r.model.isin(['', '1'])].rename(columns={'resseq': 'resnum', 'icode': 'ins'})
r['resnum'] = r.resnum.astype(int); r = r[key].drop_duplicates()
m = d.merge(r, on=key, how='left', indicator=True)
link = (m.d_prevC_N < 2.0) & (m.d_C_nextN < 2.0)
gap = m[link & (m._merge == 'left_only')].copy()
aw = gap.omega.abs()
best = ((gap.own_mc_max < 30) & (gap.own_mc_max > 0) & (gap.own_minocc >= 1) & ~(gap.H_minocc < 1) & gap.resname.isin(OK)
        & np.where(gap.rama6 == 'CisPro', aw < 30, aw.isna() | (aw >= 150)))
print(f'{kind}: residues with link both sides but not reported by ramalyze: {len(gap)}; passing the best rule: {int(best.sum())}')
feats = {'own non-standard (not 20+MSE)': ~gap.resname.isin(OK),
         'prev non-standard': ~gap.prev_rn.isin(OK), 'next non-standard': ~gap.next_resname.isin(OK),
         'own altloc': gap.own_alt > 0, 'own atom missing': gap.n_missing > 0,
         'own O missing': gap.missing.fillna('').str.split(',').apply(lambda x: 'O' in x)}
for k, v in feats.items():
    print(f'  {k:32s} all {int(v.sum()):6d}   best-rule {int((v & best).sum()):6d}')
anyf = np.logical_or.reduce(list(feats.values()))
print(f'  {"none of the above":32s} all {int((~anyf).sum()):6d}   best-rule {int((~anyf & best).sum()):6d}')
print('\nbest-rule residues omitted by ramalyze, by class:', gap[best].rama6.value_counts().to_dict())
print('next residue names (best-rule, next non-standard):', gap[best & ~gap.next_resname.isin(OK)].next_resname.value_counts().head(10).to_dict())
print('prev residue names (best-rule, prev non-standard):', gap[best & ~gap.prev_rn.isin(OK)].prev_rn.value_counts().head(10).to_dict())
rest = gap[best & ~anyf]
if len(rest):
    print('\nunexplained examples:'); print(rest[key + ['resname', 'prev_rn', 'next_resname', 'd_prevC_N', 'd_C_nextN', 'own_alt']].head(15).to_string(index=False))
