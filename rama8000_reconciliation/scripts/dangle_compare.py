"""Java Dangle 1.07 vs our parser: residue coverage, phi/psi/omega agreement, and category counts under the best rule.

Input: data/dangle_FH.tsv and data/dangle_ersatz.tsv (dangle_run.sh), data/residue_records.tsv + amideH_records.tsv.
Dangle rows restricted to model 1 and the listed top8000 chains. For each Dangle input:
  1. coverage: residues with numeric phi and psi in Dangle vs in our records; duplicate Dangle keys;
  2. angle agreement on the shared residues (|Δ| > 0.01°, wrapped);
  3. category counts: best rule (own attributes from our records) with class, trans test and phi/psi presence taken
     from Dangle, vs ours and the NDFT realcount.
Best rule: peptide-linked both sides (C-N < 2.0 Å); own mainchain 0 < B < 30; all heavy atoms occ 1; amide H occ 1;
standard residues + MSE; preceding peptide trans (|ω| ≥ 150), except CisPro (|ω| < 30).
usage: dangle_compare.py [FH] [ersatz]   Run from rama8000_reconciliation/
"""
import sys
sys.path.insert(0, 'scripts')
import numpy as np, pandas as pd
from recount3 import chains

d = pd.read_csv('data/residue_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
h = pd.read_csv('data/amideH_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
d = d.merge(h, on=['pdb', 'chain', 'resnum', 'ins'], how='left')
d['ins'] = d.ins.str.strip()
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
nd = dict(General=1061639, Gly=116589, IleVal=195083, PrePro=61568, TransPro=66364, CisPro=3654)
listed = {(p, c) for p, cs in chains.items() for c in cs}

def rama6(rn, nrn, om):
    out = np.full(len(rn), 'General', dtype=object)
    out[rn.isin(['ILE', 'VAL']).to_numpy()] = 'IleVal'
    out[(nrn == 'PRO').to_numpy()] = 'PrePro'
    out[(rn == 'PRO').to_numpy()] = np.where(om[(rn == 'PRO').to_numpy()].abs() < 30, 'CisPro', 'TransPro')
    out[(rn == 'GLY').to_numpy()] = 'Gly'
    return out

def rule(df, om, cls, res):
    aw = om.abs()
    return ((df.d_prevC_N < 2.0) & (df.d_C_nextN < 2.0) & (df.own_mc_max < 30) & (df.own_mc_max > 0) & (df.own_minocc >= 1)
            & ~(df.H_minocc < 1) & df.resname.isin(res)
            & pd.Series(np.where(cls == 'CisPro', aw < 30, aw.isna() | (aw >= 150)), index=df.index))

def counts(m, cls):
    return {c: int((m & (cls == c)).sum()) for c in nd}

print('ours (best rule):')
base_cls = d.rama6.to_numpy()
ours = counts(rule(d, d.omega, base_cls, STD | {'MSE'}), base_cls)
print('  ' + '  '.join(f'{c} {ours[c] - nd[c]:+d}' for c in nd))

for kind in (sys.argv[1:] or ['FH', 'ersatz']):
    g = pd.read_csv(f'data/dangle_{kind}.tsv', sep=':', header=None, dtype=str,
                    names=['label', 'model', 'chain', 'resnum', 'ins', 'resname', 'phi', 'psi', 'omega'])
    g = g[g.model.str.strip() == '1'].copy()
    g['pdb'] = g.label.str.replace(r'(FH|_ersatz)\.pdb$', '', regex=True)
    g['chain'] = g.chain.str.strip(); g['ins'] = g.ins.str.strip(); g['resnum'] = g.resnum.str.strip().astype(int)
    g = g[[k in listed for k in zip(g.pdb, g.chain)]]
    for a in ('phi', 'psi', 'omega'): g[a] = pd.to_numeric(g[a], errors='coerce')
    key = ['pdb', 'chain', 'resnum', 'ins']
    dup = g.duplicated(key, keep=False)
    gp = g[g.phi.notna() & g.psi.notna()].drop_duplicates(key)
    print(f'\n=== Dangle on {kind} files: {g.pdb.nunique()} structures, {len(g)} rows in listed chains, '
          f'{int(dup.sum())} rows with duplicate keys, {len(gp)} with phi and psi')
    m = d.merge(gp[key + ['phi', 'psi', 'omega', 'resname']], on=key, how='outer', suffixes=('', '_dg'), indicator=True)
    print('coverage (phi and psi computable):', m._merge.map({'both': 'both', 'left_only': 'ours only', 'right_only': 'Dangle only'}).value_counts().to_dict())
    on = m[m._merge == 'left_only']; dn = m[m._merge == 'right_only']
    if len(on): print('  ours-only examples:\n' + on[key + ['resname', 'd_prevC_N', 'd_C_nextN']].head(8).to_string(index=False))
    if len(dn): print('  Dangle-only by resname:', dn.resname_dg.value_counts().head(10).to_dict())
    b = m[m._merge == 'both'].copy()
    wrap = lambda x: (x + 180) % 360 - 180
    for a in ('phi', 'psi', 'omega'):
        dd = wrap(b[f'{a}_dg'] - b[a]).abs()
        print(f'  {a}: |Δ| > 0.01° in {int((dd > 0.01).sum())} of {int(dd.notna().sum())}; max {dd.max():.3f}°; '
              f'missing in Dangle {int(b[f"{a}_dg"].isna().sum())}, in ours {int(b[a].isna().sum())}')
    bb = d.merge(gp[key + ['omega']], on=key, how='inner', suffixes=('', '_dg'))
    cls_dg = rama6(bb.resname, bb.next_resname, bb.omega_dg)
    print(f'  class changes (our ω vs Dangle ω): {int((cls_dg != bb.rama6.to_numpy()).sum())}')
    cd = counts(rule(bb, bb.omega_dg, cls_dg, STD | {'MSE'}), cls_dg)
    print(f'{"class":10s}{"NDFT":>10s}{"ours":>10s}{"Dangle":>10s}{"ours-NDFT":>11s}{"Dangle-NDFT":>13s}{"Dangle-ours":>13s}')
    for c in nd:
        print(f'{c:10s}{nd[c]:10d}{ours[c]:10d}{cd[c]:10d}{ours[c] - nd[c]:+11d}{cd[c] - nd[c]:+13d}{cd[c] - ours[c]:+13d}')
