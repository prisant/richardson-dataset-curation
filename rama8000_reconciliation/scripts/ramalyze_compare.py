"""cctbx ramalyze vs our parser: residue coverage, class assignment, phi/psi, and category counts.

Input: data/ramalyze_{FH,ersatz}.tsv (ramalyze_run.py, cctbx-base 2026.9), our records (residue_records.tsv +
amideH_records.tsv). Rows restricted to model 1 (or blank) and the listed top8000 chains; one residue = one key
(chain, resseq, icode); ramalyze writes one row per altloc of the three-residue window.
For each input:
  1. coverage: residues ramalyze reports vs residues with phi/psi in our records; what each side lacks;
  2. class: ramalyze res_type vs our rama6 (ramalyze splits cis/trans Pro at |ω| = 90);
  3. phi/psi agreement on shared residues (first row: altloc blank or A);
  4. category counts:
     a. best rule, restricted to residues ramalyze reports;
     b. as (a) with the residue-type test (standard + MSE) replaced by ramalyze coverage;
     c. as (a) with the amide-H rule replaced by 'ramalyze reports no altloc for the residue';
     d. as (b) and (c) together (cctbx-native residue set and altloc handling).
Best rule: peptide-linked both sides (C-N < 2.0 Å); own mainchain 0 < B < 30; all heavy atoms occ 1; amide H occ 1;
standard residues + MSE; preceding peptide trans (|ω| ≥ 150), except CisPro (|ω| < 30).
usage: ramalyze_compare.py [FH] [ersatz]   Run from rama8000_reconciliation/.
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
RT = {0: 'General', 1: 'Gly', 2: 'CisPro', 3: 'TransPro', 4: 'PrePro', 5: 'IleVal'}
listed = {(p, c) for p, cs in chains.items() for c in cs}
key = ['pdb', 'chain', 'resnum', 'ins']
aw = d.omega.abs()
R = dict(link=(d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0), B=(d.own_mc_max < 30) & (d.own_mc_max > 0), occ=d.own_minocc >= 1,
         H=~(d.H_minocc < 1), res=d.resname.isin(STD | {'MSE'}),
         trans=pd.Series(np.where(d.rama6 == 'CisPro', aw < 30, aw.isna() | (aw >= 150)), index=d.index))
best = np.logical_and.reduce(list(R.values()))
cnt = lambda m: {c: int((m & (d.rama6 == c)).sum()) for c in nd}
ours = cnt(best)

for kind in (sys.argv[1:] or ['FH', 'ersatz']):
    r = pd.read_csv(f'data/ramalyze_{kind}.tsv', sep='\t', keep_default_na=False, dtype=str, comment=None)
    err = r[r.pdb == '#ERR']; r = r[r.pdb != '#ERR']
    r = r[r.model.isin(['', '1'])].rename(columns={'resseq': 'resnum', 'icode': 'ins'})
    r['resnum'] = r.resnum.astype(int)
    r = r[[k in listed for k in zip(r.pdb, r.chain)]]
    r['phi'] = r.phi.astype(float); r['psi'] = r.psi.astype(float); r['res_type'] = r.res_type.astype(int)
    r = r.sort_values(key + ['altloc'])
    g = r.groupby(key, sort=False).agg(n_rows=('altloc', 'size'), alts=('altloc', lambda x: ''.join(sorted(set(x)))),
                                       resname_r=('resname', 'first'), phi_r=('phi', 'first'), psi_r=('psi', 'first'),
                                       rt=('res_type', 'first')).reset_index()
    g['has_alt'] = g.alts != ''
    print(f'\n=== ramalyze on {kind} files: {r.pdb.nunique()} structures, {len(err)} errors, {len(r)} rows in listed chains, '
          f'{len(g)} residues ({int(g.has_alt.sum())} with altloc rows)')
    if len(err): print(err.head(10).to_string(index=False))
    m = d.merge(g, on=key, how='outer', indicator=True)
    lab = m._merge.map({'both': 'both', 'left_only': 'ours only', 'right_only': 'ramalyze only'})
    print('coverage:', lab.value_counts().to_dict())
    on = m[m._merge == 'left_only']; rn = m[m._merge == 'right_only']
    erp = set(err.model)
    print(f'  ours only: in structures ramalyze failed to read {int(on.pdb.isin(erp).sum())}')
    on = on[~on.pdb.isin(erp)]
    print('  ours only (readable structures), by resname:', on.resname.value_counts().head(12).to_dict())
    ns = d[~d.resname.isin(STD) & ~d.pdb.isin(erp)].merge(g[key], on=key, how='left', indicator=True)
    print('  non-standard residues with phi/psi in ours -> reported by ramalyze:',
          ns.groupby('resname')._merge.apply(lambda x: f'{int((x == "both").sum())}/{len(x)}').sort_values(key=lambda x: -x.str.split('/').str[1].astype(int)).head(10).to_dict())
    print(f'  ours only: link both < 2.0 Å {int(((on.d_prevC_N < 2) & (on.d_C_nextN < 2)).sum())}; '
          f'passing the best rule {int(best[d.merge(on[key], on=key, how="left", indicator=True)._merge.eq("both").to_numpy()].sum())}')
    print('  ramalyze only, by resname:', rn.resname_r.value_counts().head(12).to_dict())
    b = m[m._merge == 'both']
    cr = b.rt.map(RT)
    print(f'  class differences (ramalyze vs ours): {int((cr != b.rama6).sum())}; ',
          pd.crosstab(b.rama6[cr != b.rama6], cr[cr != b.rama6]).stack().loc[lambda x: x > 0].to_dict())
    wrap = lambda x: ((x + 180) % 360 - 180).abs()
    for a in ('phi', 'psi'):
        dd = wrap(b[f'{a}_r'] - b[a])
        print(f'  {a}: |Δ| > 0.01° in {int((dd > 0.01).sum())} of {len(dd)}; max {dd.max():.3f}°')

    inr = d.merge(g[key + ['has_alt']], on=key, how='left', indicator=True)
    cov = pd.Series(inr._merge.eq('both').to_numpy(), index=d.index)
    noalt = pd.Series(~inr.has_alt.fillna(False).astype(bool).to_numpy(), index=d.index)
    rest = lambda *drop: np.logical_and.reduce([v for k, v in R.items() if k not in drop])
    V = {'ours (best rule)': best,
         'a. best rule ∩ ramalyze coverage': best & cov,
         'b. a, residue test → ramalyze coverage': rest('res') & cov,
         'c. a, H rule → no ramalyze altloc': rest('H') & cov & noalt,
         'd. b + c (cctbx-native)': rest('res', 'H') & cov & noalt}
    print(f'\n{"variant":40s}' + ''.join(f'{c:>10s}' for c in nd) + f'{"total":>10s}')
    for name, mm in V.items():
        dl = [v - nd[c] for c, v in cnt(mm).items()]
        print(f'{name:40s}' + ''.join(f'{x:>+10,d}' for x in dl) + f'{sum(dl):>+10,d}')
