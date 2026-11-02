"""The 2011 Java tools on the lab's 2011 chain files vs the lab lists, for the residues our rule gets wrong.

Run with T8000_DATA=data/lab2011 (tables from the lab's files, see lab2011_tree.py). Rule: 'H rule + MSE B + ω defined'
from lab2013_compare.py (on these files identical to 'every atom in the file at occupancy 1').
Dangle: data/lab2011/dangle2011.tsv (dangle2011_run.sh; default alt, one row per residue).
Ramalyze: data/lab2011/ramalyze2011.tsv (ramalyze2011_run.sh; one row per backbone alternate), if present.
For ours-only and lab-only residues: is the residue in Dangle's output, with φ, ψ, ω all defined, and do Dangle's
angles agree with ours (0.01°)? Same for Ramalyze rows (count, φ/ψ).
Run from rama8000_reconciliation/.
"""
import sys, os, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    exec(open('scripts/lab2013_compare.py').read().split('# ---- Evidence for the refined rules')[0])
RULE = 'H rule + MSE B + ω defined'
R = V[RULE]
num = lambda s: pd.to_numeric(s.replace('__?__', np.nan), errors='coerce')
dg = pd.read_csv(f'{DATA}/dangle2011.tsv', sep=':', header=None, dtype=str,
                 names=['file', 'model', 'chain', 'resnum', 'ins', 'resname', 'phi', 'psi', 'omega'])
dg['pdb'] = dg.file.str[:4].str.lower(); dg['resnum'] = dg.resnum.str.strip().astype(int); dg['ins'] = dg.ins.str.strip()
dg = dg[dg.model.str.strip() == '1']
for c in ('phi', 'psi', 'omega'): dg['dg_' + c] = num(dg[c])
dg = dg.drop_duplicates(K)[K + ['dg_phi', 'dg_psi', 'dg_omega']]
X = E.merge(dg, on=K, how='left')
lab = X._merge != 'left_only'
dang = lambda a, b: np.abs((a - b + 180) % 360 - 180)
X['in_dangle'] = X.dg_phi.notna() | X.dg_psi.notna() | X.dg_omega.notna()
X['dg_all'] = X.dg_phi.notna() & X.dg_psi.notna() & X.dg_omega.notna()
X['dg_same'] = (dang(X.dg_phi, X.phi) <= 0.01) & (dang(X.dg_psi, X.psi) <= 0.01)
print(f'rule: {RULE}; Dangle rows (model 1) {len(dg):,}')
for nm, m in (('ours-only', R & ~lab), ('lab-only', lab & ~R), ('shared', R & lab)):
    Y = X[m]
    print(f'\n{nm}: {len(Y):,}; in Dangle {int(Y.in_dangle.sum()):,}; φ, ψ, ω all defined {int(Y.dg_all.sum()):,};'
          f' Dangle φ/ψ = ours {int(Y.dg_same.sum()):,}')
    if nm != 'shared':
        print('  by class:', Y.rama6.fillna(Y.lab_cls).value_counts().to_dict())
        print('  ω undefined in Dangle:', int((Y.in_dangle & Y.dg_omega.isna()).sum()), '; φ/ψ differ from ours:', int((Y.dg_all & ~Y.dg_same).sum()))
rf = f'{DATA}/ramalyze2011.tsv'
if os.path.exists(rf):
    rows = []
    for ln in open(rf):
        f, rest = ln.rstrip('\n').split(':', 1)
        lab_, *v = rest.split(':')
        rn, rt = lab_[:-4].strip(), lab_[-3:]
        rows.append((f[:4].lower(), f.split('_')[-1], int(rn.rstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ')), rn.lstrip('-0123456789'), rt,
                     float(v[1]) if v[1] else np.nan, float(v[2]) if v[2] else np.nan))
    rz = pd.DataFrame(rows, columns=K + ['rz_resname', 'rz_phi', 'rz_psi'])
    rc = rz.groupby(K).agg(rz_n=('rz_phi', 'size'), rz_phi=('rz_phi', 'first'), rz_psi=('rz_psi', 'first')).reset_index()
    Z = X.merge(rc, on=K, how='left')
    print(f'\nRamalyze rows {len(rz):,}, residues {len(rc):,}')
    for nm, m in (('ours-only', R & ~lab), ('lab-only', lab & ~R), ('shared', R & lab)):
        Y = Z[m]
        print(f'{nm}: {len(Y):,}; in Ramalyze {int(Y.rz_n.notna().sum()):,}; rows per residue',
              Y.rz_n.fillna(0).astype(int).value_counts().sort_index().to_dict())
