"""Support-mismatch diagnosis: where our Silk grid and the archived NDFT disagree on zero vs nonzero density.

usage: support_diag.py CLASS NDFT TWOPASS   (cosine 5; reconciled rule; Gly mirrored)
Reports bins nonzero only in ours / only in the NDFT, their density magnitudes, and their distance from the common support.
Run from rama8000_reconciliation/.
"""
import sys
sys.path.insert(0, 'scripts')
import numpy as np
from scipy import ndimage
src = open('scripts/silk_param_scan.py').read()
cls, path, tp = sys.argv[1], sys.argv[2], sys.argv[3]
sys.argv = ['x', cls, path, '5', tp]
exec(src[:src.index("print(f'{cls}: realcount")])   # loads d, df, rc, T, silk_ct
n, g = silk_ct(df, '5', tp, cls == 'Gly')
ours, theirs = g > 0, T > 0
both = ours & theirs
dist = ndimage.distance_transform_edt(~both)          # bins from the common support (no wrap; fine away from edges)
for name, m, vals in (('ours only', ours & ~theirs, g), ('NDFT only', theirs & ~ours, T)):
    k = int(m.sum())
    if k == 0: print(f'{name:10s} 0 bins'); continue
    v = vals[m] * rc
    print(f'{name:10s} {k:5d} bins; density×realcount min {v.min():.3g} median {np.median(v):.3g} max {v.max():.3g}; '
          f'distance from common support (bins) median {np.median(dist[m]):.1f} max {dist[m].max():.1f}')
print(f'common support {int(both.sum())} bins; ours {int(ours.sum())}, NDFT {int(theirs.sum())}')
print(f'smallest nonzero value×realcount: ours {(g[ours].min()*rc):.3g}, NDFT {(T[theirs].min()*rc):.3g}')

# residues whose (phi, psi) bin, or mirrored bin for Gly, lies in "ours only" support
import ndft as _n
only = ours & ~theirs
bi = np.array([_n.binof(x) for x in df.phi]); bj = np.array([_n.binof(x) for x in df.psi])
hit = only[bi, bj]
if cls == 'Gly': hit |= only[[_n.binof(-x) for x in df.phi], [_n.binof(-x) for x in df.psi]]
out = df[hit]
print(f'\nresidues in ours-only bins: {len(out)}')
if len(out):
    aom = out.omega.abs()
    print('preceding |omega|: <30', int((aom < 30).sum()), ' 30-150', int(((aom >= 30) & (aom < 150)).sum()), ' >=150', int((aom >= 150).sum()))
    print('own mc max B quartiles', np.percentile(out.own_mc_max, [0, 25, 50, 75, 100]).round(1).tolist())
    print(out[['pdb', 'chain', 'resnum', 'ins', 'resname', 'phi', 'psi', 'omega', 'own_mc_max']].to_string(index=False))
