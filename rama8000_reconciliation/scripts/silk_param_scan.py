"""Silk parameter scan for one non-Pro class under the reconciled rule (own mc 0 < B < 30, all atoms occ 1 incl. amide H, standard).

usage: silk_param_scan.py CLASS NDFT COSINE[,COSINE...] TWOPASS[,TWOPASS...]
Gly is mirrored explicitly ((phi, psi) and (-phi, -psi)), as Rota2D -sym does; the archived glycine.ndft is exactly
point-symmetric and titled "Gly (sym)". Run from rama8000_reconciliation/.
"""
import sys, time
sys.path.insert(0, 'scripts')
import numpy as np, pandas as pd
src = open('scripts/b30_grid_check.py').read()
exec(src[:src.index('jobs = [')])   # loads d, common, silk()
import ndft
cls, path = sys.argv[1], sys.argv[2]
cosines, twopasses = sys.argv[3].split(','), sys.argv[4].split(',')
rc, T = ndft.load(path)
df = d[common & (d.own_mc_max < 30) & (d.rama6 == cls) & ~((d.rama6 == 'TransPro') & (d.omega.abs() < 150))]

def silk_ct(df, cos, tp, sym):
    pts = list(zip(df.phi, df.psi)) + ([(-p, -q) for p, q in zip(df.phi, df.psi)] if sym else [])
    import subprocess, tempfile, os
    fd, fn = tempfile.mkstemp(suffix='.ndft'); os.close(fd)
    cmd = ['java', '-Xmx8g', '-cp', 'silk_build/classes', 'silk.SilkCmdLine', '-ndim=2', '-coords=2,3', '-bounds=-180,180,-180,180',
           '-wrap', '-gridsize=2', f'-cosine={cos}', f'-twopass={tp}', '-fraction', '-ndft', '-quiet']
    txt = '\n'.join(f'x {p} {q}' for p, q in pts) + '\n'
    open(fn, 'wb').write(subprocess.run(cmd, input=txt.encode(), capture_output=True, check=True).stdout)
    n, g = ndft.load(fn); os.remove(fn); return n, g

print(f'{cls}: realcount {rc}; n points {len(df)}{" (x2 mirrored)" if cls == "Gly" else ""}')
print(f'{"cosine":>7s}{"twopass":>8s}{"meanΔ×1e4":>11s}{"supp":>6s}{"rank-identical":>16s}{"max|Δrank|":>11s}{"sec":>6s}')
for c in cosines:
    for tp in twopasses:
        t0 = time.time(); n, g = silk_ct(df, c, tp, cls == 'Gly'); dr = np.rint(g * rc) - np.rint(T * rc)
        print(f'{c:>7s}{tp:>8s}{np.abs(g-T).mean()*1e4:11.4f}{int(((g>0)!=(T>0)).sum()):6d}'
              f'{int((dr==0).sum()):>10d}/{dr.size}{int(np.abs(dr).max()):11d}{time.time()-t0:6.0f}', flush=True)
