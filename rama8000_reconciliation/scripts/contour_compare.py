"""Compare contour polylines of our grids (best rule) with the lab's NDFT grids, per class.

NDFT values are Silk '-fraction': the fraction of input points at lower density, so the X% contour is level 1 - X.
Inner contour 98% (0.02) for all classes; outer 99.95% (0.0005) General, 99.90% (0.001) IleVal, Gly, TransPro,
PrePro, 99.85% (0.0015) CisPro.
Best rule: peptide-linked both sides (C-N < 2.0 Å); own mainchain 0 < B < 30; all heavy atoms occ 1; amide H occ 1;
standard residues + MSE; preceding peptide trans (|ω| ≥ 150), except CisPro (|ω| < 30). Silk cosine 5, twopass 8
(General 6.5); Gly mirrored.
Per class and level:
  regions   bins with value ≥ level in each grid, bins inside only one, symmetric difference as area (deg²),
            connected regions (periodic);
  polylines contourpy lines on the periodically padded grid (bin centres -179..179), clipped to [-180, 180]²;
            distance from each vertex (polylines resampled at 0.25°) to the nearest point of the other set,
            periodic in phi and psi: mean, 95th percentile, max (Hausdorff = max of both directions).
Writes results/contour_compare.txt (stdout), results/contour_compare.png (overview) and one figure per class and
contour, results/contours/<class>_<contour>.png: lines, bins inside only one grid, point of largest deviation.
Run from rama8000_reconciliation/.
"""
import sys
sys.path.insert(0, 'scripts')
import numpy as np, pandas as pd
sys.argv = ['x', 'Gly', 'data/glycine.ndft', '5', '8']
src = open('scripts/silk_param_scan.py').read()
exec(src[:src.index("print(f'{cls}: realcount")])   # d (with amide-H cols), silk_ct
import ndft, contourpy
from scipy import ndimage
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
aw = d.omega.abs()
rule = ((d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0) & (d.own_mc_max < 30) & (d.own_mc_max > 0) & (d.own_minocc >= 1)
        & ~(d.H_minocc < 1) & d.resname.isin(STD | {'MSE'})
        & pd.Series(np.where(d.rama6 == 'CisPro', aw < 30, aw.isna() | (aw >= 150)), index=d.index))
setup = [('General', 'general', '6.5', 0.0005), ('Gly', 'glycine', '8', 0.001), ('IleVal', 'ileval', '8', 0.001),
         ('PrePro', 'prepro', '8', 0.001), ('TransPro', 'transpro', '8', 0.001), ('CisPro', 'cispro', '8', 0.0015)]
P = 4
ax_c = -179 + 2 * np.arange(-P, 180 + P)

def lines(G, lev):
    # G[i, j]: i = phi bin, j = psi bin -> z[y=psi, x=phi]
    z = np.pad(G, P, mode='wrap').T
    out = []
    for L in contourpy.contour_generator(ax_c, ax_c, z, line_type='Separate').lines(lev):
        seg = []
        for p in L:                              # resample at 0.25 deg
            if seg:
                q = seg[-1]; k = max(1, int(np.ceil(np.hypot(*(p - q)) / 0.25)))
                seg += [q + (p - q) * t / k for t in range(1, k + 1)]
            else: seg.append(p)
        seg = np.array(seg); seg = seg[(np.abs(seg) <= 180).all(1)]
        if len(seg): out.append(seg)
    return np.vstack(out) if out else np.empty((0, 2)), out

def regions(mask):
    lab, n = ndimage.label(mask)
    par = list(range(n + 1))
    def f(x):
        while par[x] != x: x = par[x]
        return x
    for a, b in list(zip(lab[0], lab[-1])) + list(zip(lab[:, 0], lab[:, -1])):
        if a and b: par[f(a)] = f(b)
    return len({f(x) for x in range(1, n + 1)})

def dist(A, B):
    if not len(A) or not len(B): return np.array([np.nan])
    return cKDTree((B + 180) % 360, boxsize=360).query((A + 180) % 360)[0]

import os; os.makedirs('results/contours', exist_ok=True)
fig, axs = plt.subplots(2, 3, figsize=(15, 10))
print(f'{"class":9s}{"contour":>9s}{"level":>8s}{"bins NDFT":>10s}{"bins ours":>10s}{"only NDFT":>10s}{"only ours":>10s}{"symdiff deg²":>13s}'
      f'{"regions N/O":>12s}{"mean°":>7s}{"p95°":>7s}{"max° O→N":>9s}{"max° N→O":>9s}{"Hausdorff°":>11s}')
for ax, (cls, f, tp, outer) in zip(axs.flat, setup):
    rc, T = ndft.load(f'data/{f}.ndft'); sym = cls == 'Gly'
    n, g = silk_ct(d[rule & (d.rama6 == cls)], '5', tp, sym)
    for name, lev, ls in (('98%', 0.02, '-'), (f'{100 * (1 - outer):.2f}%', outer, '-')):
        mT, mg = T >= lev, g >= lev
        vT, sT = lines(T, lev); vg, sg = lines(g, lev)
        dOT, dTO = dist(vg, vT), dist(vT, vg); both = np.concatenate([dOT, dTO])
        print(f'{cls:9s}{name:>9s}{lev:8.4f}{int(mT.sum()):10d}{int(mg.sum()):10d}{int((mT & ~mg).sum()):10d}{int((mg & ~mT).sum()):10d}'
              f'{4 * int((mT ^ mg).sum()):13d}{f"{regions(mT)}/{regions(mg)}":>12s}{np.mean(both):7.2f}{np.percentile(both, 95):7.2f}'
              f'{np.max(dOT):9.2f}{np.max(dTO):9.2f}{max(np.max(dOT), np.max(dTO)):11.2f}', flush=True)
        for s in sT: ax.plot(s[:, 0], s[:, 1], 'k-', lw=0.9)
        for s in sg: ax.plot(s[:, 0], s[:, 1], 'r--', lw=0.9)
        f1, a1 = plt.subplots(figsize=(8.5, 8.5))
        for m, c, lab in ((mT & ~mg, 'tab:blue', f'inside NDFT only ({int((mT & ~mg).sum())} bins)'),
                          (mg & ~mT, 'tab:orange', f'inside ours only ({int((mg & ~mT).sum())} bins)')):
            for i, j in zip(*np.nonzero(m)):
                a1.add_patch(plt.Rectangle((-180 + 2 * i, -180 + 2 * j), 2, 2, color=c, alpha=0.8, lw=0))
            a1.plot([], [], 's', color=c, label=lab)
        for k, s in enumerate(sT): a1.plot(s[:, 0], s[:, 1], 'k-', lw=1.2, label='NDFT (lab)' if k == 0 else None)
        for k, s in enumerate(sg): a1.plot(s[:, 0], s[:, 1], 'r--', lw=1.2, label='ours (best rule)' if k == 0 else None)
        src_, dd = (vg, dOT) if np.max(dOT) >= np.max(dTO) else (vT, dTO)
        q = src_[np.argmax(dd)]
        a1.plot(*q, 'o', ms=14, mfc='none', mec='m', mew=1.5, label=f'max deviation {np.max(dd):.2f}° at ({q[0]:.0f}, {q[1]:.0f})')
        a1.set(xlim=(-180, 180), ylim=(-180, 180), xticks=range(-180, 181, 30), yticks=range(-180, 181, 30), aspect='equal',
               xlabel='φ (°)', ylabel='ψ (°)',
               title=f'{cls}, {name} contour (level {lev:g})\nmean deviation {np.mean(both):.2f}°, Hausdorff {max(np.max(dOT), np.max(dTO)):.2f}°')
        a1.grid(alpha=0.25); a1.legend(loc='lower right', fontsize=8, framealpha=0.9)
        f1.tight_layout(); f1.savefig(f'results/contours/{cls}_{name.rstrip("%").replace(".", "p")}.png', dpi=130); plt.close(f1)
    ax.set(xlim=(-180, 180), ylim=(-180, 180), xticks=range(-180, 181, 60), yticks=range(-180, 181, 60),
           title=f'{cls}: 98% and {100 * (1 - outer):.2f}%', xlabel='φ', ylabel='ψ', aspect='equal')
axs.flat[0].plot([], [], 'k-', label='NDFT (lab)'); axs.flat[0].plot([], [], 'r--', label='ours (best rule)'); axs.flat[0].legend(loc='lower right', fontsize=8)
fig.tight_layout(); fig.savefig('results/contour_compare.png', dpi=150)
