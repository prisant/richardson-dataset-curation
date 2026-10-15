"""Contours drawn the lab's way: Silk -fraction grid piped to kin2Dcont (rlabduke/kincontours), per class.

Lab procedure (rlabduke/reference_data Top8000/Top8000_filtered_residues_rotarama/Makefile, 2018-01-22):
  Silk  -ndim=2 -coords=2,3 -wrap -gridsize=2 -cosine=5 -twopass=8 (General 6.5) -bounds=-181,179,-181,179 -fraction
  kin2Dcont -samp -last -nosmooth -group -wrap -180 180 -180 180 -g2 -l <levels>
i.e. the contour grid has samples on even degrees (-180..178); the NDFT tables (-bounds=-180,180) on odd degrees
(-179..179). kin2Dcont places each sample on a lattice node (no smoothing), copies the wrap edge, and contours each
cell split into four triangles at a centre point (mean of the corners), linear interpolation on triangle edges
(Zyda 1988 contour tree, contree.c). Missing samples count as zero, as in the lab's sparse Silk output.
The lab's own even-degree grid needs their input points, which are not published; the lab side is therefore the
NDFT grid contoured with kin2Dcont on its odd-degree lattice (-wrap -179 181 -179 181, as in silk/old/Makefile).
Levels: inner 98% (0.02) all classes; outer 99.95% (0.0005) General, 99.90% (0.001) Gly, IleVal, PrePro, TransPro,
99.85% (0.0015) CisPro. Our data: best rule (see contour_compare.py); Gly mirrored.
Comparisons per class and contour (distances from polylines resampled at 0.25°, periodic in phi and psi):
  A  lab NDFT vs ours, both odd grid, both kin2Dcont       (replaces the contourpy comparison)
  B  ours even grid (lab's exact procedure) vs ours odd    (effect of the 1° sampling offset)
  C  NDFT: kin2Dcont vs contourpy                          (effect of the contouring method)
Writes results/contour_kin2d.txt (stdout), results/contours_kin2d/<class>_<contour>.png and the kin2Dcont output
results/contours_kin2d/<class>_{ndft,ours_odd,ours_even}.kin.
Requires kincontours_build/src/kin2Dcont (rlabduke/kincontours 3da22b4, `make kin2Dcont`).
Run from rama8000_reconciliation/.
"""
import sys, os, subprocess, tempfile
sys.path.insert(0, 'scripts')
import numpy as np, pandas as pd
sys.argv = ['x', 'Gly', 'data/glycine.ndft', '5', '8']
src = open('scripts/silk_param_scan.py').read()
exec(src[:src.index("print(f'{cls}: realcount")])   # d (with amide-H cols)
import ndft, contourpy
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

KIN2D = 'kincontours_build/src/kin2Dcont'
OUT = 'results/contours_kin2d'
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
aw = d.omega.abs()
rule = ((d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0) & (d.own_mc_max < 30) & (d.own_mc_max > 0) & (d.own_minocc >= 1)
        & ~(d.H_minocc < 1) & d.resname.isin(STD | {'MSE'})
        & pd.Series(np.where(d.rama6 == 'CisPro', aw < 30, aw.isna() | (aw >= 150)), index=d.index))
setup = [('General', 'general', '6.5', 0.0005), ('Gly', 'glycine', '8', 0.001), ('IleVal', 'ileval', '8', 0.001),
         ('PrePro', 'prepro', '8', 0.001), ('TransPro', 'transpro', '8', 0.001), ('CisPro', 'cispro', '8', 0.0015)]
WRAP = {'odd': ['-179', '181', '-179', '181'], 'even': ['-180', '180', '-180', '180']}
BOUNDS = {'odd': '-180,180,-180,180', 'even': '-181,179,-181,179'}

def silk_text(df, tp, sym, grid):
    pts = list(zip(df.phi, df.psi)) + ([(-p, -q) for p, q in zip(df.phi, df.psi)] if sym else [])
    cmd = ['java', '-Xmx8g', '-cp', 'silk_build/classes', 'silk.SilkCmdLine', '-ndim=2', '-coords=2,3', '-wrap',
           '-gridsize=2', '-cosine=5', f'-twopass={tp}', f'-bounds={BOUNDS[grid]}', '-fraction', '-quiet']
    txt = '\n'.join(f'x {p} {q}' for p, q in pts) + '\n'
    return subprocess.run(cmd, input=txt.encode(), capture_output=True, check=True).stdout.decode()

def ndft_text(T):
    return ''.join(f'{ndft.center(i)} {ndft.center(j)} {float(T[i, j])!r}\n' for i, j in zip(*np.nonzero(T > 0)))

def kin2d(text, grid, levels, kinfile):
    cmd = [KIN2D, '-samp', '-last', '-nosmooth', '-group', '-wrap', *WRAP[grid], '-g2']
    for lv in levels: cmd += ['-l', repr(lv), 'sky']
    cmd.append('-')
    kin = subprocess.run(cmd, input=text.encode(), capture_output=True, check=True).stdout.decode()
    open(kinfile, 'w').write(kin)
    return parse_kin(kin)

def parse_kin(kin):
    """{level: [polyline arrays]} from kin2Dcont vectorlists (P starts a polyline, L continues it); axis box skipped."""
    out, cur, lines = {}, None, None
    for ln in kin.splitlines():
        if ln.startswith('@vectorlist'):
            name = ln.split('{')[1].split('}')[0]
            cur = None if name == 'axis' else float(name)
            if cur is not None: lines = out.setdefault(cur, [])
            continue
        if cur is None or not ln.startswith('{'): continue
        for tok in ln.split('}')[1:]:
            tok = tok.strip().split('{')[0].split()
            if len(tok) >= 3 and tok[0] in ('P', 'L'):
                xy = (float(tok[1]), float(tok[2]))
                if tok[0] == 'P' or not lines: lines.append([xy])
                else: lines[-1].append(xy)
    return {k: [np.array(s) for s in v] for k, v in out.items()}

def resample(segs):
    out = []
    for s in segs:
        r = [s[0]]
        for p in s[1:]:
            q = r[-1]; k = max(1, int(np.ceil(np.hypot(*(p - q)) / 0.25)))
            r += [q + (p - q) * t / k for t in range(1, k + 1)]
        out.append(np.array(r))
    return np.vstack(out) if out else np.empty((0, 2))

P = 4
ax_c = -179 + 2 * np.arange(-P, 180 + P)
def contourpy_lines(G, lev):
    z = np.pad(G, P, mode='wrap').T
    out = []
    for L in contourpy.contour_generator(ax_c, ax_c, z, line_type='Separate').lines(lev):
        r = resample([L]); keep = (np.abs(r) <= 180).all(1)          # resample, then clip: no chords across the seam
        out += [piece for piece in np.split(r, np.nonzero(np.diff(keep.astype(int)))[0] + 1) if len(piece) and (np.abs(piece) <= 180).all()]
    return out

def dist(A, B):
    if not len(A) or not len(B): return np.array([np.nan])
    return cKDTree((B + 180) % 360, boxsize=360).query((A + 180) % 360)[0]

def compare(sa, sb):
    a, b = resample(sa), resample(sb); dab, dba = dist(a, b), dist(b, a); both = np.concatenate([dab, dba])
    i = np.argmax(dab) if np.max(dab) >= np.max(dba) else None
    q = a[np.argmax(dab)] if i is not None else b[np.argmax(dba)]
    return np.mean(both), np.percentile(both, 95), max(np.max(dab), np.max(dba)), (q + 180) % 360 - 180

def plot_segs(ax, segs, style, label):
    for k, s in enumerate(segs):
        w = (s + 180) % 360 - 180                       # odd lattice runs to 181: wrap, break at the seam
        cut = np.nonzero(np.abs(np.diff(w, axis=0)).max(1) > 180)[0] + 1
        for m, piece in enumerate(np.split(w, cut)):
            ax.plot(piece[:, 0], piece[:, 1], style, lw=1.1, label=label if k == 0 and m == 0 else None)

os.makedirs(OUT, exist_ok=True)
print(f'{"class":9s}{"contour":>9s}{"level":>8s}{"n lines N/O/E":>15s}'
      f'{"A mean°":>9s}{"A p95°":>8s}{"A max°":>8s}{"B mean°":>9s}{"B p95°":>8s}{"B max°":>8s}{"C mean°":>9s}{"C p95°":>8s}{"C max°":>8s}')
for cls, f, tp, outer in setup:
    rc, T = ndft.load(f'data/{f}.ndft'); sym = cls == 'Gly'; S = d[rule & (d.rama6 == cls)]
    levels = [0.02, outer]
    KN = kin2d(ndft_text(T), 'odd', levels, f'{OUT}/{cls}_ndft.kin')
    KO = kin2d(silk_text(S, tp, sym, 'odd'), 'odd', levels, f'{OUT}/{cls}_ours_odd.kin')
    KE = kin2d(silk_text(S, tp, sym, 'even'), 'even', levels, f'{OUT}/{cls}_ours_even.kin')
    for name, lev in (('98%', 0.02), (f'{100 * (1 - outer):.2f}%', outer)):
        n_, o_, e_ = KN.get(lev, []), KO.get(lev, []), KE.get(lev, [])
        A, B, C = compare(o_, n_), compare(e_, o_), compare(n_, contourpy_lines(T, lev))
        print(f'{cls:9s}{name:>9s}{lev:8.4f}{f"{len(n_)}/{len(o_)}/{len(e_)}":>15s}'
              + ''.join(f'{m:9.2f}{p:8.2f}{x:8.2f}' for m, p, x, _ in (A, B, C)), flush=True)
        fig, ax = plt.subplots(figsize=(8.5, 8.5))
        plot_segs(ax, n_, 'k-', 'lab NDFT grid (odd degrees), kin2Dcont')
        plot_segs(ax, o_, 'r--', 'ours, odd-degree grid, kin2Dcont')
        plot_segs(ax, e_, 'g:', "ours, even-degree grid (lab's procedure), kin2Dcont")
        ax.plot(*A[3], 'o', ms=14, mfc='none', mec='m', mew=1.5, label=f'max deviation lab vs ours {A[2]:.2f}° at ({A[3][0]:.0f}, {A[3][1]:.0f})')
        ax.set(xlim=(-180, 180), ylim=(-180, 180), xticks=range(-180, 181, 30), yticks=range(-180, 181, 30), aspect='equal',
               xlabel='φ (°)', ylabel='ψ (°)',
               title=f'{cls}, {name} contour (level {lev:g}), kin2Dcont\nlab vs ours: mean {A[0]:.2f}°, max {A[2]:.2f}°;'
                     f' even vs odd grid: mean {B[0]:.2f}°, max {B[2]:.2f}°')
        ax.grid(alpha=0.25); ax.legend(loc='lower right', fontsize=8, framealpha=0.9)
        fig.tight_layout(); fig.savefig(f'{OUT}/{cls}_{name.rstrip("%").replace(".", "p")}.png', dpi=130); plt.close(fig)
