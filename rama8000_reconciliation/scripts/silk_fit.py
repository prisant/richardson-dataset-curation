"""Run Silk on a point file and compare its -fraction NDFT grid to an archived rama8000 NDFT."""
import io, subprocess, sys, itertools, tempfile, numpy as np
sys.path.insert(0, 'scripts'); import ndft

CP = 'silk_build/classes'
def run_silk(pts_text, kernel, h1, k=None, lam=None, extra=()):
    cmd = ['java', '-cp', CP, 'silk.SilkCmdLine', '-ndim=2', '-coords=2,3', '-label=1',
           '-bounds=-180,180,-180,180', '-wrap', '-gridsize=2', f'-{kernel}={h1}',
           '-fraction', '-ndft', '-quiet', '-title=fit']
    if k is not None: cmd.append(f'-twopass={k}')
    if lam is not None: cmd.append(f'-lambda={lam}')
    cmd += list(extra)
    out = subprocess.run(cmd, input=pts_text.encode(), capture_output=True, check=True).stdout
    with tempfile.NamedTemporaryFile(suffix='.ndft') as tmp:
        tmp.write(out); tmp.flush()
        return ndft.load(tmp.name)

def compare(g, t, n):
    a, b = np.rint(g * n), np.rint(t * n)          # integer ranks
    supp = ((g > 0) != (t > 0)).sum()
    return dict(support_mismatch=int(supp), exact=int((a == b).sum()),
                mean_abs_rank=float(np.abs(a - b).mean()), max_abs_rank=int(np.abs(a - b).max()))

if __name__ == '__main__':
    pts = open('data/cispro_bestfit.txt').read()
    rc, tgt = ndft.load('data/cispro.ndft')
    n = pts.count('\n')
    res = []
    for kernel, h1, k in itertools.product(['cosine', 'gaussian'], [5, 8, 10, 12, 15], [None, 10, 13, 16, 20, 25]):
        try: c, g = run_silk(pts, kernel, h1, k)
        except subprocess.CalledProcessError as e: print(kernel, h1, k, 'ERR', e.stderr[:200]); continue
        m = compare(g, tgt, n); res.append((m['mean_abs_rank'], kernel, h1, k, m))
    res.sort(key=lambda x: x[0])
    for r in res[:12]: print(r[1], r[2], r[3], r[4])
