import sys, itertools, numpy as np
sys.path.insert(0, 'scripts'); import ndft
from silk_fit import run_silk
from concurrent.futures import ThreadPoolExecutor
pts = open('data/cispro_bestfit.txt').read()
rc, tgt = ndft.load('data/cispro.ndft')
def score(p):
    kernel, h1, k = p
    _, g = run_silk(pts, kernel, h1, k) if False else (None, None)
    return p
import subprocess, tempfile, os
def run(p):
    kernel, h1, k = p
    fd, fn = tempfile.mkstemp(suffix='.ndft'); os.close(fd)
    cmd = ['java','-cp','silk_build/classes','silk.SilkCmdLine','-ndim=2','-coords=2,3','-label=1',
           '-bounds=-180,180,-180,180','-wrap','-gridsize=2',f'-{kernel}={h1}',f'-twopass={k}','-fraction','-ndft','-quiet']
    open(fn,'wb').write(subprocess.run(cmd,input=pts.encode(),capture_output=True,check=True).stdout)
    _, g = ndft.load(fn); os.remove(fn)
    d = np.abs(g - tgt)
    return (float(d.mean()*1e4), int(((g>0)!=(tgt>0)).sum()), int(((g>0)&(tgt==0)).sum()), int(((g==0)&(tgt>0)).sum()), kernel, h1, k)
grid = list(itertools.product(['cosine','gaussian'], [2,3,4,5,6,7], [5,6,7,8,9,10,11,12]))
with ThreadPoolExecutor(8) as ex: res = sorted(ex.map(run, grid))
print('meanΔfrac×1e4  supp_mismatch  ours-only  ref-only  kernel h1 k')
for r in res[:15]: print(f'{r[0]:12.2f} {r[1]:13d} {r[2]:10d} {r[3]:9d}   {r[4]} {r[5]} {r[6]}')
