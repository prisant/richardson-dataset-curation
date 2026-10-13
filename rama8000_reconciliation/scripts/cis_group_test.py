"""Add groups of excluded cis-Pro candidates to the best-fit set; score Silk(c5,t8) grid vs cispro.ndft."""
import sys, json, subprocess, tempfile, os, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, 'scripts'); import ndft
_, T = ndft.load('data/cispro.ndft')
R = [json.loads(l) for l in open('data/pro_records.jsonl')]
dB = lambda r: max(r['own_mc'], r['prevC'], r['nextN'])
cis = [r for r in R if abs(r['omega']) < 30]
best = [r for r in cis if dB(r) < 30 and not r['alt']]
excl = [r for r in cis if r not in best]
def line(r): return f"{r['pdb']}_{r['chain']}{r['resnum']} {r['phi']} {r['psi']}\n"
def score(rows):
    fd, fn = tempfile.mkstemp(suffix='.ndft'); os.close(fd)
    cmd = ['java','-cp','silk_build/classes','silk.SilkCmdLine','-ndim=2','-coords=2,3','-label=1','-bounds=-180,180,-180,180',
           '-wrap','-gridsize=2','-cosine=5','-twopass=8','-fraction','-ndft','-quiet']
    open(fn,'wb').write(subprocess.run(cmd, input=''.join(map(line, rows)).encode(), capture_output=True, check=True).stdout)
    _, g = ndft.load(fn); os.remove(fn)
    return float(np.abs(g - T).mean() * 1e4), int(((g > 0) != (T > 0)).sum())
groups = {
  'baseline (best-fit, n=3609)': [],
  'alt, B ok, mainchain alt': [r for r in excl if r['alt'] and dB(r) < 30 and r['mcalt']],
  'alt, B ok, sidechain alt': [r for r in excl if r['alt'] and dB(r) < 30 and not r['mcalt']],
  'no alt, B in [30,31)': [r for r in excl if not r['alt'] and 30 <= dB(r) < 31],
  'no alt, B in [31,32)': [r for r in excl if not r['alt'] and 31 <= dB(r) < 32],
  'no alt, B in [32,35)': [r for r in excl if not r['alt'] and 32 <= dB(r) < 35],
  'no alt, B >= 35': [r for r in excl if not r['alt'] and dB(r) >= 35],
  'own mc<30 but prevC/nextN>=30': [r for r in excl if not r['alt'] and r['own_mc'] < 30 and dB(r) >= 30],
  'own mc>=30, (N,CA,C<30, O>=30)': [],
}
# B on O only: need atom-level; approximate with own_mc>=30 while prevC,nextN<30
groups['own mc>=30, neighbours ok'] = [r for r in excl if not r['alt'] and r['own_mc'] >= 30 and max(r['prevC'], r['nextN']) < 30]
del groups['own mc>=30, (N,CA,C<30, O>=30)']
items = list(groups.items())
with ThreadPoolExecutor(8) as ex:
    scores = list(ex.map(lambda kv: score(best + kv[1]), items))
print(f'{"group":34s}{"n added":>8s}{"total":>7s}{"meanΔ×1e4":>11s}{"supp mism":>10s}')
for (k, v), (s, m) in zip(items, scores):
    print(f'{k:34s}{len(v):8d}{3609+len(v):7d}{s:11.3f}{m:10d}')
