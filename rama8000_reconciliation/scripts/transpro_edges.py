"""TransPro edge tests at Makefile params: twisted-omega exclusion and single outlier removal."""
import sys, json, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.argv = ['x', 'data/pro_records_ext.jsonl', 'transpro', 'data/transpro.ndft', '8']
exec(open('scripts/class_rule_check.py').read().split("rules = {")[0])
F = [r for r in P if r['own_mc'] < 30 and r['own_minocc'] >= 1]
key = lambda r: f"{r['pdb']} {r['chain']}{r['resnum']}"
def stats(rows):
    g = silk(rows); d = np.rint(g * rc) - np.rint(T * rc)
    return len(rows), np.abs(g - T).mean() * 1e4, int(((g > 0) != (T > 0)).sum()), int((d == 0).sum()), int(np.abs(d).max())
jobs = [('rule (|ω|>=30)', F)]
for lo in (90, 120, 150, 160):
    jobs.append((f'trans = |ω|>={lo}', [r for r in F if abs(r['omega']) >= lo]))
for k in ('1eq9 B519', '2asd A160', '3pve A1574'):
    jobs.append((f'remove {k}', [r for r in F if key(r) != k]))
jobs.append(('remove all three', [r for r in F if key(r) not in ('1eq9 B519', '2asd A160', '3pve A1574')]))
with ThreadPoolExecutor(8) as ex: res = list(ex.map(lambda j: stats(j[1]), jobs))
print(f'{"variant":26s}{"n":>8s}{"Δn":>7s}{"meanΔ×1e4":>11s}{"supp":>6s}{"rank-identical":>16s}{"max|Δr|":>9s}')
for (name, _), (n, m, s, e, mx) in zip(jobs, res):
    print(f'{name:26s}{n:8d}{n-rc:+7d}{m:11.3f}{s:6d}{e:>10d}/32400{mx:9d}')
tw = sorted([r for r in F if abs(r['omega']) < 150], key=lambda r: abs(r['omega']))
print(f'\nPro in rule set with 30<=|ω|<150: {len(tw)}')
for r in tw: print(f"   {key(r):14s} ω={r['omega']:7.1f} φ={r['phi']:7.1f} ψ={r['psi']:7.1f}")
