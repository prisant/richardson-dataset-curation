"""TransPro: omega cutoff scan near 150, then edge groups on the |ω|>=150 rule set."""
import sys, json, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.argv = ['x', 'data/pro_records_ext.jsonl', 'transpro', 'data/transpro.ndft', '8']
exec(open('scripts/class_rule_check.py').read().split("rules = {")[0])
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
F0 = [r for r in P if r['own_mc'] < 30 and r['own_minocc'] >= 1]
def stats(rows):
    g = silk(rows); d = np.rint(g * rc) - np.rint(T * rc)
    return len(rows), np.abs(g - T).mean() * 1e4, int(((g > 0) != (T > 0)).sum()), int((d == 0).sum()), int(np.abs(d).max())
F = [r for r in F0 if abs(r['omega']) >= 150]
jobs = [(f'|ω|>={c}', [r for r in F0 if abs(r['omega']) >= c]) for c in (140, 145, 148, 149, 150, 151, 152, 155)]
rem = {
  'prev/next mainchain altloc': lambda r: r['prev_mcalt'] or r['next_mcalt'],
  'neighbour partial occupancy': lambda r: r['nbr_minocc'] < 1,
  'own mc B in [29,30)': lambda r: 29 <= r['own_mc'] < 30,
  'own CB B >= 30': lambda r: r['own_cb'] >= 30,
  'own any atom B >= 30': lambda r: r['own_all'] >= 30,
  'prev or next non-standard': lambda r: r['resname_prev'] not in STD or r['resname_next'] not in STD,
  'peptide C-N > 1.45 A': lambda r: max(r['d_prevC_N'], r['d_C_nextN']) > 1.45,
  'phi > -20 (implausible Pro)': lambda r: r['phi'] > -20,
}
for k, f in rem.items(): jobs.append(('[150] remove ' + k, [r for r in F if not f(r)]))
with ThreadPoolExecutor(8) as ex: res = list(ex.map(lambda j: stats(j[1]), jobs))
print(f'{"variant":44s}{"n":>8s}{"Δn":>7s}{"meanΔ×1e4":>11s}{"supp":>6s}{"rank-identical":>16s}{"max|Δr|":>9s}')
for (name, rows), (n, m, s, e, mx) in zip(jobs, res):
    print(f'{name:44s}{n:8d}{n-rc:+7d}{m:11.3f}{s:6d}{e:>10d}/32400{mx:9d}')
