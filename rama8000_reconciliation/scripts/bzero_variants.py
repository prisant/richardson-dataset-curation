"""Zero-B rule variants on TransPro (|ω|>=150) and CisPro; needs per-atom minimum B, so re-parse affected PDBs."""
import sys, json, numpy as np
sys.path.insert(0, 'scripts'); import ndft
from cis_records_ext import parse
R = [json.loads(l) for l in open('data/pro_records_ext.jsonl')]
link = lambda r: r['d_prevC_N'] < 2.0 and r['d_C_nextN'] < 2.0
base = [r for r in R if link(r) and r['own_mc'] < 30 and r['own_minocc'] >= 1]
cache = {}
def minb(r):
    if r['pdb'] not in cache: cache[r['pdb']] = parse(r['pdb'])[0]
    cur = cache[r['pdb']][(r['chain'], r['resnum'], r['ins'])]
    mc = [min(cur['b'][a]) for a in ('N','CA','C','O') if a in cur['b']]
    return min(mc), min(min(v) for v in cur['b'].values())
zero_pdbs = {'2fcr', '2yyh'}
# any-zero test only matters where some B is 0; scan all pdbs cheaply via own_mc<2 shortlist plus the two known pdbs
short = [r for r in base if r['own_mc'] < 2 or r['pdb'] in zero_pdbs]
for r in short: r['_minmc'], r['_minall'] = minb(r)
get = lambda r, k: r.get(k, 99.0)
variants = {
  'none':                          lambda r: True,
  'drop mainchain all B==0':       lambda r: r['own_mc'] > 0,
  'drop any mainchain atom B==0':  lambda r: get(r, '_minmc') > 0,
  'drop any atom B==0':            lambda r: get(r, '_minall') > 0,
}
for cls, sel, path in [('transpro', lambda r: abs(r['omega']) >= 150, 'data/transpro.ndft'),
                       ('cispro', lambda r: abs(r['omega']) < 30, 'data/cispro.ndft')]:
    rc, _ = ndft.load(path)
    for name, f in variants.items():
        n = sum(1 for r in base if sel(r) and f(r))
        print(f'{cls:9s} {name:30s} n={n} Δn={n-rc:+d}')

# grid comparison for TransPro, two surviving variants
sys.argv = ['x', 'data/pro_records_ext.jsonl', 'transpro', 'data/transpro.ndft', '8']
exec(open('scripts/class_rule_check.py').read().split("rules = {")[0])
for name in ('drop mainchain all B==0', 'drop any mainchain atom B==0'):
    rows = [r for r in base if abs(r['omega']) >= 150 and variants[name](r)]
    g = silk(rows); d = np.rint(g * rc) - np.rint(T * rc)
    print(f'grid {name:30s} n={len(rows)} meanΔ×1e4={np.abs(g-T).mean()*1e4:.4f} supp={int(((g>0)!=(T>0)).sum())} '
          f'identical={int((d==0).sum())} max|Δrank|={int(np.abs(d).max())}')
