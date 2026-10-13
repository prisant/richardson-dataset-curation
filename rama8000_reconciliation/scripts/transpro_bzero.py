"""TransPro: test exclusion of residues with zero (or near-zero) B-factors, rule |ω|>=150."""
import sys, json, numpy as np
sys.argv = ['x', 'data/pro_records_ext.jsonl', 'transpro', 'data/transpro.ndft', '8']
exec(open('scripts/class_rule_check.py').read().split("rules = {")[0])
F = [r for r in P if r['own_mc'] < 30 and r['own_minocc'] >= 1 and abs(r['omega']) >= 150]
sys.path.insert(0, 'scripts'); from cis_records_ext import parse
def minB(r):
    res, _ = parse(r['pdb']); cur = res[(r['chain'], r['resnum'], r['ins'])]
    return min(min(v) for v in cur['b'].values()), min(min(cur['b'][a]) for a in ('N','CA','C','O') if a in cur['b'])
low = [r for r in F if r['own_mc'] < 2]
print(f'TransPro rule set residues with own mainchain max B < 2: {len(low)}')
for r in low: print(f"   {r['pdb']} {r['chain']}{r['resnum']} own_mc={r['own_mc']} (min all, min mc)={minB(r)} φ={r['phi']} ψ={r['psi']}")
pdbs = sorted({r['pdb'] for r in low})
for name, rows in [('rule set', F), ('remove own mc max B == 0', [r for r in F if r['own_mc'] > 0]),
                   ('remove all residues of pdbs above', [r for r in F if r['pdb'] not in pdbs])]:
    g = silk(rows); d = np.rint(g * rc) - np.rint(T * rc)
    print(f'{name:36s} n={len(rows)} Δn={len(rows)-rc:+d} meanΔ×1e4={np.abs(g-T).mean()*1e4:.4f} supp={int(((g>0)!=(T>0)).sum())} identical={int((d==0).sum())}')
