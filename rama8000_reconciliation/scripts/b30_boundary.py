"""B-factor boundary test: own mainchain max B < 30 (NDFT titles) vs <= 30 (Williams 2018 "excluded ... B-factor > 30").

Part 1: TransPro grid at Makefile params (-cosine=5 -twopass=8) under the reconciled rule
        (0 < own mc B, all atoms occ 1, |omega| >= 150), both boundary forms.
Part 2: all six class counts under the full reconciled rule (incl. amide-H occ 1), both boundary forms.
Run from rama8000_reconciliation/.
"""
import sys, json, numpy as np, pandas as pd
sys.argv = ['x', 'data/pro_records_ext.jsonl', 'transpro', 'data/transpro.ndft', '8']
exec(open('scripts/class_rule_check.py').read().split("rules = {")[0])

base = lambda r: r['own_mc'] > 0 and r['own_minocc'] >= 1 and abs(r['omega']) >= 150
at30 = [r for r in P if base(r) and r['own_mc'] == 30]
print(f'TransPro: realcount {rc}; rule-set residues with own mc max B == 30.00: {len(at30)}')
for r in at30: print(f"   {r['pdb']} {r['chain']}{r['resnum']}{r['ins'].strip()} φ={r['phi']} ψ={r['psi']} ω={r['omega']}")
print(f'{"boundary":10s}{"n":>8s}{"Δn":>6s}{"meanΔ×1e4":>11s}{"supp":>6s}{"rank-identical":>16s}{"max|Δrank|":>11s}')
for name, f in (('B < 30', lambda r: r['own_mc'] < 30), ('B <= 30', lambda r: r['own_mc'] <= 30)):
    rows = [r for r in P if base(r) and f(r)]
    g = silk(rows); d = np.rint(g * rc) - np.rint(T * rc)
    print(f'{name:10s}{len(rows):8d}{len(rows)-rc:+6d}{np.abs(g-T).mean()*1e4:11.4f}{int(((g>0)!=(T>0)).sum()):6d}'
          f'{int((d==0).sum()):>10d}/{d.size}{int(np.abs(d).max()):11d}')

d = pd.read_csv('data/residue_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
h = pd.read_csv('data/amideH_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
d = d.merge(h, on=['pdb', 'chain', 'resnum', 'ins'], how='left')
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
common = ((d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0) & (d.own_minocc >= 1) & (d.own_mc_max > 0) & d.resname.isin(STD)
          & ~((d.rama6 == 'TransPro') & (d.omega.abs() < 150)) & ~(d.H_minocc < 1))
nd = dict(General=1061639, Gly=116589, IleVal=195083, PrePro=61568, TransPro=66364, CisPro=3654)
print(f'\nAll classes, Δ vs NDFT realcount (Gly halved), full reconciled rule')
print(f'{"boundary":10s}' + ''.join(f'{c:>10s}' for c in nd) + f'{"total":>10s}')
for name, m in (('B < 30', common & (d.own_mc_max < 30)), ('B <= 30', common & (d.own_mc_max <= 30))):
    dl = [int((m & (d.rama6 == c)).sum()) - n for c, n in nd.items()]
    print(f'{name:10s}' + ''.join(f'{x:>+10,d}' for x in dl) + f'{sum(dl):>+10,d}')
