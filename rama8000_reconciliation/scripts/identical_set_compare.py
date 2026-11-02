"""Category counts on the residues whose filtering inputs are identical in the 2011 and current files
(data/identical_inputs.tsv, identical_inputs.py): the lab's 2013 lists vs our heavy-atom rule on each file set.

On this set our rule must give the same residues on both file sets (a check of the identical-input definition); any
remaining difference from the lab lists is a rule difference. Residues outside the set are counted separately.
Rule, lab lists and loaders: current_files_reconcile.py ('heavy' occupancy test).
Optional argument: the rule name in RULES (default 'heavy'). 'all-atom' is the occupancy test we invented that also
requires the Reduce hydrogens at occupancy 1; it is not one of the lab's rules (hydrogens are not filtering inputs, so on
this set the two file sets can differ under it).
Usage (from rama8000_reconciliation/):
  python scripts/identical_set_compare.py > results/identical_set_compare.txt
  python scripts/identical_set_compare.py all-atom > results/identical_set_compare_allatom.txt
"""
import sys
import numpy as np, pandas as pd

src = open('scripts/current_files_reconcile.py').read()
exec(src[:src.index('T = {}')])          # K, CLS, L (lab lists), load, parts, RULES

I = pd.read_csv('data/identical_inputs.tsv', sep='\t', keep_default_na=False, dtype={'chain': str, 'ins': str, 'resnum': int})
I = I.set_index(K).identical.astype(bool)
RULE = sys.argv[1] if len(sys.argv) > 1 else 'heavy'
cols = RULES[RULE]
print(f'rule: {RULE} ({" + ".join(cols)})')
ours = {}
for nm, D in (('2011', 'data/lab2011'), ('current', 'data')):
    t = load(D); t['resnum'] = t.resnum.astype(int); t = t.set_index(K)
    ours[nm] = t.rama6.where(parts(t.reset_index())[cols].all(axis=1).values).dropna()
lab = L.assign(resnum=L.resnum.astype(int)).set_index(K).lab_cls

idn = set(I[I].index)
inI = lambda s: s[s.index.isin(idn)]
outI = lambda s: s[~s.index.isin(idn)]
print(f'residues with identical inputs: {len(idn):,} of {len(I):,}; not identical: {int((~I).sum()):,}')
print(f'lab-listed residues: inside the set {len(inI(lab)):,}; outside {len(outI(lab)):,}')
o11, ocur = inI(ours['2011']), inI(ours['current'])
same = set(o11.index) == set(ocur.index) and bool((o11 == ocur.reindex(o11.index)).all())
print(f'check: our rule on the identical set, 2011 vs current files: {len(o11):,} vs {len(ocur):,} residues;'
      f' identical selection and classes: {same}')

def table(o, nm):
    print(f'\ncategory counts on the identical set (lab 2013 lists vs ours, {nm} files)')
    print(f'{"category":10s}{"lab":>11s}{"ours":>11s}{"Δ":>7s}{"only ours":>11s}{"only lab":>10s}')
    j = pd.concat([inI(lab).rename('lab'), o.rename('ours')], axis=1)
    for c in CLS + ['total']:
        m = (j.lab == c) | (j.ours == c) if c != 'total' else pd.Series(True, index=j.index)
        nl = int((j.lab == c).sum()) if c != 'total' else int(j.lab.notna().sum())
        no = int((j.ours == c).sum()) if c != 'total' else int(j.ours.notna().sum())
        jo = j[m]
        print(f'{c:10s}{nl:11,d}{no:11,d}{no - nl:>+7,d}{int((jo.ours.notna() & (jo.ours != jo.lab)).sum()):11,d}'
              f'{int((jo.lab.notna() & (jo.ours != jo.lab)).sum()):10,d}')
    return j


table(o11, '2011')
j = table(ocur, 'current')
if not same:
    a, b = set(o11.index), set(ocur.index)
    print(f'\n2011 vs current selections on the identical set: only 2011 {len(a - b):,}; only current {len(b - a):,};'
          f' entries {len({k[0] for k in a ^ b}):,}')

print('\noutside the identical set (file differences): lab-listed {:,}; ours (current files) {:,}; ours (2011 files) {:,}'.format(
    len(outI(lab)), len(outI(ours['current'])), len(outI(ours['2011']))))
print('  outside-set residues by entry (top 12):', dict(pd.Series([k[0] for k in I[~I].index]).value_counts().head(12)))

# ---- the rule differences on the identical set, by category and reason ---------------------------------------------
t = load('data'); t['resnum'] = t.resnum.astype(int)
P = parts(t)[cols]
t = t.join(P.add_prefix('ok_')).set_index(K)
L2 = L.assign(resnum=L.resnum.astype(int)).set_index(K)
d = j[j.ours != j.lab].join(t, how='left').join(L2[['lab_resname']], how='left')
d['side'] = d.ours.isna().map({True: 'only lab', False: 'only ours'})
d['cat'] = d.lab.fillna(d.ours)
fail = lambda r: ' + '.join(c for c in cols if not r['ok_' + c]) or '(passes)'
d['why'] = [fail(r) if s == 'only lab' else '' for (_, r), s in zip(d.iterrows(), d.side)]
print('\nonly lab: the laboratory kept it; the parts of our rule it fails')
print(d[d.side == 'only lab'].groupby(['why', 'cat']).size().unstack(fill_value=0).to_string())
o = d[d.side == 'only ours']
sh = t[t.index.isin(j[(j.ours == j.lab)].index)]
F = {'own altloc label (occupancy 1)': lambda x: x.own_alt > 0,
     'preceding N/CA/O altloc': lambda x: (x[['prevN_alt', 'prevCA_alt', 'prevO_alt']] > 0).any(axis=1),
     'next N or CA altloc': lambda x: (x[['nextN_alt', 'nextCA_alt']] > 0).any(axis=1),
     'next N or CA occupancy < 1': lambda x: (x[['nextN_occ', 'nextCA_occ']] < 1).any(axis=1),
     'side-chain heavy B >= 30': lambda x: x.sc_maxB_heavy >= 30,
     'mainchain B 25-30': lambda x: x.own_mc_max >= 25,
     'MSE or next residue MSE': lambda x: (x.resname == 'MSE') | (x.next_resname == 'MSE')}
print(f'\nonly ours ({len(o)}): features, count (share) vs share among the {len(sh):,} residues both keep')
for k, f in F.items():
    print(f'  {k:32s} {int(f(o).sum()):4d} ({f(o).mean():6.1%})   both keep: {f(sh).mean():6.2%}')
print('  only ours with none of these:', int((~pd.concat([f(o) for f in F.values()], axis=1).any(axis=1)).sum()))
for s in ('only lab', 'only ours'):
    x = d[d.side == s]
    print(f'\n{s}: {len(x)} residues in {x.index.get_level_values(0).nunique()} entries; top entries',
          dict(pd.Series(x.index.get_level_values(0)).value_counts().head(8)))
    print(f'  residue types: {dict(x.resname.fillna(x.lab_resname).value_counts().head(8))}')

x = d[d.side == 'only lab']
print('\nonly lab, reason by MSE:'); print(pd.crosstab(x.why, (x.resname.fillna(x.lab_resname) == 'MSE').map({True: 'MSE', False: 'other'})).to_string())
y = x[x.why.str.contains('prev C')]
print(f'only lab failing the preceding-C altloc test: {len(y)}; single altloc record {int((y.prevC_n == 1).sum())},'
      f' alternates coincide {int(((y.prevC_n > 1) & (y.prevC_spread == 0)).sum())},'
      f' alternates at different positions {int(((y.prevC_n > 1) & (y.prevC_spread > 0)).sum())}')
# residue-number skip to the next residue in the residue table (peptide-linked, number not n+1)
t2 = load('data'); t2['resnum'] = t2.resnum.astype(int)
t2['skip'] = (t2.groupby(['pdb', 'chain']).resnum.shift(-1) - t2.resnum) != 1
sk = t2.set_index(K).skip
for nm, idx in (('only ours', d[d.side == 'only ours'].index), ('only lab', x.index), ('both keep', j[j.ours == j.lab].index)):
    v = sk.reindex(idx)
    print(f'next residue number not n+1: {nm} {int(v.sum()):,} / {len(v):,}')
