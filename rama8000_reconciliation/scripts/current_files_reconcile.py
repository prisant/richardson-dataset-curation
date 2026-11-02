"""How well do our current files reproduce the lab's 2013 rama8000 residue lists, and why do they differ?

Two file sets for the same 7,957 chains:
  2011:    the lab's top8000_chains_70.tgz chain files (Reduce 3.15)        tables in data/lab2011/
  current: our RCSB downloads (~/Desktop/Data/top8000_pdbs_hom70, Reduce 4.16)  tables in data/
Final rule (lab2013_compare.py 'same, mainchain B > 1'), with two forms of the occupancy test:
  all-atom: every atom of the residue in the file at occupancy 1 (includes Reduce hydrogens)
  heavy:    own heavy atoms at occupancy 1, and for non-Pro residues no altloc on the preceding C
            (heavy-atom form of the split-amide-H test; needs no hydrogens)
Each residue that differs between ours (current files) and the lab lists is assigned one cause:
  rule      - the same rule on the 2011 files also disagrees with the lab (rule not yet matched)
  file      - the 2011 files agree with the lab; the cause is the rule component that changed between the
              two file versions (numbering / residue name / class / coordinates / B / occupancy / link / altloc / H)
Run from rama8000_reconciliation/:  python scripts/current_files_reconcile.py > results/current_files_reconcile.txt
"""
import re, numpy as np, pandas as pd

K = ['pdb', 'chain', 'resnum', 'ins']
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
CLS = ['General', 'Gly', 'IleVal', 'PrePro', 'TransPro', 'CisPro']
LAB = 'data/Top8000_Rama_2013'
FILES = dict(General='general-noGPIVpreP', IleVal='ileval-nopreP', Gly='gly', PrePro='prepro-noGP',
             TransPro='transpro', CisPro='cispro')
LBL = re.compile(r'^(\w{4}) (.) (-?\d+)(\D?) (\w{3}) B<([\d.]+) wc(\S+) \(([\d.]+)A MP=([\d.]+)\)$')
rows = []
for cls, f in FILES.items():
    for ln in open(f'{LAB}/residues/rama8000-{f}.csv'):
        lab = ln.rstrip('\n').rsplit(',', 4)[0]
        p, ch, rn, ins, rt = LBL.match(lab).groups()[:5]
        rows.append((p, ch, int(rn), ins.strip(), rt, cls))
L = pd.DataFrame(rows, columns=K + ['lab_resname', 'lab_cls'])


def load(D):
    rd = lambda f: pd.read_csv(f'{D}/{f}', sep='\t', keep_default_na=False, na_values=['nan'],
                               dtype={'ins': str, 'chain': str, 'pdb': str})
    t = rd('residue_records.tsv'); t['ins'] = t.ins.str.strip()
    for f in ('amideH_records.tsv', 'allatom_records.tsv', 'prevC_spread.tsv', 'phipsi_atom_records.tsv'):
        u = rd(f); u['ins'] = u.ins.str.strip(); t = t.merge(u, on=K, how='left')
    return t


def parts(t):
    """Rule components (True = passes)."""
    ew = t.omega.abs()
    return pd.DataFrame({
        'link': (t.d_prevC_N < 2) & (t.d_C_nextN < 2),
        'B': (t.own_mc_max > 1) & (t.own_mc_max < 30) & ~((t.resname == 'MSE') & (np.fmax(t.own_mc_max, t.sc_maxB_heavy) >= 30)),
        'heavy occ': t.own_minocc >= 1,
        'prev C altloc': ~((t.prevC_alt > 0) & (t.resname != 'PRO')),
        'H occ': ~(t.all_minocc < 1),
        'omega': t.omega.notna() & np.where(t.rama6 == 'CisPro', ew < 30, ew >= 150),
        'residue type': t.resname.isin(STD | {'MSE'}) & ~(t.rama6.isin(['General', 'IleVal', 'PrePro'])
                                                          & ~t.next_resname.isin(STD | {'MSE'}))}, index=t.index)


COMMON = ['link', 'B', 'heavy occ', 'omega', 'residue type']
RULES = {'all-atom': COMMON + ['H occ'], 'heavy': COMMON + ['prev C altloc']}

T = {}
for nm, D in (('2011', 'data/lab2011'), ('current', 'data')):
    t = load(D); P = parts(t)
    t = pd.concat([t[K + ['resname', 'rama6', 'phi', 'psi']], P], axis=1)
    T[nm] = t
lab = L.set_index(K).lab_cls

print('ours vs the lab 2013 lists (1,504,897 residues); same chains, same final rule, two file sets and two occupancy tests')
print(f'{"files / occupancy test":26s}' + ''.join(f'{c:>10s}' for c in CLS) + f'{"total":>8s}{"ours-only":>11s}{"lab-only":>10s}')
nd = L.lab_cls.value_counts()
for fn in ('2011', 'current'):
    for rn, cols in RULES.items():
        t = T[fn]; v = t[cols].all(axis=1)
        s = t[v].set_index(K).rama6
        j = pd.concat([s.rename('ours'), lab.rename('lab')], axis=1)
        same = j.ours == j.lab
        dl = [int((s == c).sum()) - int(nd[c]) for c in CLS]
        print(f'{fn + ", " + rn:26s}' + ''.join(f'{x:>+10,d}' for x in dl) + f'{sum(dl):>+8,d}'
              f'{int((j.ours.notna() & ~same).sum()):11,d}{int((j.lab.notna() & ~same).sum()):10,d}')

# ---- attribution of every difference (ours on current files vs the lab) ---------------------------------------------
for rn, cols in RULES.items():
    o, c = T['2011'].set_index(K), T['current'].set_index(K)
    so = o.rama6.where(o[cols].all(axis=1)); sc = c.rama6.where(c[cols].all(axis=1))
    J = pd.concat([lab.rename('lab'), so.rename('r2011'), sc.rename('cur')], axis=1)
    J = J.join(o[['resname', 'rama6', 'phi', 'psi'] + cols].add_prefix('o_')).join(c[['resname', 'rama6', 'phi', 'psi'] + cols].add_prefix('c_'))
    J = J[J.index.isin(lab.index) | J.cur.notna()]
    diff = J.cur.fillna('-') != J.lab.fillna('-')
    D = J[diff].copy()
    D['side'] = np.select([D.lab.isna(), D.cur.isna()], ['ours-only', 'lab-only'], 'both, class differs')
    why = []
    for k, r in D.iterrows():
        if (r.r2011 if pd.notna(r.r2011) else '-') != (r.lab if pd.notna(r.lab) else '-'):
            why.append('rule: same disagreement on 2011 files'); continue
        if pd.isna(r.c_resname) or pd.isna(r.o_resname):
            why.append('file: residue number absent in one file'); continue
        if r.c_resname != r.o_resname:
            why.append('file: residue name changed'); continue
        if r.c_rama6 != r.o_rama6:
            why.append('file: class changed (next residue changed)'); continue
        flips = [x for x in cols if bool(r['o_' + x]) != bool(r['c_' + x])]
        dphi = abs((r.c_phi - r.o_phi + 180) % 360 - 180) if pd.notna(r.c_phi) and pd.notna(r.o_phi) else 0
        tag = 'file: ' + ' + '.join(flips) if flips else 'file: other'
        if dphi > 0.01: tag += ' (coordinates changed)'
        why.append(tag)
    D['why'] = why
    D['entries'] = [k[0] for k in D.index]
    print(f'\n== occupancy test: {rn}; ours on current files vs the lab: {len(D):,} residues differ '
          f'(ours-only {int((D.side == "ours-only").sum()):,}, lab-only {int((D.side == "lab-only").sum()):,},'
          f' same residue number but a different class {int((D.side == "both, class differs").sum()):,})')
    g = D.groupby(['why', 'side']).size().unstack(fill_value=0)
    g['entries'] = D.groupby('why').entries.nunique()
    g['top entries'] = D.groupby('why').entries.agg(lambda s: ', '.join(f'{a} {b}' for a, b in s.value_counts().head(4).items()))
    print(g.assign(n=g.iloc[:, :-2].sum(axis=1)).sort_values('n', ascending=False).drop(columns='n').to_string())

# ---- residues the lab kept although the preceding C has altloc records (2011 files, heavy rule otherwise passed) -----
t = load('data/lab2011'); P = parts(t)
x = t[P[COMMON].all(axis=1) & ~P['prev C altloc']].merge(L[K], on=K, how='left', indicator=True)
kept = x._merge == 'both'
grp = np.select([x.prevC_n == 1, x.prevC_spread == 0], ['single altloc record', 'alternates at the same position'],
                'alternates at different positions')
print('\n2011 files: non-Pro residues passing the heavy rule except the preceding-C altloc test, by kind of altloc')
print(pd.crosstab(grp, kept.map({True: 'lab kept', False: 'lab dropped'})).to_string())
