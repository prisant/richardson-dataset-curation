"""Per-residue differences for one class between the lab's 2013 lists and our rule (prev-C rule + MSE B + ω defined).

usage: lab2013_class_detail.py CLASS [RULE]   (General, Gly, IleVal, PrePro, TransPro, CisPro; RULE a variant name
       from lab2013_compare.py, default 'prev-C + MSE B + ω defined')
Lists every residue in the lab's CLASS list but not in ours, and vice versa, with the class the other side gives it
(if any), the first of our filters it fails, and the values behind it. Includes the four renumbered entries.
Run from rama8000_reconciliation/.
"""
import sys, io, contextlib
CLS, ARGV = sys.argv[1], list(sys.argv)
with contextlib.redirect_stdout(io.StringIO()):
    exec(open('scripts/lab2013_compare.py').read().split('# ---- Evidence for the refined rules')[0])
RULE = ARGV[2] if len(ARGV) > 2 else 'prev-C + MSE B + ω defined'
rule = V[RULE]
print(f'rule: {RULE}')
cols = K + ['resname', 'lab_resname', 'rama6', 'lab_cls', 'phi', 'lab_phi', 'psi', 'lab_psi', 'omega', 'lab_omega',
            'own_mc_max', 'lab_mcB', 'own_minocc', 'H_minocc', 'prevC_spread', 'd_prevC_N', 'd_C_nextN']
inlab = (E.lab_cls == CLS); inours = rule & (E.rama6 == CLS)
E['renumbered'] = E.pdb.isin(bad)

def why(r):
    if r['_merge'] == 'right_only': return 'not in our residue table'
    if not (r.d_prevC_N < 2.0 and r.d_C_nextN < 2.0): return 'not peptide-linked'
    if not (0 < r.own_mc_max < 30): return f'own mc B {r.own_mc_max}'
    if not r.own_minocc >= 1: return f'own heavy-atom occ {r.own_minocc}'
    if r.resname not in STD | {'MSE'}: return f'non-standard {r.resname}'
    w = abs(r.omega)
    if r.rama6 == 'CisPro' and not w < 30: return f'ω {r.omega}'
    if r.rama6 != 'CisPro' and not (w != w or w >= 150): return f'ω {r.omega}'
    if r.prevC_spread > 0 and r.resname != 'PRO': return f'prev C alternates differ by {r.prevC_spread} Å'
    if r.resname == 'MSE': return 'MSE all-heavy B >= 30'
    return f'passes; our class {r.rama6}'

pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
for title, sel, other in ((f'in the lab {CLS} list, not in ours', inlab & ~inours, 'rama6'),
                          (f'in our {CLS} set, not in the lab list', inours & ~inlab, 'lab_cls')):
    X = E[sel].copy()
    X['reason'] = [why(r) if title.startswith('in the lab') else
                   (f'lab class {r.lab_cls}' if r.lab_cls == r.lab_cls else 'lab dropped') for _, r in X.iterrows()]
    print(f'\n{title}: {len(X)}')
    if len(X): print(X[cols + ['renumbered', 'reason']].to_string(index=False))
print(f'\nlab {int(inlab.sum()):,}; ours {int(inours.sum()):,}; Δ {int(inours.sum()) - int(inlab.sum()):+d}')

# Pair lab-only with ours-only residues of the same chain by phi/psi (within 0.01°): same residue, different number.
lo, oo = E[inlab & ~inours], E[inours & ~inlab]
pairs = []
for (p, c), g in lo.groupby(['pdb', 'chain']):
    h = oo[(oo.pdb == p) & (oo.chain == c)]
    for r in g.itertuples():
        m = h[(np.abs(h.phi - r.lab_phi) <= 0.01) & (np.abs(h.psi - r.lab_psi) <= 0.01)]
        if len(m) == 1: pairs.append((p, c, r.resnum, int(m.resnum.iloc[0]), r.resnum - int(m.resnum.iloc[0])))
P = pd.DataFrame(pairs, columns=['pdb', 'chain', 'lab_resnum', 'our_resnum', 'offset'])
print(f'\nlab-only paired with an ours-only residue by φ/ψ (same residue, different number): {len(P)}')
if len(P): print(P.groupby(['pdb', 'chain']).agg(n=('offset', 'size'), offsets=('offset', lambda s: sorted(set(s)))).to_string())
up_l = lo[~lo.set_index(['pdb', 'chain', 'resnum']).index.isin(list(zip(P.pdb, P.chain, P.lab_resnum)))]
up_o = oo[~oo.set_index(['pdb', 'chain', 'resnum']).index.isin(list(zip(P.pdb, P.chain, P.our_resnum)))]
print(f'unpaired: lab-only {len(up_l)} {sorted(set(up_l.pdb))}; ours-only {len(up_o)} {sorted(set(up_o.pdb))}')
w = V['prev-C rule + MSE B'] & E.omega.isna()
print(f'\nresidues passing prev-C + MSE B with ω undefined (the lab required ω defined): {int(w.sum())};'
      f' in the lab lists {int((w & (E._merge != "left_only")).sum())}; by class', E[w].rama6.value_counts().to_dict())
