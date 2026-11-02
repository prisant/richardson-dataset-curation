"""Were the lab's 2013 residue lists computed from the 2011 chain files (top8000_chains_70.tgz)?

For lab-listed residues whose φ, ψ or mainchain max B differ between the 2011 files and our current downloads, count
which file set the lab's printed values match (φ/ψ within 0.01°; B within 0.05 of the one-decimal label).
MSE is left out of the B test (the lab's MSE B also includes SE; mse_bfactor.py).
Usage (from rama8000_reconciliation/):  python scripts/lab_file_provenance.py > results/lab_file_provenance.txt
"""
import numpy as np, pandas as pd
exec(open('scripts/current_files_reconcile.py').read().split('\ndef load')[0])   # K, LBL, FILES, LAB

rows = []
for cls, f in FILES.items():
    for ln in open(f'{LAB}/residues/rama8000-{f}.csv'):
        lab, phi, psi, om, tau = ln.rstrip('\n').rsplit(',', 4)
        m = LBL.match(lab)
        rows.append((m.group(1), m.group(2), int(m.group(3)), m.group(4).strip(), m.group(5), float(m.group(6)), float(phi), float(psi)))
L = pd.DataFrame(rows, columns=K + ['lab_resname', 'lab_B', 'lab_phi', 'lab_psi'])
rd = lambda d: pd.read_csv(f'{d}/residue_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'],
                           dtype={'ins': str, 'chain': str, 'pdb': str}).assign(ins=lambda t: t.ins.str.strip())[K + ['phi', 'psi', 'own_mc_max']]
M = L.merge(rd('data/lab2011'), on=K, how='left').merge(rd('data'), on=K, how='left', suffixes=('_11', '_cur'))
dang = lambda a, b: np.abs((a - b + 180) % 360 - 180)
ok = lambda s: (dang(M.lab_phi, M['phi' + s]) <= 0.01) & (dang(M.lab_psi, M['psi' + s]) <= 0.01)
okB = lambda s: (M.lab_B - M['own_mc_max' + s]).abs() <= 0.05 + 1e-9
print(f'lab-listed residues {len(M):,}; found in 2011 files {int(M.phi_11.notna().sum()):,}; in current files {int(M.phi_cur.notna().sum()):,}')
for nm, a, b in (('φ/ψ', ok('_11'), ok('_cur')),):
    print(f'{nm} match lab: 2011 files {int(a.sum()):,}; current files {int(b.sum()):,}')
    d = (dang(M.phi_11, M.phi_cur) > 0.01) | (dang(M.psi_11, M.psi_cur) > 0.01)
    print(f'  where the file sets differ ({int(d.sum()):,}): match 2011 only {int((d & a & ~b).sum()):,};'
          f' current only {int((d & b & ~a).sum()):,}; neither {int((d & ~a & ~b).sum()):,}')
n = M.lab_resname != 'MSE'
a, b = okB('_11') & n, okB('_cur') & n
d = n & ((M.own_mc_max_11 - M.own_mc_max_cur).abs() > 0.05)
print(f'mainchain max B (non-MSE, {int(n.sum()):,}) match lab: 2011 files {int(a.sum()):,}; current files {int(b.sum()):,}')
print(f'  where the file sets differ by > 0.05 ({int(d.sum()):,}): match 2011 only {int((d & a & ~b).sum()):,};'
      f' current only {int((d & b & ~a).sum()):,}; neither {int((d & ~a & ~b).sum()):,}')
only11 = M.phi_11.notna() & M.phi_cur.isna(); onlycur = M.phi_11.isna() & M.phi_cur.notna()
print(f'lab residue key present only in the 2011 files {int(only11.sum()):,}; only in current files {int(onlycur.sum()):,};'
      f' in neither {int((M.phi_11.isna() & M.phi_cur.isna()).sum()):,}')
