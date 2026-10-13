"""Non-Pro grid check of the B-factor boundary: own mainchain max B < 30 vs <= 30, full reconciled rule (incl. amide-H occ 1).

Silk params per silk/old/Makefile: -mask=10,13 (General) -> -twopass=6.5; -mask=10,16 (Gly, PrePro) -> -twopass=8.
IleVal has no recipe there (class is new in rama8000), so it is run at both.
Gly symmetry: SilkCmdLine has no -sym flag; Rota2D's -sym appends (-phi, -psi) for every point, done here explicitly.
usage: b30_grid_check.py [CLASS ...]   (default: all). Run from rama8000_reconciliation/.
"""
import sys, subprocess, tempfile, os, numpy as np, pandas as pd
sys.path.insert(0, 'scripts'); import ndft
d = pd.read_csv('data/residue_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
h = pd.read_csv('data/amideH_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
d = d.merge(h, on=['pdb', 'chain', 'resnum', 'ins'], how='left')
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
common = ((d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0) & (d.own_minocc >= 1) & (d.own_mc_max > 0) & d.resname.isin(STD)
          & ~(d.H_minocc < 1))

def silk(df, tp, sym):
    fd, fn = tempfile.mkstemp(suffix='.ndft'); os.close(fd)
    pts = list(zip(df.phi, df.psi)) + ([(-p, -q) for p, q in zip(df.phi, df.psi)] if sym else [])
    txt = '\n'.join(f'x {p} {q}' for p, q in pts) + '\n'
    cmd = ['java', '-Xmx8g', '-cp', 'silk_build/classes', 'silk.SilkCmdLine', '-ndim=2', '-coords=2,3', '-bounds=-180,180,-180,180',
           '-wrap', '-gridsize=2', '-cosine=5', f'-twopass={tp}', '-fraction', '-ndft', '-quiet']
    open(fn, 'wb').write(subprocess.run(cmd, input=txt.encode(), capture_output=True, check=True).stdout)
    n, g = ndft.load(fn); os.remove(fn); return n, g

jobs = [('General', 'general', '6.5'), ('PrePro', 'prepro', '8'), ('Gly', 'glycine', '8'),
        ('IleVal', 'ileval', '6.5'), ('IleVal', 'ileval', '8')]
print(f'{"class":8s}{"twopass":>8s}  {"boundary":9s}{"n":>9s}{"Δn":>7s}{"meanΔ×1e4":>11s}{"supp":>6s}{"rank-identical":>16s}{"max|Δrank|":>11s}')
for cls, f, tp in [j for j in jobs if len(sys.argv) < 2 or j[0] in sys.argv[1:]]:
    rc, T = ndft.load(f'data/{f}.ndft')
    for name, m in (('B < 30', d.own_mc_max < 30), ('B <= 30', d.own_mc_max <= 30)):
        sel = common & m & (d.rama6 == cls)
        n, g = silk(d[sel], tp, cls == 'Gly'); dr = np.rint(g * rc) - np.rint(T * rc)
        print(f'{cls:8s}{tp:>8s}  {name:9s}{n:9d}{n-rc:+7d}{np.abs(g-T).mean()*1e4:11.4f}{int(((g>0)!=(T>0)).sum()):6d}'
              f'{int((dr==0).sum()):>10d}/{dr.size}{int(np.abs(dr).max()):11d}', flush=True)
