"""Grid check for non-Pro classes: reconciled rule with / without the amide-H occupancy exclusion.

usage: nonpro_grid_check.py CLASS NDFT twopass
"""
import sys, subprocess, tempfile, os, numpy as np, pandas as pd
sys.path.insert(0, 'scripts'); import ndft
cls, path, tp = sys.argv[1], sys.argv[2], sys.argv[3]
d = pd.read_csv('data/residue_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
h = pd.read_csv('data/amideH_records.tsv', sep='\t', keep_default_na=False, na_values=['nan'], dtype={'ins': str, 'chain': str})
d = d.merge(h, on=['pdb', 'chain', 'resnum', 'ins'], how='left')
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
rule = (d.d_prevC_N < 2.0) & (d.d_C_nextN < 2.0) & (d.own_mc_max < 30) & (d.own_minocc >= 1) & (d.own_mc_max > 0) \
       & d.resname.isin(STD) & (d.rama6 == cls)
rc, T = ndft.load(path)
def silk(df):
    fd, fn = tempfile.mkstemp(suffix='.ndft'); os.close(fd)
    txt = '\n'.join(f'x {p} {q}' for p, q in zip(df.phi, df.psi)) + '\n'
    cmd = ['java', '-Xmx8g', '-cp', 'silk_build/classes', 'silk.SilkCmdLine', '-ndim=2', '-coords=2,3', '-bounds=-180,180,-180,180',
           '-wrap', '-gridsize=2', '-cosine=5', f'-twopass={tp}', '-fraction', '-ndft', '-quiet'] + (['-sym'] if cls == 'Gly' else [])
    open(fn, 'wb').write(subprocess.run(cmd, input=txt.encode(), capture_output=True, check=True).stdout)
    n, g = ndft.load(fn); os.remove(fn); return n, g
for name, m in [('rule', rule), ('rule + drop amide-H occ<1', rule & ~(d.H_minocc < 1))]:
    n, g = silk(d[m]); dr = np.rint(g * rc) - np.rint(T * rc)
    print(f'{cls} {name:28s} n={int(m.sum())} (silk realcount {n}) Δ={n-rc:+d} meanΔ×1e4={np.abs(g-T).mean()*1e4:.4f} '
          f'supp={int(((g>0)!=(T>0)).sum())} identical={int((dr==0).sum())} max|Δrank|={int(np.abs(dr).max())}')
