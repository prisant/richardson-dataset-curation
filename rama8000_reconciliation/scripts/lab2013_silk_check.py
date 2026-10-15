"""Silk on the lab's 2013 residue lists (data/Top8000_Rama_2013/residues) vs the archived rama8000 NDFT tables (data/*.ndft).

Options from the lab Makefile / run_silk.sh: -ndim=2 -label=1 -coords=2,3 -wrap -insep=, -gridsize=2 -cosine=5
-bounds=-180,180,-180,180 -fraction; -twopass=6.5 General, 8 otherwise; Gly from gly-sym.csv. Run from rama8000_reconciliation/.
"""
import sys, os, subprocess, tempfile
sys.path.insert(0, 'scripts')
import numpy as np, ndft
for csv, nd, tp in (('general-noGPIVpreP', 'general', '6.5'), ('ileval-nopreP', 'ileval', '8'), ('gly-sym', 'glycine', '8'),
                    ('prepro-noGP', 'prepro', '8'), ('transpro', 'transpro', '8'), ('cispro', 'cispro', '8')):
    fd, fn = tempfile.mkstemp(suffix='.ndft'); os.close(fd)
    cmd = ['java', '-Xmx8g', '-cp', 'silk_build/classes', 'silk.SilkCmdLine', '-ndim=2', '-label=1', '-coords=2,3', '-wrap',
           '-insep=,', '-gridsize=2', '-cosine=5', f'-twopass={tp}', '-bounds=-180,180,-180,180', '-fraction', '-ndft', '-quiet']
    with open(f'data/Top8000_Rama_2013/residues/rama8000-{csv}.csv', 'rb') as f, open(fn, 'wb') as g:
        subprocess.run(cmd, stdin=f, stdout=g, stderr=subprocess.DEVNULL, check=True)
    n, G = ndft.load(fn); rc, T = ndft.load(f'data/{nd}.ndft'); os.remove(fn)
    print(f'{csv:20s} realcount {n} vs {rc}; max |Δ| {np.abs(G - T).max():.3g}; identical bins {int((G == T).sum())}/{T.size}', flush=True)
