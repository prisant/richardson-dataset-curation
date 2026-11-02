"""Did the lab's 2013 residue lists keep residues with altloc records? (2011 chain files, data/lab2011/ tables.)
own_alt: the residue's own N, CA, C or O carries an altloc label; own_minocc: lowest heavy-atom occupancy of the residue;
prevC_alt: the preceding residue's C carries an altloc label.
Usage (from rama8000_reconciliation/):  python scripts/lab_altloc_use.py > results/lab_altloc_use.txt
"""
import pandas as pd
src = open('scripts/current_files_reconcile.py').read()
exec(src[:src.index('T = {}')])          # K, L (lab lists), load
t = load('data/lab2011'); t['resnum'] = t.resnum.astype(int)
t = t.merge(L.assign(resnum=L.resnum.astype(int))[K + ['lab_cls']], on=K, how='left')
t['listed'] = t.lab_cls.notna()
for nm, m in (('own mainchain altloc label', t.own_alt > 0),
              ('own mainchain altloc label, heavy-atom occupancy 1', (t.own_alt > 0) & (t.own_minocc >= 1)),
              ('own heavy-atom occupancy < 1', t.own_minocc < 1),
              ('preceding C altloc label', t.prevC_alt > 0),
              ('no altloc label, occupancy 1 (own or preceding C)', (t.own_alt == 0) & (t.own_minocc >= 1) & ~(t.prevC_alt > 0))):
    print(f'{nm:52s} residues {int(m.sum()):9,d}; in lab lists {int((m & t.listed).sum()):9,d}')
print(f'lab-listed residues {int(t.listed.sum()):,}')
