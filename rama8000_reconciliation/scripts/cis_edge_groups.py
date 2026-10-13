"""Score edge-rule groups (remove from / add to set S) against cispro.ndft at Makefile params."""
import sys, json, subprocess, tempfile, os, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, 'scripts'); import ndft
_, T = ndft.load('data/cispro.ndft')
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
R = [json.loads(l) for l in open('data/cis_records_ext.jsonl')]
link = lambda r, c=2.0: r['d_prevC_N'] < c and r['d_C_nextN'] < c
S = [r for r in R if abs(r['omega']) < 30 and link(r) and r['own_mc'] < 30 and not r['own_alt']]
Sid = {id(r) for r in S}
def score(rows):
    fd, fn = tempfile.mkstemp(suffix='.ndft'); os.close(fd)
    txt = ''.join(f"{r['pdb']}_{r['chain']}{r['resnum']} {r['phi']} {r['psi']}\n" for r in rows)
    cmd = ['java','-cp','silk_build/classes','silk.SilkCmdLine','-ndim=2','-coords=2,3','-label=1','-bounds=-180,180,-180,180',
           '-wrap','-gridsize=2','-cosine=5','-twopass=8','-fraction','-ndft','-quiet']
    open(fn,'wb').write(subprocess.run(cmd, input=txt.encode(), capture_output=True, check=True).stdout)
    _, g = ndft.load(fn); os.remove(fn)
    return float(np.abs(g - T).mean() * 1e4), int(((g > 0) != (T > 0)).sum())
rem = {
  'prev/next mainchain altloc':      lambda r: r['prev_mcalt'] or r['next_mcalt'],
  'prev/next any altloc':            lambda r: r['prev_alt'] or r['next_alt'],
  'own partial occupancy (<1)':      lambda r: r['own_minocc'] < 1,
  'neighbour partial occupancy':     lambda r: r['nbr_minocc'] < 1,
  'own mc B in [29,30)':             lambda r: 29 <= r['own_mc'] < 30,
  '|omega| in [25,30)':              lambda r: abs(r['omega']) >= 25,
  'own O missing':                   lambda r: not r['has_O'],
  'own any atom B >= 30':            lambda r: r['own_all'] >= 30,
  'own CB B >= 30':                  lambda r: r['own_cb'] >= 30,
  'prev or next non-standard':       lambda r: r['resname_prev'] not in STD or r['resname_next'] not in STD,
  'peptide C-N > 1.45 A':            lambda r: max(r['d_prevC_N'], r['d_C_nextN']) > 1.45,
}
add = {
  '|omega| in [30,35), rule ok':     lambda r: 30 <= abs(r['omega']) < 35 and link(r) and r['own_mc'] < 30 and not r['own_alt'],
  'C-N in [2.0,2.5), rule ok':       lambda r: abs(r['omega']) < 30 and not link(r) and link(r, 2.5) and r['own_mc'] < 30 and not r['own_alt'],
  'own mc B in [30,30.5), no alt':   lambda r: abs(r['omega']) < 30 and link(r) and 30 <= r['own_mc'] < 30.5 and not r['own_alt'],
  'altloc only on sidechain, B ok':  lambda r: abs(r['omega']) < 30 and link(r) and r['own_mc'] < 30 and r['own_alt'] and not r['own_mcalt'],
}
jobs = [('baseline S', S)]
for k, f in rem.items(): jobs.append(('remove: ' + k, [r for r in S if not f(r)]))
for k, f in add.items(): jobs.append(('add: ' + k, S + [r for r in R if id(r) not in Sid and f(r)]))
with ThreadPoolExecutor(8) as ex: sc = list(ex.map(lambda j: score(j[1]), jobs))
print(f'{"group":42s}{"Δn":>6s}{"n":>6s}{"meanΔ×1e4":>11s}{"supp":>6s}')
for (k, rows), (s, m) in zip(jobs, sc):
    print(f'{k:42s}{len(rows)-len(S):>+6d}{len(rows):6d}{s:11.3f}{m:6d}')
