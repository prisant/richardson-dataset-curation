"""Apply a residue rule to one Rama class, run Silk at fixed params, compare to the archived NDFT bin by bin.

usage: class_rule_check.py RECORDS.jsonl CLASS NDFT twopass
CLASS in {cispro, transpro}.  Prints counts under rule variants and grid agreement.
"""
import sys, json, subprocess, tempfile, os, numpy as np
sys.path.insert(0, 'scripts'); import ndft
recs, cls, ndft_path, twopass = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
R = [json.loads(l) for l in open(recs)]
link = lambda r: r['d_prevC_N'] < 2.0 and r['d_C_nextN'] < 2.0
in_cls = {'cispro': lambda r: abs(r['omega']) < 30, 'transpro': lambda r: abs(r['omega']) >= 30}[cls]
P = [r for r in R if in_cls(r) and link(r)]
rc, T = ndft.load(ndft_path)

def silk(rows):
    fd, fn = tempfile.mkstemp(suffix='.ndft'); os.close(fd)
    txt = ''.join(f"{r['pdb']}_{r['chain']}{r['resnum']}{r['ins'].strip()} {r['phi']} {r['psi']}\n" for r in rows)
    cmd = ['java','-cp','silk_build/classes','silk.SilkCmdLine','-ndim=2','-coords=2,3','-label=1','-bounds=-180,180,-180,180',
           '-wrap','-gridsize=2','-cosine=5',f'-twopass={twopass}','-fraction','-ndft','-quiet']
    open(fn,'wb').write(subprocess.run(cmd, input=txt.encode(), capture_output=True, check=True).stdout)
    _, g = ndft.load(fn); os.remove(fn); return g

rules = {
  'published-style (tri mainchain B<=30)':        lambda r: max(r['prev_mc'], r['own_mc'], r['next_mc']) <= 30,
  'own mc B<30':                                  lambda r: r['own_mc'] < 30,
  'own mc B<30, no altloc':                       lambda r: r['own_mc'] < 30 and not r['own_alt'],
  'own mc B<30, all occ=1 (CisPro rule)':         lambda r: r['own_mc'] < 30 and r['own_minocc'] >= 1,
}
print(f'{cls}: realcount {rc}; candidates with phi,psi: {len(P)}')
print(f'{"rule":42s}{"n":>8s}{"Δn":>8s}{"meanΔ×1e4":>11s}{"supp":>6s}{"rank-identical bins":>21s}{"max|Δrank|":>11s}')
for name, f in rules.items():
    rows = [r for r in P if f(r)]
    g = silk(rows); d = np.rint(g * rc) - np.rint(T * rc)
    print(f'{name:42s}{len(rows):8d}{len(rows)-rc:+8d}{np.abs(g-T).mean()*1e4:11.3f}{int(((g>0)!=(T>0)).sum()):6d}'
          f'{int((d==0).sum()):>15d}/{d.size}{int(np.abs(d).max()):11d}')
