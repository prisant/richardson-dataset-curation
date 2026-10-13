"""Final CisPro set: own mainchain B<30, no altloc, all own atoms occupancy 1. Compare grid bin-by-bin."""
import sys, json, subprocess, numpy as np
sys.path.insert(0, 'scripts'); import ndft
R = [json.loads(l) for l in open('data/cis_records_ext.jsonl')]
link = lambda r: r['d_prevC_N'] < 2.0 and r['d_C_nextN'] < 2.0
S = [r for r in R if abs(r['omega']) < 30 and link(r) and r['own_mc'] < 30 and not r['own_alt']]
F = [r for r in S if r['own_minocc'] >= 1]
print('the 9 partial-occupancy residues removed:')
for r in S:
    if r['own_minocc'] < 1:
        print(f"  {r['pdb']} {r['chain']}{r['resnum']}{r['ins'].strip():1s} phi={r['phi']:7.1f} psi={r['psi']:7.1f} omega={r['omega']:6.1f} "
              f"min occ={r['own_minocc']:.2f} ownB={r['own_mc']:.1f}")
txt = ''.join(f"{r['pdb']}_{r['chain']}{r['resnum']}{r['ins'].strip()} {r['phi']} {r['psi']}\n"
              for r in sorted(F, key=lambda r: (r['pdb'], r['chain'], r['resnum'], r['ins'])))
open('results/cispro_final_set.txt', 'w').write(txt)
cmd = ['java','-cp','silk_build/classes','silk.SilkCmdLine','-ndim=2','-coords=2,3','-label=1','-bounds=-180,180,-180,180',
       '-wrap','-gridsize=2','-cosine=5','-twopass=8','-fraction','-ndft','-quiet','-title=cis Pro, 70%, B<30']
open('results/cispro_final_c5_t8.ndft','wb').write(subprocess.run(cmd, input=txt.encode(), capture_output=True, check=True).stdout)
rc, T = ndft.load('data/cispro.ndft'); _, G = ndft.load('results/cispro_final_c5_t8.ndft')
a, b = np.rint(G * rc).astype(int), np.rint(T * rc).astype(int)
d = a - b
print(f'\nn={len(F)} (NDFT realcount {rc})')
print(f'bins with identical rank: {(d==0).sum()} / {d.size}   |Δrank|=1: {(abs(d)==1).sum()}   |Δrank|>1: {(abs(d)>1).sum()}   max |Δrank| {abs(d).max()}')
print(f'max |Δ float value| {np.abs(G-T).max():.3e}')
for i, j in zip(*np.nonzero(abs(d) > 1)):
    print(f'  bin phi={ndft.center(i)} psi={ndft.center(j)} ours={a[i,j]} ref={b[i,j]}')
