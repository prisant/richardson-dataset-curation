"""Locate reference-only support (lab has density, our best-fit set does not) and match excluded candidates."""
import sys, json, subprocess, numpy as np
from scipy import ndimage
sys.path.insert(0, 'scripts'); import ndft
pts = open('data/cispro_bestfit.txt').read()
cmd = ['java','-cp','silk_build/classes','silk.SilkCmdLine','-ndim=2','-coords=2,3','-label=1','-bounds=-180,180,-180,180',
       '-wrap','-gridsize=2','-cosine=5','-twopass=8','-fraction','-ndft','-quiet']
open('results/cispro_bestfit_c5_t8.ndft','wb').write(subprocess.run(cmd,input=pts.encode(),capture_output=True,check=True).stdout)
_, g = ndft.load('results/cispro_bestfit_c5_t8.ndft'); _, t = ndft.load('data/cispro.ndft')
R = [json.loads(l) for l in open('data/pro_records.jsonl')]
cand = [r for r in R if abs(r['omega']) < 30 and not (max(r['own_mc'], r['prevC'], r['nextN']) < 30 and not r['alt'])]
for mask, name in [((t > 0) & (g == 0), 'REF-ONLY (lab has density, we do not)'), ((g > 0) & (t == 0), 'OURS-ONLY')]:
    lab, n = ndimage.label(mask, structure=np.ones((3, 3)))
    print(f'\n{name}: {mask.sum()} bins in {n} clusters')
    for k in range(1, n + 1):
        ii, jj = np.nonzero(lab == k)
        phis, psis = [ndft.center(i) for i in ii], [ndft.center(j) for j in jj]
        print(f'  cluster {k}: {len(ii):3d} bins  phi {min(phis)}..{max(phis)}  psi {min(psis)}..{max(psis)}  max ref val {t[lab==k].max()*3654:.0f}/3654')
        lo_phi, hi_phi, lo_psi, hi_psi = min(phis) - 15, max(phis) + 15, min(psis) - 15, max(psis) + 15
        for r in cand:
            if lo_phi <= r['phi'] <= hi_phi and lo_psi <= r['psi'] <= hi_psi:
                why = []
                dB = max(r['own_mc'], r['prevC'], r['nextN'])
                if dB >= 30: why.append(f'maxB={dB:.1f}')
                if r['alt']: why.append(f"alt={''.join(r['alt'])}{'(mc)' if r['mcalt'] else '(sc)'}")
                print(f"      candidate {r['pdb']} {r['chain']}{r['resnum']}{r['ins'].strip():1s} phi={r['phi']:7.1f} psi={r['psi']:7.1f} omega={r['omega']:6.1f}  {' '.join(why)}")
