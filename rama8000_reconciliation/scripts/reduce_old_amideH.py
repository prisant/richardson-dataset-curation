"""Run Reduce 3.3 (~/bin/reduce_old -FLIP) on the original downloaded PDBs and record amide-H attributes.

Output TSV has the same columns as amideH_records.tsv. Reduce output is parsed from stdout; no files written.
usage: reduce_old_amideH.py OUT.tsv [pdb ...]   (default: all pdbs in the chain list)
"""
import sys, os, subprocess
from collections import defaultdict
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from recount3 import chains, ROOT
REDUCE = os.path.expanduser('~/bin/reduce_old')

def work(pdb):
    fn = f'{ROOT}/{pdb[:2]}/{pdb}/{pdb}.pdb'
    try:
        txt = subprocess.run([REDUCE, '-FLIP', fn], capture_output=True, timeout=1800).stdout.decode('latin-1')
    except subprocess.TimeoutExpired:
        return [f'#TIMEOUT {pdb}']
    seen, H = {}, defaultdict(list)
    for line in txt.splitlines():
        if line.startswith('ENDMDL'): break
        if not line.startswith(('ATOM', 'HETATM')): continue
        ch = line[21]
        if ch not in chains[pdb]: continue
        key = (ch, int(line[22:26]), line[26])
        seen[key] = 1
        if line[12:16].strip() == 'H':
            H[key].append((float(line[60:66]), float(line[54:60])))
    if not seen: return [f'#EMPTY {pdb}']
    return [f"{pdb}\t{k[0]}\t{k[1]}\t{k[2]}\t{int(bool(H.get(k)))}\t{max(b for b, o in H[k]) if H.get(k) else 'nan'}\t"
            f"{min(o for b, o in H[k]) if H.get(k) else 'nan'}" for k in seen]

if __name__ == '__main__':
    pdbs = sys.argv[2:] or sorted(chains)
    with Pool(8) as p, open(sys.argv[1], 'w') as f:
        f.write('pdb\tchain\tresnum\tins\thas_H\tH_maxB\tH_minocc\n')
        for recs in p.imap_unordered(work, pdbs, chunksize=4):
            for line in recs:
                (sys.stderr.write(line + '\n') if line.startswith('#') else f.write(line + '\n'))
