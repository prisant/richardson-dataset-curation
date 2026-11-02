"""Lay out the lab's 2011 Top8000 chain files in the layout of our top8000 tree, so the record scripts can run on them.

Input: top8000_chains_70.tgz (Richardson lab, 7,957 single-chain files <pdb>FH_<chain>.pdb, Reduce 3.15.091106 with
flips; server date 2012-05-29; sha256 4a6bb1b9bb340cd42c14cd9054688732ac09317dfe0af49f3fd453f77ba486e6), extracted.
Output under DST: <pdb[:2]>/<pdb>/<pdb>FH.pdb (ATOM/HETATM of model 1, all listed chains of the entry concatenated)
and <pdb>_ersatz.pdb (the same without hydrogens), plus Top8000-best_hom70_pdb_chain_list.csv.
Then run the record scripts with T8000_ROOT=DST.
usage: lab2011_tree.py SRC_DIR DST_DIR
"""
import sys, os
from collections import defaultdict

def is_h(line):
    el = line[76:78].strip()
    if el: return el in ('H', 'D')
    return line[12:16].strip().lstrip('0123456789')[:1] in ('H', 'D')

src, dst = sys.argv[1], sys.argv[2]
files = defaultdict(list)
for f in sorted(os.listdir(src)):
    p, c = f[:-4].split('FH_'); files[p.lower()].append((c, f))
with open(f'{dst}/Top8000-best_hom70_pdb_chain_list.csv', 'w') as cl:
    cl.write('pdb_id,chain\n')
    for p, cs in sorted(files.items()):
        os.makedirs(f'{dst}/{p[:2]}/{p}', exist_ok=True)
        fh, er = open(f'{dst}/{p[:2]}/{p}/{p}FH.pdb', 'w'), open(f'{dst}/{p[:2]}/{p}/{p}_ersatz.pdb', 'w')
        for c, f in cs:
            cl.write(f'{p},{c}\n')
            for line in open(f'{src}/{f}'):
                if line.startswith('ENDMDL'): break
                if not line.startswith(('ATOM', 'HETATM')): continue
                fh.write(line)
                if not is_h(line): er.write(line)
        for f in (fh, er): f.write('END\n'); f.close()
print(f'{sum(map(len, files.values()))} chain files, {len(files)} entries')
