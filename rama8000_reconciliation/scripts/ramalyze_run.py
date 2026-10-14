"""cctbx mmtbx.validation.ramalyze on every top8000 structure (model files: *FH.pdb or *_ersatz.pdb).

Runs in the cctbx environment (micromamba env 'cctbx', conda-forge cctbx-base 2026.9):
  ~/.local/bin/micromamba run -n cctbx python scripts/ramalyze_run.py FH|ersatz data/ramalyze_<kind>.tsv
Writes one row per ramalyze result (all residues, not only outliers):
  pdb model chain resseq icode altloc resname phi psi res_type
res_type: 0 General, 1 Gly, 2 cis-Pro, 3 trans-Pro, 4 pre-Pro, 5 Ile/Val (ramalyze's own assignment).
ANISOU records are dropped before parsing (Reduce output can leave an ANISOU for a hydrogen with no matching
ATOM, which iotbx rejects; ANISOU does not affect phi/psi). Errors are written as '#ERR pdb message' rows.
"""
import sys, os, glob
from multiprocessing import Pool
import iotbx.pdb
from mmtbx.validation import ramalyze

ROOT = os.path.expanduser('~/Desktop/Data/top8000_pdbs_hom70')
SUF = {'FH': 'FH.pdb', 'ersatz': '_ersatz.pdb'}

def work(fn):
    pdb = os.path.basename(fn)[:4]
    try:
        lines = [l.rstrip('\n') for l in open(fn) if not l.startswith('ANISOU')]
        h = iotbx.pdb.input(source_info=fn, lines=lines).construct_hierarchy()
        r = ramalyze.ramalyze(pdb_hierarchy=h, outliers_only=False, out=open(os.devnull, 'w'), quiet=True)
        return [f'{pdb}\t{x.model_id.strip()}\t{x.chain_id.strip()}\t{x.resseq.strip()}\t{x.icode.strip()}\t{x.altloc.strip()}\t'
                f'{x.resname}\t{x.phi:.3f}\t{x.psi:.3f}\t{x.res_type}' for x in r.results]
    except Exception as e:
        return [f'#ERR\t{pdb}\t{type(e).__name__}: {str(e).splitlines()[0] if str(e) else ""}']

if __name__ == '__main__':
    kind, out = sys.argv[1], sys.argv[2]
    files = sorted(glob.glob(f'{ROOT}/*/*/*{SUF[kind]}'))
    with Pool(8) as p, open(out, 'w') as f:
        f.write('pdb\tmodel\tchain\tresseq\ticode\taltloc\tresname\tphi\tpsi\tres_type\n')
        for i, rows in enumerate(p.imap_unordered(work, files, chunksize=4)):
            for line in rows: f.write(line + '\n')
    print(f'{kind}: {len(files)} files -> {out}')
