#!/bin/bash
# Rerun only the fixed pydangle (continuity-fix 5a22669 + 2.0 A bond cutoff
# working tree) on top100, top500, top8000; old.jsonl and java.txt from the
# first regression are unchanged and reused.  Output: new2.jsonl beside them.
set -u
PY=$HOME/Desktop/Projects/pydangle-biopython/.venv/bin/python
RUN='from pydangle_biopython.main import main; import sys; sys.exit(main())'
CMDS="phi; psi; omega; tau; chi1; chi2; chi3; chi4; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality"
D=$HOME/Desktop/Data
R=$HOME/scratch/regress
mkdir -p $R/top100
find $D/top100_pdbs -name "*_ersatz.pdb" | sort > $R/top100/filelist.txt
for ds in top100 top500 top8000; do
    cd $R/$ds
    echo "== $ds new2 $(date +%T) files $(wc -l < filelist.txt)"
    "$PY" -c "$RUN" -c "$CMDS" -o jsonl -j 4 -f filelist.txt > new2.jsonl 2> new2.log || echo "$ds new2 exit $?"
done
echo "new2 done $(date +%T)"
