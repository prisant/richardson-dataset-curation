#!/bin/bash
# Rama-precedence fix (register A2) regression on msi.  "Old" outputs are
# the continuity-regression runs of the 2.0 A code (identical to pydangle
# main 1a0e57e apart from a docstring): new2.jsonl for top100/500/8000,
# new.jsonl for top2018_full/mc.  This script runs only the fixed code:
# worktree ~/scratch/pydangle_rama = main 1a0e57e + rama.patch, via PYTHONPATH.
set -u
PY=$HOME/Desktop/Projects/pydangle-biopython/.venv/bin/python
RUN='from pydangle_biopython.main import main; import sys; sys.exit(main())'
FULL="phi; psi; omega; tau; chi1; chi2; chi3; chi4; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality"
MC="phi; psi; omega; tau; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality"
R=$HOME/scratch/regress
NEW_SRC=$HOME/scratch/pydangle_rama/src
for ds in top100 top500 top8000 top2018_full top2018_mc; do
    cmds=$FULL; [ $ds = top2018_mc ] && cmds=$MC
    cd $R/$ds
    echo "== $ds rama $(date +%T) files $(wc -l < filelist.txt)"
    PYTHONPATH=$NEW_SRC "$PY" -c "$RUN" -c "$cmds" -o jsonl -j 10 -f filelist.txt > rama.jsonl 2> rama.log || echo "$ds rama exit $?"
done
echo "rama done $(date +%T)"
