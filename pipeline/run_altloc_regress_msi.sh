#!/bin/bash
# Altloc selection regression on msi.
#
# Runs pydangle branch altloc-first (copied to ~/scratch/pydangle_alt/src,
# used via PYTHONPATH) on the same file lists as run_rama_regress_msi.sh, in
# both --altloc modes.  "Old" is rama.jsonl, the output of pydangle main
# 664e87a.  Checks:
#   1. --altloc occupancy reproduces rama.jsonl byte for byte (meta line aside);
#   2. --altloc first vs rama.jsonl, all residues (compare_altloc.py);
#   3. the same after the dataset filter (masks with polymer HETATM
#      residues, --residues std20+MSE):
#      alt_first_std.jsonl vs b2_std.jsonl (run_b2_regress_msi.sh).
# Scripts used from ~/scratch/b2 and ~/scratch/b6.
#
# Environment: MODES (default "occupancy first") selects the runs; TAG
# (default empty) is appended to the output names (alt_<mode>TAG.jsonl,
# ...) so a rerun keeps the earlier outputs.  Check 1 runs only when
# occupancy is in MODES.
set -u
MODES=${MODES:-occupancy first}
TAG=${TAG:-}
PY=$HOME/Desktop/Projects/pydangle-biopython/.venv/bin/python
RUN='from pydangle_biopython.main import main; import sys; sys.exit(main())'
FULL="phi; psi; omega; tau; chi1; chi2; chi3; chi4; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality"
MC="phi; psi; omega; tau; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality"
R=$HOME/scratch/regress
SRC=${SRC:-$HOME/scratch/pydangle_alt/src}  # pydangle source to test
for spec in top100:top100_pdbs:pruned_bfilt top500:top500_pdbs:pruned_bfilt \
            top8000:top8000_pdbs_hom70:pruned_bfilt \
            top2018_full:top2018_pdbs_full_filtered_hom70:pruned_all \
            top2018_mc:top2018_pdbs_mc_filtered_hom70:pruned_mc; do
    IFS=: read ds dir suffix <<< "$spec"
    cmds=$FULL; [ $ds = top2018_mc ] && cmds=$MC
    cd $R/$ds
    echo "== $ds $(date +%T) files $(wc -l < filelist.txt)"
    for mode in $MODES; do
        PYTHONPATH=$SRC "$PY" -W ignore -c "$RUN" --altloc $mode -c "$cmds" -o jsonl -j 10 \
            -f filelist.txt > alt_$mode$TAG.jsonl 2> alt_$mode$TAG.log || echo "  $mode exit $?"
    done
    if [[ " $MODES " == *" occupancy "* ]]; then
        if cmp -s <(tail -n +2 alt_occupancy$TAG.jsonl) <(tail -n +2 rama.jsonl); then
            echo "  occupancy == rama.jsonl (byte-identical after the meta line)"
        else
            echo "  occupancy DIFFERS from rama.jsonl"
        fi
    fi
    [[ " $MODES " == *" first "* ]] || continue
    "$PY" -I $HOME/scratch/b6/compare_altloc.py rama.jsonl alt_first$TAG.jsonl > alt_compare_all$TAG.txt 2>&1
    "$PY" -I $HOME/scratch/b2/filter_pruned_residues.py alt_first$TAG.jsonl $HOME/Desktop/Data/$dir \
        -s $suffix --residues std20+MSE > alt_first_std$TAG.jsonl 2> alt_first_std$TAG.log
    "$PY" -I $HOME/scratch/b6/compare_altloc.py b2_std.jsonl alt_first_std$TAG.jsonl > alt_compare_std$TAG.txt 2>&1
    echo "  all residues:";  sed -n '1,3p;/common records/,$p' alt_compare_all$TAG.txt | head -30 | sed 's/^/    /'
    echo "  filtered (std20+MSE):"; sed -n '1,3p;/common records/,$p' alt_compare_std$TAG.txt | head -30 | sed 's/^/    /'
done
echo "altloc regression done $(date +%T)"
