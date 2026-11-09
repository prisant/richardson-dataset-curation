#!/bin/bash
# HETATM mask fix (register B1) regression on msi.
#
# Input: ~/scratch/regress/<ds>/rama.jsonl, the unfiltered output of pydangle
# main 664e87a (run_rama_regress_msi.sh).  It is filtered twice with
# filter_pruned_residues.py: with the old masks (~/Desktop/Data/<dataset>,
# unchanged on msi, identical to the backup on nuc81) and with the new masks
# (~/scratch/b1/<dataset>: masks and synthetic pruned files rebuilt on nuc81
# by build_bfactor_masks.py / build_ersatz_pdbs.py --masks-only, copied here;
# lab pruned files symlinked).  compare_b1_filtered.py then checks the
# difference against the residues verify_hetatm_masks.py --added-tsv listed
# (~/scratch/b1/added/<ds>.tsv).
set -u
PY=$HOME/Desktop/Projects/pydangle-biopython/.venv/bin/python
FILTER=$HOME/Desktop/Projects/richardson-dataset-curation/pipeline/filter_pruned_residues.py
COMPARE=$HOME/scratch/b1/compare_b1_filtered.py
R=$HOME/scratch/regress
for spec in top100:top100_pdbs:pruned_bfilt top500:top500_pdbs:pruned_bfilt \
            top8000:top8000_pdbs_hom70:pruned_bfilt \
            top2018_full:top2018_pdbs_full_filtered_hom70:pruned_all \
            top2018_mc:top2018_pdbs_mc_filtered_hom70:pruned_mc; do
    IFS=: read ds dir suffix <<< "$spec"
    cd $R/$ds
    echo "== $ds $(date +%T)"
    "$PY" -I "$FILTER" rama.jsonl $HOME/Desktop/Data/$dir -s $suffix > b1_old_filtered.jsonl 2> b1_old_filter.log &
    "$PY" -I "$FILTER" rama.jsonl $HOME/scratch/b1/$dir -s $suffix > b1_new_filtered.jsonl 2> b1_new_filter.log &
    wait
    tail -n 2 b1_old_filter.log b1_new_filter.log
    "$PY" -I "$COMPARE" b1_old_filtered.jsonl b1_new_filtered.jsonl $HOME/scratch/b1/added/$ds.tsv > b1_compare.txt 2>&1
    cat b1_compare.txt
done
echo "b1 regression done $(date +%T)"
