#!/bin/bash
# Residue-set (--residues std20+MSE) regression on msi.
#
# Input: ~/scratch/regress/<ds>/rama.jsonl, unfiltered output of pydangle main
# 664e87a (run_rama_regress_msi.sh).  Filtered with the HETATM-fixed masks in
# ~/Desktop/Data twice, --residues all and --residues std20+MSE; the first must
# equal b1_new_filtered.jsonl (run_b1_regress_msi.sh) byte for byte, and
# check_residue_set.py checks that the second only drops records outside the set.
# filter_pruned_residues.py and check_residue_set.py are copied to ~/scratch/b2.
set -u
PY=$HOME/Desktop/Projects/pydangle-biopython/.venv/bin/python
B2=$HOME/scratch/b2
for spec in top100:top100_pdbs:pruned_bfilt top500:top500_pdbs:pruned_bfilt \
            top8000:top8000_pdbs_hom70:pruned_bfilt \
            top2018_full:top2018_pdbs_full_filtered_hom70:pruned_all \
            top2018_mc:top2018_pdbs_mc_filtered_hom70:pruned_mc; do
    IFS=: read ds dir suffix <<< "$spec"
    cd $HOME/scratch/regress/$ds
    echo "== $ds $(date +%T)"
    "$PY" -I $B2/filter_pruned_residues.py rama.jsonl $HOME/Desktop/Data/$dir -s $suffix --residues all > b2_all.jsonl 2> b2_all.log &
    "$PY" -I $B2/filter_pruned_residues.py rama.jsonl $HOME/Desktop/Data/$dir -s $suffix --residues std20+MSE > b2_std.jsonl 2> b2_std.log &
    wait
    cmp -s b2_all.jsonl b1_new_filtered.jsonl && echo "  --residues all == b1_new_filtered.jsonl (byte-identical)" || echo "  --residues all DIFFERS from b1_new_filtered.jsonl"
    tail -n 1 b2_std.log
    "$PY" -I $B2/check_residue_set.py b2_all.jsonl b2_std.jsonl > b2_check.txt 2>&1; sed 's/^/  /' b2_check.txt
done
echo "b2 regression done $(date +%T)"
