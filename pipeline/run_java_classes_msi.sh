#!/bin/bash
# Classify pydangle vs Java Dangle 1.07 phi/psi/omega disagreements
# (java_diff_classes.py) before and after the altloc change, on msi.
# Before: rama.jsonl (pydangle main 664e87a; identical to --altloc
# occupancy).  After: alt_first_v2.jsonl (--altloc first, final rule;
# run_altloc_regress_msi.sh with TAG=_v2).  java.txt: Java Dangle 1.07 on
# the same file list.  top100 has no Java Dangle run.
# Output: ~/scratch/regress/<ds>/java_classes_{before,after}.{txt,tsv}.
set -u
PY=$HOME/Desktop/Projects/pydangle-biopython/.venv/bin/python
R=$HOME/scratch/regress
for ds in top500 top8000 top2018_full top2018_mc; do
    for t in before:rama.jsonl after:alt_first_v2.jsonl; do
        (cd $R/$ds && "$PY" -I $HOME/scratch/b6/java_diff_classes.py filelist.txt ${t#*:} java.txt \
            --list 25 --tsv java_classes_${t%%:*}.tsv > java_classes_${t%%:*}.txt 2>&1) &
    done
done
wait
echo "java classes done $(date +%T)"
