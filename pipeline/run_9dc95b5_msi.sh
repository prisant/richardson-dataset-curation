#!/bin/bash
# Rerun all five datasets with pydangle 9dc95b5 (branch altloc-first) on msi.
#
#   pydangle  ~/scratch/pydangle_9dc95b5/src: git archive of commit 9dc95b5,
#             checked file by file (md5) against the commit before the run;
#             --altloc first; via PYTHONPATH
#   datasets  our current-file datasets (~/Desktop/Data/<dir>, file lists
#             ~/scratch/regress/<ds>/filelist.txt), all five
#   filter    ~/scratch/p9dc95b5/scripts/filter_pruned_residues.py (copy of
#             the b2-residue-set worktree version, recorded by md5 in the
#             log), masks of the dataset (suffix below), --residues std20+MSE,
#             in both neighbour-angle modes:
#               null  neighbour-dependent fields nulled at fragment ends
#               keep  values kept as measured (lab Top8000 rule); Ramachandran
#                     fields nulled where phi or psi is null
# Outputs per dataset in ~/scratch/regress/<ds>/p9dc95b5/ (nothing older is
# overwritten):
#   all.jsonl            pydangle, every residue
#   std_null.jsonl       filtered, neighbour angles null
#   std_keep.jsonl       filtered, neighbour angles keep
#   twisted.txt          twisted_peptides.py on both filtered files
# Waits until run_lab2011_msi.sh has finished (both use 10 workers).
set -euo pipefail
PY=$HOME/Desktop/Projects/pydangle-biopython/.venv/bin/python
SRC=$HOME/scratch/pydangle_9dc95b5/src
S=$HOME/scratch/p9dc95b5/scripts
R=$HOME/scratch/regress
RUN='from pydangle_biopython.main import main; import sys; sys.exit(main())'
FULL="phi; psi; omega; tau; chi1; chi2; chi3; chi4; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality"
MC="phi; psi; omega; tau; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality"

while pgrep -f run_lab2011_msi.sh > /dev/null; do sleep 60; done
echo "pydangle: $(cat $SRC/../COMMIT); filter md5 $(md5sum < $S/filter_pruned_residues.py | cut -c1-32)"
for spec in top100:top100_pdbs:pruned_bfilt top500:top500_pdbs:pruned_bfilt \
            top8000:top8000_pdbs_hom70:pruned_bfilt \
            top2018_full:top2018_pdbs_full_filtered_hom70:pruned_all \
            top2018_mc:top2018_pdbs_mc_filtered_hom70:pruned_mc; do
    IFS=: read ds dir suffix <<< "$spec"
    cmds=$FULL; [ $ds = top2018_mc ] && cmds=$MC
    O=$R/$ds/p9dc95b5
    mkdir -p $O
    echo "== $ds $(date +%T) files $(wc -l < $R/$ds/filelist.txt)"
    (cd $R/$ds && PYTHONPATH=$SRC "$PY" -W ignore -c "$RUN" --altloc first -c "$cmds" -o jsonl -j 10 \
        -f filelist.txt > $O/all.jsonl 2> $O/all.log) || echo "  pydangle exit $?"
    for nb in null keep; do
        "$PY" -I $S/filter_pruned_residues.py $O/all.jsonl $HOME/Desktop/Data/$dir -s $suffix \
            --residues std20+MSE --neighbour-angles $nb > $O/std_$nb.jsonl 2> $O/std_$nb.log
        sed 's/^/  '$nb': /' $O/std_$nb.log
    done
    "$PY" -I $S/twisted_peptides.py "null=$O/std_null.jsonl,$O/all.jsonl" \
        "keep=$O/std_keep.jsonl,$O/all.jsonl" > $O/twisted.txt
    sed 's/^/  /' $O/twisted.txt
done
echo "done $(date +%T)"
