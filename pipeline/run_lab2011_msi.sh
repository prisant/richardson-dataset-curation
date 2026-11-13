#!/bin/bash
# Run our v1.x top8000 pipeline on the lab's 2011 chain files and compare
# with the lab's lists (Top8000_Rama_2013), residue by residue.
#
# Same input files as the lab's lists, so every difference left is a
# difference of measurement (pydangle) or of filter rules, not of input.
#
# Layout under $W (default ~/scratch/lab2011):
#   top8000_chains_70.tgz   lab's 2011 chain files (sha256 checked)
#   chains/                 extracted
#   tree/                   FH files (lab2011_fh_tree.py), then ersatz files,
#                           masks and pruned_bfilt files (reproduce.sh steps 3-4)
#   pydangle_src/           pydangle branch altloc-first, commit 9dc95b5
#                           (git archive), used via PYTHONPATH
#   scripts/                pipeline scripts of richardson-dataset-curation
#   lab_residues/           Top8000_Rama_2013 residues/rama8000-*.csv
# Outputs in $W/out:
#   all.jsonl               pydangle, every residue (--altloc first)
#   std_keep.jsonl          filtered: B <= 30 masks, std20+MSE,
#                           --neighbour-angles keep
#   compare.txt/.tsv        rama_residue_compare.py
#   breakdown.txt/.tsv      lab_rule_breakdown.py (group "lab: other")
set -euo pipefail
W=${W:-$HOME/scratch/lab2011}
JOBS=${JOBS:-10}
PY=$HOME/Desktop/Projects/pydangle-biopython/.venv/bin/python
S=$W/scripts
CMDS="phi; psi; omega; tau; chi1; chi2; chi3; chi4; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality"
cd "$W"
mkdir -p out
log() { echo "== $(date +%H:%M:%S) $*"; }

log "check chain-file archive"
echo "4a6bb1b9bb340cd42c14cd9054688732ac09317dfe0af49f3fd453f77ba486e6  top8000_chains_70.tgz" | sha256sum -c
if [ ! -d chains ]; then
    mkdir chains && tar -xzf top8000_chains_70.tgz -C chains
fi
CH=$(find chains -name '*FH_*.pdb' -printf '%h\n' | sort -u)
echo "chain files: $(find chains -name '*FH_*.pdb' | wc -l) in $CH"

log "step 0: FH tree"
"$PY" -I $S/lab2011_fh_tree.py "$CH" tree
log "step 3: ersatz files"
"$PY" -I $S/build_ersatz_only.py tree -j $JOBS | tail -3
log "step 4: B-factor masks (mainchain B <= 30)"
"$PY" -I $S/build_bfactor_masks.py tree tree/Top8000-best_hom70_pdb_chain_list.csv -j $JOBS -b 30 | tail -3

log "step 5: pydangle $(cat pydangle_src/COMMIT) --altloc first"
find tree -name '*_ersatz.pdb' | sort > out/filelist.txt
echo "files: $(wc -l < out/filelist.txt)"
PYTHONPATH=$W/pydangle_src/src "$PY" -W ignore \
    -c 'from pydangle_biopython.main import main; import sys; sys.exit(main())' \
    --altloc first -c "$CMDS" -o jsonl -j $JOBS -f out/filelist.txt \
    > out/all.jsonl 2> out/all.log
head -c 400 out/all.jsonl; echo

log "step 6: filter (std20+MSE, neighbour angles keep)"
"$PY" -I $S/filter_pruned_residues.py out/all.jsonl tree -s pruned_bfilt \
    --residues std20+MSE --neighbour-angles keep > out/std_keep.jsonl 2> out/std_keep.log
cat out/std_keep.log

log "compare with the lab's lists"
"$PY" -I $S/rama_residue_compare.py lab_residues out/std_keep.jsonl out/all.jsonl \
    --tsv out/compare.tsv > out/compare.txt
cat out/compare.txt
"$PY" -I $S/lab_rule_breakdown.py out/compare.tsv tree --tsv out/breakdown.tsv > out/breakdown.txt
cat out/breakdown.txt
log done
