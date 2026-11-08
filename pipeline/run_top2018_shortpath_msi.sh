#!/bin/bash
# Build + sanity-check the top2018 shortpath JSONLs (full and mc) on msi.
cd ~/Desktop/Data/top2018_shortpath
PY=~/Desktop/Projects/pydangle-biopython/.venv/bin/python
R=~/Desktop/Projects/richardson-dataset-curation
S=~/shortpath_scripts
run() {
  ds=$1; pub=$2; suffix=$3
  rm -f top2018_${ds}_shortpath.jsonl.gz top2018_${ds}_shortpath_report.txt top2018_${ds}_shortpath_provenance.json
  echo "$ds start $(date +%T)" >> progress.txt
  $PY -I $S/build_top2018_shortpath.py --dataset $ds \
      --db ../williams_top2018_sqlite/top2018_build.db --published $pub \
      --pdbredo-dir ../pdbredo_dssp/top2018_$ds --dssp461-dir ../dssp461/top2018_$ds \
      --out top2018_${ds}_shortpath > build_$ds.log 2>&1
  b=$?
  echo "$ds build exit $b $(date +%T)" >> progress.txt
  [ $b -eq 0 ] && python3 $S/sanity_check_shortpath.py top2018_${ds}_shortpath.jsonl.gz $pub \
      ../top2018_pdbs_${ds}_filtered_hom70 --suffix $suffix > sanity_$ds.txt 2>&1
  echo "$ds DONE sanity exit $? $(date +%T)" >> progress.txt
}
: > progress.txt
run full $R/top2018_full/zenodo_upload/top2018_pydangle_results.jsonl.gz pruned_all &
run mc $R/top2018_mc/zenodo_upload/top2018mc_measures.jsonl.gz pruned_mc &
wait
echo "ALL DONE $(date +%T)" >> progress.txt
