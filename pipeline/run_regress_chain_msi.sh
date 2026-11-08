#!/bin/bash
# Continuity-fix regression chain on msi: top500 Java-step check, then
# top8000, top2018_full, top2018_mc.  Run from ~/scratch.
set -u
export OLD_SRC=$HOME/scratch/pydangle_main/src
D=$HOME/Desktop/Data
R=$HOME/scratch/regress
MC="phi; psi; omega; tau; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality"
# 1. parallel Java step must reproduce the serial top500 java.txt
mkdir -p $R/java_check && cd $R/java_check
sort $R/top500/java.txt > serial.sorted
cp $R/top500/filelist.txt filelist.txt
rm -rf java_parts && mkdir java_parts
export DANGLE_JAR=$HOME/bin/javadev-master/dangle/dangle.jar
xargs -P 10 -I{} sh -c 'java -jar "$DANGLE_JAR" phi psi omega tau "$1" 2>/dev/null | grep -v "^#" > "java_parts/$(basename "$1").txt" || true' _ {} < filelist.txt
cat java_parts/*.txt | sort > parallel.sorted
if cmp -s serial.sorted parallel.sorted; then echo "java check: identical $(wc -l < parallel.sorted)"; else echo "java check: DIFFER"; exit 1; fi
# 2. datasets
$HOME/scratch/regress_continuity.sh $D/top8000_pdbs_hom70 $R/top8000 10
$HOME/scratch/regress_continuity.sh $D/top2018_pdbs_full_filtered_hom70 $R/top2018_full 10
$HOME/scratch/regress_continuity.sh $D/top2018_pdbs_mc_filtered_hom70 $R/top2018_mc 10 "$MC"
echo "chain done $(date +%T)"
