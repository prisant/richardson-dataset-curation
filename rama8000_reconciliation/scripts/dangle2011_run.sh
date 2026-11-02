#!/bin/sh
# Java Dangle as of mid-2011 (rlabduke/javadev commit before 2011-07-01, Dangle 1.07) on the lab's 2011 chain files.
# Build: git clone https://github.com/rlabduke/javadev; git archive $(git rev-list -1 --before=2011-07-01 HEAD)
#   javadev/driftwood javadev/chiropraxis javadev/dangle | tar x -C jd2011; in jd2011/javadev:
#   javac -nowarn -encoding ISO-8859-1 -d classes -sourcepath driftwood/src:chiropraxis/src:dangle/src dangle/src/dangle/Dangle.java
#   and copy driftwood/, chiropraxis/, dangle/ resource/* into classes/.
# usage: sh scripts/dangle2011_run.sh JD2011_CLASSES CHAIN_DIR OUT.tsv     (from rama8000_reconciliation/)
CP=$1; SRC=$2; OUT=$3
tmp=$(mktemp -d)
ls $SRC/*.pdb | split -l 100 - $tmp/batch.
ls $tmp/batch.* | xargs -P 8 -I{} sh -c 'java -cp '$CP' dangle.Dangle "phi, psi, omega" $(cat {}) > {}.out 2> {}.err'
cat $tmp/batch.*.out | grep -v "^#" > $OUT
cat $tmp/batch.*.err > ${OUT%.tsv}.err
echo "$(ls $SRC/*.pdb | wc -l) files, $(wc -l < $OUT) rows, stderr lines $(wc -l < ${OUT%.tsv}.err)"
rm -r $tmp
