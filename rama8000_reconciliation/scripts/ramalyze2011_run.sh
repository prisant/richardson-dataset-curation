#!/bin/sh
# Java Ramalyze as of mid-2011 (rlabduke/javadev commit before 2011-07-01; Top500-era distributions) on the lab's 2011
# chain files, -raw output. Build as for dangle2011_run.sh, plus javadev/lib/iText and javadev/lib/Jama101 from the same
# commit, then: javac -nowarn -encoding ISO-8859-1 -cp lib/iText/iText-2.1.3.jar:lib/Jama101/jama.jar -d classes
#   -sourcepath driftwood/src:chiropraxis/src chiropraxis/src/chiropraxis/rotarama/Ramalyze.java
# Output rows: <file>:<Ramalyze raw row>.
# usage: sh scripts/ramalyze2011_run.sh JD2011_JAVADEV_DIR CHAIN_DIR OUT.tsv     (from rama8000_reconciliation/)
JD=$1; SRC=$2; OUT=$3
CP=$JD/classes:$JD/lib/iText/iText-2.1.3.jar:$JD/lib/Jama101/jama.jar
tmp=$(mktemp -d)
ls $SRC/*.pdb | xargs -P 8 -I{} sh -c 'f={}; b=$(basename $f .pdb); java -cp '$CP' chiropraxis.rotarama.Ramalyze -raw $f 2> '$tmp'/$b.err | sed "s/^/$b:/" > '$tmp'/$b.out'
cat $tmp/*.out > $OUT
for e in $tmp/*.err; do [ -s $e ] && { echo "== $(basename $e .err)"; cat $e; }; done > ${OUT%.tsv}.err
echo "$(ls $SRC/*.pdb | wc -l) files, $(wc -l < $OUT) rows, files with stderr $(grep -c '^==' ${OUT%.tsv}.err)"
rm -r $tmp
