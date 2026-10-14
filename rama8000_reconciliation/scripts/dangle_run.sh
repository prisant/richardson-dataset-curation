#!/bin/sh
# Java Dangle 1.07 (~/bin/dangle) phi, psi, omega on every top8000 structure, Reduce output (*FH.pdb) and ersatz.
# Output: data/dangle_FH.tsv, data/dangle_ersatz.tsv (Dangle rows, colon-delimited; gitignored, regenerable) and
# data/dangle_<kind>.err (stderr). One Dangle process per file, each writing its own output (parallel writes to a
# shared stream interleave mid-line).
# Usage: sh scripts/dangle_run.sh   (from rama8000_reconciliation/)
ROOT=$HOME/Desktop/Data/top8000_pdbs_hom70
for kind in FH ersatz; do
  suf=$kind; [ $kind = ersatz ] && suf=_ersatz
  tmp=$(mktemp -d)
  ls $ROOT/*/*/*$suf.pdb | xargs -P 8 -I{} sh -c 'f={}; ~/bin/dangle "phi, psi, omega" "$f" > '$tmp'/$(basename $f).out 2> '$tmp'/$(basename $f).err'
  cat $tmp/*.out | grep -v "^#" > data/dangle_$kind.tsv
  for e in $tmp/*.err; do [ -s $e ] && { echo "== $(basename $e .err)"; cat $e; }; done > data/dangle_$kind.err
  echo "$kind: $(ls $ROOT/*/*/*$suf.pdb | wc -l) files, $(ls $tmp/*.out | wc -l) outputs, $(grep -c '^==' data/dangle_$kind.err) with stderr, $(wc -l < data/dangle_$kind.tsv) rows"
  rm -r $tmp
done
