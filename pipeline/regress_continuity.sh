#!/bin/bash
# Continuity-fix regression driver: old pydangle (main) vs fixed (continuity-fix)
# vs Java Dangle 1.07, on one dataset's ersatz files.
#
# Usage: regress_continuity.sh DATASET_DIR OUT_DIR [JOBS] [COMMANDS]
#   DATASET_DIR  e.g. ~/Desktop/Data/top100_pdbs
#   OUT_DIR      receives filelist.txt, old.jsonl, new.jsonl, java.txt, logs
#   COMMANDS     pydangle command string; default = production string from
#                richardson-dataset-curation/top100/reproduce.sh
#
# Environment:
#   PYDANGLE_REPO  pydangle-biopython checkout on continuity-fix (venv in .venv)
#   OLD_SRC        src/ of a pydangle checkout at main (d979097)
#   DANGLE_JAR     Java Dangle 1.07 jar
set -euo pipefail

DATASET_DIR=$1
OUT_DIR=$2
JOBS=${3:-8}
COMMANDS=${4:-"phi; psi; omega; tau; chi1; chi2; chi3; chi4; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality"}
PYDANGLE_REPO=${PYDANGLE_REPO:-$HOME/Desktop/Projects/pydangle-biopython}
OLD_SRC=${OLD_SRC:?set OLD_SRC to the src/ of a pydangle checkout at main}
DANGLE_JAR=${DANGLE_JAR:-$HOME/bin/javadev-master/dangle/dangle.jar}
PY=$PYDANGLE_REPO/.venv/bin/python
RUN='from pydangle_biopython.main import main; import sys; sys.exit(main())'

mkdir -p "$OUT_DIR"
find "$DATASET_DIR" -name "*_ersatz.pdb" | sort > "$OUT_DIR/filelist.txt"
echo "files: $(wc -l < "$OUT_DIR/filelist.txt")"
echo "old pydangle: $(git -C "$OLD_SRC/.." rev-parse --short HEAD)"
echo "new pydangle: $(git -C "$PYDANGLE_REPO" rev-parse --abbrev-ref HEAD) $(git -C "$PYDANGLE_REPO" rev-parse --short HEAD) + working tree"

# cd to OUT_DIR so neither run can pick up a stray module from the cwd
cd "$OUT_DIR"
echo "== old $(date +%T)"
PYTHONPATH="$OLD_SRC" "$PY" -c "$RUN" -c "$COMMANDS" -o jsonl -j "$JOBS" \
    -f filelist.txt > old.jsonl 2> old.log || echo "old exit $?"
echo "== new $(date +%T)"
"$PY" -c "$RUN" -c "$COMMANDS" -o jsonl -j "$JOBS" \
    -f filelist.txt > new.jsonl 2> new.log || echo "new exit $?"
echo "== java $(date +%T)"
: > java.txt
while read -r f; do
    java -jar "$DANGLE_JAR" phi psi omega tau "$f" 2>>java.log | grep -v '^#' >> java.txt || true
done < filelist.txt
echo "== done $(date +%T)"
wc -l old.jsonl new.jsonl java.txt
