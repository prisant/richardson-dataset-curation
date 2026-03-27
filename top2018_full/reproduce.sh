#!/usr/bin/env bash
#
# reproduce.sh — End-to-end reproduction of the top2018 full-filtered
# pydangle measurement dataset.
#
# Prerequisites:
#   - Python 3.12+ with pydangle-biopython 0.5.1 installed
#   - Reduce (3.3+ or 4.16+) in PATH or ~/bin/
#   - mkdssp 4.x in PATH
#   - The top2018 full-filtered dataset from Zenodo:
#     https://doi.org/10.5281/zenodo.5115233
#     Unpack top2018_pdbs_full_filtered_hom70.tar.gz into this directory.
#
# Usage:
#   cd top2018_full/
#   bash reproduce.sh [JOBS]
#
# JOBS defaults to 8 (number of parallel workers).

set -euo pipefail

JOBS="${1:-8}"
DATASET_DIR="top2018_pdbs_full_filtered_hom70"
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PIPELINE="$REPO_ROOT/pipeline"
SUPPL="$REPO_ROOT/supplementary_pdbs"

REQUIRED_PYDANGLE="0.5.1"

echo "=== Step 0: Verify prerequisites ==="
command -v pydangle-biopython >/dev/null || { echo "ERROR: pydangle-biopython not found. Install with: pip install pydangle-biopython==$REQUIRED_PYDANGLE"; exit 1; }
command -v reduce >/dev/null || { echo "ERROR: reduce not found. See https://github.com/rlabduke/reduce"; exit 1; }
command -v mkdssp >/dev/null && true || command -v dssp >/dev/null || { echo "ERROR: mkdssp/dssp not found. Install with: apt install dssp"; exit 1; }

# Verify pydangle version
PYDANGLE_VER=$(pydangle-biopython --version 2>&1 | grep -oP '[\d.]+')
if [ "$PYDANGLE_VER" != "$REQUIRED_PYDANGLE" ]; then
    echo "WARNING: pydangle-biopython $PYDANGLE_VER found, but $REQUIRED_PYDANGLE was used to produce the published dataset."
    echo "  Results may differ. Install exact version with: pip install pydangle-biopython==$REQUIRED_PYDANGLE"
fi

echo "  pydangle-biopython $PYDANGLE_VER (published with $REQUIRED_PYDANGLE)"
echo "  $(python3 --version)"
echo "  Reduce $(reduce -version 2>&1 | grep 'reduce\.' || echo 'unknown')"
echo "  mkdssp $(mkdssp --version 2>&1 || dssp --version 2>&1 || echo 'unknown')"

if [ ! -d "$DATASET_DIR" ]; then
    echo "ERROR: $DATASET_DIR not found."
    echo "Download and unpack from https://doi.org/10.5281/zenodo.5115233"
    exit 1
fi

echo ""
echo "=== Step 1: Download full PDB structures and run DSSP ==="
python3 get_top2018_full.py 2>&1 | tail -5

echo ""
echo "=== Step 1b: Copy supplementary PDB files (no longer on RCSB) ==="
for pdb in 4v4m 4yza 4ztt; do
    prefix=${pdb:0:2}
    dest="$DATASET_DIR/$prefix/$pdb/$pdb.pdb"
    if [ ! -f "$dest" ]; then
        if [ -f "$SUPPL/$pdb.pdb" ]; then
            cp "$SUPPL/$pdb.pdb" "$dest"
            echo "  Copied $pdb.pdb from supplementary archive"
        else
            echo "  WARNING: $pdb.pdb not available (check supplementary_pdbs/)"
        fi
    fi
done

echo ""
echo "=== Step 2: Run Reduce (hinge flips, ~50 min at -j $JOBS) ==="
python3 "$PIPELINE/run_reduce.py" "$DATASET_DIR" -j "$JOBS"

echo ""
echo "=== Step 2b: Retry failures with Reduce 4.16 if available ==="
# 4oan and 4uzg segfault in Reduce 3.3; try 4.16 if available
REDUCE416="${REDUCE416:-}"
if [ -z "$REDUCE416" ]; then
    # Common locations
    for candidate in \
        ~/Desktop/Downloads/reduce-master/reduce_src/reduce \
        ~/bin/reduce4 \
        /usr/local/bin/reduce4; do
        if [ -x "$candidate" ]; then
            REDUCE416="$candidate"
            break
        fi
    done
fi
if [ -n "$REDUCE416" ]; then
    HETDB416="${HETDB416:-$(dirname "$REDUCE416")/../reduce_wwPDB_het_dict.txt}"
    for pdb in 4oan 4uzg; do
        prefix=${pdb:0:2}
        fh="$DATASET_DIR/$prefix/$pdb/${pdb}FH.pdb"
        if [ ! -s "$fh" ]; then
            src="$DATASET_DIR/$prefix/$pdb/$pdb.pdb"
            echo "  Retrying $pdb with Reduce 4.16..."
            timeout 300 "$REDUCE416" -db "$HETDB416" -build -norotmet "$src" > "$fh" 2>/dev/null || true
        fi
    done
fi

echo ""
echo "=== Step 3: Build ersatz PDB files (~70 sec at -j $JOBS) ==="
python3 "$PIPELINE/build_ersatz_pdbs.py" "$DATASET_DIR" -j "$JOBS"

echo ""
echo "=== Step 4: Run pydangle (~1-2 hours at -j 4) ==="
find "$DATASET_DIR" -name "*_ersatz.pdb" | sort > /tmp/ersatz_filelist.txt
echo "  $(wc -l < /tmp/ersatz_filelist.txt) ersatz files found"
pydangle-biopython \
    -c "phi; psi; omega; tau; chi1; chi2; chi3; chi4; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality" \
    -o jsonl -j 4 \
    -f /tmp/ersatz_filelist.txt \
    > top2018_ersatz_raw.jsonl 2> top2018_ersatz.log
echo "  Raw output: $(wc -l < top2018_ersatz_raw.jsonl) lines"

echo ""
echo "=== Step 5: Post-filter to quality-filtered residues ==="
python3 "$PIPELINE/filter_pruned_residues.py" \
    top2018_ersatz_raw.jsonl "$DATASET_DIR" \
    > top2018full_measures.jsonl

echo ""
echo "=== Step 6: Compress and analyze ==="
gzip -kf top2018full_measures.jsonl
python3 "$PIPELINE/analyze_jsonl.py" top2018full_measures.jsonl -o text > top2018_analysis.txt
python3 "$PIPELINE/analyze_jsonl.py" top2018full_measures.jsonl -o json > top2018_analysis.json

echo ""
echo "=== Done ==="
echo "Output: top2018full_measures.jsonl.gz"
echo "Analysis: top2018_analysis.txt"
ls -lh top2018full_measures.jsonl.gz top2018_analysis.txt
