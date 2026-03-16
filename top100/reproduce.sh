#!/bin/bash
# Reproduce the top100 pydangle dataset from scratch.
#
# Prerequisites:
#   pip install pydangle-biopython==0.5.1
#   apt install dssp
#   Reduce 4.16+ in PATH (https://github.com/rlabduke/reduce)
#
# Usage:
#   bash reproduce.sh [JOBS]    # JOBS = parallel workers (default: 8)

set -euo pipefail

JOBS=${1:-8}
DIR="top100_pdbs"
CHAIN_LIST="top100_chain_list.csv"
SUFFIX="pruned_bfilt"

echo "=== top100 reproduction pipeline ==="
echo "Workers: $JOBS"
echo ""

# 1. Download full PDB structures from RCSB
echo "Step 1: Download PDB files..."
python download_pdbs.py "$CHAIN_LIST" "$DIR" -j "$JOBS"

# 2. Run Reduce
echo ""
echo "Step 2: Run Reduce..."
python ../pipeline/run_reduce.py "$DIR" -j "$JOBS"

# 3. Build ersatz PDB files
echo ""
echo "Step 3: Build ersatz PDB files..."
python ../pipeline/build_ersatz_only.py "$DIR" -j "$JOBS"

# 4. Build B-factor masks (mainchain B <= 40, matching original Top100 methodology)
echo ""
echo "Step 4: Build B-factor masks..."
python ../pipeline/build_bfactor_masks.py "$DIR" "$CHAIN_LIST" -j "$JOBS" -b 30

# 5. Run pydangle
echo ""
echo "Step 5: Run pydangle..."
find "$DIR" -name "*_ersatz.pdb" | sort > /tmp/top100_filelist.txt
pydangle-biopython \
  -c "phi; psi; omega; tau; chi1; chi2; chi3; chi4; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality" \
  -o jsonl -j 4 \
  -f /tmp/top100_filelist.txt \
  > top100_ersatz_raw.jsonl 2> top100_ersatz.log

# 6. Post-filter
echo ""
echo "Step 6: Post-filter to quality-filtered residues..."
python ../pipeline/filter_pruned_residues.py \
  top100_ersatz_raw.jsonl "$DIR" -s "$SUFFIX" \
  > top100_measures.jsonl

# 7. Generate analysis reports
echo ""
echo "Step 7: Generate analysis reports..."
python ../pipeline/analyze_jsonl.py top100_measures.jsonl -o text -O top100_analysis.txt
python ../pipeline/analyze_jsonl.py top100_measures.jsonl -o json -O top100_analysis.json

# 8. Compress
echo ""
echo "Step 8: Compress output..."
gzip -k top100_measures.jsonl

echo ""
echo "=== Done ==="
echo "Output: top100_measures.jsonl.gz"
