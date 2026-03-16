#!/bin/bash
# Reproduce the top2018 MC-filtered pydangle dataset from scratch.
#
# Prerequisites:
#   pip install pydangle-biopython==0.5.1
#   apt install dssp
#   Reduce 4.16+ in PATH (https://github.com/rlabduke/reduce)
#
# Source data:
#   The MC-filtered pruned files must be obtained from the Richardson Lab
#   top2018 dataset (doi:10.5281/zenodo.5115233). Place the two-level
#   directory structure with {pdbid}.pdb and {pdbid}_{chain}_pruned_mc.pdb
#   files in the working directory as top2018_pdbs_mc_filtered_hom70/.
#
# Usage:
#   bash reproduce.sh [JOBS]    # JOBS = parallel workers (default: 8)

set -euo pipefail

JOBS=${1:-8}
DIR="top2018_pdbs_mc_filtered_hom70"
SUFFIX="pruned_mc"

echo "=== top2018 MC-filtered reproduction pipeline ==="
echo "Workers: $JOBS"
echo ""

# 0. Download full PDB structures from RCSB
echo "Step 0: Download PDB files from RCSB..."
echo "(Pruned MC files must be placed manually from Zenodo doi:10.5281/zenodo.5115233)"
python get_top2018_mc.py -j "$JOBS"

if [ ! -d "$DIR" ]; then
    echo "Error: $DIR not found after download."
    echo "See README.md for source data instructions."
    exit 1
fi

# 1. Run Reduce
echo ""
echo "Step 1: Run Reduce..."
python ../pipeline/run_reduce.py "$DIR" -j "$JOBS"

# 2. Build ersatz PDB files and masks from pruned files
echo ""
echo "Step 2: Build ersatz PDB files and masks..."
python ../pipeline/build_ersatz_pdbs.py "$DIR" -j "$JOBS" -s "$SUFFIX"

# 3. Run pydangle
echo ""
echo "Step 3: Run pydangle..."
find "$DIR" -name "*_ersatz.pdb" | sort > /tmp/mc_filelist.txt
pydangle-biopython \
  -c "phi; psi; omega; tau; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality" \
  -o jsonl -j 4 \
  -f /tmp/mc_filelist.txt \
  > top2018mc_ersatz_raw.jsonl 2> top2018mc_ersatz.log

# 4. Post-filter
echo ""
echo "Step 4: Post-filter to quality-filtered residues..."
python ../pipeline/filter_pruned_residues.py \
  top2018mc_ersatz_raw.jsonl "$DIR" -s "$SUFFIX" \
  > top2018mc_measures.jsonl

# 5. Generate analysis reports
echo ""
echo "Step 5: Generate analysis reports..."
python ../pipeline/analyze_jsonl.py top2018mc_measures.jsonl -o text -O top2018mc_analysis.txt
python ../pipeline/analyze_jsonl.py top2018mc_measures.jsonl -o json -O top2018mc_analysis.json

# 6. Compress
echo ""
echo "Step 6: Compress output..."
gzip -k top2018mc_measures.jsonl

echo ""
echo "=== Done ==="
echo "Output: top2018mc_measures.jsonl.gz"
