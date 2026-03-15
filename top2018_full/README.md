# pydangle top2018: Per-residue geometric measurements and classifications for 12,125 high-quality protein chains

**Zenodo DOI:** [10.5281/zenodo.19032221](https://doi.org/10.5281/zenodo.19032221) (v1.0.2) | [all versions](https://doi.org/10.5281/zenodo.19029946)

## Summary

This dataset contains per-residue backbone and sidechain geometric measurements, Ramachandran classifications, and DSSP secondary structure assignments for 2,355,814 quality-filtered protein residues from 11,843 high-resolution, low-redundancy protein chains.

The measurements were computed using [pydangle-biopython](https://github.com/prisant/pydangle-biopython) v0.5.1, a Python reimagining of the Richardson Lab's Java Dangle tool, with BioPython as the structure-parsing backend.

## Source dataset

The residues come from the top2018 dataset by the Richardson Lab at Duke University (Williams, Richardson, & Richardson, 2021). This is the full-filtered, 70% homology, 60% completeness subset.

- Williams, C. J., Richardson, D. C., & Richardson, J. S. (2021). The importance of residue-level filtering, and the Top2018 best-parts dataset of high-quality protein residues. *Protein Science*. doi:10.1002/pro.4239
- Original dataset: doi:10.5281/zenodo.5115233

## Preprocessing: "Ersatz" full-structure PDB files

The top2018 provides pruned single-chain PDB files. Running DSSP on these isolated chains gives incorrect secondary structure assignments because inter-chain hydrogen bond context is lost. We created "ersatz" full-structure PDB files that combine the complete multi-chain coordinates (for correct DSSP) with NQH sidechain flip corrections (for correct Asn/Gln/His geometry).

### Pipeline

1. **Download full PDB structures** from RCSB for each chain entry using `get_top2018_full.py`.

2. **Run Reduce** on each full PDB to add hydrogens and optimize NQH sidechain orientations using three-point hinge rotation (matching the Zenodo dataset's flip method).
   - Reduce 3.3.160422 for 11,841 entries; Reduce 4.16.250520 for 2 entries (4oan, 4uzg) that segfault in 3.3.
   - Command: `reduce -db <het_dict> -build -norotmet <input>.pdb`
   - NQH flip validation: 99.77% agreement with Zenodo pruned files (463,138/464,212 decisions match).

3. **Build ersatz PDB files** by stripping hydrogens and problematic headers from the Reduce output (`build_ersatz_pdbs.py`). This preserves multi-chain context with correct NQH corrections.

4. **Run pydangle-biopython** on all ersatz files with 11 measurements: `phi; psi; omega; tau; chi1; rama_category; dssp; is_cis; is_trans; is_left; is_right`.

5. **Post-filter** to quality-filtered residues using mask files derived from the pruned PDB files (`filter_pruned_residues.py`). Fragment boundaries from `USER INC` records in the pruned files determine where neighbor-dependent measurements (phi, psi, omega, rama_category, is_cis, is_trans) are set to null, ensuring no measurement depends on coordinates from quality-filtered residues.

### NQH flip method

Reduce was run **without** `-renameflip`, producing three-point hinge-rotation flips. This was verified to produce coordinates identical to the Zenodo pruned files (confirmed by exact coordinate match for flipped HIS B 201 in 1a4i). The alternative `-renameflip` flag swaps atom names without rotating coordinates, which produces different results.

Reference: Williams, C. J., Headd, J. J., Moriarty, N. W., et al. (2018). MolProbity: More and better reference data for improved all-atom structure validation. *Protein Science*, 27(1), 293–315. doi:10.1002/pro.3330

## Output format

Each line of the JSONL file is a JSON object with these fields:

| Field | Type | Description |
|-------|------|-------------|
| `file` | string | Source filename |
| `model` | int | Model number (1 for most X-ray structures) |
| `chain` | string | Chain identifier |
| `resnum` | int | Residue sequence number |
| `ins` | string | Insertion code (space if none) |
| `resname` | string | Three-letter residue name |
| `phi` | float/null | Phi backbone dihedral (degrees) |
| `psi` | float/null | Psi backbone dihedral (degrees) |
| `omega` | float/null | Omega peptide bond dihedral (degrees) |
| `tau` | float/null | N-CA-C bond angle (degrees) |
| `chi1` | float/null | Chi1 sidechain dihedral (degrees) |
| `rama_category` | string/null | Ramachandran category (General, Gly, IleVal, TransPro, CisPro, PrePro) |
| `dssp` | string/null | DSSP secondary structure code |
| `is_cis` | string/null | "True" if omega within ±30° of 0° |
| `is_trans` | string/null | "True" if omega within ±30° of ±180° |
| `is_left` | string/null | "True" if phi > 0° (left-handed) |
| `is_right` | string/null | "True" if phi < 0° (right-handed) |

Null values indicate the measurement could not be computed (e.g., phi requires the preceding residue, which may be absent at fragment boundaries).

## Dataset statistics

- Total residues: 2,355,814
- Unique structures: 11,842
- Standard amino acids: 20 types represented
- Ramachandran categories: General 68.6%, IleVal 14.1%, Gly 8.8%, TransPro 4.7%, PrePro 3.7%, CisPro 0.3%
- Cis peptide bonds: 5,787 (0.29%), of which 5,059 are Pro (5.0% of all Pro)
- DSSP: H 31.4%, E 23.5%, C 18.1%, T 10.9%, S 8.2%, G 3.9%, P 2.0%, B 1.4%, I 0.6%

## Files

### Zenodo deposit (data products)

| File | Description |
|------|-------------|
| `top2018_pydangle_results.jsonl.gz` | Compressed JSONL measurement output (63 MB) |
| `top2018_analysis.txt` | Human-readable analysis report |
| `top2018_analysis.json` | Machine-readable analysis report |

### GitHub repository

| File | Description |
|------|-------------|
| `top2018_pdb_chain_ids.txt` | List of 11,843 PDB chain identifiers |
| `get_top2018_full.py` | Download full PDBs from RCSB and run DSSP |
| `reproduce.sh` | End-to-end reproduction script |
| `../pipeline/` | Shared pipeline scripts (build_ersatz, filter, reduce, etc.) |
| `../supplementary_pdbs/` | 3 PDB files no longer on RCSB in PDB format |

## Reproduction

### Requirements

```bash
pip install pydangle-biopython==0.5.1    # exact version used for this dataset
sudo apt install dssp                     # mkdssp 4.x
```

Reduce must be installed separately. Version 3.3+ handles most files;
version 4.16+ is needed for 4oan and 4uzg which segfault in 3.3.
See https://github.com/rlabduke/reduce.

### Quick start

The easiest way to reproduce is with the included script:

```bash
# Clone this repository
git clone https://github.com/prisant/richardson-dataset-curation.git
cd richardson-dataset-curation/top2018_full/

# Download and unpack the top2018 full-filtered dataset from Zenodo
# https://doi.org/10.5281/zenodo.5115233
# Place top2018_pdbs_full_filtered_hom70/ in this directory

# Run the full pipeline
bash reproduce.sh 8    # 8 parallel workers; takes ~2 hours total
```

### Manual step-by-step

```bash
# 1. Download full PDB structures from RCSB and run DSSP
python get_top2018_full.py

# 2. Copy supplementary PDB files (no longer available from RCSB)
for pdb in 4v4m 4yza 4ztt; do
  cp ../supplementary_pdbs/${pdb}.pdb top2018_pdbs_full_filtered_hom70/${pdb:0:2}/${pdb}/
done

# 3. Run Reduce on full PDB files (hinge flips, ~50 min at -j 8)
python ../pipeline/run_reduce.py top2018_pdbs_full_filtered_hom70/ -j 8

# 4. Build ersatz PDB files (~70 sec at -j 8)
python ../pipeline/build_ersatz_pdbs.py top2018_pdbs_full_filtered_hom70/ -j 8

# 5. Run pydangle on ersatz files
find top2018_pdbs_full_filtered_hom70/ -name "*_ersatz.pdb" | sort > filelist.txt
pydangle-biopython \
  -c "phi; psi; omega; tau; chi1; rama_category; dssp; is_cis; is_trans; is_left; is_right" \
  -o jsonl -j 4 -f filelist.txt > raw.jsonl

# 6. Post-filter to quality-filtered residues
python ../pipeline/filter_pruned_residues.py raw.jsonl \
  top2018_pdbs_full_filtered_hom70/ > top2018_pydangle_results.jsonl
```

## Software versions used

| Software | Version |
|----------|---------|
| pydangle-biopython | 0.5.1 |
| Python | 3.12.3 |
| mkdssp | 4.2.2 |
| Reduce | 3.3.160422 (11,841 entries), 4.16.250520 (2 entries) |
| BioPython | (via pip) |
| OS | Ubuntu 24.04, Linux 6.17.0-19-generic |

## Author

Michael G. Prisant, Ph.D.
Prisant Scientific
info@prisantscientific.org

## License

CC-BY 4.0

## Related

- pydangle-biopython: https://github.com/prisant/pydangle-biopython
- top2018 dataset: doi:10.5281/zenodo.5115233
- Williams et al. (2021): doi:10.1002/pro.4239
- Williams et al. (2018): doi:10.1002/pro.3330
