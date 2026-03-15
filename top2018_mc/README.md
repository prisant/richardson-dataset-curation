# pydangle top2018 MC: Per-residue backbone geometry for 13,677 mainchain-filtered protein chains

**Zenodo DOI:** [10.5281/zenodo.19040026](https://doi.org/10.5281/zenodo.19040026)

## Summary

This dataset contains per-residue backbone geometric measurements, Ramachandran classifications at multiple granularities, DSSP secondary structure assignments, peptide bond classification, and chirality for 2,963,303 quality-filtered protein residues from 13,677 mainchain-filtered chains in 13,308 high-resolution, low-redundancy protein structures.

The measurements were computed using [pydangle-biopython](https://github.com/prisant/pydangle-biopython) v0.5.1, a Python reimagining of the Richardson Lab's Java Dangle tool, with BioPython as the structure-parsing backend.

## Source dataset

The residues come from the top2018 dataset by the Richardson Lab at Duke University (Williams, Richardson, & Richardson, 2021). This is the mainchain-filtered (MC), 70% homology, 60% completeness subset. Mainchain-only quality filtering passes more residues than the full (MC+SC) filtering because sidechain quality is not required.

- Williams, C. J., Richardson, D. C., & Richardson, J. S. (2021). The importance of residue-level filtering, and the Top2018 best-parts dataset of high-quality protein residues. *Protein Science*. doi:10.1002/pro.4239
- Original dataset: doi:10.5281/zenodo.5115233

## Comparison with full-filtered sibling

| Property | Full-filtered (MC+SC) | MC-filtered (this dataset) |
|----------|----------------------|---------------------------|
| Quality criteria | Mainchain + sidechain | Mainchain only |
| Chains | 12,125 | 13,677 |
| Structures | 11,843 | 13,308 |
| Residues | 2,405,322 | 2,963,303 |
| Includes chi1 | Yes | No |
| Rama granularity | rama_category only | rama_category + rama5/4/3 |
| Zenodo DOI | [10.5281/zenodo.19029946](https://doi.org/10.5281/zenodo.19029946) | This deposit |

## Preprocessing: "Ersatz" full-structure PDB files

The top2018 provides pruned single-chain PDB files. Running DSSP on these isolated chains gives incorrect secondary structure assignments because inter-chain hydrogen bond context is lost. We created "ersatz" full-structure PDB files that combine the complete multi-chain coordinates (for correct DSSP) with NQH sidechain flip corrections (for correct Asn/Gln/His geometry).

### Pipeline

1. **Download full PDB structures** from RCSB for each chain entry.

2. **Run Reduce 4.16.250520** on each full PDB to add hydrogens and optimize NQH sidechain orientations using three-point hinge rotation.
   - Command: `reduce -db <het_dict> -allalt -build -norotmet <input>.pdb`
   - All 13,308 structures processed with Reduce 4.16 and the current het dictionary.

3. **Build ersatz PDB files** by stripping hydrogens and problematic headers from the Reduce output (`build_ersatz_pdbs.py -s pruned_mc`). This preserves multi-chain context with correct NQH corrections. Per-chain mask files generated from the pruned MC files.

4. **Run pydangle-biopython** on all 13,308 ersatz files with 11 measurements: `phi; psi; omega; tau; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality`.

5. **Post-filter** to quality-filtered residues using mask files and `USER INC` fragment records (`filter_pruned_residues.py -s pruned_mc`). Fragment boundaries determine where neighbor-dependent measurements are set to null.

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
| `rama_category` | string/null | 6-class Ramachandran (General, Gly, IleVal, TransPro, CisPro, PrePro) |
| `rama5` | string/null | 5-class (merges CisPro into TransPro) |
| `rama4` | string/null | 4-class (also merges PrePro into General) |
| `rama3` | string/null | 3-class (General, Gly, Pro) |
| `dssp` | string/null | DSSP secondary structure code |
| `peptide_bond` | string/null | "cis", "trans", or "twisted" |
| `chirality` | string/null | "L" or "D" |

Null values indicate the measurement could not be computed (terminal residues, fragment boundaries, missing atoms).

## Dataset statistics

- Total residues: 2,963,303
- Unique structures: 13,308
- Quality-filtered chains: 13,677 (329 structures contribute multiple chains)
- Standard amino acids: 20 types represented (2 MSE mapped to MET)
- Ramachandran categories: General 70.8%, IleVal 13.5%, Gly 7.8%, TransPro 4.3%, PrePro 3.4%, CisPro 0.2%
- DSSP: H 32.4%, E 23.0%, C 17.7%, T 11.0%, S 8.1%, G 4.0%, P 2.0%, B 1.3%, I 0.6%
- Data completeness: phi 93.5%, psi 93.5%, omega 93.5%, tau 100.0%, dssp 99.3%, rama_category 87.8%

See `top2018mc_analysis.txt` for the full statistical report and `top2018mc_ISSUES.md` for known issues.

## Files

### Zenodo deposit (data products)

| File | Description |
|------|-------------|
| `top2018mc_measures.jsonl.gz` | Compressed JSONL measurement output (71 MB) |
| `top2018mc_analysis.txt` | Human-readable analysis report |
| `top2018mc_analysis.json` | Machine-readable analysis report |
| `README.md` | This file |

### GitHub repository

| File | Description |
|------|-------------|
| `../pipeline/` | Shared pipeline scripts (build_ersatz, filter, reduce, etc.) |
| `../supplementary_pdbs/` | PDB files no longer on RCSB in PDB format |

## Reproduction

### Requirements

```bash
pip install pydangle-biopython==0.5.1    # exact version used for this dataset
sudo apt install dssp                     # mkdssp 4.x
```

Reduce 4.16+ must be installed separately.
See https://github.com/rlabduke/reduce.

### Manual step-by-step

```bash
# 1. Place top2018 MC-filtered dataset in working directory
#    Source: top2018_pdbs_mc_filtered_hom70/ with pruned_mc files

# 2. Run Reduce on full PDB files (~55 min at -j 8)
python ../pipeline/run_reduce.py top2018_pdbs_mc_filtered_hom70/ -j 8

# 3. Build ersatz PDB files and masks (~2 min at -j 8)
python ../pipeline/build_ersatz_pdbs.py top2018_pdbs_mc_filtered_hom70/ -j 8 -s pruned_mc

# 4. Run pydangle on ersatz files (~35 min at -j 4)
find top2018_pdbs_mc_filtered_hom70/ -name "*_ersatz.pdb" | sort > filelist.txt
pydangle-biopython \
  -c "phi; psi; omega; tau; rama_category; rama5; rama4; rama3; dssp; peptide_bond; chirality" \
  -o jsonl -j 4 -f filelist.txt > raw.jsonl

# 5. Post-filter to quality-filtered residues
python ../pipeline/filter_pruned_residues.py raw.jsonl \
  top2018_pdbs_mc_filtered_hom70/ -s pruned_mc > top2018mc_measures.jsonl
```

## Software versions used

| Software | Version |
|----------|---------|
| pydangle-biopython | 0.5.1 |
| Python | 3.12.3 |
| mkdssp | 4.2.2 |
| Reduce | 4.16.250520 (all entries) |
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
- top2018 full-filtered sibling: doi:10.5281/zenodo.19029946
- top2018 dataset: doi:10.5281/zenodo.5115233
- Williams et al. (2021): doi:10.1002/pro.4239
- Williams et al. (2018): doi:10.1002/pro.3330
