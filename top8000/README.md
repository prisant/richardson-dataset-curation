# pydangle top8000: Per-residue backbone geometry for 7,814 quality-filtered protein chains

**Zenodo DOI:** [10.5281/zenodo.19041232](https://doi.org/10.5281/zenodo.19041232)

## Summary

This dataset contains per-residue backbone geometric measurements, Ramachandran classifications at multiple granularities, DSSP secondary structure assignments, peptide bond classification, and chirality for 1,568,561 quality-filtered protein residues from 7,814 chains in 7,607 high-resolution, low-redundancy protein structures.

The measurements were computed using [pydangle-biopython](https://github.com/prisant/pydangle-biopython) v0.5.1, a Python reimagining of the Richardson Lab's Java Dangle tool, with BioPython as the structure-parsing backend.

## Source dataset

The residues come from the Top8000 dataset by the Richardson Lab at Duke University. This is the best-quality, 70% homology chain list (Top8000-best_hom70, 7,957 chains from 7,737 structures).

- Williams, C. J., Headd, J. J., Moriarty, N. W., Prisant, M. G., et al. (2018). MolProbity: More and better reference data for improved all-atom structure validation. *Protein Science*, 27(1), 293–315. doi:10.1002/pro.3330
- Hintze, B. J., Lewis, S. M., Richardson, J. S., & Richardson, D. C. (2016). MolProbity's Ultimate Rotamer-Library Distributions for Model Validation. *Proteins*, 84(9), 1177–1189. doi:10.1002/prot.25039
- Chain lists and reference data: https://github.com/rlabduke/reference_data

## Residue-level quality filter

Residues are filtered using the mainchain B-factor criterion described in Williams et al. (2018) for Ramachandran and CaBLAM contour generation: any residue where a mainchain atom (N, CA, C, or O) has B-factor > 30 is excluded. This simple filter is applied to all residues in each quality-filtered chain.

Of the 7,957 chains in the chain list, 143 (1.8%) have zero residues passing this filter. These chains passed the chain-level quality criteria (MolProbity score < 2.0, resolution < 2.0 Å) but have elevated B-factors throughout, likely due to molecular flexibility or crystal packing effects. See `top8000_issues.md` for detailed analysis.

## Preprocessing: "Ersatz" full-structure PDB files

The Top8000 provides single-chain PDB files. Running DSSP on these isolated chains gives incorrect secondary structure assignments because inter-chain hydrogen bond context is lost. We created "ersatz" full-structure PDB files that combine the complete multi-chain coordinates (for correct DSSP) with NQH sidechain flip corrections from Reduce.

### Pipeline

1. **Download full PDB structures** from RCSB for each structure in the chain list (4,613 shared with top2018 working directories, 3,124 downloaded fresh).

2. **Run Reduce 4.16.250520** on each full PDB to add hydrogens and optimize NQH sidechain orientations.

3. **Build ersatz PDB files** by stripping hydrogens and problematic headers from the Reduce output (`build_ersatz_only.py`).

4. **Build B-factor masks** from the ersatz files (`build_bfactor_masks.py -b 30`). For each chain in the chain list, identifies residues where all mainchain atoms have B ≤ 30 and generates mask files and synthetic fragment records for contiguous passing regions.

5. **Run pydangle-biopython** on all 7,737 ersatz files with 11 measurements.

6. **Post-filter** to quality-filtered residues using mask files and synthetic fragment records (`filter_pruned_residues.py -s pruned_bfilt`).

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
| `rama_category` | string/null | 6-class Ramachandran |
| `rama5` | string/null | 5-class (merges CisPro into TransPro) |
| `rama4` | string/null | 4-class (also merges PrePro into General) |
| `rama3` | string/null | 3-class (General, Gly, Pro) |
| `dssp` | string/null | DSSP secondary structure code |
| `peptide_bond` | string/null | "cis", "trans", or "twisted" |
| `chirality` | string/null | "L" or "D" |

## Dataset statistics

- Total residues: 1,568,561
- Unique structures: 7,607
- Output chains: 7,814
- Chain list entries excluded (all fail B≤30): 143 (1.8%)
- Standard amino acids: 20 types + 1 CSD (non-standard)

See `top8000_analysis.txt` for the full statistical report and `top8000_issues.md` for known issues.

## Software versions used

| Software | Version |
|----------|---------|
| pydangle-biopython | 0.5.1 |
| Python | 3.12.3 |
| mkdssp | 4.2.2 |
| Reduce | 4.16.250520 |
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
- top2018 MC-filtered sibling: doi:10.5281/zenodo.19040026
- Richardson Lab reference_data: https://github.com/rlabduke/reference_data
- Williams et al. (2018): doi:10.1002/pro.3330
- Hintze et al. (2016): doi:10.1002/prot.25039
