# pydangle top500: Per-residue backbone geometry for 514 quality-filtered protein chains

**Zenodo DOI:** [10.5281/zenodo.19049923](https://doi.org/10.5281/zenodo.19049923) (v1.1.0) | [all versions](https://doi.org/10.5281/zenodo.19042131)
**Source:** Richardson Lab Top500 dataset

## Summary

This dataset contains per-residue backbone and sidechain geometric measurements (including chi1–chi4 torsion angles), Ramachandran classifications at multiple granularities, DSSP secondary structure assignments, peptide bond classification, and chirality for 97,519 quality-filtered protein residues from 514 chains in 492 protein structures.

The measurements were computed using [pydangle-biopython](https://github.com/prisant/pydangle-biopython) v0.5.1, a Python reimagining of the Richardson Lab's Java Dangle tool, with BioPython as the structure-parsing backend.

## Source dataset

The Top500 is the second generation of Richardson Lab quality-filtered protein structure reference datasets, curated at Duke University. It built on the original Top100 dataset of 100 high-resolution (0.9–1.7 Å) protein structures that was used for the foundational work on Asn/Gln sidechain amide flip correction using all-atom contact analysis with Probe and Reduce (Word et al., 1999a, 1999b). The Top500 expanded the dataset to 500 structures and was used to develop the "penultimate" rotamer library (Lovell et al., 2000) and to establish Ramachandran validation contours with density-dependent smoothing and the Cbeta deviation validation measure (Lovell et al., 2003). The original analysis used 81,234 non-Gly, non-Pro, non-prePro residues with B < 30 from 500 high-resolution proteins.

The validated coordinates — with computationally added and optimized hydrogen atoms and corrected Asn/Gln/His sidechain flips — became the reference standard for MolProbity's structure validation services (Davis et al., 2004).

The Top500 was subsequently superseded by the Top8000 (Hintze et al., 2016; Williams et al., 2018), which provided an order-of-magnitude increase in data and stricter residue-level filtering including real-space correlation and local electron density criteria.

### Key references

- Word, J. M., Lovell, S. C., Richardson, J. S., & Richardson, D. C. (1999a). Asparagine and glutamine: using hydrogen atom contacts in the choice of side-chain amide orientation. *J. Mol. Biol.*, 285(4), 1735–1747. doi:10.1006/jmbi.1998.2401
- Word, J. M., Lovell, S. C., Richardson, J. S., & Richardson, D. C. (1999b). Visualizing and quantifying molecular complementarity. *J. Mol. Biol.*, 285(4), 1711–1733. doi:10.1006/jmbi.1998.2400
- Lovell, S. C., Word, J. M., Richardson, J. S., & Richardson, D. C. (2000). The penultimate rotamer library. *Proteins*, 40(3), 389–408. PMID:10861930
- Lovell, S. C., Davis, I. W., Arendall, W. B. III, de Bakker, P. I. W., Word, J. M., Prisant, M. G., Richardson, J. S., & Richardson, D. C. (2003). Structure validation by Calpha geometry: phi,psi and Cbeta deviation. *Proteins*, 50(3), 437–450. doi:10.1002/prot.10286
- Davis, I. W., Murray, L. W., Richardson, J. S., & Richardson, D. C. (2004). MOLPROBITY: structure validation and all-atom contact analysis for nucleic acids and their complexes. *Nucleic Acids Research*, 32(suppl_2), W615–W619. doi:10.1093/nar/gkh398
- Hintze, B. J., Lewis, S. M., Richardson, J. S., & Richardson, D. C. (2016). MolProbity's Ultimate Rotamer-Library Distributions for Model Validation. *Proteins*, 84(9), 1177–1189. doi:10.1002/prot.25039
- Williams, C. J., Headd, J. J., Moriarty, N. W., Prisant, M. G., et al. (2018). MolProbity: More and better reference data for improved all-atom structure validation. *Protein Science*, 27(1), 293–315. doi:10.1002/pro.3330

### Chain list provenance

The Top500 chain list was cached from the Richardson Lab kinemage webserver (kinemage.biochem.duke.edu) prior to its retirement. The original FH (Flipped + Hydrogens) single-chain files processed with Reduce v2.13.2 are preserved in the kinemage archive. For this dataset, full multi-chain PDB files were downloaded from RCSB and processed with Reduce 4.16 to provide correct multi-chain context for DSSP.

## Residue-level quality filter

Residues are filtered using the mainchain B-factor ≤ 30 criterion, matching the original Top500 analysis (Lovell et al., 2003: "residues with B < 30") and the subsequent Top8000 Ramachandran methodology (Williams et al., 2018). Any residue where a mainchain atom (N, CA, C, or O) has B-factor > 30 is excluded.

Of 523 chain entries, 9 (1.7%) have zero residues passing this filter. These chains have elevated B-factors throughout, consistent with the similar 1.8% exclusion rate observed in the top8000 dataset.

## Preprocessing: "Ersatz" full-structure PDB files

### Pipeline

1. **Full PDB structures** downloaded from RCSB for each entry (494 structures, current coordinates).
2. **Run Reduce 4.16.250520** on each full PDB for NQH sidechain flip corrections.
3. **Build ersatz PDB files** by stripping hydrogens and problematic headers.
4. **Build B-factor masks** identifying residues where all mainchain atoms have B ≤ 30, with synthetic fragment records for contiguous passing regions.
5. **Run pydangle-biopython** on all ersatz files with 15 measurements (backbone geometry, chi1–chi4, Rama classifications, DSSP, peptide bond, chirality).
6. **Post-filter** to quality-filtered residues using masks and fragment records.

## Output format

Each line of the JSONL file is a JSON object with these fields:

| Field | Type | Description |
|-------|------|-------------|
| `file` | string | Source filename |
| `model` | int | Model number |
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

- Total residues: 97,519
- Unique structures: 492
- Output chains: 514
- Chain list entries: 523 (9 excluded due to all-high-B, 1.7%)
- Compressed size: 2.4 MB

## Regression validation

Phi and psi values were validated against the Richardson Lab's pre-computed angles from the top500-angles archive (2006). Of 36,810 matched phi pairs, 99.7% agree within 1 degree; of 36,809 matched psi pairs, 99.7% agree within 1 degree. The ~100 discrepancies >1 degree per angle are attributable to alternate conformation selection and RCSB coordinate remediation over the 20+ years since the reference data was computed. See `top500_issues.md` for full details.

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
- top2018 full-filtered: doi:10.5281/zenodo.19029946
- top2018 MC-filtered: doi:10.5281/zenodo.19040026
- top8000: doi:10.5281/zenodo.19041232
- Richardson Lab reference_data: https://github.com/rlabduke/reference_data
- Word et al. (1999a) NQH flips: doi:10.1006/jmbi.1998.2401
- Word et al. (1999b) Probe: doi:10.1006/jmbi.1998.2400
- Lovell et al. (2000): PMID:10861930
- Lovell et al. (2003): doi:10.1002/prot.10286
- Davis et al. (2004): doi:10.1093/nar/gkh398
