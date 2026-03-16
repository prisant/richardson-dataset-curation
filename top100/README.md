# pydangle top100: Per-residue backbone geometry for 106 quality-filtered protein chains

**Zenodo DOI:** [10.5281/zenodo.19050466](https://doi.org/10.5281/zenodo.19050466)
**Source:** Richardson Lab Top100 dataset — the original quality-filtered reference

## Summary

This dataset contains per-residue backbone geometric measurements, Ramachandran classifications at multiple granularities, DSSP secondary structure assignments, peptide bond classification, and chirality for 16,635 quality-filtered protein residues from 106 chains in 98 protein structures.

The measurements were computed using [pydangle-biopython](https://github.com/prisant/pydangle-biopython) v0.5.1, a Python reimagining of the Richardson Lab's Java Dangle tool, with BioPython as the structure-parsing backend.

## Source dataset

The Top100 is the original Richardson Lab quality-filtered protein structure reference dataset, curated at Duke University. It consists of 100 high-resolution (0.9–1.7 Å) protein structures selected for crystallographic quality. The Top100 was used for the foundational work on Asn/Gln sidechain amide flip correction using all-atom contact analysis with Probe and Reduce (Word et al., 1999a, 1999b). It was the first dataset to establish the methodology of adding explicit hydrogen atoms and using all-atom contacts to optimize sidechain orientations — a methodology that became central to the MolProbity validation pipeline.

The Top100 was subsequently expanded to the Top500 (Lovell et al., 2000, 2003), which was used for the "penultimate" rotamer library and updated Ramachandran/Cbeta deviation validation. The Top500 was in turn superseded by the Top8000 (Hintze et al., 2016; Williams et al., 2018) and ultimately the Top2018 (Williams et al., 2021).

### Key references

- Word, J. M., Lovell, S. C., Richardson, J. S., & Richardson, D. C. (1999a). Asparagine and glutamine: using hydrogen atom contacts in the choice of side-chain amide orientation. *J. Mol. Biol.*, 285(4), 1735–1747. doi:10.1006/jmbi.1998.2401
- Word, J. M., Lovell, S. C., Richardson, J. S., & Richardson, D. C. (1999b). Visualizing and quantifying molecular complementarity. *J. Mol. Biol.*, 285(4), 1711–1733. doi:10.1006/jmbi.1998.2400
- Lovell, S. C., Word, J. M., Richardson, J. S., & Richardson, D. C. (2000). The penultimate rotamer library. *Proteins*, 40(3), 389–408. PMID:10861930
- Lovell, S. C., Davis, I. W., Arendall, W. B. III, de Bakker, P. I. W., Word, J. M., Prisant, M. G., Richardson, J. S., & Richardson, D. C. (2003). Structure validation by Calpha geometry: phi,psi and Cbeta deviation. *Proteins*, 50(3), 437–450. doi:10.1002/prot.10286
- Williams, C. J., Headd, J. J., Moriarty, N. W., Prisant, M. G., et al. (2018). MolProbity: More and better reference data for improved all-atom structure validation. *Protein Science*, 27(1), 293–315. doi:10.1002/pro.3330

### Chain list provenance

The Top100 chain list was reconstructed from the kinemage webserver archive (top100H directory, 100 FH files processed with Reduce v2.13.2). The original FH files are preserved in the archive. For this dataset, full multi-chain PDB files were downloaded from RCSB and processed with Reduce 4.16 to provide correct multi-chain context for DSSP.

## Residue-level quality filter

Residues are filtered using a mainchain B-factor ≤ 40 criterion, matching the original Top100 methodology (Word et al., 1999a). This is intentionally less strict than the B ≤ 30 used for the larger top500 and top8000 datasets — the original Top100 used the higher threshold because with only 100 structures, a stricter cutoff would have sacrificed too many residues for adequate statistics. Of 108 chain entries, 2 (1.9%) have zero passing residues even at B ≤ 40.

## Preprocessing: "Ersatz" full-structure PDB files

### Pipeline

1. **Full PDB structures** downloaded from RCSB (100 structures).
2. **Run Reduce 4.16.250520** for NQH sidechain flip corrections.
3. **Build ersatz PDB files** by stripping hydrogens and problematic headers.
4. **Build B-factor masks** with synthetic fragment records.
5. **Run pydangle-biopython** with 15 measurements (backbone geometry, sidechain torsions chi1–chi4, Ramachandran classifications, DSSP, peptide bond type, chirality).
6. **Post-filter** to quality-filtered residues.

## Output format

Same JSONL format as all sibling datasets. Fields: file, model, chain, resnum, ins, resname, phi, psi, omega, tau, chi1, chi2, chi3, chi4, rama_category, rama5, rama4, rama3, dssp, peptide_bond, chirality.

Sidechain torsion angles (chi1–chi4) are included for residues that have them. Gly and Ala have no chi angles; shorter sidechains have fewer (e.g., Val has only chi1). This enables the dataset to be used for rotamer analysis as well as Ramachandran analysis, consistent with the Top100's historical role in developing the penultimate rotamer library.

## Dataset statistics

- Total residues: 17,434
- Unique structures: 98
- Output chains: 106
- Chain list entries: 108 (2 excluded at B ≤ 40, 1.9%)
- Compressed size: ~430 KB

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
- top8000: doi:10.5281/zenodo.19049701
- top500: doi:10.5281/zenodo.19049923
- Word et al. (1999a) NQH flips: doi:10.1006/jmbi.1998.2401
- Word et al. (1999b) Probe: doi:10.1006/jmbi.1998.2400
- Lovell et al. (2000): PMID:10861930
- Lovell et al. (2003): doi:10.1002/prot.10286
- Williams et al. (2018): doi:10.1002/pro.3330
