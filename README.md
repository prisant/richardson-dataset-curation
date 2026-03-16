# Richardson Dataset Curation

Curated per-residue geometric measurements and classifications for
high-quality protein structure datasets from the Richardson Laboratory
at Duke University.

These datasets are computed using
[pydangle-biopython](https://github.com/prisant/pydangle-biopython)
on "ersatz" full-structure PDB files that preserve multi-chain context
for correct DSSP secondary structure assignments while incorporating
NQH sidechain flip corrections from Reduce.

## Datasets

| Dataset | Source | Status | Zenodo DOI |
|---------|--------|--------|------------|
| [top2018_full](top2018_full/) | Full-filtered, 70% homology | Complete | [10.5281/zenodo.19029946](https://doi.org/10.5281/zenodo.19029946) |
| [top2018_mc](top2018_mc/) | Mainchain-filtered, 70% homology | Complete | [10.5281/zenodo.19040026](https://doi.org/10.5281/zenodo.19040026) |
| [top8000](top8000/) | Top 8000 chains | Complete | [10.5281/zenodo.19049701](https://doi.org/10.5281/zenodo.19049701) (v1.1.0) |
| [top500](top500/) | Top 500 chains | Complete | [10.5281/zenodo.19042131](https://doi.org/10.5281/zenodo.19042131) |
| [top100](top100/) | Top 100 chains | Planned | — |

## Pipeline overview

Each dataset is produced by the same pipeline:

1. **Download full PDB structures** from RCSB for each chain in the
   Richardson Lab source dataset.
2. **Run Reduce** on each full PDB to add hydrogens and optimize NQH
   sidechain orientations using three-point hinge rotation (matching the
   published flip method).
3. **Build ersatz PDB files** by stripping hydrogens and problematic
   headers from the Reduce output.  This preserves multi-chain context
   with correct NQH corrections.
4. **Run pydangle-biopython** on all ersatz files to compute backbone
   geometry, sidechain torsions, Ramachandran classifications, and DSSP
   secondary structure.
5. **Post-filter** to quality-filtered residues using mask files derived
   from the Richardson Lab pruned files.  Fragment boundaries from
   `USER INC` records ensure no measurement depends on coordinates from
   quality-filtered residues.

The shared pipeline scripts are in [`pipeline/`](pipeline/).  Each
dataset directory contains a `reproduce.sh` script for end-to-end
reproduction and a `README.md` with dataset-specific documentation.

## Why ersatz files?

The Richardson Lab datasets provide pruned single-chain PDB files with
hydrogens and NQH flips applied by Reduce.  Running DSSP on these
isolated chains gives incorrect secondary structure assignments because
inter-chain hydrogen bond context is lost.  The ersatz files solve this
by providing the complete multi-chain structure (for correct DSSP) with
NQH corrections applied to all chains (for correct sidechain geometry),
stripped of hydrogens (for mkdssp compatibility).

## Supplementary PDB files

Three PDB files (4v4m, 4yza, 4ztt) are no longer available from RCSB in
PDB format.  They are included in [`supplementary_pdbs/`](supplementary_pdbs/)
to ensure full reproducibility.

## References

- Williams, C. J., Richardson, D. C., & Richardson, J. S. (2021). The
  importance of residue-level filtering, and the Top2018 best-parts
  dataset of high-quality protein residues. *Protein Science*.
  doi:10.1002/pro.4239

- Williams, C. J., Headd, J. J., Moriarty, N. W., et al. (2018).
  MolProbity: More and better reference data for improved all-atom
  structure validation. *Protein Science*, 27(1), 293–315.
  doi:10.1002/pro.3330

## Author

Michael G. Prisant, Ph.D.
Prisant Scientific
info@prisantscientific.org

## License

Code: MIT
Data products: CC-BY 4.0
