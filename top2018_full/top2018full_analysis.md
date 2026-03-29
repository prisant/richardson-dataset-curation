# top2018 Full-Filtered Dataset: Summary Statistics

**Zenodo DOI:** [10.5281/zenodo.19029946](https://doi.org/10.5281/zenodo.19029946) (all versions)
**pydangle-biopython version:** 0.5.1
**Date:** 2026-03-15

## Dataset Overview

| Metric | Value |
|--------|-------|
| Total residues | 2,405,322 |
| Unique structures | 11,843 |
| Quality-filtered chains | 12,125 \[1\] |
| Output chain IDs | 12,126 \[1\] |
| Multi-chain structures | 256 (283 extra chains) |
| Standard amino acids | 2,405,322 |
| Non-standard mapped | MSE → MET (2) |

The source chain list (`top2018_chains_hom70_fullfiltered_60pct_complete.txt`)
contains **12,125 chain entries** from **11,843 unique PDB structures**. Of
these, 256 structures contribute more than one quality-filtered chain
(e.g., 5t5i contributes chains A, C, G, J, L). Each structure is processed
as a single ersatz PDB file containing all chains; the output is then
post-filtered to include only residues from the quality-filtered chains
using per-chain mask files.

\[1\] The output contains 12,126 unique (file, chain) combinations rather
than 12,125 because RCSB PDB remediation split the original chain B of
entry 1xmz into RCSB chains B (residues ≤62) and C (residues ≥66). Two
additional entries required remapping: 4ct3 chain D was renamed to RCSB
chain E, and 5tvo chain A residues were renumbered with a +315 offset.
Three entries (5kve\_L, 4oe9\_B, 3wwL\_A) have 7 total residues absent
from the current RCSB structures due to remediation (0.0003% of the
dataset). See [top2018full_issues.md](top2018full_issues.md) §8 for full details.

## Residue Composition

| Residue | Count | Percent |
|---------|------:|--------:|
| ALA | 244,156 | 10.15% |
| LEU | 233,506 | 9.71% |
| GLY | 217,173 | 9.03% |
| VAL | 193,740 | 8.05% |
| THR | 146,896 | 6.11% |
| SER | 145,239 | 6.04% |
| ILE | 145,068 | 6.03% |
| ASP | 137,869 | 5.73% |
| PRO | 123,415 | 5.13% |
| PHE | 109,667 | 4.56% |
| GLU | 108,527 | 4.51% |
| ASN | 103,144 | 4.29% |
| TYR | 94,603 | 3.93% |
| ARG | 92,125 | 3.83% |
| LYS | 72,286 | 3.00% |
| GLN | 71,944 | 2.99% |
| HIS | 57,298 | 2.38% |
| TRP | 39,555 | 1.65% |
| MET | 37,283 | 1.55% |
| CYS | 31,828 | 1.32% |

## Ramachandran Categories

| Category | Count | Percent |
|----------|------:|--------:|
| General | 1,185,578 | 68.56% |
| IleVal | 242,781 | 14.04% |
| Gly | 151,712 | 8.77% |
| TransPro | 80,908 | 4.68% |
| PrePro | 63,787 | 3.69% |
| CisPro | 4,574 | 0.26% |

## Cis Peptide Bonds

| Category | Cis | Total | Percent |
|----------|----:|------:|--------:|
| Overall | 5,902 | 2,027,346 | 0.291% |
| Proline | 5,154 | 103,564 | 4.977% |
| Non-proline | 748 | 1,923,782 | 0.039% |

## Chirality

| Metric | Value |
|--------|-------|
| L-amino acids | 2,187,890 |
| D-amino acids | 0 |

## DSSP Secondary Structure

| State | Count | Percent | Description |
|-------|------:|--------:|-------------|
| H | 750,497 | 31.39% | Alpha helix |
| E | 562,254 | 23.52% | Beta strand |
| C | 433,516 | 18.13% | Coil |
| T | 261,395 | 10.93% | Turn |
| S | 195,507 | 8.18% | Bend |
| G | 94,499 | 3.95% | 3-10 helix |
| P | 47,117 | 1.97% | Polyproline II |
| B | 32,210 | 1.35% | Beta bridge |
| I | 14,697 | 0.61% | Pi helix |

## Geometric Distributions (Overall)

| Measure | N | Mean | Std |
|---------|--:|-----:|----:|
| phi | 2,027,343 | -85.842 | 42.650 |
| psi | 2,027,390 | 21.398 | 117.893 |
| omega | 2,027,346 | 179.186 | 8.315 |
| tau | 2,405,312 | 111.152 | 2.530 |
| chi1 | 1,936,652 | -85.699 | 87.319 |

## Data Completeness

| Field | Present | Missing | % Present |
|-------|--------:|--------:|----------:|
| phi | 2,027,343 | 377,979 | 84.29% |
| psi | 2,027,390 | 377,932 | 84.29% |
| omega | 2,027,346 | 377,976 | 84.29% |
| tau | 2,405,312 | 10 | 100.00% |
| chi1 | 1,936,652 | 468,670 | 80.52% |
| is_cis | 2,027,346 | 377,976 | 84.29% |
| is_left | 2,187,890 | 217,432 | 90.96% |
| dssp | 2,390,692 | 14,630 | 99.39% |
| rama_category | 1,729,340 | 675,982 | 71.90% |

Missing values for phi/psi/omega/is_cis reflect fragment boundaries where
the neighboring residue was quality-filtered and its coordinates are not
used (see [top2018full_issues.md](top2018full_issues.md) §9). Missing chi1 values are for glycine
(no sidechain) and residues with incomplete sidechain atoms. Missing
rama_category values occur where either phi or psi is unavailable. The 10
missing tau values and 14,630 missing dssp values are from structures
where mkdssp failed or residues fell outside polypeptide fragments.
