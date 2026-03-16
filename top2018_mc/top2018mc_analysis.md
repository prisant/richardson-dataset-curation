# top2018 MC-Filtered Dataset: Summary Statistics

**pydangle-biopython version:** 0.5.1
**Date:** 2026-03-15

## Dataset Overview

| Metric | Value |
|--------|-------|
| Total residues | 2,963,303 |
| Unique structures | 13,308 |
| Quality-filtered chains | 13,677 |
| Multi-chain structures | 329 (369 extra chains) |
| Standard amino acids | 2,963,303 |
| Non-standard mapped | MSE → MET (2) |

The source chain list (`top2018_chains_hom70_mcfilter_60pct_complete.txt`)
contains **13,677 chain entries** from **13,308 unique PDB structures**. Of
these, 329 structures contribute more than one quality-filtered chain.
Each structure is processed as a single ersatz PDB file containing all
chains; the output is then post-filtered to include only residues from
the quality-filtered chains using per-chain mask files.

## Ramachandran Category Distribution

| Category | Count | Percent |
|----------|------:|--------:|
| General | 1,841,337 | 70.75% |
| IleVal | 350,387 | 13.46% |
| Gly | 202,952 | 7.80% |
| TransPro | 112,742 | 4.33% |
| PrePro | 89,121 | 3.42% |
| CisPro | 5,949 | 0.23% |

## DSSP Secondary Structure

| State | Count | Percent |
|-------|------:|--------:|
| H (alpha helix) | 952,968 | 32.38% |
| E (beta strand) | 677,384 | 23.02% |
| C (coil) | 520,455 | 17.69% |
| T (turn) | 323,358 | 10.99% |
| S (bend) | 237,704 | 8.08% |
| G (3-10 helix) | 117,542 | 3.99% |
| P (polyproline II) | 57,622 | 1.96% |
| B (beta bridge) | 38,209 | 1.30% |
| I (pi helix) | 17,515 | 0.60% |

## Geometric Distributions

| Measure | N | Mean | Std |
|---------|--:|-----:|----:|
| phi | 2,771,857 | -76.891 | 50.295 |
| psi | 2,771,715 | 39.892 | 89.074 |
| omega | 2,771,861 | 26.061 | 173.832 |
| tau | 2,963,289 | 111.140 | 2.506 |

## Data Completeness

| Field | Present | Missing | % Present |
|-------|--------:|--------:|----------:|
| phi | 2,771,857 | 191,446 | 93.54% |
| psi | 2,771,715 | 191,588 | 93.53% |
| omega | 2,771,861 | 191,442 | 93.54% |
| tau | 2,963,289 | 14 | 100.00% |
| dssp | 2,942,757 | 20,546 | 99.31% |
| rama_category | 2,602,488 | 360,815 | 87.82% |

Missing phi/psi/omega values are due to fragment-boundary nulling (191K
fragment starts lack a valid i-1 neighbor; 191K fragment ends lack a
valid i+1 neighbor). Missing rama_category values include all residues
where either phi or psi is null (both are needed for classification).

## Residue Composition (top 20)

| Residue | Count | Percent |
|---------|------:|--------:|
| LEU | 278,872 | 9.41% |
| ALA | 262,777 | 8.87% |
| GLY | 232,691 | 7.85% |
| VAL | 223,447 | 7.54% |
| GLU | 178,620 | 6.03% |
| ILE | 175,347 | 5.92% |
| ASP | 171,682 | 5.79% |
| THR | 167,699 | 5.66% |
| SER | 163,093 | 5.50% |
| LYS | 155,314 | 5.24% |
| ARG | 143,215 | 4.83% |
| PRO | 139,649 | 4.71% |
| PHE | 126,640 | 4.27% |
| ASN | 124,440 | 4.20% |
| TYR | 112,200 | 3.79% |
| GLN | 103,773 | 3.50% |
| HIS | 70,358 | 2.37% |
| MET | 50,790 | 1.71% |
| TRP | 46,657 | 1.57% |
| CYS | 36,039 | 1.22% |

## Known Issues

4 chains out of 13,677 (0.029%) have RCSB remediation mismatches
affecting 47 mask residues total (0.0016% of dataset). These residues
are silently excluded during post-filtering. See `top2018mc_issues.md`
for details.

## Software Versions

| Software | Version |
|----------|---------|
| pydangle-biopython | 0.5.1 |
| Python | 3.12.3 |
| mkdssp | 4.2.2 |
| Reduce | 4.16.250520 |
| BioPython | (via pip) |
| OS | Ubuntu 24.04, Linux 6.17.0-19-generic |
