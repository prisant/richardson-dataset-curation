# top500 Dataset: Summary Statistics

**pydangle-biopython version:** 0.5.1
**Date:** 2026-03-15
**Quality filter:** Mainchain B-factor ≤ 30 (Lovell et al. 2003; Williams et al. 2018)

## Dataset Overview

| Metric | Value |
|--------|-------|
| Total residues | 97,519 |
| Unique structures | 492 |
| Output chains | 514 |
| Chain list entries | 523 |
| Chains excluded (all fail B≤30) | 9 (1.7%) |
| Multi-chain structures | 21 (22 extra chains) |
| Standard amino acids | 97,518 |
| Non-standard | 1 (CSD) |

## Ramachandran Category Distribution

| Category | Count | Percent |
|----------|------:|--------:|
| General | 65,472 | 71.11% |
| IleVal | 12,041 | 13.08% |
| Gly | 7,275 | 7.90% |
| TransPro | 3,927 | 4.26% |
| PrePro | 3,143 | 3.41% |
| CisPro | 219 | 0.24% |

## DSSP Secondary Structure

| State | Count | Percent |
|-------|------:|--------:|
| H (alpha helix) | 30,756 | 32.41% |
| E (beta strand) | 21,745 | 22.92% |
| C (coil) | 16,462 | 17.35% |
| T (turn) | 10,357 | 10.91% |
| S (bend) | 8,042 | 8.47% |
| G (3-10 helix) | 3,872 | 4.08% |
| P (polyproline II) | 1,743 | 1.84% |
| B (beta bridge) | 1,322 | 1.39% |
| I (pi helix) | 592 | 0.62% |

## Geometric Distributions

| Measure | N | Mean | Std |
|---------|--:|-----:|----:|
| phi | 94,594 | -76.985 | 50.991 |
| psi | 94,594 | 40.569 | 89.442 |
| omega | 94,594 | 24.881 | 175.447 |
| tau | 97,519 | 111.000 | 3.190 |

## Data Completeness

| Field | Present | Missing | % Present |
|-------|--------:|--------:|----------:|
| phi | 94,594 | 2,925 | 97.00% |
| psi | 94,594 | 2,925 | 97.00% |
| omega | 94,594 | 2,925 | 97.00% |
| tau | 97,519 | 0 | 100.00% |
| dssp | 94,891 | 2,628 | 97.31% |
| rama_category | 92,077 | 5,442 | 94.42% |

## Residue Composition

| Residue | Count | Percent |
|---------|------:|--------:|
| ALA | 8,600 | 8.82% |
| LEU | 8,421 | 8.64% |
| GLY | 7,704 | 7.90% |
| VAL | 7,243 | 7.43% |
| SER | 5,811 | 5.96% |
| THR | 5,711 | 5.86% |
| ASP | 5,643 | 5.79% |
| GLU | 5,533 | 5.67% |
| ILE | 5,415 | 5.55% |
| LYS | 5,380 | 5.52% |
| ASN | 4,431 | 4.54% |
| PRO | 4,429 | 4.54% |
| ARG | 4,349 | 4.46% |
| PHE | 4,096 | 4.20% |
| TYR | 3,654 | 3.75% |
| GLN | 3,601 | 3.69% |
| HIS | 2,316 | 2.37% |
| MET | 2,004 | 2.06% |
| CYS | 1,645 | 1.69% |
| TRP | 1,532 | 1.57% |

## Comparison Across Datasets

| Metric | top500 | top8000 | top2018 MC |
|--------|-------:|--------:|-----------:|
| Residues | 97,519 | 1,568,561 | 2,963,303 |
| Structures | 492 | 7,607 | 13,308 |
| Chains | 514 | 7,814 | 13,677 |
| General % | 71.1% | 71.1% | 70.8% |
| IleVal % | 13.1% | 13.5% | 13.5% |
| Gly % | 7.9% | 7.6% | 7.8% |
| Helix (H) % | 32.4% | 32.7% | 32.4% |
| Strand (E) % | 22.9% | 23.3% | 23.0% |
| phi completeness | 97.0% | 97.2% | 93.5% |

The remarkable consistency of Ramachandran category and DSSP distributions
across three datasets spanning 20 years and two orders of magnitude in size
validates both the curation methodology and the measurement pipeline.

## Regression Validation

Phi and psi values validated against the Richardson Lab top500-angles
archive (2006; dipeptide format with psi on residue i, phi on residue i+1).

| Angle | Matched pairs | < 0.01° | < 1.0° | > 1.0° |
|-------|-------------:|--------:|-------:|-------:|
| phi | 36,810 | 99.0% | 99.7% | 103 (0.3%) |
| psi | 36,809 | 99.0% | 99.7% | 108 (0.3%) |

Large discrepancies concentrated in 1tph (trypsinogen) — alternate
conformation selection and/or RCSB coordinate remediation since 2005.
See `top500_issues.md` for details.

## Software Versions

| Software | Version |
|----------|---------|
| pydangle-biopython | 0.5.1 |
| Python | 3.12.3 |
| mkdssp | 4.2.2 |
| Reduce | 4.16.250520 |
| BioPython | (via pip) |
| OS | Ubuntu 24.04, Linux 6.17.0-19-generic |
