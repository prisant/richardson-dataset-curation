# top8000 Dataset: Summary Statistics

**pydangle-biopython version:** 0.5.1
**Date:** 2026-03-15
**Quality filter:** Mainchain B-factor ≤ 30 (Williams et al. 2018)

## Dataset Overview

| Metric | Value |
|--------|-------|
| Total residues | 1,568,561 |
| Unique structures | 7,607 |
| Output chains | 7,814 |
| Chain list entries | 7,957 |
| Chains excluded (all fail B≤30) | 143 (1.8%) |
| Multi-chain structures | 182 (207 extra chains) |
| Standard amino acids | 1,568,560 |
| Non-standard | 1 (CSD) |

## Ramachandran Category Distribution

| Category | Count | Percent |
|----------|------:|--------:|
| General | 1,056,115 | 71.06% |
| IleVal | 200,462 | 13.49% |
| Gly | 112,610 | 7.58% |
| TransPro | 63,261 | 4.26% |
| PrePro | 50,230 | 3.38% |
| CisPro | 3,524 | 0.24% |

## DSSP Secondary Structure

| State | Count | Percent |
|-------|------:|--------:|
| H (alpha helix) | 509,873 | 32.66% |
| E (beta strand) | 363,619 | 23.29% |
| C (coil) | 275,923 | 17.67% |
| T (turn) | 167,858 | 10.75% |
| S (bend) | 123,456 | 7.91% |
| G (3-10 helix) | 62,206 | 3.98% |
| P (polyproline II) | 28,842 | 1.85% |
| B (beta bridge) | 20,370 | 1.30% |
| I (pi helix) | 9,043 | 0.58% |

## Geometric Distributions

| Measure | N | Mean | Std |
|---------|--:|-----:|----:|
| phi | 1,523,961 | -77.128 | 50.062 |
| psi | 1,523,951 | 39.810 | 89.195 |
| omega | 1,523,961 | 25.966 | 174.286 |
| tau | 1,568,552 | 111.121 | 2.715 |

## Data Completeness

| Field | Present | Missing | % Present |
|-------|--------:|--------:|----------:|
| phi | 1,523,961 | 44,600 | 97.16% |
| psi | 1,523,951 | 44,610 | 97.16% |
| omega | 1,523,961 | 44,600 | 97.16% |
| tau | 1,568,552 | 9 | 100.00% |
| dssp | 1,561,190 | 7,371 | 99.53% |
| rama_category | 1,486,202 | 82,359 | 94.75% |

Missing phi/psi/omega values are due to fragment-boundary nulling
at B-factor filter gaps (44K fragment starts and ends).
The higher completeness compared to top2018 MC (97.2% vs 93.5%)
reflects the simpler B-factor filter producing longer contiguous
fragments than the multi-criteria Richardson Lab filter.

## Software Versions

| Software | Version |
|----------|---------|
| pydangle-biopython | 0.5.1 |
| Python | 3.12.3 |
| mkdssp | 4.2.2 |
| Reduce | 4.16.250520 |
| BioPython | (via pip) |
| OS | Ubuntu 24.04, Linux 6.17.0-19-generic |
