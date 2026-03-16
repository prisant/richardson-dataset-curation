# top100 Dataset: Summary Statistics

**pydangle-biopython version:** 0.5.1
**Date:** 2026-03-16
**Quality filter:** Mainchain B-factor ≤ 40 (Word et al. 1999; see Issues for rationale)
**Measurements:** 15 (backbone geometry, chi1–chi4, Rama classifications, DSSP, peptide bond, chirality)

## Dataset Overview

| Metric | Value |
|--------|-------|
| Total residues | 17,434 |
| Unique structures | 98 |
| Output chains | 106 |
| Chain list entries | 108 |
| Chains excluded (all fail B≤40) | 2 (1.9%) |
| Multi-chain structures | 8 (8 extra chains) |
| Standard amino acids | 17,434 |
| Non-standard | 0 |

## Ramachandran Category Distribution

| Category | Count | Percent |
|----------|------:|--------:|
| General | 12,187 | 71.63% |
| IleVal | 2,052 | 12.06% |
| Gly | 1,437 | 8.45% |
| TransPro | 720 | 4.23% |
| PrePro | 567 | 3.33% |
| CisPro | 50 | 0.29% |

## DSSP Secondary Structure

| State | Count | Percent |
|-------|------:|--------:|
| H (alpha helix) | 4,964 | 29.15% |
| E (beta strand) | 3,884 | 22.81% |
| C (coil) | 3,115 | 18.29% |
| T (turn) | 2,093 | 12.29% |
| S (bend) | 1,591 | 9.34% |
| G (3-10 helix) | 686 | 4.03% |
| B (beta bridge) | 316 | 1.86% |
| P (polyproline II) | 303 | 1.78% |
| I (pi helix) | 77 | 0.45% |

## Geometric Distributions

| Measure | N | Mean | Std |
|---------|--:|-----:|----:|
| phi | 17,218 | -77.138 | 52.745 |
| psi | 17,218 | 43.801 | 90.295 |
| omega | 17,218 | 19.502 | 175.810 |
| tau | 17,434 | 110.729 | 3.269 |

## Data Completeness

| Field | Present | Missing | % Present |
|-------|--------:|--------:|----------:|
| phi | 17,218 | 216 | 98.76% |
| psi | 17,218 | 216 | 98.76% |
| omega | 17,218 | 216 | 98.76% |
| tau | 17,434 | 0 | 100.00% |
| chi1 | 14,380 | 3,054 | 82.48% |
| dssp | 17,029 | 405 | 97.68% |
| rama_category | 17,013 | 421 | 97.59% |

The higher phi/psi completeness (98.8%) compared to top500/top8000
(97.0–97.2%) reflects the looser B ≤ 40 threshold producing fewer
fragment boundaries (214 vs 399 with B ≤ 30 on the same data).

## Comparison Across All Five Datasets

| Metric | top100 | top500 | top8000 | top2018 MC | top2018 full |
|--------|-------:|-------:|--------:|-----------:|-----------:|
| Residues | 17,434 | 97,519 | 1,568,561 | 2,963,303 | 2,405,322 |
| Structures | 98 | 492 | 7,607 | 13,308 | 11,843 |
| B-factor threshold | ≤ 40 | ≤ 30 | ≤ 30 | pruned | pruned |
| General % | 71.6% | 71.1% | 71.1% | 70.8% | 68.6% |
| IleVal % | 12.1% | 13.1% | 13.5% | 13.5% | 14.0% |
| Gly % | 8.5% | 7.9% | 7.6% | 7.8% | 8.8% |
| Helix (H) % | 29.2% | 32.4% | 32.7% | 32.4% | 31.4% |
| Strand (E) % | 22.8% | 22.9% | 23.3% | 23.0% | 23.5% |

The consistency of Ramachandran and DSSP distributions across five
datasets spanning 25 years and three orders of magnitude in size
validates both the curation methodology and the measurement pipeline.

## Software Versions

| Software | Version |
|----------|---------|
| pydangle-biopython | 0.5.1 |
| Python | 3.12.3 |
| mkdssp | 4.2.2 |
| Reduce | 4.16.250520 |
| BioPython | (via pip) |
| OS | Ubuntu 24.04, Linux 6.17.0-19-generic |
