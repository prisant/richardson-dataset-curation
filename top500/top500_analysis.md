# top500 Dataset: Summary Statistics

**pydangle-biopython version:** 0.5.1 (post-DSSP-fix tree as of 2026-04-26: includes commits `03e5511` ANISOU/SIGATM/SIGUIJ strip, `5e79f71` SEQRES blank-chain fill, `5eaed6b` selective HETATM strip)
**Date:** 2026-04-26 (regenerated after Phase 1 SEQRES-rebuild ersatz; supersedes 2026-03-15 numbers)
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

Percentages computed against records with non-null DSSP value (n = 80,694).
Loop residues are recorded as `dssp = null` in the JSONL (mkdssp emits a
space character for loop; the older `'C'` literal that appeared in
2026-03-15 statistics was a parse_dssp_output artefact, corrected by
commit `03e5511`).

| State | Count | Percent |
|-------|------:|--------:|
| H (alpha helix) | 31,421 | 38.94% |
| E (beta strand) | 22,536 | 27.93% |
| T (turn) | 10,718 | 13.28% |
| S (bend) | 8,296 | 10.28% |
| G (3-10 helix) | 3,954 | 4.90% |
| P (polyproline II) | 1,792 | 2.22% |
| B (beta bridge) | 1,374 | 1.70% |
| I (pi helix) | 603 | 0.75% |
| (loop, dssp = null) | 16,825 | — |

## Geometric Distributions

| Measure | N | Mean | Std |
|---------|--:|-----:|----:|
| phi | 94,594 | -85.589 | 42.247 |
| psi | 94,594 | 13.349 | 115.731 |
| omega | 94,594 | 179.505 | 7.300 |
| tau | 97,519 | 111.000 | 3.190 |

## Data Completeness

| Field | Present | Missing | % Present |
|-------|--------:|--------:|----------:|
| phi | 94,594 | 2,925 | 97.00% |
| psi | 94,594 | 2,925 | 97.00% |
| omega | 94,594 | 2,925 | 97.00% |
| tau | 97,519 | 0 | 100.00% |
| dssp | 80,694 | 16,825 | 82.74% |
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

> **Caveat (2026-04-26):** top500 numbers below are **post-DSSP-fix / post-SEQRES-rebuild** (this dataset). top8000 / top2018 MC numbers are **pre-fix** as published; they will be regenerated under the corrected pipeline (Phase B propagation) before v2.0.0 ships. The H/E/etc. percentages will shift modestly when those datasets are reprocessed (top500's helix percentage rose from 32.4% → 38.9% under the corrected DSSP labeling).

| Metric | top500 (post-fix) | top8000 (pre-fix) | top2018 MC (pre-fix) |
|--------|-------:|--------:|-----------:|
| Residues | 97,519 | 1,568,561 | 2,963,303 |
| Structures | 492 | 7,607 | 13,308 |
| Chains | 514 | 7,814 | 13,677 |
| General % | 71.1% | 71.1% | 70.8% |
| IleVal % | 13.1% | 13.5% | 13.5% |
| Gly % | 7.9% | 7.6% | 7.8% |
| Helix (H) % | **38.9%** | 32.7% (pre-fix) | 32.4% (pre-fix) |
| Strand (E) % | **27.9%** | 23.3% (pre-fix) | 23.0% (pre-fix) |
| phi completeness | 97.0% | 97.2% | 93.5% |

The Ramachandran category distributions remain remarkably consistent across
the three datasets, validating the curation methodology. The DSSP secondary-
structure percentages are temporarily inconsistent because top8000 and
top2018 MC have not yet been reprocessed under the corrected pipeline (the
older 32.x% helix and 23.x% strand reflect the old `'C'`-literal-as-loop-marker
artefact and the silent failures on HETATM-affected structures); cross-
dataset DSSP comparisons should be treated cautiously until Phase B
propagation completes.

## Regression Validation

### Backbone geometry (vs Richardson Lab top500-angles archive, 2006)

Phi and psi values validated against the Richardson Lab top500-angles
archive (dipeptide format with psi on residue i, phi on residue i+1).

| Angle | Matched pairs | < 0.01° | < 1.0° | > 1.0° |
|-------|-------------:|--------:|-------:|-------:|
| phi | 36,810 | 99.0% | 99.7% | 103 (0.3%) |
| psi | 36,809 | 99.0% | 99.7% | 108 (0.3%) |

Large discrepancies concentrated in 1tph (trypsinogen) — alternate
conformation selection and/or RCSB coordinate remediation since 2005.
See `top500_issues.md` for details.

### DSSP labels (vs `pipeline/run_dssp.py` ground truth)

DSSP values in the JSONL validated against independently-generated `.dssp`
files (mkdssp 4.2.2 on cleaned original PDBs via `pipeline/run_dssp.py`)
using `pipeline/validate_dssp.py`. The comparator reuses
`pydangle_biopython.dssp.parse_dssp_output` to keep the validation parser
in sync with the production label code.

| Metric | Value |
|--------|-------|
| JSONL records (post-filter) | 97,519 |
| `.dssp` ground-truth residues | 168,114 (across 494 structures) |
| Agreement (JSONL ↔ .dssp) | **97,511** |
| mkdssp-only (B-filter excluded from JSONL) | 70,603 |
| jsonl-only (mkdssp drops; pydangle keeps) | 8 |
| **Non-trivial mismatches** | **0** |
| Missing-dssp pdbids | 0 |

The full validation report is at `top500/top500_dssp_validation.md`,
regenerable via `python pipeline/validate_dssp.py
top500/top500_measures.jsonl ~/Desktop/Data/top500_pdbs/ -o
top500/top500_dssp_validation.md`.

This 0-mismatch result was achieved 2026-04-26 after Phase 1 propagation
of three DSSP fixes (commits `03e5511`, `5e79f71`, `5eaed6b`) and
rebuilding the top500 ersatz tree to preserve SEQRES records. Prior to
this work the published top500 dataset shipped 18,608 mismatches against
ground truth (most from the parse_dssp_output `'C'` literal artefact;
569 from silent pydangle nulling of structures with HETATM/SEQRES
collisions on 1a1y, 1gdo, 1tax). See `top500_issues.md` and
`tda-project-notes/SessionSummaries/` for the diagnostic chain.

## Software Versions

| Software | Version |
|----------|---------|
| pydangle-biopython | 0.5.1 (post-fix tree: includes commits `03e5511` + `5e79f71` + `5eaed6b`, 2026-04-26) |
| Python | 3.12.3 |
| mkdssp | 4.2.2 |
| Reduce | 4.16.250520 |
| BioPython | (via pip) |
| OS | Ubuntu 24.04, Linux 6.17.0-19-generic |

## Caveats / pending work

- **Figures (`top500_analysis_plots.{ipynb,html,pdf}`) were generated against the pre-fix JSONL** and have not yet been regenerated. They will be refreshed as part of v2.0.0 release prep, when the notebook is re-run against the corrected `top500_measures.jsonl`. The DSSP-related plots in particular will shift (helix percentage 32.4% → 38.9%, etc.). The pre-fix figures are retained as a diagnostic reference for the bug discovery; do not rely on them for current numbers.
- **Cross-dataset comparison** rows for top8000 and top2018 MC are pre-fix and will update during Phase B propagation.
