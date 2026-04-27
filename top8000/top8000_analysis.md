# top8000 Dataset: Summary Statistics

**pydangle-biopython version:** 0.5.1 (post-DSSP-fix tree as of 2026-04-27: includes commits `03e5511` ANISOU/SIGATM/SIGUIJ strip, `5e79f71` SEQRES blank-chain fill, `5eaed6b` selective HETATM strip, `c45effa` mkdssp-overflow returncode handling)
**Date:** 2026-04-27 (regenerated after Phase B/C SEQRES-rebuild ersatz; supersedes 2026-03-15 numbers)
**Quality filter:** Mainchain B-factor ≤ 30 (Williams et al. 2018)
**Measurements:** 15 (backbone geometry, chi1–chi4, Rama classifications, DSSP, peptide bond, chirality)

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

Percentages computed against records with non-null DSSP value (n = 1,310,345).
Loop residues are recorded as `dssp = null` in the JSONL (mkdssp emits a
space character for loop; the older `'C'` literal that appeared in
2026-03-15 statistics was a `parse_dssp_output` artefact, corrected by
commit `03e5511`).

| State | Count | Percent |
|-------|------:|--------:|
| H (alpha helix) | 521,563 | 39.80% |
| E (beta strand) | 374,045 | 28.55% |
| T (turn) | 168,774 | 12.88% |
| S (bend) | 124,942 | 9.54% |
| G (3-10 helix) | 63,039 | 4.81% |
| P (polyproline II) | 29,262 | 2.23% |
| B (beta bridge) | 19,432 | 1.48% |
| I (pi helix) | 9,288 | 0.71% |
| (loop, dssp = null) | 258,216 | — |

## Geometric Distributions

| Measure | N | Mean | Std |
|---------|--:|-----:|----:|
| phi | 1,523,961 | -85.225 | 42.270 |
| psi | 1,523,951 | 14.803 | 115.125 |
| omega | 1,523,961 | 179.273 | 7.902 |
| tau | 1,568,552 | 111.121 | 2.715 |

## Data Completeness

| Field | Present | Missing | % Present |
|-------|--------:|--------:|----------:|
| phi | 1,523,961 | 44,600 | 97.16% |
| psi | 1,523,951 | 44,610 | 97.16% |
| omega | 1,523,961 | 44,600 | 97.16% |
| tau | 1,568,552 | 9 | 100.00% |
| chi1 | 1,309,084 | 259,477 | 83.46% |
| dssp | 1,310,345 | 258,216 | 83.54% |
| rama_category | 1,486,202 | 82,359 | 94.75% |

Missing phi/psi/omega values are due to fragment-boundary nulling
at B-factor filter gaps (44K fragment starts and ends).
The higher completeness compared to top2018 MC (97.2% vs 93.5%)
reflects the simpler B-factor filter producing longer contiguous
fragments than the multi-criteria Richardson Lab filter.

## Regression Validation

### Backbone geometry (vs Richardson Lab Top8000 rotamer reference)

Phi and psi values validated against the Richardson Lab Top8000 rotamer
reference dataset (Hintze et al., 2016; strict RSCC/map/B-factor filter).

| Angle | Matched pairs | < 0.01° | < 1.0° | > 1.0° |
|-------|-------------:|--------:|-------:|-------:|
| phi | 850,224 | 99.94% | 99.95% | 469 (0.06%) |
| psi | 848,096 | 99.90% | 99.93% | 594 (0.07%) |

869,849 residues matched (88.3% of the 985,200 in the reference set).
Large discrepancies (>1°) are dominated by ~180° differences from
alternate conformation selection. See `top8000_issues.md` for details.

### DSSP labels (vs `pipeline/run_dssp.py` ground truth)

DSSP values in the JSONL validated against independently-generated `.dssp`
files (mkdssp 4.2.2 on cleaned original PDBs via `pipeline/run_dssp.py`)
using `pipeline/validate_dssp.py`.

| Metric | Value |
|--------|-------|
| JSONL records (post-filter) | 1,568,561 |
| `.dssp` ground-truth residues | (across 7,737 structures) |
| Agreement (JSONL ↔ .dssp) | **1,567,955** |
| mkdssp-only (B-filter excluded from JSONL) | 2,146,658 |
| jsonl-only (mkdssp drops; pydangle keeps) | 60 |
| **Non-trivial mismatches** | **0** |
| Missing-dssp pdbids | 0 |
| Non-model-1 (alternate models, skipped) | 546 |

Full validation report: `top8000/top8000_dssp_validation.md`.

This 0-mismatch result was achieved 2026-04-27 after Phase B + C
propagation of four DSSP fixes:
- `03e5511` (parse_dssp_output: space → None, ANISOU/SIGATM/SIGUIJ strip)
- `5e79f71` (SEQRES blank-chain fill)
- `5eaed6b` (selective HETATM strip — addresses 1wte, 2gn0)
- `c45effa` (mkdssp returncode-1 with overflow warning — addresses 1h64,
  1ryp, 2zzs, 3mt6 which exceed legacy DSSP header field widths)

Plus rebuild of the top8000 ersatz tree under the post-`8802e0f`
build_ersatz_only.py to preserve SEQRES records.

The previously-published top8000 dataset on Zenodo shipped with
substantial DSSP errors: every structure's loop residues carried the
literal `'C'` instead of `null` (mkdssp parsing artefact), and the 6
HETATM-survey-flagged structures with surviving residues (1h64, 1ryp,
1wte, 2gn0, 2zzs, 3mt6) carried `dssp = null` for *every* residue
because pydangle's per-file mkdssp invocation silently failed (paths
A and B). Approximate scope of silent-null records in the published
top8000:

| Structure | Records (post-filter) | All previously null? |
|---|---:|---|
| 1h64 | 24 | yes |
| 1ryp | 406 | yes |
| 1wte | 526 | yes |
| 2gn0 | 566 | yes |
| 2zzs | 80 | yes |
| 3mt6 | 169 | yes |
| 3cmy | 0 | n/a (B-filtered) |
| **Subtotal** | **1,771** | all 1,771 silently nulled |

Plus ~270k loop residues that carried `'C'` instead of `null`. All
corrected under v2.0.0; see `top8000_issues.md` for the fix history.

## Software Versions

| Software | Version |
|----------|---------|
| pydangle-biopython | 0.5.1 (post-fix tree: includes commits `03e5511` + `5e79f71` + `5eaed6b` + `c45effa`, 2026-04-27) |
| Python | 3.12.3 |
| mkdssp | 4.2.2 |
| Reduce | 4.16.250520 |
| BioPython | (via pip) |
| OS | Ubuntu 24.04, Linux 6.17.0-19-generic |

## Caveats / pending work

- **Figures (`top8000_analysis_plots.{ipynb,html,pdf}`) were generated against the pre-fix JSONL** and have not yet been regenerated. They will be refreshed as part of v2.0.0 release prep when the notebook is re-run against the corrected `top8000_measures.jsonl`. DSSP-related plots will shift (helix percentage 32.66% → 39.80% under proper labeling; the 'C' coil category disappears, reabsorbed into null/loop). The pre-fix figures should not be relied on for current numbers.
- top2018_full and top2018_mc remain pre-fix until their Phase B propagation runs.
