# top8000 — Known Issues

## 143 Chains with Zero Passing Residues (1.8%)

Of the 7,957 chains in the Top8000-best_hom70 chain list, 143 chains
(1.8%) have zero residues passing the mainchain B-factor ≤ 30 filter.
These chains are present in the current RCSB PDB files with the correct
chain IDs (no remediation issues), but every residue has at least one
mainchain atom (N, CA, C, or O) with B-factor exceeding 30.

### Breakdown by minimum mainchain B-factor

| Min MC B-factor range | Chains | Notes |
|----------------------|--------|-------|
| < 25 | 6 | Some good atoms, but every residue has ≥1 atom > 30 |
| 25–28 | 15 | Most atoms near threshold |
| 28–30 | 41 | Borderline — best atoms just under 30 |
| 30–35 | 53 | Uniformly above threshold |
| > 35 | 28 | Clearly disordered throughout |

62 of the 143 chains (43%) have their lowest mainchain B-factor below 30,
meaning some individual atoms are well-ordered but every residue has at
least one mainchain atom above the threshold. Testing with N, CA, C only
(excluding O) recovers only 16 additional residues across 13 chains,
confirming that the issue is not specific to the more mobile O atom.

### Interpretation

These chains passed the Top8000 chain-level quality filter (MolProbity
score < 2.0, resolution < 2.0 Å, ≤5% geometry outliers) but fail the
residue-level B-factor filter uniformly. This represents a gap between
chain-level and residue-level quality criteria: a chain can have good
overall geometry while still having elevated B-factors throughout,
particularly for mobile domains, flexible loops, or chains at crystal
contacts with high thermal motion.

The Williams et al. (2018) paper does not discuss this chain-level
dropout explicitly, likely because the Ramachandran contour generation
pooled all passing residues across the entire dataset and did not
track per-chain contributions.

### Impact

- 143 chains excluded from 7,957 total (1.8%)
- 7,814 chains retained with ≥1 passing residue
- 7,607 structures contribute residues (130 structures have all chains excluded)
- The B ≤ 30 threshold matches the published methodology exactly

## Regression Validation Against Richardson Lab Reference Values

The phi and psi values in this dataset were validated against the
Richardson Lab's pre-computed angle data from the Top8000 rotamer
reference dataset (Top8000_filtered_residues_rotarama, available at
https://github.com/rlabduke/reference_data). That dataset uses the
stricter rotamer-level filtering (RSCC ≥ 0.7, map σ ≥ 1.1, B < 40)
from Hintze et al. (2016), so the matched residue set is a subset of
ours.

### Match statistics

| Metric | Value |
|--------|-------|
| Richardson Lab reference residues | 985,200 |
| Pydangle top8000 residues | 1,567,394 |
| Matched (both datasets) | 869,849 (88.3% of reference) |
| Reference-only | 115,351 |

Of the 115,351 reference-only residues: 55,014 are from 449 structures
in the SFbest chain list that are not in our best_hom70 list (different
chain lists); 60,337 are from structures we have but the residue did not
appear in our output (different filtering criteria or fragment boundary
definitions).

### Phi comparison (850,224 pairs with both values non-null)

| Threshold | Count | Percent |
|-----------|------:|--------:|
| \|diff\| < 0.01° | 849,696 | 99.94% |
| \|diff\| < 0.1° | 849,704 | 99.94% |
| \|diff\| < 1.0° | 849,755 | 99.95% |
| \|diff\| > 1.0° | 469 | 0.06% |

Mean |diff|: 0.0209°. Median |diff|: 0.0000°. Max |diff|: 178.47°.

### Psi comparison (848,096 pairs with both values non-null)

| Threshold | Count | Percent |
|-----------|------:|--------:|
| \|diff\| < 0.01° | 847,273 | 99.90% |
| \|diff\| < 0.1° | 847,295 | 99.91% |
| \|diff\| < 1.0° | 847,502 | 99.93% |
| \|diff\| > 1.0° | 594 | 0.07% |

Mean |diff|: 0.0280°. Median |diff|: 0.0000°. Max |diff|: 179.86°.

### Interpretation of discrepancies

The ~500 large discrepancies (>1°) for each angle are dominated by
differences of ~180°, the signature of alternate conformation selection.
The Richardson Lab data and pydangle may select different alternate
conformations for residues with multiple modeled positions. This is
expected behavior, not an error.

The 19,625 phi null mismatches and 21,753 psi null mismatches reflect
fragment boundary differences between our B-factor-derived fragments
and the Richardson Lab's fragment definitions.

### Conclusion

The pipeline produces values that are identical to the Richardson Lab
reference to within numerical precision (99.9% within 0.01°) for the
matched residue set. The small number of large discrepancies are
attributable to alternate conformation selection differences.

## RCSB Remediation

No chain ID or residue numbering mismatches were detected. All 7,957
chain IDs from the 2011 chain list exist in the current RCSB PDB files.
This is notable given that the chain list is 15 years old — RCSB
remediation has not affected chain IDs for any structure in this dataset.

## Reduce Processing

All 7,737 structures processed successfully with Reduce 4.16.250520.
Zero errors, zero timeouts.

## Non-standard Residues

One non-standard residue (CSD, 3-sulfinoalanine) in the filtered output.
No MSE→MET mapping was needed (unlike top2018 which had 2 MSE residues).

## Filtering Methodology Difference from top2018

The top2018 datasets use Richardson Lab pruned PDB files with explicit
USER INC fragment records from a multi-criteria quality filter (including
clash analysis, real-space correlation, and local map values). The top8000
dataset predates this infrastructure and uses a simpler B-factor ≤ 30
mainchain filter as described in Williams et al. (2018). Fragment
boundaries in top8000 are inferred from gaps in the B-factor-passing
residue sequence rather than from explicit USER INC records.

## DSSP File Generation

DSSP output files (`.dssp`) were generated for all 7,737 structures
using `pipeline/run_dssp.py` with mkdssp 4.2.2. The script builds a
minimal cleaned PDB (HEADER + CRYST1 + SEQRES + ATOM + TER +
protein-only HETATM) before invoking mkdssp to work around four
mkdssp 4.x issues (also encountered in top500 and top100):

1. **ANISOU records** cause mkdssp to produce incomplete residue
   assignments.
2. **SIGATM/SIGUIJ records** cause mkdssp to produce zero residues.
3. **Blank chain IDs** (column 22 of ATOM/TER, column 12 of SEQRES)
   in old depositions cause parse failures or silent zero-residue
   output. Top8000 has 0 ersatz files with blank-chain SEQRES (verified
   by survey 2026-04-26).
4. **Non-protein HETATM records** (water, sulfate, glycerol, ions,
   sugars, ...) collide with the SEQRES protein alignment and cause
   mkdssp to produce zero residues for the entire structure.
   **Affected top8000 structures by HETATM-survey hard-fail (path-A):**
   1h64, 1ryp, 1wte, 2gn0, 2zzs, 3cmy, 3mt6 (7 of 7,737). HETATM
   lines whose residue name is in the file's SEQRES or ATOM record
   set (i.e., modified amino-acid residues like MSE, PCA, MEN, PTR,
   LOV that some old PDBs deposit as HETATM) are kept; other HETATM
   are stripped.

Additionally, **mkdssp 4.x exits with returncode 1** (with the
warning "This file contains data that won't fit in the original
DSSP format") when a structure's HEADER / COMPND fields exceed
legacy fixed-width column widths. The residue-by-residue secondary-
structure output is correct in this case; pydangle's wrapper
(`pydangle.dssp.run_dssp()`) trusts stdout when it contains the
residue-table marker (commit `c45effa`, see fix history below).
**Affected top8000 structures (returncode-1):** 1h64, 1ryp, 2zzs,
3mt6 (4 of the 7 originally-flagged HETATM hard-fails; 1wte / 2gn0
work via the HETATM strip alone; 3cmy is fully B-filtered).

SEQRES records are included so mkdssp correctly identifies modified
amino acids (MEN, PTR, PCA, LOV, MSE, etc.) as part of the
polypeptide chain. The `top8000_pdbs_hom70/*/<pdbid>/_ersatz.pdb`
files were rebuilt 2026-04-26 with `build_ersatz_only.py` post-
commit `8802e0f` to include SEQRES; previously they did not.

All 7,737 structures produce `.dssp` files with zero hard errors
under the post-fix pipeline.

## DSSP Semantics

The `dssp` column in the JSONL output contains the mkdssp secondary
structure code (H, B, E, G, I, T, S, P) or `null` for residues
classified as loop by mkdssp (space in the native DSSP output).
DSSP never outputs a literal 'C' character — the 'C' convention
is a downstream simplification used only in the 3-state `dssp3`
reduction (H/E/C). Note: the previously-published top8000 JSONL
(2026-03-15 vintage) did emit literal `'C'` for loop residues; this
was a `parse_dssp_output` artefact corrected in pydangle commit
`03e5511`.

## DSSP Regression Validation

DSSP values in the JSONL were validated against independently
generated `.dssp` files (from `pipeline/run_dssp.py` on cleaned
original PDB files) using `pipeline/validate_dssp.py`. The
comparator reuses `pydangle_biopython.dssp.parse_dssp_output` so
the validation parser stays in sync with the production label code.

| Metric | Value |
|--------|-------|
| JSONL records (post-filter) | 1,568,561 |
| `.dssp` ground-truth files | 7,737 |
| Agreement (JSONL ↔ .dssp) | 1,567,955 |
| mkdssp-only (B-filter excluded from JSONL) | 2,146,658 |
| jsonl-only (mkdssp drops; pydangle keeps) | 60 |
| **Non-trivial mismatches** | **0** |
| Missing-dssp pdbids | 0 |
| Non-model-1 (alternate models, skipped) | 546 |

The 60 jsonl-only records are chain-terminal residues that mkdssp
drops from `.dssp` output but pydangle keeps in the JSONL. All have
`dssp = null` in the JSONL — no value disagreement, just a residue-
set boundary difference between pydangle and mkdssp's chain-end
behavior. Same pattern as the 2 jsonl-only in top100 and 8 in top500.

Canonical validation report: `top8000/top8000_dssp_validation.md`.
Regenerable via:

```
python pipeline/validate_dssp.py top8000/top8000_measures.jsonl \
    ~/Desktop/Data/top8000_pdbs_hom70/ \
    -o top8000/top8000_dssp_validation.md
```

## DSSP Fix History (2026-04-26 / 2026-04-27)

The 0-mismatch result above is the conclusion of a multi-stage
forensic investigation that ran across two days. A condensed history
is recorded here because the previously-published top8000 dataset
contains substantial DSSP errors that the v2.0.0 release will correct.

### Pre-fix state (published top8000 dataset, Zenodo)

The originally-published top8000 JSONL (`top8000_measures.jsonl.gz`,
generated 2026-03-15 with pydangle 0.5.1 pre-fix) shipped with two
overlapping classes of DSSP error:

1. **~270k loop-residue mislabels: `parse_dssp_output` `'C'` literal
   artefact.** Pre-fix `parse_dssp_output` mapped mkdssp space-character
   output to literal `'C'` in the JSONL. The corrected JSONL emits
   `null` for loop residues (`null = 258,216` post-fix vs ~7k pre-fix
   non-null-but-real-null ambiguity). The mkdssp output itself is
   identical; this is purely a label-encoding correction. Fixed by
   pydangle commit `03e5511`.

2. **1,771 records carrying `dssp = null` from silent pydangle
   failures.** For 6 structures with surviving B-filtered residues
   (1h64: 24 records, 1ryp: 406, 1wte: 526, 2gn0: 566, 2zzs: 80,
   3mt6: 169 — total 1,771; plus 3cmy: 0 leak because fully
   B-filtered), every published JSONL record carries `dssp = null`
   because pydangle's per-file `run_dssp()` returned no SS
   assignments. The cause was actually two distinct bugs working in
   sequence:

   **a. HETATM-passthrough silent failure** (resolved by pydangle
   `5eaed6b` — selective HETATM stripping). Affected 1wte and 2gn0:
   ligand HETATM in those structures collided with SEQRES protein
   alignment and caused mkdssp to silently produce zero residues.

   **b. mkdssp returncode-1 with valid stdout** (resolved by
   pydangle `c45effa`). Affected 1h64, 1ryp, 2zzs, 3mt6: these
   multi-chain structures (1h64 has 28 chains, etc.) overflow legacy
   DSSP HEADER/COMPND field widths. mkdssp emits the warning "This
   file contains data that won't fit in the original DSSP format"
   and exits with returncode 1, but the residue table in stdout is
   complete and correct. pydangle's wrapper previously discarded
   the stdout; now it trusts stdout when the residue-table marker
   is present. Verification: captured stdout produces byte-identical
   residue assignments to mkdssp's file-output mode (which exits 0
   for the same data).

### Auxiliary fix that mattered for top8000

3. **Stale ersatz files lacked SEQRES.** The previously-built top8000
   ersatz tree was generated with a pre-`8802e0f` `build_ersatz_only.py`
   that did not preserve SEQRES records. Without SEQRES, mkdssp on
   the ersatz lacks the canonical-sequence reference. Resolved
   2026-04-26 by rebuilding the ersatz tree with the current
   `build_ersatz_only.py`; old stale ersatz preserved at
   `/mnt/backups/nuc81/snapshots/prisant/ersatz_pre-seqres-rebuild/top8000/`.

### Per-structure correction

| Structure | Records | Published `null` | Corrected `null` | Records gaining real DSSP code |
|---|---:|---:|---:|---:|
| 1h64 | 24 | 24 (100%) | 4 (17%) | 20 |
| 1ryp | 406 | 406 (100%) | 50 (12%) | 356 |
| 1wte | 526 | 526 (100%) | 84 (16%) | 442 |
| 2gn0 | 566 | 566 (100%) | 90 (16%) | 476 |
| 2zzs | 80 | 80 (100%) | 13 (16%) | 67 |
| 3cmy | 0 | 0 | 0 | 0 (fully B-filtered) |
| 3mt6 | 169 | 169 (100%) | 18 (11%) | 151 |
| **Subtotal** | **1,771** | **1,771** | **259** | **1,512** |

Every one of these 1,771 records flips from "silently nulled" to
either a real DSSP code (1,512 records) or a *legitimately* null
loop residue (259 records).

### Disclosure plan

The published top8000 dataset on Zenodo continues to ship with the
1,771 silently-nulled records and ~270k `'C'`-literal records
described above. Per project decision (revised 2026-04-27, supersedes
the 2026-04-26 "defer until v2.0.0" decision), the corrections will
ship as a **v1.x bugfix release** (proposed top8000 patch version:
next patch from current, under the "Phase 1 DSSP corrections" release
umbrella covering all five datasets simultaneously). The v1.x bugfix
release ships before v2.0.0; v2.0.0 follows separately for format
rectification (augmented ersatz, footer ordering, schema unification
— qualitatively different work).

The shift from the original "defer all to v2.0.0" plan addresses three
concerns: (a) holding back known-correctable data while colleagues are
informed of the bugs creates a credibility gap counter to the project's
external data-integrity commitments; (b) Phase 1 (bugfix) and Phase 2
(format change) are qualitatively different changes, and decoupling
them gives each release a clean narrative; (c) the v1.x bugfix
release shrinks the "tainted data in the wild" window from months to
days/weeks.

### References

- Pydangle commits: `03e5511` (ANISOU strip + parse_dssp_output fix),
  `5e79f71` (SEQRES blank-chain fill), `5eaed6b` (selective HETATM
  strip), `c45effa` (mkdssp returncode-1 with valid stdout).
- Richardson commits: `8802e0f` (DSSP ground-truth pipeline +
  build_ersatz_only.py SEQRES preservation), `a6f0f0d` (run_dssp.py
  HETATM strip).
