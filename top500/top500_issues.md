# top500 — Known Issues

## 9 Chains with Zero Passing Residues (1.7%)

Of the 523 chain entries derived from the Top500 chain list, 9 chains
(1.7%) have zero residues passing the mainchain B-factor ≤ 30 filter.
All 9 chains exist in the current RCSB PDB files with correct chain IDs
(no remediation issues). Every residue in these chains has at least one
mainchain atom (N, CA, C, or O) with B-factor exceeding 30.

This exclusion rate is consistent with the 1.8% (143/7,957) observed
in the top8000 dataset using the same filter. See the top8000 issues
documentation for detailed characterization of this phenomenon.

2 structures are entirely excluded (all chains fail), reducing the
output from 494 to 492 structures.

## Chain List Provenance

The Top500 chain list was reconstructed from the kinemage webserver
archive. The original list format is tab-separated (PDB ID + chain
string), with 327 entries having no explicit chain ID (single-chain
structures) and 167 entries with one or more chain IDs. Multi-chain
strings were expanded to individual chain entries, and single-chain
structures were resolved by inspecting the PDB files, producing 523
total chain entries from 494 structures.

There is no formal published chain list for the Top500 comparable to
the Top8000's CSV files on GitHub. The kinemage webserver FH file
collection (500 files) serves as the authoritative source.

## Regression Validation Against Richardson Lab Reference Values

Phi and psi values were validated against the Richardson Lab's
pre-computed angle data from the top500-angles archive (2006 vintage,
from kinemage.biochem.duke.edu). The reference data uses a dipeptide
format where each line reports the phi/psi pair at a peptide bond
between residues i and i+1: phi belongs to residue i+1, psi to
residue i.

### Match statistics

| Metric | Value |
|--------|-------|
| Reference phi entries | 37,433 |
| Reference psi entries | 37,433 |
| Matched phi pairs | 36,810 |
| Matched psi pairs | 36,809 |

### Phi comparison (36,810 pairs)

| Threshold | Count | Percent |
|-----------|------:|--------:|
| \|diff\| < 0.01° | 36,443 | 99.0% |
| \|diff\| < 1.0° | 36,707 | 99.7% |
| \|diff\| > 1.0° | 103 | 0.3% |

Mean |diff|: 0.130°. Median |diff|: 0.003°.

### Psi comparison (36,809 pairs)

| Threshold | Count | Percent |
|-----------|------:|--------:|
| \|diff\| < 0.01° | 36,453 | 99.0% |
| \|diff\| < 1.0° | 36,701 | 99.7% |
| \|diff\| > 1.0° | 108 | 0.3% |

Mean |diff|: 0.241°. Median |diff|: 0.003°.

### Interpretation

99.7% of both phi and psi values agree within 1 degree. The ~100
large discrepancies per angle are concentrated in a few structures
(particularly 1tph, trypsinogen) and show ~180° differences
characteristic of alternate conformation selection. The reference
data was computed from circa-2005 PDB coordinates; our data uses
current RCSB coordinates processed with Reduce 4.16.

The slightly lower exact-match rate compared to the top8000 regression
(99.0% vs 99.9% within 0.01°) likely reflects RCSB coordinate
remediation over the intervening 20 years between the two dataset
vintages (top500 from ~2000 vs top8000 from 2011).

## RCSB Remediation

No chain ID or residue numbering mismatches were detected. All chain
IDs from the reconstructed chain list exist in the current RCSB PDB
files. This is notable for a dataset of structures deposited primarily
in the late 1990s.

## Reduce Processing

All 494 structures processed successfully with Reduce 4.16.250520.
Zero errors, zero timeouts.

## Non-standard Residues

One non-standard residue (CSD, 3-sulfinoalanine) in the filtered
output. The same residue appears in the top8000 dataset.

## Filtering Methodology

Same as top8000: mainchain B-factor ≤ 30 filter with synthetic
USER INC fragment records inferred from contiguous passing residues.
No Richardson Lab pruned PDB files are available for this dataset.

## DSSP File Generation

DSSP output files (`.dssp`) were generated for all 494 structures
using `pipeline/run_dssp.py` with mkdssp 4.2.2. The script builds a
minimal cleaned PDB (HEADER + CRYST1 + SEQRES + ATOM + TER +
protein-only HETATM) before invoking mkdssp to work around four
mkdssp 4.x issues:

1. **ANISOU records** cause mkdssp to produce incomplete residue
   assignments. Affected top500 structures: 1iro, 1lkk, 1rge, 2erl
   (carried over from the top100 fix).
2. **SIGATM/SIGUIJ records** cause mkdssp to produce zero residues.
   No top500 structures rely solely on this fix path.
3. **Blank chain IDs** (column 22 of ATOM/TER, column 12 of SEQRES)
   in old depositions cause parse failures or silent zero-residue
   output. Top500 has 0 structures with blank-chain SEQRES (verified
   by survey 2026-04-26).
4. **Non-protein HETATM records** (water, sulfate, glycerol, ions,
   sugars, ...) collide with the SEQRES protein alignment and cause
   mkdssp to produce zero residues for the entire structure.
   **Affected top500 structures: 1a1y, 1gdo, 1tax, 2myr** (4/494 by
   path-A hard error; same 4 by path-B silent failure under the
   pre-fix pydangle code). HETATM lines whose residue name is in the
   file's SEQRES or ATOM record set (i.e., modified amino-acid
   residues like MSE, PCA, MEN, PTR, LOV that some old PDBs deposit
   as HETATM) are kept; other HETATM are stripped.

SEQRES records are included so that mkdssp correctly identifies
modified amino acids (MEN, PTR, PCA, LOV, MSE, etc.) as part of the
polypeptide chain, which affects DSSP assignments for neighboring
residues. The `top500_pdbs/*/<pdbid>/_ersatz.pdb` files were rebuilt
2026-04-26 with `build_ersatz_only.py` post-commit `8802e0f` to
include SEQRES; previously they did not.

All 494 structures produce non-empty `.dssp` files with zero errors
under the post-fix pipeline.

## DSSP Semantics

The `dssp` column in the JSONL output contains the mkdssp secondary
structure code (H, B, E, G, I, T, S, P) or `null` for residues
classified as loop by mkdssp (space in the native DSSP output).
DSSP never outputs a literal 'C' character — the 'C' convention
is a downstream simplification used only in the 3-state `dssp3`
reduction (H/E/C). Note: the previously-published top500 JSONL
(2026-03-15 vintage) did emit literal `'C'` for loop residues; this
was a `parse_dssp_output` artefact corrected in pydangle commit
`03e5511`.

## DSSP Regression Validation

DSSP values in the JSONL were validated against independently
generated `.dssp` files (from `pipeline/run_dssp.py` on original
PDB files) using `pipeline/validate_dssp.py`. The comparator reuses
`pydangle_biopython.dssp.parse_dssp_output` so the validation
parser stays in sync with the production label code.

| Metric | Value |
|--------|-------|
| JSONL records (post-filter) | 97,519 |
| `.dssp` ground-truth residues | 168,114 (across 494 structures) |
| Agreement (JSONL ↔ .dssp) | 97,511 |
| mkdssp-only (B-filter excluded from JSONL) | 70,603 |
| jsonl-only (mkdssp drops; pydangle keeps) | 8 |
| **Non-trivial mismatches** | **0** |
| Missing-dssp pdbids | 0 |

The 70,603 mkdssp-only residues are residues present in the full
structure but excluded by the mainchain B≤30 filter — expected and
correct; they reflect the dataset's pruning criterion, not a DSSP
discrepancy.

The 8 jsonl-only records are chain-terminal residues that mkdssp
drops from `.dssp` output but pydangle keeps in the JSONL. All 8
have `dssp = null` in the JSONL, so there is no value disagreement —
only the residue-set boundary differs between pydangle and mkdssp's
chain-end behavior. Documented for the record; not a bug.

The canonical validation report is at `top500/top500_dssp_validation.md`,
regenerable via:

```
python pipeline/validate_dssp.py top500/top500_measures.jsonl \
    ~/Desktop/Data/top500_pdbs/ -o top500/top500_dssp_validation.md
```

## DSSP Fix History (2026-04-26)

The 0-mismatch result above is the conclusion of a multi-stage
forensic investigation that began with the top500 propagation pass.
A condensed history is recorded here because the previously-published
top500 JSONL contains a substantial number of pre-fix DSSP errors.

### Pre-fix state (published top500 dataset, Zenodo)

The originally-published top500 JSONL (`top500_measures.jsonl.gz`,
generated 2026-03-15 with pydangle 0.5.1 pre-fix) showed **18,608
non-trivial mismatches** when validated against the corrected
ground-truth `.dssp` files generated by `pipeline/run_dssp.py`.
The mismatches fell into two categories:

1. **~18,000+ mismatches: parse_dssp_output `'C'` literal artefact.**
   Pre-fix `parse_dssp_output` mapped mkdssp space-character output
   to literal `'C'` in the JSONL. `validate_dssp.py` correctly
   reports `'C'` ≠ any other DSSP code as a mismatch. Fixed by
   pydangle commit `03e5511` (parse_dssp_output now returns `None`
   for space; DSSP never emits `'C'`).

2. **569 silent-null records: HETATM-passthrough silent failure.**
   For 3 structures (1a1y: 62 records, 1gdo: 205, 1tax: 302; plus
   2myr fully B-filtered → 0 leak), every published JSONL record
   carries `dssp = null` because pydangle's per-file `run_dssp()`
   call returned empty mkdssp output. Root cause: the input ersatz
   files contained HETATM (water / sulfate / glycerol / ...) which
   collided with the SEQRES protein alignment and caused mkdssp 4.x
   to silently produce zero residues. Fixed by pydangle commit
   `5eaed6b` (selective HETATM stripping; protein-residue HETATM
   like MSE/PCA/MEN/PTR/LOV are kept, non-protein HETATM are
   removed).

### Auxiliary fix that mattered for top500

3. **Stale ersatz files lacked SEQRES.** The previously-built top500
   ersatz tree (`top500_pdbs/*/<pdbid>/_ersatz.pdb`) was generated
   with a pre-`8802e0f` `build_ersatz_only.py` that did not preserve
   SEQRES records. Without SEQRES, mkdssp on the ersatz lacks the
   canonical-sequence reference and labels modified amino acids
   inconsistently. Resolved 2026-04-26 by rebuilding the ersatz tree
   with the current `build_ersatz_only.py`, which preserves SEQRES.
   Top100 ersatz files were already built with the post-`8802e0f`
   script (have SEQRES); top500/top8000/top2018 ersatz trees were
   stale and required rebuild.

### Affected residue counts (published vs corrected)

| Structure | Records | Published `dssp = null` | Corrected `dssp = null` | Records gaining real DSSP code |
|---|---:|---:|---:|---:|
| 1a1y | 62 | 62 (100%) | 14 (23%) | 48 |
| 1gdo | 205 | 205 (100%) | 35 (17%) | 170 |
| 1tax | 302 | 302 (100%) | 44 (15%) | 258 |
| 2myr | 0 | 0 | 0 | 0 (fully B-filtered) |
| **Subtotal (HETATM bug)** | **569** | **569** | **93** | **476** |
| All other top500 records | 96,950 | 2,059 | 16,732 | n/a* |
| **Total** | **97,519** | **2,628** | **16,825** | **n/a** |

\* The "All other" row reflects the parse_dssp_output `'C'`-as-loop
artefact: in the published JSONL, ~14,000 records had literal `'C'`
(treated as a value) where the corrected JSONL has `null` (loop).
This is a label correction, not a value-discovery; the underlying
mkdssp output is identical.

### Disclosure plan

The published top500 dataset on Zenodo continues to ship with the
569 silently-nulled records and ~14,000 `'C'`-literal records
described above. Per project decision (revised 2026-04-27, supersedes
the 2026-04-26 "defer until v2.0.0" decision), the corrections will
ship as a **v1.x bugfix release** (proposed top500 patch version: v2.0.3,
under the "Phase 1 DSSP corrections" release umbrella covering all five
datasets simultaneously). The v1.x bugfix release ships before v2.0.0;
v2.0.0 follows separately for format rectification (augmented ersatz,
footer ordering, schema unification — qualitatively different work).

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
  strip).
- Richardson commits: `8802e0f` (DSSP ground-truth pipeline +
  build_ersatz_only.py SEQRES preservation), `a6f0f0d` (run_dssp.py
  HETATM strip).
- Diagnostic provenance for the original `'C'` discovery is preserved
  as `top100/top100_analysis_plots_files/*.png` in the corresponding
  top100 docs (committed 2026-04-26).
