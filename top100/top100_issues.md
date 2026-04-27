# top100 — Known Issues

## B-factor Threshold: 40 (Different from top500/top8000)

This dataset uses a mainchain B-factor ≤ 40 filter, in contrast to the
B ≤ 30 used for the top500 and top8000 datasets. This matches the
original published methodology for the Top100 era (Word et al., 1999a).

The original Top100 used B ≤ 40 because with only 100 structures, a
stricter threshold would have sacrificed too many residues for adequate
statistics. The later, larger datasets (Top500 with 500 structures,
Top8000 with ~8000 chains) had sufficient data to afford the stricter
B ≤ 30 threshold introduced by Williams et al. (2018).

This is an intentional methodological difference that respects the
published practice for each dataset's era, rather than imposing a
uniform modern threshold retroactively.

## 2 Chains with Zero Passing Residues (1.9%)

Of the 108 chain entries derived from the Top100 chain list, 2 chains
(1.9%) have zero residues passing even the B ≤ 40 filter. Both chains
exist in the current RCSB PDB files with correct chain IDs. This
exclusion rate is consistent with the 1.7–1.8% observed in the
top500 and top8000 datasets (which use the stricter B ≤ 30).

2 structures are entirely excluded, reducing the output from 100 to
98 structures.

## Chain List Provenance

The Top100 chain list was reconstructed from the kinemage webserver
archive (top100H directory, 100 FH files). The Top100 is the original
Richardson Lab quality-filtered reference dataset, used for the
foundational Asn/Gln sidechain amide flip correction work (Word et al.,
1999a, doi:10.1006/jmbi.1998.2401).

The FH filename format encodes PDB ID and chain(s): `{pdbid}{chains}H`.
Three entries use a `bio` prefix (bio1rpo, bio2wrp) indicating
biological assembly files. 80 entries have no explicit chain ID
(single-chain structures), 17 have one chain, and 1 has two chains
(1benAB). Expansion and PDB inspection produces 108 chain entries.

## RCSB Remediation

No chain ID or residue numbering mismatches were detected.

## Reduce Processing

All 100 structures processed successfully with Reduce 4.16.250520.
Zero errors, zero timeouts.

## Filtering Methodology

Same pipeline as top8000 and top500: B-factor filter (here B ≤ 40)
with synthetic USER INC fragment records inferred from contiguous
passing residues. No Richardson Lab pruned PDB files are available
for this dataset.

## DSSP File Generation

DSSP output files (`.dssp`) were generated for all 100 structures
using `pipeline/run_dssp.py` with mkdssp 4.2.2. The script builds a
minimal cleaned PDB (HEADER + CRYST1 + SEQRES + ATOM + TER +
protein-only HETATM) before invoking mkdssp to work around four
mkdssp 4.x issues:

1. **ANISOU records** cause mkdssp to produce incomplete residue
   assignments (4 structures affected: 1iro, 1lkk, 1rge, 2erl).
2. **SIGATM/SIGUIJ records** cause mkdssp to produce zero residues
   (1 structure: 1etm).
3. **Blank chain IDs** (column 22 of ATOM/TER, column 12 of SEQRES)
   in old depositions cause parse failures or silent zero-residue
   output (3 structures: 2wrp, 3b5c, 4ptp). Blanks are filled
   with the first non-blank ATOM chain ID.
4. **Non-protein HETATM records** (water, sulfate, glycerol, ions,
   sugars, ...) collide with the SEQRES protein alignment and cause
   mkdssp to produce zero residues for the entire structure
   (0 structures affected in top100 — verified by backward-compat
   regression run; 4 structures affected in top500: 1a1y, 1gdo, 1tax,
   2myr). HETATM lines whose residue name is in the file's SEQRES or
   ATOM record set (i.e., modified amino-acid residues like MSE, PCA,
   MEN, PTR, LOV that some old PDBs deposit as HETATM) are kept;
   other HETATM are stripped.

SEQRES records are included so that mkdssp correctly identifies
modified amino acids (MEN, PTR, PCA, LOV) as part of the polypeptide
chain, which affects DSSP assignments for neighboring residues.

All 100 structures produce `.dssp` files with zero errors.

## DSSP Semantics

The `dssp` column in the JSONL output contains the mkdssp secondary
structure code (H, B, E, G, I, T, S, P) or `null` for residues
classified as loop by mkdssp (space in the native DSSP output).
DSSP never outputs a literal 'C' character — the 'C' convention
is a downstream simplification used only in the 3-state `dssp3`
reduction (H/E/C).

## Regression Validation: DSSP

DSSP values in the JSONL were validated against independently
generated `.dssp` files (from `pipeline/run_dssp.py` on original
PDB files) using `pipeline/validate_dssp.py`. The comparator reuses
`pydangle_biopython.dssp.parse_dssp_output` to ensure the validation
parser stays in sync with the production label code.

| Metric | Value |
|--------|-------|
| JSONL records (post-filter) | 17,434 |
| `.dssp` ground-truth residues | 21,010 |
| Agreement (JSONL ↔ .dssp) | 17,432 |
| mkdssp-only (B-filter excluded from JSONL) | 3,578 |
| jsonl-only (mkdssp drops; pydangle keeps) | 2 |
| Non-trivial mismatches | 0 |
| Missing-dssp pdbids | 0 |

The 3,578 mkdssp-only residues are residues present in the full
structure but excluded by B-factor filtering — expected and correct.

The 2 jsonl-only records are 1cka chain A residue 190 and chain B
residue 9 (both with `dssp = null`). These are chain-terminal residues
that mkdssp drops from `.dssp` output but pydangle keeps in the
JSONL. There is no value disagreement (both are null), only a
residue-set boundary difference between the two tools' chain-end
behavior. Documented for the record; not a Phase 1 bug.

The canonical validation report is at
`top100/top100_dssp_validation.md`, regenerable on demand by
`python pipeline/validate_dssp.py top100/top100_measures.jsonl.gz
~/Desktop/Data/top100_pdbs/ -o top100/top100_dssp_validation.md`.

## Regression Validation: Backbone Geometry

Phi, psi, omega, and tau values were validated against Java Dangle
1.07 (`/home/prisant/bin/dangle`) run on the same ersatz PDB files.
Full details in pydangle-biopython `BENCHMARK.md`.

### Summary (top100pdb, 21,003 matched rows)

| Metric | Value |
|--------|-------|
| Rows agreeing (≤ 0.01°) | 20,976 (99.87%) |
| Coverage diffs | 14 |
| Value mismatches | 13 |
| Rows only in pydangle | 4 |

### Value mismatches (13)

All 13 trace to two root causes:

- **Chain break boundary differences (4 rows):** 1cpc residues
  71–75 (adjacent to missing residues 72–74) and 1fxd residue 12
  (adjacent to missing residue 11). Both tools agree on the break
  but select different atoms at fragment boundaries.

- **Alternate conformation selection (9 rows):** 1ctj residues 1–3,
  1ifc residues 30/31/47/48, 1fxd residue 10. BioPython and Java
  Dangle select different altloc atoms when multiple conformations
  have equal or near-equal occupancy.

Neither represents a computational error.
