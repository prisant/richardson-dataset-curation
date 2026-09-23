# top2018 MC-filtered — Known Issues

## RCSB Remediation Mismatches

4 chains out of 13,677 (0.029%) have residue numbering or chain ID
mismatches between the Richardson Lab pruned files and the current RCSB
PDB coordinates. These affect 47 mask residues out of 2,963,303 total
(0.0016%).

The mismatched residues are silently excluded by the post-filtering step
(mask residues not found in the ersatz PDB are skipped). No remediation
was applied — the effort/benefit ratio is unfavorable for 47 residues.

### 4oe9_B — 3 missing residues

Residues B:66, B:101, B:116 are in the mask but absent from the current
RCSB PDB. The current PDB has gaps at these positions (65→67, 100→102,
115→117). Likely RCSB remediation removed these residues due to model
quality or renumbered surrounding residues.

### 4z54_A — 22 missing residues

Mask contains residues A:23–A:52, but the current RCSB PDB chain A starts
at residue 53. The entire N-terminal segment was renumbered by approximately
+30, making 22 mask residues unmatchable.

### 5e7x_A — 17 missing residues

Mask contains residues A:10–A:26, but the current RCSB PDB chain A starts
at residue 27. The N-terminal segment was renumbered by approximately +17.

### 5kve_L — 5 missing residues

Mask contains 5 entries for L:30 with different insertion codes. The
current RCSB PDB has residues at L:28–L:33 without insertion codes.
Insertion code convention changed during remediation.

## Reduce Timeout

Structure 4zw9 (membrane transporter, 3,638 ATOM lines) timed out at the
default 120s limit. Successfully processed with an extended 300s timeout.
The FH file was generated manually and the ersatz/mask were built
successfully. No data lost.

## Non-standard Residues

Only 2 MSE (selenomethionine) residues in the entire dataset, mapped to
MET for analysis. No other non-standard residues survived quality filtering.

## DSSP Coverage

20,546 residues (0.69%) have null DSSP assignments. These are typically
residues in structures where mkdssp encounters issues (non-standard
connectivity, missing backbone atoms in neighboring chains, etc.).

## DSSP — mkdssp 4.x legacy-format chain-count cap (>26 chains)

**Problem.** The classic 80-column DSSP output format reserves a single
character (PDB column 22) for the chain identifier.  mkdssp 4.x silently
caps its legacy DSSP output at the first 26 distinct case-distinguished
chain IDs encountered (uppercase processed before lowercase) and emits
`This file contains data that won't fit in the original DSSP format` to
stderr with returncode 1 — but the partial .dssp file is written to disk
**with the surplus chains entirely absent.**  This is an mkdssp 4.x
format limitation, not a chain-naming issue.

**Affected structures (9 in top2018_mc).** Identified by direct test —
every structure with > 26 case-distinguished chain IDs is affected, and
all `mkdssp` errors in this dataset are from this cap (verified by
biconditional: candidates = errors = 9).

| pdbid | chains | curated chain(s) in mc |
|---|---:|---|
| 5Le5 |  28 | 6 chains |
| 2zzs |  32 | 1 chain |
| 5b5e |  38 | 1 chain |
| 5b66 |  38 | F |
| 5v2c |  38 | C, b, d |
| 4ub6 |  39 | 1 chain |
| 4ub8 |  39 | 2 chains |
| 3uoi |  48 | 1 chain |
| 4v4m |  60 | w |

These are large multi-monomer assemblies (photosynthesis reaction
centres, GroEL-class chaperonins, etc.) deliberately retained in
top2018_mc to sample rare main-chain conformations.

**Resolution — chain-batch fallback.** The pipeline now detects the
"won't fit" stderr and automatically falls back to
`pipeline/dssp_chain_batch_fix.py`:

1. Split the structure's chains into ordered batches of ≤ 26.
2. Run mkdssp on each batch independently (each fits the 26-chain cap;
   each returns RC=0 with no warnings).
3. Concatenate the per-batch residue tables into a single `.dssp` file
   with renumbered global sequential indices and BP1/BP2 references
   shifted accordingly.

After fallback, **all input chains and all residues are present** in the
ground-truth `.dssp` file (no chain dropping).

**Information loss disclosure — cross-batch H-bonds.** Per-batch mkdssp
runs cannot detect hydrogen bonds whose donor and acceptor are in
different batches.  Empirical bound, on a worst-case test structure
(20-chain inter-chain β-sheet, 13/7 split):

- Per-residue **secondary-structure code (H/E/B/G/I/T/S/P)** preserved
  with **>99.5% fidelity** vs full-run reference.
- Per-residue **ACC, TCO, KAPPA, ALPHA, PHI, PSI, X-CA/Y-CA/Z-CA** are
  computed locally per residue and are **fully preserved.**
- Per-residue **4 H-bond partner offset+energy columns** lose the small
  fraction of partners that span the batch boundary (~4% of residues
  affected on the worst-case test; less on homo-multimers where
  inter-batch contact is small).

**For the 9 affected structures specifically:** the curated chain in
each (the chain that contributes residues to top2018_mc) was empirically
verified to be **inside the kept-26 set even before the fallback** —
i.e., the JSONL DSSP assignments consumed by downstream analysis were
already correct for these structures.  The fallback corrects the
ground-truth `.dssp` file used by `validate_dssp.py`, not the JSONL
science output.

**Validation gate.**  `validate_dssp.py` 0-mismatch result is preserved
after the fallback; counts are unchanged for residues already present.
The fallback eliminates the "missing-dssp" failure mode that would
otherwise be triggered for surplus chains in any future dataset.
