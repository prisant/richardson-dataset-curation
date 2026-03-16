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
