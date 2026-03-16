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
