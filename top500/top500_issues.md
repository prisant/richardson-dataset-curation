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
