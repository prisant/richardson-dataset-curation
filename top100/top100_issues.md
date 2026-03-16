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
