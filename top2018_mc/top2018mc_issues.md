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
