# top2018 Full-Filtered Dataset: Preparation Issues and Resolutions

This document records every issue encountered during preparation of the
top2018 full-filtered pydangle dataset, how each was resolved, and what
impact (if any) it has on the final published data.

## 1. DSSP on single-chain files gives incorrect assignments

**Problem:** The Richardson Lab pruned files are single-chain extracts.
Running mkdssp on these isolated chains loses inter-chain hydrogen bond
context, producing incorrect secondary structure assignments for
residues near chain interfaces.

**Resolution:** Created "ersatz" full-structure PDB files that preserve
all chains from the original crystal structure. mkdssp runs on the
complete multi-chain ersatz file, giving correct H-bond context. Output
is then post-filtered to the quality-filtered residues.

**Impact:** 34% of DSSP assignments differ from the single-chain
approach. This is the intended improvement.

## 2. NQH flip method: renameflip vs hinge rotation

**Problem:** Reduce's `-renameflip` flag swaps atom names without
changing coordinates, while the default behavior rotates coordinates
around the CB-CG hinge axis. The Zenodo pruned files use the hinge
rotation method, producing coordinates that differ by up to 0.16 A
from the renameflip approach.

**Discovery:** Comparing HIS B 201 in 1a4i between our initial
`-renameflip` output and the Zenodo pruned file showed exact coordinate
identity for renameflip (name swap only) but different coordinates in
the Zenodo file (hinge rotation with CG displacement).

**Resolution:** Re-ran Reduce without `-renameflip` to match the Zenodo
approach. Verified exact coordinate match on test cases.

**Impact:** chi1 values for ~993 NQH-flipped residues differ slightly
from the renameflip approach. The hinge rotation values match the
published Zenodo coordinates.

**Reference:** Williams, C. J., Headd, J. J., Moriarty, N. W., et al.
(2018). *Protein Science*, 27(1), 293-315. doi:10.1002/pro.3330

## 3. Reduce version differences: 35 borderline flip decisions

**Problem:** Our local Reduce (3.3.160422) makes different flip/keep
decisions than the Zenodo version (3.3.160602) for 35 residues across
22 structures. These are all borderline scoring cases where the flip
vs keep energy difference is small.

**Resolution:** Accepted as inherent version variation. The 35
disagreements out of 464,212 decisions (0.008%) are within the noise
of the scoring function.

**Impact:** For the 35 affected residues, sidechain atom positions may
differ slightly from the Zenodo pruned files. All are borderline cases
by definition.

## 4. Three PDB files no longer available from RCSB in PDB format

**Problem:** RCSB no longer serves PDB-format files for 4v4m, 4yza,
and 4ztt (only mmCIF is available). Michael's original download script
successfully obtained these files in an earlier run, but they cannot
be re-downloaded from RCSB today.

**Resolution:** Located copies in local archives:
- 4v4m.pdb from `FromMolpy01/FilteringPdbs/MiscStrucFiles/`
- 4yza.pdb and 4ztt.pdb from `ClevelandFilteredResidues/top2018_pdbs_mc_filtered_hom70/`

These are included in the `supplementary_pdbs/` directory of the
GitHub repository to ensure full reproducibility.

**Impact:** None on the final data. All three structures are fully
processed.

## 5. Two Reduce segfaults: 4oan and 4uzg

**Problem:** Reduce 3.3.160422 crashes with a segmentation fault on
4oan.pdb and 4uzg.pdb, regardless of flags (`-allalt`, `-norotmet`,
different het dictionaries). The crash occurs in the dot-scoring
engine during sidechain optimization.

**Resolution:** Used Reduce 4.16.250520 (`reduce-master` from the
Richardson Lab GitHub) with its companion het dictionary. Both files
processed successfully.

**Impact:** None on the final data. Both structures are fully processed.
The two versions of Reduce produce equivalent results for all other
structures (validated by the NQH flip comparison).

## 6. One Reduce timeout: 4zw9

**Problem:** Reduce on 4zw9.pdb exceeded the 120-second timeout in the
parallel processing pool.

**Resolution:** Re-ran with a 300-second timeout. Completed successfully
in ~180 seconds.

**Impact:** None.

## 7. Multi-chain entries: 282 chains initially missed

**Problem:** 255 PDB structures in the dataset contribute multiple
chains (e.g., 5t5i has chains A, C, G, J, L). The initial
`build_ersatz_pdbs.py` script processed only the first pruned file per
directory (`pruned_files[0]`), missing 282 chains. This was not caught
before the initial Zenodo publication because the title said "12,125
chains" (from the chain list) while the data contained only 11,843
(unique PDB IDs).

**Root cause:** Confusion between "structure" count (11,843 unique PDB
IDs) and "chain" count (12,125 entries in the chain list). The
discrepancy was visible from the start but was not flagged or
investigated.

**Resolution:** Fixed `build_ersatz_pdbs.py` to loop over ALL pruned
files per entry, generating a mask for each chain. Fixed
`filter_pruned_residues.py` to merge masks for the same PDB ID
instead of overwriting. Verified: 12,125 mask files, 12,123 chains
in final output (2 chains excluded due to issue #8).

**Impact:** Zenodo versions 1.0.0 through 1.0.2 contain data for only
11,843 chains. Version 1.0.3 contains the corrected 12,123-chain
dataset.

## 8. RCSB vs Zenodo chain ID and residue numbering mismatches

**Problem:** For 6 entries, the full PDB file downloaded from RCSB has
different chain identifiers or residue numbering than the Zenodo pruned
files. This prevents the post-filter from matching mask residues to
pydangle output.

| Entry | Chain | Missing | Issue | Resolution |
|-------|-------|---------|-------|------------|
| 4ct3 | D | 137/137 (100%) | Chain ID renamed: RCSB has E,G,I,K; Zenodo has D | Remapped mask D→E |
| 5tvo | A | 194/194 (100%) | Renumbered: RCSB 402-671, Zenodo 87-353 | Remapped mask +315 offset |
| 1xmz | B | 112/159 (70%) | RCSB split chain B into B (≤62) and C (≥66) | Remapped mask residues ≥66 from B→C |
| 5kve | L | 4/88 (5%) | 4 insertion-code residues (30B,30D,30E,30F) removed by RCSB remediation | Unfixable — residues absent |
| 3wwL | A | 1/48 (2%) | Residue 54 (GLU) removed; RCSB chain ends at 53 | Unfixable — residue absent |
| 4oe9 | B | 2/99 (2%) | Residues 66 (GLN) and 101 (LEU) removed by RCSB remediation | Unfixable — residues absent |

**Root cause:** The RCSB PDB archive is periodically remediated.
Structures may be renumbered, re-chained, or have coordinates revised
between the time the Richardson Lab prepared the Zenodo dataset and our
download. The Zenodo pruned files were prepared from an earlier PDB
archive snapshot.

**Resolution:** Three entries (4ct3, 5tvo, 1xmz) were fixed by
remapping the mask and USER INC fragment records to match the RCSB
chain IDs and residue numbering. The remaining three entries (5kve,
3wwL, 4oe9) have 7 residues that were genuinely removed from the
current RCSB structures and cannot be recovered.

**Impact:** 7 residues out of 2,405,329 mask residues unfixable
(0.0003%). Final dataset contains all 12,125 chains from all 11,843
structures.

## 9. Fragment-aware filtering: quality boundary measurements

**Problem:** The ersatz files contain all residues from the original
structure, including those that failed quality filtering. Geometric
measurements that depend on neighbor atoms (phi, psi, omega) could
silently use coordinates from filtered-out residues if the measurement
reaches across a quality-filter boundary.

**Discovery:** The old pruned-file approach had the same issue in
reverse: BioPython's CaPPBuilder uses CA-CA distance (< 4.6 A) to
define polypeptide fragments, not quality filter status. Filtered-out
residues with nearby CA atoms were silently included in fragments,
allowing tau and other measurements to use their coordinates.

**Resolution:** Post-filter uses `USER INC` fragment records from the
Richardson Lab pruned files as the authoritative source for fragment
boundaries. Measurements at fragment boundaries are nulled:
- phi, omega, is_cis, is_trans: null when i-1 was quality-filtered
- psi: null when i+1 was quality-filtered
- rama_category: null when either neighbor was quality-filtered
- tau, chi1, dssp: unaffected (use only current residue's atoms, or
  full-structure context)

**Impact:** ~370K residues at fragment starts have phi/omega nulled;
~370K at fragment ends have psi nulled. rama_category is null for ~662K
residues (where either phi or psi is unavailable). This is more
conservative than the old approach but more correct.

## 10. Stale intermediate deleted during cleanup

**Problem:** The 1.6 GB raw pydangle output
(`top2018_ersatz_raw.jsonl`) was deleted during a working directory
cleanup, requiring a full re-run (~2 hours) when the multi-chain mask
bug was discovered afterward.

**Resolution:** Re-ran pydangle. Output was identical (6,150,347 lines)
since the ersatz files hadn't changed.

**Lesson:** Expensive intermediate outputs should be preserved until
the final dataset is confirmed complete and all validation checks pass.

## Summary

| Issue | Entries affected | Residues affected | Resolution |
|-------|-----------------|-------------------|------------|
| Multi-chain masking | 255 entries, 282 chains | ~49K residues recovered | Fixed in scripts |
| RCSB/Zenodo mismatch | 6 entries | 3 fixed by remapping, 7 residues unfixable (0.0003%) | Mask remapping + documented |
| Reduce segfaults | 2 entries | 0 | Used Reduce 4.16 |
| Reduce timeout | 1 entry | 0 | Extended timeout |
| Missing RCSB PDB files | 3 entries | 0 | Local archive copies |
| NQH flip method | All entries | ~993 chi1 values | Matched Zenodo method |
| Flip version diffs | 22 entries, 35 residues | Negligible | Accepted |
| Fragment boundary | ~378K residues | Measurements nulled | USER INC filtering |

## Notes for Future Work

Several findings from this preparation effort have implications for
processing other Richardson Lab datasets (top2018_mc, top8000, top500,
top100) and for the longer-term goal of building a custom filtration
pipeline.

**RCSB archive remediation is an ongoing problem.** The PDB archive is
periodically updated with chain renames, residue renumbering, and
coordinate revisions. Any pipeline that combines current RCSB downloads
with a historical dataset (like the Richardson Lab Zenodo deposits) must
detect and handle these discrepancies. In this dataset, 6 of 11,843
entries were affected. The rate may be higher for older structures. A
robust approach would be to validate every entry by comparing the
sequence and numbering between the downloaded PDB and the pruned file,
and apply remapping automatically where the alignment is unambiguous.

**Multi-chain entries must be handled explicitly.** The Richardson Lab
chain lists contain more chain entries than unique PDB IDs because some
structures contribute multiple quality-filtered chains. Any script that
iterates over entry directories (one per PDB ID) must process all pruned
files and mask files in each directory, not just the first one found.
The chain list file is the authoritative count.

**Reduce version matters.** Reduce 3.3 segfaults on 2 structures that
Reduce 4.16 handles correctly. For future datasets, prefer Reduce 4.16
(`~/Desktop/Downloads/reduce-master/reduce_src/reduce`) with its
companion het dictionary. The NQH flip decisions are nearly identical
between versions (99.77% agreement on 464K decisions).

**Fragment-aware filtering is essential.** Using CA-CA distance alone
(as BioPython's CaPPBuilder does) to define polypeptide fragments allows
measurements to silently use coordinates from quality-filtered residues.
The `USER INC` fragment records in the Richardson Lab pruned files are
the authoritative source for fragment boundaries and should always be
used for neighbor-dependent measurement validation.

**Three PDB files are no longer available from RCSB.** The structures
4v4m, 4yza, and 4ztt are only available in mmCIF format from RCSB.
PDB-format copies are preserved in `supplementary_pdbs/`. This problem
will likely grow as RCSB continues migrating large structures away from
PDB format. Future pipelines may need to support mmCIF input or maintain
a local PDB-format archive.

**Count validation should be built into the pipeline.** Every stage
should validate its output count against the chain list: number of
mask files must equal number of chain entries, number of chains in
the final output must equal the same. Discrepancies should halt the
pipeline, not be silently absorbed.
