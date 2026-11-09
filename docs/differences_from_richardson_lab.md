# Differences between these datasets and the Richardson Lab reference data

*File: differences_from_richardson_lab.md · Git: Michael G. Prisant e6c33d0 · Txd: 11/09/26 08:09 with local modifications!*

**Date:** 2026-11-09.
**Status:** DRAFT for review by M. G. Prisant. Statements about the Richardson Lab's procedures must be confirmed before this file is merged or cited.
**Scope:** how the pydangle measures databases in this repository (top100, top500, top8000, top2018_full, top2018_mc; v1.x pipeline) differ from the Richardson Lab's own reference data, as far as we understand that work today.

## How to read this file

Each statement about the lab's work carries its evidence:

| Tag | Meaning |
|---|---|
| **[M]** | Measured from the lab's own data or files by a committed script (script and output named) |
| **[P]** | Stated in a publication (page given) |
| **[L]** | Said by Richardson Lab members (online meeting 2026-10-14: J. S. Richardson, C. J. Williams, V. B. Chen; or the personal communication dated) |
| **[I]** | Inferred by us from counts; not confirmed by the lab |
| **[U]** | Unknown |

Statements about our own pipeline are from the code in this repository and in pydangle-biopython, and are not tagged.

Reference material used:

- **top8000:** the lab's per-residue Ramachandran lists `Top8000_Rama_2013.zip` (written 2013-02-16 by `get_srcdata.py` from the lab's top8000 MySQL database; Silk on these lists reproduces all six archived rama8000 `.ndft` tables bin for bin), and the lab's 2011 chain files `top8000_chains_70.tgz`, which are the input of those lists. Analysis on branch `rama8000-reconciliation` (`e00fa7b`), directory `rama8000_reconciliation/` (README sections "The lab's 2013 residue lists" and "The lab's 2011 chain files and 2011 tools").
- **top2018:** the lab's pruned files (Zenodo 10.5281/zenodo.5115233) and C. J. Williams' Top2018 SQLite database of validation and geometry measures (Zenodo 10.5281/zenodo.22085722).
- **top100, top500:** no per-residue reference list has been recovered; only the publications cited in each README.
- **Publications:** Williams et al. 2018, *Protein Sci* 27:293–315 (WHM+18); Read et al. 2011, *Structure* 19:1395–1412 (RAA+11).

## 1. Which residue types get records

| | |
|---|---|
| **These datasets (v1.x)** | The 20 standard amino acids and selenomethionine (MSE) only (`filter_pruned_residues.py --residues std20+MSE`, used by every `reproduce.sh`). MSE keeps the name MSE in the records. Other modified residues (MLY, CSD, PCA, CSO, KCX, LLP, ...) get no record. The v1.0 datasets kept whatever pydangle reports (its list, copied from Java Dangle, also has MLY, M3L, CSD, CYO, PCA, AIB, ABU, DGN, ASX, GLX), but the masks dropped every residue written as `HETATM` (register B1), so in practice v1.0 had almost no modified residues. |
| **top8000 lab lists** | Standard residues and MSE. [M] 3,440 non-standard residues pass the best-fitting rule; 2,985 are MSE; adding the other 455 (mostly MLY) worsens the fit to the lab's grids (`best_rule_summary.py`, `nonstandard_types.py`, c406ccd). [L] "Best guess: only MSE included; not certain." |
| **Top2018 database** | Standard residues and MSE. [M] MSE: 11,116 of 11,532 residues have a Ramachandran row. All 110 other non-standard types together: 13 of 1,045 have one, and none passes the Ramachandran filter (`pipeline/lab_modified_residues.py`). |
| **cctbx ramalyze 2026.9, run by us on the top8000 files** | [M] Reports MSE (8,422 of 8,458) and essentially no other non-standard residue (MLY 0 of 277, CSO 0 of 48) (`ramalyze_compare.py`). |
| **Agreement** | Same residue-type set as the lab in both generations. |

**Check:** the Top2018 database has 8,852 MSE passing its mainchain filter; our top2018_mc masks (B1 fix) hold 8,852 MSE. [M] (`lab_modified_residues.py`; `verify_hetatm_masks.py`) The records the v1.x residue set drops from top2018_mc (MLY 256, CSD 44, PCA 39, M3L 7) equal the Top2018 database's mainchain-passing counts for those types (MLY 256, CSD 44, PCA 39, M3L 7).

Records dropped by the v1.x residue set (pydangle 664e87a output with the B1 masks; `run_b2_regress_msi.sh`, `check_residue_set.py`; no other record or field changes):

| Dataset | Dropped | By type |
|---|--:|---|
| top100 | 2 | PCA 2 |
| top500 | 12 | CSD 7, PCA 5 |
| top8000 | 279 | MLY 212, CSD 36, PCA 30, M3L 1 |
| top2018_full | 100 | MLY 75, PCA 19, CSD 4, M3L 2 |
| top2018_mc | 346 | MLY 256, CSD 44, PCA 39, M3L 7 |

**Microheterogeneity:** at 1din A123 (top500) and 1ysl B111 (top8000) the file holds CSD and Cys as alternates at one position, and pydangle writes two records with the same chain, number and insertion code, one with φ and one with ψ. These two CSD records are the "one non-standard residue (CSD)" of the v1.0 top500 and top8000 documentation. The v1.x residue set keeps only the Cys record; how such sites should be reported belongs to the altloc decision (section 5).

## 2. φ and ψ of residues next to a modified residue (other than MSE)

| | |
|---|---|
| **These datasets (v1.x)** | Measured through the bonded modified residue, using its atoms, whenever the peptide bond is present (C–N ≤ 2.0 Å). Java Dangle 1.07 gives the same values (pydangle continuity-fix regression, A3). Records of the modified residue itself are dropped (section 1), not its neighbours. |
| **Top2018 database** | Null. [M] Standard residues bonded (table `connectivity`) to a non-MSE non-standard residue: φ present for 13 of 1,011 residues after one, ψ for 13 of 937 before one. Next to MSE: 11,129 of 11,198 and 10,637 of 10,900 (`lab_modified_residues.py`). |
| **top8000 lab lists** | Kept, as far as we can tell. [I] Restricting our best rule to the residues current cctbx ramalyze reports (which, like the Top2018 database, gives no φ/ψ next to non-MSE modified residues) breaks the Gly and Pro classes, which otherwise match the lab's lists exactly (`ramalyze_gaps.py`). So the 2013 lists include residues next to non-standard residues. |
| **Difference** | Ours agrees with the top8000 lists (inferred) and with Java Dangle, and differs from the Top2018 database, where these neighbour values are null. [U] Why the Top2018 database leaves them null; we have not read the cctbx source for this. |

## 3. Terminal caps (ACE, FOR, NH2, NME)

| | |
|---|---|
| **These datasets** | No records. pydangle's residue list names the caps, but they have no CA atom and BioPython's peptide builder leaves them out, so pydangle never reports them (B1 regression: 19 caps in top2018_full and 29 in top2018_mc masks, none in the output). |
| **Top2018 database** | [M] No cap residues in the `residues` table (`lab_modified_residues.py`). |
| **Agreement** | Same. |

## 4. Residue-level B-factor filter

| | |
|---|---|
| **These datasets** | top8000 and top500: every mainchain atom (N, CA, C, O) of the residue has B ≤ 30, and at least three of the four are present; top100: B ≤ 40 (`build_bfactor_masks.py`). The rule is applied to the residue itself; a residue whose i−1 (i+1) neighbour fails keeps its record but φ, ω, `peptide_bond` (ψ) and the Ramachandran category are set to null (`filter_pruned_residues.py`). MSE is tested on N, CA, C, O like any residue. top2018: no B rule of ours; the residue set is the lab's pruned files (section 7). |
| **top8000 lab lists** | [M] `get_srcdata.py`: `1 < max_mc_b < 30`, and φ, ψ, ω, τ defined, so a residue with an unusable neighbour has no row at all. [M] For MSE the B tested is the maximum over N, CA, C, O and SE (1,820 of 1,820; `mse_bfactor.py`). [P] WHM+18 p.295 and p.312: "we excluded residues having any mainchain atom with a B-factor > 30". [M] The archived `.ndft` titles say "B<30"; the CisPro/TransPro fits support < 30 (`b30_boundary.py`). [L] Cutoff 30 on mainchain N, CA, C, O; later versions (top2018 filtered) add CB. |
| **Differences (top8000)** | B = 30.00 is kept by us and excluded by the lab (340 residues in the lab's classes have own mainchain max B = 30.00; reconciliation README, "Published methodology"); B ≤ 1 is kept by us and excluded by the lab; MSE side-chain SE is tested by the lab, not by us; a residue next to a failing residue is a record with null φ or ψ for us and absent for the lab. |
| **top100, top500** | The READMEs cite Word et al. 1999a (top100, B ≤ 40) and Lovell et al. 2003 (top500, "residues with B < 30"). [U] Neither rule has been checked against a per-residue reference list; the quotations have not been re-checked against the papers for this note. Note "B < 30" (Lovell 2003 as quoted) versus our B ≤ 30. |

## 5. Alternate conformations (altlocs)

| | |
|---|---|
| **These datasets** | pydangle uses BioPython's default: for each atom separately, the alternate with the highest occupancy (ties: first listed). No residue is excluded for having alternates. Decision for v1.x still open (register B6). |
| **top8000 lab lists** | [M] Almost every residue with altloc records was left out: 1,504,264 of the 1,504,897 listed residues have no altloc label and occupancy 1 (`lab_altloc_use.py`). [M] `get_srcdata.py` requires `min_res_q = 1`. [M] Non-Pro residues whose preceding C alternates lie at different positions were dropped (10,228 of 10,263), those whose preceding C alternates coincide were kept (14,172 of 14,285) (`lab2013_compare.py`, `prevC_spread.py`). [L] Altloc handling on adjoining residues was added to ramalyze after the top8000 analysis; the old version probably took the first alternate of each neighbour. [L] (J. S. Richardson, pers. comm. 2026-10-02) residues with alternates were dropped. |
| **Top2018 database** | [M] φ/ψ from altloc A / the first conformer: in 552 values where the rules disagree (300 random entries with alternates), 551 match altloc A and none highest occupancy (`altloc_vs_williams.py`, a80f002). [M] Conformer kept consistent across neighbours (6b6w B301, 5j1s B431). |
| **Effect** | top2018 filtered sets: 57 (full) and 78 (mc) passing residues have alternates; the rules differ for 2 ψ values. top8000: we keep residues the lab dropped (no occupancy or altloc filter in our masks). |

## 6. ω, cis/trans and twisted peptides

| | |
|---|---|
| **These datasets** | No ω filter. Pro is CisPro if \|ω\| < 30°, otherwise TransPro, so twisted Pro (30° ≤ \|ω\| < 150°) are labelled TransPro. Non-Pro residues after a cis or twisted peptide stay in their class. Decision for v1.x still open (register B5). |
| **top8000 lab lists** | [M] `get_srcdata.py`: every class except CisPro requires \|ω\| > 150°; CisPro requires \|ω\| < 30°; twisted peptides are in no class. [P] WHM+18 p.305: cis −30° to +30°, trans beyond ±150°, twisted otherwise. [L] Not known for the original build; the lab now thinks twisted peptides should be excluded from filtered sets. |
| **cctbx ramalyze 2026.9** | [M] Splits cis/trans Pro at \|ω\| = 90° (30 twisted Pro in top8000 classed differently; `ramalyze_compare.py`). |

## 7. Ramachandran category assignment

| | |
|---|---|
| **These datasets** | Gly, Pro > pre-Pro > Ile/Val > General (pydangle 664e87a, register A2). |
| **Lab** | [P] RAA+11 p.1399: "Gly, Pro > pre-Pro > Ile/Val > general". [M] `get_srcdata.py`: pre-Pro Ile/Val are PrePro; Gly includes pre-Pro Gly. |
| **Agreement** | Same since 664e87a. The v1.0 datasets put Ile/Val before Pro in IleVal. |

## 8. Coordinates and residue selection

| | |
|---|---|
| **top8000 (ours)** | Current RCSB files, Reduce for Asn/Gln/His flips; the B filter of section 4. |
| **top8000 lab lists** | [M] Built from the lab's 2011 chain files (Reduce 3.15.091106; flips are name swaps, heavy-atom coordinates unchanged): where the 2011 files and our downloads differ, the lab's φ/ψ match the 2011 files in 1,039 of 1,039 residues (`lab_file_provenance.py`). [M] Mainchain records: 7,543,539 of 7,549,697 paired records identical; 4 entries (2oyp, 3h87, 3hkx, 3me5) carry different residue types at the same numbers (`mc_atom_compare.py`, `mc_unmatched_pairing.py`). [L] (Oct 2026) Filtering used deposited heavy-atom data only, not Reduce hydrogens. |
| **top2018 (ours)** | Residue selection is the lab's: our masks are the residues kept in the lab's pruned files, and fragment boundaries are the lab's `USER INC` records. Measurement is on current RCSB files with Reduce. Where RCSB has renamed or renumbered since the lab's files, we remap (4ct3 chain D → E, 5tvo +315, 1xmz B ≥ 66 → C; `build_ersatz_pdbs.py --rcsb-remap`). 5 residues have a different type in the current RCSB file (2omk A39, 3esm A150: lab ACY, RCSB Gly; 4oe9 B71 MSE → Ala, B106 MSE → Ser). |
| **top2018 lab filter** | [L] Mainchain criterion adds CB where present. We do not re-apply any lab filter; we use its output. |

## 9. Chain breaks and isolated residues

| | |
|---|---|
| **These datasets** | A peptide bond exists when C(i)–N(i+1) ≤ 2.0 Å (Java Dangle's `AtomSpec.checkConnection` test). A residue with no bond to either neighbour is not reported (decision 2026-11-08). |
| **Top2018 database** | [M] Table `connectivity` gives prev/next residues (documented as sequence adjacency plus peptide-bond distance). [U] Its distance threshold. |
| **top8000 lab lists** | [U] The bond test of the tool that filled the lab's `dangle` table. [M] Java Dangle 1.07 and our parser give the same φ/ψ coverage on the top8000 files except at chain breaks with C–N > 2 Å (`dangle_compare.py`). |

## 10. Measurement tool

| | |
|---|---|
| **These datasets** | pydangle-biopython (BioPython). Agrees with Java Dangle 1.07 except where the two pick different alternates (section 5; over all residues in the top8000 ersatz files, 15,994 φ/ψ/ω values, `altloc_java_diffs.py`) and where altloc labels do not match across a bond (3rdo A−5; pydangle does not link it, Java Dangle does). |
| **Lab** | [L] top8000: probably an early ramalyze; earlier datasets possibly Dangle or Prekin. [M] The 2013 lists' SQL reads a table named `dangle` (`get_srcdata.py`). [L] (M. G. Prisant, 2026-11-08) reference angles from Top8000 onward come from cctbx. [U] Exact versions. |

## 11. Side-chain χ angles

| | |
|---|---|
| **These datasets** | Raw dihedral to the first-named atom (e.g. Phe/Tyr χ2 to CD1), −180° to 180°, not folded for symmetric side chains (register B4). |
| **Top2018 database** | [M] The same raw dihedrals on 0° to 360°, also not folded (Asp χ2: 135,617 of 200,806 values ≥ 180°). |
| **Difference** | Range only. [U] Whether MolProbity rotamer scoring folds internally (not checked against source). |

## 12. DSSP

| | |
|---|---|
| **These datasets** | mkdssp 4.2.2 on the ersatz files. |
| **Lab** | [M] The Top2018 database has no DSSP. No lab DSSP reference is used here. |

## Open questions for the lab

1. The φ/ψ handling next to non-MSE modified residues (section 2): is the null in the Top2018 database intended, and which convention should v1.x follow? v1.x currently keeps the values.
2. Twisted peptides (section 6) and altloc choice (section 5): decisions pending (register B5, B6).
3. The top8000 residue-level differences (section 4) are documented, not adopted: v1.x top8000 keeps our filter (register B3).

## Scripts

This repository, `main`: `pipeline/lab_modified_residues.py`, `pipeline/verify_hetatm_masks.py`, `pipeline/filter_pruned_residues.py`, `pipeline/check_residue_set.py`, `pipeline/altloc_vs_williams.py`, `pipeline/altloc_java_diffs.py`. Branch `rama8000-reconciliation` (`e00fa7b`), `rama8000_reconciliation/scripts/`: `best_rule_summary.py`, `nonstandard_types.py`, `ramalyze_compare.py`, `ramalyze_gaps.py`, `b30_boundary.py`, `mse_bfactor.py`, `lab_altloc_use.py`, `lab2013_compare.py`, `prevC_spread.py`, `lab_file_provenance.py`, `mc_atom_compare.py`, `mc_unmatched_pairing.py`, `dangle_compare.py`.
