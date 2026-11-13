# Design choices and differences from the Richardson Lab methodology

*File: differences_from_richardson_lab.md · Git: Michael G. Prisant 038e7f7 · Txd: 11/11/26 08:28 with local modifications!*

**Date:** 2026-11-07 onward; restructured as a design-choice register 2026-11-09.
**Status:** DRAFT for review by M. G. Prisant. Statements about the Richardson Lab's procedures must be confirmed before this file is merged or cited.
**Scope:** every design choice in the pydangle measures databases of this repository (top100, top500, top8000, top2018_full, top2018_mc; v1.x pipeline) that bears on comparability with the Richardson Lab's original analyses, with what we understand the original work did, the difference, and its measured consequence. A running list: add an entry whenever a choice is made or changed.

## How to read this file

Each entry gives the **status** (decided, with date and by whom, or open), **our choice**, the **original** as we understand it, the **difference**, and the **consequence** for our datasets.

Each statement about the lab's work carries its evidence:

| Tag | Meaning |
|---|---|
| **[M]** | Measured from the lab's own data or files by a committed script (script named) |
| **[P]** | Stated in a publication (page given) |
| **[L]** | Said by Richardson Lab members (online meeting 2026-10-14: J. S. Richardson, C. J. Williams, V. B. Chen; or the dated personal communication) |
| **[I]** | Inferred by us from counts; not confirmed by the lab |
| **[U]** | Unknown |

Statements about our own pipeline come from the code in this repository and in pydangle-biopython and are not tagged.

Reference material:

- **top8000:** the lab's per-residue Ramachandran lists `Top8000_Rama_2013.zip` (written 2013-02-16 by `get_srcdata.py` from the lab's top8000 MySQL database; Silk on these lists reproduces all six archived rama8000 `.ndft` tables bin for bin), and the lab's 2011 chain files `top8000_chains_70.tgz`, the input of those lists. Analysis on branch `rama8000-reconciliation` (`e00fa7b`), directory `rama8000_reconciliation/`.
- **top2018:** the lab's pruned files (all-atom v2.0, Zenodo 10.5281/zenodo.5773255; mainchain v2.01, Zenodo 10.5281/zenodo.5777651) and C. J. Williams' Top2018 SQLite database of validation and geometry measures (Zenodo 10.5281/zenodo.22085722).
- **top100, top500:** no per-residue reference list recovered; only the publications cited in each README.
- **Publications:** Williams et al. 2018, *Protein Sci* 27:293–315 (WHM+18); Read et al. 2011, *Structure* 19:1395–1412 (RAA+11).

## Summary

| # | Choice | Original (as understood) | Same? | Status |
|---|---|---|---|---|
| DC-1 | Current RCSB coordinates, Reduce flips | top8000: lab's 2011 chain files | Inputs differ in a few entries | decided (v1.0) |
| DC-2 | top2018 residues = the lab's pruned files | the lab's pruned files | Same selection | decided (v1.0) |
| DC-3 | Mainchain B ≤ 30 (top100: ≤ 40), own residue | top8000: 1 < B < 30, MSE incl. SE | Differs at the edges | open |
| DC-4 | Neighbour fields nulled at fragment ends, record kept | top8000: residue absent | Differs | decided (v1.0) |
| DC-5 | pydangle-biopython | cctbx (probably via ramalyze or omegalyze); Prekin or Dangle before Top8000 | Different tool | decided |
| DC-6 | Peptide bond C–N ≤ 2.0 Å; isolated residues dropped | threshold unknown | Unknown | decided 2026-11-08 |
| DC-7 | Altloc: first conformer (per-residue file order), one per residue | Top2018: first conformer, consistent; top8000: altloc residues dropped | Top2018: 12 values in 7 entries differ | decided 2026-11-09, revised 2026-11-10 |
| DC-8 | Microheterogeneity: first standard-residue conformer | Top2018: standard conformer measured, modified name kept | Name only | decided 2026-11-09, revised 2026-11-10 |
| DC-9 | Records for std20 + MSE only | std20 + MSE | Same | decided 2026-11-09 |
| DC-10 | φ/ψ next to modified residues measured through them | top8000: kept [I]; Top2018: null | Differs from Top2018 | decided 2026-11-09 |
| DC-11 | Caps never reported | Top2018: no caps | Same | — |
| DC-12 | Gly, Pro > pre-Pro > Ile/Val > General | same (RAA+11) | Same since v1.x | decided 2026-11-08 |
| DC-13 | No ω filter; twisted Pro labelled TransPro | top8000: twisted excluded | Differs | open |
| DC-14 | χ raw dihedral, −180..180, not folded | Top2018: raw, 0..360, not folded | Range only | open |
| DC-15 | mkdssp 4.2.2, PDB input | no lab DSSP reference | — | decided 2026-09-25 |

---

### DC-1. Coordinates

- **Status:** decided (v1.0).
- **Our choice:** current RCSB files; Reduce for Asn/Gln/His flips only (heavy atoms otherwise as deposited).
- **Original:** [M] the top8000 lists were built from the lab's 2011 chain files (Reduce 3.15.091106; flips are name swaps, heavy-atom coordinates unchanged): where the 2011 files and our downloads differ, the lab's φ/ψ match the 2011 files in 1,039 of 1,039 residues (`lab_file_provenance.py`). [L] (Oct 2026) filtering used deposited heavy-atom data only, not Reduce hydrogens.
- **Difference:** [M] mainchain records, 2011 files vs our downloads: 7,543,539 of 7,549,697 paired records identical; 4 entries (2oyp, 3h87, 3hkx, 3me5) carry different residue types at the same numbers (`mc_atom_compare.py`, `mc_unmatched_pairing.py`).
- **Consequence:** residue-level comparisons with the lab's top8000 lists differ in those entries and in entries re-refined since 2011 (e.g. 3kda, 1gug: φ/ψ up to 49°).

### DC-2. top2018 residue selection

- **Status:** decided (v1.0); remaps scripted 2026-11-08.
- **Our choice:** the residue set is the lab's: masks are the residues kept in the lab's pruned files (ATOM and, since 2026-11-08, polymer HETATM residues), and fragment boundaries are the lab's `USER INC` records. Measurement is on current RCSB files. Where RCSB has renamed or renumbered since the lab's files, we remap (4ct3 chain D → E, 5tvo +315, 1xmz B ≥ 66 → C; `build_ersatz_pdbs.py --rcsb-remap`).
- **Original:** [L] the lab's mainchain filter adds CB where present. We do not re-apply any lab filter; we use its output.
- **Difference:** none in selection. 5 residues have a different type in the current RCSB file (2omk A39, 3esm A150: lab ACY, RCSB Gly; 4oe9 B71 MSE → Ala, B106 MSE → Ser).
- **Consequence:** [M] the Top2018 database's 8,852 MSE passing its mainchain filter equal the 8,852 MSE in our top2018_mc masks (`lab_modified_residues.py`, `verify_hetatm_masks.py`).

### DC-3. Residue-level B-factor filter (top100, top500, top8000)

- **Status:** open; v1.x keeps our rule unless decided otherwise.
- **Our choice:** every mainchain atom (N, CA, C, O) of the residue has B ≤ 30 (top100: ≤ 40), at least three of the four present (`build_bfactor_masks.py`). MSE is tested on N, CA, C, O like any residue. No occupancy or altloc filter.
- **Original:** [M] top8000 lists, `get_srcdata.py`: `1 < max_mc_b < 30`; `min_res_q = 1`. [M] MSE: B is the maximum over N, CA, C, O and SE (1,820 of 1,820; `mse_bfactor.py`). [P] WHM+18 p.295, p.312: "we excluded residues having any mainchain atom with a B-factor > 30". [M] archived `.ndft` titles "B<30"; CisPro/TransPro fits support < 30 (`b30_boundary.py`). [L] cutoff 30 on mainchain N, CA, C, O. [M] Almost every residue with altloc records was left out (1,504,264 of 1,504,897 listed residues have no altloc label and occupancy 1; `lab_altloc_use.py`), and non-Pro residues whose preceding C alternates lie at different positions were dropped (10,228 of 10,263; `lab2013_compare.py`, `prevC_spread.py`). top100, top500: the READMEs cite Word et al. 1999a (B ≤ 40) and Lovell et al. 2003 ("residues with B < 30"); [U] not checked against a per-residue list, quotations not re-checked.
- **Difference:** B = 30.00 kept by us, excluded by the lab (340 residues in the lab's classes); B ≤ 1 kept by us; MSE SE not tested by us; residues with alternates kept by us.
- **Consequence:** top8000 keeps residues the lab's lists exclude (B at the edges, residues with alternates, MSE with high side-chain B); per-class counts against the lab's lists: rama8000 reconciliation README. The angles of residues with alternates are DC-7.

### DC-4. Residues next to a residue that fails the filter

- **Status:** decided (v1.0).
- **Our choice:** the residue keeps its record; φ, ω, `peptide_bond` (failing i−1) or ψ (failing i+1), and the Ramachandran category, are set to null (`filter_pruned_residues.py`). τ, χ, chirality, DSSP are kept.
- **Original:** [M] top8000 lists: φ, ψ, ω, τ must be defined (`get_srcdata.py`), so such a residue has no row.
- **Difference:** a record with null neighbour fields versus no row.
- **Consequence:** fewer φ/ψ pairs than records; any comparison with the lab's lists must count residues with both φ and ψ.

### DC-5. Measurement tool

- **Status:** decided.
- **Our choice:** pydangle-biopython (BioPython parser). Agrees with Java Dangle 1.07 except at altloc choice (DC-7), in the third decimal where BioPython's single-precision coordinates matter (DC-7), and at mismatched altloc labels across a bond (3rdo A−5: pydangle does not link it, Java Dangle does).
- **Original:** [L] (M. G. Prisant, 2026-11-12) Top8000 and later: cctbx, probably called by ramalyze or omegalyze; before Top8000 (Top100, Top500): probably Prekin or Dangle. [M] the 2013 lists' SQL reads a table named `dangle` (`get_srcdata.py`). [U] exact versions.
- **Difference / consequence:** see DC-6 to DC-10.

### DC-6. Peptide bonds and isolated residues

- **Status:** decided 2026-11-08 (isolated residues: M. G. Prisant).
- **Our choice:** a peptide bond exists when C(i)–N(i+1) ≤ 2.0 Å (Java Dangle's `AtomSpec.checkConnection`), for any compatible altloc pair. A residue with no bond to either neighbour is not reported.
- **Original:** [M] Top2018 database: table `connectivity` gives prev/next residues (documented as sequence adjacency plus peptide-bond distance). [U] its threshold, and the bond test of the tool that filled the top8000 `dangle` table. [M] Java Dangle 1.07 and our parser give the same φ/ψ coverage on the top8000 files except at chain breaks with C–N > 2 Å (`dangle_compare.py`).
- **Consequence:** 6–10 isolated residues per dataset not reported (Java Dangle reports them with τ only).

### DC-7. Alternate conformations (altlocs)

- **Status:** decided 2026-11-09 (M. G. Prisant), ordering rule revised 2026-11-10 (M. G. Prisant); pydangle branch `altloc-first`, default `--altloc first`; `--altloc occupancy` keeps the earlier behaviour.
- **Our choice:** the first conformer. "First" is per residue: among the alternates of the atoms concerned, the label written first in the file (lowest atom serial; normally A). Each residue is measured in the first conformer shared by every atom of its φ/ψ window C(i−1), N, CA, C, N(i+1); a measurement whose atoms lack that conformer uses the first conformer its own atoms share, else each atom's own first conformer. So a residue present only as B, and the residue after it, are measured in B.
- **Rejected (2026-11-10):** ranking labels by first appearance anywhere in the file. In 3zeu the first alternates in the file are A1 Met conformers B and C (no A), so that rule measured B throughout the entry, where the lab measured A (3zeu B168: lab φ −67.634, ψ −44.071 = conformer A). [M] Over every entry below it left 311 values differing near alternates; alphabetical order left 38; the per-residue order 36, before DC-8's revision (`pydangle_vs_williams.py`).
- **Before v1.x:** BioPython's per-atom highest occupancy (ties: first listed), so one dihedral could mix two conformers.
- **Original:** [M] Top2018 database: the first conformer, the same across the φ/ψ window (6b6w B301–B302, 5j1s B431–B432, 4qb6 A57, 3zeu B168). [M] top8000 lists: residues with alternates dropped (DC-3). [L] the old ramalyze probably took the first alternate of each neighbour.
- **Difference:** [M] every Top2018 entry with alternates (9,095; `pydangle_vs_williams.py --altloc first`, with DC-8): φ/ψ of residues with an alternate in the window agree within 0.01° for 451,148 of 451,186 values. Of the other 38, 26 are in two entries the lab numbers differently from our files (4oe9 chain B offset 5, 4z54 chain A offset 38; DC-1), and 12 are in 7 entries where the lab used conformer B and we used A: 4hvy A163, 6fj7 A18, 5ice A10 (only CA has alternates), 3zoj A194, 2q8k A40 (only the preceding C has alternates), 5ahn A304, 3eoj A212 (whole window; in 5ahn A is written first and both have occupancy 0.50). Each lab value equals the all-B value in our coordinates (`altloc_dihedrals.py`); why B is chosen there is [U]. Residues without alternates: 4,249,028 values agree; the 740 that do not are all in 4oe9 and 4z54.
- **Java Dangle 1.07:** [M] φ/ψ/ω disagreements (`java_diff_classes.py`, `run_java_classes_msi.sh`). Before = `--altloc occupancy` (the earlier behaviour), after = `--altloc first`; values differing by more than 0.002° (both programs print three decimals), every residue, model 1; top100 has no Java Dangle run. Each disagreement is recomputed from the coordinates for every combination of its four atoms' alternates: *precision* = the same atoms give both values within 0.05°; *alternate choice* = the two values come from different alternates. No disagreement falls outside these two classes.

  **Table 1.** φ/ψ/ω values where pydangle and Java Dangle 1.07 differ by more than 0.002° (model 1). Before = `--altloc occupancy`, the earlier behaviour; after = `--altloc first`.

  | Dataset | Values compared | Before | After |
  |---|--:|--:|--:|
  | top500 | 501,352 | 391 | 10 |
  | top8000 | 11,080,899 | 15,994 | 1,579 |
  | top2018_full | 18,349,291 | 65,015 | 2,420 |
  | top2018_mc | 20,325,102 | 73,864 | 3,283 |

  **Table 2.** Classification of the disagreements in Table 1. Each cell is *before / after*. Precision: the same atoms give both values within 0.05°. Alternate choice: the two values come from different alternates. Total = precision + alternate choice; no disagreement falls outside these two classes.

  | Dataset | Precision | Alternate choice | Total |
  |---|--:|--:|--:|
  | top500 | 8 / 4 | 383 / 6 | 391 / 10 |
  | top8000 | 1,743 / 1,464 | 14,251 / 115 | 15,994 / 1,579 |
  | top2018_full | 3,725 / 2,115 | 61,290 / 305 | 65,015 / 2,420 |
  | top2018_mc | 4,751 / 2,967 | 69,113 / 316 | 73,864 / 3,283 |

  After the change, precision: pydangle is the program off in 4 / 1,450 / 2,097 / 2,949 (top500 / top8000 / top2018_full / top2018_mc): recomputing from coordinates rounded to single precision, as BioPython stores them, reproduces pydangle's value, and Java's is the double-precision value (1gdj A3, C–N–CA ≈ 177°: Java φ 47.971, coordinates 47.9705, pydangle 47.965); Java is off in 0 / 5 / 9 / 9; 0 / 9 / 9 / 9 are within 0.05° with cause not established. Alternate choice: microheterogeneity 2 / 48 / 185 / 191 (including positions where pydangle keeps the standard type, DC-8, and Java the first, e.g. 2vqr A56 next to A57, conformer A DDZ; 1ejg A23 next to A22 Pro/Ser); Java mixing conformers within one dihedral 3 / 37 / 55 / 56 (3sil A5: A4 exists only as B; pydangle φ −115.14, all B, the 6b6w rule; Java 120.86 from B, A, A, A); a different single conformer 1 / 30 / 65 / 69. Before, alternate choice is almost all the occupancy rule taking B where Java takes A (378 / 14,229 / 61,230 / 69,052).
- **Consequence:** [M] `--altloc occupancy` reproduces the pre-v1.x output byte for byte on all five datasets. Records changed by `--altloc first` (any field; `run_altloc_regress_msi.sh`, `compare_altloc.py`, TAG `_v2`), all residues / after the filter (masks with polymer HETATM residues, std20+MSE): top100 68 / 62; top500 481 / 336; top8000 8,464 / 4,674; top2018_full 32,559 / 11; top2018_mc 35,165 / 17. φ values changed, all residues: top100 7, top500 147, top8000 4,394, top2018_full 18,523, top2018_mc 20,941. In the filtered top2018 sets almost nothing changes (mostly ω by under 5°); [I] the lab's top2018 residue selection already excludes most residues with mainchain alternates.

### DC-8. Microheterogeneity (two residue types at one position)

- **Status:** decided 2026-11-09 (with DC-7), revised 2026-11-10 (M. G. Prisant): standard residue type preferred.
- **Our choice:** keep one residue type per position: the first conformer whose type is a standard residue (20 amino acids or MSE, the DC-9 set), or the first conformer if none is. Handles BioPython's `DisorderedResidue` (both types ATOM, e.g. 2zs0 D31 Ile/Val) and split pairs (one type HETATM, e.g. 1din A123 CSD/Cys, 5k0a A157 CSO/Cys, 1xpm D111 SCY/Cys, 3zeu D194 MSE/Met). Where the first type has an atom with no altloc label, BioPython's parser keeps only that type in both modes (1ejg A22 Pro/Ser keeps the Pro).
- **Rejected (2026-11-10):** the strict first conformer's type. Where that is a modified residue, pydangle does not report it (DC-9) and its neighbours were measured through the modified conformer.
- **Before v1.x:** split pairs gave two records with the same chain, number and insertion code, each with half the neighbour angles (1din A123, 1ysl B111 in top500/top8000; these were the "one non-standard residue (CSD)" of the v1.0 docs); `DisorderedResidue` sites gave the higher-occupancy type (2zs0 D31: Val; Java Dangle: Ile).
- **Original:** [M] Top2018 database: at the 13 split sites where conformer A is a modified residue and B a standard one (1xpm, 3fcr, 3o2r, 3piu, 3sqz, 4dpz, 4kuk, 4upi, 4uup, 4z0y, 5jfn, 5k0a, 2vqr), φ/ψ are those of the standard conformer, but the residue is named as the modified type (e.g. 5k0a A157 "CSO", φ −137.641, ψ 111.071 = the Cys conformer; the CSO conformer gives −122.152, 129.179). [U] top8000 and earlier.
- **Difference:** [M] values: none, 26 of 26 φ/ψ at those sites agree within 0.01° (`pydangle_vs_williams.py`). Name: we label the residue with the type measured (CYS, LYS, ...); the lab's database labels it with the modified type. [M] These 13 are exactly the non-MSE modified residues with Ramachandran rows in the Top2018 database (DC-9): the database gives φ/ψ to a modified residue other than MSE only at such sites, and then measures the standard conformer.
- **Consequence:** [M] residue types that appear or disappear at microheterogeneity positions (`compare_altloc.py`, all residues): top100 none; top500 6 / 5 (e.g. Val→Ile, 1din A123 CSD→Cys); top8000 26 / 26; top2018_full 88 / 61; top2018_mc 88 / 61. After the filter: top500 2 / 2, top8000 10 / 11, top2018 none.

### DC-9. Which residue types get records

- **Status:** decided 2026-11-09 (M. G. Prisant).
- **Our choice:** the 20 standard amino acids and MSE (`filter_pruned_residues.py --residues std20+MSE`, every `reproduce.sh`). MSE keeps the name MSE. pydangle itself still reports its Java-Dangle list (also MLY, M3L, CSD, CYO, PCA, AIB, ABU, DGN, ASX, GLX); the restriction is a dataset choice.
- **Original:** [M] top8000 lists: std + MSE (3,440 non-standard residues pass the best-fitting rule, 2,985 MSE; adding the other 455 worsens the fit; `best_rule_summary.py`, `nonstandard_types.py`). [L] "best guess: only MSE included; not certain." [M] Top2018 database: MSE 11,116 of 11,532 with a Ramachandran row; the other 110 non-standard types 13 of 1,045, none passing the Ramachandran filter (`lab_modified_residues.py`). [M] cctbx ramalyze 2026.9 on the top8000 files: MSE 8,422 of 8,458, MLY 0 of 277, CSO 0 of 48 (`ramalyze_compare.py`).
- **Difference:** none.
- **Consequence:** [M] records dropped (`run_b2_regress_msi.sh`, `check_residue_set.py`): top100 2, top500 12, top8000 279 (MLY 212, CSD 36, PCA 30, M3L 1), top2018_full 100, top2018_mc 346 (MLY 256, CSD 44, PCA 39, M3L 7 = the Top2018 database's mainchain-passing counts of those types). v1.0 had almost no modified residues because of the ATOM-only masks (fixed 2026-11-08).

### DC-10. φ and ψ next to a modified residue (other than MSE)

- **Status:** decided 2026-11-09 (with DC-9); open question to the lab.
- **Our choice:** measured through the bonded modified residue, using its atoms (Java Dangle 1.07 gives the same values). Only the modified residue's own record is dropped (DC-9).
- **Original:** [M] Top2018 database: null. Standard residues bonded (`connectivity`) to a non-MSE non-standard residue have φ in 13 of 1,011 cases (residue after) and ψ in 13 of 937 (residue before); next to MSE 11,129 of 11,198 and 10,637 of 10,900 (`lab_modified_residues.py`). [I] top8000 lists: kept: restricting our best rule to the residues current ramalyze reports breaks the Gly and Pro classes, which otherwise match the lab's lists exactly (`ramalyze_gaps.py`). [U] why the Top2018 database leaves them null (cctbx source not read).
- **Difference:** agrees with the top8000 lists (inferred) and Java Dangle; differs from the Top2018 database.
- **Consequence:** over all Top2018 chains before filtering, 998 φ and 924 ψ values that the Top2018 database leaves null are measured by us (1,011 − 13 and 937 − 13, `lab_modified_residues.py`).

### DC-11. Terminal caps (ACE, FOR, NH2, NME)

- **Status:** no choice needed.
- **Our choice:** never reported: BioPython's peptide builder accepts a residue only if `is_aa` knows the name or it has a CA atom; caps pass neither (HETATM mask fix regression: 19 caps in top2018_full and 29 in top2018_mc masks, none in the output).
- **Original:** [M] the Top2018 database has no cap residues (`lab_modified_residues.py`).
- **Difference:** none.

### DC-12. Ramachandran category precedence

- **Status:** decided 2026-11-08 (pydangle 664e87a).
- **Our choice:** Gly, Pro > pre-Pro > Ile/Val > General.
- **Original:** [P] RAA+11 p.1399: "Gly, Pro > pre-Pro > Ile/Val > general". [M] `get_srcdata.py`: pre-Pro Ile/Val are PrePro; Gly includes pre-Pro Gly.
- **Difference:** none since v1.x. v1.0 put Ile/Val before Pro in IleVal (filtered top8000: 9,283 residues).

### DC-13. ω, cis/trans and twisted peptides

- **Status:** open.
- **Our choice:** no ω filter. Pro is CisPro if |ω| < 30°, otherwise TransPro, so twisted Pro (30° ≤ |ω| < 150°) are TransPro. Non-Pro residues after a cis or twisted peptide stay in their class.
- **Original:** [M] top8000 lists, `get_srcdata.py`: every class except CisPro requires |ω| > 150°, CisPro |ω| < 30°; twisted peptides in no class. [P] WHM+18 p.305: cis −30° to +30°, trans beyond ±150°, twisted otherwise. [L] not known for the original build; the lab now thinks twisted peptides should be excluded. [M] cctbx ramalyze 2026.9 splits cis/trans Pro at |ω| = 90° (`ramalyze_compare.py`).
- **Difference:** twisted Pro labelled TransPro; residues after twisted or cis peptides kept.

### DC-14. Side-chain χ angles

- **Status:** open (stay unfolded for now, 2026-11-07).
- **Our choice:** raw dihedral to the first-named atom (e.g. Phe/Tyr χ2 to CD1), −180° to 180°, not folded for symmetric side chains.
- **Original:** [M] Top2018 database: the same raw dihedrals on 0° to 360°, also not folded (Asp χ2: 135,617 of 200,806 values ≥ 180°). [U] whether MolProbity rotamer scoring folds internally.
- **Difference:** range only.

### DC-15. DSSP

- **Status:** decided 2026-09-25 (v1.x on mkdssp 4.2.2, PDB input).
- **Our choice:** mkdssp 4.2.2 on the ersatz files.
- **Original:** [M] the Top2018 database has no DSSP; no lab DSSP reference is used.

## Open questions for the lab

1. DC-10: is the null φ/ψ next to non-MSE modified residues in the Top2018 database intended, and which convention should v1.x follow? v1.x keeps the values.
2. DC-13 and DC-3: twisted peptides and the top8000 residue-level rules.
3. DC-8: at split sites the Top2018 database measures the standard conformer but names the residue as the modified type; is that intended? And how did earlier analyses (top500, top8000) treat these positions?
4. DC-7: why the Top2018 database uses conformer B at 4hvy A163, 6fj7 A18, 5ice A10, 3zoj A194, 2q8k A40, 5ahn A304 and 3eoj A212.
5. DC-1: 4oe9 (chain B) and 4z54 (chain A) are numbered differently in the Top2018 database from the coordinate files we use; which coordinate version did the lab use?

## Scripts

This repository, `main`: `pipeline/lab_modified_residues.py`, `pipeline/verify_hetatm_masks.py`, `pipeline/filter_pruned_residues.py`, `pipeline/check_residue_set.py`, `pipeline/altloc_vs_williams.py`, `pipeline/altloc_java_diffs.py`; altloc checks: `pipeline/pydangle_vs_williams.py`, `pipeline/altloc_dihedrals.py`, `pipeline/compare_altloc.py`, `pipeline/run_altloc_regress_msi.sh`, `pipeline/java_diff_classes.py`, `pipeline/run_java_classes_msi.sh`. Branch `rama8000-reconciliation` (`e00fa7b`), `rama8000_reconciliation/scripts/`: `best_rule_summary.py`, `nonstandard_types.py`, `ramalyze_compare.py`, `ramalyze_gaps.py`, `b30_boundary.py`, `mse_bfactor.py`, `lab_altloc_use.py`, `lab2013_compare.py`, `prevC_spread.py`, `lab_file_provenance.py`, `mc_atom_compare.py`, `mc_unmatched_pairing.py`, `dangle_compare.py`.
