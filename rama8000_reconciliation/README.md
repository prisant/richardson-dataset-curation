# rama8000 reconciliation

Goal: explain the residue-count differences between the published top8000
JSONL (Zenodo 10.5281/zenodo.19051741, v2.0.0) and the Richardson Lab
rama8000 NDFT contour tables
(`rlabduke/javadev@24e0666`, `chiropraxis/resource/chiropraxis/rotarama/rama8000/*.ndft`),
first raised by Erik Novak (2026-10-03).

Branch: `rama8000-reconciliation`.

## Layout

| Path | Contents |
|---|---|
| `silk_build/javadev/` | javadev source at `24e0666` (silk, driftwood, rama8000 resources); tarball sha256 `5d6dfc42…ebac5` |
| `silk_build/classes/` | Silk compiled with OpenJDK 21 (`javac -d classes -sourcepath silk/src:driftwood/src`) |
| `data/` | the six archived `.ndft` files; `pro_records.jsonl` (all Pro with φ,ψ in listed chains, parsed from ersatz PDBs); `cispro_bestfit.txt` |
| `scripts/` | `ndft.py` (reader), `recount3.py` (class recount under filter variants), `pro_records.py`, `silk_fit*.py`, `cis_*.py` |
| `results/` | Silk outputs and per-residue test results |

## Findings so far (CisPro)

1. NDFT `realcount` is the observation count: grid values are exact multiples
   of 1/realcount (Silk `-fraction`). Gly realcount is doubled by `-sym`.
2. Silk parameters that reproduce `cispro.ndft`: `-cosine=5 -twopass=8`
   (2° grid, wrap, `-fraction`). Neighbouring values fit clearly worse.
3. Residue filter that best reproduces the grid: **own mainchain N, CA, C, O
   B < 30; neighbours not checked; residues with any altloc excluded**
   (n = 3,663 vs realcount 3,654). Adding the 54 residues that this rule
   admits but the φ/ψ-atom rule rejects cuts the grid error 1.15 → 0.39
   (×1e-4 mean |Δfraction|) and the support mismatch 73 → 2 bins. Adding
   residues with B just over 30, or with altlocs, makes the fit worse.
   The altloc exclusion is not stated in Williams et al. (2018); Jane
   Richardson confirmed it as part of the rama8000 filtering (personal
   communication to M. Prisant, 2026-10-02).
5. **CisPro reconciled (2026-10-04).** Final rule: own mainchain N, CA, C, O
   B < 30 (neighbours not checked) **and every atom of the residue has
   occupancy 1** (`scripts/cis_edge_groups.py`, `scripts/cis_final_check.py`).
   n = 3,654 = realcount exactly. The occupancy test subsumes the altloc
   drop (no altloc residue has all occupancies 1; the two rules give
   identical sets). It additionally removes 9 residues with partial or zero
   occupancy but no altloc label: 1gtv A37, 1k3i A52/A163/A350, 2gw3 B49/B84,
   2gw4 D84, 2qod A694, 3o0p A221 (sidechain-only partial occupancy, so the
   test covers all atoms, not just mainchain). Silk at `-cosine=5
   -twopass=8` on this set reproduces `cispro.ndft` with identical integer
   rank in 32,090 / 32,400 bins; the rest differ by 1–3 ranks (of 3,654).
   No other edge group tested (neighbour altlocs or occupancy, B in
   [29,30) or [30,30.5), |ω| near 30°, CB or all-atom B, non-standard
   neighbours, long peptide bonds) improves the fit; all worsen it.
   Caveat: set identity is established by count plus grid; swapping two
   residues at nearly identical φ,ψ would be invisible.
6. **TransPro (2026-10-04, in progress).** Same Makefile params
   (`-cosine=5 -twopass=8`; `-twopass=6.5` fits worse). The CisPro rule
   gives n = 66,404 (+40) but leaves 759 bins of support mismatch, caused by
   a few isolated outliers with twisted ω. Defining trans as **|ω| ≥ 150°**
   (twisted 30–150° excluded; plateau 149–152°) gives n = 66,377 (**+13**),
   support mismatch 0, rank-identical 30,482 / 32,400 bins, max |Δrank| 60
   (of 66,364). Edge groups (neighbour altloc/occupancy, B near 30, CB or
   all-atom B, non-standard neighbours, long peptide bonds) do not explain
   the +13. Scripts: `class_rule_check.py`, `transpro_islands.py`,
   `transpro_edges*.py`, `transpro_loo_sparse.py`.
   The +13 is explained by **zero B-factors**: 14 residues (2fcr A ×10,
   2yyh C ×4) have all mainchain B = 0.00 (unrefined/unreported B).
   Dropping them gives n = 66,363 (**−1**), grid error 0.040 → 0.025,
   rank-identical 30,851 bins, max |Δrank| 48 (`transpro_bzero.py`,
   `bzero_variants.py`). "Any mainchain atom B = 0" gives −3 with an
   indistinguishable grid; "any atom B = 0" is ruled out (breaks CisPro,
   3,653). Effective B rule is therefore 0 < B < 30 on mainchain, exact
   form at zero unresolved. Residual: grid still differs by up to 48 ranks
   (of 66,364) in ~1,550 bins, i.e. some same-count membership
   differences remain that the count cannot show.
7. **Non-Pro classes (2026-10-04).** Under the Pro rule (own mainchain
   0 < B < 30, all heavy atoms occ 1, standard residues) General/Gly/IleVal/
   PrePro were each +0.6–0.8 % over. The excess is the **amide H**: with
   Reduce hydrogens present, "every atom occupancy 1" also drops residues
   whose H has occ < 1 (Reduce builds partial-occupancy alternate H next to
   altloc residues: 11,272 residues; plus 179 residues in 21 structures
   deposited with occ-0 riding H). Pro has no amide H, hence exact.
   Result: General −772, Gly +76, IleVal +167, PrePro −39 (all within
   0.1 %). Grids improve (PrePro 1.25 → 0.55, General 0.71 → 0.56 ×1e-4 at
   Makefile params) but remain ~20× worse than TransPro.
   Reduce 3.3 (2016-04-22, `~/bin/reduce_old`) gives identical amide-H
   occupancies to 4.16 on all 7,737 structures (0 of 10,219 flagged
   rule-set residues differ), so the Reduce version is not the residual.
   Scripts: `residue_records.py`, `amideH_records.py`, `amideH_counts.py`,
   `nonpro_grid_check.py`, `reduce_old_amideH.py`.
4. Single-residue add/remove tests do not identify individuals: within a
   dense cluster every point has the same effect. They show a regional
   imbalance instead — S has slightly too many points near (φ −70, ψ −17)
   and too few near (φ −77, ψ +162) relative to the lab's set.
