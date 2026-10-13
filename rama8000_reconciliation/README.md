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

## Published methodology: Williams et al. (2018), checked 2026-10-13

Source: Williams CJ et al., "MolProbity: More and better reference data
for improved all-atom structure validation", *Protein Sci* 27:293–315
(2018), doi:10.1002/pro.3330 (local copy
`~/Desktop/UnFiled261002_RamaRichardson/Refs/WHM+18.pdf`). The residue
filter for the Ramachandran contours is stated twice, in Results (p.295)
and Methods (p.312):

> In preparation of Ramachandran and CaBLAM contours, we excluded residues
> having any mainchain atom with a B-factor > 30. This simple filter
> provides a proxy for local model fit to the electron density, and it
> indeed produced quite clean, reproducible contours.

Comparison with the inferred rules:

| Inferred rule | Paper | Status |
|---|---|---|
| Chain set: Top8000-best_hom70, 7,957 chains | "filtered at the 70% homology level and contains 7957 protein chains" (p.295, p.312) | Agrees (our list is byte-identical to rlabduke's) |
| B on the residue's **own** mainchain only | "residues having any mainchain atom" — the residue's own atoms | Agrees. Our pipeline's i−1/i+1 neighbour test is **stricter than published** |
| B **< 30** (30.00 excluded) | "> 30" excluded, i.e. 30.00 kept | **Conflicts.** NDFT titles say "B<30". 340 residues have own mainchain max B = 30.00 (General 235, IleVal 47, Gly 33, TransPro 15, PrePro 10, CisPro 0). The earlier edge tests did not cover it (no CisPro residue has B = 30.00; `transpro_edges2.py` only removed B in [29, 30)). **Tested 2026-10-13** (`scripts/b30_boundary.py`, `results/b30_boundary.txt`): TransPro supports **< 30**; see "B = 30.00 boundary test" below |
| TransPro = \|ω\| ≥ 150°; twisted Pro excluded | Taxonomy (p.305), "matching the PDB definition": cis −30° to +30°, trans > +150° or < −150°, "twisted" > 30° from planar; per-residue categories cis Pro / cis nonPro / twisted Pro / twisted nonPro; ω assigned to the following residue | Consistent: the "trans Pro" plot (Fig. 8) holds trans-defined Pro, so twisted Pro are not in it. The paper does not state the exclusion as a filter. The boundary at exactly 150° is not distinguishable in our data (plateau 149–152°) |
| Occupancy 1 on all atoms, incl. Reduce amide H | Not stated | Undocumented (altloc part confirmed by J. Richardson, pers. comm. 2026-10-02) |
| Mainchain B = 0.00 excluded | Not stated | Undocumented |
| Standard residues only | Not stated | Undocumented; possibly a Dangle side effect |
| Class precedence Gly, Pro > pre-Pro > Ile/Val > General | Six plots named (Fig. 8) but no precedence given; precedence is from Read et al. (2011) p.1399 | No conflict |

Other relevant points in the paper:

- **Two Top8000 residue filters.** In 2015 the Top8000 was revised to
  require deposited structure factors, and a new residue filter was built:
  RSCC ≥ 0.7, 2mFo−DFc ρ ≥ 1.1 and B < 40 on all relevant atoms (p.311–312,
  ref. 44 = Hintze et al. 2016, local `Refs/HLR+16.pdf`). It was used for
  the rotamers (983,574 residues) and "future validations". The rama8000
  contours use the earlier B > 30 filter. The rotamer filter writes "< 40",
  so the paper uses both inequality forms.
- **CaBLAM** (p.311) excludes residues where "any of the atoms used in the
  calculation were missing or had B-factor >30". That covers the atoms
  needed for the measurement, across neighbours. The Ramachandran wording
  says "mainchain atom" instead, and our CisPro/TransPro fits support
  own-residue atoms only.
- **Hand checks near thresholds** (Methods, Reference datasets, p.310):
  "we examine several dozen individual examples near the proposed
  thresholds to check for distortions, unconvincing electron density, or
  other circumstances". This describes how thresholds were chosen, not a
  per-residue curation. It is a possible but unconfirmed source of the
  General/PrePro grid residual.
- **Smoothing** (p.310): "iterative kernel-density methods", ref. 11
  (Lovell et al. 2003). No Silk parameters are given; the Makefile in
  javadev is still the only source for them.
- **Size**: "a million quality-filtered residues" (Fig. 8 caption, p.294).
  The NDFT realcounts total 1,504,897 (Gly counted once), so the figure is
  approximate.

### B = 30.00 boundary test (2026-10-13)

`scripts/b30_boundary.py` (output `results/b30_boundary.txt`) reruns the
reconciled rule with own mainchain max B < 30 and with ≤ 30.

TransPro (Silk `-cosine=5 -twopass=8`; 15 rule-set residues at exactly
30.00):

| Boundary | n | Δn | meanΔ×1e4 | support mismatch | rank-identical bins | max \|Δrank\| |
|---|--:|--:|--:|--:|--:|--:|
| B < 30 | 66,363 | −1 | 0.0253 | 0 | 30,851 / 32,400 | 48 |
| B ≤ 30 | 66,378 | +14 | 0.0370 | 0 | 30,718 / 32,400 | 59 |

All classes, Δ vs NDFT realcount under the full reconciled rule:

| Boundary | General | Gly | IleVal | PrePro | TransPro | CisPro | Total |
|---|--:|--:|--:|--:|--:|--:|--:|
| B < 30 | −772 | +76 | +167 | −39 | −1 | 0 | −569 |
| B ≤ 30 | −550 | +108 | +214 | −29 | +14 | 0 | −243 |

TransPro, the class reconciled to within one residue, rejects ≤ 30 on
count and on grid fit. The non-Pro counts are mixed (General and PrePro
closer, Gly and IleVal further), but adding residues moves every count up,
so counts cannot decide the boundary for classes with other residuals.

Non-Pro grid check (`scripts/b30_grid_check.py`, `results/b30_grid_check.txt`;
full reconciled rule; Gly mirrored explicitly because `SilkCmdLine` has no
`-sym`, Gly n is doubled):

| Class | twopass | Boundary | n | Δn | meanΔ×1e4 | supp | rank-identical | max \|Δrank\| |
|---|--:|---|--:|--:|--:|--:|--:|--:|
| General | 6.5 | B < 30 | 1,060,867 | −772 | 0.5576 | 256 | 8,696 | 2,576 |
| General | 6.5 | B ≤ 30 | 1,061,089 | −550 | 0.5659 | 256 | 8,690 | 2,588 |
| PrePro | 8 | B < 30 | 61,529 | −39 | 0.5499 | 7 | 26,795 | 337 |
| PrePro | 8 | B ≤ 30 | 61,539 | −29 | 0.5695 | 7 | 26,761 | 390 |
| Gly | 8 | B < 30 | 233,330 | +152 | 1.2133 | 2,084 | 2,328 | 1,220 |
| Gly | 8 | B ≤ 30 | 233,394 | +216 | 1.3020 | 2,084 | 2,320 | 1,252 |
| IleVal | 6.5 | B < 30 | 195,250 | +167 | 2.9909 | 2,845 | 12,390 | 4,605 |
| IleVal | 6.5 | B ≤ 30 | 195,297 | +214 | 2.9975 | 2,845 | 12,383 | 4,644 |
| IleVal | 8 | B < 30 | 195,250 | +167 | 0.2497 | 98 | 18,435 | 337 |
| IleVal | 8 | B ≤ 30 | 195,297 | +214 | 0.2591 | 98 | 18,468 | 352 |

Every class has a larger mean grid error with ≤ 30, including General and
PrePro, whose counts moved closer. Grid fit, unlike count, does not
improve by adding points in general, so this is a real cross-class result.

Side findings from this run:
- **IleVal smoothing:** `-twopass=8` (Makefile `-mask=10,16`) fits 12×
  better than 6.5 (0.25 vs 2.99). IleVal has no recipe in
  `silk/old/Makefile` (the class is new in rama8000); it was evidently
  smoothed like Gly/Pro/PrePro, not like General.
- **Gly** fits poorly at `-twopass=8` (1.21, 2,084 support-mismatch bins)
  even with mirroring, so its smoothing or symmetry handling still differs.
  Open.
- General was only run at 6.5 here; whether it, too, used 8 is untested.

### Further sources on the residue filter (2026-10-13)

1. **Lab code.** `silk/old/Rota2D.java` (javadev `24e0666`)
   `readDataPoints`: a point is kept only `if(te.Bfactor < cutoff)`; the
   rama targets in `silk/old/Makefile` pass `-blt=30`. `RamaStat.java:106`
   uses `B >= minB && B < maxB`. The lab's tools exclude B = 30.00.
2. **Read et al. (2011)**, *Structure* 19:1395–1412 (`Refs/RAA+11.pdf`),
   Fig. 2 caption: "about 825,000 residues after quality-filtering by
   resolution (<2Å), alternate conformations, and backbone B-factor
   (<30Å2)". General: 582,000 residues from 4,400 nonhomologous files, so
   this is the pre-Top8000 dataset. Also gives the class precedence (p.1399).
3. **Richardson, Keedy & Richardson (2013)**, "The plot thickens"
   (`Refs/RKR13.pdf`), sect. 2.1, describing the Top8000 itself:
   > For purposes of Ramachandran-plot analysis, we omit any residue that
   > has a backbone atom with B-factor ≥30 or an alternate conformation.
   > That leaves 1,573,349 total residues for analysis
   Fig. 3 is labelled "Top8000 B<30 distributions". Other details: PDB
   release as of March 29, 2011 (WHM+18 says March 25); Reduce run on all
   files; chains > 37 residues; 7,957 chains.
   RKR13's total does not reproduce from our φ,ψ table
   (`scripts/rkr13_total.py`, `results/rkr13_total.txt`): B < 30 alone
   gives 1,569,549 and B < 30 + no altloc 1,530,608 (−42,741). The total
   is larger than any φ,ψ-residue count, so it probably includes residues
   without both φ and ψ (chain ends, breaks). It is 68,452 more than the
   NDFT total (1,504,897), consistent with further exclusions (occupancy
   incl. H, zero B, twisted Pro, non-standard) applied when the NDFTs were
   built. Not usable as a check without a table of all residues.

4. **Read et al. (2011) Supplemental Experimental Procedures**
   (`Refs/RAA+11_supp.pdf`, Document S1 via Europe PMC, PMC3195755),
   "Categories for Ramachandran validation":
   > These distributions were made from a MySQL database containing PDB
   > and validation data for over a million residues from 4400
   > nonhomologous chains (at the 70% sequence identity level), chosen for
   > resolution < 2.0Å and an average of resolution and MolProbity score
   > of < 2.0. Individual residues were omitted if they had occupancy < 1.0
   > or any backbone B-factor > 30.
   This is the first published statement of the **occupancy rule**
   (inferred rule 2). It does not say which atoms; our fits require all
   atoms including Reduce amide H. Same text: "> 30" again, and the
   precedence example "a Gly that is also a pre-Pro is judged on the Gly
   distribution". Applies to the pre-Top8000 (4,400-chain) dataset.
   The Williams 2018 supplement (PMC5734394, `PRO-27-293-s001.pdf`) and
   Hintze 2016 (PMC4983197) are not open access in Europe PMC; PMC serves a
   CAPTCHA to scripted downloads, so they need a browser download.

**Conclusion on the boundary:** B < 30 (30.00 excluded). It is stated by
RKR13 ("≥30" omitted) and Read 2011 ("<30Å2"), implemented by `Rota2D`
(`Bfactor < cutoff`), written in the NDFT titles ("B<30"), and supported
by the TransPro count and by every class's grid fit. WHM+18's "> 30" is
the outlier wording. The altloc exclusion is published in Read 2011 and
RKR13, and confirmed by J. Richardson. The occupancy < 1.0 exclusion is
published in the Read 2011 supplement. The "> 30" wording appears in the
Read 2011 supplement as well as WHM+18, against "≥ 30 omitted" (RKR13),
"<30" (Read 2011 Fig. 2, NDFT titles) and `Bfactor < cutoff` (Rota2D);
the data decide for < 30.

Conclusion: the paper documents the chain set and the B-factor test
(own mainchain atoms) but not occupancy, zero-B or residue-type rules,
and its "> 30" disagrees with RKR13, Read 2011, the lab code, the NDFT
titles and our fits, all of which give < 30 (sections above). Questions for
J. Richardson: (1) confirm B = 30.00 was excluded; (2) were the occupancy (including H), zero-B and twisted-Pro
exclusions deliberate, or side effects of Reduce/Dangle; (3) the exact
Dangle and Silk commands used for rama8000.

## Silk parameter scans and support diagnosis (2026-10-13)

`scripts/silk_param_scan.py` (outputs `results/silk_scan_{gly,general}.txt`)
runs Silk over a parameter grid on the reconciled-rule set (B < 30).
`SilkCmdLine` has no `-sym` flag (that is `silk/old/Rota2D`), so Gly points
are mirrored explicitly to (φ, ψ) and (−φ, −ψ). The archived `glycine.ndft`
is titled "Gly (sym)" and is exactly point-symmetric, so this is the right
model.

| Class | Scanned | Best (meanΔ×1e4) | Next best |
|---|---|---|---|
| General | cosine 5; twopass 6, 6.5, 7, 7.5, 8 | twopass 6.5 (0.558) | 1.16 |
| Gly | cosine 4, 5, 6 × twopass 6.5, 7, 8, 9, 10 | cosine 5, twopass 8 (1.213) | 2.40 |

Both minima are sharp and equal the `silk/old/Makefile` values, so the
residual misfit is in residue membership, not smoothing.

`scripts/support_diag.py CLASS NDFT TWOPASS` (outputs
`results/support_diag.txt`, `results/support_diag_nonpro.txt`) splits the
support mismatch into bins nonzero only in our grid and only in the NDFT,
with their densities and distance from the common support, and lists
residues whose own bin lies in our-only support.

| Class | ours-only bins | NDFT-only bins | ours-only density (counts) | distance from common support |
|---|--:|--:|---|---|
| TransPro (twisted excluded) | 0 | 0 | – | – |
| PrePro | 5 | 2 | 1 | 1 bin |
| IleVal | 89 | 9 | 1–2 | median 1, max 3 |
| General | 239 | 17 | median 3, max 6 | median 2.2, max 7.6 |
| Gly | 2,084 | 0 | median 2, max 12 | median 3, max 17 |

Reading:
- **Gly** (+76): all mismatch is ours, a low-density fringe outside the
  lab's support. Our set has extra points in sparse regions that the lab's
  set lacks. Only 2 residues (1w0p A729, 3k26 A190) sit inside our-only
  bins; the others lie inside the common support and widen it through the
  density-dependent second pass. Next: rank Gly by local density and
  characterise the sparse ones.
- **General** (−772): fewer residues overall but more extra fringe, so
  membership differs in both directions.
- **IleVal** (+167): thin fringe, mostly ours; consistent with a small
  excess of edge points.
- **PrePro** (−39): support almost exact; the misfit is interior density.

## Preceding-peptide (non-trans) exclusion for non-Pro classes (2026-10-13)

`scripts/gly_sparse.py` (`results/gly_sparse.txt`) profiled the sparsest
Gly residues and tested candidate exclusions. Only one helps: dropping Gly
whose **preceding peptide is non-trans** (|ω| < 150°, cis or twisted
non-Pro). Neighbour-ω, neighbour-B and neighbour-occupancy rules all make
the fit much worse. `scripts/nontrans_nonpro.py`
(`results/nontrans_nonpro.txt`) splits cis from twisted and extends the
test to all non-Pro classes (Makefile smoothing; reconciled rule, B < 30):

| Class | Variant | n | Δn | meanΔ×1e4 | supp (ours/NDFT) | rank-identical |
|---|---|--:|--:|--:|---|--:|
| Gly | reconciled rule | 116,665 | +76 | 1.2133 | 2,084 / 0 | 2,328 |
| Gly | drop cis | 116,616 | +27 | 0.7399 | 2,084 / 0 | 2,220 |
| Gly | drop twisted | 116,645 | +56 | 0.8775 | 2 / 0 | 10,194 |
| Gly | **drop non-trans** | 116,596 | **+7** | **0.4750** | **2 / 0** | **14,748** |
| PrePro | reconciled rule | 61,529 | −39 | 0.5499 | 5 / 2 | 26,795 |
| PrePro | **drop non-trans** | 61,500 | −68 | **0.4460** | 5 / 1 | 26,933 |
| IleVal | reconciled rule | 195,250 | +167 | 0.2497 | 89 / 9 | 18,435 |
| IleVal | **drop non-trans** | 195,192 | **+109** | **0.1474** | **0 / 9** | 20,618 |
| General | reconciled rule | 1,060,867 | −772 | 0.5576 | 239 / 17 | 8,696 |
| General | **drop non-trans** | 1,060,356 | −1,283 | **0.3630** | **10 / 12** | **15,899** |

(Full split by cis/twisted for every class in the results file.) Dropping
non-trans improves the grid in **every** class; both cis and twisted
contribute, and the combination is best. With the Pro rule (TransPro =
|ω| ≥ 150°, CisPro = |ω| < 30°) this generalises to: **a residue's
preceding peptide must be trans, except for the CisPro class.** Counts
after the rule: Gly +7, IleVal +109, PrePro −68, General −1,283 — the
grids now say General and PrePro lack residues that the lab included,
i.e. some other exclusion of ours is too broad (candidates: amide-H
occupancy scope, peptide-link test, zero-B, non-standard neighbours).

Literature: no paper states this exclusion. WHM+18 (p.305) treats
cis-nonPro and twisted peptides as near-certain errors ("twisted peptides
should virtually always be considered modeling errors"); LDA+03 notes the
Top500 Pro plot "includes prolines preceded by either trans or cis
peptide bonds", i.e. before the cis/trans Pro split.

## Further literature review (2026-10-13)

Sources in `~/Desktop/UnFiled261002_RamaRichardson/Refs/`:
- **LDA+03** (Lovell 2003, Top500): residues kept only if no main-chain
  atom has B ≥ 30 (strict); files screened for "large numbers of B-factors
  ≤1, which is an indication of the use of U rather than B, or of
  unrefined B-factors" — precedent for the zero-B exclusion (2fcr, 2yyh);
  dihedrals from DANG.
- **HLR+16** (Hintze 2016) and supplement: Top8000 chain selection —
  PDB as of March 29, 2011; chains ≥ 38 residues at < 2.0 Å; Reduce with
  N/Q/H flips; chain MolProbity score < 2.0; ≤ 5% bond-length, bond-angle
  (> 4σ) and Cβ (> 0.25 Å) outliers; best average of resolution and
  MolProbity score per homology cluster, ties broken by PDB ID + chain ID.
  Rotamer residue filter (B < 40, RSCC, map value, clashes, no altlocs) is
  not the rama8000 filter. Top500 rotamer filter: "no clashes, no
  alternates, and no B factors ≥ 40".
- **WRR22** (Top2018) and supplement: later dataset, B < 40 + RSCC > 0.7 +
  map > 1.2σ + no geometry outliers / clashes / altlocs; non-standard
  residues kept if part of the chain. Not the rama8000 regime.
- **RVW+25** (cis-nonPro): Top8000 50% SF list with B ≤ 40 backbone + Cβ
  and density filters; does not describe the rama contours.
- **WHM+18 supplement**: H-parameter methods only; Top8000 used with B
  cutoffs on atom pairs (< 10, raised to < 30), not residue filtering.
- Not reviewed in depth (not Richardson-lab datasets): SD13, TWS+10.
