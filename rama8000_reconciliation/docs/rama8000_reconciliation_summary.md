# rama8000 count reconciliation: top8000 vs the Richardson Lab NDFT tables

Michael G. Prisant · 2026-10-04 · Prepared for Erik Novak

The 19,273-observation shortfall Erik found comes from residue selection, not smoothing. Under five inferred selection rules, a recount from coordinates matches every rama8000 class to within 0.1%, and CisPro exactly.

## Results

All counts are residues with both φ and ψ, in the 7,957 chains of Top8000-best\_hom70. The published column is Zenodo v2.0.0 after Erik's Ile/Val-before-Pro correction (9,281 residues moved to PrePro). Gly's NDFT realcount (233,178) is halved because Silk's `-sym` mirrors each point.

| Class | NDFT realcount | Published (corrected) | Inferred rules | Δ inferred − NDFT | Grid check (Silk, Makefile params) |
| --- | --: | --: | --: | --: | --- |
| General | 1,061,639 | 1,055,696 | 1,060,867 | −772 | mean \|Δfraction\| 0.56×10⁻⁴; 256 support-mismatch bins |
| Gly | 116,589 | 112,542 | 116,665 | +76 | not yet run |
| IleVal | 195,083 | 191,138 | 195,250 | +167 | not yet run |
| PrePro | 61,568 | 59,498 | 61,529 | −39 | 0.55×10⁻⁴; 7 bins |
| TransPro | 66,364 | 63,227 | 66,363 | −1 | 0.025×10⁻⁴; 0 bins; 30,851 / 32,400 bins rank-identical, max \|Δrank\| 48 |
| CisPro | 3,654 | 3,523 | 3,654 | 0 | 0.028×10⁻⁴; 1 bin; 32,090 / 32,400 bins rank-identical, max \|Δrank\| 3 |
| **Total** | **1,504,897** | **1,485,624** | **1,504,328** | **−569** |  |

The two Pro classes reproduce essentially exactly, grid included. The four non-Pro classes match in count to within 0.1%, but their grids fit roughly 20× worse than TransPro's; see Open questions.

## Inferred rama8000 selection rules

A residue enters its class only if it passes all five rules. Rules 1 and 2 account for almost all of the gap; 3–5 fix the last few hundred.

| # | Rule (rama8000, inferred) | Published top8000 pipeline | Evidence |
| --- | --- | --- | --- |
| 1 | B < 30 on the residue's **own** mainchain N, CA, C, O only | B ≤ 30 required on the mainchain of i−1, i and i+1 | Literal reading of Williams 2018 ("residues having any mainchain atom with a B-factor > 30" excluded); NDFT titles say "B<30". Adding the 54 cis-Pro residues this admits cut CisPro support mismatch from 73 to 2 bins |
| 2 | **Every atom** of the residue at occupancy 1, **including Reduce hydrogens** | Not tested; altloc residues kept (one conformation chosen) | Jane Richardson confirmed altloc residues were dropped (personal communication to MP, 2026-10-02). Occupancy covers all altlocs, plus 9 cis-Pro residues with partial or zero occupancy and no altloc label (CisPro then exactly 3,654). Reduce gives the amide H of a residue next to an altloc residue partial occupancy; dropping those (\~10k residues) took General/Gly/IleVal/PrePro from +0.6–0.8% to within 0.1%. Pro has no amide H, so the Pro classes were unaffected |
| 3 | TransPro means \|ω\| ≥ 150°; Pro with 30° ≤ \|ω\| < 150° (twisted) is in no class | TransPro = \|ω\| ≥ 30° | TransPro support mismatch 759 → 0 bins; results identical for cutoffs 149–152°, worse at 145° and 155° |
| 4 | Residues whose mainchain B-factors are all 0.00 excluded | Kept | 14 Pro residues (2fcr A ×10, 2yyh C ×4); TransPro +13 → −1. "Any mainchain atom = 0" gives −3 with an indistinguishable grid; "any atom = 0" is ruled out (breaks CisPro) |
| 5 | Standard amino acids only (no MSE or other modified residues) | Same | Adding MSE makes General worse (+2,854) |

Class precedence is Read et al. 2011's: Gly, Pro > pre-Pro > Ile/Val > General. CisPro is |ω| < 30°, with ω the peptide bond preceding the Pro.

## Inputs and provenance

The published JSONL cannot be used for this: it carries no B-factors or occupancies. Start from coordinates.

| Input | Where | Version / check |
| --- | --- | --- |
| Chain list (7,957 chains, 7,737 entries) | [Top8000-best\_hom70\_pdb\_chain\_list.csv](https://raw.githubusercontent.com/rlabduke/reference_data/master/Top8000/Top8000-best_hom70_pdb_chain_list.csv) | Byte-identical to our copy (checked 2026-10-04) |
| Coordinates | `https://files.rcsb.org/download/<ID>.pdb` | RCSB files: 4,613 downloaded 2022-05-29 (shared with top2018), 3,124 on 2026-03-15; 4v4m, 4yza, 4ztt no longer in PDB format at RCSB |
| Hydrogens | Reduce 4.16.250520: `reduce -db reduce_wwPDB_het_dict.txt -allalt -build -norotmet <ID>.pdb` | Reduce 3.3 (2016-04-22) with `-FLIP` gives identical amide-H occupancies on all 7,737 entries |
| NDFT tables | [rlabduke/javadev @ 24e0666](https://github.com/rlabduke/javadev/tree/24e0666307a376fcb1b7017a0706d9069fadd93a/chiropraxis/resource/chiropraxis/rotarama/rama8000): general, glycine, ileval, prepro, transpro, cispro `.ndft` | Commit `24e0666307a376fcb1b7017a0706d9069fadd93a` |
| Silk (lab's KDE tool) | Same commit, `silk/src` + `driftwood/src`; source tarball `codeload.github.com/rlabduke/javadev/tar.gz/24e0666…` | Tarball sha256 `5d6dfc42d483e7c02aa765f3410619ff9328d408cce30520b576f044f50ebac5`; built with OpenJDK 21 |
| Smoothing recipe | `silk/old/Makefile` (old top500 recipe) in the same commit | `-mask=10,16` (Pro, Gly, PrePro); `-mask=10,13` (General) |
| Published dataset (for the baseline only) | [Zenodo 10.5281/zenodo.19051741](https://doi.org/10.5281/zenodo.19051741), v2.0.0 | `top8000_measures.jsonl.gz` md5 `7592037eff515404e43a0b21d8bdf0d9` |

Method reference: Williams et al. 2018, Protein Science 27:293–315, [doi:10.1002/pro.3330](https://doi.org/10.1002/pro.3330).

## Workflow to reproduce

Following these steps should give the "Inferred rules" counts above to within a few residues. Our own run parsed Reduce output with hydrogens stripped for the heavy atoms, plus the Reduce file itself for amide H; parsing the Reduce file directly is equivalent.

1. **Download.** For each PDB ID in the chain list, fetch `<ID>.pdb` from RCSB. Use model 1 only (stop at the first `ENDMDL`). Keep only chains named in the list.
2. **Add hydrogens.** Run Reduce with the flags in the provenance table. Use its output for every later step; Reduce does not move backbone atoms, so φ/ψ are unchanged.
3. **Parse.** Read `ATOM` and `HETATM` records. Group atoms into residues by (chain, resSeq, iCode), in file order. For each atom name keep the first coordinates seen (first altloc), **every** B value, and **every** occupancy across altloc copies. Include H atoms.
4. **Eligibility.** Residue i needs N, CA, C; residues i−1 and i+1 adjacent in file order in the same chain; C(i−1)–N(i) < 2.0 Å and C(i)–N(i+1) < 2.0 Å; CA on i−1.
5. **Angles.** φ = C(i−1)–N–CA–C; ψ = N–CA–C–N(i+1); ω = CA(i−1)–C(i−1)–N(i)–CA(i). Degrees, (−180, 180\].
6. **Class.** Residue i must be one of the 20 standard amino acids. GLY → Gly; PRO → CisPro if |ω| < 30, else TransPro; next residue PRO → PrePro; ILE/VAL → IleVal; otherwise General.
7. **Filter.** Keep i only if all hold:
   1. max B over all copies of its own N, CA, C, O < 30;
   2. that max B > 0;
   3. min occupancy over **all** its atoms, hydrogens included, = 1.00;
   4. for TransPro, |ω| ≥ 150.
8. **Count** per class and compare with the realcounts (Gly's NDFT realcount is 2× the residue count).
9. **Build Silk.** From javadev `24e0666`: `javac -d classes -encoding ISO-8859-1 -sourcepath silk/src:driftwood/src $(find silk/src -name '*.java')`, then copy `silk/resource/*` into `classes/`.
10. **Smooth.** Write one line per residue, `label phi psi`, space-separated, then run:
    - `java -cp classes silk.SilkCmdLine -ndim=2 -coords=2,3 -label=1 -bounds=-180,180,-180,180 -wrap -gridsize=2 -cosine=5 -twopass=8 -fraction -ndft < points.txt > out.ndft`
    - Gly: add `-sym`. General: `-twopass=6.5`. (Makefile `-mask=a,b` = Silk `-cosine=a/2 -twopass=b/2`, because Silk's cosine mask radius is 2 × halfwidth.)
11. **Read NDFT** (big-endian): NUL-terminated title; int32 nDim; float32 min\[nDim\]; float32 max\[nDim\]; int32 bins\[nDim\]; byte wrap\[nDim\]; int32 realcount; then float32 values. Reshape to 180 × 180: axis 0 = φ, axis 1 = ψ, bin centre = −179 + 2i degrees.
12. **Compare** your grid with the archived one (metrics in the next section).

## Testing a rule against the grid

The NDFT holds no residue identities, but its grid is a fingerprint of the input set. Smoothing settings stay fixed at the Makefile values; only the residue filter varies.

- **Why realcount is a true count.** Silk's `-fraction` value at a bin is the fraction of input points with lower density, so every grid value is an exact multiple of 1/realcount (the smallest nonzero CisPro value is 1/3,654).
- **Metrics.** Mean |Δfraction| over all 32,400 bins; support mismatch (bins zero in one grid, nonzero in the other); and rank-identical bins (round(value × realcount) equal).
- **Group test.** Add or remove a whole group of residues defined by one attribute (B in a range, altloc, occupancy, ω range), rerun Silk and keep the change only if the grid improves. Single-residue tests rarely work: inside a dense region every point has the same effect.
- **Islands.** In sparse regions one isolated point paints many bins, so clusters of support mismatch point straight at outliers. That is how the twisted cis/trans Pro were found (for example 3pve A1574, ω = +65°; 2asd A160, ω = −117°).

The 9 cis-Pro residues removed by the occupancy rule are 1gtv A37, 1k3i A52, A163, A350, 2gw3 B49, B84, 2gw4 D84, 2qod A694 and 3o0p A221. The last has partial occupancy only on CB and CG, which shows the test covers every atom, not just the mainchain.

## Open questions and caveats

- **Non-Pro grids.** General and PrePro match in count but their grids fit about 20× worse than TransPro's. Candidate causes: the real rama8000 build used different smoothing than the old top500 Makefile; coordinates changed since about 2011; or a filter detail not yet found. Reduce version is ruled out (3.3 = 4.16).
- **Gly and IleVal** have not been grid-tested yet.
- **Zero-B rule form** (all mainchain atoms at 0 versus any one) cannot be told apart from the data.
- **Identity.** A match rests on count plus grid. Two residues swapped at nearly the same φ/ψ would be invisible.
- **Smoothing recipe** is inferred from the old top500 Makefile; the rama8000 build script itself has not been found.

Questions for the Richardson Lab:

1. The rama8000 build commands: the Dangle invocation and the Silk options for each class.
2. Was the altloc drop implemented as "every atom, hydrogens included, at occupancy 1"?
3. Was Pro with 30° < |ω| < 150° excluded?
4. Were residues with zero B-factors excluded?
