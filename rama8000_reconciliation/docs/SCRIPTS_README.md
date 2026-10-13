# rama8000 reconciliation — script bundle

Companion to `rama8000_reconciliation_summary.md` (the method and results).
These are working research scripts, written around Michael Prisant's local
file layout. They are shared so the workflow can be checked or reproduced;
they are not a packaged tool.

## Requirements

- Python 3.10+ with `numpy`, `pandas`, `scipy`
- Java (OpenJDK 21 used) to build and run Silk
- Reduce 4.16 (`reduce` on PATH) and its het dictionary
- `pandoc` is not needed

## Local layout the scripts expect

Edit `ROOT` in `scripts/recount3.py` to point at the data tree. Per PDB entry:

```
ROOT/<id[0:2]>/<id>/<id>.pdb          # RCSB download
ROOT/<id[0:2]>/<id>/<id>FH.pdb        # Reduce output (hydrogens added)
ROOT/<id[0:2]>/<id>/<id>_ersatz.pdb   # Reduce output, hydrogens and some headers stripped
ROOT/Top8000-best_hom70_pdb_chain_list.csv
```

`pipeline/download_pdbs.py`, `pipeline/run_reduce.py` and
`pipeline/build_ersatz_only.py` build that tree. Parsing the `FH.pdb` files
directly is equivalent (see the summary, Workflow step 3).

Scripts that call Silk run from the bundle root and expect the compiled
classes in `silk_build/classes` (summary, Workflow step 9).
`scripts/reduce_old_amideH.py` expects an older Reduce at `~/bin/reduce_old`;
it is only needed to repeat the Reduce 3.3 vs 4.16 check.

## Run order (main path)

| Step | Script | Output | Expected result |
|---|---|---|---|
| 1 | `scripts/residue_records.py data/residue_records.tsv` | one row per residue with φ, ψ, ω, class, B and occupancy | ~1.85 M rows |
| 2 | `scripts/amideH_records.py data/amideH_records.tsv` | amide-H presence, B, occupancy from `FH.pdb` | one row per residue |
| 3 | `scripts/amideH_counts.py` | class counts vs NDFT under the rule variants | last-but-one row: General −772, Gly +76, IleVal +167, PrePro −39, TransPro −1, CisPro 0 |
| 4 | `scripts/cis_records_ext.py data/pro_records_ext.jsonl 999` | every Pro with extended attributes | 85,481 records |
| 5 | `scripts/class_rule_check.py data/pro_records_ext.jsonl cispro data/cispro.ndft 8` | CisPro counts + grid agreement | occupancy rule: n = 3,654, 32,090 rank-identical bins |
| 6 | `scripts/bzero_variants.py` | TransPro zero-B variants + grid | 66,363 with all-zero mainchain B dropped |
| 7 | `scripts/nonpro_grid_check.py General data/general.ndft 6.5` (and `PrePro … 8`) | non-Pro grid agreement | see summary, Results |

`data/*.ndft` are the six archived tables from javadev `24e0666`
(`chiropraxis/resource/chiropraxis/rotarama/rama8000/`).

## Exploratory scripts (history of the investigation)

`recount3.py` (first recount, also holds shared helpers), `pro_records.py`,
`silk_fit.py` / `silk_fit_refine.py` (the early parameter search — superseded:
parameters are now fixed at the Makefile values), `cis_group_test.py`,
`cis_single_test.py`, `cis_missing_islands.py`, `cis_edge_groups.py`,
`cis_final_check.py`, `transpro_islands.py`, `transpro_edges.py`,
`transpro_edges2.py`, `transpro_loo_sparse.py` (≈10 min),
`transpro_bzero.py`, `general_counts.py`, `reduce_old_amideH.py`.
Several import helpers from one another, so keep them together in `scripts/`.
