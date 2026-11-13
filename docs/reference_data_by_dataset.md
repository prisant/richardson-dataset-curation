# Reference data available for each dataset

*File: reference_data_by_dataset.md · Git: Michael G. Prisant 4c7b2f1 · Txd: 11/12/26 with local modifications!*

Pipeline diagram: `pipeline_block_diagram.pdf` (source `pipeline_block_diagram.tex`).

"Lab" means the Richardson Lab's own files. "Ours" means what our pipeline
uses now. "Not located" means a source document says the files exist but I
have not found them on disk.

| Dataset | Chain list | Lab's coordinate files (as the lab used them) | Our coordinate files | Lab's filter output | Lab's filter rules | Lab's reference results |
|---|---|---|---|---|---|---|
| top100 | Reconstructed from the lab's `top100H` files (100 entries, 108 chains) | 100 FH files, Reduce 2.13.2: `~/Desktop/MyProgs/PolarWork/top100H/` | RCSB downloads, Reduce 4.16 | None | Word et al. 1999b: atoms with B ≥ 40 omitted (our filter: mainchain B ≤ 40) | None located |
| top500 | Reconstructed from the lab's kinemage-server list (494 entries, 523 chains) | 500 FH files, Reduce 2.13.2 (per `top500/README.md`): not located | RCSB downloads, Reduce 4.16 | None | Lovell et al. 2003: "residues with B < 30" | φ/ψ archive, 2006 (`top500_issues.md`): not located |
| top8000 | `Top8000-best_hom70_pdb_chain_list.csv` (7,737 entries, 7,957 chains) | 2011 chain files, Reduce 3.15.091106: `~/Desktop/Data/top8000_chains_70_2011/` (sha256 4a6bb1b9…) | RCSB downloads, Reduce 4.16 | The filtered residue lists themselves (`Top8000_Rama_2013.zip`); no masks | `get_srcdata.py` in `Top8000_Rama_2013.zip` | Filtered residue lists (`Top8000_Rama_2013.zip`); `.ndft` tables |
| top2018_full | Top2018 all-atom chain list v2.0 (Zenodo 10.5281/zenodo.5773255) | Pruned single-chain files (`*_pruned_all.pdb`, kept residues only, Reduce 3.3 hinge flips); whole-entry files the lab started from: not held | RCSB downloads, Reduce 3.3.160422 (4.16 for 4oan, 4uzg) | Pruned files: our `.mask` files are derived from them | Williams et al. 2021; Top2018 database tables `allatom_filters`, `mc_filters` | Top2018 database: `~/Desktop/Data/williams_top2018_sqlite/top2018_build.db` |
| top2018_mc | Top2018 mainchain chain list v2.01 (Zenodo 10.5281/zenodo.5777651) | Pruned single-chain files (`*_pruned_mc.pdb`); whole-entry files: not held | RCSB downloads, Reduce 4.16 | Pruned files: our `.mask` files are derived from them | Williams et al. 2021; Top2018 database table `mc_filters` | Top2018 database (as above) |

## Masks for top2018

The lab did not supply mask files; it supplied the pruned chain files.
`build_ersatz_pdbs.py` derives our `.mask` file for each chain from the
residues present in the lab's pruned file, and keeps the file's `USER INC`
records as the fragment boundaries. So for top2018 the lab's filter is
already applied; we do not apply filter rules ourselves. The rules (Williams
et al. 2021, and the database's per-residue filter flags) serve only to check
the result.
