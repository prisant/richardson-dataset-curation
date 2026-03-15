# Supplementary PDB Files

These three PDB files are no longer available from RCSB in PDB format
(only mmCIF is offered). They are included here to ensure full
reproducibility of the pipeline, which requires PDB-format input for
Reduce.

| File | Chains | Atoms | Notes |
|------|--------|-------|-------|
| 4v4m.pdb | w (target) | 86,696 | Large ribosome structure |
| 4yza.pdb | A (target) | 3,369 | |
| 4ztt.pdb | C (target) | 10,729 | |

These files were originally downloaded from RCSB (https://files.rcsb.org)
in PDB format before RCSB discontinued PDB format for these entries.

To use: copy these files into the appropriate entry directories before
running the pipeline:

```bash
for pdb in 4v4m 4yza 4ztt; do
  prefix=${pdb:0:2}
  cp supplementary_pdbs/${pdb}.pdb top2018_pdbs_full_filtered_hom70/${prefix}/${pdb}/${pdb}.pdb
done
```
