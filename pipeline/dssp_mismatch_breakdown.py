"""Tabulate (JSONL dssp, .dssp dssp) pairs among validate_dssp mismatches.

Scratch diagnostic; reuses pipeline/validate_dssp.py helpers so the
matching logic is identical to the committed validator.
Usage: python mismatch_breakdown.py MEASURES.jsonl[.gz] SRC_DIR
"""
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, "/home/prisant/Desktop/Projects/richardson-dataset-curation/pipeline")
import validate_dssp as v  # noqa: E402

jsonl_path, src_dir = Path(sys.argv[1]), Path(sys.argv[2])
index = v._scan_dssp_files(src_dir)
cache = {}
pairs = Counter()
per_pdb_allnull = Counter()
per_pdb_total = Counter()
with v._open_jsonl(jsonl_path) as fh:
    for raw in fh:
        rec = json.loads(raw)
        if rec.get("_meta") or rec.get("model", 1) != 1:
            continue
        pdbid = v._pdbid_from_file_field(rec["file"])
        if pdbid not in cache:
            p = index.get(pdbid)
            cache[pdbid] = v._load_assignments(p) if p else None
        a = cache[pdbid]
        if a is None:
            continue
        key = (rec["chain"], int(rec["resnum"]), (rec.get("ins") or "").strip())
        if key not in a:
            continue
        j = rec.get("dssp") or None
        r = a[key]
        per_pdb_total[pdbid] += 1
        if j is None and r is not None:
            per_pdb_allnull[pdbid] += 1
        if j != r:
            pairs[(j, r)] += 1

total = sum(pairs.values())
enc = pairs[("C", None)]
null_vs_code = sum(n for (j, r), n in pairs.items() if j is None and r is not None)
other = total - enc - null_vs_code
print(f"mismatches total           {total:>9,}")
print(f"  'C' vs blank (encoding)  {enc:>9,}")
print(f"  null vs real code        {null_vs_code:>9,}")
print(f"  other code disagreements {other:>9,}")
print("\ntop 15 (JSONL, .dssp) pairs:")
for (j, r), n in pairs.most_common(15):
    print(f"  {j!r:>6} vs {r!r:<6} {n:>9,}")
