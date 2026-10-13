#!/usr/bin/env python3
"""Independent arithmetic audit of top8000 Ramachandran 6-category counts.

Streams top8000_measures.jsonl(.gz), tallies rama_category values, and diffs
against the numbers reported in top8000_analysis.md (regenerated 2026-04-27,
commit a02b65e). Deliberately independent of pipeline/analyze_jsonl.py: uses
only json + Counter so any shared bug would not hide behind a shared code path.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

EXPECTED = {
    "General":  1_056_115,
    "IleVal":     200_462,
    "Gly":        112_610,
    "TransPro":    63_261,
    "PrePro":      50_230,
    "CisPro":       3_524,
}
EXPECTED_PRESENT = 1_486_202
EXPECTED_NULL    =    82_359
EXPECTED_TOTAL   = 1_568_561


def open_jsonl(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "jsonl",
        type=Path,
        nargs="?",
        default=Path(__file__).resolve().parent.parent
                / "top8000" / "top8000_measures.jsonl.gz",
    )
    args = ap.parse_args()

    counts: Counter[str] = Counter()
    null_count = 0
    total = 0
    meta_seen = 0

    with open_jsonl(args.jsonl) as fh:
        for line in fh:
            rec = json.loads(line)
            if rec.get("_meta"):
                meta_seen += 1
                continue
            total += 1
            cat = rec.get("rama_category")
            if cat is None:
                null_count += 1
            else:
                counts[cat] += 1

    present = sum(counts.values())

    print(f"File: {args.jsonl}")
    print(f"_meta records skipped: {meta_seen}")
    print(f"Data records:          {total:>10,}")
    print(f"rama_category present: {present:>10,}")
    print(f"rama_category null:    {null_count:>10,}")
    print()
    print(f"{'Category':<10} {'Observed':>10} {'Expected':>10} {'Delta':>7}")
    ok = True
    for cat, exp in EXPECTED.items():
        obs = counts.get(cat, 0)
        delta = obs - exp
        flag = "" if delta == 0 else "  MISMATCH"
        if delta != 0:
            ok = False
        print(f"{cat:<10} {obs:>10,} {exp:>10,} {delta:>+7,d}{flag}")

    extras = set(counts) - set(EXPECTED)
    if extras:
        ok = False
        print("\nUnexpected categories in JSONL (not in the 6-bucket schema):")
        for cat in sorted(extras):
            print(f"  {cat}: {counts[cat]:,}")

    print()
    checks = [
        ("sum(6 buckets) == present",         present == sum(EXPECTED.values())),
        ("present matches doc",               present == EXPECTED_PRESENT),
        ("null matches doc",                  null_count == EXPECTED_NULL),
        ("present + null == total",           present + null_count == total),
        ("total matches doc",                 total == EXPECTED_TOTAL),
    ]
    for name, passed in checks:
        print(f"  [{'OK' if passed else 'FAIL'}] {name}")
        if not passed:
            ok = False

    print("\nAUDIT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
