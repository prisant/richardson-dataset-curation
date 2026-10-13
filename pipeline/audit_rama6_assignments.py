#!/usr/bin/env python3
"""Independent assignment-correctness audit of top8000 rama_category.

Companion to audit_rama6.py, which only checks bucket totals.
This script verifies that each residue's stored rama_category value is
consistent with the pydangle decision tree in
`pydangle_biopython/labels.py::label_rama_category`:

    1. GLY                                    -> "Gly"
    2. ILE or VAL                             -> "IleVal"   (priority over PrePro)
    3. PRO with |omega| < 30                  -> "CisPro"
       PRO with |omega| >= 30 or omega unknown -> "TransPro" (default)
    4. Next residue is PRO (non Gly/Ile/Val/Pro) -> "PrePro"
    5. Otherwise                              -> "General"

Reimplemented from scratch here so any bug in pydangle would not hide
behind a shared code path. Uses only json + stdlib.

Two check tiers:
  Backward (strong) - the label constrains what resname must be.
  Forward  (partial) - PrePro requires the *physical* next residue to be
    PRO. The pipeline B-filter can drop residues from the JSONL, so
    forward checks are only made when JSONL adjacency is unambiguous
    (same file/model/chain, contiguous resnum, matching ins code).
    Ambiguous cases are counted as "unverifiable", not as mismatches.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

GLY = {"GLY"}
ILEVAL = {"ILE", "VAL"}
PRO = {"PRO"}
EXCLUDED_FROM_PREPRO = GLY | ILEVAL | PRO
ALL_SIX = {"Gly", "IleVal", "CisPro", "TransPro", "PrePro", "General"}


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
    ap.add_argument(
        "--sample-mismatches",
        type=int,
        default=10,
        help="Print up to N example mismatches per rule (default: 10)",
    )
    args = ap.parse_args()

    # First pass: load all records, keyed by (file, model, chain), preserving order.
    # Needed so we can inspect the "next residue" for forward PrePro checks.
    chains: dict[tuple[str, int, str], list[dict[str, Any]]] = defaultdict(list)
    total = 0
    with open_jsonl(args.jsonl) as fh:
        for line in fh:
            rec = json.loads(line)
            if rec.get("_meta"):
                continue
            total += 1
            key = (rec["file"], rec["model"], rec["chain"])
            chains[key].append(rec)

    # Tally structures
    bucket_counts: Counter[str] = Counter()
    resname_by_bucket: dict[str, Counter[str]] = defaultdict(Counter)

    # Backward-check mismatches: label => list of example records
    backward_mismatches: dict[str, list[dict[str, Any]]] = defaultdict(list)
    backward_mismatch_counts: Counter[str] = Counter()

    # peptide_bond / rama_category consistency
    pep_mismatches: list[dict[str, Any]] = []
    pep_mismatch_count = 0

    # Forward PrePro checks
    prepro_verified_correct = 0
    prepro_verified_wrong: list[dict[str, Any]] = []
    prepro_verified_wrong_count = 0
    prepro_unverifiable = 0

    # Also check: labeled PrePro but JSONL-next is clearly NOT Pro
    prepro_labeled_but_next_not_pro: list[dict[str, Any]] = []
    prepro_labeled_but_next_not_pro_count = 0

    # Also check: NOT labeled PrePro but JSONL-next is clearly Pro
    should_be_prepro: list[dict[str, Any]] = []
    should_be_prepro_count = 0

    unknown_labels: Counter[str] = Counter()
    null_count = 0
    nonprotein_null_count = 0

    def record_backward(rule: str, rec: dict[str, Any]) -> None:
        nonlocal backward_mismatch_counts
        backward_mismatch_counts[rule] += 1
        if len(backward_mismatches[rule]) < args.sample_mismatches:
            backward_mismatches[rule].append(rec)

    for key, records in chains.items():
        for i, rec in enumerate(records):
            cat = rec.get("rama_category")
            resname = rec.get("resname", "")
            if cat is None:
                null_count += 1
                # Sanity: null rama_category typically comes from missing phi/psi.
                # We do not require a specific resname for nulls.
                continue
            if cat not in ALL_SIX:
                unknown_labels[cat] += 1
                continue

            bucket_counts[cat] += 1
            resname_by_bucket[cat][resname] += 1

            # ---- Backward checks (label -> resname constraint) ----
            if cat == "Gly" and resname not in GLY:
                record_backward("Gly_but_not_GLY", rec)
            if cat == "IleVal" and resname not in ILEVAL:
                record_backward("IleVal_but_not_ILE_or_VAL", rec)
            if cat in ("CisPro", "TransPro") and resname not in PRO:
                record_backward(f"{cat}_but_not_PRO", rec)
            if cat == "PrePro" and resname in EXCLUDED_FROM_PREPRO:
                record_backward("PrePro_but_resname_excluded", rec)
            if cat == "General" and resname in EXCLUDED_FROM_PREPRO:
                record_backward("General_but_resname_excluded", rec)

            # ---- Forward check for Gly bidirection: GLY must be Gly (unless null) ----
            if resname == "GLY" and cat != "Gly":
                record_backward("GLY_resname_not_labeled_Gly", rec)

            # ---- peptide_bond consistency for CisPro/TransPro ----
            pb = rec.get("peptide_bond")
            if cat == "CisPro" and pb != "cis":
                pep_mismatch_count += 1
                if len(pep_mismatches) < args.sample_mismatches:
                    pep_mismatches.append(rec)
            if cat == "TransPro" and pb not in ("trans", None):
                pep_mismatch_count += 1
                if len(pep_mismatches) < args.sample_mismatches:
                    pep_mismatches.append(rec)

            # ---- Forward PrePro check ----
            # Only meaningful when resname is a candidate (not GLY/ILE/VAL/PRO).
            if resname not in EXCLUDED_FROM_PREPRO:
                # Determine "true next residue is PRO?" from JSONL adjacency.
                next_rec = records[i + 1] if i + 1 < len(records) else None
                next_is_verified_pro = False
                next_is_verified_non_pro = False
                if next_rec is not None:
                    contiguous = (
                        next_rec.get("resnum") == rec.get("resnum", 0) + 1
                        and next_rec.get("ins") == rec.get("ins")
                    )
                    # Even at contiguous resnums, a fragment-boundary
                    # signature (psi(i) null AND phi(i+1)/omega(i+1) null)
                    # means the pre-filter data had a physical gap, so
                    # pydangle's residue_list did not treat them as
                    # neighbours. Treat as unverifiable.
                    frag_boundary = (
                        rec.get("psi") is None
                        and next_rec.get("phi") is None
                        and next_rec.get("omega") is None
                    )
                    if contiguous and not frag_boundary:
                        if next_rec.get("resname") == "PRO":
                            next_is_verified_pro = True
                        else:
                            next_is_verified_non_pro = True

                if cat == "PrePro":
                    if next_is_verified_pro:
                        prepro_verified_correct += 1
                    elif next_is_verified_non_pro:
                        # Definitely wrong: visible-and-contiguous next isn't PRO
                        prepro_labeled_but_next_not_pro_count += 1
                        if len(prepro_labeled_but_next_not_pro) < args.sample_mismatches:
                            prepro_labeled_but_next_not_pro.append(rec)
                    else:
                        prepro_unverifiable += 1
                else:
                    # cat is General (already filtered non-excluded resnames)
                    if next_is_verified_pro:
                        # Should have been PrePro
                        should_be_prepro_count += 1
                        if len(should_be_prepro) < args.sample_mismatches:
                            should_be_prepro.append(rec)

    # ------------ Report ------------
    def pr_hdr(s: str) -> None:
        print(f"\n=== {s} ===")

    print(f"File: {args.jsonl}")
    print(f"Records:       {total:>10,}")
    print(f"Chains:        {len(chains):>10,}")
    print(f"rama null:     {null_count:>10,}")
    print(f"rama non-null: {sum(bucket_counts.values()):>10,}")

    if unknown_labels:
        pr_hdr("UNKNOWN rama_category labels (outside 6-bucket schema)")
        for lab, n in unknown_labels.most_common():
            print(f"  {lab!r}: {n:,}")

    pr_hdr("Backward checks (label -> resname constraint)")
    if not backward_mismatch_counts:
        print("  [OK] all six labels satisfy resname constraints for every record")
    else:
        for rule, n in backward_mismatch_counts.most_common():
            print(f"  [FAIL] {rule}: {n:,} record(s)")
            for ex in backward_mismatches[rule]:
                print(f"    e.g. {ex['file']} model={ex['model']} chain={ex['chain']} "
                      f"resnum={ex['resnum']}{ex['ins']!r} resname={ex['resname']} "
                      f"rama={ex['rama_category']!r}")

    pr_hdr("Resname composition per bucket (sanity)")
    for bucket in ("Gly", "IleVal", "TransPro", "CisPro", "PrePro", "General"):
        cnt = resname_by_bucket[bucket]
        total_b = sum(cnt.values())
        top = cnt.most_common(6)
        top_str = ", ".join(f"{r}={n:,}" for r, n in top)
        print(f"  {bucket:<10} n={total_b:>9,} : {top_str}")

    pr_hdr("peptide_bond consistency (CisPro <-> cis, TransPro <-> trans|null)")
    if pep_mismatch_count == 0:
        print("  [OK] no CisPro/TransPro record disagrees with peptide_bond field")
    else:
        print(f"  [FAIL] {pep_mismatch_count:,} record(s) with cis/trans field disagreement")
        for ex in pep_mismatches:
            print(f"    e.g. {ex['file']} model={ex['model']} chain={ex['chain']} "
                  f"resnum={ex['resnum']} resname={ex['resname']} "
                  f"rama={ex['rama_category']!r} peptide_bond={ex.get('peptide_bond')!r} "
                  f"omega={ex.get('omega')!r}")

    pr_hdr("Forward PrePro check (verifiable via JSONL adjacency)")
    print(f"  PrePro labeled and next residue verified PRO:  {prepro_verified_correct:,}")
    print(f"  PrePro labeled but visible-contiguous next != PRO: {prepro_labeled_but_next_not_pro_count:,}")
    print(f"  Not PrePro but visible-contiguous next == PRO: {should_be_prepro_count:,}")
    print(f"  Unverifiable (next residue absent or non-contiguous in JSONL): {prepro_unverifiable:,}")
    if prepro_labeled_but_next_not_pro_count:
        print("  Sample PrePro-labeled-but-next-not-PRO records:")
        for ex in prepro_labeled_but_next_not_pro:
            print(f"    {ex['file']} model={ex['model']} chain={ex['chain']} "
                  f"resnum={ex['resnum']} resname={ex['resname']}")
    if should_be_prepro_count:
        print("  Sample should-be-PrePro records:")
        for ex in should_be_prepro:
            print(f"    {ex['file']} model={ex['model']} chain={ex['chain']} "
                  f"resnum={ex['resnum']} resname={ex['resname']}")

    pr_hdr("Summary")
    strong_ok = (
        not backward_mismatch_counts
        and pep_mismatch_count == 0
        and not unknown_labels
    )
    forward_ok = (
        prepro_labeled_but_next_not_pro_count == 0
        and should_be_prepro_count == 0
    )
    print(f"  Backward checks (strong):  {'PASS' if strong_ok else 'FAIL'}")
    print(f"  Forward PrePro check (verifiable subset only): {'PASS' if forward_ok else 'FAIL'}")

    return 0 if (strong_ok and forward_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
