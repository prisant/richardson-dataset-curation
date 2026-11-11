#!/usr/bin/env python3
"""Check filter_pruned_residues.py --residues against the unrestricted output.

ALL and SUBSET are the same raw JSONL filtered with --residues all and with
a restricted set.  Expected: SUBSET is a sub-multiset of ALL (so every kept
record is identical in every field), every dropped record's resname is outside the
set, and every kept record's resname is inside it.  Records are compared as
exact JSON lines, as a multiset.

Usage:
    check_residue_set.py ALL.jsonl SUBSET.jsonl [--set std20+MSE]
"""

import argparse
import importlib.util
import json
from collections import Counter
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "filter_pruned_residues", Path(__file__).with_name("filter_pruned_residues.py"))
_filter = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_filter)


def _load(path):
    """Records as a multiset of exact JSON lines (meta lines left out).

    A multiset, not a dict keyed on residue: at microheterogeneity sites
    pydangle writes two records with the same chain, number and insertion
    code (e.g. 1din A123 CSD/CYS, 1ysl B111 CYS/CSD).
    """
    out = Counter()
    with open(path) as fh:
        for line in fh:
            if not json.loads(line).get("_meta"):
                out[line.rstrip("\n")] += 1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("all")
    ap.add_argument("subset")
    ap.add_argument("--set", default="std20+MSE", choices=[k for k, v in _filter.RESIDUE_SETS.items() if v])
    a = ap.parse_args()
    names = _filter.RESIDUE_SETS[a.set]
    full, sub = _load(a.all), _load(a.subset)
    extra = sub - full          # records in SUBSET not identical to one in ALL
    dropped = full - sub
    resname = {ln: json.loads(ln)["resname"] for ln in set(full) | set(sub)}
    by_name = Counter()
    for ln, n in dropped.items():
        by_name[resname[ln]] += n
    wrong_drop = sum(n for ln, n in dropped.items() if resname[ln] in names)
    wrong_keep = sum(n for ln, n in sub.items() if resname[ln] not in names)
    n_extra = sum(extra.values())
    print(f"all {sum(full.values())}  subset {sum(sub.values())}  dropped {sum(dropped.values())}")
    print(f"dropped by resname: {dict(by_name.most_common())}")
    print(f"subset records not identical to an ALL record: {n_extra}")
    print(f"dropped although in set: {wrong_drop}")
    print(f"kept although outside set: {wrong_keep}")
    ok = not (n_extra or wrong_drop or wrong_keep)
    print("OK" if ok else "FAIL")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
