#!/usr/bin/env python3
"""Compare pydangle JSONL outputs before and after the altloc change.

OLD and NEW are pydangle outputs (raw or filtered) of the same files.
Records are keyed on (file, model, chain, resnum, ins, resname), so the two
residue types at a microheterogeneity position stay separate.  Reports
records only in OLD / only in NEW (by resname), and for common records the
number of changed values per field (|delta| > --tol degrees for numbers,
any change for labels; value <-> null counted separately), with the largest
angle changes.

Usage:
    compare_altloc.py OLD.jsonl NEW.jsonl [--tol 0.0005] [--list 15]
"""

import argparse
import json
from collections import Counter

KEY = ("file", "model", "chain", "resnum", "ins", "resname")


def _load(path):
    out = {}
    dup = 0
    with open(path) as fh:
        for line in fh:
            r = json.loads(line)
            if r.get("_meta"):
                continue
            k = tuple(r[f] for f in KEY)
            dup += k in out
            out[k] = r
    return out, dup


def _delta(a, b):
    return abs(((a - b + 180.0) % 360.0) - 180.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--tol", type=float, default=0.0005)
    ap.add_argument("--list", type=int, default=15)
    a = ap.parse_args()
    old, dup_old = _load(a.old)
    new, dup_new = _load(a.new)
    only_old = set(old) - set(new)
    only_new = set(new) - set(old)
    print(f"old {len(old)}  new {len(new)}  duplicate keys old {dup_old} new {dup_new}")
    print(f"only old {len(only_old)} {dict(Counter(k[5] for k in only_old).most_common(10))}")
    print(f"only new {len(only_new)} {dict(Counter(k[5] for k in only_new).most_common(10))}")
    for k in sorted(only_old)[: a.list]:
        print("   only old", k)
    for k in sorted(only_new)[: a.list]:
        print("   only new", k)
    changes = Counter()
    records_changed = 0
    big = []
    for k in set(old) & set(new):
        o, n = old[k], new[k]
        changed = False
        for f in sorted(set(o) | set(n)):
            if f in KEY:
                continue
            ov, nv = o.get(f), n.get(f)
            if ov == nv:
                continue
            if ov is None or nv is None:
                changes[(f, "null<->value")] += 1
                changed = True
            elif isinstance(ov, (int, float)) and isinstance(nv, (int, float)):
                d = _delta(float(ov), float(nv))
                if d > a.tol:
                    changes[(f, "value")] += 1
                    changed = True
                    big.append((d, k, f, ov, nv))
            else:
                changes[(f, "label")] += 1
                changed = True
        records_changed += changed
    print(f"common records with any change: {records_changed}")
    for (f, kind), n in sorted(changes.items()):
        print(f"  {f:14s} {kind:13s} {n}")
    big.sort(reverse=True)
    for d, k, f, ov, nv in big[: a.list]:
        print(f"   {d:8.3f}  {k[0]} {k[2]}{k[3]}{k[4].strip()} {k[5]} {f} {ov} -> {nv}")


if __name__ == "__main__":
    main()
