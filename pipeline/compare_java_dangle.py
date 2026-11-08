#!/usr/bin/env python3
"""Compare pydangle JSONL angles against Java Dangle 1.07 text output.

Java Dangle is run as ``dangle.jar phi psi omega tau FILE`` (header lines
starting with ``#`` removed).  Records are matched on (file, chain, resnum,
icode).  For each of phi, psi, omega, tau the script counts:

    agree        both values present and within TOL degrees (wrapped)
    differ       both present, |difference| > TOL
    null_py      Java has a value, pydangle null
    null_java    pydangle has a value, Java null

and lists residues present in only one source.  Differing rows can be
written to a TSV for inspection.

Usage:
    compare_java_dangle.py JAVA.txt PYDANGLE.jsonl [--tol 0.002]
                           [--tsv DIFFS.tsv] [--sample N]
"""

import argparse
import csv
import json
from collections import Counter

FIELDS = ("phi", "psi", "omega", "tau")


def _num(text):
    text = text.strip()
    if text in ("", "__?__", "null", "NaN"):
        return None
    return float(text)


def load_java(path):
    records = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split(":")
            # file:model:chain:resnum:icode:resname:phi:psi:omega:tau
            key = (f[0], f[2].strip(), int(f[3]), f[4] if f[4].strip() else " ")
            records[key] = (f[5], dict(zip(FIELDS, map(_num, f[6:10]), strict=True)))
    return records


def load_pydangle(path):
    records = {}
    with open(path) as fh:
        for line in fh:
            rec = json.loads(line)
            if rec.get("_meta"):
                continue
            key = (rec["file"], rec["chain"], rec["resnum"], rec["ins"] or " ")
            records[key] = (rec["resname"], {k: rec.get(k) for k in FIELDS})
    return records


def wrapped(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("java")
    ap.add_argument("pydangle")
    ap.add_argument("--tol", type=float, default=0.002)
    ap.add_argument("--tsv")
    ap.add_argument("--sample", type=int, default=10)
    args = ap.parse_args()

    java = load_java(args.java)
    py = load_pydangle(args.pydangle)
    common = sorted(set(java) & set(py))
    only_java = sorted(set(java) - set(py))
    only_py = sorted(set(py) - set(java))
    print(f"java residues {len(java)}; pydangle residues {len(py)}; "
          f"common {len(common)}; only java {len(only_java)}; only pydangle {len(only_py)}")
    for label, keys, src in (("only java", only_java, java), ("only pydangle", only_py, py)):
        names = Counter(src[k][0] for k in keys)
        if keys:
            print(f"  {label} by resname: {names.most_common(8)}")
            for k in keys[: args.sample]:
                print(f"    {k} {src[k][0]}")

    counts = {f: Counter() for f in FIELDS}
    diffs = []
    for k in common:
        for f in FIELDS:
            a, b = java[k][1][f], py[k][1][f]
            if a is None and b is None:
                counts[f]["both_null"] += 1
            elif b is None:
                counts[f]["null_py"] += 1
                diffs.append((k, py[k][0], f, a, b))
            elif a is None:
                counts[f]["null_java"] += 1
                diffs.append((k, py[k][0], f, a, b))
            elif wrapped(a, b) <= args.tol:
                counts[f]["agree"] += 1
            else:
                counts[f]["differ"] += 1
                diffs.append((k, py[k][0], f, a, b))
    for f in FIELDS:
        c = counts[f]
        print(f"  {f:5s} agree {c['agree']}; differ {c['differ']}; null_py {c['null_py']}; "
              f"null_java {c['null_java']}; both_null {c['both_null']}")
    for d in diffs[: args.sample]:
        print(f"    {d}")
    if args.tsv:
        with open(args.tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["file", "chain", "resnum", "icode", "resname", "field", "java", "pydangle"])
            for (fi, ch, rn, ic), rname, f, a, b in diffs:
                w.writerow([fi, ch, rn, ic, rname, f, a, b])


if __name__ == "__main__":
    main()
