#!/usr/bin/env python3
"""Triage an old-vs-new pydangle regression against Java Dangle 1.07.

Inputs are the three outputs of ``regress_continuity.sh``: old.jsonl (pydangle
main), new.jsonl (continuity-fix) and java.txt (``dangle.jar phi psi omega tau``).

For every residue in both pydangle runs and every field in phi, psi, omega:

* changed rows (old != new) are classified by how new compares with Java
  (agree / differ / java_null / new_null / both_null / not_in_java);
* unchanged rows that disagree with Java are counted as pre-existing.

Also lists records present in only one pydangle run and every change of
rama_category / peptide_bond, so that each one can be accounted for.

Residues are keyed by (file, model, chain, resnum, icode); Java Dangle
reports model 1 only, so later models count as not_in_java.

Usage:
    triage_continuity_regression.py OUT_DIR [--new new.jsonl] [--tol 0.002] [--list N]
"""

import argparse
import json
from collections import Counter
from pathlib import Path

FIELDS = ("phi", "psi", "omega")
JAVA_FIELDS = ("phi", "psi", "omega", "tau")


def _num(text):
    text = text.strip()
    if text in ("", "__?__", "null", "NaN"):
        return None
    return float(text)


def load_java(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split(":")
            key = (f[0], int(f[1]), f[2].strip(), int(f[3]), f[4] if f[4].strip() else " ")
            out[key] = dict(zip(JAVA_FIELDS, map(_num, f[6:10]), strict=True))
    return out


def load_jsonl(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            rec = json.loads(line)
            if rec.get("_meta"):
                continue
            out[(rec["file"], rec["model"], rec["chain"], rec["resnum"], rec["ins"] or " ")] = rec
    return out


def same(a, b, tol):
    if a is None or b is None:
        return a is None and b is None
    d = abs(a - b) % 360.0
    return min(d, 360.0 - d) <= tol


def classify(new, java, tol):
    """Compare a pydangle value with the Java value of a residue Java reports."""
    if new is None and java is None:
        return "both_null"
    if new is None:
        return "new_null"
    if java is None:
        return "java_null"
    return "agree" if same(new, java, tol) else "differ"


def fmt(x):
    return "   null" if x is None else f"{x:8.2f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--new", default="new.jsonl", help="fixed-pydangle JSONL in OUT_DIR")
    ap.add_argument("--tol", type=float, default=0.002)
    ap.add_argument("--list", type=int, default=200, help="max rows to list per section")
    args = ap.parse_args()

    old = load_jsonl(args.out_dir / "old.jsonl")
    new = load_jsonl(args.out_dir / args.new)
    java = load_java(args.out_dir / "java.txt")

    only_old = sorted(set(old) - set(new))
    only_new = sorted(set(new) - set(old))
    print(f"old {len(old)}  new {len(new)}  java {len(java)}")
    print(f"only old {len(only_old)}  only new {len(only_new)}")
    for label, keys, src in (("only old", only_old, old), ("only new", only_new, new)):
        for k in keys[: args.list]:
            r = src[k]
            j = java.get(k)
            print(f"  {label:8s} {k[0]} m{k[1]} {k[2]}{k[3]}{k[4].strip()} {r['resname']}"
                  f"  phi{fmt(r.get('phi'))} psi{fmt(r.get('psi'))} omega{fmt(r.get('omega'))}"
                  f"  java={'yes' if j else 'no'}")

    changed = Counter()
    preexisting = Counter()
    changed_rows = []
    for k in sorted(set(old) & set(new)):
        j = java.get(k)
        for f in FIELDS:
            o, n = old[k].get(f), new[k].get(f)
            jv = j[f] if j else None
            if not same(o, n, args.tol):
                c = classify(n, jv, args.tol) if j else "not_in_java"
                changed[(f, c)] += 1
                changed_rows.append((k, f, o, n, jv, c, old[k]["resname"]))
            elif j and not same(n, jv, args.tol):
                preexisting[(f, classify(n, jv, args.tol))] += 1

    print("\n=== changed values (old != new), new vs java ===")
    for f in FIELDS:
        row = {c: n for (ff, c), n in changed.items() if ff == f}
        print(f"  {f:6s} {dict(sorted(row.items()))}")
    print("\n=== unchanged values that disagree with java (pre-existing) ===")
    for f in FIELDS:
        row = {c: n for (ff, c), n in preexisting.items() if ff == f}
        print(f"  {f:6s} {dict(sorted(row.items()))}")

    bad = [r for r in changed_rows if r[5] not in ("agree", "both_null")]
    print(f"\n=== changed values where new does not agree with java: {len(bad)} ===")
    for k, f, o, n, jv, c, rn in bad[: args.list]:
        print(f"  {k[0]} m{k[1]} {k[2]}{k[3]}{k[4].strip()} {rn} {f:5s} old{fmt(o)} new{fmt(n)}"
              f" java{fmt(jv)}  {c}")

    print("\n=== categorical changes ===")
    for k in sorted(set(old) & set(new)):
        for f in ("rama_category", "peptide_bond"):
            if old[k].get(f) != new[k].get(f):
                print(f"  {k[0]} m{k[1]} {k[2]}{k[3]}{k[4].strip()} {old[k]['resname']} {f}:"
                      f" {old[k].get(f)} -> {new[k].get(f)}"
                      f"  omega old{fmt(old[k].get('omega'))} new{fmt(new[k].get('omega'))}")


if __name__ == "__main__":
    main()
