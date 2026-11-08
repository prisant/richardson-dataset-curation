#!/usr/bin/env python3
"""Explain pydangle vs Java Dangle phi/psi/omega disagreements by altloc choice.

For each residue where pydangle (new.jsonl) and Java Dangle (java.txt) both
report a value and differ by more than TOL, the angle is recomputed from the
ersatz PDB under three atom-selection rules:

    maxocc   per atom, the altloc with the highest occupancy (BioPython default;
             ties keep the first altloc listed)
    first    per atom, the first altloc listed in the file
    altA     altloc 'A' (or blank) where present, else the only/first altloc

and the script reports which rule reproduces pydangle and which reproduces Java.

Only model 1 is compared (Java Dangle reports model 1 only).

Usage:
    altloc_java_diffs.py OUT_DIR [--tol 0.002] [--list N]
OUT_DIR holds filelist.txt, new.jsonl and java.txt from regress_continuity.sh.
"""

import argparse
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

ATOMS = {
    "phi": (("-", "C"), ("0", "N"), ("0", "CA"), ("0", "C")),
    "psi": (("0", "N"), ("0", "CA"), ("0", "C"), ("+", "N")),
    "omega": (("-", "CA"), ("-", "C"), ("0", "N"), ("0", "CA")),
}


def dihedral(p0, p1, p2, p3):
    b0 = [a - b for a, b in zip(p0, p1, strict=True)]
    b1 = [a - b for a, b in zip(p2, p1, strict=True)]
    b2 = [a - b for a, b in zip(p3, p2, strict=True)]
    n1 = math.sqrt(sum(x * x for x in b1))
    b1 = [x / n1 for x in b1]
    d0 = sum(a * b for a, b in zip(b0, b1, strict=True))
    d2 = sum(a * b for a, b in zip(b2, b1, strict=True))
    v = [a - d0 * b for a, b in zip(b0, b1, strict=True)]
    w = [a - d2 * b for a, b in zip(b2, b1, strict=True)]
    x = sum(a * b for a, b in zip(v, w, strict=True))
    c = (b1[1] * v[2] - b1[2] * v[1], b1[2] * v[0] - b1[0] * v[2], b1[0] * v[1] - b1[1] * v[0])
    y = sum(a * b for a, b in zip(c, w, strict=True))
    return math.degrees(math.atan2(y, x))


def read_pdb(path):
    """Return {(chain, resnum, icode): {atom: [(altloc, occ, xyz), ...]}} and residue order."""
    res = defaultdict(lambda: defaultdict(list))
    order = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if line[:6] not in ("ATOM  ", "HETATM"):
                continue
            key = (line[21], int(line[22:26]), line[26])
            if key not in res:
                order.append(key)
            res[key][line[12:16].strip()].append(
                (line[16], float(line[54:60]), (float(line[30:38]), float(line[38:46]), float(line[46:54])))
            )
    return res, order


def pick(alts, rule):
    if rule == "first":
        return alts[0][2]
    if rule == "maxocc":
        best = alts[0]
        for a in alts[1:]:
            if a[1] > best[1]:
                best = a
        return best[2]
    for a in alts:  # altA
        if a[0] in ("A", " "):
            return a[2]
    return alts[0][2]


def angle(res, order, key, field, rule, index=None):
    i = index[key] if index is not None else order.index(key)
    pos = {"-": i - 1, "0": i, "+": i + 1}
    pts = []
    for which, name in ATOMS[field]:
        j = pos[which]
        if j < 0 or j >= len(order) or name not in res[order[j]]:
            return None
        pts.append(pick(res[order[j]][name], rule))
    return dihedral(*pts)


def close(a, b, tol):
    if a is None or b is None:
        return False
    d = abs(a - b) % 360.0
    return min(d, 360.0 - d) <= tol


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--tol", type=float, default=0.002)
    ap.add_argument("--list", type=int, default=0)
    args = ap.parse_args()
    paths = {os.path.basename(p.strip()): p.strip() for p in open(args.out_dir / "filelist.txt")}

    java = {}
    for line in open(args.out_dir / "java.txt"):
        f = line.rstrip("\n").split(":")
        num = [None if t.strip() in ("", "__?__") else float(t) for t in f[6:9]]
        java[(f[0], int(f[1]), f[2].strip(), int(f[3]), f[4] if f[4].strip() else " ")] = dict(
            zip(("phi", "psi", "omega"), num, strict=True))

    diffs = defaultdict(list)
    for line in open(args.out_dir / "new.jsonl"):
        r = json.loads(line)
        if r.get("_meta"):
            continue
        k = (r["file"], r["model"], r["chain"], r["resnum"], r["ins"] or " ")
        if r["model"] != 1:  # read_pdb reads model 1 only; Java reports model 1 only
            continue
        j = java.get(k)
        if not j:
            continue
        for f in ("phi", "psi", "omega"):
            if r.get(f) is not None and j[f] is not None and not close(r[f], j[f], args.tol):
                diffs[r["file"]].append((k, f, r[f], j[f]))

    rules = ("maxocc", "first", "altA")
    tally = Counter()
    listed = 0
    for fn, rows in sorted(diffs.items()):
        res, order = read_pdb(paths[fn])
        for k, f, py, jv in rows:
            key = (k[2], k[3], k[4])
            vals = {ru: angle(res, order, key, f, ru) for ru in rules}
            py_rule = "+".join(ru for ru in rules if close(vals[ru], py, 0.01)) or "none"
            j_rule = "+".join(ru for ru in rules if close(vals[ru], jv, 0.01)) or "none"
            tally[(py_rule, j_rule)] += 1
            if listed < args.list:
                listed += 1
                print(f"  {fn} {k[2]}{k[3]}{k[4].strip()} {f:5s} py {py:8.2f} java {jv:8.2f}"
                      f"  py~{py_rule}  java~{j_rule}")
    total = sum(tally.values())
    print(f"differing values: {total}")
    print("pydangle reproduced by / Java reproduced by -> count")
    for (pr, jr), n in tally.most_common():
        print(f"  {pr:20s} / {jr:20s} {n}")


if __name__ == "__main__":
    main()
