#!/usr/bin/env python3
"""Classify pydangle vs Java Dangle phi/psi/omega disagreements.

A value disagrees when pydangle (JSONL) and Java Dangle (java.txt, model 1)
both report it and differ by more than --tol degrees (default 0.002; both
programs print three decimals).  Each disagreement is recomputed from the
PDB file for every combination of alternates of its four atoms (all
alternates of each atom, from every residue type at that position) and
classed:

  precision      one combination of atoms reproduces both values within
                 --match-tol: same atoms, different arithmetic or rounding.
                 Split by which program is off from the dihedral computed
                 in double precision from the file's coordinates (off =
                 more than --exact-tol): Java, pydangle or both.  (BioPython
                 holds coordinates in single precision; where the geometry
                 is ill-conditioned, e.g. a bond angle near 180 degrees,
                 that alone moves the dihedral by 0.005 degrees, 1gdj A3.
                 Each precision case is also recomputed from coordinates
                 rounded to single precision; the report counts how many
                 pydangle values that reproduces.)
  microhet       the two values come from different residue types at a
                 microheterogeneity position (two residue names) of any
                 residue of the dihedral
  java_mixed     Java's value needs atoms of two different conformers
                 (labels A and B in one dihedral), pydangle's does not
  conformer      both values come from one conformer each, but different
                 ones (e.g. pydangle B, Java A)
  other          anything else (no combination reproduces one of the
                 values: e.g. a different neighbour residue)

Neighbours are taken in file order of first appearance (model 1), as in
altloc_java_diffs.py.  Reports totals per class and per angle, the |delta|
distribution of the precision class, and with --list N the first N
values of every class other than precision.

Usage:
    java_diff_classes.py FILELIST NEW.jsonl JAVA.txt [--tol 0.002]
        [--match-tol 0.05] [--exact-tol 0.0006] [--list 10] [--tsv OUT.tsv]
"""

import argparse
import itertools
import json
import math
import os
import struct
from collections import Counter, defaultdict

ATOMS = {
    "phi": (("-", "C"), ("0", "N"), ("0", "CA"), ("0", "C")),
    "psi": (("0", "N"), ("0", "CA"), ("0", "C"), ("+", "N")),
    "omega": (("-", "CA"), ("-", "C"), ("0", "N"), ("0", "CA")),
}
FIELDS = ("phi", "psi", "omega")
CLASSES = ("precision", "microhet", "java_mixed", "conformer", "other")
OFF = ("java", "pydangle", "both")
BINS = (0.005, 0.01, 0.02, 0.05)


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


def f32(xyz):
    """Coordinates rounded to single precision, as BioPython stores them."""
    return tuple(struct.unpack("f", struct.pack("f", v))[0] for v in xyz)


def delta(a, b):
    d = abs(a - b) % 360.0
    return min(d, 360.0 - d)


def read_pdb(path):
    """{(chain, resnum, icode): {atom: [(altloc, resname, xyz)]}}, residue order (model 1)."""
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
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            res[key][line[12:16].strip()].append((line[16], line[17:20], xyz))
    return res, order


def combos(res, order, index, key, field):
    """[(value, labels, resnames, points)] for every combination of alternates, or None."""
    i = index[key]
    pos = {"-": i - 1, "0": i, "+": i + 1}
    choices = []
    for which, name in ATOMS[field]:
        j = pos[which]
        if j < 0 or j >= len(order) or name not in res[order[j]]:
            return None
        choices.append([(alt, rn, xyz, which) for alt, rn, xyz in res[order[j]][name]])
    out = []
    for combo in itertools.product(*choices):
        labels = frozenset(a for a, _, _, _ in combo if a != " ")
        names = frozenset((w, rn) for _, rn, _, w in combo)
        pts = tuple(c[2] for c in combo)
        out.append((dihedral(*pts), labels, names, pts))
    return out


def classify(cands, py, jv, tol, exact):
    """(class, which program is off, single precision reproduces pydangle).

    The last two are None except for the precision class.
    """
    if cands is None:
        return "other", None, None
    py_m = [c for c in cands if delta(c[0], py) <= tol]
    j_m = [c for c in cands if delta(c[0], jv) <= tol]
    if not py_m or not j_m:
        return "other", None, None
    both = [c for c in py_m if delta(c[0], jv) <= tol]
    if both:
        c = min(both, key=lambda c: delta(c[0], py) + delta(c[0], jv))
        py_off = delta(c[0], py) > exact
        j_off = delta(c[0], jv) > exact
        single = delta(dihedral(*(f32(p) for p in c[3])), py) <= exact
        return ("precision", "both" if py_off and j_off else "pydangle" if py_off else "java",
                single)
    py_names = set().union(*(c[2] for c in py_m))
    j_names = set().union(*(c[2] for c in j_m))
    if py_names != j_names and not any(c[2] in {d[2] for d in j_m} for c in py_m):
        return "microhet", None, None
    if all(len(c[1]) > 1 for c in j_m) and any(len(c[1]) <= 1 for c in py_m):
        return "java_mixed", None, None
    return "conformer", None, None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("filelist")
    ap.add_argument("new")
    ap.add_argument("java")
    ap.add_argument("--tol", type=float, default=0.002)
    ap.add_argument("--match-tol", type=float, default=0.05)
    ap.add_argument("--exact-tol", type=float, default=0.0006)
    ap.add_argument("--list", type=int, default=10)
    ap.add_argument("--tsv")
    a = ap.parse_args()
    paths = {os.path.basename(p.strip()): p.strip() for p in open(a.filelist) if p.strip()}

    java = {}
    with open(a.java) as fh:
        for line in fh:
            f = line.rstrip("\n").split(":")
            num = [None if t.strip() in ("", "__?__") else float(t) for t in f[6:9]]
            java[(f[0], int(f[1]), f[2].strip(), int(f[3]), f[4] if f[4].strip() else " ")] = num

    diffs = defaultdict(list)
    compared = Counter()
    with open(a.new) as fh:
        for line in fh:
            r = json.loads(line)
            if r.get("_meta") or r["model"] != 1:
                continue
            k = (r["file"], 1, r["chain"], r["resnum"], r["ins"] or " ")
            j = java.get(k)
            if not j:
                continue
            for f, jv in zip(FIELDS, j, strict=True):
                py = r.get(f)
                if py is None or jv is None:
                    continue
                compared[f] += 1
                if delta(py, jv) > a.tol:
                    diffs[r["file"]].append((k, f, py, jv))

    tally = Counter()
    by_field = Counter()
    prec_bins = Counter()
    prec_max = 0.0
    prec_off = Counter()
    prec_single = Counter()
    examples = defaultdict(list)
    entries = defaultdict(set)
    tsv = open(a.tsv, "w") if a.tsv else None
    if tsv:
        tsv.write("file\tchain\tresnum\tins\tfield\tpydangle\tjava\tdelta\tclass\toff\tsingle\n")
    for fn in sorted(diffs):
        res, order = read_pdb(paths[fn])
        index = {key: i for i, key in enumerate(order)}
        for k, f, py, jv in diffs[fn]:
            key = (k[2], k[3], k[4])
            cands = combos(res, order, index, key, f) if key in index else None
            cls, off, single = classify(cands, py, jv, a.match_tol, a.exact_tol)
            d = delta(py, jv)
            tally[cls] += 1
            by_field[(cls, f)] += 1
            entries[cls].add(fn)
            if cls == "precision":
                prec_max = max(prec_max, d)
                prec_bins[next((b for b in BINS if d <= b), "more")] += 1
                prec_off[off] += 1
                prec_single[(off, single)] += 1
            elif len(examples[cls]) < a.list:
                examples[cls].append(f"{fn} {k[2]}{k[3]}{k[4].strip()} {f} py {py:.3f} java {jv:.3f}")
            if tsv:
                tsv.write(f"{fn}\t{k[2]}\t{k[3]}\t{k[4].strip()}\t{f}\t{py}\t{jv}\t{d:.4f}\t{cls}"
                          f"\t{off or ''}\t{'' if single is None else int(single)}\n")
    if tsv:
        tsv.close()

    total = sum(tally.values())
    print(f"values compared: {sum(compared.values())} "
          f"({', '.join(f'{f} {compared[f]}' for f in FIELDS)})")
    print(f"differing values (> {a.tol}): {total}")
    print(f"{'class':12s} {'total':>7s} {'phi':>7s} {'psi':>7s} {'omega':>7s} {'entries':>8s}")
    for cls in CLASSES:
        print(f"{cls:12s} {tally[cls]:7d} " + " ".join(f"{by_field[(cls, f)]:7d}" for f in FIELDS)
              + f" {len(entries[cls]):8d}")
    print(f"precision |delta|: max {prec_max:.3f}; "
          + ", ".join(f"<= {b}: {prec_bins[b]}" for b in BINS) + f", more: {prec_bins['more']}")
    print("precision, off from the coordinates: " + ", ".join(f"{o} {prec_off[o]}" for o in OFF))
    print("precision, pydangle reproduced from single-precision coordinates: "
          + ", ".join(f"{o} {prec_single[(o, True)]}/{prec_off[o]}" for o in OFF))
    for cls in CLASSES[1:]:
        for line in examples[cls]:
            print(f"  {cls:10s} {line}")


if __name__ == "__main__":
    main()
