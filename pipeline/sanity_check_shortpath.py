#!/usr/bin/env python3
"""Sanity checks on a top2018 shortpath JSONL before release.

Checks the file on its own terms (structure, ranges, internal label
consistency, Ramachandran plausibility per DSSP class, the 1gk8-type
island) against the published file, and recomputes phi/psi/omega with
Java Dangle 1.07 on the lab's Top2018 pruned files for a random sample of
entries (the coordinates the Williams database was built from).

Usage:
    sanity_check_shortpath.py NEW.jsonl.gz PUBLISHED.jsonl.gz DATASET_DIR
        [--dangle JAR] [--sample 150] [--seed 1] [--suffix pruned_all]
"""

import argparse
import gzip
import json
import random
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

RAMA5 = {"General": "General", "Gly": "Gly", "IleVal": "IleVal", "TransPro": "Pro",
         "CisPro": "Pro", "PrePro": "PrePro"}
RAMA4 = {"General": "General", "Gly": "Gly", "IleVal": "General", "TransPro": "Pro",
         "CisPro": "Pro", "PrePro": "PrePro"}
RAMA3 = {"General": "General", "Gly": "Gly", "IleVal": "General", "TransPro": "Pro",
         "CisPro": "Pro", "PrePro": "General"}


def load(path):
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            rows.append(json.loads(line))
    return rows


def key(r):
    return (r["file"], r["chain"], r["resnum"], r["ins"])


def wrap(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def report(name, bad, total, examples=()):
    flag = "OK  " if bad == 0 else "FAIL"
    print(f"[{flag}] {name}: {bad:,} of {total:,}")
    for e in list(examples)[:5]:
        print(f"         e.g. {e}")


def structure_and_ranges(rows, pub_rows):
    keys = Counter(key(r) for r in rows)
    dups = [k for k, n in keys.items() if n > 1]
    report("duplicate residue keys", len(dups), len(rows), dups)
    chains = {(r["file"], r["chain"]) for r in rows}
    pub_chains = {(r["file"], r["chain"]) for r in pub_rows}
    print(f"[info] chains: new {len(chains):,}; published {len(pub_chains):,}; "
          f"only new {len(chains - pub_chains)}; only published {len(pub_chains - chains)}")
    for f in ("phi", "psi", "omega", "chi1", "chi2", "chi3", "chi4"):
        bad = [(key(r), r[f]) for r in rows if r.get(f) is not None and not -180.0 <= r[f] <= 180.0]
        report(f"{f} outside [-180, 180]", len(bad), len(rows), bad)
    bad = [(key(r), r["tau"]) for r in rows if r["tau"] is not None and not 90.0 <= r["tau"] <= 130.0]
    report("tau outside [90, 130]", len(bad), len(rows), bad)


def label_consistency(rows):
    n = len(rows)
    checks = {
        "Gly category on non-GLY / GLY not Gly":
            lambda r: (r["rama_category"] == "Gly") != (r["resname"] == "GLY") and r["rama_category"],
        "Pro category on non-PRO / PRO not Pro":
            lambda r: (r["rama_category"] in ("TransPro", "CisPro")) != (r["resname"] == "PRO")
            and r["rama_category"],
        "IleVal on non-ILE/VAL":
            lambda r: r["rama_category"] == "IleVal" and r["resname"] not in ("ILE", "VAL"),
        "CisPro/TransPro inconsistent with omega (|w|<30 cis)":
            lambda r: r["rama_category"] in ("TransPro", "CisPro") and r["omega"] is not None
            and (r["rama_category"] == "CisPro") != (abs(r["omega"]) < 30.0),
        "peptide_bond inconsistent with omega":
            lambda r: (r["omega"] is None) != (r["peptide_bond"] is None) or (
                r["omega"] is not None and (r["peptide_bond"] == "cis") != (abs(r["omega"]) < 30.0)),
        "rama5/4/3 inconsistent with rama_category":
            lambda r: (r["rama_category"] is None and any(r[f] for f in ("rama5", "rama4", "rama3")))
            or (r["rama_category"] is not None and (
                r["rama5"] != RAMA5[r["rama_category"]] or r["rama4"] != RAMA4[r["rama_category"]]
                or r["rama3"] != RAMA3[r["rama_category"]])),
        "rama_category set but phi or psi null":
            lambda r: r["rama_category"] is not None and (r["phi"] is None or r["psi"] is None),
    }
    for name, test in checks.items():
        bad = [(key(r), r["resname"], r["rama_category"], r.get("omega"), r.get("peptide_bond"))
               for r in rows if test(r)]
        report(name, len(bad), n, bad)
    # The category fix: Ile/Val immediately followed (same chain, resnum+1, no icode) by PRO
    by_key = {key(r): r for r in rows}
    hits = bad = 0
    examples = []
    for r in rows:
        if r["resname"] in ("ILE", "VAL") and r["ins"] == " " and r["rama_category"]:
            nxt = by_key.get((r["file"], r["chain"], r["resnum"] + 1, " "))
            if nxt and nxt["resname"] == "PRO":
                hits += 1
                if r["rama_category"] != "PrePro":
                    bad += 1
                    examples.append((key(r), r["rama_category"]))
    report(f"Ile/Val followed by Pro not labelled PrePro (checked {hits:,})", bad, hits, examples)


def plausibility(rows, label):
    """Fraction of H in the alpha basin and E in the beta basin; island count."""
    h = [r for r in rows if r["dssp"] == "H" and r["phi"] is not None and r["psi"] is not None]
    e = [r for r in rows if r["dssp"] == "E" and r["phi"] is not None and r["psi"] is not None]
    h_ok = sum(1 for r in h if -100 <= r["phi"] <= -30 and -80 <= r["psi"] <= -10)
    e_ok = sum(1 for r in e if -180 <= r["phi"] <= -45 and (r["psi"] >= 90 or r["psi"] <= -150))
    island = sum(1 for r in rows if r["phi"] is not None and r["psi"] is not None
                 and r["rama_category"] == "General" and 10 <= r["phi"] <= 40 and -65 <= r["psi"] <= -30)
    island_h = sum(1 for r in h if 10 <= r["phi"] <= 40 and -65 <= r["psi"] <= -30)
    print(f"[info] {label}: H in alpha basin {h_ok / len(h):.4f} ({len(h):,}); "
          f"E in beta basin {e_ok / len(e):.4f} ({len(e):,}); "
          f"General near (25,-47): {island:,} (of which DSSP H: {island_h:,})")


def java_sample(rows, dataset_dir, jar, n, seed, suffix):
    entries = sorted({r["file"][:4] for r in rows})
    random.Random(seed).shuffle(entries)
    sample = set(entries[:n])
    new = {(r["file"][:4], r["chain"], r["resnum"], r["ins"]): r for r in rows if r["file"][:4] in sample}
    java = {}
    for pid in sorted(sample):
        for f in Path(dataset_dir).glob(f"*/{pid}/{pid}_*_{suffix}.pdb"):
            out = subprocess.run(["java", "-jar", jar, "phi", "psi", "omega", str(f)],
                                 capture_output=True, text=True, check=False).stdout
            for line in out.splitlines():
                if line.startswith("#"):
                    continue
                p = line.split(":")
                java[(pid, p[2].strip(), int(p[3]), p[4] if p[4].strip() else " ")] = [
                    None if x.strip() in ("", "__?__") else float(x) for x in p[6:9]]
    stats = Counter()
    worst = []
    for k, r in new.items():
        if k not in java:
            stats["row_not_in_java"] += 1
            continue
        for f, jv in zip(("phi", "psi", "omega"), java[k]):
            v = r[f]
            if v is None:
                stats[f"{f}_null_in_new"] += 1
            elif jv is None:
                stats[f"{f}_null_in_java"] += 1
            elif wrap(v, jv) <= 0.01:
                stats[f"{f}_agree"] += 1
            else:
                stats[f"{f}_differ"] += 1
                worst.append((round(wrap(v, jv), 3), k, r["resname"], f, v, jv))
    print(f"[info] Java Dangle on lab {suffix} files, {len(sample)} random entries, "
          f"{len(new):,} rows: {dict(sorted(stats.items()))}")
    worst.sort(reverse=True)
    total = sum(v for k, v in stats.items() if k.endswith(("_agree", "_differ")))
    report("phi/psi/omega differing from Java Dangle by > 0.01 deg", len(worst), total, worst)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("new")
    ap.add_argument("published")
    ap.add_argument("dataset_dir")
    ap.add_argument("--dangle", default=str(Path.home() / "bin/javadev-master/dangle/dangle.jar"))
    ap.add_argument("--sample", type=int, default=150)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--suffix", default="pruned_all")
    args = ap.parse_args()
    rows = load(args.new)
    pub = load(args.published)
    print(f"[info] new {len(rows):,} rows; published {len(pub):,} rows")
    structure_and_ranges(rows, pub)
    label_consistency(rows)
    plausibility(rows, "new")
    plausibility(pub, "published")
    java_sample(rows, args.dataset_dir, args.dangle, args.sample, args.seed, args.suffix)


if __name__ == "__main__":
    main()
