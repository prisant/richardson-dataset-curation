#!/usr/bin/env python3
"""Verify the HETATM mask fix before and after rebuilding masks.

Runs the old mask code (a copy of the pre-fix build_ersatz_pdbs.py or
build_bfactor_masks.py, passed with --old) and the new code in this
directory, on every mask in a dataset, without writing anything into it.

Checks
  1. old code reproduces the masks on disk (and, for pruned_bfilt, the
     USER INC fragment records of the synthetic pruned files); for
     top2018_full this includes the three RCSB remaps (--rcsb-remap);
  2. new mask = old mask plus added residues, ATOM residues and their
     order unchanged, and every added residue comes from a HETATM line;
  3. pruned_bfilt: fragment records old vs new, with the change tallied.

With --after, the masks on disk are expected to equal the NEW code's output
instead (run after build_ersatz_pdbs.py --masks-only / build_bfactor_masks.py).

Usage:
    verify_hetatm_masks.py DATASET_DIR SUFFIX --old OLD_SCRIPT
        [--rcsb-remap] [--after] [--added-tsv OUT.tsv] [--bmax 30] [-j 4]
"""

from __future__ import annotations

import argparse
import importlib.util
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


NEW_E = _load("new_build_ersatz_pdbs", HERE / "build_ersatz_pdbs.py")
NEW_B = _load("new_build_bfactor_masks", HERE / "build_bfactor_masks.py")
POLY = _load("polymer_residues", HERE / "polymer_residues.py")
OLD = None  # set in main (inherited by forked workers)


def read_mask(path):
    keys = []
    for ln in open(path):
        toks = ln[1:].split()
        if toks:
            keys.append((ln[0], int(toks[0]), toks[1] if len(toks) > 1 else " "))
    return keys


def read_inc(path):
    out = []
    for ln in open(path):
        if ln.startswith("USER  INC:"):
            p = ln[10:].split(":")
            out.append((p[0].strip(), int(p[1]), p[2].strip() or " ",
                        int(p[4]), p[5].strip() or " "))
    return out


def resnames(lines, chains=None):
    return {(ln[21], int(ln[22:26]), ln[26]): (ln[17:20].strip(), ln[:6])
            for ln in POLY.mask_lines(lines, chains)}


def check_pruned(args):
    entry_dir, pdb_id, chain, pruned, suffix, rcsb, after, _ = args
    res = Counter()
    added = []
    remap = NEW_E._RCSB_REMAPS.get((pdb_id, chain)) if rcsb else None
    old = OLD._build_mask(pruned)
    new = NEW_E._build_mask(pruned)
    if remap:
        old = [remap(k) for k in old]
        new = [remap(k) for k in new]
    disk = read_mask(entry_dir / f"{pdb_id}_{chain}.mask")
    res["disk_eq_" + ("new" if after else "old")] += disk == (new if after else old)
    oldset = set(old)
    res["atom_order_kept"] += [k for k in new if k in oldset] == old
    names = resnames(open(pruned).read().splitlines())
    inv = {(remap(k) if remap else k): v for k, v in names.items()}
    for k in new:
        if k not in oldset:
            rn, rec = inv[k]
            res["added_from_HETATM" if rec == "HETATM" else "added_from_ATOM"] += 1
            added.append((pdb_id, *k, rn, ""))
    res["masks"] += 1
    return res, added, []


def check_bfilt(args):
    entry_dir, pdb_id, chain, _, suffix, rcsb, after, bmax = args
    res = Counter()
    added = []
    ersatz = entry_dir / f"{pdb_id}_ersatz.pdb"
    old_all, old_pass = OLD._parse_residues(ersatz, chain, bmax)
    new_all, new_pass = NEW_B._parse_residues(ersatz, chain, bmax)
    old_fr = OLD._find_fragments(old_all, set(old_pass))
    new_fr = NEW_B._find_fragments(new_all, set(new_pass))
    disk_mask = read_mask(entry_dir / f"{pdb_id}_{chain}.mask")
    disk_fr = read_inc(entry_dir / f"{pdb_id}_{chain}_pruned_bfilt.pdb")
    want_mask, want_fr = (new_pass, new_fr) if after else (old_pass, old_fr)
    tag = "new" if after else "old"
    res[f"disk_mask_eq_{tag}"] += disk_mask == want_mask
    res[f"disk_frag_eq_{tag}"] += [(c.strip() or c, *r) for c, *r in want_fr] == disk_fr
    oldset = set(old_all)
    res["atom_order_kept"] += [k for k in new_all if k in oldset] == old_all
    res["passing_atom_unchanged"] += [k for k in new_pass if k in set(old_pass)] == old_pass
    names = resnames(open(ersatz).read().splitlines(), [chain])
    npass = set(new_pass)
    for k in new_all:
        if k not in oldset:
            rn, rec = names[k]
            res["added_from_HETATM" if rec == "HETATM" else "added_from_ATOM"] += 1
            res["added_pass" if k in npass else "added_fail"] += 1
            added.append((pdb_id, *k, rn, "pass" if k in npass else "fail"))
    res["fragments_old"] += len(old_fr)
    res["fragments_new"] += len(new_fr)
    res["chains_fragments_changed"] += old_fr != new_fr
    res["masks"] += 1
    return res, added, ([(pdb_id, chain)] if old_fr != new_fr else [])


def main():
    global OLD
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dataset", type=Path)
    ap.add_argument("suffix")
    ap.add_argument("--old", required=True, type=Path)
    ap.add_argument("--rcsb-remap", action="store_true")
    ap.add_argument("--after", action="store_true")
    ap.add_argument("--added-tsv")
    ap.add_argument("-j", type=int, default=4)
    ap.add_argument("--bmax", type=float, default=30.0,
                    help="pruned_bfilt B-factor cutoff (top100 uses 40; see its reproduce.sh)")
    a = ap.parse_args()
    OLD = _load("old_builder", a.old)
    bfilt = a.suffix == "pruned_bfilt"

    jobs = []
    for mask in sorted(a.dataset.rglob("*.mask")):
        pdb_id, chain = mask.stem.rsplit("_", 1)
        pruned = mask.parent / f"{pdb_id}_{chain}_{a.suffix}.pdb"
        jobs.append((mask.parent, pdb_id, chain, pruned, a.suffix, a.rcsb_remap, a.after, a.bmax))
    # pruned files without a mask (expected: only the superseded 4ct3_D)
    if not bfilt:
        for pruned in sorted(a.dataset.rglob(f"*_{a.suffix}.pdb")):
            pdb_id, chain = pruned.name[: -len(f"_{a.suffix}.pdb")].rsplit("_", 1)
            if not (pruned.parent / f"{pdb_id}_{chain}.mask").exists():
                print(f"pruned file without mask: {pruned.name}"
                      f"{' (in _RCSB_SKIP)' if (pdb_id, chain) in NEW_E._RCSB_SKIP else ''}")

    if not jobs:
        raise SystemExit(f"no .mask files under {a.dataset}")

    total = Counter()
    added = []
    changed = []
    fn = check_bfilt if bfilt else check_pruned
    with ProcessPoolExecutor(max_workers=a.j) as pool:
        for res, add, ch in pool.map(fn, jobs, chunksize=16):
            total.update(res)
            added += add
            changed += ch
    print(f"{a.dataset} suffix={a.suffix} rcsb_remap={a.rcsb_remap} after={a.after}"
          + (f" bmax={a.bmax}" if bfilt else ""))
    for k in sorted(total):
        print(f"  {k:28s} {total[k]}")
    by = Counter(r[4] for r in added)
    print("  added by resname:", ", ".join(f"{n} {c}" for n, c in by.most_common(12)))
    if changed:
        print("  fragment changes, examples:", changed[:6])
    if a.added_tsv:
        with open(a.added_tsv, "w") as fh:
            fh.write("pdb\tchain\tresnum\ticode\tresname\tbfilt\n")
            for r in added:
                fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
