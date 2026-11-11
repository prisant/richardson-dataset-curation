#!/usr/bin/env python3
"""Compare pydangle phi/psi with Williams' Top2018 SQLite, per altloc mode.

Runs pydangle-biopython (whatever version is importable, e.g. via
PYTHONPATH=<worktree>/src) on the ersatz files of a random sample of
entries that have alternates (or of the entries given with --pids), in the
given --altloc mode, and compares every phi and psi with the database
(model 1, matched on chain, number, insertion code; values within 0.01 deg
agree).  A residue whose name differs between the database and our
output (the entry is numbered differently, e.g. 4oe9 chain B is offset by
5, or a microheterogeneity site where we keep the standard residue type
and the database names the modified one, 1xpm D111 CYS/SCY) is counted as
"resname differs" and its values are compared separately ("resname
differs, agree" / "resname differs, differ").  Residues are split by whether the database lists alternates for
the residue or a bonded neighbour (`connectivity` prev/next).

Prints the tally, then per entry the number of values that differ or are
null in our output, then the first --list rows.  With --misses FILE, every
row that does not agree is written to FILE as TSV
(pdbid chain resseq icode group field db ours).

Usage:
    PYTHONPATH=... pydangle_vs_williams.py DB DATASET_DIR --altloc first
        [--entries 300] [--seed 1] [--pids 4oe9,3zeu] [-j 1]
        [--list 30] [--misses FILE]
"""

import argparse
import json
import random
import sqlite3
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

from Bio.PDB import PDBParser

from pydangle_biopython import measure

_CFG = {}


def _entry(pid):
    """Compare one entry; return (pid, tally, rows)."""
    tally = Counter()
    rows = []
    files = list(Path(_CFG["dataset_dir"]).glob(f"*/{pid}/{pid}_ersatz.pdb"))
    if not files:
        tally["no_file"] += 1
        return pid, tally, rows
    con = sqlite3.connect(f"file:{_CFG['db']}?mode=ro", uri=True)
    db = {}
    for chain, resseq, icode, resname, alts, phi, psi, p_alt, n_alt in con.execute(
            "select r.chain, r.resseq, r.icode, r.resname, r.altlocs, a.phi, a.psi, "
            "rp.altlocs, rn.altlocs "
            "from residues r join ramachandran a using (chain_i, residue_i) "
            "join connectivity c using (chain_i, residue_i) "
            "left join residues rp on rp.chain_i = c.chain_i and rp.residue_i = c.prev_i "
            "left join residues rn on rn.chain_i = c.chain_i and rn.residue_i = c.next_i "
            "where r.pdbid = ?", (pid,)):
        near = bool(alts) or bool(p_alt) or bool(n_alt)
        db[(chain, int(resseq), (icode or " ").strip() or " ")] = (resname, phi, psi, near)
    con.close()
    structure = PDBParser(QUIET=True).get_structure(pid, str(files[0]))
    lines = measure.process_measurement_commands(
        files[0].name, structure, "phi; psi", output_format="jsonl",
        altloc_mode=_CFG["altloc"])
    ours = {}
    for ln in lines:
        r = json.loads(ln)
        if r["model"] == 1:
            ours.setdefault((r["chain"], r["resnum"], r["ins"]), r)
    for key, (resname, phi, psi, near) in db.items():
        grp = "altloc in window" if near else "no altloc"
        r = ours.get(key)
        renamed = bool(r) and r["resname"] != resname
        if renamed:
            tally[(grp, "resname differs")] += 1
            rows.append((pid, *key, grp, "resname", resname, r["resname"]))
        for field, dbv in (("phi", phi), ("psi", psi)):
            if dbv is None:
                continue
            v = r.get(field) if r else None
            agree = v is not None and abs(((v - float(dbv) + 180) % 360) - 180) <= 0.01
            if renamed:
                tally[(grp, "resname differs, " + ("agree" if agree else "differ"))] += 1
            elif v is None:
                tally[(grp, "ours null")] += 1
                rows.append((pid, *key, grp, field, float(dbv), None))
            elif agree:
                tally[(grp, "agree")] += 1
            else:
                tally[(grp, "differ")] += 1
                rows.append((pid, *key, grp, field, float(dbv), v))
    return pid, tally, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("db")
    ap.add_argument("dataset_dir")
    ap.add_argument("--altloc", required=True)
    ap.add_argument("--entries", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pids", help="comma-separated entries instead of a sample")
    ap.add_argument("-j", type=int, default=1)
    ap.add_argument("--list", type=int, default=30)
    ap.add_argument("--misses", help="write every non-agreeing row here (TSV)")
    a = ap.parse_args()
    _CFG.update(db=a.db, dataset_dir=a.dataset_dir, altloc=a.altloc)

    con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    with_alts = sorted({r[0] for r in con.execute(
        "select distinct pdbid from residues where altlocs != ''")})
    con.close()
    if a.pids:
        sample = a.pids.split(",")
    else:
        random.Random(a.seed).shuffle(with_alts)
        sample = with_alts[: a.entries]

    print(f"pydangle {measure.__file__}")
    print(f"altloc mode {a.altloc}; {len(sample)} entries "
          f"({len(with_alts)} in the database have alternates)")
    if a.j > 1:
        with Pool(a.j) as pool:
            results = pool.map(_entry, sample, chunksize=4)
    else:
        results = [_entry(pid) for pid in sample]
    tally = Counter()
    rows = []
    per_entry = []
    for pid, t, r in results:
        tally.update(t)
        rows += r
        if r:
            vals = [x for x in r if x[5] != "resname"]
            per_entry.append((pid, sum(1 for x in vals if x[4] == "altloc in window"),
                              sum(1 for x in vals if x[4] == "no altloc"),
                              len(r) - len(vals)))
    for k in sorted(tally, key=str):
        print(f"  {str(k):40s} {tally[k]}")
    per_entry.sort(key=lambda x: -(x[1] + x[2]))
    print(f"entries with any non-agreeing value: {len(per_entry)}"
          "  (pdbid, values in altloc window, values elsewhere, resname differs)")
    for row in per_entry:
        print("  ", *row)
    for row in rows[: a.list]:
        print("   row", row)
    if a.misses:
        with open(a.misses, "w") as fh:
            fh.write("pdbid\tchain\tresseq\ticode\tgroup\tfield\tdb\tours\n")
            for row in rows:
                fh.write("\t".join("" if x is None else str(x) for x in row) + "\n")


if __name__ == "__main__":
    main()
