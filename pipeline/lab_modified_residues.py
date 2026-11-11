#!/usr/bin/env python3
"""How the Richardson Lab Top2018 database treats modified residues.

Reads Christopher Williams' Top2018 SQLite database (Zenodo
10.5281/zenodo.22085722) read-only and reports, for residue types other
than the 20 standard amino acids:

  1. per residue name: residues, rows in `ramachandran`, rows with both
     phi and psi, residues passing the mainchain filter (mc_filters
     .mc_pass_filters), rows passing the Ramachandran filter;
  2. for standard residues bonded to such a residue (the `connectivity`
     table's prev_i / next_i, so chain breaks are not counted): how often
     the neighbour's phi (residue after) or psi (residue before) is
     present, for MSE and for all other non-standard types;
  3. how many cap residues (ACE, FOR, NH2, NME) the `residues` table holds.

Usage:
    lab_modified_residues.py DB [--top 20]
"""

import argparse
import sqlite3

STD20 = ("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
         "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")
CAPS = ("ACE", "FOR", "NH2", "NME")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("db")
    ap.add_argument("--top", type=int, default=20)
    a = ap.parse_args()
    con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    std = ",".join(f"'{n}'" for n in STD20)

    print("1. Non-standard residue types")
    print(f"{'resname':8} {'residues':>9} {'rama_rows':>10} {'phi_psi':>8} {'mc_pass':>8} {'rama_pass':>10}")
    rows = con.execute(f"""
        SELECT r.resname, COUNT(*),
               SUM(ra.chain_i IS NOT NULL),
               SUM(ra.phi IS NOT NULL AND ra.psi IS NOT NULL),
               SUM(m.mc_pass_filters = 1),
               SUM(ra.ramachandran_pass_filters = 1)
        FROM residues r
        LEFT JOIN ramachandran ra USING (chain_i, residue_i)
        LEFT JOIN mc_filters m USING (chain_i, residue_i)
        WHERE r.resname NOT IN ({std})
        GROUP BY r.resname ORDER BY COUNT(*) DESC""").fetchall()
    widths = (9, 10, 8, 8, 10)

    def fmt(vals):
        return " ".join(f"{v or 0:>{w}}" for v, w in zip(vals, widths))

    for name, *vals in rows[: a.top]:
        print(f"{name:8} {fmt(vals)}")
    non_mse = [sum(v or 0 for v in col) for col in zip(*[r[1:] for r in rows if r[0] != "MSE"])]
    print(f"all {len(rows) - 1} non-MSE types: {fmt(non_mse)}")

    print("\n2. Standard residues bonded to a non-standard residue (connectivity)")
    print(f"{'group':16} {'side':14} {'neighbours':>10} {'angle_present':>14}")
    for side, link, col in (("after (phi)", "prev_i", "phi"), ("before (psi)", "next_i", "psi")):
        q = f"""
            SELECT CASE WHEN m.resname = 'MSE' THEN 'MSE' ELSE 'other non-std' END AS grp,
                   COUNT(*), SUM(ra.{col} IS NOT NULL)
            FROM connectivity c
            JOIN residues nb ON nb.chain_i = c.chain_i AND nb.residue_i = c.residue_i
            JOIN residues m ON m.chain_i = c.chain_i AND m.residue_i = c.{link}
            LEFT JOIN ramachandran ra ON ra.chain_i = nb.chain_i AND ra.residue_i = nb.residue_i
            WHERE nb.resname IN ({std}) AND m.resname NOT IN ({std})
            GROUP BY grp ORDER BY grp"""
        for grp, n, present in con.execute(q):
            print(f"{grp:16} {side:14} {n:>10} {present or 0:>14}")

    caps = ",".join(f"'{n}'" for n in CAPS)
    n_caps = con.execute(f"SELECT COUNT(*) FROM residues WHERE resname IN ({caps})").fetchone()[0]
    print(f"\n3. Cap residues ({', '.join(CAPS)}) in the residues table: {n_caps}")


if __name__ == "__main__":
    main()
