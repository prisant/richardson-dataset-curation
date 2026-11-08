#!/usr/bin/env python3
"""Build an interim top2018 JSONL from Williams' Top2018 SQLite database.

"Shortpath" rebuild of the published pydangle top2018 JSONL (Zenodo
10.5281/zenodo.19053293 full v2.0.0; 10.5281/zenodo.19040026 MC v1.0.0)
without rerunning pydangle.  Same record schema, field order and residue
naming as the published files, so existing readers work unchanged.  It
corrects three defects of the published files:

  * continuity: phi/psi/omega reached across dropped modified residues
    (here: MolProbity values; ramalyze leaves phi/psi null next to
    non-standard residues, so those become null, not wrong);
  * rama_category: Ile/Val checked before pre-Pro (here: MolProbity
    ramalyze rama_res_type, Read et al. 2011 precedence);
  * DSSP: ~1.7-1.8% wrong assignments (here: PDB-REDO DSSP databank,
    DSSP 4.6.1 on deposited coordinates; fallback to our own DSSP 4.6.1
    run where PDB-REDO lacks the residue).

and admits HETATM residues (MSE etc.) that our ATOM-only masks dropped.

Field sources
-------------
  row set        residues passing the dataset filter in the database:
                 full = allatom_filters.all_pass_filters AND chain in
                 allatom_filters_chainlist; mc = mc_filters.mc_pass_filters.
                 Restricted to pydangle's output whitelist (_PROTEIN_NAMES).
  file           "{pdbid}_ersatz.pdb" (published naming)
  phi, psi       ramachandran.phi / .psi
  omega          omega.omega
  tau, chi1-4,   within-residue values: copied from the published record
  chirality      when it is the same residue on the same coordinates (same
                 resname, tau within 0.002 deg) -- unaffected by the fixes;
                 keeps pydangle/Java Dangle conventions incl. aromatic ring
                 chi3/chi4.  Otherwise (new residues; entries whose RCSB
                 coordinates or numbering changed since Top2018):
                 angles.N_CA_C and rotamer.chi1-4 converted from 0..360 to
                 -180..180; chirality kept only if resname matches.
  rama_category  ramachandran.rama_res_type renamed to pydangle labels;
                 Pro relabelled by omega (|omega| < 30 -> CisPro, else
                 TransPro) as pydangle does; rama5/4/3 by pydangle's
                 collapse tables
  peptide_bond   from omega, pydangle rule: |omega| < 30 -> cis, else trans
  dssp           PDB-REDO classic DSSP, else our DSSP 4.6.1 run; residue
                 accepted only if the DSSP one-letter AA matches resname.
                 Loop written as 'C' (published notation); null = no
                 assignment available.

Neighbour rule (mirrors pipeline/filter_pruned_residues.py): phi, omega,
peptide_bond need the peptide-bonded previous residue (connectivity.prev_i)
to pass the same filter; psi needs the next one; rama_* need both.

Outputs: OUT.jsonl.gz, OUT_provenance.json, OUT_report.txt (comparison
with the published file).

Usage:
    build_top2018_shortpath.py --dataset full|mc --db DB --published PUB.jsonl.gz
        --pdbredo-dir DIR --dssp461-dir DIR --out PREFIX
"""

import argparse
import gzip
import hashlib
import json
import math
import sqlite3
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from Bio.Data.IUPACData import protein_letters_3to1
from pydangle_biopython.dssp import parse_dssp_output
from pydangle_biopython.labels import (
    _RAMA6_TO_RAMA3,
    _RAMA6_TO_RAMA4,
    _RAMA6_TO_RAMA5,
)
from pydangle_biopython.measure import _PROTEIN_NAMES

RAMA_NAME = {
    "General": "General",
    "Gly": "Gly",
    "IleVal": "IleVal",
    "pre-Pro": "PrePro",
    "trans-Pro": "TransPro",
    "cis-Pro": "CisPro",
}
LOOP_CODE = "C"  # published top2018 files write DSSP loop as 'C'
WITHIN_RESIDUE = ("tau", "chi1", "chi2", "chi3", "chi4", "chirality")
ONE_LETTER = {k.upper(): v for k, v in protein_letters_3to1.items()}
ONE_LETTER["MSE"] = "M"

FIELDS_FULL = ["file", "model", "chain", "resnum", "ins", "resname", "phi", "psi",
               "omega", "tau", "chi1", "chi2", "chi3", "chi4", "rama_category",
               "rama5", "rama4", "rama3", "dssp", "peptide_bond", "chirality"]
FIELDS_MC = [f for f in FIELDS_FULL if not f.startswith("chi") or f == "chirality"]
ANGLES = ("phi", "psi", "omega", "tau", "chi1", "chi2", "chi3", "chi4")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def signed(angle):
    """0..360 -> (-180, 180]; None passes through."""
    if angle is None:
        return None
    a = float(angle)
    return round(a - 360.0, 3) if a > 180.0 else round(a, 3)


def key_of(rec):
    return (rec["file"], rec["chain"], rec["resnum"], rec["ins"])


# ---------------------------------------------------------------- DSSP
def load_dssp(path):
    """Return {(chain, resnum, icode_stripped): (aa, code_or_None)}."""
    if not path.exists():
        return None
    text = path.read_text()
    codes = parse_dssp_output(text)
    aas = {}
    in_data = False
    for line in text.splitlines():
        if line.strip().startswith("#  RESIDUE AA"):
            in_data = True
            continue
        if not in_data or len(line) < 17 or line[13:14] == "!":
            continue
        try:
            resnum = int(line[5:10].strip())
        except ValueError:
            continue
        aas[(line[11:12], resnum, line[10:11].strip())] = line[13:14]
    return {k: (aas.get(k), v) for k, v in codes.items()}


def aa_matches(resname, aa):
    if aa is None:
        return False
    if aa.islower():  # DSSP marks disulfide-bonded Cys with lowercase letters
        aa = "C"
    expected = ONE_LETTER.get(resname)
    return expected is None or aa == expected or aa == "X"


class DsspSource:
    def __init__(self, pdbredo_dir, ours_dir):
        self.pdbredo_dir = Path(pdbredo_dir)
        self.ours_dir = Path(ours_dir)
        self.pdbid = None
        self.redo = self.ours = None
        self.stats = Counter()

    def lookup(self, pdbid, chain, resnum, icode, resname):
        if pdbid != self.pdbid:
            self.pdbid = pdbid
            self.redo = load_dssp(self.pdbredo_dir / f"{pdbid}.dssp")
            self.ours = load_dssp(self.ours_dir / f"{pdbid}.dssp")
        key = (chain, resnum, icode.strip())
        for name, table in (("pdbredo", self.redo), ("ours_4.6.1", self.ours)):
            if table is not None and key in table:
                aa, code = table[key]
                if aa_matches(resname, aa):
                    self.stats[name] += 1
                    # Published-file notation: loop is 'C'; null means
                    # no assignment available.
                    return LOOP_CODE if code is None else code
                self.stats[f"{name}_aa_mismatch"] += 1
        self.stats["none"] += 1
        return None


# ---------------------------------------------------------------- build
SQL = """
SELECT r.pdbid, r.chain, r.resname, r.resseq, r.icode,
       fp.{flag} AS pass_prev, fn.{flag} AS pass_next,
       ra.phi, ra.psi, ra.rama_res_type, o.omega, a.N_CA_C,
       ro.chi1, ro.chi2, ro.chi3, ro.chi4
FROM residues r
JOIN {table} f ON f.chain_i = r.chain_i AND f.residue_i = r.residue_i
LEFT JOIN connectivity c ON c.chain_i = r.chain_i AND c.residue_i = r.residue_i
LEFT JOIN {table} fp ON fp.chain_i = r.chain_i AND fp.residue_i = c.prev_i
LEFT JOIN {table} fn ON fn.chain_i = r.chain_i AND fn.residue_i = c.next_i
LEFT JOIN ramachandran ra ON ra.chain_i = r.chain_i AND ra.residue_i = r.residue_i
LEFT JOIN omega o ON o.chain_i = r.chain_i AND o.residue_i = r.residue_i
LEFT JOIN angles a ON a.chain_i = r.chain_i AND a.residue_i = r.residue_i
LEFT JOIN rotamer ro ON ro.chain_i = r.chain_i AND ro.residue_i = r.residue_i
WHERE f.{flag} = 1 {chainlist}
ORDER BY r.pdbid, r.chain_i, r.residue_i
"""


def build_records(args, published):
    if args.dataset == "full":
        table, flag = "allatom_filters", "all_pass_filters"
        chainlist = "AND r.chain_i IN (SELECT chain_i FROM allatom_filters_chainlist)"
        fields = FIELDS_FULL
    else:
        table, flag, chainlist, fields = "mc_filters", "mc_pass_filters", "", FIELDS_MC
    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    dssp = DsspSource(args.pdbredo_dir, args.dssp461_dir)
    stats = Counter()
    records = []
    for (pdbid, chain, resname, resseq, icode, pass_prev, pass_next, phi, psi,
         rama, omega, tau, chi1, chi2, chi3, chi4) in con.execute(
            SQL.format(table=table, flag=flag, chainlist=chainlist)):
        stats["db_rows_passing_filter"] += 1
        if resname not in _PROTEIN_NAMES:
            stats[f"skipped_not_whitelisted:{resname}"] += 1
            continue
        prev_ok, next_ok = pass_prev == 1, pass_next == 1
        rec = {
            "file": f"{pdbid}_ersatz.pdb",
            "model": 1,
            "chain": chain,
            "resnum": int(resseq),
            "ins": icode if icode else " ",
            "resname": resname,
            "phi": signed(phi) if prev_ok else None,
            "psi": signed(psi) if next_ok else None,
            "omega": signed(omega) if prev_ok else None,
            "tau": signed(tau),
            "chi1": signed(chi1), "chi2": signed(chi2),
            "chi3": signed(chi3), "chi4": signed(chi4),
        }
        cat = RAMA_NAME.get(rama) if (prev_ok and next_ok and rec["phi"] is not None
                                       and rec["psi"] is not None) else None
        # pydangle / published rule: Pro is CisPro only when |omega| < 30.
        # ramalyze labels some twisted Pro (30 <= |omega| < 150) cis-Pro.
        if cat in ("CisPro", "TransPro") and rec["omega"] is not None:
            rule = "CisPro" if abs(rec["omega"]) < 30.0 else "TransPro"
            if rule != cat:
                stats[f"pro_relabelled_{cat}_to_{rule}"] += 1
                cat = rule
        rec["rama_category"] = cat
        rec["rama5"] = _RAMA6_TO_RAMA5.get(cat) if cat else None
        rec["rama4"] = _RAMA6_TO_RAMA4.get(cat) if cat else None
        rec["rama3"] = _RAMA6_TO_RAMA3.get(cat) if cat else None
        rec["dssp"] = dssp.lookup(pdbid, chain, rec["resnum"], rec["ins"], resname)
        om = rec["omega"]
        rec["peptide_bond"] = None if om is None else ("cis" if abs(om) < 30.0 else "trans")
        # Within-residue values are untouched by the corrected defects:
        # keep the published values (pydangle / Java Dangle conventions,
        # e.g. aromatic ring chi3/chi4) wherever the residue was published;
        # new residues take the database values (chirality unavailable).
        # Copy only when the published row is the same residue on the same
        # coordinates: same resname and tau within 0.002 deg.  Entries whose
        # RCSB coordinates or numbering changed since Top2018 (e.g. 5u97,
        # 4oe9, 5kve) fail this test and keep the database values.
        pub = published.get(key_of(rec))
        same_residue = pub is not None and pub["resname"] == resname
        same_coords = (same_residue and rec["tau"] is not None and pub["tau"] is not None
                       and angle_equal(rec["tau"], pub["tau"]))
        if same_coords:
            for f in WITHIN_RESIDUE:
                if f in pub:
                    rec[f] = pub[f]
            stats["within_residue_from_published"] += 1
        else:
            rec["chirality"] = pub["chirality"] if same_residue else None
            stats["within_residue_from_db_"
                  + ("coords_differ" if same_residue else
                     "resname_differs" if pub else "new_row")] += 1
        records.append({f: rec[f] for f in fields})
    return records, stats, dssp.stats, fields


# ---------------------------------------------------------------- compare
def angle_equal(a, b, tol=0.002):
    if a is None or b is None:
        return a is None and b is None
    return abs((a - b + 180.0) % 360.0 - 180.0) <= tol


def compare(records, published, fields, out):
    new = {key_of(r): r for r in records}
    common = set(new) & set(published)
    only_new = set(new) - set(published)
    only_pub = set(published) - set(new)
    print(f"rows: shortpath {len(new):,}; published {len(published):,}; common {len(common):,}; "
          f"only shortpath {len(only_new):,}; only published {len(only_pub):,}", file=out)
    print("  only shortpath by resname: "
          f"{Counter(new[k]['resname'] for k in only_new).most_common(10)}", file=out)
    print("  only published by resname: "
          f"{Counter(published[k]['resname'] for k in only_pub).most_common(10)}", file=out)
    print("\nper-field agreement on common rows (angles: |diff| <= 0.002 deg):", file=out)
    for f in fields[6:]:
        same = diff = new_null = pub_null = 0
        for k in common:
            a, b = new[k][f], published[k][f]
            ok = angle_equal(a, b) if f in ANGLES else a == b
            if ok:
                same += 1
            elif a is None:
                new_null += 1
            elif b is None:
                pub_null += 1
            else:
                diff += 1
        print(f"  {f:13s} same {same:>9,}  differ {diff:>7,}  null_in_shortpath {new_null:>7,}  "
              f"null_in_published {pub_null:>7,}", file=out)
    trans = Counter((published[k]["rama_category"], new[k]["rama_category"])
                    for k in common if published[k]["rama_category"] != new[k]["rama_category"])
    print("\nrama_category changes published -> shortpath (top 15):", file=out)
    for (a, b), n in trans.most_common(15):
        print(f"  {str(a):9s} -> {str(b):9s} {n:>7,}", file=out)
    return {"common": len(common), "only_shortpath": len(only_new), "only_published": len(only_pub)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dataset", choices=("full", "mc"), required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--published", required=True)
    ap.add_argument("--pdbredo-dir", required=True)
    ap.add_argument("--dssp461-dir", required=True)
    ap.add_argument("--out", required=True, help="output prefix")
    args = ap.parse_args()

    published = {}
    with gzip.open(args.published, "rt") as fh:
        for line in fh:
            rec = json.loads(line)
            if not rec.get("_meta"):
                published[key_of(rec)] = rec
    print(f"published records: {len(published):,}", file=sys.stderr)

    records, stats, dstats, fields = build_records(args, published)
    out_jsonl = Path(f"{args.out}.jsonl.gz")
    with gzip.open(out_jsonl, "wt") as fh:
        for rec in records:
            fh.write(json.dumps(rec, separators=(",", ":")) + "\n")

    report_path = Path(f"{args.out}_report.txt")
    with open(report_path, "w") as rep:
        print(f"top2018 {args.dataset} shortpath build — {datetime.now(timezone.utc).isoformat()}",
              file=rep)
        print(f"records written: {len(records):,}", file=rep)
        print(f"build stats: {dict(sorted(stats.items()))}", file=rep)
        print(f"dssp sources: {dict(sorted(dstats.items()))}", file=rep)
        dssp_null = sum(1 for r in records if r["dssp"] is None)
        print(f"dssp null (no assignment available): {dssp_null:,}\n", file=rep)
        summary = compare(records, published, fields, rep)

    provenance = {
        "description": f"Interim top2018 {args.dataset} JSONL rebuilt from the Williams Top2018 "
                       "SQLite database (shortpath); same schema as the published pydangle file.",
        "created": datetime.now(timezone.utc).isoformat(),
        "script": Path(__file__).name,
        "dataset": args.dataset,
        "records": len(records),
        "inputs": {
            "williams_db": {"path": str(args.db), "zenodo_doi": "10.5281/zenodo.22085722",
                            "md5": md5(args.db)},
            "published": {"path": str(args.published), "md5": md5(args.published)},
            "pdbredo_dssp_dir": str(args.pdbredo_dir),
            "ours_dssp461_dir": str(args.dssp461_dir),
        },
        "dssp_sources": dict(dstats),
        "build_stats": dict(stats),
        "comparison_with_published": summary,
        "output_md5": md5(out_jsonl),
    }
    Path(f"{args.out}_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(report_path.read_text())


if __name__ == "__main__":
    main()
