#!/usr/bin/env python3
"""List residues present in one DSSP output directory but not the other.

Question (2026-09-25): top8000 run through our pipeline gives 1,699
residues in the DSSP 4.6.1 output that are absent from the 4.2.2 output
(and 22 the other way).  Which entries, and what distinguishes them?

For each entry present in both A and B, collects residue keys
(chain, resnum, icode) found in only one side and reports per entry:
count each way, number of chains in the original PDB (ATOM records,
model 1), whether the entry exceeds 26 chains (chain-batch fallback),
whether every one-sided key has a same-(chain,resnum) partner on the
other side differing only in insertion code, and sample keys.

Usage:
    python pipeline/diagnose_dssp_keysets.py A_DIR B_DIR SRC_DIR [LABEL_A LABEL_B] [-o REPORT.md]
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from pydangle_biopython.dssp import parse_dssp_output  # noqa: E402
from run_dssp import collect_entries  # noqa: E402


def n_chains(pdb: Path) -> int:
    chains = set()
    with open(pdb, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM  ") and len(line) > 21:
                chains.add(line[21])
    return len(chains)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("a_dir", type=Path)
    ap.add_argument("b_dir", type=Path)
    ap.add_argument("src_dir", type=Path)
    ap.add_argument("labels", nargs="*", default=["A", "B"])
    ap.add_argument("-o", "--output", type=Path)
    x = ap.parse_args()
    la, lb = (x.labels + ["A", "B"])[:2]

    entries = {e.name: e for e in collect_entries(x.src_dir)}
    ia = {p.stem: p for p in x.a_dir.rglob("*.dssp")}
    ib = {p.stem: p for p in x.b_dir.rglob("*.dssp")}
    rows = []
    tot: Counter[str] = Counter()
    for pid in sorted(set(ia) & set(ib)):
        ka = set(parse_dssp_output(ia[pid].read_text(errors="replace")))
        kb = set(parse_dssp_output(ib[pid].read_text(errors="replace")))
        only_a, only_b = ka - kb, kb - ka
        if not only_a and not only_b:
            continue
        pdb = entries[pid] / f"{pid}.pdb" if pid in entries else None
        nch = n_chains(pdb) if pdb and pdb.exists() else -1
        pa = {(k[0], k[1]) for k in only_a}
        pb = {(k[0], k[1]) for k in only_b}
        icode_only = pa == pb and len(only_a) == len(only_b)
        chains_b = "".join(sorted({k[0] for k in only_b}))
        chains_a = "".join(sorted({k[0] for k in only_a}))
        tot[f"only-{la}"] += len(only_a)
        tot[f"only-{lb}"] += len(only_b)
        tot["entries"] += 1
        tot["entries >26 chains"] += nch > 26
        tot["entries icode-only"] += icode_only
        sample = ", ".join(f"{k[0]}{k[1]}{k[2]}" for k in sorted(only_b | only_a)[:6])
        rows.append((pid, len(only_a), len(only_b), nch, nch > 26, icode_only,
                     chains_a or "-", chains_b or "-", sample))

    lines = [f"# DSSP residue-key differences: {la} vs {lb}", "",
             f"- {la}: `{x.a_dir}`", f"- {lb}: `{x.b_dir}`", "",
             "| Count | Value |", "|---|--:|"]
    lines += [f"| {k} | {v:,} |" for k, v in tot.items()]
    lines += ["", f"| Entry | only-{la} | only-{lb} | Chains | >26 | icode-only | chains only-{la} | chains only-{lb} | Sample keys |",
              "|---|--:|--:|--:|---|---|---|---|---|"]
    for row in sorted(rows, key=lambda r: -(r[1] + r[2])):
        lines.append("| " + " | ".join(str(v) for v in row) + " |")
    text = "\n".join(lines) + "\n"
    print(text[:4000])
    if x.output:
        x.output.write_text(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
