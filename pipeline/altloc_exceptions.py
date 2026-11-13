#!/usr/bin/env python3
"""Count the two cases pydangle --altloc first does not resolve.

pydangle --altloc first (branch altloc-first) measures each residue in one
conformation and, at microheterogeneity (two amino acids at one position),
keeps one amino acid.  Two cases fall outside that rule:

  Exception 1 (one angle): the atoms of an angle do not all have the
      residue's chosen conformation.
        other  the angle is measured in another conformation, the first
               one all its atoms share (e.g. a side-chain atom only in B
               and C while the residue is measured in A);
        mixed  no conformation is shared by all its atoms, so the angle
               mixes atoms of two conformations.
      Counted with pydangle's own code (imported from PYTHONPATH): every
      geometric measurement of a polymer residue that gives a value, model
      1, with a hook on pydangle_biopython.altloc.vectors.

  Exception 2 (one position): the file models a position as two or more
      amino acids (each with a CA atom, at least one atom with an altloc
      label, same chain, residue number and insertion code, model 1), but
      BioPython's parser keeps fewer of them, so pydangle never sees the
      others.  Each case is reported with the amino acid kept and the one
      pydangle's rule would have chosen (first standard type in file
      order, altloc.STANDARD; else the first).

Only files with at least one altloc label are measured (others cannot have
either exception).  Writes a summary to stdout and, with --tsv, one line
per case.

Usage:
    PYTHONPATH=<pydangle altloc-first>/src python altloc_exceptions.py \
        FILELIST [-c COMMANDS] [-j JOBS] [--tsv OUT.tsv]
"""

import argparse
import multiprocessing
import os
import sys
from collections import Counter

from pydangle_biopython import altloc, main as pmain, measure

DEFAULT_COMMANDS = "phi; psi; omega; tau; chi1; chi2; chi3; chi4"

_events = []          # exception 1 cases of the current measurement
_current = {}         # the measurement being computed


def _vectors_hook(atoms):
    """altloc.vectors, recording exception 1 cases."""
    dis = [a for a in atoms if a.is_disordered()]
    if dis:
        common = altloc._common_labels(atoms)
        preferred = altloc._state["residue"]
        if not common:
            kind = "mixed"
        elif preferred is not None and preferred not in common:
            kind = "other"
        else:
            kind = None
        if kind:
            used = sorted({a.get_altloc() for a in dis}) if kind == "mixed" else []
            _events.append((kind, preferred, "".join(sorted(common or ())), "".join(used)))
        else:
            _events.append(None)
    else:
        _events.append(None)
    return _orig_vectors(atoms)


def _compute_hook(command, residue_list, residue_index, unknown_str):
    """measure.compute_measurement, keeping the case of the value returned."""
    _events.clear()
    out = _orig_compute(command, residue_list, residue_index, unknown_str)
    if command[0] != "label" and out != unknown_str and _events and _events[-1]:
        res = residue_list[residue_index]
        kind, preferred, common, used = _events[-1]
        _current["found"].append(("exc1_" + kind, res.get_parent().id, res.id[1],
                                  res.id[2].strip(), res.get_resname(), command[1],
                                  f"residue {preferred or '-'}; shared {common or '-'}"
                                  + (f"; used {used}" if used else "")))
    return out


_orig_vectors = altloc.vectors
_orig_compute = measure.compute_measurement
altloc.vectors = _vectors_hook
measure.compute_measurement = _compute_hook


def _file_positions(path):
    """{(chain, resnum, icode): {resname: (first serial, has CA, labelled)}}, model 1."""
    pos = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if line[:6] not in ("ATOM  ", "HETATM"):
                continue
            key = (line[21], int(line[22:26]), line[26])
            name = line[17:20].strip()
            try:
                serial = int(line[6:11])
            except ValueError:
                serial = 10**12
            first, ca, lab = pos.setdefault(key, {}).get(name, (serial, False, False))
            pos[key][name] = (min(first, serial), ca or line[12:16].strip() == "CA",
                              lab or line[16] != " ")
    return pos


def _exception2(path, model, found):
    for (chain, num, ins), types in _file_positions(path).items():
        if len(types) < 2 or not all(t[1] for t in types.values()):
            continue
        if not any(t[2] for t in types.values()):
            continue
        kept = set()
        if chain in model:
            for res in model[chain]:
                if res.id[1] == num and res.id[2] == ins:
                    if res.is_disordered() == 2:
                        kept.update(res.disordered_get_id_list())
                    else:
                        kept.add(res.get_resname())
        if len(kept) >= len(types):
            continue
        order = sorted(types, key=lambda n: types[n][0])
        std = [n for n in order if n in altloc.STANDARD]
        rule = (std or order)[0]
        found.append(("exc2_dropped", chain, num, ins.strip(), "/".join(order), "",
                      f"kept {'/'.join(sorted(kept)) or '-'}; rule {rule}; "
                      f"{'same' if rule in kept else 'DIFFERENT'}"))


def _has_altloc(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                return False
            if line[:6] in ("ATOM  ", "HETATM") and line[16] != " ":
                return True
    return False


def process(args):
    path, commands = args
    found = []
    if not _has_altloc(path):
        return os.path.basename(path), False, found, None
    try:
        structure = pmain._parse_pdb_resilient(path)
        for m in list(structure)[1:]:
            structure.detach_child(m.id)
        model = next(iter(structure))
        _exception2(path, model, found)
        _current["found"] = found
        measure.process_measurement_commands(os.path.basename(path), structure, commands,
                                             output_format="jsonl", altloc_mode="first")
    except Exception as exc:  # report, do not hide
        return os.path.basename(path), True, found, repr(exc)
    return os.path.basename(path), True, found, None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("filelist")
    ap.add_argument("-c", "--commands", default=DEFAULT_COMMANDS)
    ap.add_argument("-j", "--jobs", type=int, default=1)
    ap.add_argument("--tsv")
    a = ap.parse_args()
    paths = [p.strip() for p in open(a.filelist) if p.strip()]
    work = [(p, a.commands) for p in paths]
    with multiprocessing.Pool(a.jobs) as pool:
        results = pool.map(process, work, chunksize=4)

    kinds = Counter()
    by_field = Counter()
    entries = {}
    rule_diff = 0
    with_alt = 0
    errors = []
    tsv = open(a.tsv, "w") if a.tsv else None
    if tsv:
        tsv.write("case\tfile\tchain\tresnum\tins\tresname\tfield\tdetail\n")
    for fn, alt, found, err in sorted(results):
        with_alt += alt
        if err:
            errors.append(f"{fn}: {err}")
        for case in found:
            kind, *rest = case
            kinds[kind] += 1
            entries.setdefault(kind, set()).add(fn)
            if kind != "exc2_dropped":
                by_field[(kind, rest[4])] += 1
            elif "DIFFERENT" in rest[5]:
                rule_diff += 1
            if tsv:
                tsv.write("\t".join(map(str, (kind, fn, *rest))) + "\n")
    if tsv:
        tsv.close()

    fields = [c.strip() for c in a.commands.split(";") if c.strip()]
    print(f"files {len(paths)}; with altloc labels {with_alt}; errors {len(errors)}")
    print("exception 1, angles measured in another conformation than the residue's (other)")
    print("            or mixing two conformations (mixed), model 1:")
    for kind in ("exc1_other", "exc1_mixed"):
        print(f"  {kind[5:]:6s} {kinds[kind]:6d} in {len(entries.get(kind, ())):4d} entries; "
              + ", ".join(f"{f} {by_field[(kind, f)]}" for f in fields))
    print(f"exception 2, positions where BioPython keeps fewer amino acids than the file "
          f"models: {kinds['exc2_dropped']} in {len(entries.get('exc2_dropped', ()))} entries;"
          f" amino acid kept differs from pydangle's rule: {rule_diff}")
    for e in errors:
        print("  error", e, file=sys.stderr)


if __name__ == "__main__":
    main()
