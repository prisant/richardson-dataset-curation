"""Select the coordinate lines that belong in a residue mask.

Shared by build_ersatz_pdbs.py (top2018 masks from the lab's pruned files)
and build_bfactor_masks.py (top100/500/8000 masks from the ersatz files).

Before 2026-11-08 both builders read ATOM lines only, so a modified residue
written as HETATM (MSE, SMC, KCX, ...) never entered a mask, even when it
was part of the kept chain (defect register B1).  This module keeps every
ATOM line in the requested chains, exactly as before, and adds the HETATM
lines of polymer residues.

A HETATM line is taken as a polymer residue when all of these hold:

* its chain is one of the kept chains;
* its residue name is not water (HOH, DOD, WAT);
* its residue name appears in that chain's SEQRES records (the same
  self-evidencing rule the DSSP cleanup uses: each file declares its own
  protein composition, so no list of modified residues is needed);
* it comes before that chain's first TER record.

The TER test excludes free amino acids modelled as ligands (e.g. a free
MSE with OXT after the chain's TER, 1vje A 301).  The chain test excludes
HETATM records of other chains that the lab's pruned files keep while
deleting those chains' ATOM records (e.g. 1ci4_B_pruned_all.pdb carries
MSE A 1).  Which of the admitted residues pydangle measures is a separate
question (register B2); the mask only records what the lab kept.

Load from a sibling script with importlib (works under ``python -I``)::

    _spec = importlib.util.spec_from_file_location(
        "polymer_residues", Path(__file__).with_name("polymer_residues.py"))
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator

WATER = frozenset({"HOH", "DOD", "WAT"})


def seqres_names(lines: Iterable[str]) -> dict[str, set[str]]:
    """Map chain ID -> residue names listed in that chain's SEQRES records."""
    names: dict[str, set[str]] = {}
    for line in lines:
        if line.startswith("SEQRES"):
            names.setdefault(line[11], set()).update(line[19:].split())
    return names


def mask_lines(
    lines: list[str], chains: Iterable[str] | None = None,
) -> Iterator[str]:
    """Yield the ATOM and polymer HETATM lines of the kept chains, in file order.

    Parameters
    ----------
    lines
        All lines of a PDB file.
    chains
        Kept chain IDs.  ``None`` means the chains that have ATOM records
        (a lab pruned file holds ATOM records for its kept chain only).
    """
    if chains is None:
        kept = {ln[21] for ln in lines if ln.startswith("ATOM  ")}
    else:
        kept = set(chains)
    seq = seqres_names(lines)
    ended: set[str] = set()      # chains whose first TER has been seen
    last_chain = None
    for line in lines:
        if line.startswith("TER"):
            # A TER without a chain ID ends the chain of the preceding atom.
            chain = line[21] if len(line) > 21 and line[21] != " " else last_chain
            if chain is not None:
                ended.add(chain)
            continue
        if line.startswith("ATOM  "):
            last_chain = line[21]
            if line[21] in kept:
                yield line
        elif line.startswith("HETATM"):
            chain = line[21]
            last_chain = chain
            resname = line[17:20].strip()
            if (chain in kept and chain not in ended
                    and resname not in WATER
                    and resname in seq.get(chain, ())):
                yield line
