"""Repair missing heavy atoms in a receptor PDB before gmx pdb2gmx.

Crystal structures commonly have unresolved side-chain density beyond Cβ
for surface-exposed residues (Gln, Lys, Arg, Glu...). gmx pdb2gmx can only
generate missing *hydrogens* (-ignh) — it has no mechanism to rebuild a
missing heavy atom, and fails outright ("atom CG ... not found") when a
present residue doesn't match its template atom-for-atom.

This only fills in missing atoms *within residues that are present* — it
deliberately does not model in missing residues/loops (real structural
gaps are left alone) and does not add hydrogens (pdb2gmx does that itself,
per the target force field). Chain IDs and residue numbers are preserved
exactly (`keepIds=True`) since downstream stages (covalent-anchor lookup,
the Mg surrogate placement) key off the original crystal numbering.
"""
from __future__ import annotations

import argparse
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-pdb", required=True)
    ap.add_argument("--out-pdb", required=True)
    args = ap.parse_args()

    try:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer
    except ImportError:
        print(
            "ERROR: pdbfixer/openmm not installed. Add 'pdbfixer' to the md "
            "conda env (conda-forge) and retry.",
            file=sys.stderr,
        )
        return 1

    fixer = PDBFixer(filename=args.in_pdb)
    fixer.findMissingResidues()
    # Only repair atoms within residues already present in the structure;
    # do not model in missing loops/termini as new residues.
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    n_missing = sum(len(v) for v in fixer.missingAtoms.values()) + sum(
        len(v) for v in fixer.missingTerminals.values()
    )
    if n_missing:
        print(f"[fix_missing_atoms] adding {n_missing} missing heavy atom(s)")
    fixer.addMissingAtoms()

    with open(args.out_pdb, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    print(f"[fix_missing_atoms] wrote {args.out_pdb}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
