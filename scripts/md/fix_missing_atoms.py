"""Clean up a receptor PDB for gmx pdb2gmx: strip cofactor HETATMs and
repair missing heavy atoms.

02_prepare_receptors.py deliberately keeps catalytic cofactor HETATMs
(SAM, Fe, 2OG, Zn, Mg, m6A...) for Tier 1-2's docking/geometric analysis.
gmx pdb2gmx has no residue templates for those and can't process them —
worse, they show up as spurious 1-atom "chains" that make it fail outright
("This chain does not appear to contain a recognized chain molecule").
Per docs/md_tier3.md, Tier-3 MD doesn't model natural cofactor chemistry
anyway (arsenic binding is tested via a Mg2+ surrogate + distance
restraints, not the real cofactor), so cofactors are simply removed before
gmx sees the structure.

Crystal structures also commonly have unresolved side-chain density beyond
Cβ for surface-exposed residues (Gln, Lys, Arg, Glu...). gmx pdb2gmx can
only generate missing *hydrogens* (-ignh) — it has no mechanism to rebuild
a missing heavy atom, and fails outright ("atom CG ... not found") when a
present residue doesn't match its template atom-for-atom.

Heavy-atom repair only fills in atoms *within residues that are present*
— it deliberately does not model in missing residues/loops (real
structural gaps are left alone) and does not add hydrogens (pdb2gmx does
that itself, per the target force field). Chain IDs and residue numbers
are preserved exactly (`keepIds=True`) since downstream stages
(covalent-anchor lookup, the Mg surrogate placement) key off the original
crystal numbering.
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
    fixer.removeHeterogens(keepWater=False)
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
