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
        from openmm.app import Modeller, PDBFile
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

    # Belt-and-suspenders: removeHeterogens() has been observed to still
    # leave a stray cofactor residue behind (FTO's Fe2+ ion, residue name
    # FE2) that gmx pdb2gmx then rejects outright with "Residue 'FE2' not
    # found in residue topology database" / "chain does not appear to
    # contain a recognized chain molecule". Explicitly strip anything left
    # over that isn't a standard amino acid the target force field knows,
    # regardless of what removeHeterogens' own heterogen list caught.
    standard_residues = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "HIP", "HID", "HIE", "HSD", "HSE", "HSP", "CYX", "CYM", "ASH", "GLH", "LYN",
    }
    stray = [r for r in fixer.topology.residues() if r.name not in standard_residues]
    if stray:
        print(
            f"[fix_missing_atoms] stripping {len(stray)} residual non-standard "
            f"residue(s) removeHeterogens missed: {sorted({r.name for r in stray})}"
        )
        modeller = Modeller(fixer.topology, fixer.positions)
        modeller.delete(stray)
        fixer.topology, fixer.positions = modeller.topology, modeller.positions

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
