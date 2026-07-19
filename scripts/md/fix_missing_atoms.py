"""Rebuild missing heavy side-chain atoms in a receptor PDB using PDBFixer,
in place of pdb2gmx's `-missing` ideal-geometry filler.

Why this exists: pdb2gmx's `-missing` flag adds missing heavy atoms using a
simple geometric template, which is normally fine — but on a handful of
residues (seen on FTO: LYS121, ASP189, GLN499-Cterm) it produces a
degenerate placement that comes out as NaN coordinates for one hydrogen
built on top of the added atom later, corrupting everything downstream
(editconf/solvate fatal errors: "does not contain a '.'"). PDBFixer's
atom-placement (via OpenMM's Modeller machinery) is more robust for these
same cases. Use this as an opt-in pre-pass for genes that hit that failure
mode, rather than switching every gene off pdb2gmx's simpler/faster path.

Only fills missing *atoms* within residues already present — does not add
missing whole residues (loop building is a separate, much bigger modeling
decision) and does not touch non-standard residues (cofactors like FE2 are
left exactly as upstream prep — scripts/02_prepare_receptors.py — wrote
them).

Usage:
    fix_missing_atoms.py --in-pdb <in> --out-pdb <out>
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from openmm.app import PDBFile
from pdbfixer import PDBFixer


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in-pdb", type=Path, required=True)
    p.add_argument("--out-pdb", type=Path, required=True)
    args = p.parse_args()

    fixer = PDBFixer(filename=str(args.in_pdb))
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # no loop building
    fixer.findNonstandardResidues()
    fixer.nonstandardResidues = []  # leave cofactors (e.g. FE2) untouched
    fixer.findMissingAtoms()
    if not fixer.missingAtoms and not fixer.missingTerminals:
        print("[fix_missing_atoms] no missing heavy atoms found; nothing to do")
    else:
        for res, atoms in fixer.missingAtoms.items():
            print(f"[fix_missing_atoms] {res}: adding {[a.name for a in atoms]}")
        for res, atoms in fixer.missingTerminals.items():
            print(f"[fix_missing_atoms] {res}: adding terminal {atoms}")
    fixer.addMissingAtoms()
    with args.out_pdb.open("w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    print(f"[fix_missing_atoms] wrote {args.out_pdb}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
