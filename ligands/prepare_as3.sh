#!/usr/bin/env bash
# Convert As(OH)3 mol2 -> pdbqt for AutoDock Vina / AutoDock4.
# Vina lacks an As atom type by default; we map As -> metal-like atom in the PDBQT
# (treated as A or custom). For AutoDock4 covalent docking the receptor side
# carries the As anchor instead, so the ligand PDBQT here is only used by Vina.
set -euo pipefail
cd "$(dirname "$0")"

# OpenBabel: assign Gasteiger charges, output PDBQT
obabel as3.mol2 -O as3.pdbqt --partialcharge gasteiger -xh

# Patch atom type: OpenBabel may write 'As' which Vina rejects.
# Substitute with 'Mg' as a workaround metal atom type (vdW close enough for As;
# results are RELATIVE only — see plan notes).
sed -i 's/ As$/ Mg/' as3.pdbqt || true

echo "Wrote as3.pdbqt"
