#!/usr/bin/env bash
# Build As(OH)3 PDBQT for AutoDock Vina by writing it directly.
#
# OpenBabel can't reliably handle arsenic:
#   - SMILES route: Gasteiger charge model has no parameters for As
#   - MOL2 route: 'As' isn't a SYBYL atom type
# so we skip conversion and emit the PDBQT verbatim. As is remapped to
# Mg in the AutoDock atom-type column because Vina's default force field
# has no As parameters; Vina ΔG becomes a relative-ranking proxy only.
# Primary readout for inhibition inference is the geometric covalent
# docking in scripts/06_dock_covalent.py.
#
# Geometry: As-O 1.78 A, O-As-O 94 deg (trigonal pyramidal, C3v).
set -euo pipefail
cd "$(dirname "$0")"

cat > as3.pdbqt <<'EOF'
REMARK  As(OH)3 arsenous acid; As remapped to Mg for Vina FF compatibility
REMARK  Geometry: As-O 1.78 A, O-As-O 94 deg (trigonal pyramidal)
ROOT
HETATM    1  Mg  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.450 Mg
HETATM    2  O1  LIG A   1       1.301   0.000  -1.214  1.00  0.00    -0.500 OA
HETATM    3  O2  LIG A   1      -0.651   1.127  -1.214  1.00  0.00    -0.500 OA
HETATM    4  O3  LIG A   1      -0.651  -1.127  -1.214  1.00  0.00    -0.500 OA
HETATM    5  H1  LIG A   1       2.002   0.000  -1.869  1.00  0.00     0.350 HD
HETATM    6  H2  LIG A   1      -1.001   1.733  -1.869  1.00  0.00     0.350 HD
HETATM    7  H3  LIG A   1      -1.001  -1.733  -1.869  1.00  0.00     0.350 HD
ENDROOT
TORSDOF 0
EOF

if [ ! -s as3.pdbqt ]; then
  echo "[prepare_as3] ERROR: write failed"
  exit 1
fi

echo "[prepare_as3] wrote as3.pdbqt ($(wc -l < as3.pdbqt) lines)"
head -20 as3.pdbqt
