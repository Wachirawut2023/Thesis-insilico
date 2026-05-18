#!/usr/bin/env bash
# Build As(OH)3 PDBQT for AutoDock Vina. Tries SMILES first (most reliable),
# falls back to as3.mol2 if SMILES path fails.
#
# Vina has no native As atom type — we remap to Mg in the PDBQT (metal-like
# vdW radius). Vina ΔG values are RELATIVE ranking only; the primary readout
# for inhibition inference is geometric covalent docking in scripts/06.
set -euo pipefail
cd "$(dirname "$0")"

OUT=as3.pdbqt
rm -f "$OUT"

# Route 1: SMILES -> 3D -> PDBQT
if obabel -:"O[As](O)O" -O "$OUT" --gen3d --partialcharge gasteiger 2>/tmp/obabel.log; then
  if [ -s "$OUT" ]; then
    echo "[prepare_as3] generated $OUT from SMILES"
  fi
fi

# Route 2: fallback to mol2 if SMILES route produced nothing
if [ ! -s "$OUT" ]; then
  echo "[prepare_as3] SMILES route empty, trying as3.mol2"
  cat /tmp/obabel.log 2>/dev/null || true
  obabel as3.mol2 -O "$OUT" --partialcharge gasteiger
fi

# Patch As -> Mg for Vina compatibility (Vina parameter file lacks As).
sed -i 's/\bAs\b/Mg/g' "$OUT"

if [ -s "$OUT" ]; then
  echo "[prepare_as3] wrote $OUT ($(wc -l < "$OUT") lines)"
  grep -E '^(ATOM|HETATM)' "$OUT" | head -10
else
  echo "[prepare_as3] ERROR: $OUT is empty after both routes"
  exit 1
fi
