#!/usr/bin/env bash
# Run one Tier-3 MD simulation (apo or As-bound) for a single protein.
#
# Usage:
#   bash run_md.sh <gene> <mode> [anchor_chain] [anchor_resi] [partner_resi]
#
# Examples:
#   bash run_md.sh METTL3 apo
#   bash run_md.sh METTL3 bound A 375
#   bash run_md.sh TXN1 bound A 32 35       # bidentate with partner
#
# Requires:
#   - GROMACS 2024+ (gmx in PATH, built with CUDA)
#   - Input PDB at /opt/Thesis-insilico/data/prepared/<gene>.clean.pdb
#   - The mdp/ files in this directory
#
# Bound mode notes:
#   As(III) is not in standard force fields. We approximate by placing a
#   Mg2+ ion at the predicted bonding distance from the target Cys SG
#   (2.25 Å from a single SG for monodentate, or at the bidentate apex
#   for two SGs) and adding distance restraints to maintain the geometry.
#   This tests pose stability, not absolute binding energetics.

set -euo pipefail

GENE="${1:?usage: $0 <gene> <mode=apo|bound> [chain] [resi] [partner_resi]}"
MODE="${2:?mode required (apo|bound)}"
CHAIN="${3:-A}"
ANCHOR_RESI="${4:-}"
PARTNER_RESI="${5:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MDP_DIR="$SCRIPT_DIR/mdp"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
PDB_IN="${PDB_IN:-$REPO_DIR/data/prepared/$GENE.clean.pdb}"
OUT_BASE="${OUT_BASE:-$REPO_DIR/results/md}"
RUN_DIR="$OUT_BASE/${GENE}_${MODE}"

if [ "$MODE" = "bound" ] && [ -z "$ANCHOR_RESI" ]; then
  echo "ERROR: bound mode requires anchor_resi"
  exit 1
fi
if [ ! -f "$PDB_IN" ]; then
  echo "ERROR: input PDB not found: $PDB_IN"
  echo "Run stages 1-2 of the static pipeline first (or pull data/prepared/ from your previous run)."
  exit 1
fi
if ! command -v gmx >/dev/null 2>&1; then
  echo "ERROR: gmx (GROMACS) not in PATH"
  exit 1
fi

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

LOG="$RUN_DIR/run.log"
echo "[md] $(date -Is)  gene=$GENE  mode=$MODE  chain=$CHAIN  anchor=$ANCHOR_RESI partner=$PARTNER_RESI" | tee -a "$LOG"

# Genes needing PDBFixer instead of pdb2gmx's `-missing` for rebuilding
# missing heavy side-chain atoms: on FTO, `-missing`'s ideal-geometry
# filler produced NaN coordinates for one hydrogen on 3 disordered
# residues (LYS121, ASP189, GLN499-Cterm), which corrupted every step
# downstream. PDBFixer's placement is more robust for these same cases.
# Gene-scoped opt-in (like FORCE_HIE_GENES below) since pdb2gmx's
# `-missing` is fine — faster, no extra dependency — everywhere else.
PDBFIXER_GENES=" FTO "
if [[ "$PDBFIXER_GENES" == *" $GENE "* ]]; then
  echo "[md] step 0: PDBFixer pre-pass (rebuild missing heavy atoms)" | tee -a "$LOG"
  python3 "$SCRIPT_DIR/fix_missing_atoms.py" --in-pdb "$PDB_IN" --out-pdb "$RUN_DIR/${GENE}.pdbfixer.pdb" \
    >> "$LOG" 2>&1
  PDB_IN="$RUN_DIR/${GENE}.pdbfixer.pdb"
fi

# Speed: use all CPU threads, offload non-bonded to GPU.
NT=$(nproc)
# No "-update gpu": OPC water uses a virtual site (massless MW), and
# GROMACS's GPU update explicitly doesn't support virtual sites at all
# ("Virtual sites are not supported") — this applies to every dynamical
# step (NVT/NPT/production), not just EM, so update stays on CPU throughout.
GPU_FLAGS="-nb gpu -pme gpu -bonded cpu"
COMMON_RUN="-v -nt $NT $GPU_FLAGS"
# EM uses a non-dynamical integrator (steep/cg) — GROMACS rejects GPU PME
# for that (on top of the update restriction above), so EM can only
# offload nonbonded to the GPU.
EM_RUN="-v -nt $NT -nb gpu -bonded cpu"

# ── 1. Topology generation ──
# AMBER ff19SB protein + OPC water (ff19SB's authors recommend OPC/OPC3,
# not TIP3P — GROMACS's own force-field docs say the same; ff19SB's CMAP
# backbone corrections weren't validated against TIP3P).
# -ignh strips existing hydrogens; pdb2gmx rebuilds them per the force field.
# -missing: several targets (e.g. METTL3, FTO, ALKBH5) have surface
# side chains unresolved past CB in the crystal structure (high B-factor
# disorder, not a data-prep bug — checked the raw PDBs directly). Without
# -missing, pdb2gmx fatal-errors instead of rebuilding those atoms; with
# it, the missing heavy atoms are added by ideal geometry and relaxed by
# the EM step that already runs right after this.
#
# Genes needing an explicit histidine protonation state: pdb2gmx's
# automatic detector (hizzie) infers HID/HIE/HIP from each histidine's
# local H-bonding, but needs the ring atoms present to do it — checked
# ring completeness only, before -missing gets a chance to rebuild them.
# A histidine with the whole ring unresolved in the crystal structure
# (not just CB, like the -missing cases above) fails "Incomplete ring"
# even with -missing. -his forces interactive selection for *every*
# histidine in the structure (pdb2gmx has no per-residue override), so
# only opt genes into this when they actually hit that error — for
# everything else, automatic per-residue detection is the better default.
# User-approved default for the current case (METTL3/HIS116, ring fully
# absent, not obviously catalytic): HIE for all of METTL3's histidines.
FORCE_HIE_GENES=" METTL3 "
PDB2GMX_STDIN="1"
if [[ "$FORCE_HIE_GENES" == *" $GENE "* ]]; then
  N_HIS=$(awk '$1=="ATOM" && $4=="HIS" {print substr($0,22,1) substr($0,23,4)}' "$PDB_IN" | sort -u | wc -l)
  PDB2GMX_HIS_FLAG="-his"
  PDB2GMX_STDIN="$(printf '1\n%.0s' $(seq 1 $((N_HIS + 1))))"
else
  PDB2GMX_HIS_FLAG=""
fi
echo "[md] step 1: pdb2gmx" | tee -a "$LOG"
echo "$PDB2GMX_STDIN" | gmx pdb2gmx -f "$PDB_IN" -o protein.gro -p topol.top \
  -i posre.itp -ff amber19sb -water opc -ignh -missing $PDB2GMX_HIS_FLAG \
  >> "$LOG" 2>&1

# ── 2. Define box ──
echo "[md] step 2: editconf (1.0 nm cubic padding)" | tee -a "$LOG"
gmx editconf -f protein.gro -o protein_box.gro -c -d 1.0 -bt cubic \
  >> "$LOG" 2>&1

# ── 3. Solvate ──
# OPC is a 4-point model (OW/HW1/HW2/MW) like TIP4P, not 3-point like SPC —
# the solvent box coordinate file's per-molecule atom count has to match,
# so this must be tip4p.gro, not the default spc216.gro (3-point).
echo "[md] step 3: solvate" | tee -a "$LOG"
gmx solvate -cp protein_box.gro -cs tip4p.gro -o protein_solv.gro -p topol.top \
  >> "$LOG" 2>&1

# ── 4. Add ions for neutrality + 150 mM NaCl ──
echo "[md] step 4: add ions" | tee -a "$LOG"
gmx grompp -f "$MDP_DIR/ions.mdp" -c protein_solv.gro -p topol.top -o ions.tpr -maxwarn 2 \
  >> "$LOG" 2>&1
echo "SOL" | gmx genion -s ions.tpr -o protein_ions.gro -p topol.top \
  -pname NA -nname CL -neutral -conc 0.15 \
  >> "$LOG" 2>&1

# ── 4b. (bound mode only) inject Mg2+ surrogate at predicted As position + restraints ──
if [ "$MODE" = "bound" ]; then
  echo "[md] step 4b: place Mg surrogate for As(III) + add distance restraints" | tee -a "$LOG"
  python3 "$SCRIPT_DIR/place_as_surrogate.py" \
    --in-gro protein_ions.gro \
    --out-gro protein_ions_as.gro \
    --topol topol.top \
    --chain "$CHAIN" \
    --anchor-resi "$ANCHOR_RESI" \
    --partner-resi "$PARTNER_RESI" \
    >> "$LOG" 2>&1
  START_GRO="protein_ions_as.gro"
else
  START_GRO="protein_ions.gro"
fi

# ── 5. Energy minimization ──
echo "[md] step 5: energy minimization" | tee -a "$LOG"
gmx grompp -f "$MDP_DIR/em.mdp" -c "$START_GRO" -p topol.top -o em.tpr -maxwarn 2 \
  >> "$LOG" 2>&1
gmx mdrun -deffnm em $EM_RUN \
  >> "$LOG" 2>&1

# ── 6. NVT equilibration ──
echo "[md] step 6: NVT 100 ps @ 310 K" | tee -a "$LOG"
gmx grompp -f "$MDP_DIR/nvt.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2 \
  >> "$LOG" 2>&1
gmx mdrun -deffnm nvt $COMMON_RUN \
  >> "$LOG" 2>&1

# ── 7. NPT equilibration ──
echo "[md] step 7: NPT 100 ps @ 310 K, 1 bar" | tee -a "$LOG"
gmx grompp -f "$MDP_DIR/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2 \
  >> "$LOG" 2>&1
gmx mdrun -deffnm npt $COMMON_RUN \
  >> "$LOG" 2>&1

# ── 8. Production ──
echo "[md] step 8: production 50 ns (this is the long step, ~3-4 hr on RTX 4090)" | tee -a "$LOG"
gmx grompp -f "$MDP_DIR/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2 \
  >> "$LOG" 2>&1
gmx mdrun -deffnm md $COMMON_RUN \
  >> "$LOG" 2>&1

echo "[md] DONE $(date -Is)  $GENE/$MODE" | tee -a "$LOG"
ls -lh md.* | tee -a "$LOG"
