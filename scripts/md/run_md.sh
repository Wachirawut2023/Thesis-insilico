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
#
# Resumability (important on Colab, where sessions crash/disconnect):
#   Every step is skipped if its output already exists, and the production
#   step (the multi-hour one) resumes from its last GROMACS checkpoint via
#   `gmx mdrun -cpi md.cpt -append` instead of restarting the 50 ns run from
#   t=0. Just re-run this script (or run_all.sh) after a crash — it picks up
#   wherever it left off. Set MD_MAXH to make production stop cleanly (and
#   write a clean checkpoint) before a Colab time limit hits, e.g.:
#     MD_MAXH=5 bash run_md.sh METTL3 apo
#   then re-run the same command to continue; repeat until "DONE" is logged.

set -euo pipefail

GENE="${1:?usage: $0 <gene> <mode=apo|bound> [chain] [resi] [partner_resi]}"
MODE="${2:?mode required (apo|bound)}"
CHAIN="${3:-A}"
ANCHOR_RESI="${4:-}"
PARTNER_RESI="${5:-}"

# Max wall-clock hours for the production mdrun before it stops cleanly and
# writes a checkpoint (gmx -maxh). Unset/0 = run until GROMACS finishes or
# the process is killed. On Colab, set this comfortably below your expected
# session length so you always get a clean, resumable checkpoint.
MD_MAXH="${MD_MAXH:-0}"
# Checkpoint-write interval in minutes (gmx mdrun -cpt). Shorter = less lost
# work on an abrupt kill, at the cost of a bit more I/O.
MD_CPT_MIN="${MD_CPT_MIN:-5}"

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

if [ -f md.gro ]; then
  echo "[md] $GENE/$MODE already complete (md.gro exists) — nothing to do"
  exit 0
fi

LOG="$RUN_DIR/run.log"
echo "[md] $(date -Is)  gene=$GENE  mode=$MODE  chain=$CHAIN  anchor=$ANCHOR_RESI partner=$PARTNER_RESI" | tee -a "$LOG"

# Speed: use all CPU threads, offload non-bonded to GPU.
NT=$(nproc)
GPU_FLAGS="-nb gpu -pme gpu -bonded cpu -update gpu"
COMMON_RUN="-v -nt $NT $GPU_FLAGS -cpt $MD_CPT_MIN"

# ── 1. Topology generation ──
# AMBER ff19SB protein + TIP3P water.
# -ignh strips existing hydrogens; pdb2gmx rebuilds them per the force field.
if [ -f protein.gro ]; then
  echo "[md] step 1: pdb2gmx — already done, skipping" | tee -a "$LOG"
else
  echo "[md] step 1: pdb2gmx" | tee -a "$LOG"
  echo "1" | gmx pdb2gmx -f "$PDB_IN" -o protein.gro -p topol.top \
    -i posre.itp -ff amber19sb -water tip3p -ignh \
    >> "$LOG" 2>&1
fi

# ── 2. Define box ──
if [ -f protein_box.gro ]; then
  echo "[md] step 2: editconf — already done, skipping" | tee -a "$LOG"
else
  echo "[md] step 2: editconf (1.0 nm cubic padding)" | tee -a "$LOG"
  gmx editconf -f protein.gro -o protein_box.gro -c -d 1.0 -bt cubic \
    >> "$LOG" 2>&1
fi

# ── 3. Solvate ──
if [ -f protein_solv.gro ]; then
  echo "[md] step 3: solvate — already done, skipping" | tee -a "$LOG"
else
  echo "[md] step 3: solvate" | tee -a "$LOG"
  gmx solvate -cp protein_box.gro -cs spc216.gro -o protein_solv.gro -p topol.top \
    >> "$LOG" 2>&1
fi

# ── 4. Add ions for neutrality + 150 mM NaCl ──
if [ -f protein_ions.gro ]; then
  echo "[md] step 4: add ions — already done, skipping" | tee -a "$LOG"
else
  echo "[md] step 4: add ions" | tee -a "$LOG"
  gmx grompp -f "$MDP_DIR/ions.mdp" -c protein_solv.gro -p topol.top -o ions.tpr -maxwarn 2 \
    >> "$LOG" 2>&1
  echo "SOL" | gmx genion -s ions.tpr -o protein_ions.gro -p topol.top \
    -pname NA -nname CL -neutral -conc 0.15 \
    >> "$LOG" 2>&1
fi

# ── 4b. (bound mode only) inject Mg2+ surrogate at predicted As position + restraints ──
if [ "$MODE" = "bound" ]; then
  if [ -f protein_ions_as.gro ]; then
    echo "[md] step 4b: place Mg surrogate — already done, skipping" | tee -a "$LOG"
  else
    echo "[md] step 4b: place Mg surrogate for As(III) + add distance restraints" | tee -a "$LOG"
    python3 "$SCRIPT_DIR/place_as_surrogate.py" \
      --in-gro protein_ions.gro \
      --out-gro protein_ions_as.gro \
      --topol topol.top \
      --chain "$CHAIN" \
      --anchor-resi "$ANCHOR_RESI" \
      --partner-resi "$PARTNER_RESI" \
      >> "$LOG" 2>&1
  fi
  START_GRO="protein_ions_as.gro"
else
  START_GRO="protein_ions.gro"
fi

# ── 5. Energy minimization ──
if [ -f em.gro ]; then
  echo "[md] step 5: energy minimization — already done, skipping" | tee -a "$LOG"
else
  echo "[md] step 5: energy minimization" | tee -a "$LOG"
  gmx grompp -f "$MDP_DIR/em.mdp" -c "$START_GRO" -p topol.top -o em.tpr -maxwarn 2 \
    >> "$LOG" 2>&1
  gmx mdrun -deffnm em $COMMON_RUN \
    >> "$LOG" 2>&1
fi

# ── 6. NVT equilibration ──
if [ -f nvt.gro ]; then
  echo "[md] step 6: NVT — already done, skipping" | tee -a "$LOG"
else
  echo "[md] step 6: NVT 100 ps @ 310 K" | tee -a "$LOG"
  if [ ! -f nvt.tpr ]; then
    gmx grompp -f "$MDP_DIR/nvt.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2 \
      >> "$LOG" 2>&1
  fi
  if [ -f nvt.cpt ]; then
    echo "[md] resuming NVT from checkpoint" | tee -a "$LOG"
    gmx mdrun -deffnm nvt -cpi nvt.cpt -append $COMMON_RUN >> "$LOG" 2>&1
  else
    gmx mdrun -deffnm nvt $COMMON_RUN >> "$LOG" 2>&1
  fi
fi

# ── 7. NPT equilibration ──
if [ -f npt.gro ]; then
  echo "[md] step 7: NPT — already done, skipping" | tee -a "$LOG"
else
  echo "[md] step 7: NPT 100 ps @ 310 K, 1 bar" | tee -a "$LOG"
  if [ ! -f npt.tpr ]; then
    gmx grompp -f "$MDP_DIR/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2 \
      >> "$LOG" 2>&1
  fi
  if [ -f npt.cpt ]; then
    echo "[md] resuming NPT from checkpoint" | tee -a "$LOG"
    gmx mdrun -deffnm npt -cpi npt.cpt -append $COMMON_RUN >> "$LOG" 2>&1
  else
    gmx mdrun -deffnm npt $COMMON_RUN >> "$LOG" 2>&1
  fi
fi

# ── 8. Production ──
echo "[md] step 8: production 50 ns (this is the long step, ~3-4 hr on RTX 4090; much longer on a Colab GPU)" | tee -a "$LOG"
MAXH_FLAG=""
if [ "$MD_MAXH" != "0" ]; then
  MAXH_FLAG="-maxh $MD_MAXH"
fi
if [ ! -f md.tpr ]; then
  gmx grompp -f "$MDP_DIR/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2 \
    >> "$LOG" 2>&1
fi
if [ -f md.cpt ]; then
  echo "[md] resuming production from checkpoint (md.cpt)" | tee -a "$LOG"
  gmx mdrun -deffnm md -cpi md.cpt -append $COMMON_RUN $MAXH_FLAG >> "$LOG" 2>&1
else
  gmx mdrun -deffnm md $COMMON_RUN $MAXH_FLAG >> "$LOG" 2>&1
fi

if [ -f md.gro ]; then
  echo "[md] DONE $(date -Is)  $GENE/$MODE" | tee -a "$LOG"
  ls -lh md.* | tee -a "$LOG"
else
  echo "[md] production stopped before completion (maxh limit or interrupted) — re-run this script to resume from md.cpt" | tee -a "$LOG"
fi
