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
#   - GROMACS >= 2026.3 (gmx in PATH), built with GPU offload support —
#     either a CUDA build (NVIDIA: conda-forge's `nompi_cuda*` build string,
#     see infra/gcp-gpu/ and infra/runpod/) or a HIP build (AMD: built from
#     source against ROCm, see infra/amd-gpu/). GROMACS's `-nb gpu -pme gpu`
#     offload is backend-agnostic, so this script runs unchanged on either.
#     2026.3+ specifically because that's the release that added a native
#     port of AMBER's ff19SB (this pipeline's protein force field) — earlier
#     releases fail pdb2gmx with "force field 'amber19sb' not found".
#   - Input PDB at <repo>/data/prepared/<gene>.clean.pdb
#   - The mdp/ files in this directory
#
# Bound mode notes:
#   As(III) is not in standard force fields. We approximate by placing a
#   Mg2+ ion at the predicted bonding distance from the target Cys SG
#   (2.25 Å from a single SG for monodentate, or at the bidentate apex
#   for two SGs) and adding distance restraints to maintain the geometry.
#   This tests pose stability, not absolute binding energetics.
#
# Resumability (important on Colab, where sessions crash/disconnect, and on
# any long unattended cloud run in general):
#   Every step is skipped if its output already exists, and production (the
#   multi-hour step) is segmented into MD_PROD_MAXH_HOURS-long `gmx mdrun`
#   processes that resume from the last checkpoint via `-cpi md.cpt -append`
#   instead of restarting the 50 ns run from t=0. Just re-run this script
#   (or run_all.sh) after a crash — it picks up wherever it left off.
#
# Device selection:
#   MD_DEVICE=auto (default) — use the GPU if one is visible via
#   `nvidia-smi -L` (NVIDIA) or `rocminfo` (AMD), otherwise fall back to
#   CPU-only offload. Override with MD_DEVICE=gpu or MD_DEVICE=cpu to force
#   a path (e.g. to force GPU and fail loudly if one isn't actually
#   available, or to force CPU-only as a last resort if the GPU backend
#   can't see a device in your environment, e.g. MD_DEVICE=cpu bash run_md.sh TXN1 apo).

set -euo pipefail

GENE="${1:?usage: $0 <gene> <mode=apo|bound> [chain] [resi] [partner_resi]}"
MODE="${2:?mode required (apo|bound)}"
CHAIN="${3:-A}"
ANCHOR_RESI="${4:-}"
PARTNER_RESI="${5:-}"

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

# ── 0. Repair missing heavy atoms + strip cofactor HETATMs ──
# Crystal structures (e.g. METTL3/METTL14, FTO, ALKBH5) commonly have
# unresolved side-chain density beyond Cβ for surface residues, and
# 02_prepare_receptors.py deliberately keeps catalytic cofactor HETATMs
# (SAM, Fe, 2OG, Zn, Mg, m6A...) for Tier 1-2's docking analysis that gmx
# pdb2gmx has no template for. Run every receptor through PDBFixer first
# (scripts/md/fix_missing_atoms.py) so neither of those varies by protein —
# see that script's docstring for the full rationale.
FIXED_PDB="protein_fixed.pdb"
if [ -f "$FIXED_PDB" ]; then
  echo "[md] step 0: repair missing heavy atoms — already done, skipping" | tee -a "$LOG"
else
  echo "[md] step 0: repair missing heavy atoms + strip cofactors (PDBFixer)" | tee -a "$LOG"
  python3 "$SCRIPT_DIR/fix_missing_atoms.py" --in-pdb "$PDB_IN" --out-pdb "$FIXED_PDB" \
    >> "$LOG" 2>&1
fi
PDB_FOR_GMX="$RUN_DIR/$FIXED_PDB"

# Speed: use all CPU threads, offload non-bonded (+ PME, on GPU) work. NT is
# overridable so a second, smaller job can share the box with a larger
# already-running one without fully oversubscribing the CPU.
NT="${NT:-$(nproc)}"

MD_DEVICE="${MD_DEVICE:-auto}"
case "$MD_DEVICE" in
  auto)
    if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L 2>/dev/null | grep -q '^GPU'; then
      USE_GPU=1
    elif command -v rocminfo >/dev/null 2>&1 && rocminfo 2>/dev/null | grep -qE 'gfx[0-9]{2,4}[a-z]?'; then
      USE_GPU=1
    else
      USE_GPU=0
    fi
    ;;
  gpu) USE_GPU=1 ;;
  cpu) USE_GPU=0 ;;
  *)
    echo "ERROR: MD_DEVICE must be auto, gpu, or cpu (got: $MD_DEVICE)"
    exit 1
    ;;
esac

if [ "$USE_GPU" = "1" ]; then
  echo "[md] device: GPU (nb+pme offloaded, bonded on CPU)" | tee -a "$LOG"
  # No "-update gpu": OPC water uses a virtual site (massless MW), and
  # GROMACS's GPU update explicitly doesn't support virtual sites at all
  # ("Virtual sites are not supported") — this applies to every dynamical
  # step (NVT/NPT/production), not just EM, so update stays on CPU throughout.
  GPU_FLAGS="-nb gpu -pme gpu -bonded cpu"
  COMMON_RUN="-v -nt $NT $GPU_FLAGS -cpt $MD_CPT_MIN"
  # EM uses a non-dynamical integrator (steep/cg) — GROMACS rejects GPU PME
  # for that (on top of the update restriction above), so EM can only
  # offload nonbonded to the GPU.
  EM_RUN="-v -nt $NT -nb gpu -bonded cpu -cpt $MD_CPT_MIN"
else
  echo "[md] device: CPU-only (no GPU detected / MD_DEVICE=cpu)" | tee -a "$LOG"
  # CPU-only: no GPU flags at all. GROMACS's Verlet scheme (the only
  # scheme it supports since 2020) runs the full nonbonded+PME+bonded
  # stack on CPU threads fine — just slower than GPU offload.
  COMMON_RUN="-v -nt $NT -cpt $MD_CPT_MIN"
  EM_RUN="-v -nt $NT -cpt $MD_CPT_MIN"
fi

# ── 1. Topology generation ──
# AMBER ff19SB protein + OPC water (ff19SB's authors recommend OPC/OPC3,
# not TIP3P — GROMACS's own force-field docs say the same; ff19SB's CMAP
# backbone corrections weren't validated against TIP3P). Needs GROMACS
# >= 2026.3 (see header) — earlier releases don't ship ff19SB at all.
# -ignh strips existing hydrogens; pdb2gmx rebuilds them per the force
# field. Missing heavy atoms and cofactors are already handled by step 0's
# PDBFixer pre-pass, so no -missing flag is needed here.
#
# Genes needing an explicit histidine protonation state: pdb2gmx's
# automatic detector (hizzie) infers HID/HIE/HIP from each histidine's
# local H-bonding, but needs the ring atoms present to do it. A histidine
# with the whole ring unresolved in the crystal structure fails "Incomplete
# ring" even after step 0's heavy-atom repair. -his forces interactive
# selection for *every* histidine in the structure (pdb2gmx has no
# per-residue override), so only opt genes into this when they actually hit
# that error — for everything else, automatic per-residue detection is the
# better default. User-approved default for the current case
# (METTL3/HIS116, ring fully absent, not obviously catalytic): HIE for all
# of METTL3's histidines.
FORCE_HIE_GENES=" METTL3 "
PDB2GMX_STDIN="1"
if [[ "$FORCE_HIE_GENES" == *" $GENE "* ]]; then
  N_HIS=$(awk '$1=="ATOM" && $4=="HIS" {print substr($0,22,1) substr($0,23,4)}' "$PDB_FOR_GMX" | sort -u | wc -l)
  PDB2GMX_HIS_FLAG="-his"
  PDB2GMX_STDIN="$(printf '1\n%.0s' $(seq 1 $((N_HIS + 1))))"
else
  PDB2GMX_HIS_FLAG=""
fi
if [ -f protein.gro ]; then
  echo "[md] step 1: pdb2gmx — already done, skipping" | tee -a "$LOG"
else
  echo "[md] step 1: pdb2gmx" | tee -a "$LOG"
  echo "$PDB2GMX_STDIN" | gmx pdb2gmx -f "$PDB_FOR_GMX" -o protein.gro -p topol.top \
    -i posre.itp -ff amber19sb -water opc -ignh $PDB2GMX_HIS_FLAG \
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
# OPC is a 4-point model (OW/HW1/HW2/MW) like TIP4P, not 3-point like SPC —
# the solvent box coordinate file's per-molecule atom count has to match,
# so this must be tip4p.gro, not the default spc216.gro (3-point).
if [ -f protein_solv.gro ]; then
  echo "[md] step 3: solvate — already done, skipping" | tee -a "$LOG"
else
  echo "[md] step 3: solvate" | tee -a "$LOG"
  gmx solvate -cp protein_box.gro -cs tip4p.gro -o protein_solv.gro -p topol.top \
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
  gmx mdrun -deffnm em $EM_RUN \
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
# Segmented via periodic -maxh restarts + checkpoint resume, for two
# independent reasons: it's what makes a Colab-session crash (or any other
# unattended-cloud interruption) resumable instead of restarting the 50 ns
# run from t=0, and it works around an in-process memory leak observed in
# the GROMACS 2026.3 build that scales with simulation progress — three
# independent OOM kills on a 235GB-RAM droplet (PIN1_bound, GAPDH/apo,
# CBLL1/apo) all died late in their respective runs (83-91% of steps)
# regardless of system size (27k-609k atoms), consistent with a leak rather
# than genuinely needing that much memory. A fresh mdrun process starts
# with a clean heap, so capping each process to MD_PROD_MAXH_HOURS and
# resuming from checkpoint avoids both problems without changing the
# resulting trajectory at all.
MD_PROD_MAXH_HOURS="${MD_PROD_MAXH_HOURS:-2}"
echo "[md] step 8: production 50 ns (segmented, ${MD_PROD_MAXH_HOURS}h per mdrun process)" | tee -a "$LOG"
if [ ! -f md.tpr ]; then
  gmx grompp -f "$MDP_DIR/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2 \
    >> "$LOG" 2>&1
fi

PROD_MAX_SEGMENTS=50
seg=0
while [ ! -f md.gro ]; do
  seg=$((seg + 1))
  if [ "$seg" -gt "$PROD_MAX_SEGMENTS" ]; then
    echo "[md] production aborted: exceeded $PROD_MAX_SEGMENTS restart segments without reaching md.gro" | tee -a "$LOG"
    exit 1
  fi
  if [ -f md.cpt ]; then
    echo "[md] production segment $seg (resuming from checkpoint)" | tee -a "$LOG"
    gmx mdrun -deffnm md -cpi md.cpt -append -maxh "$MD_PROD_MAXH_HOURS" $COMMON_RUN \
      >> "$LOG" 2>&1 || true
  else
    echo "[md] production segment $seg (fresh start)" | tee -a "$LOG"
    gmx mdrun -deffnm md -maxh "$MD_PROD_MAXH_HOURS" $COMMON_RUN \
      >> "$LOG" 2>&1 || true
  fi
  if [ ! -f md.gro ] && [ ! -f md.cpt ]; then
    echo "[md] production segment $seg failed before writing any checkpoint — aborting" | tee -a "$LOG"
    exit 1
  fi
done

if [ -f md.gro ]; then
  echo "[md] DONE $(date -Is)  $GENE/$MODE" | tee -a "$LOG"
  ls -lh md.* | tee -a "$LOG"
else
  echo "[md] production stopped before completion (maxh limit or interrupted) — re-run this script to resume from md.cpt" | tee -a "$LOG"
fi
