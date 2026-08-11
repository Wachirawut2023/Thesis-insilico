#!/usr/bin/env bash
# Run on the GCP instance (inside tmux, as root). Executes the full pipeline with logging.
# Identical to the Hetzner runner except for the WORK_DIR path.
set -euo pipefail

WORK_DIR="${WORK_DIR:-/opt/Thesis-insilico}"
cd "$WORK_DIR"
source /opt/miniforge/etc/profile.d/conda.sh
conda activate arsenic-m6a

LOG_DIR=results/logs
mkdir -p "$LOG_DIR"
ts="$(date +%Y%m%d_%H%M%S)"
log="$LOG_DIR/run_${ts}.log"

echo "Pipeline start $(date -Is)" | tee "$log"
{
  make ligand
  make fetch
  make prep
  make cys
  make pockets
  make dock-nc
  make dock-cov
  make rank
} 2>&1 | tee -a "$log"
echo "Pipeline done $(date -Is)" | tee -a "$log"

echo
echo "Key outputs:"
ls -la results/*.tsv results/REPORT.md results/figures/*.png 2>/dev/null || true
