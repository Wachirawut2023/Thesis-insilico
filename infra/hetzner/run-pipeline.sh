#!/usr/bin/env bash
# Run on the Hetzner server (inside tmux). Executes the full pipeline with logging.
set -euo pipefail

cd /root/Thesis-insilico
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
