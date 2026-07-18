#!/usr/bin/env bash
# Orchestrate Tier-3 MD across all top-8 hits, both modes.
#
# Reads scripts/md/data/md_targets.tsv and runs run_md.sh + analyze.py per row.
# Logs per-protein wall time; resumable at two levels:
#   - whole gene/mode runs are skipped once md.gro exists
#   - a gene/mode run interrupted mid-production (crash, disconnect, or an
#     MD_MAXH time limit inside run_md.sh) resumes from its GROMACS
#     checkpoint (md.cpt) instead of restarting the 50 ns trajectory
# So on Colab, just re-run this script (e.g. from a fresh notebook cell)
# after any crash/disconnect — it picks up exactly where it stopped.
#
# Usage:   bash run_all.sh
#          GENES_FILTER="TXN1,PIN1" bash run_all.sh     # subset
#          MODES_FILTER="apo" bash run_all.sh           # apo only
#          MD_MAXH=5 bash run_all.sh                    # cap each run's production step at 5h wall time
#          MD_GPU=0 bash run_all.sh                      # CPU-only (GPU backend can't see a device)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGETS_TSV="$SCRIPT_DIR/data/md_targets.tsv"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
COV_TSV="$REPO_DIR/results/docking_covalent.tsv"
OUT_BASE="${OUT_BASE:-$REPO_DIR/results/md}"

GENES_FILTER="${GENES_FILTER:-}"
MODES_FILTER="${MODES_FILTER:-}"

mkdir -p "$OUT_BASE"
overall_log="$OUT_BASE/run_all.log"
echo "=== Tier-3 MD orchestration $(date -Is) ===" | tee -a "$overall_log"

resolve_anchor() {
  # If anchor_resi is "-" or empty, look up the best-covalent-score anchor
  # for this gene from results/docking_covalent.tsv.
  local gene="$1" anchor_in="$2" partner_in="$3"
  if [ "$anchor_in" != "-" ] && [ -n "$anchor_in" ]; then
    echo "$anchor_in $partner_in"
    return
  fi
  if [ ! -f "$COV_TSV" ]; then
    echo "-" "-"
    return
  fi
  # Pick the row with highest covalent_score for this gene
  local row
  row=$(awk -F'\t' -v g="$gene" 'NR>1 && $1==g {print}' "$COV_TSV" \
        | sort -t$'\t' -k11,11 -gr | head -1)
  if [ -z "$row" ]; then
    echo "-" "-"
    return
  fi
  local resi partner
  resi=$(echo "$row" | awk -F'\t' '{print $3}')
  partner=$(echo "$row" | awk -F'\t' '{print $8}' | awk -F'[;:]' '{print $1}' | sed 's/[A-Za-z]//g')
  echo "$resi" "${partner:--}"
}

# Read TSV (skip header)
tail -n +2 "$TARGETS_TSV" | while IFS=$'\t' read -r gene mode chain anchor_in partner_in binding_mode source_pdb notes; do
  if [ -n "$GENES_FILTER" ]; then
    [[ ",$GENES_FILTER," == *",$gene,"* ]] || continue
  fi
  if [ -n "$MODES_FILTER" ]; then
    [[ ",$MODES_FILTER," == *",$mode,"* ]] || continue
  fi

  run_dir="$OUT_BASE/${gene}_${mode}"
  if [ -f "$run_dir/md.gro" ]; then
    echo "[skip] $gene/$mode — already completed" | tee -a "$overall_log"
    continue
  fi

  echo "" | tee -a "$overall_log"
  echo "=== $gene / $mode  (started $(date -Is)) ===" | tee -a "$overall_log"
  start=$(date +%s)

  if [ "$mode" = "bound" ]; then
    read resi partner <<< "$(resolve_anchor "$gene" "$anchor_in" "$partner_in")"
    echo "[orchestrator] $gene bound anchor: chain=$chain resi=$resi partner=$partner" | tee -a "$overall_log"
    if [ "$resi" = "-" ] || [ -z "$resi" ]; then
      echo "[skip] $gene/$mode — no covalent anchor available" | tee -a "$overall_log"
      continue
    fi
    bash "$SCRIPT_DIR/run_md.sh" "$gene" bound "$chain" "$resi" "$partner" \
        2>&1 | tee -a "$overall_log" || {
      echo "[FAIL] $gene/$mode — see $run_dir/run.log" | tee -a "$overall_log"
      continue
    }
  else
    bash "$SCRIPT_DIR/run_md.sh" "$gene" apo \
        2>&1 | tee -a "$overall_log" || {
      echo "[FAIL] $gene/$mode — see $run_dir/run.log" | tee -a "$overall_log"
      continue
    }
  fi

  end=$(date +%s)
  wall_min=$(( (end - start) / 60 ))

  if [ ! -f "$run_dir/md.gro" ]; then
    # run_md.sh stopped early (MD_MAXH limit, disconnect, or crash) without
    # finishing production. It checkpointed via md.cpt, so re-running this
    # script later will resume this exact target where it left off.
    echo "[partial] $gene/$mode — production not finished after ${wall_min} min, will resume next run" | tee -a "$overall_log"
    continue
  fi

  echo "[done] $gene/$mode in ${wall_min} min" | tee -a "$overall_log"

  # Analyse immediately
  python3 "$SCRIPT_DIR/analyze.py" --gene "$gene" --mode "$mode" --md-root "$OUT_BASE" \
      2>&1 | tee -a "$overall_log" || true
done

echo "" | tee -a "$overall_log"
echo "=== Tier-3 MD complete $(date -Is) ===" | tee -a "$overall_log"
echo "Tarball results: tar czf md_results.tar.gz -C $(dirname $OUT_BASE) md/" | tee -a "$overall_log"
