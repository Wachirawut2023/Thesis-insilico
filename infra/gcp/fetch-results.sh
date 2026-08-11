#!/usr/bin/env bash
# Pull results from the GCP instance to ./results-pulled/<timestamp>/ on your laptop.
# Uses gcloud compute scp so no manual SSH key wrangling is needed.
#
# Usage:
#   ./fetch-results.sh
#   INSTANCE_NAME=as-m6a-2 ZONE=us-east1-b ./fetch-results.sh

set -euo pipefail

INSTANCE_NAME="${INSTANCE_NAME:-as-m6a}"
ZONE="${ZONE:-us-central1-a}"
WORK_DIR="${WORK_DIR:-/opt/Thesis-insilico}"

ts="$(date +%Y%m%d_%H%M%S)"
dest="results-pulled/$ts"
mkdir -p "$dest"

# Use rsync over gcloud's SSH wrapper for incremental, resumable transfer.
gcloud compute scp --recurse --zone="$ZONE" \
  "$INSTANCE_NAME:$WORK_DIR/results/" "$dest/"

# Drop the heavy docking_poses tree to keep the local copy small.
rm -rf "$dest/results/docking_poses" 2>/dev/null || true

echo
echo "Pulled to $dest"
find "$dest" -maxdepth 3 -type f | head -30
