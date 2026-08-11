#!/usr/bin/env bash
# Pull MD results from the GPU instance and back them up to GCS.
# Run from Cloud Shell or local gcloud.

set -euo pipefail

INSTANCE_NAME="${INSTANCE_NAME:-as-m6a-gpu}"
ZONE="${ZONE:-us-central1-a}"
WORK_DIR="${WORK_DIR:-/opt/Thesis-insilico}"

PROJECT="$(gcloud config get-value project 2>/dev/null || true)"
if [ -z "$PROJECT" ] || [ "$PROJECT" = "(unset)" ]; then
  echo "No active project. Run:  gcloud config set project PROJECT_ID"
  exit 1
fi

ts="$(date +%Y%m%d_%H%M%S)"
local_dir="md-results-$ts"
mkdir -p "$local_dir"

echo "1) Tarball on the GPU instance ..."
gcloud compute ssh "$INSTANCE_NAME" --zone="$ZONE" -- \
  "cd $WORK_DIR && sudo tar czf /tmp/md_results.tar.gz results/md/ && ls -lh /tmp/md_results.tar.gz"

echo "2) SCP to here ..."
gcloud compute scp --zone="$ZONE" \
  "$INSTANCE_NAME:/tmp/md_results.tar.gz" "$local_dir/"

echo "3) Push to GCS for permanence ..."
gsutil cp "$local_dir/md_results.tar.gz" \
  "gs://${PROJECT}-thesis-as-m6a-data/md_results-${ts}.tar.gz"
gsutil cp "$local_dir/md_results.tar.gz" \
  "gs://${PROJECT}-thesis-as-m6a-data/md_results_latest.tar.gz"

echo
echo "Pulled to:  $local_dir/md_results.tar.gz"
echo "Backed up:  gs://${PROJECT}-thesis-as-m6a-data/md_results_latest.tar.gz"
echo
echo "Download to your iPad (Cloud Shell):"
echo "  cloudshell download $local_dir/md_results.tar.gz"
echo
echo "Now STOP THE GPU INSTANCE — it bills \$0.82/hr while running:"
echo "  gcloud compute instances delete $INSTANCE_NAME --zone=$ZONE --quiet"
