#!/usr/bin/env bash
# Backup the pipeline outputs to a GCS bucket so they survive VM lifecycle.
# Run on the *old* VM before deleting it.
#
# Usage:
#   bash infra/gcp/backup-to-gcs.sh
#   BUCKET=my-bucket bash infra/gcp/backup-to-gcs.sh
#
# Cost: GCS standard storage is ~$0.02/GB-month. Typical pipeline output
# is <200 MB so monthly cost is rounding error.

set -euo pipefail

WORK_DIR="${WORK_DIR:-/opt/Thesis-insilico}"
BUCKET="${BUCKET:-thesis-as-m6a-data}"
LOCATION="${LOCATION:-us-central1}"
PROJECT="$(gcloud config get-value project 2>/dev/null || true)"

if [ -z "$PROJECT" ] || [ "$PROJECT" = "(unset)" ]; then
  echo "No active project. Run:  gcloud config set project PROJECT_ID"
  exit 1
fi

# Bucket names are globally unique — prefix with project ID so they're stable.
FULL_BUCKET="${PROJECT}-${BUCKET}"

echo "Project: $PROJECT"
echo "Bucket:  gs://$FULL_BUCKET"
echo

# Create bucket if it doesn't exist (idempotent).
if ! gsutil ls -b "gs://$FULL_BUCKET" >/dev/null 2>&1; then
  echo "Creating bucket gs://$FULL_BUCKET in $LOCATION ..."
  gsutil mb -l "$LOCATION" "gs://$FULL_BUCKET"
fi

cd "$WORK_DIR"

ts="$(date +%Y%m%d_%H%M%S)"
archive="/tmp/thesis-backup-${ts}.tar.gz"

echo "Packing results, data, and ligand artifacts ..."
tar czf "$archive" \
  results \
  data/prepared \
  ligands/as3.pdbqt \
  --warning=no-file-changed 2>/dev/null || true

size=$(du -h "$archive" | awk '{print $1}')
echo "Uploading $archive ($size) to gs://$FULL_BUCKET/ ..."
gsutil cp "$archive" "gs://$FULL_BUCKET/"

# Also push as "latest.tar.gz" for restore-from-gcs.sh to grab without
# needing to know the timestamp.
gsutil cp "$archive" "gs://$FULL_BUCKET/latest.tar.gz"

rm -f "$archive"

echo
echo "Backed up to:"
echo "  gs://$FULL_BUCKET/thesis-backup-${ts}.tar.gz   (timestamped)"
echo "  gs://$FULL_BUCKET/latest.tar.gz                (alias for restore)"
echo
echo "Now safe to delete the VM:"
echo "  gcloud compute instances delete as-m6a --zone=us-central1-a --quiet"
