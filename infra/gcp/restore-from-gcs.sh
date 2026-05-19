#!/usr/bin/env bash
# Restore pipeline outputs from a GCS bucket onto a fresh VM.
# Run on the *new* VM after cloud-init finishes.
#
# Usage:
#   bash infra/gcp/restore-from-gcs.sh                 # latest backup
#   ARCHIVE=thesis-backup-20260518_153022.tar.gz \
#     bash infra/gcp/restore-from-gcs.sh               # specific timestamp

set -euo pipefail

WORK_DIR="${WORK_DIR:-/opt/Thesis-insilico}"
BUCKET="${BUCKET:-thesis-as-m6a-data}"
ARCHIVE="${ARCHIVE:-latest.tar.gz}"
PROJECT="$(gcloud config get-value project 2>/dev/null || true)"

if [ -z "$PROJECT" ] || [ "$PROJECT" = "(unset)" ]; then
  echo "No active project. Run:  gcloud config set project PROJECT_ID"
  exit 1
fi

FULL_BUCKET="${PROJECT}-${BUCKET}"
echo "Project: $PROJECT"
echo "Bucket:  gs://$FULL_BUCKET"
echo "Archive: $ARCHIVE"
echo

if ! gsutil ls "gs://$FULL_BUCKET/$ARCHIVE" >/dev/null 2>&1; then
  echo "Archive not found. Available backups:"
  gsutil ls "gs://$FULL_BUCKET/" || true
  exit 1
fi

tmp="/tmp/$(basename "$ARCHIVE")"
echo "Downloading gs://$FULL_BUCKET/$ARCHIVE -> $tmp ..."
gsutil cp "gs://$FULL_BUCKET/$ARCHIVE" "$tmp"

cd "$WORK_DIR"
echo "Extracting into $WORK_DIR ..."
tar xzf "$tmp"
rm -f "$tmp"

echo
echo "Restored:"
ls -la results/ 2>/dev/null | head -15 || true
echo
echo "You can now re-run the pipeline from where you left off, or just"
echo "inspect the results. To continue from the ranking stage:"
echo "  source /opt/miniforge/etc/profile.d/conda.sh && conda activate arsenic-m6a"
echo "  python scripts/07_score_and_rank.py"
