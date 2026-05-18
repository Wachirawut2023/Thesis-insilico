#!/usr/bin/env bash
# Provision a GCP Compute Engine instance for the As(III)-m6A pipeline (on-demand).
#
# Prereqs:
#   1. GCP account: https://console.cloud.google.com (new users get $300 / 90-day free trial)
#   2. Create or select a Project. Note the PROJECT_ID.
#   3. Install gcloud CLI:
#        macOS:  brew install --cask google-cloud-sdk
#        Linux:  https://cloud.google.com/sdk/docs/install-sdk
#   4. Authenticate + select project:
#        gcloud auth login
#        gcloud config set project PROJECT_ID
#        gcloud auth application-default login   # optional, for SDKs later
#   5. Enable the Compute Engine API (one-time):
#        gcloud services enable compute.googleapis.com
#      (Billing must be enabled on the project — happens automatically with free trial.)
#
# Usage:
#   ./provision.sh                                    # n2d-standard-8 in us-central1-a
#   MACHINE_TYPE=c3-standard-8 ZONE=us-east1-b ./provision.sh
#
# Tear down (STOPS BILLING):
#   gcloud compute instances delete as-m6a --zone=us-central1-a

set -euo pipefail

INSTANCE_NAME="${INSTANCE_NAME:-as-m6a}"
MACHINE_TYPE="${MACHINE_TYPE:-n2d-standard-8}"   # 8 AMD EPYC vCPU, 32 GB RAM ($0.388/hr us-central1)
ZONE="${ZONE:-us-central1-a}"                    # cheapest region; alt: us-east1-b, europe-west4-a
IMAGE_FAMILY="${IMAGE_FAMILY:-ubuntu-2404-lts-amd64}"
IMAGE_PROJECT="${IMAGE_PROJECT:-ubuntu-os-cloud}"
DISK_SIZE_GB="${DISK_SIZE_GB:-50}"
DISK_TYPE="${DISK_TYPE:-pd-balanced}"            # cheaper than pd-ssd; faster than pd-standard
REPO_URL="${REPO_URL:-https://github.com/Wachirawut2023/Thesis-insilico.git}"
BRANCH="${BRANCH:-claude/arsenic-m6a-inhibition-model-6pwWS}"

cd "$(dirname "$0")"

# 0. Sanity checks
if ! command -v gcloud >/dev/null 2>&1; then
  echo "gcloud CLI not found. Install: https://cloud.google.com/sdk/docs/install-sdk"
  exit 1
fi
PROJECT="$(gcloud config get-value project 2>/dev/null || true)"
if [ -z "$PROJECT" ] || [ "$PROJECT" = "(unset)" ]; then
  echo "No active project. Run:  gcloud config set project PROJECT_ID"
  exit 1
fi
echo "Project: $PROJECT"
echo "Zone:    $ZONE"
echo "Type:    $MACHINE_TYPE"
echo

# 1. Ensure Compute Engine API is enabled
if ! gcloud services list --enabled --filter='config.name:compute.googleapis.com' --format='value(config.name)' | grep -q compute; then
  echo "Enabling Compute Engine API ..."
  gcloud services enable compute.googleapis.com
fi

# 2. Render cloud-init with REPO_URL / BRANCH baked in (so cloud-init picks the right branch)
tmp_user_data="$(mktemp)"
trap 'rm -f "$tmp_user_data"' EXIT
sed -e "s|\${REPO_URL:-https://github.com/Wachirawut2023/Thesis-insilico.git}|$REPO_URL|g" \
    -e "s|\${BRANCH:-claude/arsenic-m6a-inhibition-model-6pwWS}|$BRANCH|g" \
    cloud-init.yaml > "$tmp_user_data"

# 3. Create or reuse the instance
if gcloud compute instances describe "$INSTANCE_NAME" --zone="$ZONE" >/dev/null 2>&1; then
  echo "Instance $INSTANCE_NAME already exists in $ZONE; reusing."
else
  echo "Creating $INSTANCE_NAME ($MACHINE_TYPE in $ZONE) ..."
  gcloud compute instances create "$INSTANCE_NAME" \
    --zone="$ZONE" \
    --machine-type="$MACHINE_TYPE" \
    --image-family="$IMAGE_FAMILY" \
    --image-project="$IMAGE_PROJECT" \
    --boot-disk-size="${DISK_SIZE_GB}GB" \
    --boot-disk-type="$DISK_TYPE" \
    --metadata-from-file=user-data="$tmp_user_data" \
    --tags=as-m6a \
    --scopes=cloud-platform
fi

ip="$(gcloud compute instances describe "$INSTANCE_NAME" --zone="$ZONE" --format='value(networkInterfaces[0].accessConfigs[0].natIP)')"

echo
echo "Instance ready at $ip (zone $ZONE)"
echo
echo "Watch cloud-init progress (~5-10 min):"
echo "  gcloud compute ssh $INSTANCE_NAME --zone=$ZONE -- tail -f /var/log/bootstrap.log"
echo
echo "Once you see '[bootstrap] done':"
echo "  gcloud compute ssh $INSTANCE_NAME --zone=$ZONE"
echo "  sudo -i"
echo "  cd /opt/Thesis-insilico && tmux new -s pipeline"
echo "  bash infra/gcp/run-pipeline.sh"
echo
echo "Pull results back from your laptop:"
echo "  ./fetch-results.sh"
echo
echo "Tear down when done (STOP BILLING):"
echo "  gcloud compute instances delete $INSTANCE_NAME --zone=$ZONE --quiet"
