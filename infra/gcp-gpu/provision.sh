#!/usr/bin/env bash
# Provision a GCP GPU instance for Tier-3 MD work.
#
# Defaults to g2-standard-8 with NVIDIA L4 in us-central1-a.
# Uses the Deep Learning VM image (CUDA 12.x + drivers pre-installed) so
# bootstrap only needs to install GROMACS via conda.
#
# Prereqs:
#   1. Same GCP account/project you used for the static pipeline.
#   2. GPU quota in the chosen region. Most projects have 0 GPU quota by
#      default. Request an increase:
#         Console -> IAM & Admin -> Quotas -> filter "NVIDIA L4 GPUs"
#         in us-central1 -> Edit Quotas -> request 1 (or more)
#      Approval typically takes minutes to a few hours.
#   3. Compute Engine API enabled (was done by the CPU provisioning).
#
# Usage:
#   ./provision.sh
#   INSTANCE_NAME=as-m6a-gpu-2 ZONE=us-east1-c ./provision.sh
#
# Tear down (STOPS BILLING):
#   gcloud compute instances delete as-m6a-gpu --zone=us-central1-a --quiet

set -euo pipefail

INSTANCE_NAME="${INSTANCE_NAME:-as-m6a-gpu}"
MACHINE_TYPE="${MACHINE_TYPE:-g2-standard-8}"    # 8 vCPU, 32 GB, 1× L4
GPU_TYPE="${GPU_TYPE:-nvidia-l4}"
GPU_COUNT="${GPU_COUNT:-1}"
ZONE="${ZONE:-us-central1-a}"
IMAGE_FAMILY="${IMAGE_FAMILY:-common-cu123-debian-11}"
IMAGE_PROJECT="${IMAGE_PROJECT:-deeplearning-platform-release}"
DISK_SIZE_GB="${DISK_SIZE_GB:-100}"
DISK_TYPE="${DISK_TYPE:-pd-balanced}"
REPO_URL="${REPO_URL:-https://github.com/Wachirawut2023/Thesis-insilico.git}"
BRANCH="${BRANCH:-claude/arsenic-m6a-inhibition-model-6pwWS}"

cd "$(dirname "$0")"

if ! command -v gcloud >/dev/null 2>&1; then
  echo "gcloud not in PATH. Install: https://cloud.google.com/sdk/docs/install-sdk"
  exit 1
fi
PROJECT="$(gcloud config get-value project 2>/dev/null || true)"
if [ -z "$PROJECT" ] || [ "$PROJECT" = "(unset)" ]; then
  echo "No active project. Run:  gcloud config set project PROJECT_ID"
  exit 1
fi
echo "Project: $PROJECT"
echo "Zone:    $ZONE"
echo "Type:    $MACHINE_TYPE + ${GPU_COUNT}× ${GPU_TYPE}"
echo

# Render cloud-init with REPO_URL / BRANCH baked in
tmp_user_data="$(mktemp)"
trap 'rm -f "$tmp_user_data"' EXIT
sed -e "s|\${REPO_URL:-https://github.com/Wachirawut2023/Thesis-insilico.git}|$REPO_URL|g" \
    -e "s|\${BRANCH:-claude/arsenic-m6a-inhibition-model-6pwWS}|$BRANCH|g" \
    cloud-init.yaml > "$tmp_user_data"

if gcloud compute instances describe "$INSTANCE_NAME" --zone="$ZONE" >/dev/null 2>&1; then
  echo "Instance $INSTANCE_NAME exists in $ZONE; reusing."
else
  echo "Creating $INSTANCE_NAME ($MACHINE_TYPE with ${GPU_COUNT}× ${GPU_TYPE}) ..."
  gcloud compute instances create "$INSTANCE_NAME" \
    --zone="$ZONE" \
    --machine-type="$MACHINE_TYPE" \
    --accelerator="type=$GPU_TYPE,count=$GPU_COUNT" \
    --image-family="$IMAGE_FAMILY" \
    --image-project="$IMAGE_PROJECT" \
    --boot-disk-size="${DISK_SIZE_GB}GB" \
    --boot-disk-type="$DISK_TYPE" \
    --metadata="install-nvidia-driver=True" \
    --metadata-from-file=user-data="$tmp_user_data" \
    --maintenance-policy=TERMINATE \
    --restart-on-failure \
    --tags=as-m6a-gpu \
    --scopes=cloud-platform
fi

ip="$(gcloud compute instances describe "$INSTANCE_NAME" --zone="$ZONE" --format='value(networkInterfaces[0].accessConfigs[0].natIP)')"

echo
echo "GPU instance ready at $ip (zone $ZONE)"
echo
echo "Cloud-init takes ~10-15 min (NVIDIA driver init + conda env build)."
echo "Watch progress:"
echo "  gcloud compute ssh $INSTANCE_NAME --zone=$ZONE -- tail -f /var/log/bootstrap.log"
echo
echo "Then connect:"
echo "  gcloud compute ssh $INSTANCE_NAME --zone=$ZONE"
echo "  sudo -i && nvidia-smi"
echo
echo "Run MD (after pulling data/prepared/ from GCS):"
echo "  cd /opt/Thesis-insilico"
echo "  source /opt/miniforge/etc/profile.d/conda.sh && conda activate md"
echo "  tmux new -s md"
echo "  bash scripts/md/run_all.sh"
echo
echo "Tear down (STOP BILLING — \$0.82/hr while running):"
echo "  gcloud compute instances delete $INSTANCE_NAME --zone=$ZONE --quiet"
