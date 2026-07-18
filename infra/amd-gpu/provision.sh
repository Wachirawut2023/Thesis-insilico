#!/usr/bin/env bash
# Provision a DigitalOcean AMD GPU Droplet (MI300X) for Tier-3 MD work.
#
# Prereqs:
#   1. `doctl` installed and authenticated: https://docs.digitalocean.com/reference/doctl/
#        doctl auth init
#   2. An SSH key already added to your DO account (doctl compute ssh-key list).
#   3. GPU Droplets are region-limited — check availability with:
#        doctl compute size list | grep -i gpu
#      (typically NYC2/TOR1/ATL1; region names shift as DO expands capacity).
#
# Usage:
#   SSH_KEY_NAME=my-key ./provision.sh
#   DROPLET_NAME=as-m6a-amd-gpu REGION=nyc2 ./provision.sh
#
# Tear down (STOPS BILLING):
#   doctl compute droplet delete as-m6a-amd-gpu

set -euo pipefail

DROPLET_NAME="${DROPLET_NAME:-as-m6a-amd-gpu}"
REGION="${REGION:-nyc2}"
SIZE="${SIZE:-gpu-mi300x1-192gb}"   # 1x MI300X GPU Droplet slug (verify with `doctl compute size list`)
IMAGE="${IMAGE:-gpu-amd-base}"      # AI/ML-ready AMD GPU image slug (verify with `doctl compute image list-distribution --public | grep -i amd`)
REPO_URL="${REPO_URL:-https://github.com/Wachirawut2023/Thesis-insilico.git}"
BRANCH="${BRANCH:-claude/tier3-md-amd-gpu-oydfoj}"
SSH_KEY_NAME="${SSH_KEY_NAME:?set SSH_KEY_NAME to the name/ID of an SSH key already registered with doctl (doctl compute ssh-key list)}"

cd "$(dirname "$0")"

if ! command -v doctl >/dev/null 2>&1; then
  echo "doctl not in PATH. Install: https://docs.digitalocean.com/reference/doctl/how-to/install/"
  exit 1
fi
doctl account get >/dev/null || { echo "doctl not authenticated. Run: doctl auth init"; exit 1; }

SSH_KEY_ID="$(doctl compute ssh-key list --format ID,Name --no-header | awk -v n="$SSH_KEY_NAME" '$2==n{print $1}')"
if [ -z "$SSH_KEY_ID" ]; then
  echo "ERROR: SSH key '$SSH_KEY_NAME' not found in your DO account."
  echo "List available keys: doctl compute ssh-key list"
  exit 1
fi

echo "Region:  $REGION"
echo "Size:    $SIZE"
echo "Image:   $IMAGE"
echo

tmp_user_data="$(mktemp)"
trap 'rm -f "$tmp_user_data"' EXIT
sed -e "s|\${REPO_URL:-https://github.com/Wachirawut2023/Thesis-insilico.git}|$REPO_URL|g" \
    -e "s|\${BRANCH:-claude/tier3-md-amd-gpu-oydfoj}|$BRANCH|g" \
    cloud-init.yaml > "$tmp_user_data"

if doctl compute droplet get "$DROPLET_NAME" >/dev/null 2>&1; then
  echo "Droplet $DROPLET_NAME already exists; reusing."
else
  echo "Creating $DROPLET_NAME ($SIZE in $REGION) ..."
  doctl compute droplet create "$DROPLET_NAME" \
    --region "$REGION" \
    --size "$SIZE" \
    --image "$IMAGE" \
    --ssh-keys "$SSH_KEY_ID" \
    --user-data-file "$tmp_user_data" \
    --wait
fi

ip="$(doctl compute droplet get "$DROPLET_NAME" --format PublicIPv4 --no-header)"

echo
echo "GPU droplet ready at $ip"
echo
echo "Cloud-init takes ~30-45 min (ROCm verify + GROMACS built from source)."
echo "Watch progress:"
echo "  ssh root@$ip 'tail -f /var/log/bootstrap.log'"
echo
echo "Then connect:"
echo "  ssh root@$ip"
echo "  rocm-smi"
echo
echo "Run MD (after pulling data/prepared/):"
echo "  cd /opt/Thesis-insilico"
echo "  source /etc/profile.d/md-env.sh"
echo "  tmux new -s md"
echo "  bash scripts/md/run_all.sh"
echo
echo "Tear down (STOP BILLING):"
echo "  doctl compute droplet delete $DROPLET_NAME"
