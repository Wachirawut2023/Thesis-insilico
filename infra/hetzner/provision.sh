#!/usr/bin/env bash
# Provision a Hetzner Cloud server for the As(III)-m6A pipeline.
#
# Prereqs:
#   1. Hetzner Cloud account: https://console.hetzner.cloud
#   2. Create a Project, then a "Read & Write" API token under
#      Project -> Security -> API Tokens. Save the token; you only see it once.
#   3. Install hcloud CLI:
#        macOS:  brew install hcloud
#        Linux:  curl -L https://github.com/hetznercloud/cli/releases/latest/download/hcloud-linux-amd64.tar.gz | tar -xz && sudo mv hcloud /usr/local/bin/
#   4. hcloud context create thesis  # paste token when prompted
#
# Usage:
#   ./provision.sh              # creates server "as-m6a" with CCX33
#   SERVER_TYPE=ccx43 ./provision.sh
#   SERVER_NAME=as-m6a-2 ./provision.sh
#
# To tear down at the end:
#   hcloud server delete as-m6a

set -euo pipefail

SERVER_NAME="${SERVER_NAME:-as-m6a}"
SERVER_TYPE="${SERVER_TYPE:-ccx33}"        # 8 dedicated vCPU, 32GB RAM, 240GB NVMe, EUR 0.073/h
IMAGE="${IMAGE:-ubuntu-24.04}"
LOCATION="${LOCATION:-nbg1}"               # Nuremberg, DE; alt: fsn1, hel1, ash, hil, sin
SSH_KEY_NAME="${SSH_KEY_NAME:-thesis}"
REPO_URL="${REPO_URL:-https://github.com/Wachirawut2023/Thesis-insilico.git}"
BRANCH="${BRANCH:-claude/arsenic-m6a-inhibition-model-6pwWS}"

cd "$(dirname "$0")"

# 1. SSH key
if ! hcloud ssh-key describe "$SSH_KEY_NAME" >/dev/null 2>&1; then
  if [ ! -f ~/.ssh/id_ed25519.pub ] && [ ! -f ~/.ssh/id_rsa.pub ]; then
    echo "No SSH public key found at ~/.ssh/id_ed25519.pub or ~/.ssh/id_rsa.pub."
    echo "Generate one with:  ssh-keygen -t ed25519 -C thesis"
    exit 1
  fi
  pubkey="$([ -f ~/.ssh/id_ed25519.pub ] && cat ~/.ssh/id_ed25519.pub || cat ~/.ssh/id_rsa.pub)"
  hcloud ssh-key create --name "$SSH_KEY_NAME" --public-key "$pubkey"
fi

# 2. Render user-data with REPO_URL / BRANCH baked in
tmp_user_data="$(mktemp)"
trap 'rm -f "$tmp_user_data"' EXIT
sed -e "s|\${REPO_URL:-https://github.com/Wachirawut2023/Thesis-insilico.git}|$REPO_URL|g" \
    -e "s|\${BRANCH:-claude/arsenic-m6a-inhibition-model-6pwWS}|$BRANCH|g" \
    cloud-init.yaml > "$tmp_user_data"

# 3. Create or reuse the server
if hcloud server describe "$SERVER_NAME" >/dev/null 2>&1; then
  echo "Server $SERVER_NAME already exists; reusing."
else
  echo "Creating $SERVER_NAME ($SERVER_TYPE in $LOCATION) ..."
  hcloud server create \
    --name "$SERVER_NAME" \
    --type "$SERVER_TYPE" \
    --image "$IMAGE" \
    --location "$LOCATION" \
    --ssh-key "$SSH_KEY_NAME" \
    --user-data-from-file "$tmp_user_data" \
    --start-after-create
fi

ip="$(hcloud server ip "$SERVER_NAME")"
echo
echo "Server ready at $ip"
echo
echo "Next steps:"
echo "  ssh root@$ip"
echo "  tail -f /var/log/bootstrap.log    # watch cloud-init progress (~5-10 min)"
echo "  ls /root/.bootstrap-complete       # appears when bootstrap finishes"
echo
echo "  cd /root/Thesis-insilico"
echo "  conda activate arsenic-m6a"
echo "  tmux new -s pipeline"
echo "  make smoke                         # 3-target smoke test first"
echo "  make all                           # full panel"
echo
echo "Pull results back from your laptop:"
echo "  ./fetch-results.sh $ip"
echo
echo "Tear down when done (stops billing):"
echo "  hcloud server delete $SERVER_NAME"
