#!/usr/bin/env bash
# Pull MD results from the AMD GPU droplet back to your local machine.

set -euo pipefail

DROPLET_NAME="${DROPLET_NAME:-as-m6a-amd-gpu}"
WORK_DIR="${WORK_DIR:-/opt/Thesis-insilico}"

if ! command -v doctl >/dev/null 2>&1; then
  echo "doctl not in PATH. Install: https://docs.digitalocean.com/reference/doctl/how-to/install/"
  exit 1
fi

ip="$(doctl compute droplet get "$DROPLET_NAME" --format PublicIPv4 --no-header)"
if [ -z "$ip" ]; then
  echo "ERROR: could not resolve IP for droplet $DROPLET_NAME"
  exit 1
fi

ts="$(date +%Y%m%d_%H%M%S)"
local_dir="md-results-$ts"
mkdir -p "$local_dir"

echo "1) Tarball on the droplet ..."
ssh "root@$ip" "cd $WORK_DIR && tar czf /tmp/md_results.tar.gz results/md/ && ls -lh /tmp/md_results.tar.gz"

echo "2) SCP to here ..."
scp "root@$ip:/tmp/md_results.tar.gz" "$local_dir/"

echo
echo "Pulled to: $local_dir/md_results.tar.gz"
echo
echo "Now DELETE THE DROPLET — it bills per hour while running:"
echo "  doctl compute droplet delete $DROPLET_NAME"
