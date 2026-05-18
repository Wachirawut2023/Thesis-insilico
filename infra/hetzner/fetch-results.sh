#!/usr/bin/env bash
# Pull results from the Hetzner server to ./results-pulled/<timestamp>/ on your laptop.
# Usage:  ./fetch-results.sh <server-ip>

set -euo pipefail
ip="${1:?usage: $0 <server-ip>}"
ts="$(date +%Y%m%d_%H%M%S)"
dest="results-pulled/$ts"
mkdir -p "$dest"

rsync -avz --partial \
  --include='*.tsv' \
  --include='*.md' \
  --include='figures/' \
  --include='figures/*' \
  --include='logs/' \
  --include='logs/*' \
  --exclude='docking_poses/' \
  --exclude='tmp/' \
  "root@$ip:/root/Thesis-insilico/results/" "$dest/"

echo
echo "Pulled to $dest"
ls -la "$dest"
