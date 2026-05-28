#!/usr/bin/env bash
# One-shot bootstrap for a fresh RunPod GPU pod.
# Installs Miniforge + GROMACS (GPU) + analysis dependencies, then clones the
# repo. Designed to be pasted directly into the RunPod web terminal.

set -euo pipefail

apt-get update
apt-get install -y git tmux htop rsync build-essential curl

# Miniforge
if [ ! -d /opt/miniforge ]; then
  curl -L -o /tmp/miniforge.sh \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
  bash /tmp/miniforge.sh -b -p /opt/miniforge
  rm /tmp/miniforge.sh
fi
source /opt/miniforge/etc/profile.d/conda.sh

# MD environment
if ! conda env list | grep -q '^md'; then
  mamba create -n md -c conda-forge -y \
    gromacs=2024 \
    python=3.11 \
    mdanalysis \
    freesasa \
    biopython \
    matplotlib \
    pandas \
    numpy \
    scipy \
    tqdm \
    pyyaml
fi
conda activate md

# Repo
if [ ! -d /opt/Thesis-insilico ]; then
  git clone https://github.com/Wachirawut2023/Thesis-insilico.git /opt/Thesis-insilico
fi
cd /opt/Thesis-insilico
git checkout claude/arsenic-m6a-inhibition-model-6pwWS
git pull

chmod +x scripts/md/run_md.sh scripts/md/run_all.sh

# Verification
echo
echo "=== bootstrap done ==="
gmx --version | grep -E "GROMACS|GPU support|CUDA" || true
nvidia-smi --query-gpu=name,memory.total --format=csv
echo
echo "Next:"
echo "  1) Pull data/prepared/ from your GCS bucket OR rerun static stages 1-2"
echo "  2) tmux new -s md"
echo "  3) bash scripts/md/run_all.sh"
