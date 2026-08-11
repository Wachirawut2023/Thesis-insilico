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

# MD environment. IMPORTANT: conda-forge's *plain* "nompi" gromacs build
# offloads to GPU via OpenCL, and GROMACS' OpenCL backend does not support
# Volta/Turing/Ampere-or-newer NVIDIA GPUs for compute -- gmx enumerates the
# card, reports it "incompatible", and silently drops to CPU-only with no
# error. Pin the real CUDA-enabled build (build string "nompi_cuda*")
# explicitly, since the solver won't prefer it on its own.
CUDA_VER=$(nvidia-smi 2>/dev/null | grep -oP 'CUDA Version:\s*\K[0-9]+\.[0-9]+' | head -1)
export CONDA_OVERRIDE_CUDA="${CUDA_VER:-12.4}"
# >= 2026.3 specifically: that's the GROMACS release that added native
# ff19SB support (this pipeline's protein force field, see
# docs/md_tier3.md) -- earlier releases fail pdb2gmx with "force field
# 'amber19sb' not found".
GMX_SPEC="gromacs=2026.3=nompi_cuda*"
if ! conda env list | grep -q '^md'; then
  mamba create -n md -c conda-forge -y \
    "$GMX_SPEC" \
    python=3.11 \
    mdanalysis \
    freesasa \
    biopython \
    pdbfixer \
    matplotlib \
    pandas \
    numpy \
    scipy \
    tqdm \
    pyyaml
fi
conda activate md
gmx --version 2>&1 | grep -qi "GPU support:.*CUDA" || {
  echo "ERROR: gromacs env is not CUDA-enabled (spec was $GMX_SPEC)."
  echo "Check https://anaconda.org/conda-forge/gromacs/files for the current nompi_cuda* build and bump GMX_SPEC."
  exit 1
}

# Repo
if [ ! -d /opt/Thesis-insilico ]; then
  git clone https://github.com/Wachirawut2023/Thesis-insilico.git /opt/Thesis-insilico
fi
cd /opt/Thesis-insilico
git checkout main
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
