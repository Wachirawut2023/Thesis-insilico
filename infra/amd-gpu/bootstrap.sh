#!/usr/bin/env bash
# One-shot bootstrap for a DigitalOcean AMD GPU Droplet (MI300X, ROCm).
#
# GROMACS on AMD needs a HIP build — conda-forge's `gromacs` package is
# CUDA-only, so unlike the GCP/RunPod NVIDIA paths we build GROMACS from
# source against the droplet's preinstalled ROCm stack. Everything else
# (analysis env, repo, run scripts) is identical to the NVIDIA runbooks.
#
# GMX_GPU=HIP requires GROMACS >= 2025.2 (added in that release; earlier
# versions only accept OFF/CUDA/OpenCL/SYCL for GMX_GPU and will fail
# cmake configure). As of 2026.1 the HIP backend still only offloads
# NBNxM (nonbonded) kernels — PME/bonded stay on CPU.
#
# Needs GROMACS >= 2026.3, not just >= 2025.2: that's when GROMACS added a
# *native* port of AMBER's ff19SB (this pipeline's protein force field,
# see docs/md_tier3.md), including per-amino-acid CMAP support — a feature
# ff19SB needs that earlier GROMACS releases don't support at all (no
# working workaround short of a patched build). Earlier we tried pinning
# to 2025.4 instead because GROMACS bumps its .tpr binary format (tpx
# version) across major releases — 2025.x writes tpx 137, 2026.x writes
# tpx 138, and MDAnalysis tops out at 137 as of its latest release (2.10.0)
# — but correct force-field physics wins over that; scripts/md/analyze.py
# was changed to build its MDAnalysis Universe from the .gro output instead
# of the .tpr, which sidesteps the tpx-version dependency entirely.
#
# Designed to be pasted into the droplet's web console / SSH session, or
# run non-interactively via cloud-init (see cloud-init.yaml).

set -euo pipefail
exec > >(tee -a /var/log/bootstrap.log) 2>&1
echo "[amd-bootstrap] starting $(date -Is)"

GROMACS_VERSION="${GROMACS_VERSION:-2026.3}"
WORK_DIR="${WORK_DIR:-/opt/Thesis-insilico}"
REPO_URL="${REPO_URL:-https://github.com/Wachirawut2023/Thesis-insilico.git}"
BRANCH="${BRANCH:-claude/tier3-md-amd-gpu-oydfoj}"

# 1. Verify ROCm + GPU visible.
# Note: rocminfo also lists generic ISA family names like "gfx9-4-generic"
# (the dash breaks a loose gfx[0-9a-z]* match into a bogus "gfx9"), so require
# 2-4 digits to land on the real arch code (e.g. gfx942), not the family stub.
rocminfo | grep -i "gfx" || { echo "ERROR: no AMD GPU visible via rocminfo"; exit 1; }
AMDGPU_TARGET="$(rocminfo | grep -oE 'gfx[0-9]{2,4}[a-z]?' | sort -u | head -1)"
echo "[amd-bootstrap] detected GPU target: $AMDGPU_TARGET"

apt-get update
apt-get install -y --no-install-recommends \
  build-essential cmake git curl tmux htop rsync unzip jq \
  libfftw3-dev ninja-build

# 2. Miniforge for the Python analysis env (not used to build GROMACS).
if [ ! -d /opt/miniforge ]; then
  curl -L -o /tmp/miniforge.sh \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
  bash /tmp/miniforge.sh -b -p /opt/miniforge
  rm /tmp/miniforge.sh
fi
export PATH="/opt/miniforge/bin:$PATH"
source /opt/miniforge/etc/profile.d/conda.sh

if ! conda env list | grep -q '^md'; then
  mamba create -n md -c conda-forge -y \
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

# 3. Build GROMACS with HIP (AMD GPU) support from source.
# Uses the system ROCm install. Preinstalled on DigitalOcean's AI/ML-ready
# AMD GPU Droplet image at /opt/rocm; on a bare Ubuntu image with ROCm from
# apt (this box), it lands under /usr instead (see `hipconfig`'s ROCM_PATH).
ROCM_PREFIX="$(hipconfig -R 2>/dev/null || echo /opt/rocm)"
if ! command -v gmx >/dev/null 2>&1 || ! gmx --version 2>&1 | grep -qi "HIP"; then
  echo "[amd-bootstrap] building GROMACS $GROMACS_VERSION with HIP for $AMDGPU_TARGET"
  cd /opt
  if [ ! -d "gromacs-${GROMACS_VERSION}" ]; then
    curl -L -o gromacs.tar.gz \
      "https://ftp.gromacs.org/gromacs/gromacs-${GROMACS_VERSION}.tar.gz"
    tar xzf gromacs.tar.gz
    rm gromacs.tar.gz
  fi
  cd "gromacs-${GROMACS_VERSION}"
  mkdir -p build && cd build
  # Ubuntu's packaged hip-config.cmake has a templating bug (`if("")` around
  # the block that would set these) so find_package(hip) never populates
  # them, and GROMACS's arch-support check silently invokes an empty
  # COMMAND and reports every architecture as rejected. Set explicitly.
  cmake .. \
    -DGMX_GPU=HIP \
    -DGMX_HIP_TARGET_ARCH="$AMDGPU_TARGET" \
    -DCMAKE_PREFIX_PATH="$ROCM_PREFIX" \
    -DHIP_HIPCC_EXECUTABLE="$(command -v hipcc)" \
    -DHIP_HIPCONFIG_EXECUTABLE="$(command -v hipconfig)" \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DGMX_MPI=OFF \
    -DCMAKE_INSTALL_PREFIX=/opt/gromacs \
    -DCMAKE_BUILD_TYPE=Release
  make -j"$(nproc)"
  make install
  cd /
fi
# GMXRC references a couple of variables ($shell, $GMXLDLIB) without
# defaults, which is fatal under this script's `set -u` — relax nounset
# just for the source.
set +u
# shellcheck disable=SC1091
source /opt/gromacs/bin/GMXRC
set -u

echo "[amd-bootstrap] gmx version check:"
gmx --version 2>&1 | grep -E "GROMACS|GPU support|HIP" || true

# 4. Repo.
if [ ! -d "$WORK_DIR" ]; then
  git clone "$REPO_URL" "$WORK_DIR"
fi
cd "$WORK_DIR"
git fetch origin
git checkout "$BRANCH"
git pull --ff-only
chmod +x scripts/md/run_md.sh scripts/md/run_all.sh

# 5. Persist environment for future shells (SSH, tmux, etc.)
cat > /etc/profile.d/md-env.sh <<'EOF'
export PATH="/opt/miniforge/bin:$PATH"
source /opt/miniforge/etc/profile.d/conda.sh
conda activate md 2>/dev/null || true
source /opt/gromacs/bin/GMXRC
EOF

echo
echo "=== bootstrap done $(date -Is) ==="
rocm-smi --showproductname || true
echo
echo "Next:"
echo "  1) Pull data/prepared/ (receptor structures from the static run)"
echo "  2) tmux new -s md"
echo "  3) bash scripts/md/run_all.sh"
