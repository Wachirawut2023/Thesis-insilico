# RunPod GPU setup for Tier-3 MD

End-to-end runbook for running the 8-protein × 2-mode MD panel on a single
RunPod RTX 4090. Total cost: **~$20–30**. Wall time: **~24–30 hours
running serially**; faster if you spin up two pods in parallel.

## 1. Create a RunPod account (~5 min)

1. Sign up at https://www.runpod.io
2. Add **$25 credit** via card or crypto (you'll spend ~$20)
3. (Optional) install the CLI: https://docs.runpod.io/cli/install — not
   required if you only use the web Console.

## 2. Create a pod with GPU + persistent storage

In the RunPod Console:

1. **Deploy → GPU Cloud → Secure Cloud**.
2. Filter by **RTX 4090** (≈ $0.40/hr) or **A40** (≈ $0.50/hr, more VRAM).
3. Template: **"PyTorch 2.4"** or **"Ubuntu 22.04 + CUDA 12.x"** — anything
   with CUDA pre-installed.
4. **Container Disk**: 50 GB minimum.
5. **Volume Disk**: 50 GB (persistent — survives if you stop the pod).
6. **Expose SSH port** (port 22) so you can connect from your laptop /
   Cloud Shell.
7. Click **Deploy**. Pod starts in ~1 minute.

## 3. Connect to the pod

Once running, the RunPod Console shows an **SSH command** like:

```
ssh root@<pod-ip> -p <port> -i ~/.ssh/id_ed25519
```

Use this from your laptop terminal or Cloud Shell. RunPod also has a
**Web Terminal** button — works from iPad Safari.

## 4. One-time bootstrap (~10 min)

Inside the pod:

```bash
# 1. Install GROMACS with GPU support + dependencies
apt-get update
apt-get install -y git tmux htop rsync build-essential

# 2. Install Miniforge
curl -L -o /tmp/miniforge.sh \
  "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash /tmp/miniforge.sh -b -p /opt/miniforge
source /opt/miniforge/etc/profile.d/conda.sh

# 3. Create the MD environment (GROMACS with GPU, MDAnalysis, freesasa).
# IMPORTANT: conda-forge's *plain* "nompi" gromacs build offloads to GPU
# via OpenCL, and GROMACS' OpenCL backend does not support Volta/Turing/
# Ampere-or-newer NVIDIA GPUs (the RTX 4090 included) for compute -- gmx
# enumerates the card, reports it "incompatible", and silently drops to
# CPU-only with no error at all. conda-forge also ships a real
# CUDA-enabled build (build string "nompi_cuda*"); pin it explicitly, since
# the solver won't prefer it on its own.
CUDA_VER=$(nvidia-smi 2>/dev/null | grep -oP 'CUDA Version:\s*\K[0-9]+\.[0-9]+' | head -1)
export CONDA_OVERRIDE_CUDA="${CUDA_VER:-12.4}"
GMX_SPEC="gromacs=2024.5=nompi_cuda*"
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
conda activate md

# Verify GROMACS actually got the CUDA build -- fail here rather than
# discover 24-30 GPU-hours later that it ran CPU-only the whole time.
gmx --version | grep -i gpu
gmx --version 2>&1 | grep -qi "GPU support:.*CUDA" || {
  echo "ERROR: not a CUDA build (spec was $GMX_SPEC)."
  echo "Check https://anaconda.org/conda-forge/gromacs/files for the current nompi_cuda* build and bump GMX_SPEC."
  exit 1
}
nvidia-smi
```

`gmx --version` should print `GPU support: CUDA`. `nvidia-smi` should show
your RTX 4090.

## 5. Pull repo + prepared structures

```bash
# Clone the repo
git clone https://github.com/Wachirawut2023/Thesis-insilico.git /opt/Thesis-insilico
cd /opt/Thesis-insilico
git checkout claude/arsenic-m6a-inhibition-model-6pwWS

# Pull data/prepared/ from your GCS backup
# (alternative: re-run static pipeline stages 1-2 here — ~10 min)
gsutil cp gs://YOUR_PROJECT-thesis-as-m6a-data/latest.tar.gz /tmp/
tar xzf /tmp/latest.tar.gz -C /tmp/
cp -r /tmp/thesis-results*/data/prepared data/

# Or, easiest if gsutil isn't set up on the pod:
# from your laptop or Cloud Shell, scp the prepared structures over:
#   gcloud compute scp --recurse data/prepared/ root@<runpod-ip>:/opt/Thesis-insilico/data/
```

If you'd rather regenerate the prepared structures from scratch:

```bash
mamba install -c conda-forge -c bioconda -y \
    pdb2pqr propka meeko vina fpocket
conda activate md   # same env
python scripts/01_fetch_structures.py
python scripts/02_prepare_receptors.py
```

## 6. Run the MD panel (~24–30 hr)

```bash
chmod +x scripts/md/run_md.sh scripts/md/run_all.sh
tmux new -s md

# Full panel: all 8 proteins, both modes (apo + bound)
bash scripts/md/run_all.sh

# Or filter:
GENES_FILTER="TXN1,PIN1,METTL3,FTO" bash scripts/md/run_all.sh
MODES_FILTER="apo" bash scripts/md/run_all.sh
```

`Ctrl-b d` to detach tmux and close your iPad. Reattach later with
`tmux a -t md`. Each protein logs to `results/md/<gene>_<mode>/run.log`;
the orchestrator logs to `results/md/run_all.log`.

**Monitor**:
```bash
nvidia-smi -l 5            # GPU utilisation (refresh every 5 s)
tail -f results/md/run_all.log
```

A healthy production run shows GPU at 80–100% utilisation. If
utilisation is below 30%, GROMACS is CPU-bottlenecked and you need
more CPU threads (try a pod with more vCPUs).

## 7. Pull results back

When `run_all.sh` finishes (last line of `results/md/run_all.log` says
"complete"):

```bash
cd /opt/Thesis-insilico
tar czf md_results.tar.gz results/md/

# Push to your GCS bucket (recommended)
gsutil cp md_results.tar.gz gs://YOUR_PROJECT-thesis-as-m6a-data/md_results.tar.gz

# Or use runpodctl / scp to your laptop directly
```

## 8. STOP THE POD — IMPORTANT

```bash
# In the RunPod Console: Stop or Terminate the pod
# "Stop" preserves the volume (you still pay ~$0.10/day for 50 GB storage)
# "Terminate" destroys everything (zero ongoing cost)
```

**An idle RTX 4090 still bills at $0.40/hr (~$10/day).** Always stop or
terminate when not actively running MD.

## Cost & time summary

| Item | Per protein | Total (8 × 2 = 16 sims) |
|---|---|---|
| GPU wall time | ~1.5–2 hr | ~24–32 GPU-hr |
| GPU cost (RTX 4090 @ $0.40/hr) | ~$0.60–0.80 | **~$10–13** |
| Pod overhead (CPU/storage idle) | — | ~$2–5 |
| Storage (50 GB volume, ~1 day) | — | ~$0.10 |
| **Estimated total** | | **~$15–20** |

Add ~$5 buffer for re-runs of any failed simulations and you stay
under $25 for the whole Tier-3 phase.

## Common issues

- **"Cannot allocate memory" during solvation**: pod has too little RAM
  for a large protein (RBM15B is 977 aa). Scale up to a pod with 32+ GB
  RAM, or skip RBM15B in `GENES_FILTER`.
- **"No suitable GPU found"**: `nvidia-smi` returns nothing → pod was
  provisioned without GPU passthrough. Re-deploy with a different
  template that explicitly has CUDA.
- **GROMACS warns about charge mismatch after As surrogate placement**:
  expected — the Mg²⁺ adds +2 charge that `gmx genion` already
  neutralised; either rerun ion placement after the Mg insertion or
  accept the residual charge (it doesn't affect dynamics meaningfully).
- **`gmx pdb2gmx` complains about unrecognised residues**: usually means
  a HETATM cofactor wasn't whitelisted in stage 2. Strip the offending
  residue from the input PDB and re-run.

## Sanity check: a 5-ns smoke run before the full panel

Before committing 30 hours of GPU, run just TXN1 apo with a shorter
production run to confirm the whole pipeline works:

```bash
# Quick edit: 5 ns instead of 50
sed 's/nsteps          = 25000000/nsteps          = 2500000/' \
    scripts/md/mdp/md.mdp > scripts/md/mdp/md_smoke.mdp
cp scripts/md/mdp/md.mdp scripts/md/mdp/md_full.mdp
cp scripts/md/mdp/md_smoke.mdp scripts/md/mdp/md.mdp

GENES_FILTER=TXN1 MODES_FILTER=apo bash scripts/md/run_all.sh
# ~15 min on RTX 4090

# If TXN1 produced md.gro + analysis/ outputs, restore and run the full panel:
cp scripts/md/mdp/md_full.mdp scripts/md/mdp/md.mdp
bash scripts/md/run_all.sh
```
