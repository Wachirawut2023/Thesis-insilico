# AMD GPU Droplet setup for Tier-3 MD

Runbook to run the 8-protein × 2-mode MD panel on a **DigitalOcean AMD GPU
Droplet (MI300X)**, as an alternative to Google Colab (which hit resource
/session-length limits for the 50 ns × 16-run panel) and to the NVIDIA
paths in `infra/gcp-gpu/` and `infra/runpod/`.

**Wall time: ~24–32 hours unattended.** Cost depends on DO's current
GPU Droplet pricing for MI300X (check your DO console — this changes
more often than NVIDIA on-demand pricing).

## Why this is different from the NVIDIA runbooks

GROMACS's GPU offload (`-nb gpu -pme gpu` in `scripts/md/run_md.sh`) is
backend-agnostic — the same run scripts work unchanged on AMD, and
`MD_DEVICE=auto` (the default) detects an AMD GPU via `rocminfo` the same
way it detects an NVIDIA one via `nvidia-smi`. The only real difference is
**how GROMACS gets built**:

- NVIDIA path: `mamba install gromacs` from conda-forge (prebuilt CUDA binary).
- AMD path: conda-forge does not ship a HIP/ROCm GROMACS build, so
  `bootstrap.sh` / `cloud-init.yaml` here **compile GROMACS from source**
  with `-DGMX_GPU=HIP` against the droplet's ROCm install. This adds
  ~20-30 min to first boot but is otherwise a one-time cost.

## 1. Prerequisites

1. A DigitalOcean account with GPU Droplet access enabled (MI300X
   capacity is region-limited; check availability in your console under
   **Droplets → GPU Droplets**, or `doctl compute size list | grep -i gpu`).
2. `doctl` installed and authenticated:
   ```bash
   doctl auth init
   ```
3. An SSH key added to your DO account:
   ```bash
   doctl compute ssh-key list
   # if empty, add one:
   doctl compute ssh-key import my-key --public-key-file ~/.ssh/id_ed25519.pub
   ```

## 2. Provision the droplet (~1 min to create, ~30-45 min bootstrap)

```bash
cd infra/amd-gpu
SSH_KEY_NAME=my-key ./provision.sh
```

Defaults: 1× MI300X (`gpu-mi300x1-192gb`), DO's AI/ML-ready AMD GPU image
(Ubuntu 22.04 + ROCm preinstalled), region `nyc2`. **Verify the size and
image slugs before running** — DO's GPU Droplet catalogue is newer and
slug names may differ by region/account tier:

```bash
doctl compute size list | grep -i gpu
doctl compute image list-distribution --public | grep -i amd
```

Override with `SIZE=... IMAGE=... REGION=... ./provision.sh` if the
defaults don't match what's available to you.

Watch bootstrap progress:

```bash
ssh root@<droplet-ip> 'tail -f /var/log/bootstrap.log'
```

Watch for `[amd-bootstrap] done`. This step builds GROMACS from source,
so it's slower than the NVIDIA conda-install path — budget 30-45 min.

## 3. Get the receptor structures onto the droplet

The MD pipeline needs `data/prepared/*.clean.pdb` from the static
pipeline (stages 1-2). Two options:

**A. Copy from wherever you ran the static pipeline** (e.g. your GCS
bucket used in `infra/gcp/`):

```bash
ssh root@<droplet-ip>
cd /opt/Thesis-insilico
# if you have gsutil creds, or just scp from your machine:
#   scp -r data/prepared root@<droplet-ip>:/opt/Thesis-insilico/data/
ls data/prepared/ | head
```

**B. Regenerate stages 1-2 on the droplet** (~10 min, needs the
non-GPU tools too):

```bash
source /etc/profile.d/md-env.sh
mamba install -n md -c conda-forge -c bioconda -y \
    pdb2pqr propka meeko vina fpocket
python scripts/01_fetch_structures.py
python scripts/02_prepare_receptors.py
```

## 4. Smoke test (~15 min) — STRONGLY recommended

Same as the NVIDIA runbooks: verify GROMACS-HIP + the As-surrogate
pipeline work end-to-end on one protein before committing ~30 hours.

```bash
cd /opt/Thesis-insilico
source /etc/profile.d/md-env.sh
tmux new -s md

cp scripts/md/mdp/md.mdp scripts/md/mdp/md_full.mdp
sed -i 's/nsteps          = 25000000/nsteps          = 2500000/' scripts/md/mdp/md.mdp

GENES_FILTER=TXN1 MODES_FILTER=apo bash scripts/md/run_all.sh

ls -lh results/md/TXN1_apo/md.gro results/md/TXN1_apo/analysis/rmsd.png
```

Check `results/md/TXN1_apo/run.log` for `GPU support: enabled` and a HIP
device listed. If both output files exist, restore the full-length
config:

```bash
cp scripts/md/mdp/md_full.mdp scripts/md/mdp/md.mdp
```

## 5. Run the full panel (~24-30 hr unattended)

```bash
tmux a -t md      # if not still attached
bash scripts/md/run_all.sh
```

Detach with `Ctrl-b d`; reattach later with `tmux a -t md`. Per-protein
logs: `results/md/<gene>_<mode>/run.log`; orchestrator log:
`results/md/run_all.log`.

**Monitor while running**:

```bash
rocm-smi -l 5                            # GPU utilisation (ROCm equivalent of nvidia-smi -l)
tail -f results/md/run_all.log
ls -ltr results/md/*/md.gro 2>/dev/null  # which proteins are done
```

A healthy run shows the MI300X at high (>70%) utilisation during the
production step.

## 6. Pull results back

```bash
cd infra/amd-gpu
./fetch-results.sh
```

This tarballs `results/md/` on the droplet and `scp`s it to a local
`md-results-<timestamp>/` directory.

## 7. TEAR DOWN THE DROPLET — IMPORTANT

```bash
doctl compute droplet delete as-m6a-amd-gpu
```

**An idle GPU droplet keeps billing per hour.** Always delete when done;
there is no "stop but keep GPU allocation" option worth using here since
data has already been pulled off in step 6 — if you want to keep the
disk without the GPU cost, snapshot the droplet first
(`doctl compute droplet-action snapshot`) instead of leaving it running.

## Common issues

- **"size not found" / "image not found" on provision**: MI300X
  availability and slug names vary by region and account. Run
  `doctl compute size list | grep -i gpu` and
  `doctl compute image list-distribution --public | grep -i amd` and
  override `SIZE`/`IMAGE`/`REGION`.
- **`rocminfo` shows no GPU / GROMACS build fails to detect `gfx*`
  target**: confirm you're on an actual GPU Droplet (not a CPU droplet)
  and that ROCm is installed (`dpkg -l | grep rocm`). If ROCm is
  missing, the "AI/ML ready" image wasn't used — reprovision with the
  correct `IMAGE` slug.
- **GROMACS build takes forever / OOMs**: the source build uses
  `-j$(nproc)`; on droplets with less RAM per core than expected, cap
  parallelism with `make -jN` for a smaller N (edit `bootstrap.sh`).
- **`gmx --version` doesn't show HIP**: the `cmake` step likely fell
  back to CPU-only because `-DGMX_HIP_TARGET_ARCH` didn't match your
  card. Confirm with `rocminfo | grep gfx` and re-run the build with
  the correct target set explicitly.
- **Out of memory on RBM15B** (977 aa, large solvated box): same
  caveat as the NVIDIA runbooks — MI300X's large VRAM/HBM should absorb
  this easily, but if host RAM (not GPU memory) is the constraint,
  reduce concurrent runs.

## Why Colab didn't work for Tier 3

Colab's free/Pro tiers cap session length (12-24 hr) and can silently
disconnect or reclaim the GPU mid-run, which is fatal for a 50 ns
production run taking several hours per protein × mode × 16 runs. A
dedicated cloud GPU instance (this droplet, or the GCP/RunPod
alternatives) gives an uninterrupted, resumable, tmux-persisted session
for the full ~24-30 hr unattended panel.
