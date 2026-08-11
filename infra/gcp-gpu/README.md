# GCP GPU setup for Tier-3 MD

Full runbook to run the 8-protein × 2-mode MD panel on a GCP `g2-standard-8`
with NVIDIA L4. **Estimated cost: ~$25–30**. **Wall time: ~24–32 hours
unattended**. Stays in the same GCP project / GCS bucket you already use.

## 0. Request GPU quota (one-time, may take minutes to hours)

Most GCP projects have **zero GPU quota** by default. You need at least
**1 NVIDIA L4 in us-central1**.

1. **Console → IAM & Admin → Quotas & System Limits**
2. Filter: **NVIDIA L4 GPUs** in **us-central1**
3. Check the box → **Edit Quotas** → request **1** → submit
4. Wait for approval email (typically minutes for small requests; can
   take hours during high demand).

While you wait, you can complete steps 1-3 below.

## 1. From Cloud Shell (or local gcloud) — one-time

```bash
gcloud config set project YOUR_PROJECT_ID
gcloud services enable compute.googleapis.com
git clone https://github.com/Wachirawut2023/Thesis-insilico.git ~/Thesis-insilico
cd ~/Thesis-insilico
git checkout main
```

## 2. Provision the GPU instance (~30 sec, then 10-15 min bootstrap)

```bash
cd infra/gcp-gpu
./provision.sh
```

Defaults: `g2-standard-8` (8 vCPU, 32 GB RAM, 1× NVIDIA L4) in
us-central1-a, 100 GB pd-balanced boot disk, Deep Learning VM image
(CUDA 12.3 + NVIDIA drivers preinstalled).

Wait for cloud-init to finish:

```bash
gcloud compute ssh as-m6a-gpu --zone=us-central1-a -- tail -f /var/log/bootstrap.log
```

Watch for `[gpu-bootstrap] done`.

## 3. Pull data/prepared/ from your GCS backup

The MD pipeline needs the receptor PDB files from the static run.

```bash
gcloud compute ssh as-m6a-gpu --zone=us-central1-a
sudo -i
cd /opt/Thesis-insilico

# Pull from GCS (faster, intra-cloud)
PROJECT=$(gcloud config get-value project)
gsutil cp gs://${PROJECT}-thesis-as-m6a-data/latest.tar.gz /tmp/
tar xzf /tmp/latest.tar.gz -C /tmp/

# Copy the prepared structures into place
cp -r /tmp/thesis-results*/data/prepared data/  2>/dev/null || \
  cp -r /tmp/results/data/prepared data/  2>/dev/null || \
  echo "Manual copy needed — check /tmp/ structure"

ls data/prepared/ | head
```

You should see 18+ `.clean.pdb` and `.pdbqt` files. If not, list /tmp/
and copy the prepared/ directory by hand.

**Alternative**: regenerate stages 1–2 of the static pipeline on this
GPU node (~10 min):

```bash
mamba install -n md -c conda-forge -c bioconda -y \
    pdb2pqr propka meeko vina fpocket
conda activate md
python scripts/01_fetch_structures.py
python scripts/02_prepare_receptors.py
```

## 4. Smoke test (~15 min) — STRONGLY recommended

Verify GROMACS + GPU + As-surrogate pipeline work end-to-end on TXN1
before committing 30 hours of compute.

```bash
cd /opt/Thesis-insilico
source /opt/miniforge/etc/profile.d/conda.sh && conda activate md
tmux new -s md

# Edit md.mdp to 5 ns instead of 50 (one-line sed)
cp scripts/md/mdp/md.mdp scripts/md/mdp/md_full.mdp
sed -i 's/nsteps          = 25000000/nsteps          = 2500000/' scripts/md/mdp/md.mdp

# Run TXN1 apo
GENES_FILTER=TXN1 MODES_FILTER=apo bash scripts/md/run_all.sh

# Check it produced output
ls -lh results/md/TXN1_apo/md.gro results/md/TXN1_apo/analysis/rmsd.png
```

If both files exist, the pipeline works. Restore the full-length config:

```bash
cp scripts/md/mdp/md_full.mdp scripts/md/mdp/md.mdp
```

## 5. Run the full panel (~24–30 hr unattended)

```bash
tmux a -t md      # if not still attached
bash scripts/md/run_all.sh
```

Detach tmux with `Ctrl-b d` and close your iPad. Reattach with
`tmux a -t md` later. Per-protein logs go to
`results/md/<gene>_<mode>/run.log`; orchestrator log to
`results/md/run_all.log`.

**Monitor while running**:

```bash
nvidia-smi -l 5                           # GPU utilisation (refresh 5s)
tail -f results/md/run_all.log
ls -ltr results/md/*/md.gro 2>/dev/null   # which proteins are done
```

A healthy run shows L4 GPU at 70–95% utilisation.

## 6. Pull results back

When complete:

```bash
cd /opt/Thesis-insilico
tar czf md_results.tar.gz results/md/

PROJECT=$(gcloud config get-value project)
gsutil cp md_results.tar.gz gs://${PROJECT}-thesis-as-m6a-data/
ls -lh md_results.tar.gz
```

Then from Cloud Shell to your iPad:

```bash
gsutil cp gs://${PROJECT}-thesis-as-m6a-data/md_results.tar.gz ~/
cloudshell download md_results.tar.gz
```

## 7. STOP THE INSTANCE — IMPORTANT

```bash
gcloud compute instances delete as-m6a-gpu --zone=us-central1-a --quiet
```

**An idle `g2-standard-8` still bills $0.82/hr (~$590/month).** Always
delete or stop the instance when done.

Use **stop** (not delete) if you want to keep the data on the disk for
later — you still pay ~$4/mo for the 100 GB volume but no vCPU/GPU charges.

## Cost summary

| Item | Rate | Duration | Total |
|---|---|---|---|
| g2-standard-8 + L4 | $0.82/hr | ~30 hr | **~$25** |
| 100 GB pd-balanced | $0.10/GB-mo | prorated 2 days | ~$0.65 |
| External IPv4 | $0.005/hr | ~30 hr | ~$0.15 |
| Egress (1-3 GB results) | $0.12/GB after 1 GB free | small | ~$0.25 |
| **Estimated total** | | | **~$26** |

If you still have GCP $300 free trial credit: this run is $0 out of pocket.

## Side-by-side with RunPod (the alternative I built earlier)

| Aspect | GCP L4 (this) | RunPod RTX 4090 (alternative) |
|---|---|---|
| Cost | ~$26 | ~$15-20 |
| Setup friction | Reuses existing GCP project/bucket/CLI | New account, separate CLI |
| GPU speed | L4 — modern Ada generation, well-suited to GROMACS | RTX 4090 — faster (consumer-grade) |
| Storage workflow | Direct gsutil to/from your GCS bucket | Need to bridge RunPod ↔ GCS |
| Existing data | One gsutil cp away | Need to push from GCS to RunPod |
| Tear-down | Familiar `gcloud instances delete` | RunPod web console |

If you have free trial credit left, **GCP is the better choice**. If not
and budget matters more than convenience, **RunPod is ~$10 cheaper**.

## Common issues

- **"Cannot meet GPU quota" on provision**: see step 0 above; quota
  request takes minutes to hours.
- **"Image not available" error**: the deep learning image family name
  occasionally changes. List current ones with:
  `gcloud compute images list --project=deeplearning-platform-release | grep -i cu`
  Update `IMAGE_FAMILY` in provision.sh if needed.
- **NVIDIA driver not installed**: `--metadata=install-nvidia-driver=True`
  in provision.sh triggers a first-boot driver install. If `nvidia-smi`
  fails, run `sudo /opt/deeplearning/install-driver.sh` manually and
  retry.
- **GROMACS slow / GPU idle**: `nvidia-smi -l 5` shows utilisation. If
  it's below 30%, GROMACS is running CPU-only. Check `gmx --version` shows
  `GPU support: CUDA` (not `OpenCL` — conda-forge's plain `gromacs` build
  offloads via OpenCL, which doesn't support Volta/Turing/Ampere-or-newer
  NVIDIA GPUs for compute and silently falls back). Fix: recreate the `md`
  env pinning the CUDA build explicitly,
  `mamba create -n md -c conda-forge -y "gromacs=2026.3=nompi_cuda*" ...`
  (see `cloud-init.yaml` step 4 for the full command) — a bare
  `cudatoolkit` dependency isn't enough, the build string itself must be
  pinned.
- **Out of memory on RBM15B**: 977 aa protein in a solvated box can
  need >32 GB. Either skip RBM15B from the panel or upgrade to
  `g2-standard-12` ($0.94/hr).
