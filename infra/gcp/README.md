# GCP Compute Engine setup — As(III) m6A pipeline (on-demand)

End-to-end runbook for GCP. Total cost for a full-panel run on `n2d-standard-8`:
**~$2.50** at list price, **$0** if you have the $300 free-trial credit.

## TL;DR

```bash
# Local machine, one-time:
gcloud auth login
gcloud config set project YOUR_PROJECT_ID
gcloud services enable compute.googleapis.com

# Provision:
cd infra/gcp
./provision.sh                               # ~30 s; prints the external IP

# Wait for cloud-init (~5–10 min):
gcloud compute ssh as-m6a --zone=us-central1-a -- tail -f /var/log/bootstrap.log
# wait for "[bootstrap] done"

# Run:
gcloud compute ssh as-m6a --zone=us-central1-a
sudo -i
cd /opt/Thesis-insilico
tmux new -s pipeline
bash infra/gcp/run-pipeline.sh               # ~3–6 hr

# Local: pull results, tear down.
./fetch-results.sh
gcloud compute instances delete as-m6a --zone=us-central1-a --quiet
```

## 1. GCP account + billing (one-time, ~5 min)

1. https://console.cloud.google.com — sign in with a Google account.
2. **New users**: accept the **$300 / 90-day free trial** when prompted.
   A credit card is required for identity check but you won't be charged
   while in the free tier.
3. Create or pick a Project (sidebar **Project selector → New Project**).
   Note the **Project ID** (not the display name).
4. Make sure billing is linked: **Billing → Account management → Link account**.

## 2. Pick a connection method

You have three options. Pick whichever feels easiest.

### 2A. Cloud Shell (no local install — recommended)

A free Linux terminal inside the browser, with `gcloud` pre-installed and a
5 GB persistent home directory.

1. https://console.cloud.google.com → click the terminal icon (`>_`) top-right.
2. From Cloud Shell:
   ```bash
   gcloud config set project YOUR_PROJECT_ID
   gcloud services enable compute.googleapis.com
   git clone https://github.com/Wachirawut2023/Thesis-insilico.git
   cd Thesis-insilico && git checkout claude/arsenic-m6a-inhibition-model-6pwWS
   cd infra/gcp
   ```

Everything below works identically from Cloud Shell or from a laptop with
`gcloud` installed.

### 2B. Local `gcloud` CLI

```bash
# macOS
brew install --cask google-cloud-sdk
# Linux: https://cloud.google.com/sdk/docs/install-sdk

gcloud auth login
gcloud config set project YOUR_PROJECT_ID
gcloud services enable compute.googleapis.com
```

### 2C. Pure Console UI (no terminal at all)

You can create the VM by clicking through the Console — see
[Manual VM creation](#manual-vm-creation-console-only) at the bottom of this
file. Then use **SSH-in-browser** (Option 5A below) to connect.

## 3. Provision the instance

```bash
cd infra/gcp
./provision.sh
```

Defaults:

| Setting | Value | Why |
|---|---|---|
| machine type | `n2d-standard-8` | 8 dedicated AMD EPYC vCPU, 32 GB RAM, $0.388/hr us-central1 — same chip family as the Hetzner CCX33 we validated against |
| zone | `us-central1-a` | Iowa; among the cheapest regions globally |
| image | `ubuntu-2404-lts-amd64` | Cloud-init compatible |
| disk | 50 GB `pd-balanced` | Cheaper than pd-ssd; enough for pipeline |
| network | default VPC | SSH (port 22) open via IAP/gcloud by default |

Override via env vars:

```bash
MACHINE_TYPE=c3-standard-8 ZONE=us-east1-b ./provision.sh
# c3-standard-8 ($0.418/hr) is newer Intel Sapphire Rapids, ~10% faster per core
```

## 4. Wait for cloud-init

The startup script installs Miniforge, clones the repo, and builds the conda
env. Track progress:

```bash
gcloud compute ssh as-m6a --zone=us-central1-a -- tail -f /var/log/bootstrap.log
```

Done when you see `[bootstrap] done <timestamp>` (~5–10 min).

If `tail` exits with "file not found" the script hasn't started yet — wait
30 s and retry.

## 5. Connect and run the pipeline

Two equivalent ways in:

### 5A. SSH-in-browser (no CLI needed)

1. Console → **Compute Engine → VM instances**
2. Click the **SSH** button on the `as-m6a` row → opens a browser terminal
3. You're logged in as your Google username

### 5B. `gcloud compute ssh` (from Cloud Shell or local laptop)

```bash
gcloud compute ssh as-m6a --zone=us-central1-a
```

### Then, regardless of how you connected:

```bash
sudo -i                                        # bootstrap installs to /opt/
cd /opt/Thesis-insilico
tmux new -s pipeline                           # Ctrl-b d to detach; `tmux a` to reattach

# Recommended: smoke test first (3 targets, ~30 min)
source /opt/miniforge/etc/profile.d/conda.sh && conda activate arsenic-m6a
make smoke

# Then full panel
bash infra/gcp/run-pipeline.sh
```

Detach tmux and close the browser tab / laptop — the pipeline keeps running.
Re-attach later with the SSH button or `gcloud compute ssh`, then `tmux a -t pipeline`.

## 6. Pull results back

Three ways, pick whichever fits your connection method:

### 6A. SSH-in-browser — click to download

In the browser SSH window, click the **gear icon (⚙)** at top-right →
**Download file** → paste the path, e.g.:

```
/opt/Thesis-insilico/results/REPORT.md
/opt/Thesis-insilico/results/composite_ranking.tsv
/opt/Thesis-insilico/results/figures/family_heatmap.png
```

Each download goes straight to your browser's downloads folder. Best for
grabbing a handful of small files.

### 6B. Cloud Shell — `fetch-results.sh` + `cloudshell download`

```bash
# In Cloud Shell:
cd ~/Thesis-insilico/infra/gcp
./fetch-results.sh
tar czf results.tar.gz results-pulled/
cloudshell download results.tar.gz       # opens a download dialog in your browser
```

### 6C. Local gcloud — same `fetch-results.sh`

```bash
cd infra/gcp
./fetch-results.sh
# -> ./results-pulled/<timestamp>/results/ on your laptop
```

## 7. Destroy — STOP BILLING

```bash
gcloud compute instances delete as-m6a --zone=us-central1-a --quiet

# Verify nothing else is running:
gcloud compute instances list
```

**An idle n2d-standard-8 still bills $0.388/hr (~$280/month).** Always delete
or stop the instance when done. Stopping (`instances stop`) halts vCPU billing
but you still pay for the boot disk (~$2/mo).

## Re-installing the VM while preserving your data

If you need a fresh VM (different machine type, fresh OS, recovering from a
broken state) but want to keep your `results/`, `data/prepared/`, and ligand
artifacts, use a GCS bucket as the persistent layer between VM lifecycles.

The pipeline regenerates everything in `data/` from public PDB/AlphaFold in
~10 min, so the only thing you can't recreate cheaply is `results/` (your
actual scientific output). The backup script saves all three.

**Cost**: ~$0.02/GB-month for GCS standard storage. A typical run produces
<200 MB so the persistent storage is rounding error.

### Step 1 — back up the current VM (before deleting it)

On the existing VM:

```bash
cd /opt/Thesis-insilico
git pull                                          # make sure backup scripts are present
bash infra/gcp/backup-to-gcs.sh
```

The script:
- creates a GCS bucket `gs://<PROJECT>-thesis-as-m6a-data` if it doesn't exist
- tars `results/`, `data/prepared/`, and `ligands/as3.pdbqt`
- uploads as both `thesis-backup-<timestamp>.tar.gz` and `latest.tar.gz`

### Step 2 — delete the old VM

```bash
gcloud compute instances delete as-m6a --zone=us-central1-a --quiet
```

Billing for the VM stops immediately. The GCS bucket persists.

### Step 3 — create the new VM

Same `provision.sh` as before, or use the Console UI:

```bash
cd infra/gcp
./provision.sh                                    # ~30 s
# wait ~5–10 min for cloud-init
gcloud compute ssh as-m6a --zone=us-central1-a -- tail -f /var/log/bootstrap.log
# wait for [bootstrap] done
```

### Step 4 — restore your data on the new VM

```bash
gcloud compute ssh as-m6a --zone=us-central1-a
sudo -i
cd /opt/Thesis-insilico
bash infra/gcp/restore-from-gcs.sh
```

You'll find `results/` and `data/prepared/` populated as they were before the
old VM was deleted. The conda env was rebuilt fresh by cloud-init, so you have
a clean OS + your existing outputs.

### Optional — pick a specific backup instead of `latest`

```bash
# list available backups
gsutil ls gs://$(gcloud config get-value project)-thesis-as-m6a-data/

# restore a specific one
ARCHIVE=thesis-backup-20260518_153022.tar.gz bash infra/gcp/restore-from-gcs.sh
```

### Tear down the GCS bucket if you don't need it any more

```bash
gsutil -m rm -r gs://$(gcloud config get-value project)-thesis-as-m6a-data
```

## Machine type cheat-sheet (8 vCPU / 32 GB tier)

| Type | CPU | $/hr (us-central1) | 6 hr cost | Notes |
|---|---|---|---|---|
| **n2d-standard-8** | AMD EPYC Milan | **$0.388** | **$2.33** | Default — proven for our workload |
| n2-standard-8 | Intel Cascade Lake | $0.388 | $2.33 | Older Intel; similar perf to n2d |
| c3d-standard-8 | AMD EPYC Genoa | $0.398 | $2.39 | Newer than n2d, ~10% faster IPC |
| c3-standard-8 | Intel Sapphire Rapids | $0.418 | $2.51 | Fastest per core; ~10% better than n2d |
| e2-standard-8 | shared / burstable | $0.268 | $1.61 | **Avoid** — burstable CPU makes Vina wall time unpredictable |

For Vina/AutoDock workloads, dedicated CPU (n2d / c3d / c3) beats shared-core
e2 by a wide margin in *real* wall time, often making e2 more expensive once
you account for actual run length.

## Cost estimate breakdown

| Item | Rate | 6 hr |
|---|---|---|
| n2d-standard-8 vCPU+RAM | $0.388/hr | $2.33 |
| 50 GB pd-balanced | $0.10/GB-mo prorated | ~$0.04 |
| External IPv4 (auto-assigned) | $0.005/hr | $0.03 |
| Egress (~1–2 GB of results) | 1 GB free, then $0.12/GB | ≤$0.12 |
| **Total** | | **~$2.50** |

With the **$300 free trial**, that's <1% of your credit. You can run the full
pipeline ~120 times before the credit is exhausted.

## Region picks

If you're not in North America and latency matters for `gcloud ssh`:

| Region | Zone example | Notes |
|---|---|---|
| us-central1 | us-central1-a | Iowa; cheapest, default |
| us-east1 | us-east1-b | South Carolina; same price |
| europe-west4 | europe-west4-a | Netherlands; +5–10% |
| asia-southeast1 | asia-southeast1-a | Singapore; +20%, but lower latency from Thailand |

Run `gcloud compute regions list` to see all options. Price differences are
small relative to your free trial budget — pick by latency, not cost.

## Common gotchas

- **`gcloud compute ssh` permission denied.** Wait for cloud-init to finish — the
  metadata SSH key sync can take 30–60 s after instance creation.
- **`API not enabled`** errors on first run — `gcloud services enable
  compute.googleapis.com` (provision.sh does this automatically).
- **Free trial expiry.** The $300 expires after 90 days, even if unspent. After
  that, you'll need to upgrade to a paid account or your instance will be
  paused.
- **Quota errors on n2d-standard-8.** Default CPU quota in a region is 24. The
  8-core instance fits. If you ever scale up: **IAM & Admin → Quotas → search
  for "CPUs (all regions)"**.
- **Conda env build is slow on first run** (~5 min, mamba downloading
  dependencies). Subsequent updates are fast.

## Why on-demand vs spot (preemptible) for this workload

- **On-demand** (what this setup uses): instance runs until you delete it.
  Guaranteed availability. Cost: $0.388/hr.
- **Spot** (preemptible): ~75% cheaper (~$0.10/hr) but GCP can kill it with
  30 s notice. Max lifetime 24 hr.

You explicitly chose **on-demand for no-babysitting**. For a 6-hour run that
costs $2.30, the spot savings (~$1.70) aren't worth managing interruptions.

If you ever want to switch to spot: add `--provisioning-model=SPOT
--instance-termination-action=DELETE` to the `gcloud compute instances create`
in `provision.sh`. The pipeline writes per-protein outputs, so it's
checkpoint-friendly anyway.

## Manual VM creation (Console only)

If you'd rather click through the GCP Console than run `provision.sh`:

1. **Compute Engine → VM instances → CREATE INSTANCE**
2. **Name**: `as-m6a`
3. **Region**: `us-central1`, **Zone**: `us-central1-a`
4. **Machine configuration**:
   - Series: **N2D**
   - Machine type: **n2d-standard-8** (8 vCPU, 32 GB memory)
5. **Boot disk → Change**:
   - Operating system: **Ubuntu**
   - Version: **Ubuntu 24.04 LTS Minimal (amd64)**
   - Boot disk type: **Balanced persistent disk**
   - Size: **50 GB**
6. **Advanced options → Networking**: leave defaults (default VPC, ephemeral
   external IP). SSH is allowed by default via IAP/Console.
7. **Advanced options → Management → Metadata → Add item**:
   - Key: `user-data`
   - Value: paste the **entire contents** of `infra/gcp/cloud-init.yaml`
     (open the file, copy everything starting from `#cloud-config`)
8. Click **CREATE**.

Wait 5–10 min for cloud-init to finish. Track progress: click **SSH** on the
instance row in the Console, then `tail -f /var/log/bootstrap.log`.

Once `[bootstrap] done` shows up, jump to [step 5](#5-connect-and-run-the-pipeline).

## Tear-down checklist

```bash
gcloud compute instances list                  # what's running
gcloud compute instances delete as-m6a --zone=us-central1-a --quiet
gcloud compute disks list                      # any orphan disks?
gcloud compute addresses list                  # any reserved static IPs?
```

The instance creation here doesn't reserve a static IP or create extra disks,
so a single `instances delete` is normally enough to stop all billing.

You can also tear down from the Console: **Compute Engine → VM instances →**
check the box next to `as-m6a` **→ DELETE** at the top.
