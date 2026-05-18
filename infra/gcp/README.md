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

## 2. Install + auth the gcloud CLI (one-time, ~2 min)

```bash
# macOS
brew install --cask google-cloud-sdk

# Linux  (Debian/Ubuntu): see https://cloud.google.com/sdk/docs/install-sdk
```

```bash
gcloud auth login                              # browser-based OAuth
gcloud config set project YOUR_PROJECT_ID
gcloud services enable compute.googleapis.com  # one-time, ~30 s
```

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

## 5. Run the pipeline

```bash
gcloud compute ssh as-m6a --zone=us-central1-a
sudo -i                                        # bootstrap installs to /opt/, so use root
cd /opt/Thesis-insilico
tmux new -s pipeline                           # detach Ctrl-b d; reattach `tmux a`

# Recommended: smoke test first (3 targets, ~30 min)
source /opt/miniforge/etc/profile.d/conda.sh && conda activate arsenic-m6a
make smoke

# Then full panel
bash infra/gcp/run-pipeline.sh
```

Detach tmux and close your laptop — it runs unattended.

## 6. Pull results back

From your laptop:

```bash
cd infra/gcp
./fetch-results.sh
# -> results-pulled/<timestamp>/results/  with TSVs, REPORT.md, figures, logs
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

## Tear-down checklist

```bash
gcloud compute instances list                  # what's running
gcloud compute instances delete as-m6a --zone=us-central1-a --quiet
gcloud compute disks list                      # any orphan disks?
gcloud compute addresses list                  # any reserved static IPs?
```

The instance creation here doesn't reserve a static IP or create extra disks,
so a single `instances delete` is normally enough to stop all billing.
