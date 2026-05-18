# Hetzner Cloud setup — As(III) m6A pipeline

End-to-end runbook. Total cost for a full-panel run: **~€0.50–€3** (on-demand,
no spot). Wall time: ~5–10 hours on CCX33.

## TL;DR

```bash
# Local machine (one-time):
hcloud context create thesis            # paste API token
./provision.sh                          # creates server, prints SSH command

# Server (once cloud-init finishes, ~5–10 min):
ssh root@<ip>
tail -f /var/log/bootstrap.log          # wait for "[bootstrap] done"
cd /root/Thesis-insilico
tmux new -s pipeline
bash infra/hetzner/run-pipeline.sh      # ~5–10 hr on CCX33

# Local machine:
./fetch-results.sh <ip>
hcloud server delete as-m6a             # STOP BILLING
```

## 1. Hetzner account + API token (one-time, ~5 min)

1. Sign up at https://console.hetzner.cloud (German address not required; credit
   card or PayPal accepted).
2. Create a Project (any name — e.g. `thesis`).
3. In the project: **Security → API Tokens → Generate API Token**.
   Permission: **Read & Write**. Copy the token immediately (only shown once).

## 2. Install hcloud CLI (one-time, ~1 min)

```bash
# macOS
brew install hcloud

# Linux
curl -L https://github.com/hetznercloud/cli/releases/latest/download/hcloud-linux-amd64.tar.gz \
  | tar -xz && sudo mv hcloud /usr/local/bin/
```

Authenticate:

```bash
hcloud context create thesis
# Paste your API token when prompted.
```

## 3. SSH key (one-time)

If you don't already have one:

```bash
ssh-keygen -t ed25519 -C thesis    # accept defaults
```

`provision.sh` auto-uploads `~/.ssh/id_ed25519.pub` to your Hetzner project.

## 4. Provision the server

```bash
cd infra/hetzner
./provision.sh
```

This creates a **CCX33** server (8 dedicated AMD EPYC vCPU, 32 GB RAM, 240 GB
NVMe) in Nuremberg, with cloud-init that installs Miniforge, clones the repo,
and builds the conda env automatically.

You'll see something like:

```
Server ready at 49.12.xxx.xxx
```

Bootstrap takes 5–10 minutes after the server boots. Track progress:

```bash
ssh root@<ip>
tail -f /var/log/bootstrap.log
```

Wait for `[bootstrap] done <timestamp>` and the existence of
`/root/.bootstrap-complete`.

## 5. Run the pipeline

```bash
ssh root@<ip>
cd /root/Thesis-insilico
tmux new -s pipeline      # detach with Ctrl-b d; reattach with `tmux a`
bash infra/hetzner/run-pipeline.sh
```

Suggested workflow:

```bash
# Smoke test first (METTL3, ALKBH5, YTHDF2 — one per family). ~30 min.
make smoke

# If smoke looks sane (controls in cys_table.tsv, no script errors), run full:
bash infra/hetzner/run-pipeline.sh
```

Detach the tmux session and let it run — no babysitting needed.

## 6. Pull results back

From your laptop:

```bash
cd infra/hetzner
./fetch-results.sh <ip>
# -> results-pulled/<timestamp>/  contains TSVs, REPORT.md, figures
```

## 7. Destroy the server — STOP BILLING

```bash
hcloud server delete as-m6a
```

You will be billed for partial hours; round up. A 6-hour run on CCX33 is
6 × €0.073 = **€0.44**. Plus ~€0.50/mo if you keep a snapshot.

## Server type cheat-sheet

| Type | vCPU | RAM | Disk | €/hr | Recommended for |
|---|---|---|---|---|---|
| CPX31 | 4 shared | 8 GB | 160 GB | 0.026 | tiny smoke test only |
| **CCX33** | **8 dedicated** | **32 GB** | **240 GB** | **0.073** | **default — full panel** |
| CCX43 | 16 dedicated | 64 GB | 360 GB | 0.146 | 2× wall-time speedup if in a hurry |

Hetzner has **no GPU offering** for CCX/CPX lines. Since this pipeline is
CPU-only, that's fine. If you ever add MD (Tier 3), switch to a GPU provider
(Lambda Labs, RunPod, Vast.ai) or AWS `g5.xlarge` for that step alone.

## Cost / runtime estimates

Pipeline is parallelised across proteins; Vina internally uses all cores.

| Stage | CCX33 wall time |
|---|---|
| fetch | 5–10 min (network-bound) |
| prep | 10–20 min (PROPKA per protein) |
| cys | <5 min |
| pockets | 10–15 min (fpocket) |
| dock-nc | 2–4 hr (Vina × boxes) |
| dock-cov | <5 min (geometric, fast) |
| rank | <1 min |
| **total** | **~3–6 hr** → **€0.22–€0.44** |

## Common gotchas

- **Region matters for cost only after egress.** Hetzner doesn't charge per-GB
  egress within the included 20 TB / month. Any location is fine.
- **`mamba env create` slow on first run** (~5 min). Subsequent updates are fast.
- **Vina ligand atom type warnings.** Expected — As is handled via the
  `Mg`-substitution hack in `prepare_as3.sh`. Treat Vina ΔG as relative ranking.
- **`pdb2pqr30` warnings about chain breaks** are also expected for AlphaFold
  models with disordered termini. Falls back to `obabel -p 7.4` automatically.
- **Disk space.** 240 GB on CCX33 is overkill — full pipeline outputs are <2 GB.

## Tear-down checklist

```bash
hcloud server list                    # confirm what's running
hcloud server delete as-m6a           # destroy server
hcloud ssh-key list                   # optional: clean up
hcloud context delete thesis          # optional
```

Without explicit deletion, **the server keeps billing**. Set a calendar
reminder if you're stepping away.
