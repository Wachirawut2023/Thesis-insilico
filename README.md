# Thesis-insilico — Arsenite (As(III)) direct inhibition of m6A machinery

In silico scoping study for the thesis question:

> Does arsenic disrupt hematopoietic stem cell differentiation by directly
> inhibiting m6A writers, erasers, and readers?

A structure-based screen ranks ~18 m6A proteins (plus controls) for biophysical
plausibility as direct As(III)/arsenite targets, to scope wet-lab assays.

See `/root/.claude/plans/your-are-expert-of-prancy-quiche.md` for the full
plan and rationale.

## Quick start

```bash
# 1. Create environment
make env
conda activate arsenic-m6a

# 2. Prepare the As(OH)3 ligand
make ligand

# 3. Smoke test on 3 targets (one writer, one eraser, one reader)
make smoke

# 4. Full run (~18 targets + 4 controls)
make all
```

Outputs land in `results/`:

| File | Stage | Content |
|---|---|---|
| `cys_table.tsv` | 2 | Per-Cys reactivity (SASA, pKa, vicinal partners, functional proximity) |
| `protein_summary.tsv` | 2 | Per-protein Cys landscape summary |
| `pocket_table.tsv` | 3 | fpocket predictions |
| `docking_boxes.tsv` | 3 | Pocket × Cys-cluster docking boxes |
| `docking_noncovalent.tsv` | 4 | Vina ΔG (relative ranking only) |
| `docking_covalent.tsv` | 5 | Geometric covalent As(III)-Cys feasibility |
| `composite_ranking.tsv` | 6 | Final ranked list |
| `figures/family_heatmap.png` | 6 | Family-level susceptibility heatmap |
| `REPORT.md` | 6 | Interpretation skeleton |

## Pipeline stages

```
00_verify_mettl3_cys.py  METTL3 Cys375-vs-Cys376 numbering + geometry diagnostic (run after fetch)
01_fetch_structures.py   PDB + AlphaFold downloads
02_prepare_receptors.py  Strip, protonate (pH 7.4, PROPKA), PDBQT
03_cys_landscape.py      Per-Cys SASA, pKa, disulfide, vicinal pairs/triads
04_pocket_detection.py   fpocket; intersect with Cys clusters -> docking boxes
05_dock_noncovalent.py   AutoDock Vina across boxes
06_dock_covalent.py      Geometric covalent docking (primary readout; incl. induced-fit bidentate)
07_score_and_rank.py     Composite score, heatmap, REPORT.md
```

> **METTL3 Cys375 vs Cys376.** The pipeline copies residue numbers verbatim from
> the deposited structure (no renumbering), so the first-pass reactive cysteine
> was reported as "Cys375" in raw 5IL0 author numbering, whereas the literature
> SAM-pocket cysteine (adjacent to the SAM-binding Asp377) is **Cys376**. Run
> `make verify` (Stage 0) after `make fetch` to reconcile the author↔UniProt
> numbering and measure each candidate Sγ's distance to the SAM/SAH cofactor and
> Asp377 — see `docs/chapter_in_silico.md` §5.2.

## What this pipeline does *not* prove

- It proves *binding plausibility*, not enzymatic *inhibition*. Inhibition is
  inferred from binding × functional proximity (see `inference_tier` column in
  `composite_ranking.tsv`).
- Vina scoring of arsenic is approximate (non-standard ligand); ΔG values are
  relative, not absolute affinities.
- Static structures: cryptic Cys clusters that open in MD are missed.
- Intracellular GSH (~mM) competes for As(III) — kinetic reachability requires
  experimental confirmation.

## Targets

See `data/targets.tsv`. Panel: METTL3, METTL14, WTAP, VIRMA, ZC3H13, RBM15(B),
CBLL1, METTL16 (writers); FTO, ALKBH5 (erasers); YTHDF1/2/3, YTHDC1/2,
IGF2BP1/2/3, HNRNPA2B1, HNRNPC (readers); plus PIN1, TXN1, GAPDH (positive
controls — known As(III) targets) and SH3-SRC (negative control).

## Compute

CPU-only workload. Two provisioning paths are wired up:

**GCP (recommended — covered by $300 free trial)** — `n2d-standard-8`,
8 AMD EPYC vCPU, 32 GB RAM, $0.388/hr → ~$2.50 / run, ~3–6 hr wall time.
Runbook: [`infra/gcp/README.md`](infra/gcp/README.md).

```bash
gcloud auth login && gcloud config set project YOUR_PROJECT_ID
cd infra/gcp
./provision.sh
gcloud compute ssh as-m6a --zone=us-central1-a
# then on the server: cd /opt/Thesis-insilico && bash infra/gcp/run-pipeline.sh
./fetch-results.sh
gcloud compute instances delete as-m6a --zone=us-central1-a --quiet
```

**Hetzner (cheapest list price, no free credit)** — `CCX33`, 8 dedicated
AMD EPYC vCPU, €0.073/hr → ~€0.50 / run. Runbook:
[`infra/hetzner/README.md`](infra/hetzner/README.md).

**GitHub Actions (free, no cloud account needed)** — Tier 1-2 (stages 00-08)
can also run on GitHub-hosted runners via
[`.github/workflows/tier1-2-pipeline.yml`](.github/workflows/tier1-2-pipeline.yml).
Trigger it manually from the Actions tab (`workflow_dispatch`), starting with
`target_set: smoke` to validate the toolchain, then `target_set: all` for the
full panel. Stages 03-06 checkpoint per-gene and results/ is committed back to
the branch after every run, so — since a single GitHub-hosted job is capped at
6h — the full run may need a few re-dispatches to finish; each one resumes
from the previous run's checkpoints instead of restarting. Once Tier 1-2's
results are on the branch, only Tier 3 (MD, see `docs/md_tier3.md`) needs a
GPU instance.
