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
01_fetch_structures.py   PDB + AlphaFold downloads
02_prepare_receptors.py  Strip, protonate (pH 7.4, PROPKA), PDBQT
03_cys_landscape.py      Per-Cys SASA, pKa, disulfide, vicinal pairs/triads
04_pocket_detection.py   fpocket; intersect with Cys clusters -> docking boxes
05_dock_noncovalent.py   AutoDock Vina across boxes
06_dock_covalent.py      Geometric covalent docking (primary readout)
07_score_and_rank.py     Composite score, heatmap, REPORT.md
```

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

CPU-only workload. Cheapest option is **Hetzner Cloud CCX33** (8 dedicated
AMD EPYC vCPU, 32 GB RAM) at €0.073/hr → total cost **~€0.50–€3** for the
full panel, ~3–6 hr wall time.

Full provisioning runbook in [`infra/hetzner/README.md`](infra/hetzner/README.md):

```bash
cd infra/hetzner
./provision.sh                 # creates server with cloud-init bootstrap
ssh root@<ip>                  # then: cd /root/Thesis-insilico && make all
./fetch-results.sh <ip>        # pull TSVs/figures back to laptop
hcloud server delete as-m6a    # stop billing
```

AWS equivalent: `c7i.4xlarge` on-demand ~$50–$80; spot ~$20–$30.
