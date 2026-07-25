# Tier 3 — Molecular Dynamics methodology

Supplementary methods to be merged into the in silico chapter (likely as
a new §10 or as a sub-section of §6 Discussion).

## What MD adds beyond the static screen

The static pipeline (chapters 1–9) ranks proteins by predicted
arsenic-binding plausibility based on **a single conformation** of each
protein. MD addresses two specific limitations of that approach:

1. **Cryptic cysteines** — buried cysteines that may become accessible
   during normal protein motion. The static SASA is a snapshot; the MD
   trajectory gives the *distribution* of SASA over time.
2. **Binding-pose stability** — does the protein conformation tolerate
   arsenic at the predicted cysteine, or does it distort / dissociate
   in dynamics?

Two MD modes are run per protein:

### A. Apo MD (50 ns)
Protein in explicit water (OPC — the water model ff19SB was
parameterised and validated against; GROMACS's own force-field docs
recommend OPC/OPC3 over TIP3P for this force field), 150 mM NaCl,
AMBER ff19SB force field. Standard biomolecular protocol: energy
minimisation → 100 ps NVT
equilibration at 310 K (position restraints on heavy atoms) → 100 ps
NPT equilibration at 1 bar → 50 ns production with 2 fs timestep,
LINCS constraints on hydrogen bonds, PME electrostatics.

**Analyses produced**:
- Cα RMSD over time (overall protein stability)
- Per-residue RMSF (mobility map)
- Per-cysteine Sγ SASA over time (cryptic-Cys detection)
- Vicinal-Cys distance over time (persistence of bidentate-capable pairs)

### B. As-bound MD (50 ns)
Same protocol as apo, but with a Mg²⁺ surrogate placed at the predicted
As(III) bonding position and harmonic distance restraints maintaining
Mg²⁺–SG distance at 2.25 ± 0.02 Å (the As–S bond length from crystal
structures).

**Why Mg²⁺ and not As(III) directly?** As(III) is not parameterised in
AMBER ff19SB or any standard biomolecular force field. The two paths
to "real" As MD are:

1. Derive partial atomic charges via QM (e.g. RESP fitting from
   Gaussian B3LYP/6-31G* calculations on As-cysteinate model
   complexes), then assign bonded parameters from analogous metalloid
   systems. ~1–2 days of QM/setup work *per protein*.
2. Use MCPB.py from AmberTools to generate parameters for As-bonded
   residues automatically. Faster but still requires QM input and
   custom validation.

For this thesis, the Mg²⁺ surrogate with distance restraints is the
defensible pragmatic choice because:

- Mg²⁺ has a similar ionic radius to As(III) (0.72 Å vs 0.58 Å).
- Both have +0 to +3 nominal charge ranges in their protein-bound forms.
- The distance restraint enforces the As–S geometry without needing
  custom bonded parameters.
- The question we want MD to answer ("does the protein binding pocket
  tolerate occupation by a metal-like atom at the predicted As
  geometry?") is exactly what this setup tests.
- The simulation does *not* claim to compute absolute As binding
  energies — those would need full QM/MM or the parameterisation work
  above. The scientific claim is restricted to **pose stability and
  protein conformational response**.

This caveat is explicit in the Methods section.

## What the analyses tell us

### RMSD plateau
- A protein whose Cα RMSD plateaus below 2 Å is stable in its starting
  fold. Drift above 3–4 Å suggests partial unfolding (uninteresting
  artefact) or a real conformational rearrangement (potentially
  interesting).
- For As-bound vs apo: if As-bound RMSD is markedly higher than apo,
  the predicted binding causes structural perturbation — possibly
  indicating the binding is destabilising rather than stable.

### RMSF map
- High RMSF (> 3 Å) regions are floppy loops/termini — usually not
  functional.
- Catalytic-site residues should be low-RMSF (~1 Å). If As-bound RMSF
  at the active site spikes, that's evidence of allosteric perturbation.

### Per-Cys SASA over time
- The headline analysis. Each reactive Cys identified in the static
  screen should remain accessible (mean SASA > 5 Å²) throughout the
  MD trajectory. If it goes buried (mean < 5) most of the simulation,
  the static prediction was a snapshot artefact.
- **Cryptic-Cys flag**: any Cys whose minimum SASA < 2 Å² (buried)
  but maximum SASA > 10 Å² (transiently exposed) is flagged. These are
  candidates that the static screen missed.

### Vicinal-Cys distance over time
- For bidentate hits (TXN1 Cys32-Cys35 at 3.92 Å in static): track
  the Sγ–Sγ distance over the trajectory. If it stays in the 3.0–4.4 Å
  bidentate window, the bidentate prediction is dynamically supported.
  If it drifts above 5–6 Å, the pair only fleetingly satisfies the As
  geometry.

## Interpretation matrix

For each protein, after both apo + bound MD complete:

| Apo result | Bound result | Interpretation |
|---|---|---|
| Cys stays accessible (mean SASA > 5) | Mg²⁺ stays bound (distance restraint satisfied), low ΔRMSD apo→bound | **Confirmed binding pose, stable.** Strongest result. |
| Cys stays accessible | Mg²⁺ pose tolerated but RMSD higher | **Binding plausible but causes conformational change.** Possibly allosteric inhibition. |
| Cys becomes buried (mean SASA < 5) | n/a (binding probably can't form) | **Static prediction overcalled.** Demote from likely_inhibitory tier. |
| Cryptic Cys flagged (was buried, opens transiently) | n/a | **New candidate target** the static screen missed. |
| Cys stays accessible | Large protein distortion in bound | **Geometric possibility but biological cost.** Mark as ambiguous; needs in vitro test. |

## Where the results plug into the thesis

In `docs/chapter_in_silico.md`:

- §3.7 (limits) — strike "static structures only" from limitations once
  Tier 3 is done.
- §5 (results) — add per-protein paragraph: "MD over 50 ns confirms
  Cys375 remains solvent-accessible (mean SASA = ... Å², range ... to
  ...) and the bound-state RMSD ... ."
- §6 (discussion) — strengthen the METTL3 / FTO inhibition story with
  pose-stability evidence.
- New §10 — "Molecular dynamics validation of top hits" — summary of
  the apo + bound analyses with the interpretation matrix above.

## Cost & runtime

**Local (your own computer, CPU-only)** — no cloud provisioning, $0 cost.
Most of Tier 3's wall time is the production `mdrun` step, and GROMACS's
CPU-only path (Verlet scheme, thread-MPI) runs the exact same protocol as
the GPU paths below — just slower. Reasonable if you don't mind the
simulations running for days instead of hours, or you're only running a
handful of genes rather than the full top-8 panel.

```bash
conda activate arsenic-m6a
conda install -c conda-forge -c bioconda "gromacs=2024" pdbfixer
bash scripts/md/run_all.sh          # or run_md.sh for a single gene/mode
```

`run_md.sh` auto-detects the absence of an NVIDIA GPU (`MD_DEVICE=auto`,
the default) and drops the `-nb gpu -pme gpu` offload flags, running
everything on CPU threads instead (`-nt $(nproc)` by default; override
with `NT=<n>`). Force a path explicitly with `MD_DEVICE=cpu` or
`MD_DEVICE=gpu` if you have an NVIDIA GPU locally and want to use it (or
confirm it's actually being picked up). No infra runbook needed — this
runs directly against your local `data/prepared/` and writes to
`results/md/` like any other pipeline stage.

Three cloud provisioning paths are also supported, if you'd rather trade
money for wall time:

**GCP `g2-standard-8` with NVIDIA L4** (recommended — same project as
static pipeline):

| Resource | Specification | Cost |
|---|---|---|
| Compute | g2-standard-8 + L4 | $0.82/hr |
| Wall time | 8 × 2 × ~1.5–2 hr each | ~24–32 GPU-hr |
| Total compute | | ~$20–26 |
| Storage + IP + egress | | ~$1 |
| **Total** | | **~$25–28** (or $0 with GCP free trial credit) |

Runbook: `infra/gcp-gpu/README.md`. Quota request needed first (most
projects have 0 GPU quota by default — see step 0 of the runbook).

**RunPod RTX 4090** (alternative — cheaper but separate account):

| Resource | Specification | Cost |
|---|---|---|
| Compute | RunPod RTX 4090 | $0.40/hr |
| Total GPU time | 24–32 GPU-hr | ~$10–13 |
| Storage + buffer | | ~$5 |
| **Total** | | **~$15–20** |

Runbook: `infra/runpod/README.md`. New account/CLI but ~$10 cheaper.

**DigitalOcean AMD GPU Droplet (MI300X)** (alternative — used to replace
a Google Colab run that hit session-length/disconnect limits mid-panel):

GROMACS's GPU offload path works unchanged on AMD; the only difference
is GROMACS is compiled from source with HIP support against the
droplet's ROCm stack instead of installed as a prebuilt CUDA conda
package. Cost depends on current DO GPU Droplet pricing for MI300X.

Runbook: `infra/amd-gpu/README.md`.

**Known issue — production-step OOM on long runs**: on the MI300X
droplet, `gmx mdrun`'s production step has OOM-killed 3 separate times
(different genes/system sizes), always late in the run (83-91% of
steps) at ~176GB RSS — consistent with an in-process memory leak in
this GROMACS 2026.3 build tied to simulation progress rather than a
one-off hardware hiccup. `run_md.sh` now runs production in
`PROD_MAXH_HOURS`-long (default 2h) segments via `-maxh` + checkpoint
resume so no single process runs long enough to reach the failure
point — same output trajectory, just restarted periodically.
