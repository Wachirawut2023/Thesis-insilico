# Methods explained, step-by-step

A plain-language reference for the pipeline. Companion to the thesis
chapter at `docs/chapter_in_silico.md`. Reading this end-to-end gives
you the full ability to interpret any number in the output TSVs.

---

## The big picture in one paragraph

We feed in 25 proteins (the m6A panel + four controls). For each protein,
we identify every cysteine residue, measure several geometric and
biochemical properties, then ask: *given arsenic's known chemistry, can
it physically attack this cysteine, and would attacking it disrupt the
protein's function?* The answers per cysteine are summarised into a
per-protein score and mapped to one of four qualitative tiers. The whole
computation is deterministic — same input, same output every time. No
machine learning, no probabilistic inference; it's geometric measurement
plus a transparent weighted sum.

---

## Stage-by-stage measurements

### Stage 1 — fetch (`scripts/01_fetch_structures.py`)

Downloads PDB (experimental) and AlphaFold (predicted) protein
structures from public databases. Each structure is a list of atom
coordinates in Ångström (Å) units.

A typical PDB line:
```
ATOM    123  SG  CYS A 375      10.234  -5.789  23.456  1.00 12.34   S
```
"Atom number 123, atom name SG (sulfur in cysteine side chain), residue
Cysteine, chain A, residue number 375, position (x, y, z) in Å."

### Stage 2 — prepare (`scripts/02_prepare_receptors.py`)

Strips waters and crystallization additives, deduplicates alt-locs, adds
hydrogen atoms at pH 7.4 with PROPKA. **The pH 7.4 step matters because
cysteine reactivity depends on protonation state** — thiolate (S⁻) is
~1000× more reactive with As(III) than thiol (–SH).

### Stage 3 — cys reactivity landscape (`scripts/03_cys_landscape.py`)

The **most informative single stage**. For every cysteine in every
protein, measures:

| Column | What it is | Threshold | Meaning |
|---|---|---|---|
| `sg_sasa` | Solvent-accessible surface area of the sulfur, in Å² | ≥ 5.0 Å² → "accessible" | 0 = buried; 5–10 = partial; 20+ = surface-exposed |
| `pka` | Predicted pKa of the thiol from PROPKA | ≤ 8.0 → "reactive thiolate likely" | low pKa = more often deprotonated = more reactive |
| `in_disulfide` | Sγ–Sγ distance to nearest other Cys | < 2.3 Å | Y = locked in S–S bond, not free to bind As |
| `vicinal_partners` | List of other Cys Sγ within 7 Å | n/a | informational |
| `functional_proximity_A` | Distance from Sγ to nearest curated active-site residue | ≤ 6 Å = at active site; 6–12 Å = nearby; >12 Å = far | the key biology signal |
| `reactive` | `accessible AND (pka ≤ 8 OR pka=NA) AND NOT in_disulfide` | n/a | the structural fingerprint of an As target candidate |

**Default chosen thresholds**:
- 5 Å² SASA — standard biophysical chemistry cutoff for "accessible"
- 8.0 pKa — conservative; permissive (catches more candidates)
- 2.3 Å Sγ–Sγ for disulfide — standard from crystal-structure literature
- 7 Å Sγ–Sγ for vicinal — slightly wider than the biophysical 4.4 Å bidentate window, to catch all candidates for Stage 6

### Stage 4 — pocket detection (`scripts/04_pocket_detection.py`)

Uses `fpocket` to find protein-surface cavities. Keeps only pockets that
overlap with reactive-Cys clusters, plus monodentate boxes around any
reactive Cys within 8 Å of a functional anchor.

### Stage 5 — Vina docking (`scripts/05_dock_noncovalent.py`)

Runs AutoDock Vina searching each docking box. Outputs a predicted
binding free energy ΔG in kcal/mol — more negative = stronger binding.

| ΔG (kcal/mol) | Loose interpretation |
|---|---|
| 0 to -2 | weak / no binding |
| -3 to -5 | moderate binding |
| -6 to -10 | strong (drug-like) binding |

**Caveat (stated in chapter)**: Vina is trained on organic drug-like
molecules, not inorganic arsenic. We treat ΔG here as a **relative**
ranking signal only.

### Stage 6 — geometric covalent docking (`scripts/06_dock_covalent.py`)

**The primary scientific readout.** For each reactive cysteine, scores
feasibility of three As(III) binding modes:

**Arsenic geometric facts** (from crystal structures):
- As–S bond length ≈ 2.25 Å
- S–As–S angle ≈ 94°

**Modes tested**:

| Mode | What it is | Geometric requirement |
|---|---|---|
| Monodentate | one As–S bond | reactive Cys exists |
| Bidentate | As bridges two Cys | second Sγ within 3.0–4.4 Å of anchor Sγ |
| Tridentate | As bonds to three Cys | all 3 pairwise Sγ–Sγ distances ≤ 4.4 Å |

**For bidentate**, As is placed at the apex of an isoceles triangle
bridging the two sulfurs. Heavy-atom clashes within 2.0 Å of the placed
As are counted; ≤ 2 clashes = bidentate feasible.

**Per-anchor score**:
```
score = 3.0 × tridentate_feasible
      + 2.0 × bidentate_feasible
      + func_bonus
      + 0.5 × thiolate_bonus
      - 0.5 × clash_penalty
```

Where:
- `func_bonus` = 2.0 if func_prox ≤ 6 Å, 1.0 if 6 < func_prox ≤ 12 Å, else 0
- `thiolate_bonus` = 1.0 if pKa ≤ 7.5 else 0
- `clash_penalty` = number of clashes in the best bidentate placement

**Binding mode**:
- tridentate_feasible > 0 → "tridentate"
- elif bidentate_feasible > 0 → "bidentate"
- else → "monodentate"

### Stage 7 — composite ranking (`scripts/07_score_and_rank.py`)

One row per protein. Combines all the above into a single composite score
and a four-category tier.

**Composite formula**:
```
composite = 0.6  × best_covalent_score
          + 0.25 × min(vicinal_score, 5)
          - 0.25 × best_vina_dG
          + 3.0  × prox_bonus
          + 2.0  × reactive_at_func
```

Where:
- `prox_bonus` = max(0, (6 − min_func_prox) / 6) — ranges 0 to 1
- `reactive_at_func` = 1 if min_func_prox ≤ 6 else 0

**Why each term**:
- +0.6 × covalent: reward geometric feasibility
- +0.25 × capped vicinal: modest credit for Cys clusters; cap prevents
  Cys-rich scaffolds (RBM15B 36 reactive Cys) from swamping small
  high-quality hits (PIN1 1 reactive Cys)
- −0.25 × Vina: small weight, kept for transparency
- +3.0 × prox_bonus: continuous reward for active-site proximity
- +2.0 × reactive_at_func: binary "yes there's a reactive Cys at the
  active site" — the headline biological signal

**Tier mapping**:
```
if n_reactive == 0:                       tier = no_binding
elif min_func_prox is NA:                 tier = binding_only
elif min_func_prox <= 6:                  tier = likely_inhibitory
elif min_func_prox <= 12:                 tier = possibly_allosteric
else:                                     tier = binding_only
```

---

## Worked examples

### METTL3 composite = 5.406

From the data: best_covalent_score = 2.000, vicinal_score = 1.25,
best_vina_dG = −2.472, min_func_proximity_A = 3.45.

| Term | Computation | Value |
|---|---|---|
| 0.6 × covalent_best | 0.6 × 2.000 | **+1.200** |
| 0.25 × min(vicinal, 5) | 0.25 × 1.25 | **+0.313** |
| −0.25 × vina_dG | −0.25 × (−2.472) | **+0.618** |
| 3.0 × prox_bonus | 3.0 × (6−3.45)/6 = 3.0 × 0.425 | **+1.275** |
| 2.0 × reactive_at_func | 2.0 × 1 | **+2.000** |
| **Total** | | **5.406** ✓ |

Tier: min_func_prox = 3.45 ≤ 6 → **likely_inhibitory** ✓

### METTL3 SAM-pocket cysteine covalent_score = 2.000

> These worked values use the first-pass numbers (raw 5IL0 author numbering,
> pocket cysteine reported as "Cys375", func_prox to the DPPW motif). The
> residue label (Cys375 vs Cys376) and the func_prox value are re-derived by
> Stage 0 (`00_verify_mettl3_cys.py`) and the Asp377 anchor added in
> `functional_sites.yaml`; the arithmetic of the formula is unchanged.

From the data: vicinal partner (Cys375/Cys376 pair) at 6.75 Å, func_prox
3.45 Å, pKa NA.

| Term | Value |
|---|---|
| 3.0 × tridentate | 0 (only 1 vicinal partner) |
| 2.0 × bidentate | 0 (static 6.75 Å is outside the 3.0–4.4 window) |
| func_bonus | 2.0 (3.45 ≤ 6) |
| 0.5 × thiolate | 0 (pKa NA) |
| −0.5 × clashes | 0 |
| **Total** | **2.000** ✓ |

Static binding mode: no bidentate, no tridentate → **monodentate**.

Induced-fit check (new): the covalent scorer also reports
`min_rotamer_sg_sg_A` and `bidentate_feasible_rotamer` — the minimum
clash-free Sγ–Sγ distance reachable by rotating the Cys χ1 of both partners.
If that value enters the 3.0–4.4 Å window, an As(III) **bidentate bridge**
across the Cys375/Cys376 dithiol is feasible despite the 6.75 Å static gap,
which strengthens the covalent-inhibition claim.

### TXN1 Cys32 covalent_score = 4.000

From the data: vicinal Cys35 at 3.92 Å, func_prox 0.0, pKa NA.

| Term | Value |
|---|---|
| 3.0 × tridentate | 0 |
| 2.0 × bidentate | 2.0 (3.92 Å in window, 0 clashes) |
| func_bonus | 2.0 (0 ≤ 6) |
| 0.5 × thiolate | 0 |
| −0.5 × clashes | 0 |
| **Total** | **4.000** ✓ |

Binding mode: bidentate feasible → **bidentate** ✓

This bidentate +2.0 is exactly why TXN1 ranks above METTL3 in composite:
both have a Cys at the active site, but TXN1's CGPC pair allows the
geometrically stronger bidentate mode.

---

## How to read the outputs

| File | Use for |
|---|---|
| `results/composite_ranking.tsv` | Headline ranking; thesis Table 1 |
| `results/cys_table.tsv` | Per-Cys deep dive; supplementary table |
| `results/docking_covalent.tsv` | Per-anchor feasibility; supplementary |
| `results/docking_noncovalent.tsv` | Vina ΔG values; supplementary |
| `results/pocket_table.tsv` | fpocket detail; supplementary |
| `results/REPORT.md` | Auto-generated interpretation skeleton |
| `results/figures/family_heatmap.png` | Family-level summary figure |
| `results/figures/top_<gene>.png` | Per-hit binding-pose figure |

**The interpretive workflow**:
1. Open `composite_ranking.tsv`. Look at the top — are the positive
   controls (TXN1, GAPDH, PIN1) in the top three or four?
2. If yes (the calibration check passes), look at which m6A proteins
   appear in `likely_inhibitory` and `possibly_allosteric`.
3. For each interesting hit, open `cys_table.tsv` and find the rows for
   that protein. Which cysteine is driving the signal?
4. Open `docking_covalent.tsv` for the per-anchor breakdown of that
   cysteine. What binding mode? What clashes?
5. Cross-reference with the per-protein PyMOL figure in
   `results/figures/top_<gene>.png` for visual verification.

---

## Tuning the thresholds

If you want to test the robustness of the predictions:

| Constant | Where | Effect of changing |
|---|---|---|
| `SASA_ACCESSIBLE` (default 5.0) | `scripts/03_cys_landscape.py` | Higher → fewer reactive Cys, more conservative |
| `REACTIVE_PKA` (default 8.0) | `scripts/03_cys_landscape.py` | Lower → only very acidic thiolates counted |
| `VICINAL_MAX` (default 7.0) | `scripts/03_cys_landscape.py` | Lower → fewer vicinal pairs detected |
| `BIDENTATE_RANGE` (default 3.0–4.4) | `scripts/06_dock_covalent.py` | The geometric window for As bidentate; tightening reduces false positives |
| `CLASH_TOLERANCE` (default 2) | `scripts/06_dock_covalent.py` | Lower → stricter steric filtering |
| Composite weights | `scripts/07_score_and_rank.py` | See the worked-example computation; balance shifts the ranking |
| Tier cutoffs (6 / 12 Å) | `scripts/07_score_and_rank.py` `_tier()` | Tighter cutoffs reduce `likely_inhibitory` count |

For thesis defense, the meaningful sensitivity test would be: rerun the
pipeline with composite weights ±30%, show that the top-4 and bottom-9
positions don't change. That demonstrates the ranking is robust to
weighting choices.
