# As(III) direct inhibition of m6A machinery — in silico scoping report

**Methodology.** Tier 1 (Cys reactivity landscape) + Tier 2 (geometric covalent docking + Vina non-covalent).
Inhibition is *inferred* from binding plausibility × functional proximity; in vitro confirmation required.

## Top 5 candidates (any family)

| Gene | Family | Composite | Best covalent | Mode | Min func.prox (Å) | Tier |
|---|---|---|---|---|---|---|
| TXN1 | control_pos | 8.708 | 4.000 | bidentate | 0.00 | likely_inhibitory |
| GAPDH | control_pos | 7.516 | 2.000 | monodentate | 0.00 | likely_inhibitory |
| PIN1 | control_pos | 7.132 | 2.000 | monodentate | 0.00 | likely_inhibitory |
| METTL3 | writer | 5.404 | 2.000 | monodentate | 3.45 | likely_inhibitory |
| VIRMA | writer | 3.154 | 2.000 | bidentate | NA | binding_only |

## Writers

- **METTL3**: composite=5.404, covalent=2.000 (monodentate), reactive Cys=1, func.prox=3.45 Å — *likely_inhibitory*
- **VIRMA**: composite=3.154, covalent=2.000 (bidentate), reactive Cys=10, func.prox=NA Å — *binding_only*
- **ZC3H13**: composite=2.904, covalent=2.000 (bidentate), reactive Cys=9, func.prox=NA Å — *binding_only*
- **CBLL1**: composite=2.676, covalent=2.500 (bidentate), reactive Cys=3, func.prox=8.68 Å — *possibly_allosteric*
- **RBM15B**: composite=1.994, covalent=0.000 (monodentate), reactive Cys=36, func.prox=NA Å — *binding_only*
- **METTL14**: composite=0.931, covalent=0.000 (monodentate), reactive Cys=1, func.prox=34.90 Å — *binding_only*
- **METTL16**: composite=0.85, covalent=1.000 (monodentate), reactive Cys=4, func.prox=9.21 Å — *possibly_allosteric*
- **RBM15**: composite=0.188, covalent=0.000 (monodentate), reactive Cys=3, func.prox=NA Å — *binding_only*
- **WTAP**: composite=0.0, covalent=NA (NA), reactive Cys=0, func.prox=NA Å — *no_binding*

## Erasers

- **FTO**: composite=3.126, covalent=3.000 (bidentate), reactive Cys=4, func.prox=9.18 Å — *possibly_allosteric*
- **ALKBH5**: composite=0.725, covalent=1.000 (monodentate), reactive Cys=2, func.prox=9.40 Å — *possibly_allosteric*

## Readers

- **YTHDF1**: composite=0.062, covalent=0.000 (monodentate), reactive Cys=1, func.prox=18.66 Å — *binding_only*
- **IGF2BP3**: composite=0.062, covalent=0.000 (monodentate), reactive Cys=1, func.prox=34.55 Å — *binding_only*
- **YTHDF2**: composite=0.0, covalent=NA (NA), reactive Cys=0, func.prox=NA Å — *no_binding*
- **YTHDF3**: composite=0.0, covalent=NA (NA), reactive Cys=0, func.prox=NA Å — *no_binding*
- **YTHDC1**: composite=0.0, covalent=NA (NA), reactive Cys=0, func.prox=NA Å — *no_binding*
- **YTHDC2**: composite=0.0, covalent=NA (NA), reactive Cys=0, func.prox=NA Å — *no_binding*
- **IGF2BP1**: composite=0.0, covalent=NA (NA), reactive Cys=0, func.prox=NA Å — *no_binding*
- **IGF2BP2**: composite=0.0, covalent=NA (NA), reactive Cys=0, func.prox=NA Å — *no_binding*
- **HNRNPA2B1**: composite=0.0, covalent=NA (NA), reactive Cys=0, func.prox=NA Å — *no_binding*
- **HNRNPC**: composite=0.0, covalent=NA (NA), reactive Cys=0, func.prox=NA Å — *no_binding*

## Controls

- TXN1 (control_pos): composite=8.708, tier=likely_inhibitory
- GAPDH (control_pos): composite=7.516, tier=likely_inhibitory
- PIN1 (control_pos): composite=7.132, tier=likely_inhibitory
- SH3_SRC (control_neg): composite=0.125, tier=binding_only

## Caveats

- Vina ΔG is *relative* — As is non-standard ligand, scoring not calibrated for inorganic.
- Covalent score is a geometric proxy; not a free energy.
- Static structures: cryptic Cys clusters revealed by dynamics are missed.
- Intracellular GSH (~1–10 mM) competes; reachable Cys may never see arsenite in vivo.
- Inhibition vs binding: only confirmed by in vitro activity assay (± DTT/GSH rescue).