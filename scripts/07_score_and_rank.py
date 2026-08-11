"""Stage 6: composite ranking + interpretation deliverable.

Combines per-protein outputs from:
  - results/protein_summary.tsv  (stage 2: vicinal_score)
  - results/docking_noncovalent.tsv (stage 4: best Vina ΔG)
  - results/docking_covalent.tsv (stage 5: best covalent score)

into a single composite inhibition-likelihood score, family heatmap (PNG),
and skeleton REPORT.md.
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

from _common import RESULTS, get_logger, load_targets, write_tsv

LOG = get_logger("rank")
SUMMARY_TSV = RESULTS / "protein_summary.tsv"
NC_TSV = RESULTS / "docking_noncovalent.tsv"
COV_TSV = RESULTS / "docking_covalent.tsv"
RANK_TSV = RESULTS / "composite_ranking.tsv"
HEATMAP_PNG = RESULTS / "figures" / "family_heatmap.png"
REPORT_MD = RESULTS / "REPORT.md"

W_COV = 0.6           # covalent score weight (down from 1.0; sign too volatile)
W_VICINAL = 0.25      # vicinal cluster score weight (down from 0.4; capped below)
VICINAL_CAP = 5.0     # avoid letting Cys-rich scaffolds dominate
W_VINA = -0.25
W_PROX_BONUS = 3.0    # functional proximity is the biology-relevant signal (up from 0.5)
W_REACTIVE_AT_FUNC = 2.0   # explicit bonus for "reactive Cys at active site"


def _load_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _safe_float(s: str) -> float | None:
    if s in (None, "", "NA"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def main() -> int:
    targets = {t.gene: t for t in load_targets()}
    summary = {r["gene"]: r for r in _load_tsv(SUMMARY_TSV)}
    nc = _load_tsv(NC_TSV)
    cov = _load_tsv(COV_TSV)

    best_vina: dict[str, float] = {}
    for r in nc:
        v = _safe_float(r["vina_dG_kcal_mol"])
        if v is None:
            continue
        if r["gene"] not in best_vina or v < best_vina[r["gene"]]:
            best_vina[r["gene"]] = v

    best_cov: dict[str, float] = defaultdict(lambda: float("-inf"))
    min_func_prox: dict[str, float] = {}
    sum_cov: dict[str, float] = defaultdict(float)
    n_cov: dict[str, int] = defaultdict(int)
    best_mode: dict[str, str] = {}
    for r in cov:
        s = _safe_float(r["covalent_score"])
        if s is None:
            continue
        if s > best_cov[r["gene"]]:
            best_cov[r["gene"]] = s
            best_mode[r["gene"]] = r.get("binding_mode", "")
        sum_cov[r["gene"]] += s
        n_cov[r["gene"]] += 1
        fp = _safe_float(r["functional_proximity_A"])
        if fp is not None:
            if r["gene"] not in min_func_prox or fp < min_func_prox[r["gene"]]:
                min_func_prox[r["gene"]] = fp

    rows: list[dict] = []
    for gene, t in targets.items():
        s = summary.get(gene, {})
        cov_best = best_cov.get(gene)
        if cov_best == float("-inf"):
            cov_best = None
        vina = best_vina.get(gene)
        vicinal = _safe_float(s.get("vicinal_score", "0")) or 0.0
        prox = min_func_prox.get(gene)
        prox_bonus = 0.0
        if prox is not None and prox <= 6.0:
            prox_bonus = (6.0 - prox) / 6.0  # 0..1
        vicinal_capped = min(vicinal, VICINAL_CAP)
        reactive_at_func = 1.0 if (prox is not None and prox <= 6.0) else 0.0
        composite = (
            W_COV * (cov_best or 0.0)
            + W_VICINAL * vicinal_capped
            + (W_VINA * vina if vina is not None else 0.0)
            + W_PROX_BONUS * prox_bonus
            + W_REACTIVE_AT_FUNC * reactive_at_func
        )
        rows.append(
            {
                "gene": gene,
                "family": t.family,
                "subfamily": t.subfamily,
                "n_cys": s.get("n_cys", ""),
                "n_reactive": s.get("n_reactive", ""),
                "n_vicinal_pairs": s.get("n_vicinal_pairs", ""),
                "vicinal_score": s.get("vicinal_score", ""),
                "best_vina_dG": f"{vina:.3f}" if vina is not None else "NA",
                "best_covalent_score": f"{cov_best:.3f}" if cov_best is not None else "NA",
                "best_binding_mode": best_mode.get(gene, "NA"),
                "min_func_proximity_A": f"{prox:.2f}" if prox is not None else "NA",
                "composite_score": round(composite, 3),
                "inference_tier": _tier(prox, cov_best, int(s.get("n_reactive", 0) or 0)),
            }
        )

    rows.sort(key=lambda r: r["composite_score"], reverse=True)
    write_tsv(
        rows,
        RANK_TSV,
        [
            "gene",
            "family",
            "subfamily",
            "n_cys",
            "n_reactive",
            "n_vicinal_pairs",
            "vicinal_score",
            "best_vina_dG",
            "best_covalent_score",
            "best_binding_mode",
            "min_func_proximity_A",
            "composite_score",
            "inference_tier",
        ],
    )
    LOG.info("wrote %s", RANK_TSV)

    _write_heatmap(rows)
    _write_report(rows)
    return 0


def _tier(prox: float | None, cov: float | None, n_reactive: int = 0) -> str:
    """Qualitative inference tier.

    Anchored on functional proximity (the biologically interpretable signal)
    rather than covalent_score sign (which is volatile near zero — a slight
    func-proximity penalty pushes it negative even when binding is plausible).

      likely_inhibitory:    reactive Cys exists with Sγ ≤ 6 Å from active site
      possibly_allosteric:  reactive Cys exists 6 < Sγ ≤ 12 Å from active site
      binding_only:         reactive Cys present but none near a curated anchor
                            (either no anchor curated, or all far from active)
      no_binding:           no reactive Cys at all
    """
    if n_reactive <= 0:
        return "no_binding"
    if prox is None:
        return "binding_only"
    if prox <= 6.0:
        return "likely_inhibitory"
    if prox <= 12.0:
        return "possibly_allosteric"
    return "binding_only"


def _write_heatmap(rows: list[dict]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        LOG.warning("matplotlib not available; skipping heatmap")
        return
    HEATMAP_PNG.parent.mkdir(parents=True, exist_ok=True)
    families = ["writer", "eraser", "reader", "control_pos", "control_neg"]
    by_family: dict[str, list[tuple[str, float]]] = {f: [] for f in families}
    for r in rows:
        if r["family"] in by_family:
            by_family[r["family"]].append((r["gene"], r["composite_score"]))
    max_len = max((len(v) for v in by_family.values()), default=1)
    if max_len == 0:
        return
    grid = np.full((len(families), max_len), np.nan)
    labels: list[list[str]] = []
    for i, fam in enumerate(families):
        items = sorted(by_family[fam], key=lambda x: x[1], reverse=True)
        row_labels: list[str] = []
        for j, (g, s) in enumerate(items):
            grid[i, j] = s
            row_labels.append(g)
        while len(row_labels) < max_len:
            row_labels.append("")
        labels.append(row_labels)
    fig, ax = plt.subplots(figsize=(max(8, max_len * 0.7), 4))
    im = ax.imshow(grid, aspect="auto", cmap="viridis")
    ax.set_yticks(range(len(families)))
    ax.set_yticklabels(families)
    ax.set_xticks([])
    for i in range(len(families)):
        for j in range(max_len):
            if labels[i][j]:
                ax.text(j, i, labels[i][j], ha="center", va="center", color="white", fontsize=8)
    fig.colorbar(im, ax=ax, label="composite score")
    ax.set_title("As(III) direct-inhibition likelihood across m6A families")
    fig.tight_layout()
    fig.savefig(HEATMAP_PNG, dpi=150)
    plt.close(fig)
    LOG.info("wrote %s", HEATMAP_PNG)


def _write_report(rows: list[dict]) -> None:
    top = rows[:5]
    by_family: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_family[r["family"]].append(r)
    out: list[str] = []
    out.append("# As(III) direct inhibition of m6A machinery — in silico scoping report")
    out.append("")
    out.append("**Methodology.** Tier 1 (Cys reactivity landscape) + Tier 2 (geometric covalent docking + Vina non-covalent).")
    out.append("Inhibition is *inferred* from binding plausibility × functional proximity; in vitro confirmation required.")
    out.append("")
    out.append("## Top 5 candidates (any family)")
    out.append("")
    out.append("| Gene | Family | Composite | Best covalent | Mode | Min func.prox (Å) | Tier |")
    out.append("|---|---|---|---|---|---|---|")
    for r in top:
        out.append(
            f"| {r['gene']} | {r['family']} | {r['composite_score']} | "
            f"{r['best_covalent_score']} | {r['best_binding_mode']} | "
            f"{r['min_func_proximity_A']} | {r['inference_tier']} |"
        )
    out.append("")
    for fam in ("writer", "eraser", "reader"):
        out.append(f"## {fam.capitalize()}s")
        out.append("")
        for r in by_family.get(fam, []):
            out.append(
                f"- **{r['gene']}**: composite={r['composite_score']}, "
                f"covalent={r['best_covalent_score']} ({r['best_binding_mode']}), "
                f"reactive Cys={r['n_reactive']}, "
                f"func.prox={r['min_func_proximity_A']} Å — *{r['inference_tier']}*"
            )
        out.append("")
    out.append("## Controls")
    out.append("")
    for fam in ("control_pos", "control_neg"):
        for r in by_family.get(fam, []):
            out.append(
                f"- {r['gene']} ({fam}): composite={r['composite_score']}, tier={r['inference_tier']}"
            )
    out.append("")
    out.append("## Caveats")
    out.append("")
    out.append("- Vina ΔG is *relative* — As is non-standard ligand, scoring not calibrated for inorganic.")
    out.append("- Covalent score is a geometric proxy; not a free energy.")
    out.append("- Static structures: cryptic Cys clusters revealed by dynamics are missed.")
    out.append("- Intracellular GSH (~1–10 mM) competes; reachable Cys may never see arsenite in vivo.")
    out.append("- Inhibition vs binding: only confirmed by in vitro activity assay (± DTT/GSH rescue).")
    REPORT_MD.write_text("\n".join(out))
    LOG.info("wrote %s", REPORT_MD)


if __name__ == "__main__":
    sys.exit(main())
