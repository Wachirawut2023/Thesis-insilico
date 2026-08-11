"""Stage 8: PyMOL renders of top inhibition candidates.

For each of the top N proteins by composite score, produces:
  - figures/top_<gene>.pse   — PyMOL session (open locally to rotate/zoom)
  - figures/top_<gene>.png   — static ray-traced render (for thesis chapter)
  - figures/top_<gene>.pml   — generating script (for reproducibility)

Highlights:
  yellow sticks   reactive anchor Cys (the best-scoring covalent anchor)
  orange sticks   vicinal Cys partner(s) for bidentate/tridentate cases
  cyan sticks     curated functional anchor residues (catalytic / RNA-binding)
  gray cartoon    rest of the protein at 30% transparency

Headless render via `pymol -cq <script.pml>` — needs pymol-open-source
from the conda env.
"""
from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

from _common import DATA, PREP_DIR, RESULTS, get_logger

LOG = get_logger("viz")
TOP_N = 5
RANK_TSV = RESULTS / "composite_ranking.tsv"
COV_TSV = RESULTS / "docking_covalent.tsv"
SITES_YAML = DATA / "functional_sites.yaml"
FIGURES = RESULTS / "figures"


def _build_pml(gene: str, anchor: dict, sites_cfg: dict, pdb: Path, pse_out: Path, png_out: Path) -> str:
    anchor_chain = anchor["anchor_chain"]
    anchor_resi = anchor["anchor_resi"]
    partner_tags = (anchor.get("partner_tags") or "").strip()
    site_anchors = (sites_cfg.get(gene) or {}).get("anchors", [])
    binding_mode = anchor.get("binding_mode", "")

    lines = [
        f"# Top hit: {gene}",
        f"# anchor: chain {anchor_chain} Cys{anchor_resi}",
        f"# binding_mode: {binding_mode}",
        f"# covalent_score: {anchor.get('covalent_score', '')}",
        f"# functional_proximity_A: {anchor.get('functional_proximity_A', '')}",
        "",
        f"load {pdb}, prot",
        "bg_color white",
        "hide everything",
        "show cartoon, prot",
        "color gray80, prot",
        "set cartoon_transparency, 0.3",
        "",
        "# Reactive anchor Cys (yellow)",
        f"select anchor_cys, chain {anchor_chain} and resi {anchor_resi} and resn CYS",
        "show sticks, anchor_cys",
        "color yellow, anchor_cys",
        "label first (anchor_cys and name CA), \"C%s\" % resi",
        "",
    ]

    if partner_tags:
        lines.append("# Bidentate/tridentate partner Cys (orange)")
        for tag in partner_tags.split(";"):
            if not tag:
                continue
            p_chain = tag[0]
            p_resi = tag[1:]
            sel = f"partner_{p_chain}{p_resi}"
            lines.extend([
                f"select {sel}, chain {p_chain} and resi {p_resi} and resn CYS",
                f"show sticks, {sel}",
                f"color orange, {sel}",
                f"label first ({sel} and name CA), \"C%s\" % resi",
            ])
        lines.append("")

    if site_anchors:
        lines.append("# Functional anchor residues (cyan)")
        for site in site_anchors[:20]:
            if isinstance(site, int):
                resi = site
                sel_name = f"funcanc_{resi}"
                lines.extend([
                    f"select {sel_name}, resi {resi} and (not anchor_cys)",
                    f"show sticks, {sel_name}",
                    f"color cyan, {sel_name}",
                ])
            elif isinstance(site, dict):
                resi = site.get("resi")
                chain = site.get("chain", "")
                if resi is None:
                    continue
                sel_name = f"funcanc_{resi}"
                sel = f"resi {resi}"
                if chain:
                    sel = f"chain {chain} and " + sel
                lines.extend([
                    f"select {sel_name}, ({sel}) and (not anchor_cys)",
                    f"show sticks, {sel_name}",
                    f"color cyan, {sel_name}",
                ])
        lines.append("")

    lines.extend([
        "# Camera + render",
        "orient anchor_cys",
        "zoom anchor_cys, 12",
        "set ray_shadows, 0",
        "set ray_opaque_background, 1",
        "set label_size, 14",
        "set label_color, black",
        f"ray 1200, 900",
        f"png {png_out}, dpi=150",
        f"save {pse_out}",
    ])
    return "\n".join(lines) + "\n"


def _best_anchor_per_gene(rows: list[dict]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for r in rows:
        try:
            s = float(r["covalent_score"])
        except (KeyError, ValueError, TypeError):
            continue
        cur = out.get(r["gene"])
        if cur is None or s > float(cur["covalent_score"]):
            out[r["gene"]] = r
    return out


def main() -> int:
    if not RANK_TSV.exists() or not COV_TSV.exists():
        LOG.error("Run pipeline (stages 1-7) first")
        return 1
    if not shutil.which("pymol"):
        LOG.error("pymol not in PATH; install pymol-open-source via the conda env")
        return 1

    sites_cfg = yaml.safe_load(SITES_YAML.read_text()) if SITES_YAML.exists() else {}
    rank_rows = list(csv.DictReader(RANK_TSV.open(), delimiter="\t"))
    rank_rows.sort(key=lambda r: float(r.get("composite_score") or 0), reverse=True)

    cov_rows = list(csv.DictReader(COV_TSV.open(), delimiter="\t"))
    best_anchor = _best_anchor_per_gene(cov_rows)

    FIGURES.mkdir(parents=True, exist_ok=True)
    rendered = 0
    for rank in rank_rows[:TOP_N]:
        gene = rank["gene"]
        anchor = best_anchor.get(gene)
        if anchor is None:
            LOG.info("no covalent anchor for %s; skip viz", gene)
            continue
        pdb = PREP_DIR / f"{gene}.clean.pdb"
        if not pdb.exists():
            LOG.info("no PDB for %s; skip viz", gene)
            continue
        pml_path = FIGURES / f"top_{gene}.pml"
        pse_path = FIGURES / f"top_{gene}.pse"
        png_path = FIGURES / f"top_{gene}.png"
        pml_path.write_text(_build_pml(gene, anchor, sites_cfg, pdb, pse_path, png_path))
        try:
            r = subprocess.run(
                ["pymol", "-cq", str(pml_path)],
                capture_output=True,
                text=True,
                timeout=240,
            )
            if r.returncode != 0:
                LOG.warning("pymol failed for %s: %s", gene, r.stderr[-300:])
                continue
        except (subprocess.TimeoutExpired, FileNotFoundError) as exc:
            LOG.warning("pymol error for %s: %s", gene, exc)
            continue
        rendered += 1
        LOG.info("rendered figures/top_%s.png (+ .pse, .pml)", gene)

    LOG.info("rendered %d / %d top hits", rendered, min(TOP_N, len(rank_rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
