"""Stage 3: pocket detection with fpocket; intersect pockets with vicinal-Cys clusters.

Defines docking boxes for Stage 4/5 by finding pockets whose centroid is within
8 A of a vicinal Cys cluster centre. Writes results/pocket_table.tsv and
results/docking_boxes.tsv (one row per docking site to evaluate).
"""
from __future__ import annotations

import csv
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

from _common import PREP_DIR, RESULTS, get_logger, load_targets

LOG = get_logger("pockets")
POCKET_TSV = RESULTS / "pocket_table.tsv"
BOXES_TSV = RESULTS / "docking_boxes.tsv"
CYS_TSV = RESULTS / "cys_table.tsv"
CLUSTER_RADIUS = 8.0
BOX_PAD = 8.0


def _run_fpocket(pdb: Path) -> Path | None:
    if not shutil.which("fpocket"):
        LOG.error("fpocket not in PATH")
        return None
    out_dir = pdb.parent / f"{pdb.stem}_out"
    try:
        subprocess.run(["fpocket", "-f", str(pdb)], check=True, capture_output=True, timeout=600)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        LOG.warning("fpocket failed on %s: %s", pdb.name, exc)
        return None
    return out_dir if out_dir.exists() else None


POCKET_HDR = re.compile(r"Pocket\s+(\d+)\s*:")


def _parse_fpocket_info(info_file: Path) -> list[dict]:
    if not info_file.exists():
        return []
    pockets: list[dict] = []
    cur: dict | None = None
    for line in info_file.read_text().splitlines():
        m = POCKET_HDR.match(line.strip())
        if m:
            if cur:
                pockets.append(cur)
            cur = {"pocket": int(m.group(1))}
            continue
        if cur is None:
            continue
        if "Druggability Score" in line:
            cur["druggability"] = float(line.split(":")[-1].strip())
        elif "Score :" in line and "Druggability" not in line:
            try:
                cur["score"] = float(line.split(":")[-1].strip())
            except ValueError:
                pass
        elif "Number of Alpha Spheres" in line:
            cur["n_alpha"] = int(line.split(":")[-1].strip())
    if cur:
        pockets.append(cur)
    return pockets


def _pocket_centroid(pocket_pdb: Path) -> tuple[float, float, float] | None:
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    if not pocket_pdb.exists():
        return None
    for line in pocket_pdb.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
    if not xs:
        return None
    n = len(xs)
    return (sum(xs) / n, sum(ys) / n, sum(zs) / n)


def _read_cys_clusters(gene: str) -> list[tuple[float, float, float, str]]:
    """Build cluster centroids from reactive Cys with >=1 vicinal partner."""
    if not CYS_TSV.exists():
        return []
    pdb = PREP_DIR / f"{gene}.clean.pdb"
    if not pdb.exists():
        return []
    coords: dict[tuple[str, int], tuple[float, float, float]] = {}
    for line in pdb.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[17:20].strip() != "CYS" or line[12:16].strip() != "SG":
            continue
        coords[(line[21], int(line[22:26]))] = (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        )
    clusters: list[tuple[float, float, float, str]] = []
    with CYS_TSV.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["gene"] != gene or row["reactive"] != "Y":
                continue
            if int(row["n_vicinal_partners"]) < 1:
                continue
            key = (row["chain"], int(row["resi"]))
            if key not in coords:
                continue
            x, y, z = coords[key]
            tag = f"{row['chain']}{row['resi']}"
            clusters.append((x, y, z, tag))
    return clusters


def _box_around(centre: tuple[float, float, float], pad: float = BOX_PAD) -> tuple[float, float, float]:
    return (pad * 2, pad * 2, pad * 2)


def main() -> int:
    pocket_rows: list[dict] = []
    box_rows: list[dict] = []
    for t in load_targets():
        pdb = PREP_DIR / f"{t.gene}.clean.pdb"
        if not pdb.exists():
            continue
        out_dir = _run_fpocket(pdb)
        pockets: list[dict] = []
        if out_dir is not None:
            info = out_dir / f"{pdb.stem}_info.txt"
            pockets = _parse_fpocket_info(info)
            for p in pockets:
                pp = out_dir / "pockets" / f"pocket{p['pocket']}_atm.pdb"
                p["centroid"] = _pocket_centroid(pp)
        clusters = _read_cys_clusters(t.gene)
        # log pockets
        for p in pockets:
            cx, cy, cz = p.get("centroid") or (None, None, None)
            pocket_rows.append(
                {
                    "gene": t.gene,
                    "pocket_id": p["pocket"],
                    "score": p.get("score", ""),
                    "druggability": p.get("druggability", ""),
                    "n_alpha": p.get("n_alpha", ""),
                    "cx": f"{cx:.3f}" if cx is not None else "",
                    "cy": f"{cy:.3f}" if cy is not None else "",
                    "cz": f"{cz:.3f}" if cz is not None else "",
                }
            )
        # boxes: pocket+cluster intersections
        seen_box_ids: set[str] = set()
        for p in pockets:
            cen = p.get("centroid")
            if cen is None:
                continue
            for cx, cy, cz, tag in clusters:
                d = math.sqrt((cen[0] - cx) ** 2 + (cen[1] - cy) ** 2 + (cen[2] - cz) ** 2)
                if d <= CLUSTER_RADIUS:
                    box_id = f"P{p['pocket']}_C{tag}"
                    if box_id in seen_box_ids:
                        continue
                    seen_box_ids.add(box_id)
                    sx, sy, sz = _box_around(cen)
                    box_rows.append(
                        {
                            "gene": t.gene,
                            "box_id": box_id,
                            "origin": "pocket+cys",
                            "cx": f"{cen[0]:.3f}",
                            "cy": f"{cen[1]:.3f}",
                            "cz": f"{cen[2]:.3f}",
                            "sx": f"{sx:.1f}",
                            "sy": f"{sy:.1f}",
                            "sz": f"{sz:.1f}",
                            "anchor_cys": tag,
                            "pocket_id": p["pocket"],
                        }
                    )
        # fallback: cys-only boxes if no pocket overlap found
        if not any(r["gene"] == t.gene for r in box_rows):
            for cx, cy, cz, tag in clusters:
                sx, sy, sz = _box_around((cx, cy, cz))
                box_rows.append(
                    {
                        "gene": t.gene,
                        "box_id": f"C{tag}",
                        "origin": "cys_only",
                        "cx": f"{cx:.3f}",
                        "cy": f"{cy:.3f}",
                        "cz": f"{cz:.3f}",
                        "sx": f"{sx:.1f}",
                        "sy": f"{sy:.1f}",
                        "sz": f"{sz:.1f}",
                        "anchor_cys": tag,
                        "pocket_id": "",
                    }
                )

    from _common import write_tsv

    write_tsv(
        pocket_rows,
        POCKET_TSV,
        ["gene", "pocket_id", "score", "druggability", "n_alpha", "cx", "cy", "cz"],
    )
    write_tsv(
        box_rows,
        BOXES_TSV,
        ["gene", "box_id", "origin", "cx", "cy", "cz", "sx", "sy", "sz", "anchor_cys", "pocket_id"],
    )
    LOG.info("wrote %d pockets, %d docking boxes", len(pocket_rows), len(box_rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
