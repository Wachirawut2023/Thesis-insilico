"""Stage 4: non-covalent docking with AutoDock Vina across all docking boxes.

Reads results/docking_boxes.tsv (from stage 3) and runs Vina for each
gene x box combination. As(III) is treated approximately (As atom remapped via
ligand prep script — see ligands/prepare_as3.sh). ΔG values here are RELATIVE
ranking only; the primary readout for inhibition inference is covalent docking
in stage 5.
"""
from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path

from _common import (
    CHECKPOINT_DIR,
    LIGAND_DIR,
    PREP_DIR,
    RESULTS,
    append_tsv_rows,
    get_logger,
    load_checkpoint,
    mark_done,
)

LOG = get_logger("dock_nc")
BOXES_TSV = RESULTS / "docking_boxes.tsv"
OUT_TSV = RESULTS / "docking_noncovalent.tsv"
POSES_DIR = RESULTS / "docking_poses" / "noncovalent"
EXHAUSTIVENESS = 32
N_POSES = 20


def _have_vina() -> bool:
    return shutil.which("vina") is not None


def _run_vina(receptor: Path, ligand: Path, box: dict, out: Path) -> float | None:
    cmd = [
        "vina",
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", box["cx"], "--center_y", box["cy"], "--center_z", box["cz"],
        "--size_x", box["sx"], "--size_y", box["sy"], "--size_z", box["sz"],
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(N_POSES),
        "--out", str(out),
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    except subprocess.TimeoutExpired:
        LOG.warning("vina timeout for %s", out.name)
        return None
    if r.returncode != 0:
        LOG.warning("vina failed for %s: %s", out.name, r.stderr[-200:])
        return None
    # Parse best score from output log
    for line in r.stdout.splitlines():
        s = line.strip().split()
        if len(s) >= 4 and s[0] == "1":
            try:
                return float(s[1])
            except ValueError:
                pass
    return None


OUT_FIELDS = ["gene", "box_id", "anchor_cys", "vina_dG_kcal_mol", "pose_file"]
STAGE = "05_dock_nc"


def _box_key(gene: str, box_id: str) -> str:
    return f"{gene}::{box_id}"


def main() -> int:
    if not _have_vina():
        LOG.error("vina not in PATH; install via conda (vina>=1.2)")
        return 1
    ligand = LIGAND_DIR / "as3.pdbqt"
    if not ligand.exists():
        LOG.error("ligand missing: %s — run ligands/prepare_as3.sh first", ligand)
        return 1
    if not BOXES_TSV.exists():
        LOG.error("missing %s — run 04_pocket_detection.py first", BOXES_TSV)
        return 1
    POSES_DIR.mkdir(parents=True, exist_ok=True)
    done = load_checkpoint(STAGE)
    with BOXES_TSV.open() as fh:
        boxes = list(csv.DictReader(fh, delimiter="\t"))
    n_new = 0
    for box in boxes:
        key = _box_key(box["gene"], box["box_id"])
        if key in done:
            continue
        recep = PREP_DIR / f"{box['gene']}.pdbqt"
        if not recep.exists():
            LOG.warning("no receptor for %s; skipping box %s (not checkpointed, retry after prep)", box["gene"], box["box_id"])
            continue
        pose = POSES_DIR / f"{box['gene']}_{box['box_id']}.pdbqt"
        score = _run_vina(recep, ligand, box, pose)
        row = {
            "gene": box["gene"],
            "box_id": box["box_id"],
            "anchor_cys": box["anchor_cys"],
            "vina_dG_kcal_mol": f"{score:.3f}" if score is not None else "NA",
            "pose_file": str(pose.relative_to(RESULTS)) if pose.exists() else "",
        }
        append_tsv_rows([row], OUT_TSV, OUT_FIELDS)
        mark_done(STAGE, key)
        n_new += 1
        LOG.info("[%s/%s] vina_dG=%s (%d/%d done)", box["gene"], box["box_id"], row["vina_dG_kcal_mol"], len(done) + n_new, len(boxes))

    LOG.info(
        "processed %d new box(es), %d already checkpointed (of %d total) -> %s",
        n_new, len(done), len(boxes), OUT_TSV,
    )
    LOG.info("checkpoints in %s — delete %s/%s.done to force a full recompute", CHECKPOINT_DIR, CHECKPOINT_DIR, STAGE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
