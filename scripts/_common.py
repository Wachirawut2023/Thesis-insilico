"""Shared utilities for the As(III) m6A in silico pipeline."""
from __future__ import annotations

import csv
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
PDB_DIR = DATA / "pdb"
AF_DIR = DATA / "af"
PREP_DIR = DATA / "prepared"
LIGAND_DIR = ROOT / "ligands"
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
TARGETS_TSV = DATA / "targets.tsv"
CHECKPOINT_DIR = RESULTS / ".checkpoints"


@dataclass(frozen=True)
class Target:
    family: str
    subfamily: str
    gene: str
    uniprot: str
    pdb_primary: str
    pdb_alt: tuple[str, ...]
    role: str
    notes: str

    @property
    def has_experimental(self) -> bool:
        return bool(self.pdb_primary) and self.pdb_primary != "-"

    @property
    def all_pdbs(self) -> tuple[str, ...]:
        out: list[str] = []
        if self.has_experimental:
            out.append(self.pdb_primary)
        out.extend(p for p in self.pdb_alt if p and p != "-")
        return tuple(out)


def load_targets(path: Path = TARGETS_TSV) -> list[Target]:
    targets: list[Target] = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            alts = tuple(
                p.strip() for p in (row.get("pdb_alt") or "").split(",") if p.strip() and p.strip() != "-"
            )
            targets.append(
                Target(
                    family=row["family"],
                    subfamily=row["subfamily"],
                    gene=row["gene"],
                    uniprot=row["uniprot"],
                    pdb_primary=row["pdb_primary"].strip(),
                    pdb_alt=alts,
                    role=row["role"],
                    notes=row.get("notes", ""),
                )
            )
    # Honour SMOKE_GENES env var: comma-separated gene names to keep.
    smoke = os.environ.get("SMOKE_GENES", "").strip()
    if smoke:
        keep = {g.strip() for g in smoke.split(",") if g.strip()}
        targets = [t for t in targets if t.gene in keep]
    return targets


def get_logger(name: str) -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(name)s] %(levelname)s %(message)s",
    )
    return logging.getLogger(name)


def write_tsv(rows: Iterable[dict], path: Path, fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ── Checkpoint/resume helpers ──────────────────────────────────────────────
#
# Long stages (pocket detection, docking) process one target at a time and
# can run for hours; a Colab crash/disconnect mid-stage should not lose
# already-computed targets. Each stage keeps a `results/.checkpoints/<stage>.done`
# file listing genes it has fully finished (one per line, appended as each
# gene completes) and appends that gene's output rows to the TSV immediately
# rather than buffering everything in memory until the end. On restart, a
# stage skips any gene already listed in its checkpoint file.
#
# `make clean` removes results/.checkpoints/ along with the TSVs it guards,
# so a fresh run never skips work it hasn't actually produced.


def load_checkpoint(stage: str) -> set[str]:
    """Return the set of keys (usually gene symbols) already completed for `stage`."""
    p = CHECKPOINT_DIR / f"{stage}.done"
    if not p.exists():
        return set()
    return {line.strip() for line in p.read_text().splitlines() if line.strip()}


def mark_done(stage: str, key: str) -> None:
    """Record that `key` has fully completed for `stage` (call after its rows are on disk)."""
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    p = CHECKPOINT_DIR / f"{stage}.done"
    with p.open("a") as fh:
        fh.write(key + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def append_tsv_rows(rows: Iterable[dict], path: Path, fieldnames: list[str]) -> None:
    """Append rows to a TSV, creating it (with header) if it doesn't exist yet."""
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    is_new = not path.exists()
    with path.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        if is_new:
            w.writeheader()
        for r in rows:
            w.writerow(r)
        fh.flush()
        os.fsync(fh.fileno())
