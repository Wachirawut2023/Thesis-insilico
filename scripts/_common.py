"""Shared utilities for the As(III) m6A in silico pipeline."""
from __future__ import annotations

import csv
import logging
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
