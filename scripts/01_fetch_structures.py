"""Stage 1a: download PDB structures and AlphaFold models for the target panel.

Sources:
  - RCSB: https://files.rcsb.org/download/<PDBID>.pdb
  - AlphaFold DB: https://alphafold.ebi.ac.uk/files/AF-<UniProt>-F1-model_v4.pdb
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import requests

from _common import AF_DIR, PDB_DIR, get_logger, load_targets

LOG = get_logger("fetch")
RCSB_URL = "https://files.rcsb.org/download/{pdb}.pdb"
AF_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"
RETRIES = 3
BACKOFF = 2.0


def _download(url: str, dest: Path) -> bool:
    if dest.exists() and dest.stat().st_size > 0:
        LOG.info("cached %s", dest.name)
        return True
    for attempt in range(1, RETRIES + 1):
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 200 and r.content:
                dest.write_bytes(r.content)
                LOG.info("fetched %s (%d bytes)", dest.name, len(r.content))
                return True
            LOG.warning("HTTP %d for %s", r.status_code, url)
        except requests.RequestException as exc:
            LOG.warning("attempt %d failed for %s: %s", attempt, url, exc)
        time.sleep(BACKOFF * attempt)
    return False


def main() -> int:
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    AF_DIR.mkdir(parents=True, exist_ok=True)
    targets = load_targets()
    failed: list[str] = []
    for t in targets:
        for pdb in t.all_pdbs:
            dest = PDB_DIR / f"{pdb}.pdb"
            if not _download(RCSB_URL.format(pdb=pdb), dest):
                failed.append(f"{t.gene}:{pdb}")
        if t.uniprot and t.uniprot != "-":
            dest = AF_DIR / f"AF-{t.uniprot}.pdb"
            if not _download(AF_URL.format(uniprot=t.uniprot), dest):
                LOG.info("no AF model for %s (%s)", t.gene, t.uniprot)
    if failed:
        LOG.warning("missing experimental structures: %s", ", ".join(failed))
    LOG.info("done: %d targets processed", len(targets))
    return 0


if __name__ == "__main__":
    sys.exit(main())
