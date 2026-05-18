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
AF_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
RETRIES = 3
BACKOFF = 2.0


def _download(url: str, dest: Path, retries: int = RETRIES) -> bool:
    if dest.exists() and dest.stat().st_size > 0:
        LOG.info("cached %s", dest.name)
        return True
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 200 and r.content:
                dest.write_bytes(r.content)
                LOG.info("fetched %s (%d bytes)", dest.name, len(r.content))
                return True
            if r.status_code == 404:
                return False  # don't retry 404s
            LOG.warning("HTTP %d for %s", r.status_code, url)
        except requests.RequestException as exc:
            LOG.warning("attempt %d failed for %s: %s", attempt, url, exc)
        time.sleep(BACKOFF * attempt)
    return False


def _fetch_af(uniprot: str, dest: Path) -> bool:
    """Query the AlphaFold DB API for the current model URL, then download.

    The hard-coded v4 URL pattern 404s for many entries; the API returns
    the active URL whatever the version suffix.
    """
    if dest.exists() and dest.stat().st_size > 0:
        LOG.info("cached %s", dest.name)
        return True
    try:
        r = requests.get(AF_API.format(uniprot=uniprot), timeout=30)
    except requests.RequestException as exc:
        LOG.info("AF API unreachable for %s: %s", uniprot, exc)
        return False
    if r.status_code != 200:
        LOG.info("no AF model for %s (API HTTP %d)", uniprot, r.status_code)
        return False
    try:
        entries = r.json()
        pdb_url = entries[0]["pdbUrl"]
    except (ValueError, IndexError, KeyError, TypeError) as exc:
        LOG.info("AF API response unparseable for %s: %s", uniprot, exc)
        return False
    return _download(pdb_url, dest, retries=2)


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
            if not _fetch_af(t.uniprot, dest):
                LOG.info("no AF model for %s (%s)", t.gene, t.uniprot)
    if failed:
        LOG.warning("missing experimental structures: %s", ", ".join(failed))
    LOG.info("done: %d targets processed", len(targets))
    return 0


if __name__ == "__main__":
    sys.exit(main())
