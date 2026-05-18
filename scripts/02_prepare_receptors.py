"""Stage 1b: clean and protonate each receptor at pH 7.4, then export PDBQT.

Pipeline per structure:
  1. Strip waters, alt-locs (keep highest occupancy), non-cofactor HETATMs.
     Keep SAM/SAH (writers), Fe/2OG (erasers), m6A nucleotides (readers) as references.
  2. pdb2pqr30 --ff=AMBER --with-ph=7.4 --titration-state-method=propka -> .pqr
  3. prepare_receptor4.py -> .pdbqt
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

from _common import AF_DIR, PDB_DIR, PREP_DIR, get_logger, load_targets

LOG = get_logger("prepare")
KEEP_HETATMS = {"SAM", "SAH", "FE", "FE2", "AKG", "2OG", "M6A", "6MA", "MG", "ZN"}


def _strip_pdb(src: Path, dst: Path) -> None:
    """Keep ATOM records and whitelisted HETATM cofactors."""
    seen_alts: dict[tuple[str, str, int, str], str] = {}
    lines_out: list[str] = []
    with src.open() as fh:
        for line in fh:
            tag = line[:6]
            if tag not in ("ATOM  ", "HETATM"):
                if tag in ("HEADER", "TITLE ", "REMARK", "SEQRES", "TER   ", "END   "):
                    lines_out.append(line)
                continue
            altloc = line[16]
            resn = line[17:20].strip()
            chain = line[21]
            resi = int(line[22:26])
            atom = line[12:16].strip()
            if tag == "HETATM" and resn not in KEEP_HETATMS and resn != "MSE":
                continue
            if altloc not in (" ", "A"):
                key = (chain, resn, resi, atom)
                if key in seen_alts:
                    continue
                seen_alts[key] = altloc
            line = line[:16] + " " + line[17:]
            lines_out.append(line)
    dst.write_text("".join(lines_out))


def _run(cmd: list[str]) -> bool:
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if r.returncode != 0:
            LOG.warning("cmd failed: %s\nstderr: %s", " ".join(cmd), r.stderr[-500:])
            return False
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
        LOG.error("cmd error: %s: %s", " ".join(cmd), exc)
        return False


def _protonate(stripped: Path, pqr: Path) -> bool:
    return _run(
        [
            "pdb2pqr30",
            "--ff=AMBER",
            "--with-ph=7.4",
            "--titration-state-method=propka",
            "--keep-chain",
            str(stripped),
            str(pqr),
        ]
    )


def _to_pdbqt(pqr: Path, pdbqt: Path) -> bool:
    if shutil.which("prepare_receptor"):
        return _run(["prepare_receptor", "-r", str(pqr), "-o", str(pdbqt), "-A", "checkhydrogens"])
    return _run(["obabel", str(pqr), "-O", str(pdbqt), "-xr"])


def _pick_input(gene: str, uniprot: str, primary_pdb: str) -> Path | None:
    if primary_pdb and primary_pdb != "-":
        cand = PDB_DIR / f"{primary_pdb}.pdb"
        if cand.exists():
            return cand
    af = AF_DIR / f"AF-{uniprot}.pdb"
    return af if af.exists() else None


def main() -> int:
    PREP_DIR.mkdir(parents=True, exist_ok=True)
    targets = load_targets()
    ok = 0
    for t in targets:
        src = _pick_input(t.gene, t.uniprot, t.pdb_primary)
        if src is None:
            LOG.warning("no input for %s (%s)", t.gene, t.uniprot)
            continue
        stripped = PREP_DIR / f"{t.gene}.clean.pdb"
        pqr = PREP_DIR / f"{t.gene}.pqr"
        pdbqt = PREP_DIR / f"{t.gene}.pdbqt"
        _strip_pdb(src, stripped)
        if not _protonate(stripped, pqr):
            LOG.warning("protonation failed for %s; falling back to obabel -p", t.gene)
            _run(["obabel", str(stripped), "-O", str(pqr), "-p", "7.4"])
        if not _to_pdbqt(pqr if pqr.exists() else stripped, pdbqt):
            LOG.error("PDBQT conversion failed for %s", t.gene)
            continue
        ok += 1
    LOG.info("prepared %d / %d", ok, len(targets))
    return 0


if __name__ == "__main__":
    sys.exit(main())
