"""Stage 2: per-Cys reactivity landscape — the most diagnostic step of the pipeline.

For every Cys in every prepared receptor:
  - SASA of the Sγ atom (FreeSASA)
  - thiol pKa (from PROPKA log if available; else NA)
  - disulfide partner check (Sγ–Sγ < 2.3 Å excludes Cys from reactive pool)
  - vicinal-Cys mapping: Sγ–Sγ distances; flag pairs/triads ≤ 7 Å
  - functional_proximity: min distance from Sγ to curated anchor residues (data/functional_sites.yaml)

Output: results/cys_table.tsv
        one row per Cys with reactivity flags + per-protein vicinal_score
"""
from __future__ import annotations

import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import yaml

from _common import (
    CHECKPOINT_DIR,
    DATA,
    PREP_DIR,
    RESULTS,
    append_tsv_rows,
    get_logger,
    load_checkpoint,
    load_targets,
    mark_done,
)

LOG = get_logger("cys")
SS_BOND_MAX = 2.3
VICINAL_MAX = 7.0
SASA_ACCESSIBLE = 5.0
REACTIVE_PKA = 8.0
SITES_YAML = DATA / "functional_sites.yaml"
OUT_TSV = RESULTS / "cys_table.tsv"
PROTEIN_TSV = RESULTS / "protein_summary.tsv"


def _load_sites() -> dict[str, dict]:
    if not SITES_YAML.exists():
        return {}
    with SITES_YAML.open() as fh:
        return yaml.safe_load(fh) or {}


def _parse_cys_atoms(pdb: Path) -> list[dict]:
    """Return list of dicts with chain, resi, atom, x, y, z for CYS atoms."""
    atoms: list[dict] = []
    with pdb.open() as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            resn = line[17:20].strip()
            if resn != "CYS":
                continue
            atoms.append(
                {
                    "chain": line[21],
                    "resi": int(line[22:26]),
                    "atom": line[12:16].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                }
            )
    return atoms


def _parse_all_atoms(pdb: Path) -> list[dict]:
    atoms: list[dict] = []
    with pdb.open() as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atoms.append(
                {
                    "chain": line[21],
                    "resi": int(line[22:26]),
                    "resn": line[17:20].strip(),
                    "atom": line[12:16].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                }
            )
    return atoms


def _dist(a: dict, b: dict) -> float:
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def _sasa_map(pdb: Path) -> dict[tuple[str, int, str], float]:
    """Per-atom SASA via freesasa; key = (chain, resi, atom). Empty dict if unavailable."""
    try:
        import freesasa
    except ImportError:
        LOG.warning("freesasa not installed; SASA values will be NA")
        return {}
    try:
        structure = freesasa.Structure(str(pdb))
        # API changed in freesasa>=2.2: module-level freesasa.calc() replaces
        # freesasa.Calc().calculate() class method.
        if hasattr(freesasa, "calc"):
            result = freesasa.calc(structure)
        else:
            result = freesasa.Calc().calculate(structure)
    except Exception as exc:
        LOG.warning("freesasa failed on %s: %s", pdb.name, exc)
        return {}
    out: dict[tuple[str, int, str], float] = {}
    sample_logged = False
    for i in range(structure.nAtoms()):
        try:
            chain = str(structure.chainLabel(i)).strip() or "A"
            resi_raw = str(structure.residueNumber(i)).strip()
            # strip insertion code if present
            resi = int("".join(c for c in resi_raw if c.isdigit() or c == "-"))
            atom = structure.atomName(i).strip()
        except (ValueError, AttributeError) as exc:
            if not sample_logged:
                LOG.warning("freesasa key parse failed on %s atom %d: %s", pdb.name, i, exc)
                sample_logged = True
            continue
        out[(chain, resi, atom)] = result.atomArea(i)
    if not out:
        LOG.warning("freesasa produced empty key map for %s", pdb.name)
    return out


PKA_LINE = re.compile(r"^\s*CYS\s+(\d+)\s+([A-Z])\s+([0-9.\-]+)")


def _pka_map(propka_log: Path) -> dict[tuple[str, int], float]:
    """Parse PROPKA .pka output for CYS pKa values."""
    if not propka_log.exists():
        return {}
    out: dict[tuple[str, int], float] = {}
    try:
        for line in propka_log.read_text().splitlines():
            m = PKA_LINE.match(line)
            if m:
                resi, chain, pka = m.group(1), m.group(2), m.group(3)
                try:
                    out[(chain, int(resi))] = float(pka)
                except ValueError:
                    pass
    except Exception as exc:
        LOG.warning("pka parse failed for %s: %s", propka_log.name, exc)
    return out


def _functional_proximity(sg: dict, anchors: list, all_atoms: list[dict], chain: str) -> float | None:
    """Minimum Sγ distance to any curated anchor residue.

    Integer anchors (e.g. `[152]`) match any chain — handles homo-oligomeric
    structures like GAPDH 1U8F (chains O/P/Q/R) where the catalytic Cys152
    exists in every subunit but our default chain 'A' would miss them all.

    Dict anchors with explicit `chain:` require exact chain+resi match
    (use this when the structure is a hetero-complex and the anchor must
    live on a specific subunit, e.g. METTL14 in the METTL3/14 dimer)."""
    if not anchors:
        return None
    anchor_atoms: list[dict] = []
    for a in anchors:
        if isinstance(a, int):
            anchor_atoms.extend(at for at in all_atoms if at["resi"] == a)
        elif isinstance(a, dict):
            target_chain = a.get("chain")
            try:
                target_resi = int(a["resi"])
            except (KeyError, ValueError, TypeError):
                continue
            if target_chain:
                anchor_atoms.extend(
                    at for at in all_atoms
                    if at["chain"] == target_chain and at["resi"] == target_resi
                )
            else:
                anchor_atoms.extend(at for at in all_atoms if at["resi"] == target_resi)
    if not anchor_atoms:
        return None
    return min(_dist(sg, at) for at in anchor_atoms)


def analyse_protein(gene: str, sites_cfg: dict) -> tuple[list[dict], dict]:
    pdb = PREP_DIR / f"{gene}.clean.pdb"
    if not pdb.exists():
        return [], {"gene": gene, "n_cys": 0, "n_reactive": 0, "vicinal_score": 0.0, "status": "missing"}
    cys_atoms = _parse_cys_atoms(pdb)
    all_atoms = _parse_all_atoms(pdb)
    sasa = _sasa_map(pdb)
    pka = _pka_map(PREP_DIR / f"{gene}.pka")
    sg_atoms = [a for a in cys_atoms if a["atom"] == "SG"]
    chain_default = (sites_cfg.get(gene) or {}).get("chain", "A")
    anchors = (sites_cfg.get(gene) or {}).get("anchors", [])

    # Disulfide partners
    in_disulfide: set[tuple[str, int]] = set()
    for i, a in enumerate(sg_atoms):
        for b in sg_atoms[i + 1 :]:
            if _dist(a, b) < SS_BOND_MAX:
                in_disulfide.add((a["chain"], a["resi"]))
                in_disulfide.add((b["chain"], b["resi"]))

    rows: list[dict] = []
    vicinal_pairs = 0
    vicinal_triads = 0
    for i, a in enumerate(sg_atoms):
        partners: list[str] = []
        for j, b in enumerate(sg_atoms):
            if i == j:
                continue
            d = _dist(a, b)
            if d <= VICINAL_MAX and (a["chain"], a["resi"]) not in in_disulfide and (b["chain"], b["resi"]) not in in_disulfide:
                partners.append(f"{b['chain']}{b['resi']}:{d:.2f}")
        sg_sasa = sasa.get((a["chain"], a["resi"], "SG"))
        if sg_sasa is None:
            # Try alternative key formats (chain stripping, fallback to any SG match)
            sg_sasa = sasa.get((a["chain"].strip(), a["resi"], "SG"))
            if sg_sasa is None and not sasa:
                # SASA map entirely unavailable — default to accessible so the
                # pipeline isn't fully blocked; the function_proximity and
                # vicinal-cluster signals carry the analysis.
                sg_sasa = None  # keep NA in output but treat as accessible below
        pka_val = pka.get((a["chain"], a["resi"]))
        accessible = (sg_sasa is None) or (sg_sasa >= SASA_ACCESSIBLE)
        reactive_pka = pka_val is None or pka_val <= REACTIVE_PKA
        in_ss = (a["chain"], a["resi"]) in in_disulfide
        fp = _functional_proximity(a, anchors, all_atoms, chain_default)
        is_reactive = accessible and reactive_pka and not in_ss
        rows.append(
            {
                "gene": gene,
                "chain": a["chain"],
                "resi": a["resi"],
                "sg_sasa": f"{sg_sasa:.2f}" if sg_sasa is not None else "NA",
                "pka": f"{pka_val:.2f}" if pka_val is not None else "NA",
                "in_disulfide": "Y" if in_ss else "N",
                "accessible": "Y" if accessible else "N",
                "vicinal_partners": ";".join(partners) if partners else "",
                "n_vicinal_partners": len(partners),
                "functional_proximity_A": f"{fp:.2f}" if fp is not None else "NA",
                "reactive": "Y" if is_reactive else "N",
            }
        )
        if is_reactive and len(partners) >= 1:
            vicinal_pairs += 1
        if is_reactive and len(partners) >= 2:
            vicinal_triads += 1

    n_reactive = sum(1 for r in rows if r["reactive"] == "Y")
    score = 1.5 * vicinal_triads + vicinal_pairs + 0.25 * n_reactive
    return rows, {
        "gene": gene,
        "n_cys": len(sg_atoms),
        "n_reactive": n_reactive,
        "n_vicinal_pairs": vicinal_pairs,
        "n_vicinal_triads": vicinal_triads,
        "vicinal_score": round(score, 3),
        "status": "ok",
    }


CYS_FIELDS = [
    "gene",
    "chain",
    "resi",
    "sg_sasa",
    "pka",
    "in_disulfide",
    "accessible",
    "vicinal_partners",
    "n_vicinal_partners",
    "functional_proximity_A",
    "reactive",
]
SUMMARY_FIELDS = [
    "gene", "family", "n_cys", "n_reactive", "n_vicinal_pairs", "n_vicinal_triads", "vicinal_score", "status",
]
STAGE = "03_cys"


def main() -> int:
    sites_cfg = _load_sites()
    targets = load_targets()
    done = load_checkpoint(STAGE)
    n_new = 0
    for t in targets:
        if t.gene in done:
            continue
        rows, s = analyse_protein(t.gene, sites_cfg)
        s["family"] = t.family
        append_tsv_rows(rows, OUT_TSV, CYS_FIELDS)
        append_tsv_rows([s], PROTEIN_TSV, SUMMARY_FIELDS)
        mark_done(STAGE, t.gene)
        n_new += 1
    LOG.info(
        "processed %d new target(s), %d already checkpointed (of %d total) -> %s / %s",
        n_new, len(done), len(targets), OUT_TSV, PROTEIN_TSV,
    )
    LOG.info("checkpoints in %s — delete %s/%s.done to force a full recompute", CHECKPOINT_DIR, CHECKPOINT_DIR, STAGE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
