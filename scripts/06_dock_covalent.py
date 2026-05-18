"""Stage 5: geometric covalent docking — the primary readout for inhibition inference.

Rather than running a full covalent docking suite, we exploit the fact that
As(III)–Cys coordination is geometrically highly constrained:
  - As–S bond length  ~ 2.25 A (literature)
  - S–As–S angle      ~ 94 deg (trigonal pyramidal)
  - preferred coordination: bidentate As(OH)(SR)(SR'), tridentate As(SR)3

For each reactive Cys (from cys_table.tsv), we:
  1. Build the As anchor at Sγ at the ideal As–S distance, placing the As
     position by reflecting Sγ through Cβ along the Cβ→Sγ vector.
  2. For each vicinal reactive Cys partner, check feasibility of a 2nd As–S
     bond (As–Sγ_partner distance allowed range 3.0–4.5 A after small As
     translation, S–As–S angle within 75–115 deg).
  3. Triadic coordination: check 3-way feasibility (As within ~2.4–3.5 A of
     all 3 Sγ centroids).
  4. Score each anchor:
        S_anchor =  w_tri * tridentate_feasible
                  + w_di  * bidentate_feasible
                  - w_clash * heavy_atom_clashes
                  - w_func * func_proximity_penalty
                  + w_pka  * (1 if reactive thiolate else 0)
  5. Per-protein covalent score = max anchor score + sum of anchor scores * 0.05

Writes results/docking_covalent.tsv (one row per Cys anchor evaluated)
and adds covalent_dG_proxy to the protein summary in stage 7.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

from _common import PREP_DIR, RESULTS, get_logger, load_targets, write_tsv

LOG = get_logger("dock_cov")
CYS_TSV = RESULTS / "cys_table.tsv"
OUT_TSV = RESULTS / "docking_covalent.tsv"

AS_S_BOND = 2.25
AS_S_TOL = 0.30
S_AS_S_IDEAL = 94.0
S_AS_S_TOL = 20.0
TRI_MAX = 3.5
DI_RANGE = (3.0, 4.6)
CLASH_THR = 2.0


def _parse_atoms(pdb: Path) -> list[dict]:
    atoms: list[dict] = []
    for line in pdb.read_text().splitlines():
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


def _vec(a, b):
    return (b["x"] - a["x"], b["y"] - a["y"], b["z"] - a["z"])


def _norm(v):
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


def _unit(v):
    n = _norm(v)
    return (v[0] / n, v[1] / n, v[2] / n) if n > 0 else (0.0, 0.0, 0.0)


def _add(a, b, scale=1.0):
    return {"x": a["x"] + b[0] * scale, "y": a["y"] + b[1] * scale, "z": a["z"] + b[2] * scale}


def _dist(a, b) -> float:
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def _angle(a, b, c) -> float:
    """Angle at b in degrees."""
    ba = (a["x"] - b["x"], a["y"] - b["y"], a["z"] - b["z"])
    bc = (c["x"] - b["x"], c["y"] - b["y"], c["z"] - b["z"])
    nba = _norm(ba)
    nbc = _norm(bc)
    if nba == 0 or nbc == 0:
        return 0.0
    cos = (ba[0] * bc[0] + ba[1] * bc[1] + ba[2] * bc[2]) / (nba * nbc)
    cos = max(-1.0, min(1.0, cos))
    return math.degrees(math.acos(cos))


def _build_as_anchor(sg: dict, cb: dict | None) -> dict:
    """Place As at AS_S_BOND from Sγ, opposite Cβ→Sγ direction."""
    if cb is None:
        return {"x": sg["x"] + AS_S_BOND, "y": sg["y"], "z": sg["z"]}
    direction = _unit(_vec(cb, sg))
    return _add(sg, direction, AS_S_BOND)


def _clash_count(pos: dict, all_atoms: list[dict], skip_resi: set[tuple[str, int]]) -> int:
    n = 0
    for a in all_atoms:
        if (a["chain"], a["resi"]) in skip_resi:
            continue
        if a["atom"] == "H":
            continue
        if _dist(pos, a) < CLASH_THR:
            n += 1
    return n


def analyse(gene: str, cys_rows: list[dict]) -> list[dict]:
    pdb = PREP_DIR / f"{gene}.clean.pdb"
    if not pdb.exists():
        return []
    atoms = _parse_atoms(pdb)
    sg = {(a["chain"], a["resi"]): a for a in atoms if a["resn"] == "CYS" and a["atom"] == "SG"}
    cb = {(a["chain"], a["resi"]): a for a in atoms if a["resn"] == "CYS" and a["atom"] == "CB"}

    out: list[dict] = []
    for row in cys_rows:
        if row["reactive"] != "Y":
            continue
        key = (row["chain"], int(row["resi"]))
        if key not in sg:
            continue
        anchor_sg = sg[key]
        anchor_cb = cb.get(key)
        as_pos = _build_as_anchor(anchor_sg, anchor_cb)

        # Parse vicinal partners string: "A33:4.21;A36:6.10"
        partners: list[tuple[tuple[str, int], float]] = []
        if row["vicinal_partners"]:
            for tok in row["vicinal_partners"].split(";"):
                if not tok:
                    continue
                tag, dstr = tok.split(":")
                ch, resi = tag[0], int(tag[1:])
                partners.append(((ch, resi), float(dstr)))

        # Bidentate / tridentate feasibility
        di_feas = 0
        tri_feas = 0
        partner_tags: list[str] = []
        for pkey, _d in partners:
            psg = sg.get(pkey)
            if psg is None:
                continue
            d_as_p = _dist(as_pos, psg)
            if DI_RANGE[0] <= d_as_p <= DI_RANGE[1]:
                ang = _angle(anchor_sg, as_pos, psg)
                if abs(ang - S_AS_S_IDEAL) <= S_AS_S_TOL:
                    di_feas += 1
                    partner_tags.append(f"{pkey[0]}{pkey[1]}")
        # tridentate: any pair of partners that are also close to each other
        for i in range(len(partners)):
            for j in range(i + 1, len(partners)):
                p1 = sg.get(partners[i][0])
                p2 = sg.get(partners[j][0])
                if p1 is None or p2 is None:
                    continue
                d1 = _dist(as_pos, p1)
                d2 = _dist(as_pos, p2)
                if d1 <= TRI_MAX and d2 <= TRI_MAX:
                    tri_feas += 1

        clashes = _clash_count(as_pos, atoms, {key})
        fp_str = row.get("functional_proximity_A", "NA")
        try:
            fp = float(fp_str) if fp_str != "NA" else None
        except ValueError:
            fp = None
        func_penalty = max(0.0, (fp - 6.0) / 6.0) if fp is not None else 0.5
        pka_str = row.get("pka", "NA")
        try:
            pka_val = float(pka_str) if pka_str != "NA" else None
        except ValueError:
            pka_val = None
        thiolate_bonus = 1.0 if (pka_val is not None and pka_val <= 7.5) else 0.0

        score = (
            2.0 * tri_feas
            + 1.0 * di_feas
            - 0.5 * clashes
            - 1.0 * func_penalty
            + 0.5 * thiolate_bonus
        )

        out.append(
            {
                "gene": gene,
                "anchor_chain": row["chain"],
                "anchor_resi": row["resi"],
                "n_vicinal_partners": row["n_vicinal_partners"],
                "bidentate_feasible": di_feas,
                "tridentate_feasible": tri_feas,
                "partner_tags": ";".join(partner_tags),
                "clashes": clashes,
                "functional_proximity_A": fp_str,
                "thiolate_bonus": thiolate_bonus,
                "covalent_score": round(score, 3),
            }
        )
    return out


def main() -> int:
    if not CYS_TSV.exists():
        LOG.error("missing %s — run 03_cys_landscape.py first", CYS_TSV)
        return 1
    by_gene: dict[str, list[dict]] = {}
    with CYS_TSV.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_gene.setdefault(row["gene"], []).append(row)
    rows: list[dict] = []
    for t in load_targets():
        rows.extend(analyse(t.gene, by_gene.get(t.gene, [])))
    write_tsv(
        rows,
        OUT_TSV,
        [
            "gene",
            "anchor_chain",
            "anchor_resi",
            "n_vicinal_partners",
            "bidentate_feasible",
            "tridentate_feasible",
            "partner_tags",
            "clashes",
            "functional_proximity_A",
            "thiolate_bonus",
            "covalent_score",
        ],
    )
    LOG.info("wrote %d covalent anchor rows -> %s", len(rows), OUT_TSV)
    return 0


if __name__ == "__main__":
    sys.exit(main())
