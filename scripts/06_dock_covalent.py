"""Stage 5: geometric covalent docking — primary inhibition-inference readout.

Scores every reactive Cys for As(III) covalent binding feasibility in three modes:

  monodentate As-SR
      single Cys; always plausible if Cys is reactive. Scored via the
      functional-proximity bonus rather than as a standalone numeric feasibility.

  bidentate As(OH)(SR)(SR')
      vicinal Cys pair with Sγ-Sγ distance in [3.0, 4.4] Å. Geometric window
      derived from As-Sγ 2.25 Å and S-As-S ~94° (law of cosines: ideal Sγ-Sγ
      ≈ 3.3 Å, max feasible ≈ 4.4 Å before As-Sγ stretches beyond bond length).
      As is placed at the apex of the isoceles triangle bridging the two Sγ.
      A vicinal pair whose *static* Sγ-Sγ distance exceeds this window may still
      bridge on induced fit — As(III) has strong affinity for vicinal dithiols —
      so we additionally scan the Cys χ1 rotamer (its only side-chain rotatable)
      of both partners and report the minimum achievable, clash-free Sγ-Sγ
      distance (columns min_rotamer_sg_sg_A / bidentate_feasible_rotamer). This
      distinguishes "bidentate geometrically impossible" from "bidentate
      reachable by a modest side-chain rotation" — the latter is a much
      stronger, As-specific covalent-inhibition claim (relevant to the METTL3
      Cys375/Cys376 pair, ~6.75 Å apart in the crystal).

  tridentate As(SR)3
      Three Cys forming a triangle with all pairwise Sγ-Sγ ≤ 4.4 Å.

Anchor score combines feasibility flags, clash penalty, functional proximity
(Sγ to catalytic/binding-site residues from functional_sites.yaml), and
thiolate reactivity (pKa).
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
BIDENTATE_RANGE = (3.0, 4.4)  # Sγ-Sγ window for bidentate As(III)
TRIDENTATE_MAX = 4.4
CLASH_THR = 2.0
CLASH_TOLERANCE = 2  # accept this many heavy-atom clashes (typical for small ligands)

# Induced-fit / rotamer bidentate feasibility. Cys has a single side-chain
# rotatable (χ1 = N-Cα-Cβ-Sγ); scanning it tests whether a vicinal pair whose
# static Sγ-Sγ distance is outside BIDENTATE_RANGE can still form an As(III)
# bridge after a modest side-chain rotation.
CB_SG_BOND = 1.808            # Å, Cβ-Sγ bond length
CA_CB_SG_ANGLE = math.radians(114.0)  # Cα-Cβ-Sγ angle
ROTAMER_STEP = 10             # degrees, χ1 scan resolution


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


def _dist(a: dict, b: dict) -> float:
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def _vec(a: dict, b: dict) -> tuple[float, float, float]:
    return (b["x"] - a["x"], b["y"] - a["y"], b["z"] - a["z"])


def _unit(v: tuple[float, float, float]) -> tuple[float, float, float]:
    n = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    return (v[0] / n, v[1] / n, v[2] / n) if n > 0 else (0.0, 0.0, 0.0)


def _xyz(a: dict) -> tuple[float, float, float]:
    return (a["x"], a["y"], a["z"])


def _cross(u: tuple, v: tuple) -> tuple[float, float, float]:
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def _place_sg(n: dict, ca: dict, cb: dict, chi1: float) -> dict:
    """Place a Sγ from backbone N, Cα, Cβ at dihedral χ1 (N-Cα-Cβ-Sγ).

    Standard internal-coordinate ("NeRF") placement at fixed Cβ-Sγ bond length
    and Cα-Cβ-Sγ angle; only χ1 varies. Sign of χ1 is irrelevant here because
    callers scan the full 0-360° range.
    """
    A, B, C = _xyz(n), _xyz(ca), _xyz(cb)
    bc = _unit((C[0] - B[0], C[1] - B[1], C[2] - B[2]))
    nrm = _unit(_cross((B[0] - A[0], B[1] - A[1], B[2] - A[2]), bc))
    m2 = _cross(nrm, bc)
    theta = CA_CB_SG_ANGLE
    d = (
        -CB_SG_BOND * math.cos(theta),
        CB_SG_BOND * math.sin(theta) * math.cos(chi1),
        CB_SG_BOND * math.sin(theta) * math.sin(chi1),
    )
    return {
        "atom": "SG",
        "x": C[0] + bc[0] * d[0] + m2[0] * d[1] + nrm[0] * d[2],
        "y": C[1] + bc[1] * d[0] + m2[1] * d[1] + nrm[1] * d[2],
        "z": C[2] + bc[2] * d[0] + m2[2] * d[1] + nrm[2] * d[2],
    }


def _rotamer_sg_positions(res_atoms: dict, all_atoms: list[dict], own_key: tuple[str, int]) -> list[dict]:
    """Clash-free Sγ candidate positions over a χ1 scan; [] if backbone incomplete."""
    n, ca, cb = res_atoms.get("N"), res_atoms.get("CA"), res_atoms.get("CB")
    if not (n and ca and cb):
        return []
    out: list[dict] = []
    for deg in range(0, 360, ROTAMER_STEP):
        sg = _place_sg(n, ca, cb, math.radians(deg))
        # Allow contact with the partner (they bridge via As); forbid only
        # clashes of the rotated Sγ with the rest of the protein.
        if _clash_count(sg, all_atoms, {own_key}) <= CLASH_TOLERANCE:
            out.append(sg)
    return out


def _rotamer_bidentate(
    anchor_atoms: dict,
    partner_atoms: dict,
    anchor_key: tuple[str, int],
    partner_key: tuple[str, int],
    all_atoms: list[dict],
) -> tuple[float | None, bool]:
    """Minimum achievable Sγ-Sγ distance over χ1 rotamers of both cysteines, and
    whether some clash-free rotamer pair admits an As(III) bridge in the window."""
    a_sgs = _rotamer_sg_positions(anchor_atoms, all_atoms, anchor_key)
    p_sgs = _rotamer_sg_positions(partner_atoms, all_atoms, partner_key)
    if not a_sgs or not p_sgs:
        return None, False
    best: float | None = None
    feasible = False
    for sa in a_sgs:
        for sp in p_sgs:
            d = _dist(sa, sp)
            if best is None or d < best:
                best = d
            if BIDENTATE_RANGE[0] <= d <= BIDENTATE_RANGE[1] and not feasible:
                as_pos = _bidentate_as_position(sa, sp)
                if as_pos is not None and _clash_count(as_pos, all_atoms, {anchor_key, partner_key}) <= CLASH_TOLERANCE:
                    feasible = True
    return best, feasible


def _bidentate_as_position(sg1: dict, sg2: dict) -> dict | None:
    """As position at apex of isoceles triangle with Sγ-Sγ as base.

    Returns None if Sγ-Sγ distance > 2·AS_S_BOND (geometrically infeasible —
    As cannot reach both sulfurs at the As-S bond length).
    """
    d = _dist(sg1, sg2)
    half_d = d / 2
    if half_d >= AS_S_BOND:
        return None
    h = math.sqrt(AS_S_BOND ** 2 - half_d ** 2)
    mid = {
        "x": (sg1["x"] + sg2["x"]) / 2,
        "y": (sg1["y"] + sg2["y"]) / 2,
        "z": (sg1["z"] + sg2["z"]) / 2,
    }
    sg_axis = _unit(_vec(sg1, sg2))
    # Pick an arbitrary reference not parallel to the Sγ-Sγ axis.
    ref = (1.0, 0.0, 0.0) if abs(sg_axis[0]) < 0.9 else (0.0, 1.0, 0.0)
    perp = _unit(
        (
            sg_axis[1] * ref[2] - sg_axis[2] * ref[1],
            sg_axis[2] * ref[0] - sg_axis[0] * ref[2],
            sg_axis[0] * ref[1] - sg_axis[1] * ref[0],
        )
    )
    return {
        "x": mid["x"] + perp[0] * h,
        "y": mid["y"] + perp[1] * h,
        "z": mid["z"] + perp[2] * h,
    }


def _clash_count(pos: dict, all_atoms: list[dict], skip_resi: set[tuple[str, int]]) -> int:
    n = 0
    for a in all_atoms:
        if (a["chain"], a["resi"]) in skip_resi:
            continue
        if a["atom"].startswith("H"):
            continue
        if _dist(pos, a) < CLASH_THR:
            n += 1
    return n


def _parse_partners(vicinal_str: str) -> list[tuple[tuple[str, int], float]]:
    """Parse 'A35:3.92;A38:5.10' -> [(('A', 35), 3.92), (('A', 38), 5.10)]"""
    out: list[tuple[tuple[str, int], float]] = []
    if not vicinal_str:
        return out
    for tok in vicinal_str.split(";"):
        if not tok:
            continue
        try:
            tag, dstr = tok.split(":")
            ch, resi = tag[0], int(tag[1:])
            out.append(((ch, resi), float(dstr)))
        except (ValueError, IndexError):
            continue
    return out


def analyse(gene: str, cys_rows: list[dict]) -> list[dict]:
    pdb = PREP_DIR / f"{gene}.clean.pdb"
    if not pdb.exists():
        return []
    atoms = _parse_atoms(pdb)
    sg = {(a["chain"], a["resi"]): a for a in atoms if a["resn"] == "CYS" and a["atom"] == "SG"}
    # N/Cα/Cβ (and Sγ) per Cys, for the χ1 rotamer bidentate scan.
    cys_res: dict[tuple[str, int], dict[str, dict]] = {}
    for a in atoms:
        if a["resn"] == "CYS":
            cys_res.setdefault((a["chain"], a["resi"]), {})[a["atom"]] = a

    out: list[dict] = []
    for row in cys_rows:
        if row["reactive"] != "Y":
            continue
        key = (row["chain"], int(row["resi"]))
        if key not in sg:
            continue
        anchor_sg = sg[key]
        partners = _parse_partners(row.get("vicinal_partners", ""))

        # Bidentate feasibility: Sγ-Sγ in geometric window, As-apex placed
        # between sulfurs, clash count below tolerance.
        di_feas = 0
        di_clash_min = 0
        partner_tags: list[str] = []
        for pkey, d_sg_sg in partners:
            if not (BIDENTATE_RANGE[0] <= d_sg_sg <= BIDENTATE_RANGE[1]):
                continue
            psg = sg.get(pkey)
            if psg is None:
                continue
            as_pos = _bidentate_as_position(anchor_sg, psg)
            if as_pos is None:
                continue
            clashes = _clash_count(as_pos, atoms, {key, pkey})
            if clashes <= CLASH_TOLERANCE:
                di_feas += 1
                partner_tags.append(f"{pkey[0]}{pkey[1]}")
                di_clash_min = max(di_clash_min, clashes)

        # Tridentate: 3 Sγ forming a triangle with all edges ≤ TRIDENTATE_MAX.
        tri_feas = 0
        for i in range(len(partners)):
            for j in range(i + 1, len(partners)):
                p1_key, d1 = partners[i]
                p2_key, d2 = partners[j]
                if d1 > TRIDENTATE_MAX or d2 > TRIDENTATE_MAX:
                    continue
                p1_sg = sg.get(p1_key)
                p2_sg = sg.get(p2_key)
                if p1_sg is None or p2_sg is None:
                    continue
                if _dist(p1_sg, p2_sg) <= TRIDENTATE_MAX:
                    tri_feas += 1

        # Induced-fit bidentate: does a χ1 rotation of the anchor Cys and any
        # vicinal partner bring their Sγ into the bidentate window (clash-free)?
        # Considers ALL vicinal partners, including ones outside the static
        # window (e.g. the METTL3 375/376 pair at ~6.75 Å).
        rot_min: float | None = None
        rot_feasible = False
        anchor_res = cys_res.get(key, {})
        for pkey, _d in partners:
            partner_res = cys_res.get(pkey)
            if not partner_res:
                continue
            d_min, feas = _rotamer_bidentate(anchor_res, partner_res, key, pkey, atoms)
            if d_min is not None and (rot_min is None or d_min < rot_min):
                rot_min = d_min
            rot_feasible = rot_feasible or feas

        # Functional proximity from cys_table.
        fp_str = row.get("functional_proximity_A", "NA")
        try:
            fp = float(fp_str) if fp_str != "NA" else None
        except ValueError:
            fp = None

        pka_str = row.get("pka", "NA")
        try:
            pka_val = float(pka_str) if pka_str != "NA" else None
        except ValueError:
            pka_val = None
        thiolate_bonus = 1.0 if (pka_val is not None and pka_val <= 7.5) else 0.0

        # Active-site proximity bonus — biologically the strongest signal.
        func_bonus = 0.0
        if fp is not None:
            if fp <= 6.0:
                func_bonus = 2.0
            elif fp <= 12.0:
                func_bonus = 1.0

        score = (
            3.0 * tri_feas
            + 2.0 * di_feas
            + func_bonus
            + 0.5 * thiolate_bonus
            - 0.5 * di_clash_min
        )

        # Best binding mode for this anchor (informs the REPORT).
        if tri_feas > 0:
            mode = "tridentate"
        elif di_feas > 0:
            mode = "bidentate"
        else:
            mode = "monodentate"

        out.append(
            {
                "gene": gene,
                "anchor_chain": row["chain"],
                "anchor_resi": row["resi"],
                "n_vicinal_partners": row["n_vicinal_partners"],
                "bidentate_feasible": di_feas,
                "tridentate_feasible": tri_feas,
                "binding_mode": mode,
                "partner_tags": ";".join(partner_tags),
                "functional_proximity_A": fp_str,
                "thiolate_bonus": thiolate_bonus,
                "covalent_score": round(score, 3),
                "min_rotamer_sg_sg_A": f"{rot_min:.2f}" if rot_min is not None else "NA",
                "bidentate_feasible_rotamer": 1 if rot_feasible else 0,
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
            "binding_mode",
            "partner_tags",
            "functional_proximity_A",
            "thiolate_bonus",
            "covalent_score",
            "min_rotamer_sg_sg_A",
            "bidentate_feasible_rotamer",
        ],
    )
    LOG.info("wrote %d covalent anchor rows -> %s", len(rows), OUT_TSV)
    return 0


if __name__ == "__main__":
    sys.exit(main())
