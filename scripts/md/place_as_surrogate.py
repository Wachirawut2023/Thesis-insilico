"""Place a Mg2+ surrogate at the predicted As(III) binding position in a GRO
file, and append distance restraints to the GROMACS topology so the
surrogate stays at the predicted bonding distance to the target Cys SG
atom(s) during the simulation.

Rationale: As(III) has no parameters in standard AMBER/CHARMM force
fields. Deriving them via QM/RESP is a multi-day project per protein.
For the question we want MD to answer ("is the binding pose stable?")
the cleanest approximation is to use Mg2+ (similar ionic radius, well-
parameterised) plus harmonic distance restraints to maintain the As-S
bond length of 2.25 A. This tests pose stability and protein
conformational response; it does not give absolute binding free energy.

Usage:
    place_as_surrogate.py --in-gro <in> --out-gro <out> --topol topol.top \
        --chain A --anchor-resi 113 [--partner-resi 35]

For monodentate: omit --partner-resi. Surrogate is placed 2.25 A from
the anchor Sγ along the (-Cβ → Sγ) direction.

For bidentate: provide --partner-resi. Surrogate is placed at the apex
of the isoceles triangle bridging anchor Sγ and partner Sγ.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

AS_S_BOND = 2.25  # Å
FORCE_CONST = 1000  # kJ/mol/nm^2 — harmonic distance restraint stiffness


def parse_gro(path: Path) -> tuple[list[str], list[dict], str]:
    """Return (header_lines, atom_records, box_line). GRO format is fixed-width."""
    lines = path.read_text().splitlines()
    title = lines[0]
    n_atoms = int(lines[1].strip())
    atoms: list[dict] = []
    for i in range(2, 2 + n_atoms):
        ln = lines[i]
        atoms.append({
            "raw": ln,
            "resi": int(ln[0:5]),
            "resn": ln[5:10].strip(),
            "atom": ln[10:15].strip(),
            "serial": int(ln[15:20]),
            "x": float(ln[20:28]),
            "y": float(ln[28:36]),
            "z": float(ln[36:44]),
        })
    box = lines[2 + n_atoms]
    return [title, str(n_atoms)], atoms, box


def find_sg(atoms: list[dict], resi: int) -> dict | None:
    for a in atoms:
        if a["resn"] == "CYS" and a["resi"] == resi and a["atom"] == "SG":
            return a
    return None


def find_cb(atoms: list[dict], resi: int) -> dict | None:
    for a in atoms:
        if a["resn"] == "CYS" and a["resi"] == resi and a["atom"] == "CB":
            return a
    return None


def _vec(a: dict, b: dict) -> tuple[float, float, float]:
    return (b["x"] - a["x"], b["y"] - a["y"], b["z"] - a["z"])


def _unit(v: tuple[float, float, float]) -> tuple[float, float, float]:
    n = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    return (v[0] / n, v[1] / n, v[2] / n) if n > 0 else (0.0, 0.0, 0.0)


def _dist(a: dict, b: dict) -> float:
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def monodentate_pos(anchor_sg: dict, anchor_cb: dict) -> tuple[float, float, float]:
    """Place surrogate at 2.25 Å from Sγ, opposite the Cβ-Sγ direction (Å -> nm)."""
    d = _unit(_vec(anchor_cb, anchor_sg))
    as_pos_A = (
        anchor_sg["x"] * 10 + d[0] * AS_S_BOND,  # convert Å to position in nm space first
        anchor_sg["y"] * 10 + d[1] * AS_S_BOND,
        anchor_sg["z"] * 10 + d[2] * AS_S_BOND,
    )
    # GRO is in nm; convert back
    return (as_pos_A[0] / 10, as_pos_A[1] / 10, as_pos_A[2] / 10)


def bidentate_pos(sg1: dict, sg2: dict) -> tuple[float, float, float] | None:
    """As at apex of isoceles triangle bridging two Sγ atoms (positions in nm)."""
    d_nm = _dist(sg1, sg2)
    d_A = d_nm * 10
    half_d_A = d_A / 2
    if half_d_A >= AS_S_BOND:
        return None
    h_A = math.sqrt(AS_S_BOND ** 2 - half_d_A ** 2)
    mid = {
        "x": (sg1["x"] + sg2["x"]) / 2,
        "y": (sg1["y"] + sg2["y"]) / 2,
        "z": (sg1["z"] + sg2["z"]) / 2,
    }
    sg_axis = _unit(_vec(sg1, sg2))
    ref = (1.0, 0.0, 0.0) if abs(sg_axis[0]) < 0.9 else (0.0, 1.0, 0.0)
    perp = _unit((
        sg_axis[1] * ref[2] - sg_axis[2] * ref[1],
        sg_axis[2] * ref[0] - sg_axis[0] * ref[2],
        sg_axis[0] * ref[1] - sg_axis[1] * ref[0],
    ))
    h_nm = h_A / 10
    return (
        mid["x"] + perp[0] * h_nm,
        mid["y"] + perp[1] * h_nm,
        mid["z"] + perp[2] * h_nm,
    )


def write_gro(header: list[str], atoms: list[dict], extra: list[dict], box: str, out: Path) -> None:
    new_total = len(atoms) + len(extra)
    lines = [header[0], str(new_total)]
    for a in atoms:
        lines.append(a["raw"])
    last_resi = max((a["resi"] for a in atoms), default=0)
    last_serial = max((a["serial"] for a in atoms), default=0)
    for i, e in enumerate(extra, start=1):
        resi = last_resi + i
        serial = last_serial + i
        # GRO line: %5d%-5s%5s%5d%8.3f%8.3f%8.3f
        line = f"{resi:5d}{'MG':<5s}{'MG':>5s}{serial:5d}{e['x']:8.3f}{e['y']:8.3f}{e['z']:8.3f}"
        lines.append(line)
    lines.append(box)
    out.write_text("\n".join(lines) + "\n")


def append_topol(topol_path: Path, n_mg: int) -> None:
    """Add the MG ion to [ molecules ] in topol.top and ensure ions.itp/atomtypes
    include MG (standard AMBER ion .itp already has MG2)."""
    text = topol_path.read_text()
    if "MG" not in text or "[ molecules ]" in text:
        # Append MG entry to [ molecules ] section
        # GROMACS topol.top usually has lines like "Protein_chain_A    1" / "SOL  12345" / "NA  12"
        new_lines: list[str] = []
        in_mols = False
        added = False
        for ln in text.splitlines():
            new_lines.append(ln)
            if ln.strip().startswith("[ molecules ]"):
                in_mols = True
            elif in_mols and not added and ln.strip() and not ln.strip().startswith(";"):
                # Append after the first non-comment line in molecules section
                pass
        if not added:
            new_lines.append(f"MG               {n_mg}")
        topol_path.write_text("\n".join(new_lines) + "\n")


def write_restraints(topol_dir: Path, sg1_serial: int, sg2_serial: int | None, mg_serial: int) -> None:
    """Append a distance restraint .itp file referencing the Mg surrogate
    and the target SG atom(s). User must #include this in topol.top
    (we print the include line and instruct the caller)."""
    distance_nm = AS_S_BOND / 10
    out = ["[ distance_restraints ]",
           "; ai  aj  type  index  type'  low  up1  up2  fac"]
    out.append(
        f"  {sg1_serial}  {mg_serial}   1   0     1    {distance_nm-0.02:.3f}  {distance_nm+0.02:.3f}  {distance_nm+0.1:.3f}  1.0"
    )
    if sg2_serial is not None:
        out.append(
            f"  {sg2_serial}  {mg_serial}   1   1     1    {distance_nm-0.02:.3f}  {distance_nm+0.02:.3f}  {distance_nm+0.1:.3f}  1.0"
        )
    (topol_dir / "as_restraints.itp").write_text("\n".join(out) + "\n")
    print(f"[place_as_surrogate] wrote {topol_dir/'as_restraints.itp'} — include in topol.top after [moleculetype]")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in-gro", type=Path, required=True)
    p.add_argument("--out-gro", type=Path, required=True)
    p.add_argument("--topol", type=Path, required=True)
    p.add_argument("--chain", type=str, default="A")
    p.add_argument("--anchor-resi", type=int, required=True)
    p.add_argument("--partner-resi", type=str, default="",
                   help="Partner Cys residue number for bidentate (or empty for monodentate)")
    args = p.parse_args()

    header, atoms, box = parse_gro(args.in_gro)
    anchor_sg = find_sg(atoms, args.anchor_resi)
    if anchor_sg is None:
        print(f"ERROR: no CYS {args.anchor_resi} SG atom in {args.in_gro}", file=sys.stderr)
        return 1
    anchor_cb = find_cb(atoms, args.anchor_resi)

    partner_resi = int(args.partner_resi) if args.partner_resi and args.partner_resi != "-" else None
    if partner_resi:
        partner_sg = find_sg(atoms, partner_resi)
        if partner_sg is None:
            print(f"ERROR: no CYS {partner_resi} SG atom (partner) in {args.in_gro}", file=sys.stderr)
            return 1
        pos = bidentate_pos(anchor_sg, partner_sg)
        if pos is None:
            print(f"ERROR: Sγ-Sγ distance > 2 x As-S bond; bidentate infeasible", file=sys.stderr)
            return 1
        partner_serial = partner_sg["serial"]
        mode_str = "bidentate"
    else:
        if anchor_cb is None:
            print(f"ERROR: anchor Cys has no CB atom; cannot orient As", file=sys.stderr)
            return 1
        pos = monodentate_pos(anchor_sg, anchor_cb)
        partner_serial = None
        mode_str = "monodentate"

    print(f"[place_as_surrogate] mode={mode_str} anchor_serial={anchor_sg['serial']}"
          f"{' partner_serial=' + str(partner_serial) if partner_serial else ''}"
          f" Mg_pos=({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f}) nm")

    extra = [{"x": pos[0], "y": pos[1], "z": pos[2]}]
    write_gro(header, atoms, extra, box, args.out_gro)
    append_topol(args.topol, n_mg=1)
    mg_serial = atoms[-1]["serial"] + 1
    write_restraints(args.topol.parent, anchor_sg["serial"], partner_serial, mg_serial)
    print(f"[place_as_surrogate] wrote {args.out_gro}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
