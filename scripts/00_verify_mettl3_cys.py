"""Stage 0 (diagnostic): resolve the METTL3 Cys375-vs-Cys376 numbering question.

Motivation
----------
The pipeline's headline hit is arsenite binding a reactive cysteine in METTL3,
reported as **Cys375** (chain A, 5IL0 author numbering). The literature instead
identifies **Cys376** — adjacent to the SAM-binding Asp377 — as the functionally
decisive SAM-pocket cysteine (covalent modification there, e.g. S-palmitoylation,
lowers SAM binding and activity; a C376S mutant abolishes the effect). Because the
pipeline does **no** residue renumbering (`02_prepare_receptors.py` /
`03_cys_landscape.py` copy `line[22:26]` verbatim), "Cys375" is raw author numbering
and may be a one-residue offset from the UniProt/literature Cys376 — or 375 and 376
may be two genuinely adjacent cysteines (a vicinal dithiol, which As(III) favours).

This script answers that empirically, against the *actual* deposited coordinates,
using only the standard library so it runs without the docking conda env. It is
read-only: it downloads/reads structures and writes a diagnostic table + summary.

What it does (offset-robust — it locates motifs by sequence, not by trusting numbers)
------------------------------------------------------------------------------------
For each METTL3 structure (primary + alternates in data/targets.tsv) it:
  1. Reads DBREF records to recover the author->UniProt numbering map (detects any
     off-by-one directly) and SEQADV records for engineered mutations near the site.
  2. Locates the DPPW catalytic motif and any Cys-immediately-followed-by-Asp pair by
     scanning the modeled sequence, so the analysis holds whatever the numbering.
  3. Inventories every Cys near the site (author 370-385), reports its UniProt number,
     and measures its Sγ distance to: the bound SAM/SAH cofactor (if present), the
     adjacent residue, the Cys-adjacent Asp, and the DPPW motif — plus pairwise
     Sγ-Sγ distances between neighbouring cysteines.
  4. Emits a verdict: `off_by_one`, `distinct_adjacent`, or `single_cys`, and which
     Cys is the true SAM-pocket residue.

Outputs
-------
  results/mettl3_cys_verification.tsv  — one row per candidate Cys per structure
  results/mettl3_cys_verification.md   — human-readable summary + verdict

Run after 01_fetch_structures.py (needs data/pdb/<PDBID>.pdb). Falls back to the
prepared receptor for coordinates if raw PDBs are absent (DBREF/SEQRES then unknown).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import yaml

from _common import DATA, PDB_DIR, PREP_DIR, RESULTS, get_logger, load_targets

LOG = get_logger("verify_mettl3")
SITES_YAML = DATA / "functional_sites.yaml"
OUT_TSV = RESULTS / "mettl3_cys_verification.tsv"
OUT_MD = RESULTS / "mettl3_cys_verification.md"

GENE = "METTL3"
SITE_WINDOW = range(370, 386)  # author residue numbers to inventory
COFACTOR_RESN = {"SAM", "SAH"}
DPPW = "DPPW"

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U",
}


def _dist(a: dict, b: dict) -> float:
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def _min_dist(src: list[dict], tgt: list[dict]) -> float | None:
    """Minimum pairwise heavy-atom distance between two atom sets."""
    src = [a for a in src if not a["atom"].startswith("H")]
    tgt = [a for a in tgt if not a["atom"].startswith("H")]
    if not src or not tgt:
        return None
    return min(_dist(a, b) for a in src for b in tgt)


def _parse_pdb(path: Path) -> dict:
    """Parse ATOM/HETATM + DBREF + SEQADV. Returns dict with raw records."""
    atoms: list[dict] = []
    dbref: list[dict] = []
    seqadv: list[dict] = []
    for line in path.read_text(errors="replace").splitlines():
        rec = line[:6]
        if rec in ("ATOM  ", "HETATM"):
            try:
                atoms.append(
                    {
                        "hetatm": rec == "HETATM",
                        "atom": line[12:16].strip(),
                        "altloc": line[16].strip(),
                        "resn": line[17:20].strip(),
                        "chain": line[21],
                        "resi": int(line[22:26]),
                        "icode": line[26].strip(),
                        "x": float(line[30:38]),
                        "y": float(line[38:46]),
                        "z": float(line[46:54]),
                    }
                )
            except (ValueError, IndexError):
                continue
        elif rec == "DBREF ":
            try:
                dbref.append(
                    {
                        "chain": line[12].strip(),
                        "seq_begin": int(line[14:18]),
                        "seq_end": int(line[20:24]),
                        "database": line[26:32].strip(),
                        "db_acc": line[33:41].strip(),
                        "db_begin": int(line[55:60]),
                        "db_end": int(line[62:67]),
                    }
                )
            except (ValueError, IndexError):
                continue
        elif rec == "SEQADV":
            try:
                seqadv.append(
                    {
                        "resn": line[12:15].strip(),
                        "chain": line[16].strip(),
                        "resi": line[18:22].strip(),
                        "db_acc": line[29:38].strip(),
                        "note": line[49:].strip(),
                    }
                )
            except (ValueError, IndexError):
                continue
    return {"atoms": atoms, "dbref": dbref, "seqadv": seqadv}


def _residues(atoms: list[dict], chain: str) -> dict[int, dict]:
    """Map author resi -> {'resn', 'atoms': [...]} for one chain (protein only).

    Keeps the highest-occupancy alt-loc implicitly by taking all; distances use
    min over atoms so alt-locs do not bias the result.
    """
    out: dict[int, dict] = {}
    for a in atoms:
        if a["hetatm"] or a["chain"] != chain:
            continue
        r = out.setdefault(a["resi"], {"resn": a["resn"], "atoms": []})
        r["atoms"].append(a)
    return out


def _chain_sequence(residues: dict[int, dict]) -> list[tuple[int, str]]:
    """Ordered (resi, one-letter) for a chain, unknown residues as 'X'."""
    return [(ri, AA3TO1.get(residues[ri]["resn"], "X")) for ri in sorted(residues)]


def _find_motif(seq: list[tuple[int, str]], motif: str) -> int | None:
    """Return the author resi of the first residue of `motif`, or None."""
    letters = "".join(c for _, c in seq)
    idx = letters.find(motif)
    return seq[idx][0] if idx != -1 else None


def _dbref_for_chain(dbref: list[dict], chain: str) -> dict | None:
    for d in dbref:
        if d["chain"] == chain and d["database"].upper() in ("UNP", "UNIPROT"):
            return d
    # fall back to any dbref on the chain
    for d in dbref:
        if d["chain"] == chain:
            return d
    return None


def _uniprot_of(author_resi: int, ref: dict | None) -> int | None:
    """author resi -> UniProt resi using the DBREF linear map."""
    if not ref:
        return None
    offset = ref["db_begin"] - ref["seq_begin"]
    return author_resi + offset


def _cofactor_atoms(atoms: list[dict]) -> tuple[list[dict], str]:
    cof = [a for a in atoms if a["hetatm"] and a["resn"] in COFACTOR_RESN]
    label = ",".join(sorted({a["resn"] for a in cof})) if cof else "none"
    return cof, label


def analyse_structure(pdbid: str, chain: str, pdb_path: Path) -> tuple[list[dict], list[str]]:
    parsed = _parse_pdb(pdb_path)
    residues = _residues(parsed["atoms"], chain)
    if not residues:
        # try any single chain if the configured one is absent
        chains = sorted({a["chain"] for a in parsed["atoms"] if not a["hetatm"]})
        return [], [f"- **{pdbid}**: chain {chain} not found (chains present: {chains or 'none'})."]

    seq = _chain_sequence(residues)
    ref = _dbref_for_chain(parsed["dbref"], chain)
    offset = (ref["db_begin"] - ref["seq_begin"]) if ref else None
    dppw_author = _find_motif(seq, DPPW)
    cof_atoms, cof_label = _cofactor_atoms(parsed["atoms"])

    # DPPW anchor atoms (the 4 residues of the motif, author numbering)
    dppw_atoms: list[dict] = []
    if dppw_author is not None:
        for ri in range(dppw_author, dppw_author + 4):
            if ri in residues:
                dppw_atoms.extend(residues[ri]["atoms"])

    # Locate every Cys-immediately-followed-by-Asp pair in the window (the
    # literature Cys376/Asp377 motif — offset-robust).
    cys_asp_pairs = [
        ri for ri, c in seq
        if c == "C" and (ri + 1) in residues and residues[ri + 1]["resn"] == "ASP"
    ]

    rows: list[dict] = []
    cys_in_window = [ri for ri in residues if ri in SITE_WINDOW and residues[ri]["resn"] == "CYS"]
    for ri in sorted(cys_in_window):
        res = residues[ri]
        sg = next((a for a in res["atoms"] if a["atom"] == "SG"), None)
        sg_list = [sg] if sg else []
        adj = residues.get(ri + 1)
        adj_resn = adj["resn"] if adj else "NA"
        # nearest Cys neighbour (Sγ-Sγ) within +/- a few residues
        nn_tag, nn_d = "NA", None
        for rj in cys_in_window:
            if rj == ri:
                continue
            osg = next((a for a in residues[rj]["atoms"] if a["atom"] == "SG"), None)
            if sg and osg:
                d = _dist(sg, osg)
                if nn_d is None or d < nn_d:
                    nn_d, nn_tag = d, f"{chain}{rj}"
        rows.append(
            {
                "pdb": pdbid,
                "chain": chain,
                "author_resi": ri,
                "uniprot_resi": _uniprot_of(ri, ref) if ref else "NA",
                "resn": res["resn"],
                "has_sg": "Y" if sg else "N",
                "adjacent_resi": ri + 1,
                "adjacent_resn": adj_resn,
                "is_cys_asp_motif": "Y" if ri in cys_asp_pairs else "N",
                "d_sg_cofactor_A": f"{_min_dist(sg_list, cof_atoms):.2f}" if (sg and cof_atoms) else "NA",
                "d_sg_adjacent_A": f"{_min_dist(sg_list, adj['atoms']):.2f}" if (sg and adj) else "NA",
                "d_sg_dppw_A": f"{_min_dist(sg_list, dppw_atoms):.2f}" if (sg and dppw_atoms) else "NA",
                "nearest_cys": nn_tag,
                "nearest_cys_sg_sg_A": f"{nn_d:.2f}" if nn_d is not None else "NA",
                "cofactor_in_struct": cof_label,
                "author_to_uniprot_offset": offset if offset is not None else "NA",
                "dppw_author_start": dppw_author if dppw_author is not None else "NA",
            }
        )

    # Per-structure narrative
    notes: list[str] = []
    off_txt = (
        f"author→UniProt offset **{offset:+d}** (author {ref['seq_begin']}–{ref['seq_end']} "
        f"↦ UniProt {ref['db_begin']}–{ref['db_end']}, {ref['database']} {ref['db_acc']})"
        if ref else "no DBREF UniProt mapping in file"
    )
    notes.append(f"- **{pdbid}** (chain {chain}): {off_txt}; cofactor in structure: `{cof_label}`; "
                 f"DPPW motif author start: {dppw_author}.")
    for ri in sorted(cys_in_window):
        up = _uniprot_of(ri, ref)
        tag = f"author Cys{ri}" + (f" = UniProt Cys{up}" if up is not None else "")
        motif = " — **Cys–Asp motif (matches lit. Cys376/Asp377)**" if ri in cys_asp_pairs else ""
        notes.append(f"    - {tag}{motif}")
    seqadv_hits = [s for s in parsed["seqadv"]
                   if s["chain"] == chain and s["resi"].isdigit() and int(s["resi"]) in SITE_WINDOW]
    for s in seqadv_hits:
        notes.append(f"    - SEQADV: {s['resn']} {s['chain']}{s['resi']} — {s['note']}")
    return rows, notes


def _verdict(all_rows: list[dict]) -> list[str]:
    """Aggregate a headline verdict across structures."""
    out: list[str] = ["", "## Verdict", ""]
    if not all_rows:
        out.append("No METTL3 cysteines resolved in the site window — check that structures were fetched "
                   "(`scripts/01_fetch_structures.py`) and that the METTL3 chain id in functional_sites.yaml "
                   "matches the deposited file.")
        return out

    # Look at the primary structure's window (first pdb encountered).
    prim = all_rows[0]["pdb"]
    prim_rows = [r for r in all_rows if r["pdb"] == prim]
    cys_authors = sorted(int(r["author_resi"]) for r in prim_rows)
    motif_rows = [r for r in prim_rows if r["is_cys_asp_motif"] == "Y"]
    offsets = {r["author_to_uniprot_offset"] for r in prim_rows if r["author_to_uniprot_offset"] != "NA"}
    offset = next(iter(offsets)) if len(offsets) == 1 else None

    both_375_376 = 375 in cys_authors and 376 in cys_authors
    if both_375_376:
        out.append("**distinct_adjacent** — author positions 375 *and* 376 are both cysteines: a genuine "
                   "vicinal dithiol. The pipeline picked 375; the literature/SAM-pocket residue is 376. "
                   "Re-dock at 376 and test the bidentate As(III) bridge across the 375–376 pair "
                   "(`scripts/06_dock_covalent.py`).")
    elif offset not in (None, 0):
        out.append(f"**off_by_one** — this file's author numbering is offset from UniProt by {offset:+}. "
                   f"The pipeline's author-'Cys375' therefore corresponds to UniProt Cys{375 + int(offset)}. "
                   "If that equals 376, our result and the literature describe the *same* thiol under "
                   "different numbering — relabel to the UniProt/literature Cys376 and re-anchor to Asp377.")
    else:
        out.append("**single_cys / aligned numbering** — author and UniProt numbering agree here. Confirm "
                   "from the rows below which single cysteine (375 or 376) is actually present and which is "
                   "nearest the SAM/SAH cofactor and the Cys–Asp motif; that is the true SAM-pocket residue.")

    if motif_rows:
        m = motif_rows[0]
        out.append("")
        out.append(f"The Cys–Asp motif (literature Cys376/Asp377) sits at **author Cys{m['author_resi']} "
                   f"= UniProt Cys{m['uniprot_resi']}**, adjacent {m['adjacent_resn']}{m['adjacent_resi']}, "
                   f"Sγ→cofactor {m['d_sg_cofactor_A']} Å, Sγ→DPPW {m['d_sg_dppw_A']} Å.")
    out.append("")
    out.append("→ Use the `author_resi`/`uniprot_resi` and `d_sg_cofactor_A` columns in the TSV to set the "
               "corrected anchor residue and to state the SAM-pocket distance in the thesis.")
    return out


def main() -> int:
    targets = {t.gene: t for t in load_targets()}
    t = targets.get(GENE)
    if t is None:
        LOG.error("%s not in targets.tsv", GENE)
        return 1
    sites = {}
    if SITES_YAML.exists():
        sites = yaml.safe_load(SITES_YAML.read_text()) or {}
    chain = (sites.get(GENE) or {}).get("chain", "A")

    all_rows: list[dict] = []
    md: list[str] = ["# METTL3 Cys375-vs-Cys376 verification", "",
                     f"UniProt: {t.uniprot} · configured chain: {chain} · "
                     f"structures: {', '.join(t.all_pdbs) or '(none)'}", "",
                     "## Per-structure findings", ""]

    scanned_any = False
    for pdbid in t.all_pdbs:
        raw = PDB_DIR / f"{pdbid}.pdb"
        if not raw.exists():
            md.append(f"- **{pdbid}**: not fetched (missing {raw}).")
            continue
        scanned_any = True
        rows, notes = analyse_structure(pdbid, chain, raw)
        all_rows.extend(rows)
        md.extend(notes)

    if not scanned_any:
        prepared = PREP_DIR / f"{GENE}.clean.pdb"
        if prepared.exists():
            LOG.warning("no raw PDBs — using prepared receptor (DBREF/UniProt mapping unavailable)")
            rows, notes = analyse_structure(f"{GENE}.clean", chain, prepared)
            all_rows.extend(rows)
            md.extend(notes)
        else:
            LOG.error("no structures found. Run scripts/01_fetch_structures.py first.")

    md.extend(_verdict(all_rows))

    fieldnames = [
        "pdb", "chain", "author_resi", "uniprot_resi", "resn", "has_sg",
        "adjacent_resi", "adjacent_resn", "is_cys_asp_motif",
        "d_sg_cofactor_A", "d_sg_adjacent_A", "d_sg_dppw_A",
        "nearest_cys", "nearest_cys_sg_sg_A", "cofactor_in_struct",
        "author_to_uniprot_offset", "dppw_author_start",
    ]
    RESULTS.mkdir(parents=True, exist_ok=True)
    from _common import write_tsv

    write_tsv(all_rows, OUT_TSV, fieldnames)
    OUT_MD.write_text("\n".join(md) + "\n")
    LOG.info("wrote %s (%d rows) and %s", OUT_TSV, len(all_rows), OUT_MD)
    for line in _verdict(all_rows):
        if line.strip():
            LOG.info("%s", line.replace("**", "").replace("`", ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
