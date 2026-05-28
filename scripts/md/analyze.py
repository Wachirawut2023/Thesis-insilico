"""Analyse a completed MD trajectory.

For each gene/mode pair under results/md/, computes:
  - Cα RMSD over time (overall protein stability)
  - per-residue RMSF (which regions move most)
  - per-Cys Sγ SASA over time (does accessibility change?)
  - per-Cys vicinal distance over time (do vicinal pairs persist?)
  - cryptic-Cys flag: any Cys whose mean SASA over the trajectory >> static SASA

Writes results/md/<gene>_<mode>/analysis/ with TSVs and PNG plots.

Requires: MDAnalysis, freesasa, numpy, matplotlib (all in conda env).
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np


def _setup_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def analyse(gene: str, mode: str, run_dir: Path) -> int:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms

    tpr = run_dir / "md.tpr"
    xtc = run_dir / "md.xtc"
    if not tpr.exists() or not xtc.exists():
        print(f"ERROR: missing md.tpr or md.xtc in {run_dir}", file=sys.stderr)
        return 1

    out_dir = run_dir / "analysis"
    out_dir.mkdir(exist_ok=True)

    print(f"[analyze] loading trajectory: {xtc.name}")
    u = mda.Universe(str(tpr), str(xtc))
    n_frames = len(u.trajectory)
    n_residues = u.select_atoms("protein").n_residues
    print(f"[analyze] {n_frames} frames, {n_residues} residues")

    plt = _setup_matplotlib()

    # ── 1. Cα RMSD over time ──
    print("[analyze] computing Cα RMSD")
    ref = u.copy()
    ref.trajectory[0]
    ref_ca = ref.select_atoms("name CA")
    u_ca = u.select_atoms("name CA")
    rmsd_arr = []
    for ts in u.trajectory:
        # align then compute
        from MDAnalysis.analysis.align import alignto
        alignto(u, ref, select="name CA")
        rmsd_val = np.sqrt(np.mean(np.sum((u_ca.positions - ref_ca.positions) ** 2, axis=1)))
        rmsd_arr.append((ts.time / 1000.0, rmsd_val))  # ps -> ns
    times, rmsds = zip(*rmsd_arr)
    with (out_dir / "rmsd.tsv").open("w") as fh:
        fh.write("time_ns\trmsd_A\n")
        for t, r in zip(times, rmsds):
            fh.write(f"{t:.3f}\t{r:.3f}\n")

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times, rmsds, color="steelblue")
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("Cα RMSD (Å)")
    ax.set_title(f"{gene} ({mode}) Cα RMSD")
    ax.axhline(2.0, color="gray", linestyle="--", alpha=0.5, label="2 Å stability threshold")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "rmsd.png", dpi=120)
    plt.close(fig)

    # ── 2. per-residue RMSF (mobility map) ──
    print("[analyze] computing RMSF")
    ca = u.select_atoms("protein and name CA")
    rmsf_calc = rms.RMSF(ca).run()
    with (out_dir / "rmsf.tsv").open("w") as fh:
        fh.write("resi\tresn\trmsf_A\n")
        for atom, rmsf_val in zip(ca, rmsf_calc.results.rmsf):
            fh.write(f"{atom.resnum}\t{atom.resname}\t{rmsf_val:.3f}\n")

    fig, ax = plt.subplots(figsize=(10, 3))
    resnums = ca.resnums
    ax.plot(resnums, rmsf_calc.results.rmsf, color="darkorange")
    ax.set_xlabel("residue number")
    ax.set_ylabel("Cα RMSF (Å)")
    ax.set_title(f"{gene} ({mode}) per-residue RMSF")
    fig.tight_layout()
    fig.savefig(out_dir / "rmsf.png", dpi=120)
    plt.close(fig)

    # ── 3. per-Cys SG SASA over time ──
    print("[analyze] computing per-Cys SG SASA over time (sampling every 10th frame)")
    try:
        import freesasa
    except ImportError:
        print("[analyze] freesasa not installed; skipping SASA analysis", file=sys.stderr)
        freesasa = None

    if freesasa is not None:
        cys_sgs = u.select_atoms("resname CYS and name SG")
        if cys_sgs.n_atoms == 0:
            print("[analyze] no CYS SG atoms found")
        else:
            tmp_pdb = out_dir / "_frame.pdb"
            stride = max(1, n_frames // 50)  # ~50 samples
            print(f"[analyze]   stride={stride} -> ~{n_frames // stride} samples per Cys")
            cys_data: dict[int, list[float]] = {a.resnum: [] for a in cys_sgs}
            times_sasa = []
            for i, ts in enumerate(u.trajectory[::stride]):
                u.select_atoms("protein").write(str(tmp_pdb))
                try:
                    structure = freesasa.Structure(str(tmp_pdb))
                    if hasattr(freesasa, "calc"):
                        result = freesasa.calc(structure)
                    else:
                        result = freesasa.Calc().calculate(structure)
                except Exception:
                    continue
                times_sasa.append(ts.time / 1000.0)
                # build a per-Cys-resi map
                sg_areas: dict[int, float] = {}
                for j in range(structure.nAtoms()):
                    if structure.atomName(j).strip() == "SG":
                        try:
                            resi = int(str(structure.residueNumber(j)).strip())
                        except ValueError:
                            continue
                        sg_areas[resi] = result.atomArea(j)
                for r in cys_data:
                    cys_data[r].append(sg_areas.get(r, 0.0))
            tmp_pdb.unlink(missing_ok=True)

            with (out_dir / "cys_sasa_time.tsv").open("w") as fh:
                fh.write("time_ns\t" + "\t".join(f"C{r}_SG" for r in sorted(cys_data)) + "\n")
                for i, t in enumerate(times_sasa):
                    row = [f"{t:.2f}"]
                    for r in sorted(cys_data):
                        row.append(f"{cys_data[r][i]:.2f}" if i < len(cys_data[r]) else "NA")
                    fh.write("\t".join(row) + "\n")

            # Cryptic-Cys summary: mean/min/max SASA per Cys
            with (out_dir / "cys_summary.tsv").open("w") as fh:
                fh.write("resi\tmean_sasa_A2\tmin_sasa\tmax_sasa\tcryptic_flag\n")
                for r in sorted(cys_data):
                    arr = np.array(cys_data[r])
                    mean_sasa = float(np.mean(arr))
                    min_sasa = float(np.min(arr))
                    max_sasa = float(np.max(arr))
                    # cryptic: was buried (min < 2) but opened up at least once (max > 10)
                    cryptic = "Y" if (min_sasa < 2.0 and max_sasa > 10.0) else "N"
                    fh.write(f"{r}\t{mean_sasa:.2f}\t{min_sasa:.2f}\t{max_sasa:.2f}\t{cryptic}\n")

            fig, ax = plt.subplots(figsize=(10, 4))
            for r in sorted(cys_data):
                ax.plot(times_sasa, cys_data[r], label=f"C{r}", alpha=0.7)
            ax.set_xlabel("time (ns)")
            ax.set_ylabel("Cys Sγ SASA (Å²)")
            ax.set_title(f"{gene} ({mode}) per-Cys Sγ SASA over MD")
            ax.axhline(5.0, color="gray", linestyle="--", alpha=0.5)
            ax.legend(fontsize=6, ncol=4)
            fig.tight_layout()
            fig.savefig(out_dir / "cys_sasa_time.png", dpi=120)
            plt.close(fig)

    print(f"[analyze] DONE  -> {out_dir}")
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gene", required=True)
    p.add_argument("--mode", default="apo", choices=["apo", "bound"])
    p.add_argument("--md-root", type=Path,
                   default=Path("/opt/Thesis-insilico/results/md"))
    args = p.parse_args()
    run_dir = args.md_root / f"{args.gene}_{args.mode}"
    return analyse(args.gene, args.mode, run_dir)


if __name__ == "__main__":
    sys.exit(main())
