# METTL3 Cys375-vs-Cys376 verification

UniProt: Q86U44 · configured chain: A · structures: 5IL0, 5IL1, 5K7M, 5K7U, 5K7W, 6TTS, 6TTT

## Per-structure findings

- **5IL0** (chain A): author→UniProt offset **+0** (author 369–580 ↦ UniProt 369–580, UNP Q86U44); cofactor in structure: `none`; DPPW motif author start: 395.
    - author Cys375 = UniProt Cys375
    - author Cys376 = UniProt Cys376 — **Cys–Asp motif (matches lit. Cys376/Asp377)**
- **5IL1** (chain A): author→UniProt offset **+0** (author 369–580 ↦ UniProt 369–580, UNP Q86U44); cofactor in structure: `SAM`; DPPW motif author start: 395.
    - author Cys375 = UniProt Cys375
    - author Cys376 = UniProt Cys376 — **Cys–Asp motif (matches lit. Cys376/Asp377)**
- **5K7M** (chain A): author→UniProt offset **+0** (author 357–580 ↦ UniProt 357–580, UNP Q86U44); cofactor in structure: `none`; DPPW motif author start: 395.
    - author Cys375 = UniProt Cys375
    - author Cys376 = UniProt Cys376 — **Cys–Asp motif (matches lit. Cys376/Asp377)**
- **5K7U** (chain A): author→UniProt offset **+0** (author 357–580 ↦ UniProt 357–580, UNP Q86U44); cofactor in structure: `SAM`; DPPW motif author start: 395.
    - author Cys375 = UniProt Cys375
    - author Cys376 = UniProt Cys376 — **Cys–Asp motif (matches lit. Cys376/Asp377)**
- **5K7W** (chain A): author→UniProt offset **+0** (author 357–580 ↦ UniProt 357–580, UNP Q86U44); cofactor in structure: `SAH`; DPPW motif author start: 395.
    - author Cys375 = UniProt Cys375
    - author Cys376 = UniProt Cys376 — **Cys–Asp motif (matches lit. Cys376/Asp377)**
- **6TTS** (chain A): no DBREF UniProt mapping in file; cofactor in structure: `none`; DPPW motif author start: None.
- **6TTT** (chain A): author→UniProt offset **+0** (author 1–580 ↦ UniProt 1–580, UNP Q86U44); cofactor in structure: `none`; DPPW motif author start: 395.
    - author Cys375 = UniProt Cys375
    - author Cys376 = UniProt Cys376 — **Cys–Asp motif (matches lit. Cys376/Asp377)**

## Verdict

**distinct_adjacent** — author positions 375 *and* 376 are both cysteines: a genuine vicinal dithiol. The pipeline picked 375; the literature/SAM-pocket residue is 376. Re-dock at 376 and test the bidentate As(III) bridge across the 375–376 pair (`scripts/06_dock_covalent.py`).

The Cys–Asp motif (literature Cys376/Asp377) sits at **author Cys376 = UniProt Cys376**, adjacent ASP377, Sγ→cofactor NA Å, Sγ→DPPW 11.96 Å.

→ Use the `author_resi`/`uniprot_resi` and `d_sg_cofactor_A` columns in the TSV to set the corrected anchor residue and to state the SAM-pocket distance in the thesis.
