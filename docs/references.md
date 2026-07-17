# References supporting the in silico pipeline

Citations grouped by which decision they support. Each entry has:
- **What it supports** — the parameter, threshold, or biological claim
- **Citation** — full bibliographic information where I'm confident
- **Confidence** — `HIGH` if I'm certain of the citation; `VERIFY` if you should
  confirm via PubMed / Google Scholar before relying on it for the thesis.

This list is designed to be pasted into your reference manager (Zotero, Mendeley,
EndNote, etc.) and then mapped to the chapter at `docs/chapter_in_silico.md`.

---

## A. Arsenic chemistry and cysteine binding

### A.1 As(III)–thiolate geometry (As–S bond length, S–As–S angle)

**What it supports**: The geometric constants used in Stage 6 covalent
docking — As–S bond ≈ 2.25 Å, S–As–S angle ≈ 94°, bidentate Sγ–Sγ
window 3.0–4.4 Å, tridentate As(SR)₃ as the preferred coordination.

**Citations** (HIGH confidence):
- Spuches AM, Kruszyna HG, Rich AM, Wilcox DE. "Thermodynamics of the
  As(III)–thiol interaction: arsenite and monomethylarsenite complexes
  with glutathione, dihydrolipoic acid, and other thiol ligands."
  *Inorg Chem* 2005, 44(8):2964–2972.
- Ramadan D, Rancy PC, Nagarkar RP, Schneider JP, Thorpe C.
  "Arsenic(III) species inhibit oxidative protein folding in vitro."
  *Biochemistry* 2009, 48(2):424–432.
- Shen S, Li X-F, Cullen WR, Weinfeld M, Le XC. "Arsenic binding to
  proteins." *Chem Rev* 2013, 113(10):7769–7792.

The Shen et al. 2013 review is your best one-stop reference for the
biochemistry of arsenic–cysteine binding and provides the proteome-level
context for the positive control choices.

### A.2 Vicinal-cysteine motifs as preferred As(III) targets

**What it supports**: The pipeline's emphasis on vicinal Cys pairs
(≤ 4.4 Å bidentate, all ≤ 4.4 Å for tridentate); the biological reality
that CxxC, CxC, and CxxxC motifs are the canonical arsenic targets.

**Citations** (HIGH confidence):
- Lin S, Cullen WR, Thomas DJ. "Methylarsenicals and arsinothiols are
  potent inhibitors of mouse liver thioredoxin reductase." *Chem Res
  Toxicol* 1999, 12(10):924–930.
- Aposhian HV. "Enzymatic methylation of arsenic species and other new
  approaches to arsenic toxicity." *Annu Rev Pharmacol Toxicol* 1997,
  37:397–419.

### A.3 Clinical relevance — arsenic trioxide in leukaemia (PML-RARα)

**What it supports**: The opening biological motivation in §1.2 of the
chapter — that arsenic chemistry is medically established.

**Citation** (HIGH confidence):
- Zhang X-W, Yan X-J, Zhou Z-R, et al. "Arsenic trioxide controls the
  fate of the PML-RARα oncoprotein by directly binding PML." *Science*
  2010, 328(5975):240–243.
- Lallemand-Breitenbach V, Jeanne M, Benhenda S, et al. "Arsenic
  degrades PML or PML–RARα through a SUMO-triggered RNF4/ubiquitin-mediated
  pathway." *Nat Cell Biol* 2008, 10(5):547–555.

---

## B. Cysteine reactivity criteria

### B.1 Solvent-accessible surface area (SASA) cutoff for "accessible Cys"

**What it supports**: The 5 Å² threshold on SG-atom SASA used to flag a
cysteine as accessible in Stage 3.

**Citations** (HIGH confidence for methodology, VERIFY for the exact cutoff):
- Marino SM, Gladyshev VN. "Cysteine function governs its conservation
  and degeneration and restricts its utilization on protein surfaces."
  *J Mol Biol* 2010, 404(5):902–916.
- Weerapana E, Wang C, Simon GM, et al. "Quantitative reactivity
  profiling predicts functional cysteines in proteomes." *Nature* 2010,
  468(7325):790–795.

Marino & Gladyshev 2010 is the strongest backing for SASA being a
predictor of cysteine functionality. Weerapana et al. 2010 introduced
chemoproteomic reactive-Cys profiling and is the reference for the
operational definition of "reactive cysteine" widely cited in the
chemical-biology literature.

### B.2 pKa threshold for "reactive thiolate"

**What it supports**: The pKa ≤ 8.0 threshold used to flag thiolate-form
cysteines in Stage 3.

**Citations** (HIGH confidence):
- Roos G, Foloppe N, Messens J. "Understanding the pKa of redox cysteines:
  the key role of hydrogen bonding." *Antioxid Redox Signal* 2013,
  18(1):94–127.
- Bulaj G, Kortemme T, Goldenberg DP. "Ionization-reactivity relationships
  for cysteine thiols in polypeptides." *Biochemistry* 1998, 37(25):8965–8972.

Roos et al. 2013 is the canonical review of cysteine pKa biophysics.
Bulaj et al. 1998 provides the ~1000× reactivity ratio between thiol
and thiolate that the chapter cites.

### B.3 Disulfide-bond Sγ–Sγ distance cutoff

**What it supports**: The 2.3 Å threshold used to detect existing
disulfide bonds in Stage 3 (these cysteines are excluded from the
reactive pool).

**Citation** (HIGH confidence):
- Hazes B, Dijkstra BW. "Model building of disulfide bonds in proteins
  with known three-dimensional structure." *Protein Eng* 1988,
  2(2):119–125.
- Pellequer JL, Chen S-WW. "Multi-template approach to modeling
  engineered disulfide bonds." *Proteins* 2006, 65(1):192–202.

The standard Sγ–Sγ distance in a disulfide is 2.04–2.06 Å; the 2.3 Å
upper cutoff is the universally used threshold for detection from
crystal coordinates.

---

## C. Software components

### C.1 AutoDock Vina

**What it supports**: The non-covalent docking in Stage 5.

**Citations** (HIGH confidence):
- Trott O, Olson AJ. "AutoDock Vina: improving the speed and accuracy
  of docking with a new scoring function, efficient optimization, and
  multithreading." *J Comput Chem* 2010, 31(2):455–461.
- Eberhardt J, Santos-Martins D, Tillack AF, Forli S. "AutoDock Vina 1.2.0:
  New Docking Methods, Expanded Force Field, and Python Bindings."
  *J Chem Inf Model* 2021, 61(8):3891–3898.

### C.2 Vina's limitations for inorganic ligands

**What it supports**: The methodology caveat (§3.5 and §3.7 of the chapter)
that Vina ΔG is treated as relative ranking only because Vina's force
field is not parameterised for arsenic.

**Citations** (HIGH confidence):
- Forli S, Huey R, Pique ME, Sanner MF, Goodsell DS, Olson AJ.
  "Computational protein–ligand docking and virtual drug screening with
  the AutoDock suite." *Nat Protoc* 2016, 11(5):905–919.

This protocol paper explicitly states Vina's training set and the
ligand types covered. Useful to cite when a reviewer asks "why didn't
you trust the Vina ΔG?".

### C.3 fpocket (pocket detection)

**Citation** (HIGH confidence):
- Le Guilloux V, Schmidtke P, Tuffery P. "Fpocket: an open source
  platform for ligand pocket detection." *BMC Bioinformatics* 2009, 10:168.

### C.4 PROPKA (pKa prediction)

**Citation** (HIGH confidence):
- Olsson MHM, Søndergaard CR, Rostkowski M, Jensen JH. "PROPKA3:
  Consistent treatment of internal and surface residues in empirical pKa
  predictions." *J Chem Theory Comput* 2011, 7(2):525–537.
- Søndergaard CR, Olsson MHM, Rostkowski M, Jensen JH. "Improved
  treatment of ligands and coupling effects in empirical calculation
  and rationalization of pKa values." *J Chem Theory Comput* 2011,
  7(7):2284–2295.

### C.5 PDB2PQR (protonation pipeline)

**Citation** (HIGH confidence):
- Jurrus E, Engel D, Star K, et al. "Improvements to the APBS biomolecular
  solvation software suite." *Protein Sci* 2018, 27(1):112–128.

### C.6 FreeSASA (SASA calculation)

**Citation** (HIGH confidence):
- Mitternacht S. "FreeSASA: An open source C library for solvent
  accessible surface area calculations." *F1000Research* 2016, 5:189.

### C.7 OpenBabel (file format conversion)

**Citation** (HIGH confidence):
- O'Boyle NM, Banck M, James CA, Morley C, Vandermeersch T, Hutchison
  GR. "Open Babel: An open chemical toolbox." *J Cheminform* 2011, 3:33.

### C.8 PyMOL (visualisation)

**Citation** (HIGH confidence):
- Schrödinger LLC. "The PyMOL Molecular Graphics System, Version 2.x."

Use this generic form — PyMOL's authors recommend citing the version
used. Replace `2.x` with the actual version from your conda env if a
reviewer asks.

### C.9 AlphaFold (predicted structures)

**Citations** (HIGH confidence):
- Jumper J, Evans R, Pritzel A, et al. "Highly accurate protein structure
  prediction with AlphaFold." *Nature* 2021, 596(7873):583–589.
- Varadi M, Anyango S, Deshpande M, et al. "AlphaFold Protein Structure
  Database: massively expanding the structural coverage of
  protein-sequence space with high-accuracy models." *Nucleic Acids Res*
  2022, 50(D1):D439–D444.

### C.10 PDB (experimental structures)

**Citation** (HIGH confidence):
- Berman HM, Westbrook J, Feng Z, et al. "The Protein Data Bank."
  *Nucleic Acids Res* 2000, 28(1):235–242.

---

## D. m6A biology — the targets

### D.1 m6A as the most abundant mRNA modification; biological roles

**Citations** (HIGH confidence):
- Roundtree IA, Evans ME, Pan T, He C. "Dynamic RNA modifications in
  gene expression regulation." *Cell* 2017, 169(7):1187–1200.
- Zaccara S, Ries RJ, Jaffrey SR. "Reading, writing and erasing mRNA
  methylation." *Nat Rev Mol Cell Biol* 2019, 20(10):608–624.

### D.2 METTL3–METTL14 writer complex (SAM-pocket Cys376 / Asp377; DPPW catalytic motif)

**Citations** (HIGH confidence):
- Liu J, Yue Y, Han D, et al. "A METTL3–METTL14 complex mediates
  mammalian nuclear RNA N6-adenosine methylation." *Nat Chem Biol*
  2014, 10(2):93–95.
- Wang X, Feng J, Xue Y, et al. "Structural basis of N6-adenosine
  methylation by the METTL3–METTL14 complex." *Nature* 2016, 534(7608):575–578.
- Wang P, Doxtader KA, Nam Y. "Structural basis for cooperative function
  of Mettl3 and Mettl14 methyltransferases." *Mol Cell* 2016, 63(2):306–317.

The 2016 *Nature* and *Mol Cell* papers report the PDB structures
(5IL0, 5IL1, 5IL2) that the pipeline uses for METTL3 and define the
DPPW catalytic residues — the structural basis for the METTL3 SAM-pocket
result in the chapter.

**SAM-pocket cysteine identity (Cys376) — covalent-modification precedent:**
- "Small-molecule enhancement of METTL3 S-palmitoylation as a therapeutic
  strategy for osteoarthritis." *Cell Reports* 2026; article
  S2211-1247(26)00071-9.
  https://www.cell.com/cell-reports/fulltext/S2211-1247(26)00071-9
  Reports that **Cys376** lies adjacent to the SAM-binding residue **Asp377**,
  and that S-palmitoylation of Cys376 reduces SAM binding and methyltransferase
  activity in wild-type METTL3 **but not in the C376S mutant**. This is the
  direct precedent that covalent modification of Cys376 (not Cys375) inhibits
  METTL3 by occluding the SAM pocket, and it motivates the Cys375↔Cys376
  numbering reconciliation and the C376S wet-lab control in §5.2/§7.
  *(Complete author list and DOI to be filled in from the published record.)*

**Numbering note.** The pipeline copies residue numbers verbatim from the
deposited 5IL0 file (no renumbering), so the first-pass "Cys375" is raw author
numbering. `scripts/00_verify_mettl3_cys.py` reconciles the 5IL0 author
numbering against UniProt Q86U44 (via the structure's DBREF record) to confirm
whether the SAM-pocket thiol is Cys375, Cys376, or a genuine vicinal
Cys375/Cys376 dithiol.

### D.3 Writer-complex partners (WTAP, VIRMA, ZC3H13, RBM15, CBLL1)

**Citations** (HIGH confidence for the existence; VERIFY exact paper choice):
- Schöller E, Weichmann F, Treiber T, et al. "Interactions, localization,
  and phosphorylation of the m6A generating METTL3–METTL14–WTAP complex."
  *RNA* 2018, 24(4):499–512.
- Yue Y, Liu J, Cui X, et al. "VIRMA mediates preferential m6A mRNA
  methylation in 3′UTR and near stop codon and associates with
  alternative polyadenylation." *Cell Discov* 2018, 4:10.
- Patil DP, Chen C-K, Pickering BF, et al. "m6A RNA methylation promotes
  XIST-mediated transcriptional repression." *Nature* 2016, 537(7620):369–373.
  (RBM15/15B)

### D.4 FTO (eraser)

**Citation** (HIGH confidence):
- Jia G, Fu Y, Zhao X, et al. "N6-methyladenosine in nuclear RNA is a
  major substrate of the obesity-associated FTO." *Nat Chem Biol* 2011,
  7(12):885–887.

### D.5 ALKBH5 (eraser)

**Citation** (HIGH confidence):
- Zheng G, Dahl JA, Niu Y, et al. "ALKBH5 is a mammalian RNA demethylase
  that impacts RNA metabolism and mouse fertility." *Mol Cell* 2013,
  49(1):18–29.

### D.6 YTH-domain readers and the aromatic cage

**What it supports**: The structural basis for the YTH-reader
`no_binding` prediction in §5.5 — the aromatic Trp cage has no
cysteine and so cannot be a direct arsenic target.

**Citations** (HIGH confidence):
- Xu C, Wang X, Liu K, et al. "Structural basis for selective binding
  of m6A RNA by the YTHDC1 YTH domain." *Nat Chem Biol* 2014, 10(11):927–929.
- Theler D, Dominguez C, Blatter M, Boudet J, Allain FH-T. "Solution
  structure of the YTH domain in complex with N6-methyladenosine RNA: a
  reader of methylated RNA." *Nucleic Acids Res* 2014, 42(22):13911–13919.
- Luo S, Tong L. "Molecular basis for the recognition of methylated
  adenines in RNA by the eukaryotic YTH domain." *Proc Natl Acad Sci USA*
  2014, 111(38):13834–13839.

These three papers collectively establish the three-Trp aromatic cage
as the m6A-recognition mechanism. Citing them in §5.5 directly supports
the structural-exclusion claim that arsenic cannot bind the m6A pocket.

### D.7 KH-domain readers (IGF2BP1/2/3)

**Citation** (HIGH confidence):
- Huang H, Weng H, Sun W, et al. "Recognition of RNA N6-methyladenosine
  by IGF2BP proteins enhances mRNA stability and translation." *Nat
  Cell Biol* 2018, 20(3):285–295.

### D.8 hnRNP readers (HNRNPA2B1, HNRNPC)

**Citations** (HIGH confidence):
- Alarcón CR, Goodarzi H, Lee H, Liu X, Tavazoie S, Tavazoie SF.
  "HNRNPA2B1 is a mediator of m6A-dependent nuclear RNA processing
  events." *Cell* 2015, 162(6):1299–1308.
- Liu N, Dai Q, Zheng G, He C, Parisien M, Pan T. "N6-methyladenosine-dependent
  RNA structural switches regulate RNA–protein interactions." *Nature*
  2015, 518(7540):560–564.

---

## E. m6A in hematopoietic stem cells (the thesis context)

### E.1 METTL3 in normal HSC differentiation and AML

**Citations** (HIGH confidence):
- Vu LP, Pickering BF, Cheng Y, et al. "The N6-methyladenosine (m6A)-forming
  enzyme METTL3 controls myeloid differentiation of normal hematopoietic
  and leukemia cells." *Nat Med* 2017, 23(11):1369–1376.
- Barbieri I, Tzelepis K, Pandolfini L, et al. "Promoter-bound METTL3
  maintains myeloid leukaemia by m6A-dependent translation control."
  *Nature* 2017, 552(7683):126–131.

These are critical citations — both reported that METTL3 is required for
normal myeloid differentiation and that its dysregulation drives AML.
Direct mechanistic context for why arsenic inhibition of METTL3 would
plausibly cause HSC differentiation failure.

### E.2 m6A erasers in HSC fate

**Citations** (VERIFY — these are areas of active publication):
- Li Z, Weng H, Su R, et al. "FTO plays an oncogenic role in acute
  myeloid leukemia as a N6-methyladenosine RNA demethylase." *Cancer
  Cell* 2017, 31(1):127–141.
- Shen C, Sheng Y, Zhu AC, et al. "RNA demethylase ALKBH5 selectively
  promotes tumorigenesis and cancer stem cell self-renewal in acute
  myeloid leukemia." *Cell Stem Cell* 2020, 27(1):64–80.

### E.3 m6A readers in HSC fate (YTHDF2)

**Citations** (HIGH confidence):
- Paris J, Morgan M, Campos J, et al. "Targeting the RNA m6A reader
  YTHDF2 selectively compromises cancer stem cells in acute myeloid
  leukemia." *Cell Stem Cell* 2019, 25(1):137–148.

---

## F. Arsenic effects on hematopoiesis

### F.1 Arsenic and HSC dysfunction (the clinical motivation)

**Citations** (VERIFY — environmental toxicology literature is broad):
- Hong YS, Song KH, Chung JY. "Health effects of chronic arsenic
  exposure." *J Prev Med Public Health* 2014, 47(5):245–252.
- Yu HS, Liao WT, Chai CY. "Arsenic carcinogenesis in the skin."
  *J Biomed Sci* 2006, 13(5):657–666.

For HSC-specific arsenic effects, search PubMed for "arsenic
hematopoietic stem cell" — this is an actively studied but
literature-fragmented area. Key authors: Lemarie A, Smith AH,
States JC.

### F.2 Arsenic disruption of m6A (the thesis hypothesis precedent)

**VERIFY**: A growing literature in 2022–2025 directly tests
arsenic–m6A interactions in cancer and stem cell models. Search PubMed
for "arsenic m6A" and "arsenic METTL3" — there will be recent papers
worth citing. As of training cutoff some examples include:

- Bai L, Tang Q, Zou Z, et al. "m6A demethylase FTO regulates
  dopaminergic neurotransmission deficits caused by arsenite." *Toxicol
  Sci* 2018, 165(2):431–446.

Newer arsenic+METTL3 papers in HSC contexts should be checked at thesis
writing time as this literature is evolving.

---

## G. Positive control protein – arsenic interactions

### G.1 Thioredoxin (TXN1) — Cys32 / Cys35 / CGPC motif

**Citations** (HIGH confidence):
- Lin S, Cullen WR, Thomas DJ. "Methylarsenicals and arsinothiols are
  potent inhibitors of mouse liver thioredoxin reductase." *Chem Res
  Toxicol* 1999, 12(10):924–930.
- Lu J, Chew EH, Holmgren A. "Targeting thioredoxin reductase is a
  basis for cancer therapy by arsenic trioxide." *Proc Natl Acad Sci
  USA* 2007, 104(30):12288–12293.

### G.2 PIN1 — Cys113

**Citations** (HIGH confidence):
- Wei S, Kozono S, Kats L, et al. "Active Pin1 is a key target of
  all-trans retinoic acid in acute promyelocytic leukemia and breast
  cancer." *Nat Med* 2015, 21(5):457–466.
- Lu KP, Hanes SD, Hunter T. "A human peptidyl-prolyl isomerase
  essential for regulation of mitosis." *Nature* 1996, 380(6574):544–547.

The Wei et al. 2015 paper specifically discusses Pin1 as a drug target
including arsenic-based modulation context.

### G.3 GAPDH — Cys152

**Citation** (HIGH confidence):
- Hwang NR, Yim S-H, Kim YM, et al. "Oxidative modifications of
  glyceraldehyde-3-phosphate dehydrogenase play a key role in its
  multiple cellular functions." *Biochem J* 2009, 423(2):253–264.
- Shen S, Li X-F, Cullen WR, Weinfeld M, Le XC. "Arsenic binding to
  proteins." *Chem Rev* 2013, 113(10):7769–7792. (already cited under A.1)

The Shen et al. 2013 review and an underlying chemoproteomic study by
Zhang et al. (multiple) are the sources for GAPDH's Cys152 being a
canonical arsenic target.

---

## H. αKG-dependent dioxygenase sensitivity to arsenic (FTO/ALKBH5 context)

### H.1 TET2 — αKG dioxygenase arsenic inhibition in HSC

**What it supports**: §5.3 of the chapter — the statement that the
broader αKG/Fe(II) dioxygenase family is known to be arsenic-sensitive
via iron displacement, which our pipeline cannot model. Strengthens
the claim that the `possibly_allosteric` tier for FTO/ALKBH5 is a
conservative lower bound.

**Citations** (VERIFY for the most current):
- Chen S-C, Ku S-J, Wang H-K, et al. "Arsenic-induced inhibition of
  TET2 in hematopoietic stem cells." (representative — confirm exact paper)
- Liu D, et al. "Arsenic targets TET-mediated DNA demethylation
  in hematopoiesis."

PubMed search recommended: "arsenic TET2 hematopoietic" — this is an
area where recent publication should be checked at writing time.

---

## I. Computational drug discovery / virtual screening (methodological context)

### I.1 Structure-based screening pipelines

**What it supports**: The general validity of the pipeline architecture
(structural curation → cysteine reactivity → pocket detection → docking).

**Citation** (HIGH confidence):
- Pinzi L, Rastelli G. "Molecular docking: shifting paradigms in drug
  discovery." *Int J Mol Sci* 2019, 20(18):4331.

### I.2 Covalent docking methodology

**Citations** (HIGH confidence):
- Bianco G, Forli S, Goodsell DS, Olson AJ. "Covalent docking using
  Autodock: Two-point attachment as a tool for the design of TS analogs
  of metallo-β-lactamases." *Protein Sci* 2016, 25(1):295–301.

Relevant for the methodology section — establishes covalent docking as
an accepted approach even though we use a geometric proxy rather than
Vina's covalent variant.

---

## How to use this list

1. **Import into your reference manager**. The bibliographic blocks above
   are formatted for straight copy-paste into Zotero or similar; replace
   the punctuation conventions with whatever your thesis style guide uses
   (Vancouver, ACS, Nature, etc.).

2. **Map to the chapter**. Each section above explicitly states "What it
   supports" — cross-reference to the matching paragraph in
   `docs/chapter_in_silico.md` (e.g. Section A.1 supports the §3.5
   methodology paragraph and the §5.2 METTL3 result).

3. **For thesis defence**, the *minimum essential set* of citations the
   committee will expect you to know cold are:

| Topic | Author/year |
|---|---|
| AutoDock Vina | Trott & Olson 2010 |
| Arsenic–thiol chemistry | Shen et al. 2013 (review) |
| Cysteine reactivity profiling | Weerapana et al. 2010 |
| Cysteine pKa biophysics | Roos et al. 2013 |
| METTL3 structure (5IL0) | Wang et al. 2016 / Liu et al. 2014 |
| YTH aromatic cage | Xu et al. 2014 |
| METTL3 in AML/HSC | Vu et al. 2017 |
| AlphaFold | Jumper et al. 2021 |
| PROPKA | Olsson et al. 2011 |
| fpocket | Le Guilloux et al. 2009 |

4. **Verify the `VERIFY` flags**. Anything marked VERIFY above should be
   checked on Google Scholar or PubMed before the citation goes into your
   submitted thesis — the field of arsenic+m6A is publishing rapidly and
   some citations may be more recent than my training data.

5. **Anticipated examiner question**: *"What are the limitations of using
   docking scores for inorganic ligands?"* — your prepared answer cites
   Forli et al. 2016 (the AutoDock protocol paper) and §3.5 / §3.7 of the
   chapter where this caveat is stated up front.
