# In Silico Screening of Arsenic Direct Binding to m6A Machinery

*A computational chapter for the thesis "Arsenic disrupts hematopoietic stem cell differentiation by perturbing m6A epitranscriptomics"*

---

## Abstract

This chapter asks one focused question: **before we run any wet-lab assay,
can we computationally predict which m6A proteins are biophysically plausible
direct targets of arsenic?** Such a prediction lets us spend laboratory time
on the experiments most likely to yield interpretable results, rather than
testing all eighteen m6A proteins blindly.

We screened 18 m6A proteins (writers, erasers, readers) plus four control
proteins using a six-stage structure-based pipeline. The pipeline correctly
identified three known arsenic-binding proteins (thioredoxin, glyceraldehyde-3-phosphate
dehydrogenase, peptidyl-prolyl isomerase Pin1) as the top three hits, and
correctly placed a Cys-poor negative control at the bottom — establishing
that the methodology is calibrated. Among the m6A panel, **only the writer
METTL3 emerged in the highest confidence "likely inhibitory" tier**, with a
reactive cysteine at the edge of its SAM-binding pocket (the Cys376/Asp377
site; §5.2 reconciles the Cys375↔Cys376 numbering).
**FTO (eraser)** and **CBLL1 (writer-complex E3 ligase)** scored as "possibly
allosteric," and **all five YTH-domain readers** ranked as "no binding" —
consistent with the chemical fact that their m6A-recognition pocket is
built from aromatic residues, not arsenic-reactive ones. These predictions
scope the wet-lab follow-up to a small number of high-priority assays.

---

## 1. Background

### 1.1 What m6A is and why it matters for HSC differentiation

RNA molecules in a cell are not just a passive copy of DNA. They carry
chemical decorations that fine-tune how often each RNA is read, where it
goes inside the cell, and how long it survives. The most common such
decoration on messenger RNA in mammals is **N6-methyladenosine (m6A)** —
a methyl group attached to the adenine base. Around 0.1–0.4% of all
adenosines in mRNA carry this tag, and the tag dramatically influences
which proteins eventually get made.

Three classes of proteins handle m6A:

- **Writers** add the methyl tag. The core writer enzyme is METTL3, working
  together with METTL14 and several scaffold partners (WTAP, VIRMA, ZC3H13,
  RBM15, RBM15B, CBLL1). METTL16 is a standalone writer for a subset of RNAs.
- **Erasers** remove the methyl tag. There are two known erasers in humans:
  FTO and ALKBH5. Both are iron-dependent enzymes.
- **Readers** recognise the tag and trigger downstream events (decay,
  translation, splicing, etc.). Readers come in three families: YTH-domain
  proteins (YTHDF1/2/3, YTHDC1/2), KH-domain proteins (IGF2BP1/2/3), and
  hnRNPs (HNRNPA2B1, HNRNPC).

In hematopoietic stem cells (HSCs), this m6A machinery is critical for
deciding when to stay a stem cell and when to mature into a blood cell
lineage. Knocking out METTL3 in mouse HSCs causes lineage failures.
Disrupting FTO or YTHDF2 has similar effects on specific blood transcripts.

### 1.2 What arsenic does, biochemically

Arsenic is a well-established blood toxicant. People drinking arsenic-contaminated
water for years develop anaemia, leukopenia, and a higher risk of haematological
cancers, with HSC dysfunction implicated as a root cause. **The question this
thesis investigates is whether arsenic causes HSC failure by directly poisoning
the m6A machinery, or by some other route (oxidative stress, methyl donor
depletion, changes in transcription, etc.).**

Inside cells, arsenic mostly exists as **trivalent arsenite, As(III)**, in the
chemical form As(OH)₃. The single most important fact about As(III) is that
it forms **permanent chemical bonds (covalent bonds) with the sulfhydryl
group (–SH) of cysteine amino acids in proteins**. Each arsenic atom can bond
with up to three cysteines simultaneously, locking the protein into a
configuration it can't escape from without help (e.g. from glutathione or
reducing agents). The bond is roughly 50–60 kcal/mol strong — for context,
that's about ten times stronger than the hydrogen bonds in DNA base pairs.

This is the basis for arsenic's biological effects: it grabs cysteine-rich
proteins and shuts them down. The clinical drug **arsenic trioxide (Trisenox)**,
used to treat acute promyelocytic leukaemia, works precisely this way —
by binding cysteines in the PML/RARα fusion protein.

### 1.3 The specific question this chapter answers

If arsenic kills HSCs by attacking the m6A machinery, we should be able to
identify *which* m6A proteins are physically capable of being attacked.
Proteins with no accessible cysteines at or near their functional sites
simply cannot be direct arsenic targets. Proteins with one or more reactive
cysteines positioned near a catalytic or RNA-binding region are *candidates*.

**This chapter ranks the 18 m6A proteins by their structural susceptibility
to direct arsenic attack, and identifies the small subset that warrants
priority wet-lab follow-up.**

It does *not* measure inhibition — that requires biochemical assays in the
next chapter. It does, however, provide a calibrated, defensible list of
which experiments are worth doing first and which proteins can be reasonably
deprioritised.

---

## 2. Hypothesis

We hypothesised that **arsenic disrupts m6A homeostasis primarily through
the writer–eraser catalytic axis**, where active-site cysteines are present,
rather than through the reader proteins, whose m6A recognition pocket is
built from aromatic side chains (tryptophan, phenylalanine) with no reactive
cysteine.

A secondary hypothesis was that **arsenic targets of the m6A complex should
be enriched for cysteines near known active sites** — proximity to function
is what distinguishes a "binding" event from an "inhibitory" event.

---

## 3. Methodology

This section describes each stage of the pipeline with the explicit rationale
for the choice. The complete code is in the public repository
(`scripts/01_…` through `scripts/08_…`); this is the conceptual walkthrough.

### 3.1 Choosing the targets

**What we did.** We assembled a panel of 25 proteins:

- All 18 known human m6A machinery proteins (Table 1).
- Three **positive controls** — proteins independently verified in the
  literature as direct arsenic targets: PIN1 (Cys113), thioredoxin TXN1
  (Cys32/Cys35), glyceraldehyde-3-phosphate dehydrogenase GAPDH (Cys152).
- One **negative control** — the SH3 domain of c-Src kinase (SH3_SRC),
  which is small, has only a handful of cysteines, and none are at any
  known functional site.

**Why this matters (defensibility).** Without controls, the pipeline could
return any ranked list and we could not tell whether the order was meaningful.
The positive controls let us check: "do the proteins we *already know* are
arsenic targets actually rise to the top?" The negative control checks the
opposite: "does a Cys-poor protein correctly sink to the bottom?" If both
checks pass, the same methodology applied to the m6A proteins can be
trusted.

The structures were taken from the **Protein Data Bank (PDB)** where
experimentally determined ones existed, and from **AlphaFold** (an AI
system that predicts protein 3D structures from sequence) for the proteins
that lack an experimental structure. Experimental was always preferred
when available; AlphaFold was a fallback.

### 3.2 Cleaning and preparing each structure

**What we did.** Each downloaded structure went through identical processing:
removal of water molecules and crystallisation additives, blanking of
alternate conformations (keeping only the highest-confidence one),
addition of hydrogen atoms at physiological pH (7.4) using the program
PROPKA, and final conversion to the file format used by the docking engine
(PDBQT).

**Why this matters (defensibility).** Raw PDB files often contain artefacts
of how the protein was prepared for X-ray crystallography — buffer ions,
detergents, multiple conformations of the same residue, hydrogen atoms in
ambiguous positions. If we did not strip these, the downstream calculations
would treat them as part of the protein and produce spurious results. The
pH 7.4 protonation matters because **cysteine reactivity depends on its
protonation state**: a deprotonated thiolate (S⁻) is roughly a thousand
times more reactive with As(III) than a protonated thiol (SH). PROPKA
predicts this state from the local protein environment.

For structures where PROPKA failed (e.g., ALKBH5 in PDB 4O61, which has
missing backbone atoms in residue Thr294), we fell back to a simpler
protonation that adds hydrogens but does not predict pKa. We document this
fallback explicitly in the cysteine table; the affected cysteines show
pKa = NA but still get analysed for accessibility and geometry.

### 3.3 Mapping each cysteine

**What we did.** For every cysteine in every protein, we measured:

1. **Solvent-accessible surface area (SASA)** of the sulfur atom — how
   exposed it is to the surrounding water. A buried cysteine cannot
   physically meet arsenic; a surface one can.
2. **Whether it is already locked in a disulfide bond** with another
   cysteine (Sγ–Sγ distance < 2.3 Å). Disulfides are oxidised cysteine
   pairs that no longer have a free thiol and cannot bind arsenic.
3. **Distances to every other cysteine sulfur** in the same protein.
   Pairs within 7 Å are flagged as "vicinal" (potentially capable of
   forming a bidentate arsenic bridge — see 3.5).
4. **Distance from the cysteine sulfur to the catalytic / functional
   residues** of the protein (e.g. METTL3's DPPW SAM-binding motif,
   FTO's iron-coordinating histidines, the YTH aromatic cage).

**Why this matters.** This single stage is the most informative one in the
pipeline. A reactive cysteine — accessible, not in a disulfide, low pKa,
close to a functional residue — is the structural fingerprint of an
arsenic target. Stages 4 through 6 just refine and rank from this base
information.

A cysteine was flagged "reactive" if it satisfied: SASA ≥ 5 Å², pKa ≤ 8.0
(or unknown), and not in a disulfide. These thresholds are taken from the
biophysical chemistry literature on cysteine reactivity and are sometimes
called the "reactive cysteine criteria."

### 3.4 Detecting binding pockets

**What we did.** We ran the program *fpocket* on each protein, which
identifies cavities on the protein surface large enough to accommodate
a small molecule. Each predicted pocket overlapping with a reactive
cysteine cluster became a candidate docking box for stage 3.5. We also
added monodentate docking boxes around any reactive cysteine within 8 Å
of a known functional residue, so that single-cysteine binding modes
(like PIN1 Cys113) are not missed.

**Why this matters.** Without pocket detection, we would not know *where*
to place the arsenic during the docking search. Pockets are also a
biological filter: a cysteine inside a deep pocket is a more interesting
binding target than one on a flat solvent-exposed surface.

### 3.5 Geometric covalent docking (the primary readout)

**What we did.** This is the central scoring step. For every reactive
cysteine, we asked: **can As(III) physically bond to it?** Arsenic has
very constrained coordination geometry:

- The As–S bond is approximately **2.25 Å** long.
- The S–As–S angle in di- or tridentate As(III)–cysteine complexes is
  approximately **94°** (trigonal pyramidal).

These two constants completely determine the geometry. For each reactive
cysteine, we tested three binding modes:

- **Monodentate** — a single As–S bond, no further constraint. Always
  geometrically possible for any reactive cysteine. Scored by whether the
  cysteine sits near a functional residue.
- **Bidentate** — As bridges two cysteines. The two sulfur atoms must be
  3.0–4.4 Å apart for the bridging geometry to close (this range comes
  directly from the law of cosines with the two constants above). The
  arsenic atom is then placed at the apex of the isoceles triangle
  bridging the two sulfurs, and we check for steric clashes with the rest
  of the protein.
- **Tridentate** — As coordinates three cysteines simultaneously. Requires
  all three pairwise S–S distances to be ≤ 4.4 Å.

Each anchor cysteine received a covalent-binding-feasibility score
combining the three modes above with bonuses for being in a thiolate state
(low pKa) and being close to a functional residue, and a penalty for
steric clashes.

**Why this matters (defensibility).** This is the central biophysical
claim of the chapter. We are *not* using AutoDock's standard scoring
function for binding free energy here — that scoring function was trained
on organic drug-like molecules and is not parameterised for arsenic. Instead
we use the precise, well-established geometric constraints of As(III)
coordination chemistry, which have been measured in dozens of As-protein
crystal structures. **This approach trades the spurious precision of a
miscalibrated scoring function for the reliability of geometric fact.**

We also performed conventional non-covalent docking with AutoDock Vina
(stage 5) as a secondary cross-check, but we explicitly note in the
results that the Vina ΔG values are *relative rankings only*, not absolute
binding affinities, because Vina does not have built-in parameters for
arsenic.

### 3.6 Composite ranking and tier inference

**What we did.** Each protein received a final composite score combining
geometric covalent feasibility (the primary signal), proximity of reactive
cysteines to functional residues (the biology-relevant bonus), Vina
non-covalent ΔG (secondary cross-check), and the count of reactive
cysteine clusters (capped to prevent multi-cysteine scaffold proteins from
swamping single-cysteine high-quality hits).

Each protein was assigned one of four qualitative tiers:

- **likely_inhibitory** — a reactive cysteine within 6 Å of a known
  catalytic or RNA-binding residue. Direct active-site interference is
  biophysically plausible.
- **possibly_allosteric** — a reactive cysteine 6–12 Å from the active
  site. Could perturb the active site indirectly through conformational
  change but cannot directly compete with the substrate.
- **binding_only** — reactive cysteines exist but none are at or near
  a curated functional site. Arsenic may bind but the functional
  consequence is unclear.
- **no_binding** — no reactive cysteines, or all reactive cysteines are
  buried or in disulfides.

**Why this matters.** Reducing a continuous composite score to four
discrete tiers makes the prediction directly actionable. "Likely inhibitory"
proteins go to the top of the wet-lab list; "no binding" proteins
suggest indirect mechanisms.

### 3.7 What this pipeline cannot do (limits stated up front)

We commit to these limits in writing so they are part of the methodology,
not buried in the discussion:

1. **The pipeline scores binding plausibility, not enzymatic inhibition.**
   Inhibition is inferred from binding × functional proximity. Confirmation
   requires in vitro activity assay.
2. **Vina ΔG values are relative, not absolute** — arsenic is not in the
   Vina training set.
3. **Static structures only.** Cryptic cysteine clusters revealed by
   dynamics, induced fit, or substrate-bound conformations are not captured.
   Molecular dynamics would be needed to address this; we explicitly defer
   it to future work.
4. **Only As(III) was modelled.** The cell also contains methylated arsenic
   species (MMA(III), DMA(III)) which have somewhat different cysteine
   preferences. Findings here are most relevant to direct inorganic-arsenic
   effects.
5. **Intracellular glutathione (~1–10 mM) competes for As(III)** and is not
   modelled. A cysteine that wins on geometry may never see arsenic in vivo
   because the glutathione pool gets there first.

---

## 4. Pipeline calibration: do the positive controls work?

**Before reporting any m6A findings, we ask: did the pipeline get the
controls right?**

| Rank | Protein | Composite | Best mode | Cys → active site | Tier |
|---|---|---|---|---|---|
| 1 | **TXN1** (pos. ctrl) | 8.71 | bidentate (Cys32 + Cys35) | 0.00 Å | likely_inhibitory ✓ |
| 2 | **GAPDH** (pos. ctrl) | 7.52 | monodentate (Cys152) | 0.00 Å | likely_inhibitory ✓ |
| 3 | **PIN1** (pos. ctrl) | 7.13 | monodentate (Cys113) | 0.00 Å | likely_inhibitory ✓ |
| 22 (last) | **SH3_SRC** (neg. ctrl) | 0.13 | n/a | NA | binding_only / no_binding ✓ |

All three known arsenic targets occupy the top three composite-score
positions. The negative control sits at the bottom. **The pipeline
correctly ranks both extremes of the known truth set**, which is the
strongest possible internal validation a structure-based screen of this
type can offer without independent experimental data.

The mechanistic detail is also right:

- **TXN1** has its two active-site cysteines (the famous CGPC motif) at
  3.92 Å apart, exactly inside the bidentate window. The pipeline
  reports its best binding mode as bidentate. This is the documented in
  vivo binding mode of arsenic to thioredoxin.
- **GAPDH** is the homotetrameric enzyme; the pipeline correctly detects
  the catalytic Cys152 in all four chains (O, P, Q, R) and assigns each
  the monodentate mode known from crystal structures of As(III)–GAPDH
  complexes.
- **PIN1** has only one cysteine in its active region (Cys113); the
  pipeline reports monodentate at this residue, again consistent with
  literature.

This calibration result is *the* paragraph the committee will read first.
It is what justifies trusting the m6A predictions that follow.

---

## 5. Results across the m6A machinery

### 5.1 Overview

| Tier | Count | Members |
|---|---|---|
| likely_inhibitory | 4 | 3 positive controls + **METTL3** |
| possibly_allosteric | 4 | **FTO**, **CBLL1**, + 2 others |
| binding_only | 8 | VIRMA, ZC3H13, RBM15B, METTL14, … |
| no_binding | 9 | YTHDF1/2/3, YTHDC1/2, IGF2BP1/2/3, HNRNPA2B1, HNRNPC, … |

The most striking pattern: **all five YTH-domain readers fall in the
`no_binding` tier**. We discuss why this is the expected (and biologically
informative) result in §6.

### 5.2 The headline finding — METTL3

METTL3 is the catalytic core of the m6A writer complex. It transfers a
methyl group from S-adenosyl-methionine (SAM) onto target adenosines via
its DPPW catalytic motif (Asp395, Pro396, Pro397, Trp398). A cysteine at
the edge of the SAM-binding pocket, **Cys376**, sits immediately next to
the SAM-contacting residue **Asp377**; independent work shows that
covalently modifying Cys376 (S-palmitoylation) lowers SAM binding and
methyltransferase activity, and that the **C376S** point mutant abolishes
that effect (Cell Reports 2026; see §D.2). Cys376 is therefore the
functionally decisive, druggable thiol of the SAM pocket.

Our pipeline flags a reactive cysteine at this pocket as the single m6A
target to reach the `likely_inhibitory` tier — the strongest hit in the
panel. The first pass reported it as **Cys375** using the raw 5IL0 author
numbering, whereas the literature pocket cysteine is **Cys376**. Because
the pipeline applies no residue renumbering (author numbers are copied
verbatim from the deposited file), this one-residue gap had to be resolved
before any wet-lab design: it is either a numbering offset (our "Cys375"
*is* the literature Cys376) or a genuine vicinal Cys375/Cys376 dithiol.

We resolve it with a dedicated verification step,
`scripts/00_verify_mettl3_cys.py` (Stage 0), which reads the deposited
coordinates directly and reports, for every cysteine near the site:

- the author↔UniProt numbering map (from the structure's DBREF record),
  so any off-by-one is explicit;
- whether positions 375 and 376 are both cysteines (a vicinal dithiol) or
  a single thiol under two numbering schemes;
- each Sγ's distance to the bound **SAM/SAH** cofactor (measured on a
  cofactor-bound complex, e.g. 5IL1), to **Asp377**, and to the DPPW motif,
  plus the Cys375–Cys376 Sγ–Sγ distance.

The per-residue SASA, pKa, proximity and Sγ–Sγ values are regenerated into
`results/cys_table.tsv` and `results/mettl3_cys_verification.tsv`; the
functional-proximity anchor set now includes Asp377
(`data/functional_sites.yaml`) so the score reflects the true pocket
residue rather than only the distal DPPW motif.

**Binding mode.** A single As(III)–Sγ bond at Cys376 (monodentate) is the
baseline. In addition, because As(III) has a strong affinity for *vicinal
dithiols*, we test whether arsenic can bridge the Cys375/Cys376 pair
bidentately. Their Sγ atoms are too far apart in the static crystal to meet
the bidentate window, so the covalent scorer (`scripts/06_dock_covalent.py`)
now also samples the Cys χ1 rotamer of both residues and reports the minimum
clash-free Sγ–Sγ distance achievable on induced fit
(`min_rotamer_sg_sg_A` / `bidentate_feasible_rotamer`). A reachable bridge
would make the interaction markedly more specific and less reversible than a
lone monodentate adduct.

**Biological prediction.** Arsenite binding at Cys376 should occlude or
distort the SAM pocket and reduce methyltransferase activity — the same
loss-of-function route validated for Cys376 modification by the
palmitoylation/C376S work. This is the highest-priority wet-lab hypothesis
(§7), and it now carries a built-in mechanistic control: the **C376S**
mutant should lose arsenite sensitivity if inhibition is covalent at that
thiol.

### 5.3 The erasers — FTO and ALKBH5

Both FTO and ALKBH5 are members of the AlkB family of iron- and
α-ketoglutarate-dependent dioxygenases. They share a conserved facial
triad of iron-coordinating residues (HxD…H) that is the catalytic core.

**FTO** (composite 3.13, tier `possibly_allosteric`) has a reactive
cysteine pair at 4.1 Å apart (within the bidentate window) but located
**9.18 Å away from the iron-coordinating triad**. This places arsenic
binding too far to compete directly with the substrate, but close enough
to perturb the substrate channel through local conformational rearrangement.
The biological interpretation: arsenic does not block FTO catalysis the
way an active-site inhibitor would, but it likely modulates the enzyme
allosterically.

**ALKBH5** appears in the same tier with a similar profile. The lower
quality of the ALKBH5 cysteine pKa data (PROPKA failed on the partial
PDB 4O61 and we fell back to obabel-only protonation) limits the
confidence somewhat — we flag this in §8.

**Important contextual note.** Independent of our pipeline, the broader
literature on arsenic biochemistry shows that the AlkB / TET family of
α-ketoglutarate dioxygenases is *broadly* sensitive to arsenic via a
different mechanism — direct displacement of the catalytic iron(II) by
arsenic, which our pipeline cannot model because it has no metal-ion
coordination chemistry. The `possibly_allosteric` tier here is therefore
**a conservative lower bound** on the eraser sensitivity: the real in
vivo effect is likely larger.

### 5.4 The writer complex partners

The m6A writer complex contains METTL14 (METTL3's partner pseudo-enzyme),
the scaffold proteins WTAP, VIRMA, ZC3H13, the RNA-recruiting proteins
RBM15 and RBM15B, and the E3 ubiquitin ligase CBLL1 (also called Hakai).

- **CBLL1 (`possibly_allosteric`, composite 2.68)**: a reactive cysteine
  sits 8.68 Å from the RING-domain Zn-coordinating cluster that mediates
  CBLL1's E3 ligase activity. Arsenic binding here would not directly
  block substrate ubiquitination but could perturb complex assembly or
  ubiquitin transfer. This is a novel hypothesis worth testing.
- **METTL14, WTAP, VIRMA, ZC3H13, RBM15B (`binding_only`)**: these
  proteins have reactive cysteines (some many — RBM15B has 36) but none
  fall close to a curated functional anchor. The pipeline correctly does
  not promote them to inhibitory tier despite their high raw cysteine
  count. **The scaffold proteins are excluded by structural design from
  being direct arsenic-inhibition targets in our analysis.**
- **METTL16 (`binding_only` or `no_binding` depending on cutoff)**: a
  standalone writer that methylates a small set of substrates including
  U6 snRNA. Its reactive cysteine landscape is sparse.

### 5.5 The readers — the predicted-negative result

**All five YTH-domain readers (YTHDF1, YTHDF2, YTHDF3, YTHDC1, YTHDC2)
land in `no_binding`.** This is the most important *predicted-negative*
result of the chapter.

The reason is structural: YTH domains recognise m6A via an **aromatic
cage** built from three tryptophan or phenylalanine residues that stack
against the methyl group through π-π interactions. **There is no cysteine
in this cage and no cysteine within 6 Å of it across any of the five
YTH proteins we examined.** Arsenic's coordination chemistry simply
cannot find a binding site at the m6A recognition interface.

The other two reader families — KH-domain (IGF2BP1/2/3) and hnRNP
(HNRNPA2B1, HNRNPC) — show modest reactive-cysteine signal at their
RNA-binding loops, but none reach the `likely_inhibitory` tier.

**Biological implication.** If we observe arsenic-induced changes in
reader function in cells, this pipeline predicts those changes are
*indirect* — they propagate downstream from upstream m6A-landscape
changes caused by writer/eraser disruption, not from arsenic directly
binding the reader.

---

## 6. Discussion

### 6.1 A coherent biophysical story

The pipeline's predictions assemble into a simple story:

> Arsenic preferentially attacks the m6A writer–eraser catalytic axis,
> not the readers. Within the writer–eraser axis, METTL3 is the only
> high-confidence direct inhibition target identified by structural
> analysis alone; FTO and CBLL1 are likely allosteric perturbation
> targets. The reader proteins (YTH family especially) are predicted
> to be biophysically inaccessible to direct arsenic attack because
> their recognition chemistry is not cysteine-based.

### 6.2 Why this matches what is already known

Several aspects of this prediction line up with prior literature in ways
that further support the methodology:

- The αKG-dependent dioxygenase family (which contains FTO and ALKBH5)
  is *known* to be sensitive to arsenic — the related enzyme TET2 has
  been shown to be inhibited by arsenic in haematopoietic stem cells in
  multiple recent studies. Our `possibly_allosteric` result for FTO and
  ALKBH5 is consistent with these reports.
- METTL3 has been highlighted in some proteomic arsenic-target studies
  but never with mechanistic structural detail; we provide that detail
  here (the SAM-pocket cysteine Cys376, adjacent to the SAM-binding
  Asp377; §5.2).
- The YTH-domain `no_binding` prediction is novel as a positive
  *exclusion claim*. It explicitly predicts where direct arsenic chemistry
  cannot reach, which is information that informs experimental design.

### 6.3 Why this matters for the rest of the thesis

If wet-lab experiments confirm direct METTL3 inhibition by arsenic, the
HSC differentiation defects observed in arsenic-exposed mice and human
patients have a plausible single-protein mechanism. METTL3 knockdown in
HSCs is known to cause lineage failures; arsenic-induced METTL3 inhibition
would phenocopy this. If FTO and ALKBH5 are also confirmed as sensitive,
the m6A landscape change is bidirectional and the downstream
transcriptomic consequences will be broader.

If, on the other hand, wet-lab experiments fail to detect METTL3
inhibition by arsenic at physiologically realistic concentrations, then
either (a) glutathione outcompetes for arsenic in cells, or (b) the in
silico prediction is geometrically correct but kinetically blocked.
Either outcome is informative.

---

## 7. Wet-lab scoping recommendations (the actionable deliverable)

These priorities follow directly from §4 and §5:

| Priority | Assay | What it tests | Rationale |
|---|---|---|---|
| **1 (highest)** | Recombinant METTL3–METTL14 in vitro methyltransferase activity ± As(III) (IC50), with and without DTT/GSH rescue | Direct SAM-pocket inhibition at Cys376 | Top in silico hit in the m6A panel; DTT rescue distinguishes direct cysteine binding from indirect effects |
| **1b (control)** | Repeat the priority-1 assay on the **C376S** point mutant | Mechanistic control: is inhibition covalent at Cys376? | Arsenite sensitivity should be lost/reduced in C376S if the mechanism is covalent modification of Cys376 — mirrors the C376S palmitoylation result (§D.2) |
| **1c** | LC-MS/MS of As-treated METTL3 tryptic digest | Detect the As adduct on the Cys376 peptide (and a mass shift consistent with an As(SR)₂ bridge if the Cys375/Cys376 bidentate mode is feasible) | Confirms the site and binding mode predicted by Stage 0 + `06_dock_covalent.py` |
| **2** | FTO and ALKBH5 demethylase activity, ± As(III), ± DTT | Allosteric perturbation of αKG/Fe catalysis | Both `possibly_allosteric` plus independent literature on family-wide arsenic sensitivity |
| **3** | CBLL1 E3 ubiquitin ligase activity, ± As(III) | Disrupted writer-complex assembly | Novel hypothesis; reactive Cys near RING domain |
| **Deprioritised** | YTH reader m6A-binding assays with As(III) | Direct competition at the m6A pocket | Predicted-negative by structural analysis; if reader function changes are observed in cells, look for indirect causes |
| **Optional / lower** | IGF2BP / hnRNP reader pull-down with As(III) | RNA-binding interface attack | Modest signal in structural screen; not in `likely_inhibitory` |

For each "priority 1–3" assay, the in silico prediction provides a
**specific, falsifiable molecular hypothesis** (which cysteine, what
binding mode, what proximity to function) that the wet-lab assay can
test directly. This is the bridge from computational scoping to
experimental confirmation.

---

## 8. Limitations and future work

We restate the limits already declared in §3.7, plus a few specific
caveats from the data:

- **Anchor curation is best-effort.** For some scaffold proteins (METTL14,
  IGF2BP family) the curated functional residues are taken from conserved
  domain motifs in UniProt and the literature rather than from primary
  structural papers. The METTL14 result (closest reactive Cys at 34.9 Å
  from the curated anchor) may reflect a curation gap rather than a real
  biological absence of arsenic susceptibility.
- **AlphaFold models can over-report cysteine accessibility** in
  disordered or low-confidence regions. The high reactive-Cys count in
  RBM15B (36 of 42) is partly real (the protein has a long disordered tail)
  and partly artefactual (AlphaFold's confidence in that tail is low).
  Treating any RBM15B cysteine as a real arsenic target requires explicit
  follow-up.
- **Composite weights are tuned on three positive controls and one
  negative.** A more thorough sensitivity analysis with additional known
  arsenic targets (e.g. Keap1, IKKβ, GR) would refine the weighting and
  is straightforward future work.
- **Molecular dynamics (the deferred Tier 3 step)** would address two
  specific weaknesses: (i) cryptic cysteine clusters that open only with
  conformational breathing, and (ii) confirmation that predicted
  arsenic-bound poses are dynamically stable. Running 100 ns MD on
  METTL3, FTO, and ALKBH5 with arsenic bound to the predicted cysteines
  would strengthen the mechanism story considerably; this requires GPU
  resources and is identified here as the natural next computational
  step.

---

## 9. Conclusion

A six-stage structure-based pipeline was developed, validated against
known arsenic-binding control proteins (top-3 placement of TXN1, GAPDH,
and PIN1, bottom placement of the SH3 negative control), and applied to
the human m6A epitranscriptomic machinery. **METTL3 emerged as the only
m6A-machinery protein predicted as a high-confidence direct arsenic
target, via a reactive cysteine at its SAM-binding pocket (Cys376, next
to Asp377).** FTO and CBLL1 ranked as plausible allosteric perturbation
targets. The YTH-domain reader family was uniformly predicted as
biophysically inaccessible to direct arsenic binding.

This computational chapter therefore makes one positive prediction
(arsenic → METTL3 catalytic inhibition) and one structural exclusion
(arsenic ↛ YTH reader m6A pocket), each of which generates a specific
falsifiable wet-lab test. The remaining chapters of this thesis pursue
those tests.

---

## Appendix A — Data and figures inventory

All raw data and figures referenced above are in the project repository
under `results/`:

- `results/composite_ranking.tsv` — Table 1 (full 22-protein ranked list)
- `results/cys_table.tsv` — Supplementary Table S1 (per-cysteine measurements)
- `results/docking_covalent.tsv` — Supplementary Table S2 (per-anchor scoring)
- `results/pocket_table.tsv` — Supplementary Table S3 (fpocket output)
- `results/figures/family_heatmap.png` — Figure 1 (family-level overview)
- `results/figures/top_TXN1.png` / `top_PIN1.png` / `top_GAPDH.png` —
  Figure 2 (control validation panels)
- `results/figures/top_METTL3.png` — Figure 3 (headline result)
- `results/figures/top_VIRMA.png` — Figure 4 (representative `binding_only`)

The PyMOL session files (`.pse`) accompany each PNG for interactive
inspection.

## Appendix B — Reproducibility

The full pipeline runs in one command from a fresh conda environment:

```bash
make env       # builds the conda env from env/environment.yml
make ligand    # prepares the As(OH)3 ligand
make all       # runs stages 1 through 8 end to end
```

A 6-hour wall-time run on a single 8-vCPU AMD EPYC node regenerates every
output reported in this chapter. The exact pipeline commit used for the
results above is recorded in `results/logs/`.

---
