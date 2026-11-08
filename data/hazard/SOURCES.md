# Regulatory priority status — sources and scoring rule

Status recorded **as of 2026-11-08**. Regulatory status is time-varying; the
"as of" date is part of the datum.

## Scoring rule

`reg_score` is the **count of binding regulatory instruments in which the
compound is individually named, or falls within an explicit homologue range
defined in the instrument text.**

Two deliberate restrictions:

1. **No degradation-pathway judgements.** Several instruments extend to
   "related substances" defined functionally, as substances with the potential
   to degrade or be transformed into the named compound. Deciding whether a
   given fluorotelomer or sulfonamide qualifies is a scientific determination,
   not a documentary one, and different sources disagree. Compounds that are
   plausibly captured this way are scored 0 and flagged in the `notes` column,
   so a reader who wishes to score them differently can.
2. **Class-wide instruments are excluded from the count.** The proposed ECHA
   universal PFAS restriction (Annex XV, 2023) and REACH Annex XVII entry 82
   (PFAS in firefighting foams, Reg. (EU) 2025/1988) are structural class
   definitions that capture all 28 compounds in this panel. They carry no
   discriminating information and would add a constant to every row.

The score is therefore a conservative lower bound on regulatory attention.
It ranges 0–4 across this panel.

## Instruments scored

**`US_EPA_NPDWR_2024`** — PFAS National Primary Drinking Water Regulation,
89 FR 32532 (26 April 2024). Individual MCLs: PFOA 4.0 ng/L, PFOS 4.0 ng/L,
PFHxS 10 ng/L, PFNA 10 ng/L, HFPO-DA (GenX) 10 ng/L. A Hazard Index of 1
applies to mixtures of two or more of PFHxS, PFNA, HFPO-DA and PFBS; PFBS
appears **only** in the Hazard Index and has no individual MCL.
https://www.federalregister.gov/documents/2024/04/26/2024-07773/pfas-national-primary-drinking-water-regulation

*Status note.* On 18–20 May 2026 EPA published proposed rules to rescind the
regulatory determinations for PFHxS, PFNA, HFPO-DA and the Hazard Index, and to
extend PFOA/PFOS compliance to 2031. Comment periods closed 20 July 2026. **No
final rescission rule had been published as of 2026-11-08**, so the 2024 rule is
recorded as in force. If the rescission is finalised, PFHxS, PFNA, GenX and PFBS
each lose one instrument.

**`Stockholm_Annex_A_or_B`** — Stockholm Convention on Persistent Organic
Pollutants. PFOS and its salts, Annex B (SC-4/17, 2009). PFOA, its salts and
related compounds, Annex A (SC-9/12, 2019). PFHxS, its salts and related
compounds, Annex A (SC-10/13, 2022). Long-chain perfluorocarboxylic acids
C9–C21, their salts and related compounds, Annex A (SC-12/12, COP-12, May 2025)
— an explicit homologue range covering PFNA, PFDA, PFUnDA and PFDoDA in this
panel. https://www.pops.int/

**`EU_DWD_2020_2184_AnnexIII`** — Directive (EU) 2020/2184, Annex III Part B
point 3, "Sum of PFAS". A closed list of 20 individually named substances,
verified verbatim from the consolidated text on EUR-Lex: the C4–C13
perfluoroalkyl carboxylic acids and the C4–C13 perfluoroalkane sulfonic acids.
Fifteen of this panel's 28 compounds appear; no ether acid, fluorotelomer or
sulfonamide is on the list. https://eur-lex.europa.eu/eli/dir/2020/2184/oj/eng

**`REACH_AnnexXVII_68_or_79`** — Entry 68, C9–C14 perfluorocarboxylic acids,
their salts and related substances (Reg. (EU) 2021/1297) — an explicit homologue
range covering PFNA, PFDA, PFUnDA and PFDoDA here. Entry 79, undecafluorohexanoic
acid (PFHxA), its salts and PFHxA-related substances (Reg. (EU) 2024/2462) —
PFHxA is named individually.

## Not used, and why

- **ECHA REACH Candidate List (SVHC).** Several of these compounds are listed,
  but the live ECHA list could not be retrieved for verification at the time of
  compilation, and an unverified list would be the weakest row in the table.
  Omitted rather than guessed.
- **State-level and national limits** (US state MCLs, Michigan, Vermont, EU
  member-state values). Coverage is uneven and changes frequently; including a
  partial set would bias toward jurisdictions that happen to be well documented
  in English.
- **Predicted or modelled quantities of any kind.** The purpose of this variable
  is to be independent of molecular structure, so that correlating it against
  docking scores is not circular. Anything QSAR-derived would defeat that.

## What this variable is, and is not

It is a documentary record of regulatory attention. It is **not** a measure of
intrinsic toxicity, and it should not be read as one. Regulatory attention
reflects the accumulated weight of evidence, the economic significance of a
compound, the length of time it has been in commerce, and the availability of
analytical standards — none of which is a hazard property.

Its one methodological virtue for the present analysis is that it is not
derived from molecular structure. Composite hazard indices that incorporate
half-life, bioaccumulation or QSAR-predicted properties covary with
perfluorocarbon chain length by construction; correlating docking scores — which
over a homologous series are themselves dominated by chain length — against such
an index measures chain length against chain length. This variable does not have
that defect, and it is the reason it was chosen over the alternatives.

It does carry an evidence-availability bias in the opposite direction: a score of
0 means the compound has not attracted binding regulation, which for a recently
introduced replacement reflects absence of evidence at least as much as evidence
of absence. Twelve of the 28 compounds score 0, including four of the five
ether-linked and cyclic replacements. That pattern is reported as a result.
