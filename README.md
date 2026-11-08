# PFAS multi-receptor screening — code and data

Code, structures, and source data for:

> **Validated multi-target docking does not resolve nuclear receptor engagement by emerging PFAS replacements: limits of binding-based hazard prioritization**
>
> Jin Yoo<sup>†</sup>, Donghyeon Kim<sup>†</sup>, Jeongyun Kim, HanYeol Nam, Doo-Hee Lee, Ku Kang\*
>
> <sup>†</sup>These authors contributed equally.  \*Corresponding author: Ku Kang

Every value reported in the manuscript can be reproduced from this repository. The one exception is instrument-level luminescence output (Tecan Spark raw files), which is held by the corresponding author and provided on request; the replicate-level fold-induction values derived from it are deposited here in full.

---

## What this study asks

Structure-based screening is widely proposed as a hazard-prioritization layer for the >13,000 registered PFAS. It rests on an assumption that is rarely tested: that computed binding affinity predicts the functional consequence of the interaction.

We tested it directly. Four emerging replacement PFAS (GenX, HFPO-TA, F-53B, PFECHS) were assayed for activation of PXR and PPAR-α, and the same compounds were docked against a panel of validated human targets.

**Docked PXR/PPAR-α differences fall inside the scoring function's error, so docking does not resolve receptor engagement for these compounds. Extending to 100 ns explicit-solvent simulation does not change this. And all four compounds score highest against TR-α or ER-α — so a panel ranked by docking would not have selected either of the receptors the assays show respond.**

> **Note on assay scope.** PXR was measured in HepG2 on a native CYP3A4-XREM promoter and PPAR-α in COS-1 on a GAL4-LBD chimera. These differ in cell background, construct, promoter and expression level, so response magnitudes are **not comparable between the two receptors**, and a ratio taken within a compound does not cancel those differences. No receptor-preference assignment is made in this repository or in the manuscript. Earlier versions of this README reported such a comparison; it has been withdrawn.

---

## Quick start

```bash
pip install rdkit meeko gemmi numpy scipy pandas matplotlib
# AutoDock Vina 1.2.5 binary required for the docking step

python3 scripts/00_verify_ligand_structures.py   # structure verification
python3 scripts/01_validate_redocking.py         # redocking validation gate
python3 scripts/02_dock.py                       # production docking (28 × 7)
python3 scripts/03_statistics.py                 # correlation and leverage analysis
python3 scripts/04_transactivation_stats.py      # reporter assay statistics
python3 scripts/05_make_figures_archtox.py       # regenerate all manuscript figures
```

Scripts 00, 03, 04 run in seconds and require no external binaries.

---

## Layout

```
├── config/
│   ├── ligands_anionic.json          SMILES, charge, formula and MW for 28 compounds
│   ├── config_verified.json          8 receptors: PDB, chain, box centre and size
│   └── ligands_neutral_reference.csv neutral forms (sensitivity analysis)
├── ligands/                          PDBQT for 28 PFAS (anionic at pH 7.4)
├── receptors/                        receptor PDBQT
├── reference_ligands/                co-crystallized ligands (validation gate)
├── data/
│   ├── docking/
│   │   ├── dG_matrix.csv                       28 × 7 ΔG (primary result)
│   │   ├── dG_matrix_neutral_sensitivity.csv   recomputed with neutral forms
│   │   └── analysis_input.csv                  s_bind, nF, nC, hazard_H
│   └── experimental/
│       ├── transactivation_source_data.csv     576 replicate-level measurements
│       └── transactivation_summary.csv         64 conditions: mean, SD, SEM
├── results/
│   └── validation_gate.json          redocking RMSD and pass/fail
├── md/
│   ├── analysis/                     trajectory analysis output for 4 systems (.xvg)
│   ├── scripts/                      system build, run and analysis
│   └── mdp/                          GROMACS parameters
├── scripts/
└── figures/
```

---

## Two verification steps that this study depends on

### 1. Ligand structure verification (`scripts/00_...`)

Every SMILES is checked against the registered molecular formula and molecular weight for its CAS number, and electron parity is confirmed. An incorrectly fluorinated PFAS produces docking scores that look entirely reasonable while describing a different molecule; nothing downstream detects this.

All 28 ligands pass.

### 2. Redocking validation gate (`scripts/01_...`)

Each receptor must reproduce the crystallographic pose of its own co-crystallized ligand within 2.0 Å (symmetry-corrected heavy-atom RMSD, computed by optimal assignment within element classes) before any PFAS is docked into it.

| Target | PDB | Chain | Reference ligand | RMSD (Å) | |
|---|---|---|---|---:|---|
| HSA | 2BXD | A | warfarin | **0.45** | pass |
| PPAR-α | 1K7L | A | GW409544 | **0.52** | pass |
| TR-α | 2H79 | A | T3 | **0.70** | pass |
| ER-α | 1GWR | A | estradiol | **0.70** | pass |
| TTR | 1IE4 | A+C | thyroxine | **1.19** | pass |
| PPAR-γ | 1FM6 | **D** | rosiglitazone | **1.97** | pass |
| PXR | 1ILH | A | SR12813 | centroid 0.94 | pass |
| **L-FABP** | 3B2H | A | palmitate | **12.68** | **excluded** |

All values are symmetry-corrected heavy-atom RMSD from `results/validation_gate.json`,
regenerated by `scripts/01_validate_redocking.py`.

Notes:

- **HSA** was placed at Sudlow site I (subdomain IIA). The fatty-acid sites were not adopted because palmitate, their only reference ligand, has 14 rotatable bonds and occupies seven crystallographic sites and so does not give a well-defined pose to reproduce; warfarin at Sudlow site I does, and passes. That earlier attempt used PDB 1E7H and is **not** part of the deposited gate: `validation_gate.json` records only the eight sites actually carried into the screen.
- **PPAR-γ** requires chain D. 1FM6 is a PPARγ–RXRα heterodimer and chain A is RXRα.
- **L-FABP** could not be validated at any tested box size. Its only natural ligands are flexible fatty acids in a bilobed cavity. It was dropped, leaving **seven targets**.

---

## Principal results

### Docking does not resolve receptor engagement

| Compound | ΔG PXR | ΔG PPAR-α | ΔΔG | Strongest-binding target |
|---|---:|---:|---:|---|
| GenX | −7.74 | −7.58 | −0.17 | TR-α (−8.21) |
| PFECHS | −9.37 | −8.57 | −0.80 | ER-α (−10.16) |
| F-53B | −9.29 | −7.87 | −1.42 | TR-α (−9.76) |
| HFPO-TA | −8.70 | −8.28 | −0.42 | ER-α (−9.17) |

All four ΔΔG values fall inside the AutoDock Vina error range (~1.5 kcal/mol), so
**docking cannot discriminate** between the two receptors for these compounds.
Repeating the screen with neutral ligand forms gives the same outcome, so this is
not an artifact of protonation state.

Note the last column: none of the four compounds binds PXR or PPAR-α most strongly.
A panel assembled by ranking docking scores would have selected TR-α and ER-α, for
which no functional data exist here, and would have discarded the two receptors that
the assays below show respond. This is reported as hypothesis-generating.

### Measured activation, within each assay

| Compound | Assay | EC₅₀ (µM) | E_max |
|---|---|---:|---:|
| PFECHS | PXR (HepG2, CYP3A4-XREM) | 17.0 ± 1.5 | 10.3 ± 0.4 |
| GenX | PXR (HepG2, CYP3A4-XREM) | 38.1 ± 12.1 | 10.9 ± 1.4 |
| F-53B | PPAR-α (COS-1, GAL4-LBD) | 23.6 ± 2.5 | 7.5 ± 0.3 |
| HFPO-TA | PPAR-α (COS-1, GAL4-LBD) | 32.8 ± 10.5 | 11.6 ± 1.4 |

All four compounds activate **both** receptors concentration-dependently; the rows
above give the fit for one assay each and are **not** a preference assignment. The
two assays are not comparable to one another (see the note at the top).

Each condition is three independent biological replicates measured in technical
triplicate. **Technical triplicates are averaged within each biological replicate
before fitting and testing**, so the unit of analysis is n = 3, not N = 9 wells.
Treating technical replicates as independent inflates significance by three to six
orders of magnitude and is not done here. Uncertainties are fit standard errors;
`scripts/04_transactivation_stats.py` reports 95% confidence intervals,
Holm–Bonferroni adjusted p-values, and per-condition technical and biological RSD.

GenX–PXR and HFPO-TA–PPAR-α did not fully saturate at 100 µM (75→100 µM increments
of 9.6% and 11.5%), so those two EC₅₀ values carry extrapolation uncertainty.
Viability was not measured in parallel, so responses at the upper concentrations
cannot be separated from a cytotoxicity contribution.

### Ensemble simulation does not add resolution either

Four complexes were simulated for 100 ns: PFECHS and F-53B — the compounds with the
smallest and largest docked ΔΔG — each with PXR and PPAR-α. All four remained stably
bound (ligand RMSD 0.07–0.27 nm, minimum distance 0.19–0.20 nm), so neither receptor
rejects either ligand.

| System | Interaction energy (kJ/mol) | Backbone RMSD (nm) |
|---|---:|---:|
| PXR–PFECHS | −211.1 ± 22.9 | 0.366 |
| PPAR-α–PFECHS | −194.8 ± 13.7 | 0.186 |
| PXR–F-53B | −253.2 ± 19.7 | 0.560 |
| PPAR-α–F-53B | −199.4 ± 20.8 | 0.283 |

**These values must not be differenced across receptors.** Short-range Coulomb +
Lennard-Jones interaction energy is an enthalpic term that scales with contact
number, so a larger or more enclosed pocket yields a more negative value for reasons
unrelated to affinity. PXR and PPAR-α differ in pocket volume, and the PXR systems
are additionally truncated (see below). A between-receptor difference formed from
this column is not a receptor assignment, and earlier versions of this README
presented one; it has been withdrawn.

One trajectory was run per complex, so there is no independent replicate variance.
See `md/README.md` for the PXR chain-break issue that had to be resolved first, and
for why interaction energies should not be read as free energies.

### Aggregate binding tracks chain length, not independent hazard

| Quantity | Value |
|---|---|
| Spearman ρ (s_bind vs hazard_H) | +0.793 (p = 8.1 × 10⁻⁷) |
| ρ(s_bind, n_F) | +0.936 |
| **Partial ρ controlling for n_F** | **+0.275 (p = 0.175, n.s.)** |
| Bootstrap 95% CI | [−0.248, +0.685] |
| Fisher top-10 overlap | 7/10 (OR 10.89, p = 0.013) |

The partial correlation is not significant and is unstable: removing PFECHS alone moves it to +0.116. With predictor–control collinearity of 0.936 at n = 27 there is not enough power to resolve an independent contribution in either direction. We therefore report that binding strength provides **no demonstrable hazard information beyond chain length**, without claiming that none exists.

`scripts/03_statistics.py` prints the full leave-one-out leverage table.

---

## Methods summary

| | |
|---|---|
| Ligand preparation | RDKit ETKDGv3 + MMFF94, Meeko 0.6; carboxylates and sulfonates deprotonated at pH 7.4, sulfonamides neutral |
| Receptor preparation | Open Babel 3.1.1, hydrogens at pH 7.4, single chain unless the site spans an interface |
| Docking | AutoDock Vina 1.2.5, exhaustiveness 16, 9 modes, seed 42 |
| Reporter assays | HepG2 (pSG5-hPXR + CYP3A4-XREM-luc + pRL-TK); COS-1 (pBIND-PPAR-α-LBD GAL4 + pG5-Luc + pRL-TK); 8 concentrations 0.3–100 µM; Dual-Glo, Tecan Spark |
| Statistics | Spearman with 5,000-iteration bootstrap; partial Spearman by the three-variable formula, verified by rank-residual regression; paired tests on biological-replicate means with Holm–Bonferroni correction; leave-one-out leverage |

---

## Data provenance

Reporter assay measurements were recorded on 2026-05-15. `data/experimental/transactivation_source_data.csv` contains all 576 replicate-level fold-induction values as measured. Instrument-level output (Tecan Spark luminescence) is held by the corresponding author and can be provided to editors or reviewers on request.

Simulated dose–response datasets were generated during analysis development to build the fitting and plotting code. They are **not** included here and are not used for any reported quantity.

---

## Declaration of generative AI use

An AI assistant (Anthropic Claude) was used to write and debug the analysis, docking-pipeline, and figure-generation scripts in this repository, and to assist in drafting and editing manuscript text. No data were generated, simulated, or altered by AI for any reported result; no references were produced by AI; and no scientific interpretation was delegated to AI. All authors verified the reported values against the source records and take full responsibility for the content.

---

## License

Code: MIT (see `LICENSE`). Data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Contact

Ku Kang (corresponding author) — CBRN Defense Research Institute, Ministry of National Defense, Seoul 06796, Republic of Korea — bisu9082@gmail.com

Repository maintained by Jin Yoo — tysonyj@snu.ac.kr
