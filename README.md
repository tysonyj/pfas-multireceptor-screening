# PFAS multi-receptor screening — code and data

Code, structures, and source data for:

> **Receptor-divergent activation of emerging PFAS replacement chemicals supports multi-receptor functional screening beyond binding-based prioritization**
>
> Ku Kang, Jeongyun Kim, Doo-Hee Lee, HanYeol Nam, Donghyeon Kim, Jin Yoo\*
>
> \*Corresponding author: tysonyj@snu.ac.kr

Everything needed to reproduce every reported number is in this repository. No component is available on request only.

---

## What this study asks

Structure-based screening is widely proposed as a hazard-prioritization layer for the >13,000 registered PFAS. It rests on an assumption that is rarely tested: that computed binding affinity predicts the functional consequence of the interaction.

We tested it directly. Four emerging ether-linked PFAS were assayed for receptor-selective activation of PXR and PPAR-α, and the same compounds were docked against a panel of validated human targets. **Docking reproduces the measured receptor preference for two of four compounds — chance level — with energy differences well inside the scoring function's error. Extending to converged 100 ns explicit-solvent simulation does not recover it either.**

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
| PPAR-α | 1K7L | A | GW409544 | **0.58** | pass |
| TR-α | 2H79 | A | T3 | **0.70** | pass |
| ER-α | 1GWR | A | estradiol | **0.72** | pass |
| TTR | 1IE4 | A+C | thyroxine | **1.23** | pass |
| PPAR-γ | 1FM6 | **D** | rosiglitazone | **1.86** | pass |
| PXR | 1ILH | A | SR12813 | centroid 0.94 | pass |
| **L-FABP** | 3B2H | A | palmitate | **4.51** | **excluded** |

Notes:

- **HSA** was placed at Sudlow site I (subdomain IIA). The fatty-acid sites failed validation with palmitate (RMSD 7.26 Å) because that reference ligand has 14 rotatable bonds and seven crystallographic sites; the drug-binding site is better defined for a screening application.
- **PPAR-γ** requires chain D. 1FM6 is a PPARγ–RXRα heterodimer and chain A is RXRα.
- **L-FABP** could not be validated at any tested box size. Its only natural ligands are flexible fatty acids in a bilobed cavity. It was dropped, leaving **seven targets**.

---

## Principal results

### Docking does not resolve receptor preference

| Compound | ΔG PXR | ΔG PPAR-α | ΔΔG | Docking | Assay | |
|---|---:|---:|---:|---|---|---|
| GenX | −7.74 | −7.58 | −0.17 | PXR | PXR | ✓ |
| PFECHS | −9.37 | −8.57 | −0.80 | PXR | PXR | ✓ |
| F-53B | −9.29 | −7.87 | −1.42 | PXR | PPAR-α | ✗ |
| HFPO-TA | −8.70 | −8.28 | −0.42 | PXR | PPAR-α | ✗ |

All four ΔΔG values fall inside the AutoDock Vina error range (~1.5 kcal/mol). The correct reading is that **docking cannot discriminate**, not that it predicts PXR. Repeating the screen with neutral ligand forms gives the same 2/4 agreement, so the result is not an artifact of protonation state.

All four compounds bind TR-α or ER-α more strongly than either PXR or PPAR-α. No functional data exist for those targets; this is reported as hypothesis-generating only.

### Measured receptor selectivity

| Compound | Preferred | EC₅₀ (µM) | E_max | Ratio PXR/PPAR-α | p (biological) |
|---|---|---:|---:|---:|---:|
| PFECHS | PXR | 17.0 ± 1.5 | 10.3 ± 0.4 | 6.54 | 7.0 × 10⁻⁴ |
| GenX | PXR | 38.1 ± 12.1 | 10.9 ± 1.4 | 4.24 | 4.1 × 10⁻³ |
| F-53B | PPAR-α | 23.6 ± 2.5 | 7.5 ± 0.3 | 0.35 | 2.8 × 10⁻² |
| HFPO-TA | PPAR-α | 32.8 ± 10.5 | 11.6 ± 1.4 | 0.33 | 4.5 × 10⁻⁷ |

N = 9 (three biological × three technical). **p-values are computed on the three
biological replicates**, since technical replicates are not independent observations.
The corresponding technical-replicate values (4.3 × 10⁻⁶ to 6.5 × 10⁻¹⁴) overstate
significance by three to six orders of magnitude; the direction of every preference
is identical under both treatments.

Uncertainties on EC₅₀ and E_max are fit standard errors. GenX–PXR and HFPO-TA–PPAR-α
did not fully saturate at 100 µM (75→100 µM increments of 9.6% and 11.5%), so those
two parameters carry extrapolation uncertainty.

### Ensemble simulation does not recover the measured preference either

Four complexes were simulated for 100 ns: PFECHS and F-53B, which show opposite
functional preferences, each with PXR and PPAR-α. All four remained stably bound
(ligand RMSD 0.07–0.27 nm, minimum distance 0.19–0.20 nm), so neither receptor
rejects either ligand.

| Compound | ΔIE (PXR − PPAR-α) | Block SEM | Sign consistency | Computed | Assay | |
|---|---:|---:|---:|---|---|---|
| PFECHS | −16.3 | ±5.4 | 6/6 | PXR | PXR | ✓ |
| F-53B | **−53.8** | ±3.5 | 6/6 | PXR | PPAR-α | ✗ |

The simulation does not merely fail to discriminate — it discriminates confidently
and, for F-53B, in the wrong direction, by roughly 12.9 kcal mol⁻¹. See `md/README.md`
for the PXR chain-break issue that had to be resolved first, and for why
interaction energies should not be read as free energies.

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
| Reporter assays | HepG2 (pSG5-hPXR + CYP3A4-XREM-luc + pRL-TK); COS-1 (pBIND-PPAR-α-LBD GAL4 + pG5-Luc + pRL-TK); 8 concentrations 0.3–100 µM; Dual-Glo, Tecan Spark; CellTiter-Glo counter-screen ≥85% viability |
| Statistics | Spearman with 5,000-iteration bootstrap; partial Spearman by the three-variable formula, verified by rank-residual regression; Welch-corrected tests; leave-one-out leverage |

---

## Data provenance

Reporter assay measurements were recorded on 2026-05-15. `data/experimental/transactivation_source_data.csv` contains all 576 replicate-level fold-induction values as measured. Instrument-level output (Tecan Spark luminescence, CellTiter-Glo viability) is held by the corresponding author and can be provided to editors or reviewers on request.

Simulated dose–response datasets were generated during analysis development to build the fitting and plotting code. They are **not** included here and are not used for any reported quantity.

---

## Declaration of generative AI use

An AI assistant (Anthropic Claude) was used to write and debug the analysis, docking-pipeline, and figure-generation scripts in this repository, and to assist in drafting and editing manuscript text. No data were generated, simulated, or altered by AI for any reported result; no references were produced by AI; and no scientific interpretation was delegated to AI. All authors verified the reported values against the source records and take full responsibility for the content.

---

## License

Code: MIT (see `LICENSE`). Data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Contact

Jin Yoo — CBRN Defense Research Institute, Ministry of National Defense, Seoul 06796, Republic of Korea — tysonyj@snu.ac.kr
