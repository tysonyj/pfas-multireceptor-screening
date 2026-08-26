# Molecular dynamics

Four protein–ligand complexes simulated for 100 ns in explicit solvent:
PFECHS and F-53B — which show opposite functional receptor preferences —
each in complex with PXR and PPAR-α.

## Systems

| System | Receptor | Ligand | Backbone RMSD (nm) | Ligand RMSD (nm) | Interaction energy (kJ mol⁻¹) |
|---|---|---|---:|---:|---:|
| `PXR_PFECHS_v2` | PXR (1ILH, 198–431) | PFECHS | 0.366 | 0.163 | −211.1 ± 22.9 |
| `PPARa_PFECHS` | PPAR-α (1K7L chain A) | PFECHS | 0.186 | 0.070 | −194.8 ± 13.7 |
| `PXR_F-53B_v2` | PXR (1ILH, 198–431) | F-53B | 0.560 | 0.263 | −253.2 ± 19.7 |
| `PPARa_F-53B` | PPAR-α (1K7L chain A) | F-53B | 0.283 | 0.272 | −199.4 ± 20.8 |

Values are means over the final 30% of the trajectory. Interaction energies are
short-range Coulomb plus Lennard-Jones over the 70–100 ns plateau.

## The PXR chain break

The PXR crystal structure 1ILH lacks residues 178–197. GROMACS treats the two
resulting chain segments as separate molecules, and in preliminary simulations
they rotated relative to one another as rigid bodies until backbone RMSD reached
2.3–2.7 nm.

**This failure is invisible in the metrics normally reported.** Radius of gyration
stayed constant to within 2%, and the closest approach between the two segments
remained 0.19 nm throughout — the protein neither unfolded nor came apart. Only an
explicit backbone RMSD check against the starting structure revealed it.

The N-terminal segment (142–177) contributes none of the twelve residues within
4.5 Å of the docked ligand, and the disordered region lies 15–24 Å from the
binding site, so it was removed and residues 198–431 simulated. `scripts/rebuild_pxr.sh`
performs this by detecting the break from C–N backbone distances rather than
residue numbering, which is unreliable after format conversion.

PPAR-α has no chain break and was simulated intact.

## Result

Both compounds interact more favourably with PXR: PFECHS by 16.3 ± 5.4 and
F-53B by 53.8 ± 3.5 kJ mol⁻¹ (mean ± standard error over 5 ns blocks). The sign
is consistent across all six blocks for both compounds.

Against the assay this is correct for PFECHS and incorrect for F-53B, which is
assigned to PXR by 53.8 kJ mol⁻¹ (≈12.9 kcal mol⁻¹) while activating PPAR-α
2.9-fold over PXR in cells.

Interaction energies are enthalpic terms that omit entropy and desolvation and
scale with contact number. They are not free energies and should not be read as
affinity estimates. Their role here is to establish that extending from a rigid
docked pose to a converged explicit-solvent ensemble does not recover the measured
receptor assignment.

## Settings

| | |
|---|---|
| Force field | AMBER99SB-ILDN (protein) / GAFF2 + AM1-BCC (ligand, acpype) |
| Water | TIP3P, dodecahedral box, 1.2 nm minimum solute–boundary distance |
| Ions | 0.15 M NaCl, neutralized |
| Equilibration | NVT 100 ps → NPT 500 ps, 310 K, 1 bar (V-rescale, C-rescale) |
| Production | 100 ns, 2 fs timestep, LINCS on H-bonds, PME, 1.2 nm cutoffs |
| Ligand charge | −1 (sulfonate deprotonated at pH 7.4) |

## Reproducing

```bash
bash scripts/rebuild_pxr.sh --run   # PXR systems (segment removal + MD)
bash scripts/run_md.sh PPARa_PFECHS 0
bash scripts/analyze.sh
python3 scripts/summarize.py
```

Raw trajectories (~200 GB) are not deposited. All analysis output required to
reproduce the reported values is in `analysis/`.
