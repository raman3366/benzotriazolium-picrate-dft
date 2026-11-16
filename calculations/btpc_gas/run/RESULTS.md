# DFT results: Benzotriazolium picrate (BTZP), CCDC 835950 -- gas phase

Generated 2026-09-17 00:31 by `dft/report.py`. All numbers below are read from the files in this folder; the interpretation paragraphs are drafts in the style of comparable papers and should be edited by the chemist.
**Compound.** 1H-benzotriazol-1-ium 2,4,6-trinitrophenolate (C12H8N6O7), CCDC 835950, original structure paper doi:10.1016/j.ijleo.2012.06.018. Model system: picrate anion + benzotriazolium cation (33 atoms).

**How this model was built.** Ion pair assembled from the CIF by symmetry: picrate at 1_555, benzotriazolium at 4_666. The two moieties of the asymmetric unit are NOT hydrogen bonded to each other (their closest approach is O4...N5 2.99 A, and N5 carries no hydrogen), so the pair had to be rebuilt from a symmetry image. The resulting N4-H4N...O1 bond (1.775 A, 154.5 deg) is to the picrate phenolate oxygen, the expected charge-assisted contact in a picrate salt. The cation's second N-H (N6-H6N) donates to a further symmetry-related picrate and is therefore not represented here.
## Computational methods (draft paragraph)
All calculations were performed with density functional theory at the B3LYP/6-311++G(d,p) level using PySCF 2.14.0 (B3LYP with the VWN-RPA local correlation, i.e. the same definition as in Gaussian; density fitting with the default auxiliary basis; integration grid level 3). Geometry optimisation used the geomeTRIC 1.1.1 optimiser with Gaussian's default convergence criteria (maximum force 4.5e-4 Ha/Bohr, RMS force 3.0e-4 Ha/Bohr), in the gas phase. Harmonic frequencies were computed by finite-difference (gradients, step 0.005 Bohr, 16 workers x 2 threads) at the same level to confirm a true minimum. Frontier molecular orbitals, the molecular electrostatic potential on the 0.001 a.u. electron-density isosurface, the reduced density gradient (NCI), the electron localisation function (ELF) and the localised orbital locator (LOL) were analysed with Multiwfn . Isosurfaces were rendered with ParaView 6.1.1 and 2-D maps with matplotlib 3.11.1.
## 1. Optimised geometry and frequency check
Final energy: -1317.208231 Ha after 3 optimisation steps (14 min). Bond lengths deviate from the X-ray structure by RMSD 0.092 Å (max 0.273 Å; X-ray H positions are systematically short, so the largest deviations are X-H bonds).
Frequencies: 93 modes, **0 imaginary** → the structure is a true minimum; lowest modes 12, 22, 42, 53, 54 cm⁻¹; ZPE 5.898 eV.
Proton check after the optimisation (gas phase): H4N#4_666–N4#4_666 1.109 Å, H4N#4_666···O1 1.437 Å → **the proton stays on N4#4_666, i.e. the ionic (salt) form is retained**.
![Optimised structure with the atom numbering of the CIF](fig/structure_labels.png)

*Optimised structure with the atom numbering of the CIF*

| Hydrogen bond | D–H (Å) | H···A (Å) | D···A (Å) | D–H···A (°) |
|---|---|---|---|---|
| N4#4_666-H4N#4_666...O1 (X-ray) | 0.907 | 1.775 | 2.623 | 154.5 |
| N4#4_666-H4N#4_666...O1 (DFT) | 1.109 | 1.437 | 2.478 | 153.2 |

Selected bond lengths (Å), heavy atoms only:

| Bond | X-ray | DFT | Δ |
|---|---|---|---|
| O1-C1 | 1.232 | 1.261 | +0.029 |
| N3-O7 | 1.214 | 1.232 | +0.018 |
| N3-C6 | 1.442 | 1.461 | +0.018 |
| N3-O6 | 1.230 | 1.226 | -0.004 |
| C6-C1 | 1.450 | 1.446 | -0.004 |
| C6-C5 | 1.367 | 1.385 | +0.018 |
| C1-C2 | 1.449 | 1.447 | -0.003 |
| C2-C3 | 1.352 | 1.375 | +0.022 |
| C2-N1 | 1.456 | 1.469 | +0.013 |
| C3-C4 | 1.382 | 1.395 | +0.013 |
| C4-C5 | 1.361 | 1.387 | +0.026 |
| C4-N2 | 1.441 | 1.461 | +0.020 |
| N1-O3 | 1.207 | 1.225 | +0.018 |
| N1-O2 | 1.186 | 1.227 | +0.041 |
| N2-O4 | 1.230 | 1.228 | -0.002 |
| N2-O5 | 1.216 | 1.227 | +0.011 |
| N4#4_666-C12#4_666 | 1.377 | 1.374 | -0.004 |
| N4#4_666-N5#4_666 | 1.324 | 1.294 | -0.030 |
| N6#4_666-C7#4_666 | 1.372 | 1.374 | +0.003 |
| N6#4_666-N5#4_666 | 1.328 | 1.326 | -0.002 |
| C7#4_666-C12#4_666 | 1.381 | 1.404 | +0.022 |
| C7#4_666-C8#4_666 | 1.384 | 1.400 | +0.016 |
| C12#4_666-C11#4_666 | 1.379 | 1.402 | +0.023 |
| C11#4_666-C10#4_666 | 1.370 | 1.381 | +0.011 |
| C10#4_666-C9#4_666 | 1.385 | 1.418 | +0.033 |
| C9#4_666-C8#4_666 | 1.360 | 1.383 | +0.023 |

Full list including X–H bonds: `geometry_vs_start.csv`.
## 2. Frontier molecular orbitals and global reactivity descriptors
| Descriptor | Formula | Value |
|---|---|---|
| E(HOMO) |  | -6.8181 eV |
| E(LUMO) |  | -3.3421 eV |
| Energy gap | E(LUMO) - E(HOMO) | 3.4760 eV |
| Ionization potential (IP) | -E(HOMO) | 6.8181 eV |
| Electron affinity (EA) | -E(LUMO) | 3.3421 eV |
| Electronegativity (chi) | (IP+EA)/2 | 5.0801 eV |
| Chemical potential (mu) | -(IP+EA)/2 | -5.0801 eV |
| Chemical hardness (eta) | (IP-EA)/2 | 1.7380 eV |
| Chemical softness (S) | 1/(2 eta) | 0.2877 eV^-1 |
| Electrophilicity index (omega) | mu^2/(2 eta) | 7.4245 eV |
| Dipole moment |  | 16.3639 Debye |

Fragment composition of the frontier orbitals (Mulliken partition, %):

| Orbital | picrate anion | benzotriazolium cation |
|---|---|---|
| HOMO-1 | 98 | 2 |
| HOMO | 99 | 1 |
| LUMO | 0 | 100 |
| LUMO+1 | 100 | 0 |

![Frontier orbital energy levels](fig/mo_levels.png)

*Frontier orbital energy levels*

![HOMO (isovalue ±0.02 a.u.)](fig/homo.png)

*HOMO (isovalue ±0.02 a.u.)*

![LUMO (isovalue ±0.02 a.u.)](fig/lumo.png)

*LUMO (isovalue ±0.02 a.u.)*

*Draft interpretation.* The HOMO and LUMO energies are -6.818 and -3.342 eV, giving an energy gap of 3.476 eV. The HOMO is localised mainly on the picrate anion (99 %) whereas the LUMO resides on the benzotriazolium cation (100 %), so the HOMO→LUMO excitation transfers charge between the two partners. From these, the ionisation potential (6.818 eV) and electron affinity (3.342 eV) yield an electronegativity of 5.080 eV, a chemical hardness of 1.738 eV and softness of 0.288 eV⁻¹; the electrophilicity index is 7.424 eV. The moderate gap indicates kinetic stability combined with significant polarisability, consistent with the interest in this salt for optical applications. The dipole moment is 16.36 D.
## 3. Molecular electrostatic potential (MEP)
ESP on the 0.001 a.u. density surface: global minimum -39.7 kcal/mol (-0.0632 a.u.), global maximum 81.9 kcal/mol (0.1306 a.u.); surface area 328.8 Å² (positive 178.5, negative 150.2); molecular polarity index 0.86 eV; polar surface fraction 70 %.
| Extremum | ESP (kcal/mol) | ESP (a.u.) | Nearest atoms |
|---|---|---|---|
| minimum (global) | -39.7 | -0.0632 | O5 (1.84 A), N2 (2.53 A) |
| minimum | -39.6 | -0.0631 | O4 (1.84 A), N2 (2.52 A) |
| minimum | -33.9 | -0.0540 | O6 (1.83 A), N3 (2.62 A) |
| minimum | -32.7 | -0.0521 | O3 (1.82 A), N1 (2.67 A) |
| maximum (global) | 81.9 | 0.1306 | H6N#4_666 (1.19 A), N6#4_666 (2.18 A) |
| maximum | 38.8 | 0.0619 | N6#4_666 (2.00 A), N5#4_666 (2.04 A) |
| maximum | 37.4 | 0.0595 | H9#4_666 (1.33 A), C9#4_666 (2.40 A) |
| maximum | 37.2 | 0.0593 | N5#4_666 (2.03 A), N6#4_666 (2.03 A) |

![MEP mapped on the 0.001 a.u. electron-density isosurface (red: negative, blue: positive)](fig/mep.png)

*MEP mapped on the 0.001 a.u. electron-density isosurface (red: negative, blue: positive)*

*Draft interpretation.* The most negative potential (-39.7 kcal/mol) lies near O5 (1.84 A), identifying the oxygen lone pairs there as the preferred sites for electrophilic attack and as hydrogen-bond acceptors. The most positive potential (81.9 kcal/mol) surrounds H6N#4_666 (1.19 A), which is therefore the nucleophilic-attack and hydrogen-bond-donor site. The strongly polarised surface (polar fraction 70 %) reflects the ionic character of the assembly and the electrostatic complementarity that drives the hydrogen bonding observed in the crystal. The molecular polarity index is 0.86 eV.
## 4. Reduced density gradient (RDG) / non-covalent interactions
Spikes in the RDG vs sign(λ₂)ρ plot (minimum RDG per bin of sign(λ₂)ρ):

| sign(λ₂)ρ (a.u.) | min RDG | interaction type |
|---|---|---|
| -0.086 | 0.08 | attractive (H-bond) |
| -0.083 | 0.18 | attractive (H-bond) |
| -0.081 | 0.22 | attractive (H-bond) |
| -0.017 | 0.03 | attractive (H-bond) |
| -0.011 | 0.11 | van der Waals |
| -0.005 | 0.10 | van der Waals |
| +0.004 | 0.13 | van der Waals |
| +0.013 | 0.04 | steric repulsion |
| +0.050 | 0.25 | steric repulsion |

![RDG versus sign(λ₂)ρ](fig/rdg_scatter.png)

*RDG versus sign(λ₂)ρ*

![RDG versus sign(λ₂)ρ, extended range showing the strong hydrogen-bond spike](fig/rdg_scatter_wide.png)

*RDG versus sign(λ₂)ρ, extended range showing the strong hydrogen-bond spike*

![NCI isosurface (RDG = 0.5 a.u.) coloured by sign(λ₂)ρ from −0.035 (blue) to +0.02 a.u. (red)](fig/nci.png)

*NCI isosurface (RDG = 0.5 a.u.) coloured by sign(λ₂)ρ from −0.035 (blue) to +0.02 a.u. (red)*

Note: Multiwfn's default screening (RDG set to 100 where ρ ≥ 0.05 a.u.) would remove this very short hydrogen bond from the map, so the threshold was raised to 0.1 a.u.; covalent bond regions remain screened.

*Draft interpretation.* The RDG scatter plot shows a spike at strongly negative sign(λ₂)ρ, which corresponds in the 3-D map to the blue disc at N4#4_666-H4N#4_666···O1: strong, attractive hydrogen bond. Green, low-density regions between the two partners and around the rings mark van der Waals contacts, while the red discs at the centres of the aromatic rings indicate steric (ring-closure) repulsion. Only the hydrogen bonds inside this cluster are present; the remaining hydrogen bonds of the crystal involve further symmetry-related neighbours and cannot appear here. The hydrogen-bond spike appears at sign(λ₂)ρ = -0.086 a.u. (minimum RDG 0.08); a value below about −0.05 a.u. indicates an unusually strong, partly covalent hydrogen bond, consistent with the short H···A distance.
## 5. Electron localisation function (ELF)
![ELF colour-filled map in the picrate ring plane](fig/elf_picrate_ring.png)

*ELF colour-filled map in the picrate ring plane*

![ELF colour-filled map in the benzotriazolium ring plane](fig/elf_benzotriazolium_ring.png)

*ELF colour-filled map in the benzotriazolium ring plane*

*Draft interpretation.* ELF values close to 1 (red) appear on the X–H bonds and on the lone-pair regions of the electronegative atoms (oxygen, nitrogen, chlorine), showing highly localised electron pairs; intermediate values (green) along the ring C–N and C–C bonds reflect delocalised π bonding within the rings, and values near zero (blue) at the ring centres and in the region between the fragments indicate electron depletion.
## 6. Localised orbital locator (LOL)
![LOL colour-filled map in the picrate ring plane](fig/lol_picrate_ring.png)

*LOL colour-filled map in the picrate ring plane*

![LOL colour-filled map in the benzotriazolium ring plane](fig/lol_benzotriazolium_ring.png)

*LOL colour-filled map in the benzotriazolium ring plane*

*Draft interpretation.* The LOL maps mirror the ELF picture with sharper contrast: high LOL (red) is concentrated on the covalent bond axes and the lone pairs, the delocalised π framework appears as continuous intermediate-LOL regions around the ring perimeters, and the low-LOL regions at the ring centres and in the hydrogen-bond region between the fragments show where electrons are delocalised or depleted.
## Files

`opt.xyz` optimised geometry · `geometry_vs_start.csv` bond lengths vs X-ray · `freq.json` frequencies · `fmo.json`/`fmo_table.md` descriptors · `mep_extrema.csv`/`mep_surface.json` ESP surface analysis · `nci_spikes.csv` RDG spikes · `fig/` all figures · `STAGES.md` stage log · `logs/` raw program output.
