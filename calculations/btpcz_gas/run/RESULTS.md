# DFT results: Benzotriazolium picrate, P2_1/n polymorph (Zeng et al. 2011, COD 2230422) -- gas phase

Generated 2026-11-22 01:46 by `dft/report.py`. All numbers below are read from the files in this folder; the interpretation paragraphs are drafts in the style of comparable papers and should be edited by the chemist.
**Compound.** 1H-benzotriazol-1-ium 2,4,6-trinitrophenolate (C12H8N6O7), CCDC ?, original structure paper doi:10.1107/S1600536811018307. Model system: picrate anion + benzotriazolium cation (33 atoms).

**How this model was built.** Second polymorph of BTPC (revision item C0): Zeng, Li & Wang, Acta Cryst. E67, o1464 (2011), doi:10.1107/S1600536811018307, COD 2230422, space group P2_1/n, a 14.4113, b 3.7608, c 24.941 A, beta 90.598 deg, Z 4, 298 K, R1 0.0616. The published paper models the P2_1/c form (CCDC 835950, case btpc). Here the asymmetric unit IS the hydrogen-bonded ion pair: N6-H6...O1 (phenolate), H...A 1.870 A, N...O 2.770 A, 160.0 deg, matching the deposited _geom_hbond loop (2.770(3) A, 160(3) deg). The cation carries H on N4 and N6 and none on N5, as in btpc (1H,3H-benzotriazolium); the deposited N-H distances are refined X-ray values (N4-H4 0.85, N6-H6 0.94 A). The cation's other N-H (N4-H4) donates a bifurcated bond to O1 and the ortho-nitro O2 of a 2_656 picrate (N...O 2.799 and 2.791 A); it is not represented here. Unconstrained, build_cluster.py would have picked that bifurcated image (two contacts outscore one), so the pair was built with --require-hbond N6:O1. No CCDC number is given in the COD CIF.
## Computational methods (draft paragraph)
All calculations were performed with density functional theory at the B3LYP/6-311++G(d,p) level using PySCF 2.14.0 (B3LYP with the VWN-RPA local correlation, i.e. the same definition as in Gaussian; density fitting with the default auxiliary basis; integration grid level 3). Geometry optimisation used the geomeTRIC 1.1.1 optimiser with Gaussian's default convergence criteria (maximum force 4.5e-4 Ha/Bohr, RMS force 3.0e-4 Ha/Bohr), in the gas phase. Harmonic frequencies were computed by finite-difference (gradients, step 0.005 Bohr, dft/_fd_worker.py in 33 sequential chunks of 3 columns x 4 threads; dft/revision/c0_case.py) at the same level to confirm a true minimum. Frontier molecular orbitals, the molecular electrostatic potential on the 0.001 a.u. electron-density isosurface, the reduced density gradient (NCI), the electron localisation function (ELF) and the localised orbital locator (LOL) were analysed with Multiwfn . Isosurfaces were rendered with ParaView 6.1.1 and 2-D maps with matplotlib 3.11.2.
## 1. Optimised geometry and frequency check
Final energy: -1317.211397 Ha after 99 optimisation steps (549 min). Bond lengths deviate from the X-ray structure by RMSD 0.168 Å (max 0.885 Å; X-ray H positions are systematically short, so the largest deviations are X-H bonds).
Frequencies: 93 modes, **1 imaginary** → NOT a minimum; lowest modes -59, 5, 17, 28, 47 cm⁻¹; ZPE 5.905 eV.
Proton check after the optimisation (gas phase): H6–N6 1.823 Å, H6···O1 1.000 Å → **the proton migrated from N6 to O1, so the optimised species is the neutral hydrogen-bonded complex rather than the salt**.
Connectivity change during the optimisation: bonds broken ['N6-H6'], formed ['O1-H6'].
![Optimised structure with the atom numbering of the CIF](fig/structure_labels.png)

*Optimised structure with the atom numbering of the CIF*

| Hydrogen bond | D–H (Å) | H···A (Å) | D···A (Å) | D–H···A (°) |
|---|---|---|---|---|
| N6-H6...O1 (X-ray) | 0.938 | 1.870 | 2.770 | 160.0 |
| N6-H6...O1 (DFT) | 1.823 | 1.000 | 2.708 | 145.6 |

Selected bond lengths (Å), heavy atoms only:

| Bond | X-ray | DFT | Δ |
|---|---|---|---|
| C1-C2 | 1.446 | 1.420 | -0.026 |
| C1-C6 | 1.454 | 1.423 | -0.032 |
| C1-O1 | 1.251 | 1.308 | +0.056 |
| C2-C3 | 1.372 | 1.375 | +0.003 |
| C2-N1 | 1.462 | 1.482 | +0.020 |
| C3-C4 | 1.372 | 1.394 | +0.022 |
| C4-C5 | 1.377 | 1.381 | +0.004 |
| C4-N2 | 1.460 | 1.472 | +0.013 |
| C5-C6 | 1.371 | 1.391 | +0.020 |
| C6-N3 | 1.455 | 1.466 | +0.011 |
| N1-O2 | 1.218 | 1.217 | -0.001 |
| N1-O3 | 1.218 | 1.223 | +0.005 |
| N2-O4 | 1.222 | 1.224 | +0.002 |
| N2-O5 | 1.219 | 1.224 | +0.005 |
| N3-O6 | 1.225 | 1.230 | +0.005 |
| N3-O7 | 1.222 | 1.221 | -0.001 |
| C7-C8 | 1.399 | 1.400 | +0.002 |
| C7-C12 | 1.390 | 1.407 | +0.017 |
| C7-N4 | 1.365 | 1.367 | +0.002 |
| C8-C9 | 1.367 | 1.384 | +0.017 |
| C9-C10 | 1.414 | 1.415 | +0.001 |
| C10-C11 | 1.356 | 1.383 | +0.027 |
| C11-C12 | 1.395 | 1.404 | +0.009 |
| C12-N6 | 1.370 | 1.382 | +0.012 |
| N4-N5 | 1.307 | 1.347 | +0.040 |
| N5-N6 | 1.311 | 1.292 | -0.019 |

Full list including X–H bonds: `geometry_vs_start.csv`.
## 2. Frontier molecular orbitals and global reactivity descriptors
| Descriptor | Formula | Value |
|---|---|---|
| E(HOMO) |  | -7.4302 eV |
| E(LUMO) |  | -3.5918 eV |
| Energy gap | E(LUMO) - E(HOMO) | 3.8383 eV |
| Ionization potential (IP) | -E(HOMO) | 7.4302 eV |
| Electron affinity (EA) | -E(LUMO) | 3.5918 eV |
| Electronegativity (chi) | (IP+EA)/2 | 5.5110 eV |
| Chemical potential (mu) | -(IP+EA)/2 | -5.5110 eV |
| Chemical hardness (eta) | (IP-EA)/2 | 1.9192 eV |
| Chemical softness (S) | 1/(2 eta) | 0.2605 eV^-1 |
| Electrophilicity index (omega) | mu^2/(2 eta) | 7.9126 eV |
| Dipole moment |  | 8.9079 Debye |

Fragment composition of the frontier orbitals (Mulliken partition, %):

| Orbital | picrate anion (protonation changed) | benzotriazolium cation (protonation changed) |
|---|---|---|
| HOMO-1 | 0 | 100 |
| HOMO | 0 | 100 |
| LUMO | 100 | 0 |
| LUMO+1 | 100 | 0 |

![Frontier orbital energy levels](fig/mo_levels.png)

*Frontier orbital energy levels*

*(figure fig/homo.png not generated)*

*(figure fig/lumo.png not generated)*

*Draft interpretation.* The HOMO and LUMO energies are -7.430 and -3.592 eV, giving an energy gap of 3.838 eV. The HOMO is localised mainly on the benzotriazolium cation (protonation changed) (100 %) whereas the LUMO resides on the picrate anion (protonation changed) (100 %), so the HOMO→LUMO excitation transfers charge between the two partners. From these, the ionisation potential (7.430 eV) and electron affinity (3.592 eV) yield an electronegativity of 5.511 eV, a chemical hardness of 1.919 eV and softness of 0.261 eV⁻¹; the electrophilicity index is 7.913 eV. The moderate gap indicates kinetic stability combined with significant polarisability, consistent with the interest in this salt for optical applications. The dipole moment is 8.91 D.
## 3. Molecular electrostatic potential (MEP)
ESP on the 0.001 a.u. density surface: global minimum -30.4 kcal/mol (-0.0484 a.u.), global maximum 63.4 kcal/mol (0.1010 a.u.); surface area 337.3 Å² (positive 199.2, negative 138.1); molecular polarity index 0.65 eV; polar surface fraction 66 %.
| Extremum | ESP (kcal/mol) | ESP (a.u.) | Nearest atoms |
|---|---|---|---|
| minimum (global) | -30.4 | -0.0484 | N5 (1.93 A), O1 (2.60 A) |
| minimum | -30.2 | -0.0481 | O1 (1.91 A), O2 (1.92 A) |
| minimum | -29.3 | -0.0467 | O4 (1.82 A), N2 (2.58 A) |
| minimum | -29.2 | -0.0466 | O5 (1.82 A), N2 (2.58 A) |
| maximum (global) | 63.4 | 0.1010 | H4 (1.21 A), N4 (2.20 A) |
| maximum | 32.6 | 0.0520 | H6 (1.85 A), N3 (2.01 A) |
| maximum | 28.2 | 0.0449 | H3 (1.69 A), N1 (1.93 A) |
| maximum | 27.2 | 0.0434 | N3 (1.83 A), C6 (2.13 A) |

*(figure fig/mep.png not generated)*

*Draft interpretation.* The most negative potential (-30.4 kcal/mol) lies near N5 (1.93 A), identifying the nitrogen lone pairs there as the preferred sites for electrophilic attack and as hydrogen-bond acceptors. The most positive potential (63.4 kcal/mol) surrounds H4 (1.21 A), which is therefore the nucleophilic-attack and hydrogen-bond-donor site. After the proton transfer the surface is much less polarised than in the salt form. The molecular polarity index is 0.65 eV.
## 4. Reduced density gradient (RDG) / non-covalent interactions
Spikes in the RDG vs sign(λ₂)ρ plot (minimum RDG per bin of sign(λ₂)ρ):

| sign(λ₂)ρ (a.u.) | min RDG | interaction type |
|---|---|---|
| -0.038 | 0.01 | attractive (H-bond) |
| -0.021 | 0.04 | attractive (H-bond) |
| +0.004 | 0.15 | van der Waals |
| +0.013 | 0.10 | steric repulsion |
| +0.021 | 0.05 | steric repulsion |
| +0.050 | 0.25 | steric repulsion |

![RDG versus sign(λ₂)ρ](fig/rdg_scatter.png)

*RDG versus sign(λ₂)ρ*

*(figure fig/rdg_scatter_wide.png not generated)*

*(figure fig/nci.png not generated)*

*Draft interpretation.* The RDG scatter plot shows a spike at strongly negative sign(λ₂)ρ, which corresponds in the 3-D map to the blue disc at N6-H6···O1: strong, attractive hydrogen bond. Green, low-density regions between the two partners and around the rings mark van der Waals contacts, while the red discs at the centres of the aromatic rings indicate steric (ring-closure) repulsion. Only the hydrogen bonds inside this cluster are present; the remaining hydrogen bonds of the crystal involve further symmetry-related neighbours and cannot appear here. The hydrogen-bond spike appears at sign(λ₂)ρ = -0.038 a.u. (minimum RDG 0.01); a value below about −0.05 a.u. indicates an unusually strong, partly covalent hydrogen bond, consistent with the short H···A distance.
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
