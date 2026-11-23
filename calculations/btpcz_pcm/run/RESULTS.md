# DFT results: Benzotriazolium picrate, P2_1/n polymorph (Zeng et al. 2011, COD 2230422) -- IEF-PCM (water)

Generated 2026-11-21 19:36 by `dft/report.py`. All numbers below are read from the files in this folder; the interpretation paragraphs are drafts in the style of comparable papers and should be edited by the chemist.
**Compound.** 1H-benzotriazol-1-ium 2,4,6-trinitrophenolate (C12H8N6O7), CCDC ?, original structure paper doi:10.1107/S1600536811018307. Model system: picrate anion + benzotriazolium cation (33 atoms).

**How this model was built.** Second polymorph of BTPC (revision item C0): Zeng, Li & Wang, Acta Cryst. E67, o1464 (2011), doi:10.1107/S1600536811018307, COD 2230422, space group P2_1/n, a 14.4113, b 3.7608, c 24.941 A, beta 90.598 deg, Z 4, 298 K, R1 0.0616. The published paper models the P2_1/c form (CCDC 835950, case btpc). Here the asymmetric unit IS the hydrogen-bonded ion pair: N6-H6...O1 (phenolate), H...A 1.870 A, N...O 2.770 A, 160.0 deg, matching the deposited _geom_hbond loop (2.770(3) A, 160(3) deg). The cation carries H on N4 and N6 and none on N5, as in btpc (1H,3H-benzotriazolium); the deposited N-H distances are refined X-ray values (N4-H4 0.85, N6-H6 0.94 A). The cation's other N-H (N4-H4) donates a bifurcated bond to O1 and the ortho-nitro O2 of a 2_656 picrate (N...O 2.799 and 2.791 A); it is not represented here. Unconstrained, build_cluster.py would have picked that bifurcated image (two contacts outscore one), so the pair was built with --require-hbond N6:O1. No CCDC number is given in the COD CIF. | Rerun with an IEF-PCM (water) continuum because the continuum model is the primary result for an ionic salt.
## Computational methods (draft paragraph)
All calculations were performed with density functional theory at the B3LYP/6-311++G(d,p) level using PySCF 2.14.0 (B3LYP with the VWN-RPA local correlation, i.e. the same definition as in Gaussian; density fitting with the default auxiliary basis; integration grid level 3). Geometry optimisation used the geomeTRIC 1.1.1 optimiser with Gaussian's default convergence criteria (maximum force 4.5e-4 Ha/Bohr, RMS force 3.0e-4 Ha/Bohr), in an IEF-PCM continuum of water. Harmonic frequencies were computed analytically at the same level to confirm a true minimum. Frontier molecular orbitals, the molecular electrostatic potential on the 0.001 a.u. electron-density isosurface, the reduced density gradient (NCI), the electron localisation function (ELF) and the localised orbital locator (LOL) were analysed with Multiwfn . Isosurfaces were rendered with ParaView 6.1.1 and 2-D maps with matplotlib 3.11.2.
## 1. Optimised geometry and frequency check
Final energy: -1317.253734 Ha after 55 optimisation steps (762 min). Bond lengths deviate from the X-ray structure by RMSD 0.073 Å (max 0.162 Å; X-ray H positions are systematically short, so the largest deviations are X-H bonds).
Proton check after the optimisation (IEF-PCM water): H6–N6 1.052 Å, H6···O1 1.636 Å → **the proton stays on N6, i.e. the ionic (salt) form is retained**.
![Optimised structure with the atom numbering of the CIF](fig/structure_labels.png)

*Optimised structure with the atom numbering of the CIF*

| Hydrogen bond | D–H (Å) | H···A (Å) | D···A (Å) | D–H···A (°) |
|---|---|---|---|---|
| N6-H6...O1 (X-ray) | 0.938 | 1.870 | 2.770 | 160.0 |
| N6-H6...O1 (DFT) | 1.052 | 1.636 | 2.672 | 167.1 |

Selected bond lengths (Å), heavy atoms only:

| Bond | X-ray | DFT | Δ |
|---|---|---|---|
| C1-C2 | 1.446 | 1.459 | +0.013 |
| C1-C6 | 1.454 | 1.459 | +0.004 |
| C1-O1 | 1.251 | 1.248 | -0.004 |
| C2-C3 | 1.372 | 1.373 | +0.001 |
| C2-N1 | 1.462 | 1.460 | -0.001 |
| C3-C4 | 1.372 | 1.397 | +0.025 |
| C4-C5 | 1.377 | 1.391 | +0.014 |
| C4-N2 | 1.460 | 1.444 | -0.015 |
| C5-C6 | 1.371 | 1.380 | +0.009 |
| C6-N3 | 1.455 | 1.454 | -0.001 |
| N1-O2 | 1.218 | 1.226 | +0.008 |
| N1-O3 | 1.218 | 1.232 | +0.015 |
| N2-O4 | 1.222 | 1.234 | +0.011 |
| N2-O5 | 1.219 | 1.233 | +0.015 |
| N3-O6 | 1.225 | 1.230 | +0.005 |
| N3-O7 | 1.222 | 1.233 | +0.011 |
| C7-C8 | 1.399 | 1.401 | +0.002 |
| C7-C12 | 1.390 | 1.405 | +0.015 |
| C7-N4 | 1.365 | 1.370 | +0.005 |
| C8-C9 | 1.367 | 1.381 | +0.014 |
| C9-C10 | 1.414 | 1.420 | +0.006 |
| C10-C11 | 1.356 | 1.381 | +0.025 |
| C11-C12 | 1.395 | 1.402 | +0.007 |
| C12-N6 | 1.370 | 1.370 | +0.000 |
| N4-N5 | 1.307 | 1.312 | +0.005 |
| N5-N6 | 1.311 | 1.304 | -0.006 |

Full list including X–H bonds: `geometry_vs_start.csv`.
## 2. Frontier molecular orbitals and global reactivity descriptors
| Descriptor | Formula | Value |
|---|---|---|
| E(HOMO) |  | -6.9095 eV |
| E(LUMO) |  | -3.1743 eV |
| Energy gap | E(LUMO) - E(HOMO) | 3.7351 eV |
| Ionization potential (IP) | -E(HOMO) | 6.9095 eV |
| Electron affinity (EA) | -E(LUMO) | 3.1743 eV |
| Electronegativity (chi) | (IP+EA)/2 | 5.0419 eV |
| Chemical potential (mu) | -(IP+EA)/2 | -5.0419 eV |
| Chemical hardness (eta) | (IP-EA)/2 | 1.8676 eV |
| Chemical softness (S) | 1/(2 eta) | 0.2677 eV^-1 |
| Electrophilicity index (omega) | mu^2/(2 eta) | 6.8059 eV |
| Dipole moment |  | 23.6824 Debye |

Fragment composition of the frontier orbitals (Mulliken partition, %):

| Orbital | picrate anion | benzotriazolium cation |
|---|---|---|
| HOMO-1 | 0 | 100 |
| HOMO | 100 | 0 |
| LUMO | 100 | 0 |
| LUMO+1 | 100 | 0 |

![Frontier orbital energy levels](fig/mo_levels.png)

*Frontier orbital energy levels*

*(figure fig/homo.png not generated)*

*(figure fig/lumo.png not generated)*

*Draft interpretation.* The HOMO and LUMO energies are -6.909 and -3.174 eV, giving an energy gap of 3.735 eV. The HOMO is localised mainly on the picrate anion (100 %) whereas the LUMO resides on the picrate anion (100 %), so the HOMO→LUMO excitation is localised on that fragment. From these, the ionisation potential (6.909 eV) and electron affinity (3.174 eV) yield an electronegativity of 5.042 eV, a chemical hardness of 1.868 eV and softness of 0.268 eV⁻¹; the electrophilicity index is 6.806 eV. The moderate gap indicates kinetic stability combined with significant polarisability, consistent with the interest in this salt for optical applications. The dipole moment is 23.68 D.
## 3. Molecular electrostatic potential (MEP)
ESP on the 0.001 a.u. density surface: global minimum -52.3 kcal/mol (-0.0833 a.u.), global maximum 107.5 kcal/mol (0.1713 a.u.); surface area 343.3 Å² (positive 164.4, negative 178.8); molecular polarity index 1.21 eV; polar surface fraction 77 %.
| Extremum | ESP (kcal/mol) | ESP (a.u.) | Nearest atoms |
|---|---|---|---|
| minimum (global) | -52.3 | -0.0833 | O4 (1.90 A), O5 (2.24 A) |
| minimum | -48.0 | -0.0765 | O3 (1.88 A), O2 (2.44 A) |
| minimum | -41.6 | -0.0663 | O2 (1.86 A), O1 (2.08 A) |
| minimum | -35.5 | -0.0565 | O7 (1.83 A), N3 (2.75 A) |
| maximum (global) | 107.5 | 0.1713 | H4 (1.14 A), N4 (2.14 A) |
| maximum | 51.1 | 0.0814 | N4 (1.97 A), C7 (2.07 A) |
| maximum | 46.7 | 0.0745 | N4 (1.95 A), N5 (2.02 A) |
| maximum | 45.2 | 0.0720 | H9 (1.52 A), H10 (1.69 A) |

*(figure fig/mep.png not generated)*

*Draft interpretation.* The most negative potential (-52.3 kcal/mol) lies near O4 (1.90 A), identifying the oxygen lone pairs there as the preferred sites for electrophilic attack and as hydrogen-bond acceptors. The most positive potential (107.5 kcal/mol) surrounds H4 (1.14 A), which is therefore the nucleophilic-attack and hydrogen-bond-donor site. The strongly polarised surface (polar fraction 77 %) reflects the ionic character of the assembly and the electrostatic complementarity that drives the hydrogen bonding observed in the crystal. The molecular polarity index is 1.21 eV.
## 4. Reduced density gradient (RDG) / non-covalent interactions
Spikes in the RDG vs sign(λ₂)ρ plot (minimum RDG per bin of sign(λ₂)ρ):

| sign(λ₂)ρ (a.u.) | min RDG | interaction type |
|---|---|---|
| -0.050 | 0.07 | attractive (H-bond) |
| -0.043 | 0.30 | attractive (H-bond) |
| -0.015 | 0.03 | van der Waals |
| -0.011 | 0.11 | van der Waals |
| +0.011 | 0.04 | steric repulsion |
| +0.049 | 0.25 | steric repulsion |

![RDG versus sign(λ₂)ρ](fig/rdg_scatter.png)

*RDG versus sign(λ₂)ρ*

![RDG versus sign(λ₂)ρ, extended range showing the strong hydrogen-bond spike](fig/rdg_scatter_wide.png)

*RDG versus sign(λ₂)ρ, extended range showing the strong hydrogen-bond spike*

*(figure fig/nci.png not generated)*

Note: Multiwfn's default screening (RDG set to 100 where ρ ≥ 0.05 a.u.) would remove this very short hydrogen bond from the map, so the threshold was raised to 0.1 a.u.; covalent bond regions remain screened.

*Draft interpretation.* The RDG scatter plot shows a spike at strongly negative sign(λ₂)ρ, which corresponds in the 3-D map to the blue disc at N6-H6···O1: strong, attractive hydrogen bond. Green, low-density regions between the two partners and around the rings mark van der Waals contacts, while the red discs at the centres of the aromatic rings indicate steric (ring-closure) repulsion. Only the hydrogen bonds inside this cluster are present; the remaining hydrogen bonds of the crystal involve further symmetry-related neighbours and cannot appear here. The hydrogen-bond spike appears at sign(λ₂)ρ = -0.050 a.u. (minimum RDG 0.07); a value below about −0.05 a.u. indicates an unusually strong, partly covalent hydrogen bond, consistent with the short H···A distance.
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
