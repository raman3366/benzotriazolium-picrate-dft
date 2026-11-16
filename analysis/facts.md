# Facts dossier: benzotriazolium picrate (BTPC, paper id `btpc`)

Authoritative source of every number for the BTPC manuscript. Machine-readable twin: `facts.json` (1007 keys; each has value, unit, case, source file[:field] and a note). Both files are generated: `python scripts/analysis/build_facts.py` then `python scripts/analysis/build_facts_md.py`. Keys are quoted in backticks so an auditor can trace any number. "Derived" means computed here from a stored file (formula in the facts.json note); nothing was recomputed with quantum chemistry.

Case names: `btpc_pcm` = IEF-PCM (water) model, **primary**; `btpc_gas` = gas-phase model, comparison. Atom labels are the CIF labels; the cation labels carry an internal `#4_666` suffix in the project files, which must never appear in the paper. Never call the compound a Schiff base.

Sources that must NOT be used for numbers: `results/SUMMARY.md`, `STATUS.md` (stale, 3ATB only), `handover/` documents (context only), and the "integration grid level 3" sentence in every `RESULTS.md` (wrong; level 4 was used).

## 1. Compound and crystal structure

- Compound: benzotriazolium picrate; the CIF calls it **BTZP** (`_chemical_name_systematic` is "?"). Case-file name: "1H-benzotriazol-1-ium 2,4,6-trinitrophenolate" — see caveat C12: H atoms sit on N4 and N6 (both outer triazole N), N5 (central) carries none, i.e. the 1H,3H-benzotriazolium cation (systematic: 1H-1,2,3-benzotriazol-3-ium). Check the name against the source paper.
- Formula C12 H8 N6 O7 (moiety: C6 H2 N3 O7, C6 H6 N3), Mr 348.24 g mol⁻¹. CCDC 835950; data-block DOI 10.5517/ccx1w4y; deposited 2011-07-21. Original structure paper: doi:10.1016/j.ijleo.2012.06.018 (CIF citation year 2013).

Table 1.1 — Crystal data (all from `molecules/btpc/input/BTPC.cif`).

| Quantity | Value | key |
| --- | --- | --- |
| Crystal system | monoclinic | `xray.crystal_system` |
| Space group | P2_1/c | `xray.space_group` |
| Hall symbol | -P 2ybc | `xray.space_group_hall` |
| a (Å) | 9.2349(8) | `xray.cell_a` |
| b (Å) | 11.9272(11) | `xray.cell_b` |
| c (Å) | 12.8848(12) | `xray.cell_c` |
| α (°) | 90.00 | `xray.cell_alpha` |
| β (°) | 97.859(2) | `xray.cell_beta` |
| γ (°) | 90.00 | `xray.cell_gamma` |
| V (Å³) | 1405.9(2) | `xray.cell_volume` |
| Z / Z′ | 4 / 1 | `xray.z`, `xray.z_prime` |
| T (K) | 298 | `xray.temperature` |
| Radiation, λ (Å) | Mo Kα, 0.71073 | `xray.radiation`, `xray.wavelength` |
| D_calc (g cm⁻³) | 1.645 | `xray.density_calc` |
| μ (mm⁻¹) | 0.139 | `xray.mu` |
| F(000) | 712 | `xray.f000` |
| Crystal | yellow needle, 0.30 x 0.22 x 0.18 mm | `xray.crystal_description` |
| θ range (°) | 2.23-26.00 | `xray.theta_range` |
| Reflections measured / unique / I>2σ(I) | 14226 / 2757 / 2320 | `xray.reflections_*` |
| R_int | 0.0351 | `xray.r_int` |
| Parameters / restraints | 258 / 0 | `xray.n_parameters`, `xray.n_restraints` |
| R1 [I>2σ(I)] / R1 (all) | 0.0673 / 0.0781 | `xray.r1_gt`, `xray.r1_all` |
| wR2 [I>2σ(I)] / wR2 (all) | 0.1671 / 0.1750 | `xray.wr2_gt`, `xray.wr2_all` |
| GoF | 1.097 | `xray.gof` |
| Δρ max / min (e Å⁻³) | 0.411 / -0.485 | `xray.diff_density_*` |
| Programs | SHELXS-97 (solution), SHELXL-97 (refinement); Bruker SMART/SAINT (collection/reduction) | `xray.software` |
| Absorption correction | multi-scan, Tmin 0.9595, Tmax 0.9754 | `xray.absorption_correction` |
| H-atom treatment | mixed; every H atom has refined coordinates (with s.u.) and a refined Uiso | `xray.hydrogen_treatment` |
| Deposited H-bond table | none | `xray.hbond_table_deposited` |

**No hydrogen-bond table was deposited** (the CIF has no `_geom_hbond` loop). Table 1.2 is derived here from the deposited coordinates and the four symmetry operations (search: H···A < 2.9 Å, D–H···A > 100°, N/O acceptors; donor at x,y,z). Also written to `papers/btpc/data/derived/crystal_hbonds.csv`. X-ray H positions (N–H 0.91/0.83 Å) are short, so H···A is long and the angle uncertain; D···A is the robust number.

Table 1.2 — Crystal contacts (derived, `xray.hbNN_*`).

| # | D–H···A | A at | D–H (Å) | H···A (Å) | D···A (Å) | ∠D–H···A (°) | in DFT cluster? | comment |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 01 | N4-H4N···O1 | 4_465 (x-1, -y+1/2, z-1/2) | 0.907 | 1.775 | 2.623 | 154.5 | yes |  |
| 02 | N6-H6N···O6 | 3_666 (-x+1, -y+1, -z+1) | 0.832 | 2.103 | 2.885 | 156.5 | no |  |
| 03 | N4-H4N···O7 | 4_465 (x-1, -y+1/2, z-1/2) | 0.907 | 2.301 | 2.901 | 123.3 | yes |  |
| 04 | N6-H6N···N3 | 3_666 (-x+1, -y+1, -z+1) | 0.832 | 2.865 | 3.583 | 145.6 | no |  |
| 05 | C5-H5···O6 | 1_555 (x, y, z) | 0.808 | 2.301 | 2.637 | 105.7 | no | intramolecular S(5)-type, not an H-bond |
| 06 | C3-H3···O3 | 1_555 (x, y, z) | 0.922 | 2.318 | 2.641 | 100.2 | no | intramolecular S(5)-type, not an H-bond |
| 07 | C9-H9···O5 | 4_575 (x, -y+3/2, z-1/2) | 0.895 | 2.538 | 3.339 | 149.1 | no |  |
| 08 | C5-H5···O5 | 3_666 (-x+1, -y+1, -z+1) | 0.808 | 2.797 | 3.372 | 129.9 | no |  |
| 09 | C11-H11···O2 | 4_465 (x-1, -y+1/2, z-1/2) | 0.961 | 2.863 | 3.806 | 167.2 | yes |  |
| 10 | C11-H11···O4 | 3_565 (-x, -y+1, -z) | 0.961 | 2.875 | 3.482 | 122.0 | no |  |
| 11 | C10-H10···O3 | 3_565 (-x, -y+1, -z) | 0.959 | 2.896 | 3.732 | 146.3 | no |  |

Reading Table 1.2: the cation donates two N–H groups to two *different* picrates. N4–H4N binds the phenolate O1 of the picrate at (x−1, −y+½, z−½) with a minor bifurcated component to the ortho-nitro O7 of the same anion; N6–H6N binds nitro O6 of an inversion-related picrate at (−x+1, −y+1, −z+1). Only the first of these is inside the DFT cluster. The asymmetric unit itself is not hydrogen bonded: its closest inter-ion contact is O4...N5 2.989 Å (`xray.asu_closest_interion_*`). The only ring stacking below 4.3 Å is triazole···triazole across an inversion centre (−x, 1−y, −z): centroid–centroid 4.106 Å, interplanar separation 3.341 Å (`xray.stacking_*`).

Table 1.3 — Selected deposited geometry with s.u. (CIF `_geom_bond` / `_geom_angle`).

| Parameter | X-ray | key |
| --- | --- | --- |
| O1–C1 | 1.233(3) Å | `xray.bond_O1_C1` |
| N4–N5 | 1.324(4) Å | `xray.bond_N4_N5` |
| N6–N5 | 1.327(4) Å | `xray.bond_N6_N5` |
| N4–C12 | 1.377(3) Å | `xray.bond_N4_C12` |
| N6–C7 | 1.372(4) Å | `xray.bond_N6_C7` |
| N4–H4N | 0.91(4) Å | `xray.bond_N4_H4N` |
| N6–H6N | 0.83(4) Å | `xray.bond_N6_H6N` |
| C6–C1 | 1.450(4) Å | `xray.bond_C6_C1` |
| C1–C2 | 1.450(3) Å | `xray.bond_C1_C2` |
| N1–O2 | 1.186(3) Å | `xray.bond_N1_O2` |
| N1–O3 | 1.207(3) Å | `xray.bond_N1_O3` |
| N2–O4 | 1.229(4) Å | `xray.bond_N2_O4` |
| N2–O5 | 1.216(4) Å | `xray.bond_N2_O5` |
| N3–O6 | 1.230(3) Å | `xray.bond_N3_O6` |
| N3–O7 | 1.214(3) Å | `xray.bond_N3_O7` |
| C2–N1 | 1.456(3) Å | `xray.bond_C2_N1` |
| C4–N2 | 1.441(3) Å | `xray.bond_C4_N2` |
| N3–C6 | 1.443(3) Å | `xray.bond_N3_C6` |
| ∠C2–C1–C6 | 111.7(2)° | `xray.angle_C2_C1_C6` |
| ∠N4–N5–N6 | 108.6(3)° | `xray.angle_N4_N5_N6` |
| ∠N5–N4–C12 | 109.1(2)° | `xray.angle_N5_N4_C12` |
| ∠N5–N6–C7 | 109.7(2)° | `xray.angle_N5_N6_C7` |
| C2–C3 | 1.352(4) Å | `xray.bond_C2_C3` |
| C6–C5 | 1.367(4) Å | `xray.bond_C6_C5` |
| C7–C12 | 1.382(4) Å | `xray.bond_C7_C12` |

## 2. Structural model (how the cluster was built)

- 33 atoms, 178 electrons (89 doubly occupied orbitals), total charge 0, singlet (restricted closed shell). Fragments: picrate anion: 18 atoms (atom indices 0-17, C6H2N3O7); benzotriazolium cation: 15 atoms (indices 18-32, C6H6N3) (`model.*`).
- Built by `dft/build_cluster.py` from the CIF: the symmetry is expanded and, for each moiety, the periodic image that hydrogen-bonds to the growing cluster is selected by an exhaustive search (criteria H···A < 2.65 Å, D–H···A > 130°, heavy-atom clash < 2.40 Å). Picrate kept at 1_555 (x, y, z); benzotriazolium moved to 4_666 (x+1, -y+1/2, z+1/2). The defining contact is N4–H4N···O1 (cation N–H to phenolate O); the bifurcated N4–H4N···O7 (123°) is inside the cluster but below the builder's angle cut-off. N6–H6N is a free N–H in the model (its crystal partner, O6 of a second picrate, is absent; that picrate is the inversion image of the asymmetric-unit picrate, but relative to the cluster's own picrate it is the 2₁-screw image 2_746, see `model.missing_hbond`) — this governs the MEP maximum and the 3578/3633 cm⁻¹ band.
- Atom order everywhere is the CIF order within each ion, picrate first (atoms 0–17), then the cation (18–32); the CIF itself interleaves the two ions (O1, N3, O7, N4, C6, …). Curated by `dft/curate_cases.py` (names, ELF/LOL planes). No basis-set overrides.
- Basis size: 606 functions (C, N, O: 22 functions each (5s4p1d, spherical d); H: 7 functions (4s1p)). PySCF removed near-linear-dependent combinations of the diffuse functions: 3 at the final PCM geometry (603 MOs), 4 for gas (602 MOs); only the PCM molden was zero-padded to 606 MOs for Multiwfn (no physical effect); the gas molden holds the 602 MOs unpadded and Multiwfn read 602 orbitals. The B3LYP/6-31G(d) pre-optimisation level has 366 functions (`model.preopt_basis_functions`).

## 3. Computational protocol

Table 3.1 — Settings. Where the value comes from the PySCF 2.14.0 source rather than a run log, the key's note says so.

| Item | Setting | key |
| --- | --- | --- |
| functional | B3LYP (PySCF xc="b3lyp" = libxc HYB_GGA_XC_B3LYP, ID 402, VWN-RPA local correlation, same definition as Gaussian) | `protocol.functional` |
| basis | 6-311++G(d,p), spherical (5d) functions, all atoms | `protocol.basis` |
| density_fitting | yes (RI-JK via mf.density_fit()) | `protocol.density_fitting` |
| aux_basis | PySCF automatic default (not logged) | `protocol.aux_basis` |
| grid_level | 4 | `protocol.grid_level` |
| grid_points | H: 60 radial x 434 angular (Lebedev order 35); C, N, O: 90 radial x 590 angular (order 41); Treutler-Ahlrichs radial, Becke partition, NWChem pruning | `protocol.grid_points` |
| solvation_model | IEF-PCM | `protocol.solvation_model` |
| solvent | water | `protocol.solvent` |
| solvent_eps | 78.3553 | `protocol.solvent_eps` |
| pcm_lebedev_order | 29 | `protocol.pcm_lebedev_order` |
| pcm_radii | modified Bondi radii (H 1.10, C 1.70, N 1.55, O 1.52 Å) scaled by 1.2 -> H 1.32, C 2.04, N 1.86, O 1.824 Å | `protocol.pcm_radii` |
| pcm_discretisation | SWIG (switching-Gaussian) surface, no probe radius, electrostatic term only (no cavitation/dispersion/repulsion) | `protocol.pcm_discretisation` |
| scf_conv_tol_opt_gas | 1e-09 Ha | `protocol.scf_conv_tol_opt_gas` |
| scf_conv_tol_opt_pcm | 1e-08 Ha | `protocol.scf_conv_tol_opt_pcm` |
| scf_conv_tol_final | 1e-10 Ha | `protocol.scf_conv_tol_final` |
| scf_max_cycle | 150 | `protocol.scf_max_cycle` |
| optimizer | geomeTRIC 1.1.1 via pyscf.geomopt.geometric_solver, TRIC coordinates: 99 (= 3N) delocalised internals built from 170 primitives (34 bonds, 35 angles, 17 out-of-plane, 72 dihedrals, translation + rotation of each ion), trust radius 0.1 (max 0.3) Å, maxiter 300 | `protocol.optimizer` |
| opt_criteria | \|dE\| < 1.0e-6 Ha, RMS grad < 3.0e-4, max grad < 4.5e-4 Ha/Bohr, RMS disp < 1.2e-3, max disp < 1.8e-3 Å (all five required) | `protocol.opt_criteria` |
| preopt | B3LYP/6-31G(d) with loose criteria (dE 1e-6, grms 1.7e-3, gmax 2.5e-3, drms 6.7e-3, dmax 1.0e-2) - pipeline default | `protocol.preopt` |
| hessian_method | central finite differences of analytic gradients (6N = 198 gradients), each displaced SCF started from the converged density | `protocol.hessian_method` |
| hessian_projection | translations and rotations projected out in mass-weighted Cartesians (QR-orthonormalised projector), own routine dft/vibrations.py | `protocol.hessian_projection` |
| frequency_scaling | none (harmonic, unscaled) | `protocol.frequency_scaling` |
| thermo_T | 298.15 K | `protocol.thermo_T` |
| thermo_p | 101325 Pa | `protocol.thermo_p` |
| fmo_descriptors | IP = -E_HOMO; EA = -E_LUMO; chi = (IP+EA)/2; mu = -chi; eta = (IP-EA)/2; S = 1/(2 eta); omega = mu^2/(2 eta) | `protocol.fmo_descriptors` |
| fragment_composition_method | Mulliken partition of each MO over the two fragments: sum over fragment AOs of c_mu (S c)_mu | `protocol.fragment_composition_method` |
| cubes | HOMO/LUMO/density 80^3, ESP 64^3, margin 4.5 Bohr (2.38 Å) beyond the outermost nuclei (PySCF cubegen); render isovalues: orbitals +-0.02 a.u., MEP on 0.001 a.u. density, NCI RDG = 0.5 coloured -0.035..+0.02 a.u. | `protocol.cubes` |
| mep_surface | Multiwfn main function 12: electron-density isosurface 0.001 a.u., grid spacing 0.25 Bohr, marching tetrahedra; ESP evaluated exactly at surface vertices | `protocol.mep_surface` |
| mep_esp_code_pcm | iESPcode 2 (LIBRETA fast) | `protocol.mep_esp_code_pcm` |
| mep_esp_code_gas | iESPcode 1 (Multiwfn original ESP code) | `protocol.mep_esp_code_gas` |
| mep_note | The ESP is that of the (PCM-polarised, for btpc_pcm) molecular charge density alone; the PCM surface charges are not part of the mapped potential | `protocol.mep_note` |
| nci_settings | Multiwfn 20 -> 1 (NCI), "High-quality grid" covering the whole system with 1.5 Bohr extension; RDG_maxrho = 0.1 a.u. (points with rho >= 0.1 set to RDG 100); promolecular approximation not used (full wavefunction) | `protocol.nci_settings` |
| nci_grid_pcm | 126 x 148 x 95 = 1,771,560 points, spacing 0.158553 Bohr (0.0839 Å) | `protocol.nci_grid_pcm` |
| nci_grid_gas | 125 x 147 x 96 = 1,764,000 points, spacing 0.153832 Bohr (0.0814 Å) | `protocol.nci_grid_gas` |
| nci_spike_method | dft/mwfn.py:nci_post: sign(l2)rho binned from -0.200 to +0.050 a.u. in 0.001 bins; a bin is a candidate if it holds >20 points and its minimum RDG < 0.3; contiguous candidate bins are merged and represented by the bin with the lowest RDG; type label: < -0.015 "attractive (H-bond)", < +0.010 "van der Waals", otherwise "steric repulsion" | `protocol.nci_spike_method` |
| elf_lol_settings | Multiwfn 4 -> 9 (ELF) / 10 (LOL), plane through three atoms, 300 x 300 points, aug2D = 6.0 Bohr, Becke definition (ELFLOL_type 0, ELF_addminimal 1) | `protocol.elf_lol_settings` |
| elf_plane_picrate | picrate ring plane through C6, C4, C2 (1-based Multiwfn atoms 4, 9, 6) | `protocol.elf_plane_picrate` |
| elf_plane_cation | benzotriazolium plane through N4, N5, C7 (1-based atoms 19, 23, 21) | `protocol.elf_plane_cation` |
| host | IBM Cloud mx3d-48x480 VM (48 vCPU / 24 physical Intel Sapphire Rapids cores, 472 GB RAM, local NVMe scratch) | `protocol.host` |

Notes that the Computational-details section must reflect:

- **Grid**: level 4 (≈ Gaussian UltraFine), not level 3 as every RESULTS.md says. Verified for `btpc_gas` by re-computing the SCF at the stored geometry (`verify_grid.log`, not deposited: grid 4 reproduces the stored energy exactly, grid 3 is 1.5 µHa off). `btpc_pcm` was not in that verification run; grid 4 + Lebedev 29 for it rests only on the `run_case.sh` defaults (added because every case was launched with them). The "grid 4" in its opt.xyz comment is a literal string written by `dft/finalize_opt.py`, not a recorded setting, so it is not independent evidence.
- **PCM**: IEF-PCM, ε = 78.3553, Lebedev order 29 (302 points per sphere), modified-Bondi radii × 1.2, SWIG discretisation; electrostatics only. The dft/common.py default (order 17, the 3ATB-on-a-laptop setting) is superseded.
- **Frequencies**: central finite differences of analytic gradients (198 gradients); PCM first with step 0.005 Bohr (2026-09-17 13:56–19:42 UTC, freq_step0.005.json), then repeated with 0.01 Bohr (2026-09-18 16:05–21:51 UTC; authoritative, freq.json/hessian.npy); gas with 0.005 Bohr. For the gas case an analytic Hessian was attempted first (2026-09-16 16:31 UTC) and abandoned; the FD run replaced it. Unscaled harmonic; no IR/Raman intensities.
- **Two-stage optimisation history for BTPC**: the pipeline default is a B3LYP/6-31G(d) pre-optimisation with loose criteria, then B3LYP/6-311++G(d,p). For BTPC, both final (grid-4) optimisations started from `preopt.xyz` files carried over from earlier grid-3 runs on the (now retired) OVH VM: the gas run from the converged grid-3 geometry (3 evaluations at grid 4), the PCM run from frame 16 of the grid-3 PCM trajectory (70 steps at grid 4). The 6-31G(d) stage itself happened in those earlier runs, which are not in the archive. Reported geometries depend only on the final-level optimisation.
- **Software versions and where they were read**:
  - pyscf: 2.14.0 — data/r2/molecules/btpc/results/btpc_pcm/final.molden header "made by pyscf v[2.14.0]"; density.cube header; install log (not deposited)
  - geometric: 1.1.1 — data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log ("geomeTRIC started. Version: 1.1.1"); install log (not deposited)
  - multiwfn: 2026.9.13 (official Linux noGUI binary) — data/r2/molecules/btpc/results/btpc_pcm/logs/mwfn_mep.out banner; install log (not deposited) (Multiwfn asks to cite J. Comput. Chem. 33, 580 (2012) and J. Chem. Phys. 161, 082503 (2024); MEP module J. Mol. Graph. Model. 38, 314 (2012); LIBRETA ESP PCCP 23, 20323 (2021); MPI Carbon 171, 514 (2021); NCI original paper J. Am. Chem. Soc. 132, 6498 (2010) plus the reviews Comprehensive Computational Chemistry vol. 2, 240 (2024, DOI 10.1016/B978-0-12-821978-2.00076-3) and Angew. Chem. Int. Ed. e202504895 (2025) (as printed in the logs; the NCI entries are in logs/mwfn_nci.out))
  - paraview: 6.1.1 (pvbatch, headless) — install log (not deposited); environment/software.md (ParaView-6.1.1-MPI-Linux-Python3.12-x86_64) (The install log shows a libgomp.so.1 load error for pvbatch at install time; the renders later succeeded (logs/render.log "ok"), so it was fixed outside the log.)
  - python: 3.11.16 — environment/compute-env.yml (conda package list of the compute environment)
  - numpy_scipy_matplotlib: NumPy 2.4.6, SciPy 1.17.1, matplotlib 3.11.1 — environment/compute-env.yml (conda package list of the compute environment)
  - libxc: version not recorded — environment/compute-env.yml; install log (searched, absent; not deposited) (PySCF 2.14.0 conda-forge build; libxc version absent from the archive)
- **Multiwfn parameters in one place**: MEP on the ρ = 0.001 a.u. isosurface, surface grid spacing 0.25 Bohr; NCI with the full wavefunction, "high-quality" grid (126 x 148 x 95 = 1,771,560 points, spacing 0.158553 Bohr (0.0839 Å) for PCM; 125 x 147 x 96 = 1,764,000 points, spacing 0.153832 Bohr (0.0814 Å) for gas), 1.5 Bohr box extension, RDG_maxrho = 0.1 a.u.; NCI isosurface drawn at RDG = 0.5 coloured by sign(λ₂)ρ from −0.035 to +0.02 a.u.; ELF/LOL (Becke definitions) on 300 × 300 planes with 6 Bohr extension. Orbital isosurfaces ±0.02 a.u.

## 4. Results

### 4.1 Energies and optimisation

| Quantity | IEF-PCM (water) | gas phase | keys |
| --- | --- | --- | --- |
| E_total (Ha) | -1317.25368856 | -1317.20823094 | `*.e_total` |
| Optimisation | 70 steps (71 gradient evaluations), ≈13.4 h, 4 threads | 3 evaluations, 14.3 min | `*.opt_steps`, `*.opt_wall` |
| Converged on all five criteria? | NO — accepted on forces only (see 4.10) | yes | `*.opt_converged` |
| Final RMS / max gradient (Ha/Bohr) | 1.40e-04 / 3.85e-04 | 3.50e-05 / 6.83e-05 (log, step 2) | `pcm.opt_final_grad_*`, `gas.opt_converged` note |
| Last step ΔE (Ha); RMS/max displacement (Å) | -8.57e-06; 1.99e-02/3.76e-02 | −6.7e-07; 5.61e-04/1.18e-03 | `pcm.opt_last_step_*` |
| Starting point | grid-3 PCM trajectory frame 16 | converged grid-3 geometry | `*.opt_start` |
| Energy lowered during the final optimisation (kcal mol⁻¹) | -0.440 | — | `pcm.opt_energy_drop_total` |
| Proton transferred? | no | no | `*.proton_transferred` |

Derived: E(PCM) − E(gas) = -0.045458 Ha = -28.53 kcal mol⁻¹; G298(PCM) − G298(gas) = -27.86 kcal mol⁻¹. These mix the electrostatic solvation term with a geometry change and omit non-electrostatic terms: do not call them solvation free energies.

### 4.2 Proton position and hydrogen bonds inside the cluster

| Contact | X-ray (cluster) | IEF-PCM | gas | keys |
| --- | --- | --- | --- | --- |
| N4–H4N···O1 D–H | 0.907 | 1.053 | 1.109 | `*.hb_N4_H4N_O1_dh` |
| N4–H4N···O1 H···A | 1.775 | 1.635 | 1.437 | `*.hb_N4_H4N_O1_ha` |
| N4–H4N···O1 D···A | 2.623 | 2.672 | 2.478 | `*.hb_N4_H4N_O1_da` |
| N4–H4N···O1 ∠D–H···A | 154.5 | 167.2 | 153.2 | `*.hb_N4_H4N_O1_angle` |
| N4–H4N···O7 D–H | 0.907 | 1.053 | 1.109 | `*.hb_N4_H4N_O7_dh` |
| N4–H4N···O7 H···A | 2.301 | 2.442 | 2.268 | `*.hb_N4_H4N_O7_ha` |
| N4–H4N···O7 D···A | 2.901 | 2.983 | 2.898 | `*.hb_N4_H4N_O7_da` |
| N4–H4N···O7 ∠D–H···A | 123.3 | 110.9 | 113.7 | `*.hb_N4_H4N_O7_angle` |
| C11–H11···O2 D–H | 0.961 | 1.082 | 1.083 | `*.hb_C11_H11_O2_dh` |
| C11–H11···O2 H···A | 2.863 | 4.263 | 2.284 | `*.hb_C11_H11_O2_ha` |
| C11–H11···O2 D···A | 3.806 | 5.087 | 3.365 | `*.hb_C11_H11_O2_da` |
| C11–H11···O2 ∠D–H···A | 167.2 | 135.1 | 176.4 | `*.hb_C11_H11_O2_angle` |
| C11–H11···O1 D–H | 0.961 | 1.082 | 1.083 | `*.hb_C11_H11_O1_dh` |
| C11–H11···O1 H···A | 2.948 | 3.622 | 2.820 | `*.hb_C11_H11_O1_ha` |
| C11–H11···O1 D···A | 3.455 | 4.060 | 3.382 | `*.hb_C11_H11_O1_da` |
| C11–H11···O1 ∠D–H···A | 114.2 | 106.2 | 112.2 | `*.hb_C11_H11_O1_angle` |
| d(H4N···O1) − d(N4–H4N) (Å) | 0.868 | 0.583 | 0.328 | `*.nh_minus_ho_asym` (X-ray computed the same way) |
| O1···O7 / O1···O2 intramolecular (Å) | 2.655 / 2.606 | 2.740 / 2.753 | 2.717 / 2.763 | `*.r_O1_O7`, `*.r_O1_O2` |

Distances in Å, angles in degrees. Picrate-internal C3–H3···O3/O4 and C5–H5···O5/O6 contacts (angles 91–106°) are also listed in facts.json (`*.hb_C3_*`, `*.hb_C5_*`); they are geometric S(5)-type contacts, not hydrogen bonds. The proton stays on N4 in both models (H···O > N–H), but the gas-phase bond is much shorter (N···O 2.478 Å vs 2.672 Å PCM and 2.623 Å X-ray) with the N–H stretched to 1.109 Å. The gas model also forms a C11–H11···O2 contact (2.284 Å, 176.4°) that is absent in PCM (4.263 Å) and weak in the crystal (2.863 Å).

### 4.3 Geometry versus X-ray

| Measure | IEF-PCM | gas | keys |
| --- | --- | --- | --- |
| Bond RMSD, all 34 bonds (Å) | 0.0886 | 0.0918 | `*.bond_rmsd_all` |
| Largest bond deviation, all (Å; C5–H5) | 0.273 | 0.273 | `*.bond_max_dev_all` |
| Bond RMSD, 26 heavy-atom bonds (Å) | 0.0186 | 0.0195 | `*.bond_rmsd_heavy` |
| Bond MAE, heavy (Å) | 0.0157 | 0.0164 | `*.bond_mae_heavy` |
| Mean signed bond deviation, heavy (Å; + = DFT longer) | 0.0124 | 0.0127 | `*.bond_mean_signed_dev_heavy` |
| Angle RMSD, all 52 (°) | 1.91 | 1.82 | `*.angle_rmsd_all` |
| Angle RMSD, 36 heavy-atom angles (°) | 1.27 | 1.43 | `*.angle_rmsd_heavy` |
| Superposition RMSD, 25 heavy atoms (Å) | 0.857 | 0.460 | `*.kabsch_rmsd_heavy` |
| … picrate heavy atoms only (Å) | 0.205 | 0.219 | `*.kabsch_rmsd_picrate_heavy` |
| … cation heavy atoms only (Å) | 0.027 | 0.028 | `*.kabsch_rmsd_cation_heavy` |
| Five largest heavy-atom bond deviations | N1-O2 +0.040; C10-C9 +0.035; C4-C5 +0.031; N1-O3 +0.025; C7-C12 +0.024 | N1-O2 +0.041; C10-C9 +0.033; N4-N5 -0.030; O1-C1 +0.029; C4-C5 +0.026 | `*.bond_largest_dev_heavy` |
| Five largest heavy-atom angle deviations | C12-N4-N5 +3.10; C7-N6-N5 +2.87; N4-C12-C11 +2.41; N4-N5-N6 -2.40; O3-N1-O2 +2.36 | O3-N1-O2 +3.17; C12-N4-N5 +3.03; C7-N6-N5 +2.78; O7-N3-C6 -2.51; N6-C7-C12 -2.46 | `*.angle_largest_dev_heavy` |

The all-bond RMSD (0.089/0.092 Å, the number in RESULTS.md) is dominated by the X-ray X–H bonds, which are short by 0.12–0.27 Å as usual for X-ray H positions; quote the heavy-atom figures for the quality of the model. Each ion on its own matches the crystal well (cation 0.03 Å, picrate 0.20–0.22 Å because of nitro rotation); the whole-cluster superposition RMSD is large because the ions rotate relative to each other (Section 4.4).

Table 4.3b — Chemically important bond lengths (Å) and angles (°). X-ray from the cluster coordinates (CIF values with s.u. in Table 1.3).

| Parameter | X-ray | IEF-PCM | gas | key |
| --- | --- | --- | --- | --- |
| O1–C1 | 1.232 | 1.247 | 1.261 | `*.r_O1_C1` |
| C1–C2 | 1.449 | 1.459 | 1.447 | `*.r_C1_C2` |
| C1–C6 | 1.450 | 1.459 | 1.446 | `*.r_C1_C6` |
| C2–C3 | 1.352 | 1.373 | 1.375 | `*.r_C2_C3` |
| C3–C4 | 1.382 | 1.397 | 1.395 | `*.r_C3_C4` |
| C4–C5 | 1.361 | 1.392 | 1.387 | `*.r_C4_C5` |
| C5–C6 | 1.367 | 1.379 | 1.384 | `*.r_C5_C6` |
| C2–N1 | 1.456 | 1.460 | 1.469 | `*.r_C2_N1` |
| C4–N2 | 1.441 | 1.444 | 1.461 | `*.r_C4_N2` |
| C6–N3 | 1.443 | 1.454 | 1.461 | `*.r_C6_N3` |
| N1–O2 | 1.186 | 1.226 | 1.227 | `*.r_N1_O2` |
| N1–O3 | 1.207 | 1.232 | 1.225 | `*.r_N1_O3` |
| N2–O4 | 1.230 | 1.234 | 1.228 | `*.r_N2_O4` |
| N2–O5 | 1.216 | 1.233 | 1.227 | `*.r_N2_O5` |
| N3–O6 | 1.230 | 1.233 | 1.226 | `*.r_N3_O6` |
| N3–O7 | 1.214 | 1.230 | 1.232 | `*.r_N3_O7` |
| N4–N5 | 1.324 | 1.305 | 1.294 | `*.r_N4_N5` |
| N5–N6 | 1.328 | 1.312 | 1.326 | `*.r_N5_N6` |
| N4–C12 | 1.377 | 1.371 | 1.374 | `*.r_N4_C12` |
| N6–C7 | 1.372 | 1.370 | 1.374 | `*.r_N6_C7` |
| C7–C12 | 1.381 | 1.405 | 1.404 | `*.r_C7_C12` |
| C7–C8 | 1.384 | 1.401 | 1.400 | `*.r_C7_C8` |
| C8–C9 | 1.360 | 1.381 | 1.383 | `*.r_C8_C9` |
| C9–C10 | 1.385 | 1.420 | 1.418 | `*.r_C9_C10` |
| C10–C11 | 1.370 | 1.381 | 1.381 | `*.r_C10_C11` |
| C11–C12 | 1.379 | 1.402 | 1.402 | `*.r_C11_C12` |
| N4–H4N | 0.907 | 1.053 | 1.109 | `*.r_N4_H4N` |
| N6–H6N | 0.832 | 1.015 | 1.010 | `*.r_N6_H6N` |
| ∠C2–C1–C6 | 111.7 | 112.1 | 113.0 | `*.a_C2_C1_C6` |
| ∠N4–N5–N6 | 108.5 | 106.1 | 106.2 | `*.a_N4_N5_N6` |
| ∠N5–N4–C12 | 109.2 | 112.3 | 112.2 | `*.a_N5_N4_C12` |
| ∠N5–N6–C7 | 109.7 | 112.6 | 112.5 | `*.a_N5_N6_C7` |
| ∠C1–O1–H4N | 141.5 | 142.7 | 135.9 | `*.a_C1_O1_H4N` |
| ∠C1–O1–N4 | 150.0 | 147.7 | 144.4 | `*.a_C1_O1_N4` |

Points for the text: C1–O1 lengthens from 1.232 (X-ray) to 1.247 (PCM) and 1.261 Å (gas) as the proton approaches O1; the ipso angle C2–C1–C6 stays small (111.7/112.1/113.0°), the phenolate signature. In the cation the two N–N bonds are nearly equal in the crystal (1.324/1.328) and in PCM (1.305/1.312), but split in the gas phase (1.294/1.326): the N4 side becomes more imine-like as H4N moves towards O1.

### 4.4 Conformation: ring orientation and nitro twists

| Angle (°) | X-ray cluster | IEF-PCM | gas | keys |
| --- | --- | --- | --- | --- |
| Picrate ring / cation plane | 5.5 | 58.2 | 29.1 | `*.interplanar_picrate_cation` |
| Ring centroid distance (Å) | 6.746 | 6.936 | 6.466 | `*.centroid_distance_picrate_cation` |
| Nitro N1 twist, ortho (C2) | 10.9 | 35.9 | 40.3 | `*.nitro_twist_N1` |
| Nitro N2 twist, para (C4) | 9.3 | 0.4 | 0.6 | `*.nitro_twist_N2` |
| Nitro N3 twist, ortho (C6, next to O1) | 11.8 | 26.2 | 21.9 | `*.nitro_twist_N3` |

In the crystal the two ions are nearly coplanar (5.5°) and all three nitro groups are twisted by only 9–12°. Both isolated models keep the para-nitro group coplanar but rotate the two ortho groups by 22–40°, and open the inter-ring angle to 29° (gas) or 58° (PCM). The PCM value is not well defined (4.10).

### 4.5 Vibrational analysis and thermochemistry (unscaled harmonic)

| Quantity | IEF-PCM | gas | keys |
| --- | --- | --- | --- |
| Hessian | finite-difference (gradients, step 0.01 Bohr, 10 workers x 2 threads) | finite-difference (gradients, step 0.005 Bohr, 16 workers x 2 threads) | `*.hessian_method` |
| Hessian wall time (min) | 346 | 149 | `*.hessian_wall` |
| Imaginary modes (project's own projection / PySCF count) | 1 / 0 | 0 / 0 | `*.n_imag`, `*.n_imag_pyscf` |
| Lowest ten (cm⁻¹) | -3.22, 9.00, 24.64, 39.84, 49.95, 52.20, 55.55, 76.66, 84.79, 95.79 | 12.19, 21.54, 41.62, 52.62, 53.84, 56.62, 63.61, 69.15, 96.52, 116.60 | `*.freq_lowest10` |
| Modes below 100 cm⁻¹ (PCM count includes the imaginary mode) | 10 | 9 | `*.n_modes_below_100` |
| Highest (cm⁻¹; free N6–H6N stretch) | 3578.2 | 3632.6 | `*.freq_highest` |
| N4–H4N stretch, H-bonded (cm⁻¹) | 2889.7 | 2167.2 | `*.nu_N4_H4N_stretch` |
| N4–H4N in-plane bend (cm⁻¹) | 1561.8 | 1276.5 | `*.nu_N4_H4N_inplane_bend` |
| N4–H4N out-of-plane wag (cm⁻¹) | 959.8, 963.0 | 1172.5, 1198.3 | `*.nu_N4_H4N_oop_wag` |
| N6–H6N stretch, free (cm⁻¹) | 3578.2 | 3632.6 | `*.nu_N6_H6N_stretch` |
| N6–H6N in-plane bend / wag (cm⁻¹) | 1432.9 / 639.2, 695.4 | 1436.8 / 590.3, 689.2 | `*.nu_N6_H6N_*` |
| ν(free N6–H) − ν(bonded N4–H) (cm⁻¹) | 689 | 1465 | `diff.nh_redshift_*` |
| ZPE (Ha; kcal mol⁻¹) | 0.217329; 136.38 | 0.216747; 136.01 | `*.zpe`, `*.zpe_kcal` |
| H298 (Ha) | -1317.015341 | -1316.969856 | `*.h298` |
| G298 (Ha) | -1317.089795 | -1317.045396 | `*.g298` |
| H − E, G − E (Ha) | 0.238347, 0.163894 | 0.238375, 0.162835 | `*.thermal_corr_*` |
| S298 (cal mol⁻¹ K⁻¹) | 156.7 | 159.0 | `*.entropy` |

Mode assignments are derived from the stored Hessians (`data/r2/.../hessian.npy`): the eigenvectors were recomputed (with the masses PySCF used, H 1.008, C 12.011, N 14.007, O 15.999, the frequencies reproduce freq.json to < 1e-10 cm⁻¹), and a mode is called an N–H stretch when the H atom carries most of the mass-weighted amplitude and >50% of its motion is along the bond. Nitro stretches mix strongly with ring modes (asymmetric NO₂ character ~1495–1575 cm⁻¹, symmetric ~1286–1351 cm⁻¹ in PCM); do not assign single nitro frequencies without a PED. The three lowest modes of both models are rigid-body motions of the two ions (rigid fraction PCM 0.999, 0.998, 0.940; gas 0.991, 0.946, 0.806).

### 4.6 Frontier orbitals, descriptors, dipole, charge separation

| Descriptor (formula, dft/fmo.py) | IEF-PCM | gas | keys |
| --- | --- | --- | --- |
| E_HOMO (eV) | -6.9027 | -6.8181 | `*.homo_ev` |
| E_LUMO (eV) | -3.1672 | -3.3421 | `*.lumo_ev` |
| ΔE_gap = E_LUMO − E_HOMO (eV) | 3.7354 | 3.4760 | `*.gap_ev` |
| E_HOMO−1 (eV) | -7.5591 | -7.7375 | `*.homo_m1_ev` |
| E_LUMO+1 (eV) | -3.0201 | -2.8957 | `*.lumo_p1_ev` |
| IP = −E_HOMO (eV) | 6.9027 | 6.8181 | `*.ip_ev` |
| EA = −E_LUMO (eV) | 3.1672 | 3.3421 | `*.ea_ev` |
| χ = (IP+EA)/2 (eV) | 5.0349 | 5.0801 | `*.chi_ev` |
| μ = −(IP+EA)/2 (eV) | -5.0349 | -5.0801 | `*.mu_ev` |
| η = (IP−EA)/2 (eV) | 1.8677 | 1.7380 | `*.eta_ev` |
| S = 1/(2η) (eV⁻¹) | 0.2677 | 0.2877 | `*.softness_ev_inv` |
| ω = μ²/(2η) (eV) | 6.7866 | 7.4245 | `*.omega_ev` |
| Dipole moment (D) | 23.7898 | 16.3639 | `*.dipole_debye` |
| Dipole vector x, y, z (D) | 15.259, -16.476, 7.853 | 11.584, -10.765, 4.208 | `*.dipole_vector_debye` |
| HOMO−1 / HOMO / LUMO / LUMO+1 on picrate (%) | 0.2 / 99.9 / 99.9 / 99.9 | 97.9 / 99.4 / 0.4 / 99.6 | `*.*_picrate_pct` (Mulliken, fmo.json) |
| Highest occupied picrate / cation orbital (eV) | -6.903 / -7.559 | -6.818 / -8.401 | `*.picrate_highest_occ_ev`, `*.cation_highest_occ_ev` |
| Lowest virtual picrate / cation orbital (eV) | -3.167 / -2.571 | -2.896 / -3.342 | `*.picrate_lowest_virt_ev`, `*.cation_lowest_virt_ev` |
| Inter-ion CT gap: picrate HOMO → lowest cation virtual (eV) | 4.332 | 3.476 | `*.ct_gap_picrate_to_cation_ev` |
| Picrate-local gap (eV) | 3.735 | 3.922 | `*.local_gap_picrate_ev` |
| Mulliken charge on picrate (e) | -1.000 | -0.881 | `*.mulliken_charge_picrate` |

Orbital ladder with Mulliken fragment shares (derived; energies from fmo.json; full lists for every MO in `derived/mo_fragment_composition_<case>.csv`):

- PCM: HOMO-7 -8.977 eV (99.8% picrate); HOMO-6 -8.639 eV (100.0% picrate); HOMO-5 -8.557 eV (100.0% picrate); HOMO-4 -8.119 eV (99.7% picrate); HOMO-3 -7.950 eV (0.3% picrate); HOMO-2 -7.757 eV (98.5% picrate); HOMO-1 -7.559 eV (0.2% picrate); HOMO -6.903 eV (99.9% picrate); LUMO -3.167 eV (99.9% picrate); LUMO+1 -3.020 eV (99.9% picrate); LUMO+2 -2.672 eV (95.4% picrate); LUMO+3 -2.571 eV (4.7% picrate); LUMO+4 -0.726 eV (3.2% picrate); LUMO+5 -0.639 eV (23.7% picrate); LUMO+6 -0.612 eV (74.0% picrate); LUMO+7 -0.149 eV (39.8% picrate)
- gas: HOMO-7 -8.612 eV (99.7% picrate); HOMO-6 -8.576 eV (99.6% picrate); HOMO-5 -8.401 eV (9.2% picrate); HOMO-4 -8.335 eV (92.0% picrate); HOMO-3 -8.106 eV (100.0% picrate); HOMO-2 -7.913 eV (99.3% picrate); HOMO-1 -7.738 eV (97.9% picrate); HOMO -6.818 eV (99.4% picrate); LUMO -3.342 eV (0.4% picrate); LUMO+1 -2.896 eV (99.6% picrate); LUMO+2 -2.553 eV (99.7% picrate); LUMO+3 -2.429 eV (99.2% picrate); LUMO+4 -1.611 eV (0.4% picrate); LUMO+5 -1.439 eV (0.8% picrate); LUMO+6 -1.318 eV (-0.5% picrate); LUMO+7 -0.735 eV (-0.2% picrate)

The composition analysis rebuilds the AO overlap matrix from the molden basis; it reproduces orthonormality (max |CᵀSC − I| ≈ 3e-11) and the four fmo.json compositions to 0.05%, so the extra orbitals are on the same footing. Key reading: in the gas phase the LUMO is the cation π* (99.6% cation) and the HOMO→LUMO gap is an inter-ion charge-transfer gap; in the continuum the cation π* rises by 0.77 eV to LUMO+3 and the LUMO becomes picrate-localised, so the PCM gap (3.735 eV) is intra-picrate. The softness convention is S = 1/(2η); a paper using S = 1/η would report twice the value.

### 4.7 Molecular electrostatic potential on the 0.001 a.u. density surface

| Statistic (Multiwfn) | IEF-PCM | gas | keys |
| --- | --- | --- | --- |
| Global minimum (kcal mol⁻¹) | -52.59 | -39.65 | `*.mep_min_kcal` |
| … site | O4 | O5 | `*.mep_min_site` |
| Global maximum (kcal mol⁻¹) | 107.87 | 81.94 | `*.mep_max_kcal` |
| … site | H6N | H6N | `*.mep_max_site` |
| Global min / max (a.u.) | -0.083806 / 0.171900 | -0.063193 / 0.130579 | `*.mep_min_au`, `*.mep_max_au` |
| Minima / maxima found | 9 / 18 | 10 / 17 | `*.mep_n_minima`, `*.mep_n_maxima` |
| Surface area (Å²) | 343.91 | 328.75 | `*.mep_surface_area` |
| Positive / negative area (Å²) | 164.65 / 179.27 | 178.55 / 150.21 | `*.mep_positive_area`, `*.mep_negative_area` |
| Enclosed volume (Å³) | 360.95 | 357.26 | `*.mep_volume` |
| Mass/volume density (g cm⁻³) | 1.6020 | 1.6186 | `*.mep_density_mv` |
| Sphericity | 0.7129 | 0.7406 | `*.mep_sphericity` |
| Overall average (kcal mol⁻¹) | 2.844 | 2.767 | `*.mep_avg_all` |
| Positive / negative average (kcal mol⁻¹) | 31.801 / -23.751 | 20.798 / -18.666 | `*.mep_avg_pos`, `*.mep_avg_neg` |
| σ²_tot (kcal² mol⁻²) | 760.30 | 419.29 | `*.mep_sigma2_tot` |
| σ²₊ / σ²₋ (kcal² mol⁻²) | 543.49 / 216.80 | 303.86 / 115.43 | `*.mep_sigma2_pos`, `*.mep_sigma2_neg` |
| Balance ν | 0.2038 | 0.1995 | `*.mep_nu` |
| σ²_tot·ν (kcal² mol⁻²) | 154.98 | 83.65 | `*.mep_sigma2_nu` |
| Internal charge separation Π (kcal mol⁻¹) | 27.878 | 19.731 | `*.mep_Pi` |
| Molecular polarity index (eV; kcal mol⁻¹) | 1.1970; 27.60 | 0.8596; 19.82 | `*.mep_mpi`, `*.mep_mpi_kcal` |
| Non-polar / polar area, \|ESP\| ≤ / > 10 kcal mol⁻¹ (Å²) | 80.12 / 263.80 | 99.49 / 229.27 | `*.mep_nonpolar_area`, `*.mep_polar_area` |
| Polar area fraction (%) | 76.70 | 69.74 | `*.mep_polar_area_pct` |
| Skewness overall / + / − | 0.992 / 0.879 / 0.041 | 0.892 / 1.175 / 0.099 | `*.mep_skew_*` |
| Surface vertices | 26483 | 25344 | `*.mep_n_vertices` |
| ESP evaluation wall time (s) | 47 | 863 | `*.mep_wall` |

Table 4.7a — every surface extremum, IEF-PCM (`pcm.mep_minNN` / `pcm.mep_maxNN`; Multiwfn index; nearest atoms from opt.xyz). Note that Multiwfn lists local maxima with negative values and local minima with positive values; they are topological extrema, not sign classes.

| key | type | ESP (kcal mol⁻¹) | ESP (a.u.) | nearest atoms (Å) | ion |
| --- | --- | --- | --- | --- | --- |
| `pcm.mep_min02` | minimum (global) | -52.59 | -0.083806 | O4 (1.94), O5 (2.10) | picrate |
| `pcm.mep_min01` | minimum | -46.72 | -0.074460 | O3 (1.87), N1 (2.54) | picrate |
| `pcm.mep_min08` | minimum | -37.36 | -0.059541 | O6 (1.84), N3 (2.74) | picrate |
| `pcm.mep_min04` | minimum | -35.77 | -0.057006 | O2 (1.89), O1 (2.03) | picrate |
| `pcm.mep_max09` | maximum | -28.45 | -0.045343 | O2 (1.74), H11 (2.71) | picrate |
| `pcm.mep_min06` | minimum | -6.30 | -0.010047 | N5 (1.92), N4 (2.71) | cation |
| `pcm.mep_min03` | minimum | -0.78 | -0.001239 | C4 (2.04), C3 (2.31) | picrate |
| `pcm.mep_max01` | maximum | 0.64 | 0.001023 | N2 (1.85), C4 (2.11) | picrate |
| `pcm.mep_min05` | minimum | 0.91 | 0.001445 | C6 (2.09), C5 (2.26) | picrate |
| `pcm.mep_max03` | maximum | 1.19 | 0.001900 | N1 (1.86), C2 (2.05) | picrate |
| `pcm.mep_max06` | maximum | 3.27 | 0.005206 | N2 (1.87), C4 (2.08) | picrate |
| `pcm.mep_max02` | maximum | 3.56 | 0.005681 | H3 (1.34), C3 (2.24) | picrate |
| `pcm.mep_max04` | maximum | 3.68 | 0.005857 | C1 (2.16), C6 (2.18) | picrate |
| `pcm.mep_max08` | maximum | 4.67 | 0.007440 | C6 (2.22), C2 (2.26) | picrate |
| `pcm.mep_max05` | maximum | 7.74 | 0.012340 | N1 (1.89), H3 (1.97) | picrate |
| `pcm.mep_max13` | maximum | 11.33 | 0.018055 | N3 (1.85), C6 (2.15) | picrate |
| `pcm.mep_min07` | minimum | 14.43 | 0.022993 | H11 (1.41), C11 (2.16) | cation |
| `pcm.mep_max07` | maximum | 14.64 | 0.023331 | N3 (1.81), C6 (2.13) | picrate |
| `pcm.mep_max10` | maximum | 14.68 | 0.023391 | H4N (1.52), N5 (1.95) | cation |
| `pcm.mep_min09` | minimum | 17.32 | 0.027601 | C11 (2.07), C10 (2.29) | cation |
| `pcm.mep_max11` | maximum | 24.03 | 0.038287 | H11 (1.81), H4N (1.81) | cation |
| `pcm.mep_max14` | maximum | 30.56 | 0.048703 | H4N (1.66), H11 (1.72) | cation |
| `pcm.mep_max16` | maximum | 36.54 | 0.058236 | H10 (1.52), H11 (1.70) | cation |
| `pcm.mep_max18` | maximum | 44.88 | 0.071521 | H9 (1.50), H10 (1.73) | cation |
| `pcm.mep_max17` | maximum | 46.36 | 0.073884 | N6 (1.95), N5 (2.04) | cation |
| `pcm.mep_max12` | maximum | 52.47 | 0.083608 | N6 (1.98), N5 (2.02) | cation |
| `pcm.mep_max15` | maximum (global) | 107.87 | 0.171900 | H6N (1.14), N6 (2.13) | cation |

Table 4.7b — every surface extremum, gas phase (`gas.mep_minNN` / `gas.mep_maxNN`; Multiwfn index; nearest atoms from opt.xyz). Note that Multiwfn lists local maxima with negative values and local minima with positive values; they are topological extrema, not sign classes.

| key | type | ESP (kcal mol⁻¹) | ESP (a.u.) | nearest atoms (Å) | ion |
| --- | --- | --- | --- | --- | --- |
| `gas.mep_min03` | minimum (global) | -39.65 | -0.063193 | O5 (1.84), N2 (2.53) | picrate |
| `gas.mep_min02` | minimum | -39.61 | -0.063129 | O4 (1.84), N2 (2.52) | picrate |
| `gas.mep_min09` | minimum | -33.89 | -0.054015 | O6 (1.83), N3 (2.62) | picrate |
| `gas.mep_min05` | minimum | -32.72 | -0.052147 | O3 (1.82), N1 (2.67) | picrate |
| `gas.mep_min01` | minimum | -26.95 | -0.042942 | O3 (1.74), N1 (2.86) | picrate |
| `gas.mep_min06` | minimum | -20.26 | -0.032283 | O1 (1.85), H4N (2.61) | picrate |
| `gas.mep_min04` | minimum | -0.55 | -0.000876 | C4 (2.06), C5 (2.30) | picrate |
| `gas.mep_min07` | minimum | 2.00 | 0.003192 | C4 (2.07), C5 (2.26) | picrate |
| `gas.mep_min08` | minimum | 2.51 | 0.003994 | C6 (2.12), C5 (2.21) | picrate |
| `gas.mep_max04` | maximum | 4.77 | 0.007594 | C1 (1.95), C6 (2.29) | picrate |
| `gas.mep_max03` | maximum | 5.34 | 0.008506 | C6 (2.22), C1 (2.24) | picrate |
| `gas.mep_max10` | maximum | 6.23 | 0.009922 | C6 (2.23), C2 (2.26) | picrate |
| `gas.mep_max12` | maximum | 7.02 | 0.011180 | C1 (1.91), O1 (2.26) | picrate |
| `gas.mep_max09` | maximum | 7.21 | 0.011492 | H5 (1.31), C5 (2.33) | picrate |
| `gas.mep_max05` | maximum | 7.56 | 0.012042 | N1 (1.84), C2 (2.11) | picrate |
| `gas.mep_max01` | maximum | 7.72 | 0.012298 | N2 (1.85), C4 (2.11) | picrate |
| `gas.mep_max07` | maximum | 9.25 | 0.014748 | N2 (1.85), C4 (2.12) | picrate |
| `gas.mep_min10` | minimum | 9.70 | 0.015457 | C11 (2.08), C10 (2.21) | cation |
| `gas.mep_max02` | maximum | 11.04 | 0.017588 | H3 (1.32), C3 (2.31) | picrate |
| `gas.mep_max08` | maximum | 12.52 | 0.019955 | N3 (1.83), C6 (2.11) | picrate |
| `gas.mep_max13` | maximum | 12.85 | 0.020477 | N3 (1.83), C6 (2.17) | picrate |
| `gas.mep_max11` | maximum | 14.14 | 0.022531 | H4N (1.55), N4 (1.87) | cation |
| `gas.mep_max06` | maximum | 15.32 | 0.024414 | N1 (1.88), H3 (1.95) | picrate |
| `gas.mep_max15` | maximum | 37.21 | 0.059290 | N5 (2.03), N6 (2.03) | cation |
| `gas.mep_max17` | maximum | 37.36 | 0.059545 | H9 (1.33), C9 (2.40) | cation |
| `gas.mep_max14` | maximum | 38.83 | 0.061879 | N6 (2.00), N5 (2.04) | cation |
| `gas.mep_max16` | maximum (global) | 81.94 | 0.130579 | H6N (1.19), N6 (2.18) | cation |

Physical reading: the most negative potential sits on the para-nitro oxygens O4/O5, not on the phenolate O1, because O1 is capped by the N4–H4N hydrogen bond; the most positive region is the free N6–H6N. The ESP maps the charge density only (PCM surface charges are not included in the mapped potential).

### 4.8 NCI (RDG) analysis

Table 4.8a — RDG spikes as reported by the pipeline (`pcm.nci_spikeN`, from nci_spikes.csv) with the atoms at their lowest-RDG points (derived from rdg_scatter.txt), IEF-PCM.

| key | sign(λ₂)ρ bin (a.u.) | min RDG | points in bin | pipeline label | lowest-RDG points between |
| --- | --- | --- | --- | --- | --- |
| `pcm.nci_spike1` | -0.0495 | 0.092 | 780 | attractive (H-bond) | H4N/O1 (inter-ion) |
| `pcm.nci_spike2` | -0.0145 | 0.052 | 4095 | van der Waals | O2/O1, O7/O1 |
| `pcm.nci_spike3` | -0.0095 | 0.186 | 6707 | van der Waals | H4N/O7 (inter-ion) |
| `pcm.nci_spike4` | +0.0095 | 0.089 | 5078 | van der Waals | H4N/O7 (inter-ion) |
| `pcm.nci_spike5` | +0.0135 | 0.050 | 3574 | steric repulsion | O2/O1 |
| `pcm.nci_spike6` | +0.0475 | 0.286 | 876 | steric repulsion | N6/N4 |

Table 4.8b — RDG spikes as reported by the pipeline (`gas.nci_spikeN`, from nci_spikes.csv) with the atoms at their lowest-RDG points (derived from rdg_scatter.txt), gas phase.

| key | sign(λ₂)ρ bin (a.u.) | min RDG | points in bin | pipeline label | lowest-RDG points between |
| --- | --- | --- | --- | --- | --- |
| `gas.nci_spike1` | -0.0865 | 0.084 | 391 | attractive (H-bond) | H4N/O1 (inter-ion) |
| `gas.nci_spike2` | -0.0825 | 0.180 | 409 | attractive (H-bond) | H4N/O1 (inter-ion) |
| `gas.nci_spike3` | -0.0805 | 0.223 | 424 | attractive (H-bond) | H4N/O1 (inter-ion) |
| `gas.nci_spike4` | -0.0165 | 0.032 | 3696 | attractive (H-bond) | H4N/O7 (inter-ion) |
| `gas.nci_spike5` | -0.0105 | 0.112 | 6697 | van der Waals | H11/O2 (inter-ion) |
| `gas.nci_spike6` | -0.0055 | 0.097 | 12638 | van der Waals | H11/O1 (inter-ion) |
| `gas.nci_spike7` | +0.0045 | 0.128 | 13691 | van der Waals | H11/O1 (inter-ion), H11/H4N |
| `gas.nci_spike8` | +0.0135 | 0.039 | 3870 | steric repulsion | O2/O1 |
| `gas.nci_spike9` | +0.0495 | 0.253 | 986 | steric repulsion | N4/C12, N6/N4 |

Table 4.8c — Real-space NCI regions (derived: connected RDG < 0.5 regions of rdg.cube with ≥20 grid points; value = sign(λ₂)ρ at the lowest-RDG point), IEF-PCM.

| key | sign(λ₂)ρ at min RDG (a.u.) | min RDG | sign(λ₂)ρ range | points | located between |
| --- | --- | --- | --- | --- | --- |
| `pcm.nci_region01` | -0.0496 | 0.092 | -0.0537..-0.0357 | 78 | H4N and O1 (inter-ion) |
| `pcm.nci_region02` | -0.0138 | 0.215 | -0.0147..+0.0192 | 95 | H3 and O3 (within one ion) |
| `pcm.nci_region03` | +0.0137 | 0.050 | -0.0148..+0.0155 | 177 | O2 and O1 (within one ion) |
| `pcm.nci_region04` | +0.0142 | 0.088 | -0.0159..+0.0163 | 348 | O7 and O1 (within one ion) |
| `pcm.nci_region05` | +0.0157 | 0.152 | -0.0160..+0.0205 | 112 | H5 and O5 (within one ion) |
| `pcm.nci_region06` | +0.0158 | 0.115 | -0.0168..+0.0203 | 109 | H3 and O4 (within one ion) |
| `pcm.nci_region07` | +0.0170 | 0.148 | -0.0165..+0.0203 | 116 | H5 and O6 (within one ion) |
| `pcm.nci_region08` | +0.0199 | 0.094 | +0.0146..+0.0221 | 59 | C6 and C2 (within one ion) |
| `pcm.nci_region09` | +0.0217 | 0.156 | +0.0149..+0.0242 | 60 | C7 and C12 (within one ion) |
| `pcm.nci_region10` | +0.0554 | 0.117 | +0.0365..+0.0732 | 101 | N6 and N4 (within one ion) |

Table 4.8d — Real-space NCI regions (derived: connected RDG < 0.5 regions of rdg.cube with ≥20 grid points; value = sign(λ₂)ρ at the lowest-RDG point), gas phase.

| key | sign(λ₂)ρ at min RDG (a.u.) | min RDG | sign(λ₂)ρ range | points | located between |
| --- | --- | --- | --- | --- | --- |
| `gas.nci_region01` | -0.0863 | 0.084 | -0.0988..+0.0591 | 114 | H4N and O1 (inter-ion) |
| `gas.nci_region02` | -0.0165 | 0.032 | -0.0175..+0.0172 | 437 | H4N and O7 (inter-ion) |
| `gas.nci_region03` | -0.0155 | 0.148 | -0.0165..+0.0206 | 129 | H3 and O4 (within one ion) |
| `gas.nci_region04` | -0.0105 | 0.112 | -0.0110..-0.0085 | 58 | H11 and O2 (inter-ion) |
| `gas.nci_region05` | -0.0052 | 0.097 | -0.0053..+0.0053 | 138 | H11 and O1 (inter-ion) |
| `gas.nci_region06` | +0.0136 | 0.039 | -0.0145..+0.0154 | 197 | O2 and O1 (within one ion) |
| `gas.nci_region07` | +0.0140 | 0.247 | -0.0132..+0.0178 | 90 | H3 and O3 (within one ion) |
| `gas.nci_region08` | +0.0161 | 0.138 | -0.0165..+0.0205 | 131 | H5 and O5 (within one ion) |
| `gas.nci_region09` | +0.0171 | 0.123 | -0.0179..+0.0214 | 133 | H5 and O6 (within one ion) |
| `gas.nci_region10` | +0.0204 | 0.072 | +0.0144..+0.0229 | 70 | C6 and C2 (within one ion) |
| `gas.nci_region11` | +0.0215 | 0.107 | +0.0154..+0.0239 | 64 | C12 and C7 (within one ion) |
| `gas.nci_region12` | +0.0542 | 0.095 | +0.0350..+0.0930 | 117 | N6 and N4 (within one ion) |

Reading: the only strongly attractive feature is the N4–H4N···O1 bond, at sign(λ₂)ρ ≈ −0.050 a.u. (PCM) and −0.086 a.u. (gas; the region reaches −0.099 a.u., close to the RDG_maxrho = 0.1 screening limit, so the gas isosurface of this bond may be slightly truncated). Weak attractive/van der Waals regions: H4N···O7 (bifurcated), and in the gas phase only, C11–H11···O2 and H11···O1. Positive (repulsive/ring) regions: the three ring centres (+0.020 to +0.055 a.u.), O1···O2 and O1···O7 between the phenolate and ortho-nitro oxygens, and the picrate C–H···O(nitro) five-membered pockets. The pipeline's spike labels are a pure threshold on sign(λ₂)ρ (e.g. the PCM −0.0145 "van der Waals" spike is the intramolecular O1···O7/O1···O2 contact region); use the located-between column, not the label.

### 4.9 ELF and LOL maps

| Quantity | IEF-PCM | gas | keys |
| --- | --- | --- | --- |
| ELF range, picrate plane (C6, C4, C2) | 0.0000..1.0000 | 0.0000..1.0000 | `*.elf_picrate_ring_range` |
| ELF range, cation plane (N4, N5, C7) | 0.0000..1.0000 | 0.0000..1.0000 | `*.elf_benzotriazolium_ring_range` |
| ELF at picrate centroid | 0.022 | 0.022 | `*.elf_picrate_centroid` |
| ELF at benzo centroid | 0.022 | 0.022 | `*.elf_benzo_centroid` |
| ELF at triazole centroid | 0.048 | 0.049 | `*.elf_triazole_centroid` |
| LOL range, picrate plane (C6, C4, C2) | 0.0011..0.9664 | 0.0011..0.9694 | `*.lol_picrate_ring_range` |
| LOL range, cation plane (N4, N5, C7) | 0.0002..0.9788 | 0.0002..0.9816 | `*.lol_benzotriazolium_ring_range` |
| LOL at picrate centroid | 0.129 | 0.130 | `*.lol_picrate_centroid` |
| LOL at benzo centroid | 0.131 | 0.131 | `*.lol_benzo_centroid` |
| LOL at triazole centroid | 0.183 | 0.184 | `*.lol_triazole_centroid` |
| O1 distance from the cation map plane (Å; H4N within 0.01 Å) | -0.362 | -0.153 | `*.cation_plane_offset_O1` |

Maps are 300 × 300 points in the plane through three atoms; the cation plane contains the whole benzotriazolium (triazole and benzo rings coplanar) and, nearly, the H4N···O1 bond, so the cation map is the one to show the hydrogen bond. ELF/LOL values at ring centroids are the usual aromatic depletion (ELF ≈ 0.02–0.05). No quantitative ELF value along the H-bond is read from the plane maps because O1 lies out of the plane (0.36 Å in PCM). Superseded for the value on the H4N–O1 line itself: the revision-1 keys `*.elf_min_H4N_O1` and `*.bcp_H4N_O1_elf` give ELF along that line, which the plane offset does not affect.

### 4.10 The PCM imaginary mode: what was checked

- The authoritative PCM Hessian (step 0.01 Bohr) gives one imaginary mode, -3.22 cm⁻¹; the earlier Hessian with step 0.005 Bohr gave -2.48 cm⁻¹. Doubling the step changed no frequency by more than 1.89 cm⁻¹ (mean 0.15), ZPE by -0.0059 and G298 by -0.046 kcal mol⁻¹. A first attempt (19 workers) failed with exit code 1 and produced nothing (`saddle.first_attempt`).
- PySCF's own harmonic analysis stored the mode as 0 + 3.2i and therefore counted zero imaginary modes (`pcm.n_imag_pyscf` = 0); the project's explicit translation/rotation projection (`dft/vibrations.py`) reports it correctly.
- Character: 99.9% of the eigenvector is rigid-body motion of the two ions relative to each other; 66% of the mass-weighted amplitude is on the cation (`saddle.imag_mode_*`). It is an inter-ion libration, not a distortion of either ion.
- The optimiser was still moving along exactly this coordinate. The net displacement from step 0 to 70 projects 87% onto the imaginary mode, and from step 40 to 70 94% (max atomic displacement 1.8 Å for 0.30 kcal mol⁻¹); the last ten steps moved along the 9.0 cm⁻¹ libration (95%).
- Along the whole trajectory the hydrogen bond hardly changed while the rings rotated:

| step | E (Ha) | N4–H4N | H4N···O1 | N4···O1 | ∠ (°) | inter-ring angle (°) | H11···O2 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | -1317.25298705 | 1.048 | 1.678 | 2.663 | 154.5 | 5.8 | 3.005 |
| 25 | -1317.25304527 | 1.049 | 1.671 | 2.676 | 158.7 | 7.3 | 3.102 |
| 30 | -1317.25313663 | 1.049 | 1.666 | 2.675 | 159.7 | 18.1 | 3.163 |
| 40 | -1317.25320810 | 1.050 | 1.658 | 2.675 | 161.7 | 20.6 | 3.288 |
| 45 | -1317.25345241 | 1.052 | 1.647 | 2.678 | 165.1 | 32.9 | 3.610 |
| 50 | -1317.25361024 | 1.053 | 1.637 | 2.676 | 168.1 | 47.7 | 4.113 |
| 60 | -1317.25365541 | 1.053 | 1.633 | 2.671 | 167.6 | 55.9 | 4.283 |
| 70 | -1317.25368856 | 1.053 | 1.635 | 2.672 | 167.2 | 58.2 | 4.263 |

- Conclusion for the paper: The PCM stationary-point search ended on an almost flat inter-ion libration coordinate: over the last 30 steps the ions rotated by ~38 deg for 0.30 kcal/mol, the residual imaginary frequency (3.2i cm-1) belongs to that same rigid-body motion and survived a doubled FD step. The intra-ion geometry and the N4-H4N...O1 bond are converged; the relative orientation of the two rings is not determined to better than a few tens of degrees.
- The gas-phase model has no imaginary mode (lowest 12.2 cm⁻¹) and met every convergence criterion.

## 5. Gas phase versus IEF-PCM (derived: PCM − gas unless stated)

| Quantity | PCM − gas | key |
| --- | --- | --- |
| e_total | -28.53 kcal mol-1 | `diff.e_total` |
| e_total_ha | -0.045458 Ha | `diff.e_total_ha` |
| g298 | -27.86 kcal mol-1 | `diff.g298` |
| zpe | +0.3657 kcal mol-1 | `diff.zpe` |
| r_N4_H4N | -0.0568 Å | `diff.r_N4_H4N` |
| hb_N4_H4N_O1_ha | +0.1983 Å | `diff.hb_N4_H4N_O1_ha` |
| hb_N4_H4N_O1_da | +0.1936 Å | `diff.hb_N4_H4N_O1_da` |
| hb_N4_H4N_O1_angle | +13.95 deg | `diff.hb_N4_H4N_O1_angle` |
| hb_N4_H4N_O7_ha | +0.1731 Å | `diff.hb_N4_H4N_O7_ha` |
| hb_C11_H11_O2_ha | +1.98 Å | `diff.hb_C11_H11_O2_ha` |
| r_O1_C1 | -0.0140 Å | `diff.r_O1_C1` |
| r_N4_N5 | +0.0105 Å | `diff.r_N4_N5` |
| r_N5_N6 | -0.0136 Å | `diff.r_N5_N6` |
| r_N6_H6N | +0.0050 Å | `diff.r_N6_H6N` |
| r_C1_C2 | +0.0123 Å | `diff.r_C1_C2` |
| r_C1_C6 | +0.0123 Å | `diff.r_C1_C6` |
| a_C2_C1_C6 | -0.8355 deg | `diff.a_C2_C1_C6` |
| interplanar_picrate_cation | +29.03 deg | `diff.interplanar_picrate_cation` |
| centroid_distance_picrate_cation | +0.4697 Å | `diff.centroid_distance_picrate_cation` |
| nitro_twist_N1 | -4.35 deg | `diff.nitro_twist_N1` |
| nitro_twist_N2 | -0.2764 deg | `diff.nitro_twist_N2` |
| nitro_twist_N3 | +4.30 deg | `diff.nitro_twist_N3` |
| nh_minus_ho_asym | +0.2551 Å | `diff.nh_minus_ho_asym` |
| bond_rmsd_all | -0.0032 Å | `diff.bond_rmsd_all` |
| bond_rmsd_heavy | -0.0009 Å | `diff.bond_rmsd_heavy` |
| angle_rmsd_heavy | -0.1556 deg | `diff.angle_rmsd_heavy` |
| kabsch_rmsd_heavy | +0.3976 Å | `diff.kabsch_rmsd_heavy` |
| homo_ev | -0.0846 eV | `diff.homo_ev` |
| lumo_ev | +0.1748 eV | `diff.lumo_ev` |
| gap_ev | +0.2594 eV | `diff.gap_ev` |
| homo_m1_ev | +0.1784 eV | `diff.homo_m1_ev` |
| lumo_p1_ev | -0.1245 eV | `diff.lumo_p1_ev` |
| ip_ev | +0.0846 eV | `diff.ip_ev` |
| ea_ev | -0.1748 eV | `diff.ea_ev` |
| chi_ev | -0.0451 eV | `diff.chi_ev` |
| mu_ev | +0.0451 eV | `diff.mu_ev` |
| eta_ev | +0.1297 eV | `diff.eta_ev` |
| softness_ev_inv | -0.0200 eV-1 | `diff.softness_ev_inv` |
| omega_ev | -0.6379 eV | `diff.omega_ev` |
| dipole_debye | +7.43 D | `diff.dipole_debye` |
| cation_highest_occ_ev | +0.8418 eV | `diff.cation_highest_occ_ev` |
| cation_lowest_virt_ev | +0.7711 eV | `diff.cation_lowest_virt_ev` |
| picrate_highest_occ_ev | -0.0846 eV | `diff.picrate_highest_occ_ev` |
| picrate_lowest_virt_ev | -0.2716 eV | `diff.picrate_lowest_virt_ev` |
| ct_gap_picrate_to_cation_ev | +0.8557 eV | `diff.ct_gap_picrate_to_cation_ev` |
| local_gap_picrate_ev | -0.1870 eV | `diff.local_gap_picrate_ev` |
| mulliken_charge_picrate | -0.1195 e | `diff.mulliken_charge_picrate` |
| mep_min_kcal | -12.93 kcal mol-1 | `diff.mep_min_kcal` |
| mep_max_kcal | +25.93 kcal mol-1 | `diff.mep_max_kcal` |
| mep_surface_area | +15.16 Å^2 | `diff.mep_surface_area` |
| mep_positive_area | -13.90 Å^2 | `diff.mep_positive_area` |
| mep_negative_area | +29.06 Å^2 | `diff.mep_negative_area` |
| mep_volume | +3.69 Å^3 | `diff.mep_volume` |
| mep_avg_all | +0.0765 kcal mol-1 | `diff.mep_avg_all` |
| mep_avg_pos | +11.00 kcal mol-1 | `diff.mep_avg_pos` |
| mep_avg_neg | -5.09 kcal mol-1 | `diff.mep_avg_neg` |
| mep_sigma2_tot | +341.01 (kcal mol-1)^2 | `diff.mep_sigma2_tot` |
| mep_sigma2_pos | +239.64 (kcal mol-1)^2 | `diff.mep_sigma2_pos` |
| mep_sigma2_neg | +101.37 (kcal mol-1)^2 | `diff.mep_sigma2_neg` |
| mep_nu | +0.0043 | `diff.mep_nu` |
| mep_sigma2_nu | +71.33 (kcal mol-1)^2 | `diff.mep_sigma2_nu` |
| mep_Pi | +8.15 kcal mol-1 | `diff.mep_Pi` |
| mep_mpi | +0.3374 eV | `diff.mep_mpi` |
| mep_polar_area_pct | +6.96 % | `diff.mep_polar_area_pct` |
| nu_N4_H4N_stretch | +722.50 cm-1 | `diff.nu_N4_H4N_stretch` |
| nu_N6_H6N_stretch | -54.45 cm-1 | `diff.nu_N6_H6N_stretch` |
| freq_lowest | -15.41 cm-1 | `diff.freq_lowest` |
| freq_highest | -54.45 cm-1 | `diff.freq_highest` |
| h298 | -0.045486 Ha | `diff.h298` |
| entropy | -2.29 cal mol-1 K-1 | `diff.entropy` |
| nci_hbond_spike | +0.0370 a.u. | `diff.nci_hbond_spike` |
| nh_redshift_pcm | +688.52 cm-1 | `diff.nh_redshift_pcm` |
| nh_redshift_gas | +1465.47 cm-1 | `diff.nh_redshift_gas` |
| xray_vs_pcm_n_o | +0.0484 Å | `diff.xray_vs_pcm_n_o` |
| xray_vs_gas_n_o | -0.1452 Å | `diff.xray_vs_gas_n_o` |
| ch_stretch_max_shift | +11.66 cm-1 | `diff.ch_stretch_max_shift` |
| lowdin_charge_picrate | -0.1314 e | `diff.lowdin_charge_picrate` |
| mulliken_charge_H4N | -0.0040 e | `diff.mulliken_charge_H4N` |
| lowdin_charge_H4N | +0.0100 e | `diff.lowdin_charge_H4N` |
| bcp_H_ratio_gas_pcm | +4.94 | `diff.bcp_H_ratio_gas_pcm` |

Energy differences mix solvation and geometry (see 4.1). Orbital-energy shifts are Koopmans quantities in two different electrostatic environments.

## 6. Storyline candidates (ranked by strength of the evidence)

1. **The continuum switches the frontier orbitals from inter-ion charge transfer to intra-picrate.** Gas: HOMO on picrate (99.4%), LUMO on the cation (99.6%), gap 3.476 eV = CT gap. PCM: cation π* rises by 0.77 eV (to LUMO+3), cation HOMO rises by 0.84 eV (to HOMO−1), the picrate LUMO drops by 0.27 eV; LUMO and LUMO+1 are both picrate, gap 3.735 eV, and the CT gap opens to 4.332 eV. So "gap" means different excitations in the two models; the descriptors (ω 7.42 → 6.79 eV, η 1.74 → 1.87 eV) inherit that. Solid: validated orbital compositions.
2. **Environment-controlled proton position in the N⁺–H···O⁻ bond, salt retained in both models.** N4···O1 2.478 Å (gas) vs 2.672 Å (PCM) vs 2.623 Å (X-ray); N–H 1.109 vs 1.053 Å; ν(N4–H) 2167 vs 2890 cm⁻¹ (red shifts of 1465 and 689 cm⁻¹ against the free N6–H); NCI sign(λ₂)ρ -0.086 vs -0.050 a.u.; Mulliken charge on picrate -0.88 vs -1.00 e; C1–O1 1.261 vs 1.247 Å; N4–N5/N5–N6 asymmetry 0.032 vs 0.008 Å. Unlike 3ATB (proton migrated in the gas phase), BTPC stays ionic even in vacuo, but the gas-phase bond is a short, strong one (N···O < 2.5 Å; "partly shared" was the early wording, superseded by the critical-point data of revision 1: both bonds lie in the transit region, the gas-phase one far deeper). PCM reproduces the crystal N···O to +0.05 Å; the gas phase under-shoots by 0.15 Å.
3. **A floppy ion pair: the hydrogen bond is robust, the relative orientation is not.** Along 70 PCM steps the H-bond metrics changed by ≤0.05 Å while the inter-ring angle went from 5.8° to 58.2° for -0.44 kcal mol⁻¹; the residual -3.2 cm⁻¹ mode is 99.9% rigid inter-ion libration and survived a doubled FD step. This is a methodological point worth making honestly (crystal packing, i.e. the second N–H···O and the stacking, fixes the orientation that a single ion pair cannot).
4. **Packing flattens the picrate.** X-ray nitro twists 9–12° and a nearly coplanar pair (5.5°); both isolated models twist the ortho nitro groups to 22–40° (N1 35.9°, N3 26.2° PCM) while keeping the para group planar (0.4°). The largest bond errors are the nitro N1–O2 (+0.04 Å) where X-ray bonds are anomalously short (1.186 Å; high Ueq of O2/O3). Heavy-atom bond RMSD is only 0.019 Å.
5. **Electrostatic surface.** Most negative ESP on the para-nitro O4/O5 (-52.6 kcal mol⁻¹ PCM; -39.7 gas), not on the H-bonded phenolate; most positive on the free N6–H6N (107.9; 81.9) — precisely the site that in the crystal donates to O6 of the next picrate. The continuum amplifies the charge separation (σ²_tot 419 → 760 kcal² mol⁻², MPI 0.86 → 1.20 eV, dipole 16.4 → 23.8 D).

## 7. Caveats and limitations the paper must state

- C1. **PCM geometry is not a fully converged minimum.** Accepted on forces (rms 1.40e-4, max 3.85e-4 Ha/Bohr, inside Gaussian thresholds) but the energy (ΔE −8.6e-6 Ha) and displacement (0.02/0.04 Å) criteria were not met and the energy was still falling; one imaginary frequency of 3.2i cm⁻¹ (2.5i with the smaller FD step) remains. It is a rigid inter-ion libration; internal geometry and the H-bond are stable. Do not write "true minimum" or "all real frequencies" for the PCM model.
- C2. **Relative ion orientation in PCM is poorly defined** (inter-ring angle 58° vs 5.5° in the crystal, 29° in gas; drifted by ~38° for 0.3 kcal/mol). Any property sensitive to the orientation (dipole magnitude, MEP surface statistics, NCI C–H···O contacts) inherits this uncertainty.
- C3. **Cluster, not crystal.** One ion pair built by symmetry; the second N–H···O bond (N6–H6N···O6), C9–H9···O5, the triazole stacking and long-range lattice electrostatics are absent. N6–H6N is artificially free (MEP maximum, 3578 cm⁻¹ band).
- C4. **Continuum solvent (water) stands in for the crystal field**; electrostatic IEF-PCM only, no cavitation/dispersion terms, equilibrium solvation. E(PCM) − E(gas) is not a solvation free energy.
- C5. **No dispersion correction** (plain B3LYP), which matters for the inter-ring orientation and any stacking.
- C6. **Frequencies are harmonic, unscaled, from finite differences; no IR/Raman intensities**, so no simulated spectrum can be compared with experiment. The strongly H-bonded N4–H stretch is highly anharmonic in reality.
- C7. **Thermochemistry**: RRHO at 298.15 K/1 atm with ~10 modes below 100 cm⁻¹ treated as harmonic oscillators; entropies and G are soft.
- C8. **Koopmans descriptors** (IP = −E_HOMO etc.) from orbital energies; in PCM they refer to a solvated pair. Softness convention S = 1/(2η).
- C9. **Mulliken charges** with a diffuse (++) basis: only fragment totals are reported, and only as indicative.
- C10. **Software**: PySCF/geomeTRIC/Multiwfn/ParaView, not Gaussian/VMD; B3LYP defined as in Gaussian (VWN-RPA); grid level 4 (not 3); density fitting with PySCF's default auxiliary basis (identity not logged).
- C11. **Starting points**: both final optimisations began from geometries produced by earlier grid-3 runs (and those from a 6-31G(d) pre-optimisation) that are not archived.
- C12. **Naming**: the case file's "1H-benzotriazol-1-ium" does not match the CIF protonation pattern (H on N4 and N6); verify the cation name against the source paper.
- C13. **X-ray reference**: room-temperature structure (298 K), R1 = 0.067, H atoms refined isotropically with short X–H; nitro O2/O3 have large displacement parameters. DFT vs X-ray bond comparisons should use heavy atoms.
- C14. **NCI screening**: RDG_maxrho = 0.1 a.u.; the gas-phase H4N···O1 region reaches ρ ≈ 0.099 a.u., at the edge of the screen.

## 8. Figure data inventory (files opened and checked)

Paths: `T` = `molecules/btpc/results/<case>/` (tracked), `R` = `data/r2/molecules/btpc/results/<case>/` (raw results archive, read-only). `<case>` = btpc_pcm | btpc_gas. Tracked copies are byte-identical to the archive copies (checked with cmp).

| Figure topic | File(s) | Content checked |
| --- | --- | --- |
| Structure / labels | `T/opt.xyz` | 33 atoms, Å, CIF atom order within each ion, picrate (0–17) then cation (18–32) (labels from `molecules/btpc/cases/<case>.json:cif_labels`, strip `#4_666`). X-ray start: `molecules/btpc/input/btpc_cluster.xyz` (centred cluster, Å). Existing render: T/fig/structure_labels.png (1800×1500, 300 dpi). |
| Optimisation trajectory | `R/opt_traj.xyz` | multi-frame xyz, Å; comment "step k  E = … Ha". PCM 71 frames (0–70), gas 3 frames. Suitable for an energy/inter-ring-angle vs step plot (Section 4.10). |
| Crystal H-bond network | `molecules/btpc/input/BTPC.cif; papers/btpc/data/derived/crystal_hbonds.csv` | CIF fractional coordinates + 4 symmetry ops; CSV columns D,H,A,symmetry_code_of_A,symmetry_operation_of_A,D-H_A,H...A_A,D...A_A,D-H...A_deg,kind. |
| Frequencies | `T/freq.json (+ T/freq_step0.005.json for PCM)` | keys frequencies_cm-1 (93, sorted, imaginary negative), frequencies_pyscf_cm-1, ZPE_Ha, H298_Ha, G298_Ha, hessian_method, n_imaginary, tr_residual_cm-1. No intensities: a stick "spectrum" can only show positions. |
| Hessian / normal modes | `R/hessian.npy` | float64 array (33, 33, 3, 3), Ha/Bohr², Cartesian, same atom order and frame as opt.xyz (Bohr = Å/0.529177210903). PCM file is the step-0.01 Hessian (the 0.005 one was overwritten). Masses as used by PySCF for freq.json: H 1.008, C 12.011, N 14.007, O 15.999 (reproduce freq.json to < 1e-10 cm⁻¹). Imaginary-mode arrows can be drawn from its lowest eigenvector. |
| Orbitals (3-D) | `R/homo.cube, R/lumo.cube` | Gaussian cube, Bohr, 80×80×80; PCM origin (−12.654, −13.796, −9.456), spacing (0.3271, 0.3708, 0.2645) Bohr; gas origin (−12.138, −13.783, −9.693), spacing (0.3175, 0.3601, 0.2618). Values in a.u.; render at ±0.02. Existing renders T/fig/homo.png, lumo.png (1600×1200). |
| Orbital energies / DOS | `T/fmo.json:mo_energies_eV; papers/btpc/data/derived/mo_fragment_composition_<case>.csv` | fmo.json: 603 (PCM) / 602 (gas) energies in eV, n_occ 89. CSV (derived): index0,label,energy_eV,occupation,picrate_pct,cation_pct for every MO (Mulliken, validated) — enough for a fragment-resolved DOS/PDOS. Existing T/fig/mo_levels.png (1260×1650). |
| Wavefunction | `R/final.molden` | PySCF 2.14.0 molden, [Atoms] in Bohr, spherical [5d], 606 AOs. PCM: 603 MOs zero-padded to 606 (the 3 padded vectors have energy max+1000 a.u. and zero coefficients); gas: 602 MOs, not padded (Multiwfn read 602 orbitals). |
| Density / MEP 3-D | `R/density.cube, R/esp.cube` | density 80³ (same grid as orbitals), e/Bohr³; ESP 64³ (PCM spacing 0.4102, 0.4650, 0.3317 Bohr; gas 0.3982, 0.4515, 0.3283), a.u. Map ESP onto the 0.001 a.u. density isosurface; colour range = T/mep_range.json (vmin/vmax a.u.). Existing T/fig/mep.png. |
| MEP surface points | `R/vtx.txt` | first line = vertex count (PCM 26483, gas 25344); then 6 columns: x, y, z (Bohr), ESP (a.u.), ESP (eV), ESP (kcal/mol). Exact Multiwfn surface — best source for an ESP histogram / area-vs-ESP plot. |
| MEP extrema | `T/mep_extrema.csv, R/surfanalysis.txt` | CSV: type,index,global,ESP_au,ESP_eV,ESP_kcal_mol,x_A,y_A,z_A,nearest_atoms (27 rows each); surfanalysis.txt: Multiwfn lists (#, a.u., eV, kcal/mol, X/Y/Z Å), * = global. Statistics: T/mep_surface.json and R/logs/mwfn_mep.out (variances, skewness, sphericity only in the log). |
| NCI scatter | `R/rdg_scatter.txt, R/rdg_scatter.npz` | txt: 5 columns x, y, z (Å), sign(λ₂)ρ (a.u.), RDG (a.u.); PCM 1,771,560 rows, gas 1,764,000; points screened by RDG_maxrho (ρ ≥ 0.1 a.u.) carry RDG = 100 exactly (PCM 94,395 rows, gas 103,505); a further 104,320 (PCM) / 89,732 (gas) low-density points have RDG > 100 naturally (RDG ≥ 100 in total: PCM 198,715, gas 193,237). npz: float32 x=sign(λ₂)ρ, y=RDG, filtered to −0.25 < x < 0.06 and y < 2.5 (PCM 365,105 points, gas 403,221), n_total. Existing T/fig/rdg_scatter.png (x −0.05…0.05) and rdg_scatter_wide.png. |
| NCI 3-D | `R/rdg.cube, R/sl2r.cube` | Multiwfn cubes, Bohr; PCM 126×148×95 from (−9.654, −10.796, −6.456) spacing 0.158553; gas 125×147×96 from (−9.138, −10.783, −6.693) spacing 0.153832. Isosurface of rdg at 0.5, coloured by sl2r in −0.035…+0.02 a.u. Existing T/fig/nci.png. Spike table T/nci_spikes.csv (sign_l2_rho_au,min_RDG,n_points,type). |
| ELF / LOL maps | `R/{elf,lol}_{picrate_ring,benzotriazolium_ring}.npz and *_plane.txt` | npz keys X, Y (in-plane coordinates, Å, 300×300), F (values), atoms_xy (33×2 projected atom positions), in_plane (bool, \|offset\| < 0.6 Å), normal, center. plane.txt: 90,000 rows, 6 columns X, Y, Z (Å), graph-x, graph-y, value. Picrate plane spans ±4.36 × ±4.20 Å (steps 0.029 × 0.028 Å), cation plane ±3.81 × ±4.22 Å (gas ±4.23; steps 0.026 × 0.028 Å). Existing T/fig/{elf,lol}_*.png (1860×1500, jet colour map, 0–1). |
| ParaView inputs | `R/vtk/*.vtk` | legacy VTK (ASCII) in Å converted by dft/render.py: density, esp, homo, lumo, rdg, sl2r volumes; atoms_<El>.vtk and bonds.vtk polydata. |
| Cross-case charts | `results/charts/btpc/*.png` | bond_lengths_comparison, frequency_spectrum, hbond_geometry, mep_extrema, mo_levels_comparison, nci_spikes, optimisation_convergence, reactivity_descriptors (300 dpi, made by dft/charts.py from the tracked JSON/CSV). Rebuild for the paper rather than reuse; the SUMMARY.md title there is truncated. |
| Multiwfn inputs/logs | `R/mwfn_*.in, R/logs/mwfn_*.out, R/mwfn_home/settings.ini` | exact menu sequences (MEP 12,0,1,7,-1,q; NCI 20,1,3,2,3,0,0,q; ELF 4,9,1,300,300,4,atoms,-6,-5,q; LOL same with 10) and the settings actually used (RDG_maxrho 0.1, aug2D 6.0, iESPcode 2 PCM / 1 gas). |

## Appendix — inconsistencies found and how they were resolved

- RESULTS.md (both cases) says "integration grid level 3"; the runs used level 4. Authoritative: level 4 (`inconsistency.grid_level_text`; sources: molecules/btpc/results/*/RESULTS.md vs verify_grid.log (not deposited)).
- btpc_pcm RESULTS.md methods say step 0.01 Bohr (true for freq.json) while the dft/frequencies.py comment says 0.005 is what "every completed case used". Authoritative: freq.json hessian_method (0.01 Bohr); freq_step0.005.json is the 0.005 check (`inconsistency.pcm_hessian_step_text`; sources: molecules/btpc/results/btpc_pcm/RESULTS.md; dft/frequencies.py comment).
- handover README and opt pipeline text claim full convergence with Gaussian criteria; btpc_pcm met only the force criteria. Authoritative: optimize.log (`inconsistency.pcm_convergence_text`; sources: hand-over package README of Sep 2026 (not deposited) vs optimize.log).
- dft/common.py defaults to Lebedev order 17; runs used 29. Authoritative: 29 (not directly verified for btpc_pcm) (`inconsistency.lebedev_text`; sources: run_case.sh; verify_grid.log (not deposited)).
- btpc_pcm RESULTS.md: "1 imaginary -> NOT a minimum" and "(nan min)" for the optimisation time; SUMMARY/charts say lowest -3.2. Authoritative: consistent with freq.json; wall time recovered from STAGES.md (`inconsistency.pcm_results_md_minimum`; sources: molecules/btpc/results/btpc_pcm/RESULTS.md).
- STAGES.md prints PCM HOMO/LUMO composition as 99.9/0.1 for both while RESULTS.md table shows HOMO-1 0/100 - both consistent with fmo.json (HOMO-1 is 99.8% cation). Authoritative: no conflict; noted to avoid confusion (`inconsistency.fmo_label_threshold`; sources: STAGES.md, RESULTS.md, fmo.json).
- dft/mwfn.py comment says rdg_scatter.txt x,y,z are in Bohr; Multiwfn output states Angstrom and the values match Å. Authoritative: Angstrom (`inconsistency.nci_scatter_units`; sources: dft/mwfn.py:nci_post vs logs/mwfn_nci.out ("Column 1/2/3: X/Y/Z in Angstrom")).
- results/charts/btpc/SUMMARY.md title truncated to "1H-benzotriazol-1-ium 2". Authoritative: cosmetic (`inconsistency.chart_title`; sources: results/charts/btpc/SUMMARY.md (pipeline chart summary, not deposited)).
- The project notes (not deposited) say frequencies only on the primary (PCM) model; btpc_gas nevertheless has a full FD Hessian (0 imaginary). Authoritative: freq.json exists and is valid (`inconsistency.gas_freq_policy`; sources: molecules/btpc/results/btpc_gas/freq.json; STAGES.md).
- The auxiliary basis for density fitting is not recorded in any log; `protocol.aux_basis` gives the inference from the PySCF source and marks it as unverified.
- `btpc_pcm` was not included in `verify_grid.py`; its grid/Lebedev settings are inferred (see `protocol.grid_level`).

## Writer's additions (manuscript, 23 Sep 2026)

Added by the manuscript writer through section 7 of `build_facts.py` (OUTLINE.md section 9): arithmetic on existing keys (with an `expr` that numbers_audit.py re-evaluates), quantities recomputed from the stored trajectory, Hessians and plane files, and plotting constants read from the figure scripts.

| key | value | unit | source / formula |
| --- | --- | --- | --- |
| `xray.nh_minus_ho_asym` | 0.8677 | Å | xray.hb_N4_H4N_O1_ha - xray.hb_N4_H4N_O1_dh |
| `xray.nn_bond_asym` | 0.003653 | Å | xray.r_N5_N6 - xray.r_N4_N5 |
| `pcm.nn_bond_asym` | 0.007778 | Å | pcm.r_N5_N6 - pcm.r_N4_N5 |
| `gas.nn_bond_asym` | 0.03189 | Å | gas.r_N5_N6 - gas.r_N4_N5 |
| `gas.lumo_cation_pct` | 99.6 | % | 100 - gas.lumo_picrate_pct |
| `pcm.lumo_cation_pct` | 0.1 | % | 100 - pcm.lumo_picrate_pct |
| `pcm.homo_m1_cation_pct` | 99.8 | % | 100 - pcm.homo_m1_picrate_pct |
| `gas.homo_cation_pct` | 0.6 | % | 100 - gas.homo_picrate_pct |
| `saddle.imag_mode_rigid_pct` | 99.93 | % | 100 * saddle.imag_mode_rigid_fraction |
| `gas.mode1_rigid_pct` | 99.06 | % | 100 * gas.mode1_rigid_fraction |
| `saddle.imag_mode_cation_pct` | 66.31 | % | 100 * (1 - saddle.imag_mode_picrate_weight) |
| `saddle.drift_0_to_70_pct` | 87.45 | % | 100 * saddle.drift_0_to_70_overlap_imag |
| `saddle.drift_40_to_70_pct` | 93.67 | % | 100 * saddle.drift_40_to_70_overlap_imag |
| `saddle.drift_60_to_70_pct_9cm` | 95.3 | % | data/r2/molecules/btpc/results/btpc_pcm/opt_traj.xyz + hessian.npy; derived: 100 x squared overlap of the mass-weighted displacement frame 60 -> 70 (aligned, rigid rotation removed) with the 9.0 cm-1 mode (second lowest) |
| `saddle.hbond_max_excursion` | 0.04811 | Å | data/r2/molecules/btpc/results/btpc_pcm/opt_traj.xyz; derived: largest \|d(step k) - d(step 0)\| of N4-H4N, H4N...O1, N4...O1 over the 71 frames (H4N...O1); N4-H4N excursion 0.005, full range 0.005; H4N...O1 excursion 0.048, full range 0.051; N4...O1 excursion 0.018, full range 0.018 |
| `saddle.ring_rotation_0_to_70` | 52.32 | deg | data/r2/molecules/btpc/results/btpc_pcm/opt_traj.xyz; derived: inter-ring angle at step 70 minus step 0 (5.8 -> 58.2 deg) |
| `saddle.ring_rotation_40_to_70` | 37.55 | deg | data/r2/molecules/btpc/results/btpc_pcm/opt_traj.xyz; derived: inter-ring angle at step 70 minus step 40 (20.6 -> 58.2 deg) |
| `saddle.energy_drop_40_to_70` | -0.3015 | kcal mol-1 | data/r2/molecules/btpc/results/btpc_pcm/opt_traj.xyz; derived: E(step 70) - E(step 40) |
| `saddle.step_change_mean` | 0.1537 | cm-1 | molecules/btpc/results/btpc_pcm/freq.json vs freq_step0.005.json; derived: mean \|nu(0.01) - nu(0.005)\| over the 93 modes |
| `saddle.step_change_max_mode` | 1521 | cm-1 | molecules/btpc/results/btpc_pcm/freq.json:frequencies_cm-1; derived: wavenumber in the authoritative 0.01 Bohr Hessian of the mode with the largest change (mode 78); the same mode lies at 1522.6 cm-1 in the 0.005 Bohr Hessian; signed change -1.886 cm-1 |
| `saddle.step_change_max_signed` | -1.886 | cm-1 | molecules/btpc/results/btpc_pcm/freq.json vs freq_step0.005.json; derived: nu(0.01) - nu(0.005) of that mode |
| `pcm.nu_C3_H3_stretch` | 3230 | cm-1 | data/r2/molecules/btpc/results/btpc_pcm/hessian.npy; derived: the mode in 3000-3400 cm-1 with the largest H3 weight (picrate C-H) |
| `pcm.nu_C5_H5_stretch` | 3234 | cm-1 | data/r2/molecules/btpc/results/btpc_pcm/hessian.npy; derived: the mode in 3000-3400 cm-1 with the largest H5 weight (picrate C-H) |
| `pcm.nu_benzo_CH_stretches` | 3188.6, 3199.4, 3216.3, 3218.2 | cm-1 | data/r2/molecules/btpc/results/btpc_pcm/hessian.npy; derived: the four C-H stretches of the benzo ring (H8-H11), coupled, ascending |
| `pcm.benzo_ch_spread` | 29.58 | cm-1 | data/r2/molecules/btpc/results/btpc_pcm/hessian.npy; derived: highest minus lowest of the four benzo C-H stretches |
| `gas.nu_C3_H3_stretch` | 3226 | cm-1 | data/r2/molecules/btpc/results/btpc_gas/hessian.npy; derived: the mode in 3000-3400 cm-1 with the largest H3 weight (picrate C-H) |
| `gas.nu_C5_H5_stretch` | 3232 | cm-1 | data/r2/molecules/btpc/results/btpc_gas/hessian.npy; derived: the mode in 3000-3400 cm-1 with the largest H5 weight (picrate C-H) |
| `gas.nu_benzo_CH_stretches` | 3183.4, 3195.7, 3204.6, 3207.1 | cm-1 | data/r2/molecules/btpc/results/btpc_gas/hessian.npy; derived: the four C-H stretches of the benzo ring (H8-H11), coupled, ascending |
| `gas.benzo_ch_spread` | 23.7 | cm-1 | data/r2/molecules/btpc/results/btpc_gas/hessian.npy; derived: highest minus lowest of the four benzo C-H stretches |
| `diff.ch_stretch_max_shift` | 11.66 | cm-1 | pcm/gas nu_C3_H3_stretch, nu_C5_H5_stretch, nu_benzo_CH_stretches; derived: largest \|PCM - gas\| over the six C-H stretches (benzo modes paired in ascending order); shifts +4.0, +2.2, +5.2, +3.7, +11.7, +11.1 |
| `pcm.cation_plane_offset_O7` | 2.181 | Å | data/r2/molecules/btpc/results/btpc_pcm/elf_benzotriazolium_ring.npz + opt.xyz; derived: distance of the nitro O7 from the ELF/LOL cation plane (Fig. 12 marks its projection in the gas phase) |
| `gas.cation_plane_offset_O7` | 1.26 | Å | data/r2/molecules/btpc/results/btpc_gas/elf_benzotriazolium_ring.npz + opt.xyz; derived: distance of the nitro O7 from the ELF/LOL cation plane (Fig. 12 marks its projection in the gas phase) |
| `pcm.picrate_plane_offset_O7` | 0.5047 | Å | data/r2/molecules/btpc/results/btpc_pcm/elf_picrate_ring.npz + opt.xyz; derived: distance of the nitro O7 from the ELF/LOL picrate-ring plane (SI ELF/LOL figure) |
| `gas.picrate_plane_offset_O7` | 0.4256 | Å | data/r2/molecules/btpc/results/btpc_gas/elf_picrate_ring.npz + opt.xyz; derived: distance of the nitro O7 from the ELF/LOL picrate-ring plane (SI ELF/LOL figure) |
| `protocol.opt_threshold_energy` | 1e-06 | Ha | data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log ("Convergence Criteria"); parsed from protocol.opt_criteria (geomeTRIC/Gaussian default convergence threshold) |
| `protocol.opt_threshold_grad_rms` | 0.0003 | Ha/Bohr | data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log ("Convergence Criteria"); parsed from protocol.opt_criteria (geomeTRIC/Gaussian default convergence threshold) |
| `protocol.opt_threshold_grad_max` | 0.00045 | Ha/Bohr | data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log ("Convergence Criteria"); parsed from protocol.opt_criteria (geomeTRIC/Gaussian default convergence threshold) |
| `protocol.opt_threshold_disp_rms` | 0.0012 | Å | data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log ("Convergence Criteria"); parsed from protocol.opt_criteria (geomeTRIC/Gaussian default convergence threshold) |
| `protocol.opt_threshold_disp_max` | 0.0018 | Å | data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log ("Convergence Criteria"); parsed from protocol.opt_criteria (geomeTRIC/Gaussian default convergence threshold) |
| `protocol.pcm_lebedev_points` | 302 |  | run_case.sh (PCM_LEBEDEV_ORDER 29); PySCF 2.14.0 pyscf/dft/gen_grid.py LEBEDEV_ORDER table; number of Lebedev points of order 29 per atomic sphere of the PCM cavity (see protocol.pcm_lebedev_order) |
| `gas.nci_region01_min` | -0.0988 | a.u. | data/r2/molecules/btpc/results/btpc_gas/sl2r.cube + rdg.cube; derived: most negative sign(l2)rho inside the gas-phase N4-H4N...O1 RDG < 0.5 region (range recorded with gas.nci_region01); close to the RDG_maxrho = 0.1 a.u. screen |
| `plot.nci_colour_min` | -0.035 | a.u. | papers/_shared/figkit rdg.NCI_RANGE (fig_nci.py VRANGE); lower end of the linear sign(l2)rho colour scale (blue) |
| `plot.nci_colour_max` | 0.02 | a.u. | papers/_shared/figkit rdg.NCI_RANGE (fig_nci.py VRANGE); upper end of the linear sign(l2)rho colour scale (red) |
| `plot.nci_colour_mid` | -0.0075 | a.u. | (plot.nci_colour_min + plot.nci_colour_max) / 2 |
| `plot.mep_colour_min` | -55 | kcal mol-1 | papers/btpc/figures/scripts/fig_mep.py (vr, outer multiples of 5 of all surface values); lower end of the MEP colour scale, both models |
| `plot.mep_colour_max` | 110 | kcal mol-1 | papers/btpc/figures/scripts/fig_mep.py (vr, outer multiples of 5 of all surface values); upper end of the MEP colour scale, both models |
| `plot.mep_hist_bin` | 2.5 | kcal mol-1 | papers/btpc/figures/scripts/fig_mep.py (bins); bin width of the ESP surface-area histogram |
| `plot.mep_label_threshold` | 30 | kcal mol-1 | papers/btpc/figures/scripts/fig_mep.py LABEL_ABOVE; extrema with \|V\| above this are labelled in the extrema strip |
| `plot.pdos_fwhm` | 0.3 | eV | papers/btpc/figures/scripts/fig_levels_pdos.py FWHM; Gaussian broadening (full width at half maximum) of the DOS |
| `plot.orbital_box_margin` | 7 | bohr | papers/btpc/figures/scripts/figS_orbitals_extra.py MARGIN_BOHR; margin of the orbital evaluation box beyond the outermost nuclei (3.70 Å) |

## Revision round 1 additions (24 Sep 2026)

Added through section 8 of `build_facts.py` in answer to the round-1 audits (`papers/btpc/audit/*_r1.md`; response in `audit/response_r1.md`). All are local post-processing of stored files: CIF values, the IEF-PCM trajectory, the two Hessians, freq.json (quasi-RRHO entropy, with the RRHO entropy reproduced first as a check), Loewdin populations with the rebuilt overlap matrix, and the density topology of the N4-H4N...O1 contact evaluated from final.molden by `data/wfn_topology.py` (validated against density.cube and the Multiwfn ELF plane; sign(lambda2)rho at the critical point reproduces the NCI region 01 value). Reading: at the H4N...O1 critical point both models have lap(rho) > 0 and H < 0, the transit region between closed-shell and shared-shell interactions; the gas-phase bond is much further into it (H -0.031 against -0.006 a.u., |V|/G 1.44 against 1.14). The Loewdin picrate charge differs between the models by about the same amount as the Mulliken one (0.13 against 0.12 e), and the H4N charge is the same in both models under either partition, so the difference is not an artefact of where H4N is counted.

| key | value | unit | source / formula |
| --- | --- | --- | --- |
| `xray.bond_C2_C3` | 1.352 | Å | molecules/btpc/input/BTPC.cif:_geom_bond C2 C3; CIF value with s.u. 1.352(4) |
| `xray.bond_C6_C5` | 1.367 | Å | molecules/btpc/input/BTPC.cif:_geom_bond C6 C5; CIF value with s.u. 1.367(4) |
| `xray.bond_C7_C12` | 1.382 | Å | molecules/btpc/input/BTPC.cif:_geom_bond C7 C12; CIF value with s.u. 1.382(4) |
| `software.figures_env` | Matplotlib 3.11.2; NumPy 2.5.3; SciPy 1.18.1; PyVista 0.49.0; VTK 9.7.0; scikit-image 0.26.0; gemmi 0.7.5; RDKit 2026.3.6 |  | tools/paper-env (installed package metadata, read with importlib.metadata when facts.json is built); environment of the figure and table scripts (papers/btpc/figures/scripts, tables/make_tables.py); the calculations used the versions in software.numpy_scipy_matplotlib on the compute host |
| `xray.ueq_O2` | 0.1096 | Å^2 | molecules/btpc/input/BTPC.cif:_atom_site_U_iso_or_equiv; U(eq) 0.1096(12) Å²; the two largest of all 33 atoms (nitro group N1 on C2) |
| `xray.ueq_O3` | 0.1006 | Å^2 | molecules/btpc/input/BTPC.cif:_atom_site_U_iso_or_equiv; U(eq) 0.1006(10) Å²; the two largest of all 33 atoms (nitro group N1 on C2) |
| `xray.ueq_O_other` | 0.0544, 0.0879 | Å^2 | molecules/btpc/input/BTPC.cif:_atom_site_U_iso_or_equiv; range of U(eq) over the other five O atoms (O1, O4-O7) |
| `saddle.h4n_o7_max_excursion` | 0.1427 | Å | data/r2/molecules/btpc/results/btpc_pcm/opt_traj.xyz; derived: largest \|d(step k) - d(step 0)\| of H4N...O7 over the 71 frames (step 0 2.348, step 70 2.442, maximum 2.490 Å); N4...O7 0.065 Å |
| `saddle.hb_angle_step0` | 154.5 | deg | data/r2/molecules/btpc/results/btpc_pcm/opt_traj.xyz; derived: angle N4-H4N...O1 at step 0 (inter-ring angle 5.8 deg) |
| `saddle.hb_angle_opening` | 12.65 | deg | pcm.hb_N4_H4N_O1_angle - saddle.hb_angle_step0 |
| `saddle.nh_no_max_excursion` | 0.01788 | Å | data/r2/molecules/btpc/results/btpc_pcm/opt_traj.xyz; derived: largest \|d(step k) - d(step 0)\| of N4-H4N and N4...O1 together over the 71 frames |
| `pcm.nu_N4_H4N_stretch_h_weight` | 0.8895 |  | data/r2/molecules/btpc/results/btpc_pcm/hessian.npy; derived: share of the mass-weighted eigenvector of the N4-H4N stretch carried by H4N |
| `pcm.nu_N4_H4N_stretch_offaxis_deg` | 0.8018 | deg | data/r2/molecules/btpc/results/btpc_pcm/hessian.npy; derived: angle between the Cartesian displacement of H4N in the N4-H4N stretch and the N4-H4N bond (fig_vibrations.py check_offaxis_caption) |
| `saddle.imag_mode_picrate_rotation_share` | 0.9716 |  | data/r2/molecules/btpc/results/btpc_pcm/hessian.npy; derived: rotational part of the rigid-body motion of the picrate in the lowest mode (rotation / (rotation + translation) of its least-squares rigid fit); translation 0.010, rotation 0.327 of the whole mode |
| `saddle.imag_mode_cation_rotation_share` | 0.9726 |  | data/r2/molecules/btpc/results/btpc_pcm/hessian.npy; derived: rotational part of the rigid-body motion of the cation in the lowest mode (rotation / (rotation + translation) of its least-squares rigid fit); translation 0.018, rotation 0.645 of the whole mode |
| `saddle.imag_mode_axes_angle_deg` | 154.4 | deg | data/r2/molecules/btpc/results/btpc_pcm/hessian.npy; derived: angle between the rotation axes (angular velocities) of the two ions in the lowest mode; > 90 deg = the ions turn in opposite senses |
| `saddle.imag_mode_cation_axis_to_centroid_line_deg` | 10.77 | deg | data/r2/molecules/btpc/results/btpc_pcm/hessian.npy + opt.xyz; derived: angle between the rotation axis of the cation in the lowest mode and the line joining the picrate and cation ring centroids |
| `gas.nu_N4_H4N_stretch_h_weight` | 0.7932 |  | data/r2/molecules/btpc/results/btpc_gas/hessian.npy; derived: share of the mass-weighted eigenvector of the N4-H4N stretch carried by H4N |
| `gas.nu_N4_H4N_stretch_offaxis_deg` | 19.34 | deg | data/r2/molecules/btpc/results/btpc_gas/hessian.npy; derived: angle between the Cartesian displacement of H4N in the N4-H4N stretch and the N4-H4N bond (fig_vibrations.py check_offaxis_caption) |
| `gas.mode1_picrate_rotation_share` | 0.9503 |  | data/r2/molecules/btpc/results/btpc_gas/hessian.npy; derived: rotational part of the rigid-body motion of the picrate in the lowest mode (rotation / (rotation + translation) of its least-squares rigid fit); translation 0.022, rotation 0.412 of the whole mode |
| `gas.mode1_cation_rotation_share` | 0.9265 |  | data/r2/molecules/btpc/results/btpc_gas/hessian.npy; derived: rotational part of the rigid-body motion of the cation in the lowest mode (rotation / (rotation + translation) of its least-squares rigid fit); translation 0.041, rotation 0.516 of the whole mode |
| `gas.mode1_axes_angle_deg` | 155.5 | deg | data/r2/molecules/btpc/results/btpc_gas/hessian.npy; derived: angle between the rotation axes (angular velocities) of the two ions in the lowest mode; > 90 deg = the ions turn in opposite senses |
| `gas.mode1_cation_axis_to_centroid_line_deg` | 28.94 | deg | data/r2/molecules/btpc/results/btpc_gas/hessian.npy + opt.xyz; derived: angle between the rotation axis of the cation in the lowest mode and the line joining the picrate and cation ring centroids |
| `gas.nu_N4_H4N_stretch_h_pct` | 79.32 | % | 100 * gas.nu_N4_H4N_stretch_h_weight |
| `pcm.nu_N4_H4N_stretch_h_pct` | 88.95 | % | 100 * pcm.nu_N4_H4N_stretch_h_weight |
| `saddle.imag_mode_picrate_rotation_pct` | 97.16 | % | 100 * saddle.imag_mode_picrate_rotation_share |
| `saddle.imag_mode_cation_rotation_pct` | 97.26 | % | 100 * saddle.imag_mode_cation_rotation_share |
| `gas.mode3_rigid_pct` | 80.61 | % | 100 * gas.mode3_rigid_fraction |
| `pcm.mode3_rigid_pct` | 93.97 | % | 100 * pcm.mode3_rigid_fraction |
| `pcm.entropy_rrho_check` | 156.7 | cal mol-1 K-1 | molecules/btpc/results/btpc_pcm/freq.json + opt.xyz; derived: translational + rotational + harmonic vibrational entropy recomputed here (sigma = 1, 92 real modes); reproduces pcm.entropy (from H298 - G298) |
| `pcm.entropy_qrrho` | 145.3 | cal mol-1 K-1 | molecules/btpc/results/btpc_pcm/freq.json + opt.xyz; derived: quasi-RRHO entropy, each vibrational term interpolated between the harmonic and a free-rotor entropy with the weight 1/(1 + (100 cm-1 / nu)^4) and B_av = 1e-44 kg m^2 (Grimme 2012); translation and rotation as in the RRHO value |
| `pcm.entropy_qrrho_lowering` | 11.43 | cal mol-1 K-1 | pcm.entropy - pcm.entropy_qrrho |
| `gas.entropy_rrho_check` | 159 | cal mol-1 K-1 | molecules/btpc/results/btpc_gas/freq.json + opt.xyz; derived: translational + rotational + harmonic vibrational entropy recomputed here (sigma = 1, 93 real modes); reproduces gas.entropy (from H298 - G298) |
| `gas.entropy_qrrho` | 147.3 | cal mol-1 K-1 | molecules/btpc/results/btpc_gas/freq.json + opt.xyz; derived: quasi-RRHO entropy, each vibrational term interpolated between the harmonic and a free-rotor entropy with the weight 1/(1 + (100 cm-1 / nu)^4) and B_av = 1e-44 kg m^2 (Grimme 2012); translation and rotation as in the RRHO value |
| `gas.entropy_qrrho_lowering` | 11.65 | cal mol-1 K-1 | gas.entropy - gas.entropy_qrrho |
| `pcm.lowdin_charge_picrate` | -0.8534 | e | data/r2/molecules/btpc/results/btpc_pcm/final.molden; derived: Loewdin fragment charge, diag(S^1/2 D S^1/2), overlap rebuilt from the molden basis (smallest overlap eigenvalue 4.0e-07); cation = -this |
| `pcm.mulliken_charge_H4N` | 0.7444 | e | data/r2/molecules/btpc/results/btpc_pcm/final.molden; derived: Mulliken atomic charge of H4N (counted with the cation) |
| `pcm.lowdin_charge_H4N` | 0.1302 | e | data/r2/molecules/btpc/results/btpc_pcm/final.molden; derived: Loewdin atomic charge of H4N (counted with the cation) |
| `gas.lowdin_charge_picrate` | -0.722 | e | data/r2/molecules/btpc/results/btpc_gas/final.molden; derived: Loewdin fragment charge, diag(S^1/2 D S^1/2), overlap rebuilt from the molden basis (smallest overlap eigenvalue 2.8e-07); cation = -this |
| `gas.mulliken_charge_H4N` | 0.7484 | e | data/r2/molecules/btpc/results/btpc_gas/final.molden; derived: Mulliken atomic charge of H4N (counted with the cation) |
| `gas.lowdin_charge_H4N` | 0.1202 | e | data/r2/molecules/btpc/results/btpc_gas/final.molden; derived: Loewdin atomic charge of H4N (counted with the cation) |
| `diff.lowdin_charge_picrate` | -0.1314 | e | pcm.lowdin_charge_picrate - gas.lowdin_charge_picrate |
| `diff.mulliken_charge_H4N` | -0.003984 | e | pcm.mulliken_charge_H4N - gas.mulliken_charge_H4N |
| `diff.lowdin_charge_H4N` | 0.009957 | e | pcm.lowdin_charge_H4N - gas.lowdin_charge_H4N |
| `pcm.bcp_H4N_O1_rho` | 0.04998 | a.u. | data/r2/molecules/btpc/results/btpc_pcm/final.molden (papers/btpc/data/wfn_topology.py); derived: electron density at the (3,-1) critical point between H4N and O1 (Newton search from the density minimum on the H4N-O1 line; \|grad rho\| 3e-14); lambda = -0.0908, -0.0902, 0.3294 a.u.; sign(lambda2)rho agrees with the lowest-RDG point of NCI region 01 (-0.0496) to 0.002 |
| `pcm.bcp_H4N_O1_lap` | 0.1484 | a.u. | data/r2/molecules/btpc/results/btpc_pcm/final.molden (papers/btpc/data/wfn_topology.py); derived: Laplacian of the density at that critical point (e bohr^-5) |
| `pcm.bcp_H4N_O1_G` | 0.04327 | a.u. | data/r2/molecules/btpc/results/btpc_pcm/final.molden (papers/btpc/data/wfn_topology.py); derived: positive-definite kinetic energy density G = 1/2 sum n_i \|grad phi_i\|^2 (hartree bohr^-3) |
| `pcm.bcp_H4N_O1_V` | -0.04944 | a.u. | data/r2/molecules/btpc/results/btpc_pcm/final.molden (papers/btpc/data/wfn_topology.py); derived: potential energy density from the local virial theorem, V = lap/4 - 2G |
| `pcm.bcp_H4N_O1_H` | -0.006175 | a.u. | data/r2/molecules/btpc/results/btpc_pcm/final.molden (papers/btpc/data/wfn_topology.py); derived: total energy density H = G + V; H < 0 with lap > 0 marks a partly covalent (shared-shell) contribution |
| `pcm.bcp_H4N_O1_VG` | 1.143 |  | data/r2/molecules/btpc/results/btpc_pcm/final.molden (papers/btpc/data/wfn_topology.py); derived: \|V\|/G; 1 < \|V\|/G < 2 is the intermediate (transit) region between closed-shell and shared interactions |
| `pcm.bcp_H4N_O1_elf` | 0.1684 |  | data/r2/molecules/btpc/results/btpc_pcm/final.molden (papers/btpc/data/wfn_topology.py); derived: ELF (Becke, with the 1e-5 addition of Multiwfn ELF_addminimal) at the critical point |
| `pcm.bcp_H4N_O1_dist_H` | 0.5311 | Å | data/r2/molecules/btpc/results/btpc_pcm/final.molden (papers/btpc/data/wfn_topology.py); derived: distance of the critical point from H4N |
| `pcm.elf_min_H4N_O1` | 0.166 |  | data/r2/molecules/btpc/results/btpc_pcm/final.molden (papers/btpc/data/wfn_topology.py); derived: minimum of ELF on the straight H4N -> O1 line between the H4N and O1 valence regions, 0.519 Å from H4N |
| `pcm.wfn_rho_rel_dev_max` | 0.0002618 |  | data/r2/molecules/btpc/results/btpc_pcm/final.molden (papers/btpc/data/wfn_topology.py) vs data/r2/molecules/btpc/results/btpc_pcm/density.cube; derived: largest relative deviation of the evaluated density from the PySCF cube over 400 random points with rho > 1e-3 a.u. |
| `pcm.wfn_elf_abs_dev_max` | 1.568e-05 |  | data/r2/molecules/btpc/results/btpc_pcm/final.molden (papers/btpc/data/wfn_topology.py) vs data/r2/molecules/btpc/results/btpc_pcm/elf_benzotriazolium_ring_plane.txt; derived: largest absolute deviation of the evaluated ELF from the Multiwfn plane map (random plane points with ELF > 0.02) |
| `pcm.bcp_nci_sign_l2_rho_diff` | 0.0004343 | a.u. | data/r2/molecules/btpc/results/btpc_pcm/final.molden (papers/btpc/data/wfn_topology.py); pcm.nci_region01; derived: \|sign(lambda2)rho at the critical point - value at the lowest-RDG point of NCI region 01\| |
| `pcm.wfn_topology_validation` | density vs density.cube: max relative deviation 2.6e-04 (400 points with rho > 1e-3 a.u.); ELF vs Multiwfn elf_benzotriazolium_ring_plane.txt: max absolute deviation 1.6e-05 (304 points) |  | data/r2/molecules/btpc/results/btpc_pcm/final.molden (papers/btpc/data/wfn_topology.py); validation of the orbital evaluator used for the bcp_* and elf_min_* keys |
| `gas.bcp_H4N_O1_rho` | 0.08676 | a.u. | data/r2/molecules/btpc/results/btpc_gas/final.molden (papers/btpc/data/wfn_topology.py); derived: electron density at the (3,-1) critical point between H4N and O1 (Newton search from the density minimum on the H4N-O1 line; \|grad rho\| 7e-14); lambda = -0.1988, -0.1941, 0.5468 a.u.; sign(lambda2)rho agrees with the lowest-RDG point of NCI region 01 (-0.0863) to 0.002 |
| `gas.bcp_H4N_O1_lap` | 0.1539 | a.u. | data/r2/molecules/btpc/results/btpc_gas/final.molden (papers/btpc/data/wfn_topology.py); derived: Laplacian of the density at that critical point (e bohr^-5) |
| `gas.bcp_H4N_O1_G` | 0.06901 | a.u. | data/r2/molecules/btpc/results/btpc_gas/final.molden (papers/btpc/data/wfn_topology.py); derived: positive-definite kinetic energy density G = 1/2 sum n_i \|grad phi_i\|^2 (hartree bohr^-3) |
| `gas.bcp_H4N_O1_V` | -0.09954 | a.u. | data/r2/molecules/btpc/results/btpc_gas/final.molden (papers/btpc/data/wfn_topology.py); derived: potential energy density from the local virial theorem, V = lap/4 - 2G |
| `gas.bcp_H4N_O1_H` | -0.03053 | a.u. | data/r2/molecules/btpc/results/btpc_gas/final.molden (papers/btpc/data/wfn_topology.py); derived: total energy density H = G + V; H < 0 with lap > 0 marks a partly covalent (shared-shell) contribution |
| `gas.bcp_H4N_O1_VG` | 1.442 |  | data/r2/molecules/btpc/results/btpc_gas/final.molden (papers/btpc/data/wfn_topology.py); derived: \|V\|/G; 1 < \|V\|/G < 2 is the intermediate (transit) region between closed-shell and shared interactions |
| `gas.bcp_H4N_O1_elf` | 0.3335 |  | data/r2/molecules/btpc/results/btpc_gas/final.molden (papers/btpc/data/wfn_topology.py); derived: ELF (Becke, with the 1e-5 addition of Multiwfn ELF_addminimal) at the critical point |
| `gas.bcp_H4N_O1_dist_H` | 0.4234 | Å | data/r2/molecules/btpc/results/btpc_gas/final.molden (papers/btpc/data/wfn_topology.py); derived: distance of the critical point from H4N |
| `gas.elf_min_H4N_O1` | 0.3302 |  | data/r2/molecules/btpc/results/btpc_gas/final.molden (papers/btpc/data/wfn_topology.py); derived: minimum of ELF on the straight H4N -> O1 line between the H4N and O1 valence regions, 0.406 Å from H4N |
| `gas.wfn_rho_rel_dev_max` | 7.264e-05 |  | data/r2/molecules/btpc/results/btpc_gas/final.molden (papers/btpc/data/wfn_topology.py) vs data/r2/molecules/btpc/results/btpc_gas/density.cube; derived: largest relative deviation of the evaluated density from the PySCF cube over 400 random points with rho > 1e-3 a.u. |
| `gas.wfn_elf_abs_dev_max` | 4.212e-05 |  | data/r2/molecules/btpc/results/btpc_gas/final.molden (papers/btpc/data/wfn_topology.py) vs data/r2/molecules/btpc/results/btpc_gas/elf_benzotriazolium_ring_plane.txt; derived: largest absolute deviation of the evaluated ELF from the Multiwfn plane map (random plane points with ELF > 0.02) |
| `gas.bcp_nci_sign_l2_rho_diff` | 0.0004638 | a.u. | data/r2/molecules/btpc/results/btpc_gas/final.molden (papers/btpc/data/wfn_topology.py); gas.nci_region01; derived: \|sign(lambda2)rho at the critical point - value at the lowest-RDG point of NCI region 01\| |
| `gas.wfn_topology_validation` | density vs density.cube: max relative deviation 7.3e-05 (400 points with rho > 1e-3 a.u.); ELF vs Multiwfn elf_benzotriazolium_ring_plane.txt: max absolute deviation 4.2e-05 (316 points) |  | data/r2/molecules/btpc/results/btpc_gas/final.molden (papers/btpc/data/wfn_topology.py); validation of the orbital evaluator used for the bcp_* and elf_min_* keys |

## Revision round 2 additions (24 Sep 2026)

Added through section 9 of `build_facts.py` in answer to the round-2 audits (`papers/btpc/audit/*_r2.md`; response in `audit/response_r2.md`). Local post-processing only. (1) The curvature of the two stored, unprojected Hessians along the overall rigid rotations and translations of the pair: an exact Hessian at a stationary point gives zero, so these values measure the numerical resolution of the finite-difference Hessians (rotations up to 7.7i cm-1, translations 0.2 cm-1 or less); the 3.2i cm-1 IEF-PCM mode lies below that resolution. The analytic gradients behind the Hessians omit the grid-weight response (`protocol.fd_gradient_grid_response`). (2) H/rho at the H4N...O1 critical point and the gas/IEF-PCM ratio of H (about 4.9): both bonds are in the transit region, and they differ in degree, not in class. (3) U(eq) of the triazole nitrogens from the CIF. (4) Three plotting constants read from the figure scripts. (5) The Delta pKa values and ranges quoted in the abstract, Fig. 1 and the Introduction, computed from the pKa values in the abstract-checked literature claims (`lit.*`). Changed values: every `*.mep_minNN` / `*.mep_maxNN` now holds the unrounded kcal mol-1 value from Multiwfn's surfanalysis.txt instead of the 3-decimal value of mep_extrema.csv (six table cells had been rounded twice and were 0.01 off; QA round 4 r4-m1).

| key | value | unit | source / formula |
| --- | --- | --- | --- |
| `pcm.hessian_rot_curvature_cm` | -7.74, 1.1, 5.01 | cm-1 | data/r2/molecules/btpc/results/btpc_pcm/hessian.npy + molecules/btpc/results/btpc_pcm/opt.xyz; derived: eigenvalues of Q^T M^-1/2 H M^-1/2 Q for the unprojected Hessian, Q = the three orthonormalised rigid rotations about the centre of mass of the whole pair, as wavenumbers (negative = imaginary); zero for an exact Hessian at a stationary point, so the spread measures the numerical resolution of the Hessian (integration grid not rotated with the molecule, gradients without grid-weight response, residual gradient) |
| `pcm.hessian_trans_curvature_cm` | -0.07, -0.05, 0.19 | cm-1 | data/r2/molecules/btpc/results/btpc_pcm/hessian.npy + molecules/btpc/results/btpc_pcm/opt.xyz; derived: eigenvalues of Q^T M^-1/2 H M^-1/2 Q for the unprojected Hessian, Q = the three orthonormalised rigid translations of the whole pair, as wavenumbers (negative = imaginary); zero for an exact Hessian at a stationary point, so the spread measures the numerical resolution of the Hessian (integration grid not rotated with the molecule, gradients without grid-weight response, residual gradient) |
| `pcm.hessian_rot_curvature_max_abs` | 7.74 | cm-1 | pcm.hessian_rot_curvature_cm; derived: largest \|value\| of the three rotational curvatures |
| `gas.hessian_rot_curvature_cm` | -6.99, -4.72, 5.54 | cm-1 | data/r2/molecules/btpc/results/btpc_gas/hessian.npy + molecules/btpc/results/btpc_gas/opt.xyz; derived: eigenvalues of Q^T M^-1/2 H M^-1/2 Q for the unprojected Hessian, Q = the three orthonormalised rigid rotations about the centre of mass of the whole pair, as wavenumbers (negative = imaginary); zero for an exact Hessian at a stationary point, so the spread measures the numerical resolution of the Hessian (integration grid not rotated with the molecule, gradients without grid-weight response, residual gradient) |
| `gas.hessian_trans_curvature_cm` | -0.07, -0.01, 0.02 | cm-1 | data/r2/molecules/btpc/results/btpc_gas/hessian.npy + molecules/btpc/results/btpc_gas/opt.xyz; derived: eigenvalues of Q^T M^-1/2 H M^-1/2 Q for the unprojected Hessian, Q = the three orthonormalised rigid translations of the whole pair, as wavenumbers (negative = imaginary); zero for an exact Hessian at a stationary point, so the spread measures the numerical resolution of the Hessian (integration grid not rotated with the molecule, gradients without grid-weight response, residual gradient) |
| `gas.hessian_rot_curvature_max_abs` | 6.99 | cm-1 | gas.hessian_rot_curvature_cm; derived: largest \|value\| of the three rotational curvatures |
| `protocol.fd_gradient_grid_response` | no: nuc_grad_method().kernel() with PySCF defaults (grid_response = False), displaced SCF conv_tol 1e-10 |  | dft/_fd_worker.py; the analytic gradients behind both finite-difference Hessians omit the response of the integration-grid weights |
| `gas.bcp_H4N_O1_H_over_rho` | -0.3519 | a.u. | gas.bcp_H4N_O1_H / gas.bcp_H4N_O1_rho |
| `pcm.bcp_H4N_O1_H_over_rho` | -0.1235 | a.u. | pcm.bcp_H4N_O1_H / pcm.bcp_H4N_O1_rho |
| `diff.bcp_H_ratio_gas_pcm` | 4.945 |  | gas.bcp_H4N_O1_H / pcm.bcp_H4N_O1_H |
| `xray.ueq_N4` | 0.0434 | Å^2 | molecules/btpc/input/BTPC.cif:_atom_site_U_iso_or_equiv; U(eq) 0.0434(5) Å²; triazole nitrogen of the cation (N5 is the unprotonated central atom) |
| `xray.ueq_N6` | 0.0563 | Å^2 | molecules/btpc/input/BTPC.cif:_atom_site_U_iso_or_equiv; U(eq) 0.0563(7) Å²; triazole nitrogen of the cation (N5 is the unprotonated central atom) |
| `xray.ueq_N5` | 0.072 | Å^2 | molecules/btpc/input/BTPC.cif:_atom_site_U_iso_or_equiv; U(eq) 0.0720(8) Å²; triazole nitrogen of the cation (N5 is the unprotonated central atom) |
| `plot.levels_degen_ev` | 0.07 | eV | papers/btpc/figures/scripts/fig_levels_pdos.py DEGEN; levels closer than this are drawn side by side |
| `plot.bond_band` | 0.02 | Å | papers/btpc/figures/scripts/fig_xray_comparison.py (parity_plot band); figS_bond_angle_deviations.py (specs); half-width of the grey band in the bond-length plots |
| `plot.angle_band` | 2 | deg | papers/btpc/figures/scripts/figS_bond_angle_deviations.py (specs); half-width of the grey band in the angle plot |
| `lit.pk_picric_acid` | 0.33 |  | papers/btpc/lit/candidates.jsonl (claims of wang2000photochemistry, davis1963spectrophotometric, alkorta2012theoretical); thermodynamic pK of picric acid in water at 25 C (davis1963spectrophotometric) |
| `lit.pk_picric_acid_spread` | 0.2, 0.8 |  | papers/btpc/lit/candidates.jsonl (claims of wang2000photochemistry, davis1963spectrophotometric, alkorta2012theoretical); spread of earlier picric acid pK determinations (davis1963spectrophotometric) |
| `lit.pka_benzotriazolium_flash` | 0.42 |  | papers/btpc/lit/candidates.jsonl (claims of wang2000photochemistry, davis1963spectrophotometric, alkorta2012theoretical); pKa of the N-protonated benzotriazolium cation, flash photolysis, I = 1.0 M (wang2000photochemistry) |
| `lit.pka_benzotriazolium_older` | 1.6 |  | papers/btpc/lit/candidates.jsonl (claims of wang2000photochemistry, davis1963spectrophotometric, alkorta2012theoretical); basicity pKa of benzotriazole from the older review quoted by Alkorta et al. (alkorta2012theoretical) |
| `lit.dpka_picric_acid_flash` | 0.09 |  | papers/btpc/lit/candidates.jsonl (claims of wang2000photochemistry, davis1963spectrophotometric, alkorta2012theoretical); derived: pKa(BTAH+) 0.42 (flash photolysis, I = 1.0 M) - pK(picric acid) 0.33 (25 C) |
| `lit.dpka_picric_acid_flash_range` | -0.38, 0.22 |  | papers/btpc/lit/candidates.jsonl (claims of wang2000photochemistry, davis1963spectrophotometric, alkorta2012theoretical); derived: 0.42 minus the spread 0.2-0.8 of published picric acid values |
| `lit.dpka_picric_acid_older` | 1.27 |  | papers/btpc/lit/candidates.jsonl (claims of wang2000photochemistry, davis1963spectrophotometric, alkorta2012theoretical); derived: basicity pKa 1.6 of the older review quoted by Alkorta et al. - pK(picric acid) 0.33 |
| `lit.dpka_picric_acid_older_range` | 0.8, 1.4 |  | papers/btpc/lit/candidates.jsonl (claims of wang2000photochemistry, davis1963spectrophotometric, alkorta2012theoretical); derived: 1.6 minus the spread 0.2-0.8 |
| `lit.dpka_picric_acid_range` | -0.38, 1.4 |  | papers/btpc/lit/candidates.jsonl (claims of wang2000photochemistry, davis1963spectrophotometric, alkorta2012theoretical); derived: lowest (flash-photolysis value, highest picric acid value) to highest (older value, lowest picric acid value) Delta pKa over the published values; printed as "about -0.4 to 1.4" |
| `pcm.hessian_step_bohr` | 0.01 | bohr | pcm.hessian_method; finite-difference displacement step of the reported Hessian |
| `gas.hessian_step_bohr` | 0.005 | bohr | gas.hessian_method; finite-difference displacement step of the reported Hessian |
| `saddle.hessian_step_first_bohr` | 0.005 | bohr | saddle.freq_imag_step0p005; step of the first IEF-PCM Hessian (the check) |

## Fact-check log

Independent fact-check, 2026-09-23. Every value below was re-derived from the raw or tracked result files, never from facts.*. Sources used: BTPC.cif, btpc_cluster.xyz, opt.xyz, opt_traj.xyz, optimize.log, freq.json, freq_step0.005.json, hessian.npy, fmo.json, final.molden (with an AO overlap matrix coded separately), surfanalysis.txt, mwfn_*.out, vtx.txt, nci_spikes.csv, rdg_scatter.txt, rdg.cube, sl2r.cube, the ELF/LOL npz and plane.txt files, cube headers, STAGES.md, dft/*.py, run_case.sh, verify_grid.log, setup_linux.log and docs/environment.

**Result:** about 1,100 values or claims were re-derived or opened. In addition, 345 table cells and all 293 backticked keys were cross-checked between facts.md and facts.json. 10 were wrong and all 10 are fixed. 6 are open questions.

### What was confirmed
- **CIF crystal data:** Table 1.1 (32 items) and Table 1.3 (22 bonds and angles with s.u.) match.
- **Table 1.2:** an independent symmetry search reproduces exactly the 11 contacts, with the same symmetry codes and geometry. Search settings: H···A < 2.9 Å, angle > 100°, N/O acceptors, all images within ±2 cells. Also confirmed: the closest contact in the asymmetric unit (O4···N5 2.989 Å) and the triazole stacking (4.106 / 3.341 Å). No other ring pair is closer than 4.5 Å.
- **Geometry (Tables 4.2, 4.3b, 4.4):** match. This covers every H-bond, 28 bonds, 6 angles, O···O distances, the N–H/H···O asymmetry, the inter-ring angle, the centroid distance and the nitro twists. Plane definitions: the 6-C picrate ring, the 9-atom cation, and the O–N–O plane for the nitro groups.
- **RMSD block:** matches (34/26 bonds, 52/36 angles, Kabsch fits on 25/16/9 atoms, and the top-5 lists).
- **Optimisation:** matches. This covers the energies, step counts, the last-step gradient, displacement, ΔE and trust radius from optimize.log, the −0.440 kcal mol⁻¹ energy drop and the 8-row trajectory table. It also covers the drift overlaps with the low modes (0.875, 0.937 and 0.953) and the 1.79 Å maximum displacement.
- **Frequencies and thermochemistry:** match. This covers the lowest ten modes, the highest mode, n_imag (PCM 1/0, gas 0/0) and the counts below 100 cm⁻¹. It also covers ZPE, H, G, the thermal corrections and S, with 627.5095 kcal mol⁻¹ per Ha.
  - Doubling the FD step changes frequencies by at most 1.886 cm⁻¹ (mean 0.154), ZPE by −0.0059 kcal mol⁻¹ and G by −0.046 kcal mol⁻¹.
  - The N–H stretch, bend and wag assignments were rechecked, with in-plane versus out-of-plane character tested against the cation plane.
  - Also confirmed: the rigid-body fractions and the 66 % cation weight of the 3.2i mode.
- **FMO:** matches. This covers HOMO, LUMO, HOMO−1, LUMO+1, the gap, IP, EA, χ, μ, η, S = 1/(2η) (the dft/fmo.py convention), ω, and the dipole with its vector.
  - The Mulliken shares of HOMO−8 to LUMO+7 and the fragment charges (−1.0001 PCM, −0.8806 gas) were recomputed with the separately coded overlap matrix.
  - That matrix reproduces orthonormality: max |CᵀSC − I| is 3.3e-11 (PCM) and 2.8e-11 (gas).
- **MEP:** matches. All 27 + 27 extrema convert correctly (a.u. × 627.5095 = kcal mol⁻¹ within 0.001), and every nearest-atom site and global site is correct. Every surface statistic matches mwfn_mep.out.
- **NCI spikes and regions:** match.
  - The 15 spikes (bin, minimum RDG, number of points) were recomputed from rdg_scatter.txt with the dft/mwfn.py binning, together with the atom pairs at the lowest-RDG points.
  - The 22 real-space regions (6-connected RDG < 0.5 with at least 20 points) were rebuilt from rdg.cube and sl2r.cube.
- **ELF/LOL:** the ranges, the ring-centroid values (bilinear interpolation on plane.txt) and the O1 offset match.
- **PCM − gas table:** about 40 deltas were recomputed, and all match.
- **Protocol, checked against code and logs:** all match.
  - Grid 4 for btpc_gas (grid 3 would differ by 1.48 µHa) and the Lebedev 29 default in run_case.sh.
  - SCF settings: conv_tol 1e-9 / 1e-8 / 1e-10 and max_cycle 150.
  - geomeTRIC criteria and trust radius; FD step, workers and wall times.
  - Linear-dependency removals: 3 (PCM) and 4 (gas), and 41 × 4 + 30 × 3 during the PCM optimisation.
  - Basis: 606 = 25 × 22 + 8 × 7 functions (290 shells).
  - Multiwfn settings: isovalue 0.001, spacing 0.25 Bohr, RDG_maxrho 0.1, aug2D 6.0, ELFLOL_type 0, iESPcode 2/1.
  - NCI grids with the 1.5 Bohr extension, and the software versions.
- **Figure-data inventory:** every path was opened, and its columns, units, shapes and row counts were checked. The errors found are listed below.

### Corrections (old → new; evidence)
1. **Cube margin** (`protocol.cubes`, Table 3.1)
   - Old: "margin 4.5 Å".
   - New: "margin 4.5 Bohr (2.38 Å) beyond the outermost nuclei".
   - Evidence: every cube origin equals the smallest nuclear coordinate minus 4.5 Bohr, on each axis and in both cases (for example PCM x: −8.1537 − 4.5 = −12.6537 Bohr). PySCF cubegen takes the margin in Bohr. The 0.001 a.u. surface reaches at most 1.8 Å beyond the nuclei (vtx.txt), so it is not clipped. The "4 A" in dft/fmo.py and the "4.5 Å" in the project notes (not deposited) are the same mistake.
2. **Screened NCI points** (inventory row "NCI scatter")
   - Old: "screened points have RDG ≥ 100 (PCM 198,849, gas 193,330 rows)".
   - New: points screened by RDG_maxrho (ρ ≥ 0.1 a.u.) have RDG = 100 exactly: 94,395 rows (PCM) and 103,505 rows (gas). Rows with RDG ≥ 100 total 198,715 (PCM) and 193,237 (gas). The difference, 104,320 (PCM) and 89,732 (gas), is low-density points where RDG exceeds 100 naturally.
   - Evidence: rdg_scatter.txt counted with numpy and awk. The count with RDG = 100 equals the count with |sign(λ₂)ρ| ≥ 0.1.
3. **Atomic masses** (inventory row "Hessian / normal modes"; the text of Section 4.5; `*.freq_recomputed_check`)
   - Old: "PySCF table (H 1.00794, C 12.0107, N 14.0067, O 15.9994)", "frequencies reproduce freq.json to 0.10 cm⁻¹", check values 0.101 / 0.103 cm⁻¹.
   - New: PySCF used H 1.008, C 12.011, N 14.007, O 15.999. With these masses the frequencies reproduce freq.json to 6.5e-11 cm⁻¹ (PCM) and 1.3e-11 cm⁻¹ (gas).
   - Evidence: hessian.npy was projected with each mass set. Standard atomic weights leave a 0.10 cm⁻¹ residual and isotopic masses 0.66 cm⁻¹.
   - Knock-on: the rigid fractions, the imaginary-mode picrate weight and the drift overlaps change only in the 5th–6th significant figure. No displayed value changes.
4. **Molden padding** (Section 2; inventory row "Wavefunction")
   - Old: "the molden files were zero-padded to 606" / "MOs zero-padded to 606".
   - New: only the PCM molden was padded (603 → 606 MOs; the padded energies are max + 1000 a.u.). The gas molden holds 602 MOs and is not padded.
   - Evidence: the files contain 606 and 602 `Sym=` blocks, and Multiwfn reports "actual number of orbitals read" as 606 and 602. The btpc_gas fmo.log has no padding line.
5. **Hessian order** (Section 3, "Frequencies")
   - Old: "PCM authoritative Hessian with step 0.01 Bohr, repeated with 0.005 Bohr".
   - New: the 0.005 Bohr Hessian came first (2026-09-17 13:56–19:42 UTC, freq_step0.005.json). The 0.01 Bohr Hessian was the repeat (2026-09-18 16:05–21:51 UTC) and is the authoritative one (freq.json, hessian.npy).
   - Evidence: STAGES.md, file times, and hessian_method in both JSON files.
6. **Grid evidence for btpc_pcm** (`protocol.grid_level` note; Section 3, "Grid")
   - Old: the "grid 4" in the btpc_pcm opt.xyz comment was cited as evidence.
   - New: withdrawn as evidence. That text is a literal string that dft/finalize_opt.py writes whatever GRID_LEVEL was. Grid level 4 for btpc_pcm now rests on the run_case.sh defaults alone (see O2).
7. **Atom order** (`model.atom_order`; Section 2; inventory row "Structure / labels")
   - Old: "CIF order".
   - New: CIF order within each ion, with the picrate first (atoms 0–17) and then the cation (atoms 18–32).
   - Evidence: the CIF atom_site loop interleaves the two ions (O1, N3, O7, N4, C6, …), whereas the case file's cif_labels list the picrate first.
8. **Optimiser coordinates** (`protocol.optimizer`)
   - Old: "99 internals incl. per-fragment translations/rotations".
   - New: 99 (= 3N) delocalised TRIC coordinates built from 170 primitives: 34 bonds, 35 angles, 17 out-of-plane terms, 72 dihedrals, and a translation and rotation for each ion.
   - Evidence: optimize.log prints "99 internal coordinates being used (instead of 99 Cartesians)" followed by a list of 170 primitives.
9. **ELF/LOL grid steps** (inventory row "ELF / LOL maps")
   - Old: "step 0.029 Å" (picrate plane) and "step 0.026 Å" (cation plane).
   - New: the steps are 0.029 × 0.028 Å (picrate) and 0.026 × 0.028 Å (cation); the gas cation plane spans ±4.23 Å.
   - Evidence: the grid vectors of *_plane.txt and the X/Y ranges in the npz files.
10. **Typo** (Table 4.5 label): "project projection" → "project's own projection".

### Clarifications and additions (not errors)
- **Second picrate** (`model.missing_hbond` note; Section 2): the picrate that H6N binds is the inversion image (3_666) of the picrate in the asymmetric unit. Relative to the picrate in the cluster, it is the 2₁-screw image 2_746 (2−x, y−½, 3/2−z), with H6N···O6 at 2.103 Å (derived from the CIF). Do not call it "inversion-related" when describing the modelled pair.
- **Modes below 100 cm⁻¹** (`*.n_modes_below_100` note; Table 4.5 label): the PCM count of 10 includes the imaginary mode, so 9 are real. The thermochemistry uses only the 92 real modes; the ZPE equals ½Σν over those 92 modes to within 2e-11 Ha.
- **Multiwfn citations** (`software.multiwfn` note): added the NCI citations that Multiwfn prints in logs/mwfn_nci.out: J. Am. Chem. Soc. 132, 6498 (2010), plus the reviews in Comprehensive Computational Chemistry vol. 2 (2024) and Angew. Chem. Int. Ed. e202504895 (2025). Register them through papers/_shared/tools/lit.py before citing.

### Open questions (not changed; the archive cannot confirm them)
- **O1. Thread count.** "4 threads" for both final optimisations (`pcm.opt_wall` value, `gas.opt_wall` note, Table 4.1) is not recorded in any log. It is only the run_case.sh default.
- **O2. btpc_pcm grid and cavity.** Grid level 4 and Lebedev order 29 for btpc_pcm are not verified: verify_grid.log (not deposited) does not cover this case. Confirming them would need an SCF recomputation, which is not authorised in this phase.
- **O3. Auxiliary basis.** The auxiliary basis (`protocol.aux_basis`) is not logged. There is no PySCF source on this machine to check the aug-cc-pVTZ-JKFIT inference.
- **O4. Values taken from the PySCF source.** `protocol.grid_points`, `protocol.pcm_radii`, `protocol.pcm_discretisation` and the libxc ID / VWN-RPA definition of B3LYP could not be re-verified on this machine.
- **O5. Cation name.** The CIF protonation (H on N4 and N6) contradicts the case-file name "1H-benzotriazol-1-ium". Check it against doi:10.1016/j.ijleo.2012.06.018.
- **O6. Pre-optimisation history.** The grid-3 OVH runs and the 6-31G(d) stage are not archived. Section 3 relies only on the preopt.xyz comment lines and the project notes (not deposited) for them.

### Generators
build_facts.py and build_facts_md.py carry the same corrections, and both files were regenerated from them. facts.json changed only in the entries listed above, and derived/*.csv are byte-identical. build_facts_md.py now re-appends this log whenever the files are rebuilt.
