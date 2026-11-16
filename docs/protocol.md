# Protocol

Every setting of the two calculations, as recorded in `analysis/facts.json` (keys `model.*`, `protocol.*` and
`software.*`; the table below is a copy of their values and sources). `scripts/analysis/build_facts.py` reads each
value from the run files or the code that ran; the sources are workspace paths, which `scripts/paths.py` and
`docs/provenance.md` map to this repository.

## Summary

1. **Model.** The 33-atom ion pair (picrate at 1_555, benzotriazolium at 4_666 of CCDC 835950; `structure/README.md`),
   charge 0, singlet, 178 electrons, 606 basis functions.
2. **Optimisation.** B3LYP/6-311++G(d,p), density fitting, PySCF grid level 4; geomeTRIC 1.1.1 in TRIC coordinates
   with Gaussian default thresholds. Gas phase: 3 steps from `preopt.xyz`, all five criteria met. IEF-PCM (water,
   ε = 78.3553, Lebedev order 29, SWIG, scaled modified Bondi radii): 70 steps from `preopt.xyz`, accepted on
   negligible forces (RMS 1.40 × 10⁻⁴, maximum 3.85 × 10⁻⁴ Ha Bohr⁻¹) by `scripts/dft/finalize_opt.py`.
3. **Single points** at the optimised geometries (SCF tolerance 10⁻¹⁰ Ha): orbitals, Mulliken fragment shares, cubes
   (`scripts/dft/fmo.py`).
4. **Hessians** by central finite differences of analytic gradients (198 displaced gradients): gas phase 0.005 Bohr;
   IEF-PCM 0.01 Bohr (reported) and 0.005 Bohr (check, wavenumbers only). Translations and rotations projected out
   (`scripts/dft/vibrations.py`); harmonic, unscaled; thermochemistry at 298.15 K and 1 atm.
5. **Multiwfn 2026.9.13** on `final.molden`: MEP on the 0.001 a.u. density surface (function 12), NCI/RDG on a
   high-quality grid (function 20), ELF and LOL in the picrate and cation ring planes (function 4), with the
   settings in `environment/software.md`.
6. **Facts, tables, figures.** `scripts/analysis/build_facts.py` derives every number of the article from the files
   above (no quantum-chemistry calculation); `scripts/tables/make_tables.py` and `scripts/figures/` read only
   `facts.json` and these files.

Wall times on the 48-vCPU compute VM: gas optimisation 14 min (4 threads), IEF-PCM optimisation about 13.4 h
(4 threads), gas Hessian 149 min (16 workers × 2 threads), IEF-PCM Hessian 346 min (10 workers × 2 threads).

## All recorded settings

| Key | Value | Unit | Source (workspace path) |
|---|---|---|---|
| `model.n_atoms` | 33 |  | molecules/btpc/input/btpc_cluster.xyz |
| `model.n_electrons` | 178 |  | data/r2/molecules/btpc/results/btpc_pcm/logs/mwfn_mep.out |
| `model.charge` | 0 | e | molecules/btpc/cases/btpc_pcm.json:charge |
| `model.multiplicity` | 1 |  | molecules/btpc/cases/btpc_pcm.json:spin |
| `model.fragments` | picrate anion: 18 atoms (atom indices 0-17, C6H2N3O7); benzotriazolium cation: 15 atoms (indices 18-32, C6H6N3) |  | molecules/btpc/cases/btpc_pcm.json:fragments |
| `model.atom_order` | CIF atom order within each ion, picrate first (indices 0-17), then cation (18-32); the CIF itself interleaves the two ions (O1, N3, O7, N4, C6, ...). Cation labels carry the internal suffix #4_666 (never show it in the paper) |  | molecules/btpc/cases/btpc_pcm.json:cif_labels |
| `model.symmetry_picrate` | 1_555 (x, y, z) |  | molecules/btpc/cases/btpc_pcm.json:notes |
| `model.symmetry_cation` | 4_666 (x+1, -y+1/2, z+1/2) |  | molecules/btpc/cases/btpc_pcm.json:cif_labels / dft/build_cluster.py |
| `model.defining_hbond` | N4-H4N...O1 (cation N-H to picrate phenolate O) |  | molecules/btpc/cases/btpc_pcm.json:hbonds |
| `model.missing_hbond` | N6-H6N...O6 (to a second, inversion-related picrate at 1-x,1-y,1-z) is not represented |  | molecules/btpc/cases/btpc_pcm.json:notes |
| `model.n_basis` | 606 |  | data/r2/molecules/btpc/results/btpc_pcm/logs/frequencies.log |
| `model.basis_per_element` | C, N, O: 22 functions each (5s4p1d, spherical d); H: 7 functions (4s1p) |  | data/r2/molecules/btpc/results/btpc_pcm/final.molden:[GTO] |
| `model.basis_overrides` | none |  | molecules/btpc/cases/btpc_pcm.json |
| `model.preopt_basis_functions` | 366 |  | basis-set count: 6-31G(d) with spherical d functions, 25 heavy atoms x 14 + 8 H x 2 |
| `protocol.functional` | B3LYP (PySCF xc="b3lyp" = libxc HYB_GGA_XC_B3LYP, ID 402, VWN-RPA local correlation, same definition as Gaussian) |  | dft/common.py:XC; PySCF 2.14.0 pyscf/dft/libxc.py (XC_CODES "B3LYP": 402, "B3LYPG": 402 "used by Gaussian") |
| `protocol.basis` | 6-311++G(d,p), spherical (5d) functions, all atoms |  | dft/common.py:BASIS; final.molden [5d] |
| `protocol.density_fitting` | yes (RI-JK via mf.density_fit()) |  | dft/common.py:make_mf |
| `protocol.aux_basis` | PySCF automatic default (not logged) |  | dft/common.py:make_mf; PySCF 2.14.0 pyscf/df/addons.py:make_auxbasis |
| `protocol.grid_level` | 4 |  | verify_grid.log (dft/verify_grid.py run, not deposited); run_case.sh (GRID_LEVEL default 4) |
| `protocol.grid_points` | H: 60 radial x 434 angular (Lebedev order 35); C, N, O: 90 radial x 590 angular (order 41); Treutler-Ahlrichs radial, Becke partition, NWChem pruning |  | PySCF 2.14.0 pyscf/dft/gen_grid.py (RAD_GRIDS/ANG_ORDER row for level 4) |
| `protocol.solvation_model` | IEF-PCM |  | dft/common.py:make_mf (with_solvent.method = "IEF-PCM") |
| `protocol.solvent` | water |  | molecules/btpc/cases/btpc_pcm.json:solvent |
| `protocol.solvent_eps` | 78.3553 |  | dft/common.py:SOLVENT_EPS["water"] |
| `protocol.pcm_lebedev_order` | 29 |  | run_case.sh (PCM_LEBEDEV_ORDER default 29) |
| `protocol.pcm_radii` | modified Bondi radii (H 1.10, C 1.70, N 1.55, O 1.52 Å) scaled by 1.2 -> H 1.32, C 2.04, N 1.86, O 1.824 Å |  | PySCF 2.14.0 pyscf/solvent/pcm.py (vdw_scale 1.2, modified_Bondi) and pyscf/data/radii.py (VDW) |
| `protocol.pcm_discretisation` | SWIG (switching-Gaussian) surface, no probe radius, electrostatic term only (no cavitation/dispersion/repulsion) |  | PySCF 2.14.0 pyscf/solvent/pcm.py (surface_discretization_method "SWIG", r_probe 0.0) |
| `protocol.scf_conv_tol_opt_gas` | 1e-09 | Ha | dft/optimize.py (conv_tol 1e-9 gas / 1e-8 solvent) |
| `protocol.scf_conv_tol_opt_pcm` | 1e-08 | Ha | dft/optimize.py |
| `protocol.scf_conv_tol_final` | 1e-10 | Ha | dft/frequencies.py, dft/fmo.py, dft/_fd_worker.py |
| `protocol.scf_max_cycle` | 150 |  | dft/common.py:make_mf |
| `protocol.optimizer` | geomeTRIC 1.1.1 via pyscf.geomopt.geometric_solver, TRIC coordinates: 99 (= 3N) delocalised internals built from 170 primitives (34 bonds, 35 angles, 17 out-of-plane, 72 dihedrals, translation + rotation of each ion), trust radius 0.1 (max 0.3) Å, maxiter 300 |  | data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log (banner and "Optimization Info") |
| `protocol.opt_criteria` | \|dE\| < 1.0e-6 Ha, RMS grad < 3.0e-4, max grad < 4.5e-4 Ha/Bohr, RMS disp < 1.2e-3, max disp < 1.8e-3 Å (all five required) |  | data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log ("Convergence Criteria") |
| `protocol.preopt` | B3LYP/6-31G(d) with loose criteria (dE 1e-6, grms 1.7e-3, gmax 2.5e-3, drms 6.7e-3, dmax 1.0e-2) - pipeline default |  | dft/optimize.py (PREOPT_BASIS, PRE_CONV) |
| `protocol.hessian_method` | central finite differences of analytic gradients (6N = 198 gradients), each displaced SCF started from the converged density |  | dft/frequencies.py; dft/_fd_worker.py |
| `protocol.hessian_projection` | translations and rotations projected out in mass-weighted Cartesians (QR-orthonormalised projector), own routine dft/vibrations.py |  | dft/vibrations.py; freq.json tr_residual_cm-1 |
| `protocol.frequency_scaling` | none (harmonic, unscaled) |  | dft/frequencies.py |
| `protocol.thermo_T` | 298.15 | K | dft/frequencies.py (thermo.thermo(mf, freq, 298.15, 101325)) |
| `protocol.thermo_p` | 101325 | Pa | dft/frequencies.py |
| `protocol.fmo_descriptors` | IP = -E_HOMO; EA = -E_LUMO; chi = (IP+EA)/2; mu = -chi; eta = (IP-EA)/2; S = 1/(2 eta); omega = mu^2/(2 eta) |  | dft/fmo.py |
| `protocol.fragment_composition_method` | Mulliken partition of each MO over the two fragments: sum over fragment AOs of c_mu (S c)_mu |  | dft/fmo.py:frag_share |
| `protocol.cubes` | HOMO/LUMO/density 80^3, ESP 64^3, margin 4.5 Bohr (2.38 Å) beyond the outermost nuclei (PySCF cubegen); render isovalues: orbitals +-0.02 a.u., MEP on 0.001 a.u. density, NCI RDG = 0.5 coloured -0.035..+0.02 a.u. |  | dft/fmo.py; dft/render.py; cube headers |
| `protocol.mep_surface` | Multiwfn main function 12: electron-density isosurface 0.001 a.u., grid spacing 0.25 Bohr, marching tetrahedra; ESP evaluated exactly at surface vertices |  | data/r2/molecules/btpc/results/btpc_pcm/logs/mwfn_mep.out (menu lines "iso: 0.00100", "Spacing of grid points ... 0.250000") |
| `protocol.mep_esp_code_pcm` | iESPcode 2 (LIBRETA fast) |  | data/r2/molecules/btpc/results/btpc_pcm/mwfn_home/settings.ini; logs/mwfn_mep.out |
| `protocol.mep_esp_code_gas` | iESPcode 1 (Multiwfn original ESP code) |  | data/r2/molecules/btpc/results/btpc_gas/mwfn_home/settings.ini; data/r2/logs/btpc_mep_retry.out |
| `protocol.mep_note` | The ESP is that of the (PCM-polarised, for btpc_pcm) molecular charge density alone; the PCM surface charges are not part of the mapped potential |  | dft/mwfn.py (Multiwfn reads only final.molden) |
| `protocol.nci_settings` | Multiwfn 20 -> 1 (NCI), "High-quality grid" covering the whole system with 1.5 Bohr extension; RDG_maxrho = 0.1 a.u. (points with rho >= 0.1 set to RDG 100); promolecular approximation not used (full wavefunction) |  | data/r2/molecules/btpc/results/btpc_pcm/logs/mwfn_nci.out; mwfn_home/settings.ini |
| `protocol.nci_grid_pcm` | 126 x 148 x 95 = 1,771,560 points, spacing 0.158553 Bohr (0.0839 Å) |  | data/r2/molecules/btpc/results/btpc_pcm/logs/mwfn_nci.out; rdg.cube header |
| `protocol.nci_grid_gas` | 125 x 147 x 96 = 1,764,000 points, spacing 0.153832 Bohr (0.0814 Å) |  | data/r2/molecules/btpc/results/btpc_gas/logs/mwfn_nci.out; rdg.cube header |
| `protocol.nci_spike_method` | dft/mwfn.py:nci_post: sign(l2)rho binned from -0.200 to +0.050 a.u. in 0.001 bins; a bin is a candidate if it holds >20 points and its minimum RDG < 0.3; contiguous candidate bins are merged and represented by the bin with the lowest RDG; type label: < -0.015 "attractive (H-bond)", < +0.010 "van der Waals", otherwise "steric repulsion" |  | dft/mwfn.py |
| `protocol.elf_lol_settings` | Multiwfn 4 -> 9 (ELF) / 10 (LOL), plane through three atoms, 300 x 300 points, aug2D = 6.0 Bohr, Becke definition (ELFLOL_type 0, ELF_addminimal 1) |  | data/r2/molecules/btpc/results/btpc_pcm/mwfn_elf_*.in; mwfn_home/settings.ini; logs/mwfn.log |
| `protocol.elf_plane_picrate` | picrate ring plane through C6, C4, C2 (1-based Multiwfn atoms 4, 9, 6) |  | dft/curate_cases.py; data/r2/.../mwfn_elf_picrate_ring.in |
| `protocol.elf_plane_cation` | benzotriazolium plane through N4, N5, C7 (1-based atoms 19, 23, 21) |  | dft/curate_cases.py; data/r2/.../mwfn_elf_benzotriazolium_ring.in |
| `protocol.host` | IBM Cloud mx3d-48x480 VM (48 vCPU / 24 physical Intel Sapphire Rapids cores, 472 GB RAM, local NVMe scratch) |  | IBM Cloud instance profile mx3d-48x480 (public profile specification); run_case.sh comment on the scratch volume; run paths <run-dir> and <scratch-dir> in the logs; install log "threads 48" (not deposited) |
| `software.pyscf` | 2.14.0 |  | data/r2/molecules/btpc/results/btpc_pcm/final.molden header "made by pyscf v[2.14.0]"; density.cube header; install log (not deposited) |
| `software.geometric` | 1.1.1 |  | data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log ("geomeTRIC started. Version: 1.1.1"); install log (not deposited) |
| `software.multiwfn` | 2026.9.13 (official Linux noGUI binary) |  | data/r2/molecules/btpc/results/btpc_pcm/logs/mwfn_mep.out banner; install log (not deposited) |
| `software.paraview` | 6.1.1 (pvbatch, headless) |  | install log (not deposited); environment/software.md (ParaView-6.1.1-MPI-Linux-Python3.12-x86_64) |
| `software.python` | 3.11.16 |  | environment/compute-env.yml (conda package list of the compute environment) |
| `software.numpy_scipy_matplotlib` | NumPy 2.4.6, SciPy 1.17.1, matplotlib 3.11.1 |  | environment/compute-env.yml (conda package list of the compute environment) |
| `software.libxc` | version not recorded |  | environment/compute-env.yml; install log (searched, absent; not deposited) |
| `protocol.opt_threshold_energy` | 1e-06 | Ha | data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log ("Convergence Criteria") |
| `protocol.opt_threshold_grad_rms` | 0.0003 | Ha/Bohr | data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log ("Convergence Criteria") |
| `protocol.opt_threshold_grad_max` | 0.00045 | Ha/Bohr | data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log ("Convergence Criteria") |
| `protocol.opt_threshold_disp_rms` | 0.0012 | Å | data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log ("Convergence Criteria") |
| `protocol.opt_threshold_disp_max` | 0.0018 | Å | data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log ("Convergence Criteria") |
| `protocol.pcm_lebedev_points` | 302 |  | run_case.sh (PCM_LEBEDEV_ORDER 29); PySCF 2.14.0 pyscf/dft/gen_grid.py LEBEDEV_ORDER table |
| `software.figures_env` | Matplotlib 3.11.2; NumPy 2.5.3; SciPy 1.18.1; PyVista 0.49.0; VTK 9.7.0; scikit-image 0.26.0; gemmi 0.7.5; RDKit 2026.3.6 |  | tools/paper-env (installed package metadata, read with importlib.metadata when facts.json is built) |
| `protocol.fd_gradient_grid_response` | no: nuc_grad_method().kernel() with PySCF defaults (grid_response = False), displaced SCF conv_tol 1e-10 |  | dft/_fd_worker.py |
