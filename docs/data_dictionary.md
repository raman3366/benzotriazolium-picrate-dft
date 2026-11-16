# Data dictionary

## Conventions

- **Cases.** `btpc_gas` (gas phase) and `btpc_pcm` (IEF-PCM, water). The folder `calculations/<case>/` holds one case.
- **Atoms.** 33 atoms in the order of `cif_labels` in `cases/<case>.json`: picrate O1 … H3 (indices 0–17), then the
  benzotriazolium cation (18–32). Labels ending in `#4_666` belong to the symmetry-generated cation; the article
  drops the suffix. JSON files and the Python scripts use 0-based indices; Multiwfn inputs (`mwfn_*.in`) use 1-based
  indices.
- **Units.** Energies in hartree (Ha) unless the key says eV or kcal (`*_eV`, `*_kcal`); coordinates in Å in `.xyz`,
  `.csv` and Multiwfn text outputs; cube files in Bohr (Gaussian cube convention); Hessian in Ha Bohr⁻²; wavenumbers
  in cm⁻¹ (an imaginary mode is stored as a negative number); dipoles in debye.
- **Text files** are UTF-8 with LF line ends. Grids are xz-compressed; `scripts/unpack_grids.sh` unpacks them and
  checks them against `docs/grids.sha256`.

## cases/

`btpc_gas.json`, `btpc_pcm.json`: the input of each calculation. `xyz` (starting geometry, workspace path of the
rebuilt X-ray pair), `charge`, `spin`, `source_cif`, `cif_labels`, `fragments` (atom indices per ion), `hbonds` and
`distances` (watched contacts), `planes` (three atoms per ELF/LOL plane), `cell` (Å, degrees), `proton` (H, donor,
acceptor indices), `chemical_name`, `formula`, `ccdc`, `source_doi`, `notes`, and `solvent` (PCM case only).

## calculations/<case>/geometry/

| File | Contents |
|---|---|
| `preopt.xyz` | Starting geometry of the final optimisation (Å); the comment line says where it came from (see README, Provenance) |
| `opt.xyz` | Optimised geometry (Å); the comment line gives the level of theory and the final energy |
| `opt_traj.xyz` | Every optimisation step (multi-frame xyz, Å); comment lines `step k  E = … Ha`. Gas phase 3 frames, IEF-PCM 71 |
| `opt_summary.json` | `E_final_Ha`, `steps`, `wall_min` (NaN for IEF-PCM; the wall time is in `run/STAGES.md`), `n_basis`, `bond_rmsd_A`, `bond_maxdev_A` (against the X-ray pair), `converged_on` (IEF-PCM only), `distances_A`, `hbonds` |
| `geometry_vs_start.csv` | Every covalent bond: `bond`, `start_A` (X-ray), `dft_A`, `diff_A` |
| `proton_check.json` | Whether H4N moved from N4 to O1: `d_H_N_A`, `d_H_O_A`, `proton_transferred`, `bonds_broken`, `bonds_formed`, `connectivity_changed` |

## calculations/<case>/wavefunction/

| File | Contents |
|---|---|
| `final.molden` | Kohn–Sham orbitals at the optimised geometry (Molden format, spherical d, written by PySCF 2.14.0). PySCF removed near-linear dependencies (gas phase 4, IEF-PCM 3 of 606 functions), leaving 602 and 603 orbitals; the IEF-PCM file is padded to 606 with zero vectors for Multiwfn |
| `final.chk` | PySCF checkpoint (HDF5): `mol` (the molecule as PySCF JSON) and `scf/e_tot`, `scf/mo_coeff`, `scf/mo_energy`, `scf/mo_occ` (Ha). Load with `pyscf.lib.chkfile.load_mol` and `load('final.chk', 'scf')`. One text field of `mol` was edited for deposit (docs/provenance.md, rule L3) |

## calculations/<case>/frequencies/

| File | Contents |
|---|---|
| `hessian.npy` | Cartesian Hessian, shape (33, 33, 3, 3), Ha Bohr⁻², unprojected and symmetrised; `H.transpose(0, 2, 1, 3).reshape(99, 99)` gives the 99 × 99 matrix |
| `freq.json` | `frequencies_cm-1` (99 − 6 = 93 modes after projection; imaginary as negative), `n_imaginary`, `hessian_method` (step and workers), `hessian_wall_min`, `frequencies_pyscf_cm-1` and `n_imaginary_pyscf` (PySCF's own harmonic analysis, for comparison), `tr_residual_cm-1`, `ZPE_Ha`, `H298_Ha`, `G298_Ha` |
| `freq_step0.005.json` | IEF-PCM only: the same analysis with a 0.005 Bohr step (the Hessian itself is not archived) |

## calculations/<case>/orbitals/

| File | Contents |
|---|---|
| `fmo.json` | `E_total_Ha`; HOMO, LUMO, gap and HOMO−1/LUMO+1 (eV); Koopmans descriptors `IP_eV`, `EA_eV`, `electronegativity_chi_eV`, `chemical_potential_mu_eV`, `hardness_eta_eV`, `softness_S_eV-1` (S = 1/(2η)), `electrophilicity_omega_eV`; `dipole_Debye` and `dipole_vector_Debye`; `n_basis`, `n_electrons`, `n_occ`; `mo_energies_eV` (all orbitals); `fragment_composition` (Mulliken shares of the frontier orbitals) |
| `fmo_table.md` | The frontier orbitals as a small Markdown table |

## calculations/<case>/multiwfn/

| File | Contents |
|---|---|
| `mwfn_*.in` | Keyboard input given to Multiwfn for each analysis (`mep`, `nci`, `elf_*`, `lol_*`) |
| `settings.ini` | Multiwfn settings of the run (environment/software.md lists the four changed parameters) |
| `surfanalysis.txt` | Multiwfn list of the MEP surface minima and maxima (a.u., eV, kcal mol⁻¹, coordinates in Å) |
| `vtx.txt` | Every vertex of the 0.001 a.u. density surface: count, then x, y, z (Å) and ESP (a.u., eV, kcal mol⁻¹) |
| `mep_extrema.csv` | The extrema with `type`, `index`, `global`, ESP in a.u., eV and kcal mol⁻¹, position (Å) and nearest atoms |
| `mep_range.json` | `vmin`, `vmax` (a.u.) of the surface ESP |
| `mep_surface.json` | Surface statistics: area, volume, positive and negative areas (Å², Å³), averages, σ²tot, ν, Π, MPI, polar area, numbers of extrema |
| `nci_spikes.csv` | NCI interaction regions from the RDG scatter: `sign_l2_rho_au`, `min_RDG`, `n_points`, `type` (method: `protocol.nci_spike_method` in facts.json) |
| `rdg_scatter.npz` | The plotted part of the RDG scatter: `x` = sign(λ₂)ρ (a.u.), `y` = RDG, for −0.25 < x < 0.06 and y < 2.5; `n_total` = all points |
| `{elf,lol}_{picrate,benzotriazolium}_ring.npz` | ELF or LOL in a ring plane on 300 × 300 points: `X`, `Y` (in-plane coordinates, Å, origin `center`), `F` (function value), `atoms_xy` (atoms projected on the plane, Å), `in_plane` (atoms within 0.6 Å of the plane), `normal`, `center` (Å) |

## calculations/<case>/grids/ (xz-compressed)

| File | Contents |
|---|---|
| `density.cube`, `homo.cube`, `lumo.cube` | Electron density and the frontier orbitals, 80³ points, Bohr (PySCF cubegen, margin 4.5 Bohr) |
| `esp.cube` | Electrostatic potential, 64³ points (used only for rendering) |
| `rdg.cube`, `sl2r.cube` | RDG and sign(λ₂)ρ on the Multiwfn NCI grid |
| `rdg_scatter.txt` | Every NCI grid point: x, y, z (Å; the code comment in `scripts/dft/mwfn.py` says Bohr, which is wrong), sign(λ₂)ρ (a.u.), RDG |
| `*_plane.txt` | Multiwfn plane output behind the `.npz` maps: X, Y, Z (Å), two in-plane plot coordinates, function value |

## calculations/<case>/run/

| File | Contents |
|---|---|
| `RESULTS.md` | Report written by the pipeline at the end of the run. Known errors are listed under `inconsistency.*` in `analysis/facts.json` (for example it says integration grid level 3; the runs used level 4) |
| `STAGES.md` | Time-stamped notes of every pipeline stage (UTC) |
| `status.json` | Machine-readable stage status (stale for the IEF-PCM optimisation, see README) |
| `logs/` | Output of each stage: `optimize.log` (geomeTRIC), `frequencies.log` and `fd_worker*.log`, `fmo.log`, `mwfn.log` and `mwfn_*.out` (Multiwfn), `plots.log`, `render.log`, `report.log` |
| `btpc_*.out` | Driver logs of `run_case.sh` (threads, memory, stage times); `btpc_mep_retry.out` records the gas-phase MEP retry |

## analysis/

| File | Contents |
|---|---|
| `facts.json` | 1007 entries `key: {value, unit, case, source, note}` (some also `expr`, the arithmetic on other keys, and `added`, the revision round). `source` names the file (workspace path) and field the value was read from |
| `facts.md` | The same register as a readable dossier, made by `scripts/analysis/build_facts_md.py`; its last section is a hand-written fact-check log kept across rebuilds |
| `derived/mo_fragment_composition_<case>.csv` | Every orbital: `index0`, `label`, `energy_eV`, `occupation`, `picrate_pct`, `cation_pct` (Mulliken) |
| `derived/crystal_hbonds.csv` | Hydrogen bonds of the crystal computed from the CIF by symmetry: donor, H, acceptor, symmetry code and operation, D–H, H···A, D···A (Å), angle (°), kind |
| `literature_claims.jsonl` | Five literature records (DOI, bibliographic data, `claims`) from which the scripts parse the pKa values and earlier H···O distances quoted in the article |

## tables/

`tbl_*.md` (main text) and `tblS_*.md` (SI): one pandoc pipe table each, with the caption line
`Table: … {#tbl:<stem>}` after it. `tbl*.csv` are plain copies (scripts/tables/tables_to_csv.py). `INDEX.md` lists, per
table, the facts keys, the other input files and the cross-checks made before writing.
