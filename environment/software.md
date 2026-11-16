# Software

## Programs that produced the data

| Program | Version | Used for | Where the version is recorded |
|---|---|---|---|
| PySCF | 2.14.0 (conda-forge) | SCF, gradients, IEF-PCM, finite-difference Hessians, cubes | Molden header of `final.molden`; cube headers |
| geomeTRIC | 1.1.1 | geometry optimisation (TRIC coordinates) | banner of `run/logs/optimize.log` |
| Python | 3.11.16 | driver of the above | conda package list of the run environment |
| NumPy / SciPy / matplotlib | 2.4.6 / 1.17.1 / 3.11.1 | numerics, pipeline plots | conda package list |
| libxc | not recorded | B3LYP (libxc ID 402, VWN-RPA local correlation) inside PySCF | `analysis/facts.json` key `software.libxc` |
| Multiwfn | 2026.9.13, official Linux noGUI binary | MEP surface analysis, NCI (RDG), ELF and LOL planes | banner of `run/logs/mwfn_*.out` |
| ParaView | 6.1.1 (ParaView-6.1.1-MPI-Linux-Python3.12-x86_64, pvbatch, headless) | 3-D isosurface renders of the pipeline (not deposited) | `run/logs/render.log` |

`environment/compute-env.yml` pins the Python packages. `scripts/setup_linux.sh` installs them, Multiwfn and ParaView
into `tools/` (not part of the dataset).

## Multiwfn

- Download: `Multiwfn_2026.9.13_bin_Linux_noGUI.zip` from the Multiwfn website (http://sobereva.com/multiwfn/).
  The binary is not redistributed here. SHA-256 of the `Multiwfn_noGUI` binary that produced the data:

  ```
  bff6c2621a4454f38a4f5c6ced09c78ac80ad1a82d6d7c7429ec258b3fe8460d  Multiwfn_noGUI
  ```

  If this exact build is no longer offered, the grids it produced are deposited (`calculations/<case>/grids/`), so
  every number derived from them can still be checked.
- `settings.ini` of each run is deposited (`calculations/<case>/multiwfn/settings.ini`). Four parameters were set
  for these runs (`scripts/setup_linux.sh` sets the first three, `scripts/dft/mwfn.py` the fourth):

  | Parameter | Gas phase | IEF-PCM | Why |
  |---|---|---|---|
  | `nthreads` | 12 | 12 | thread count |
  | `aug2D` | 6.0 Bohr | 6.0 Bohr | wider plane maps for ELF and LOL |
  | `RDG_maxrho` | 0.1 a.u. | 0.1 a.u. | raised from 0.05 so that the short N4–H4N···O1 bond stays in the NCI map (its critical-point density is 0.087 a.u. in the gas phase) |
  | `iESPcode` | **1** (original ESP code) | 2 (LIBRETA, fast; the Multiwfn default) | LIBRETA (codes 2 and 3) crashed on the gas-phase wavefunction (exit 174, `calculations/btpc_gas/run/btpc_mep_retry.out`) |

- Licence: Multiwfn is free and may be redistributed; its licence requires citing T. Lu, F. Chen,
  *J. Comput. Chem.* 33 (2012) 580–592 and T. Lu, *J. Chem. Phys.* 161 (2024) 082503.

## Environment variables that set the level of theory

`scripts/dft/common.py` defaults to integration grid level 3 and Lebedev order 17 for the PCM cavity. Every reported
run used **grid level 4 and Lebedev order 29**, which `scripts/run_case.sh` sets (`GRID_LEVEL=4`,
`PCM_LEBEDEV_ORDER=29`). Running `scripts/dft/pipeline.py` directly therefore gives a different level of theory.
The IEF-PCM Hessian in `freq.json`/`hessian.npy` used `FD_STEP=0.01` Bohr; the code default, used for the gas
phase and for the IEF-PCM check in `freq_step0.005.json`, is 0.005 Bohr.

## Analysis, tables and figures

| Program | Version |
|---|---|
| Python | 3.12.3 |
| NumPy, SciPy, matplotlib | 2.5.3, 1.18.1, 3.11.2 |
| PyVista, VTK | 0.49.0, 9.7.0 |
| gemmi, RDKit, scikit-image, pandas, Pillow | 0.7.5, 2026.3.6, 0.26.0, 3.0.6, 12.3.0 |
| pandoc (only for `make_tables.py --check`) | 3.11; the check needs 3.11 or newer (older versions, such as 3.1.3, do not write the table ids it looks for) |
| Liberation Sans (figure font; metric-compatible with Arial) | system package `fonts-liberation` |

The full list is `environment/paper-requirements.txt`. The figure scripts look for the font files in
`tools/fonts/extracted/usr/share/fonts/truetype/liberation/`; without them they fall back to DejaVu Sans, whose wider
glyphs make the strict layout checks of `figkit.qa` fail for some figures (docs/REPRODUCE.md).
