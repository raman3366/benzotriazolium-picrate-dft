# Revision calculations (November 2026)

During the revision of the article (first round of review, November 2026) further calculations were run on the
benzotriazolium picrate ion pair. They extend the two published calculations (`calculations/btpc_gas/`,
`calculations/btpc_pcm/`) and add a second crystal form of the salt. This page says what was computed, where the
files are, how the jobs were run, and how to run one again. `calculations/revision1/README.md` lists every job.

The folder names end in `1` for the first revision of the article.

## Where the files are

| Folder | Contents |
|---|---|
| `calculations/revision1/<job>/` | Outputs of the 57 revision jobs, as the jobs wrote them (checkpoints left out), and `queue.log`, the job queue's log of the job |
| `calculations/revision1/runner/vm1/`, `vm2/` | Job-queue records of the two compute VMs: `status.json` (every job's expanded command, host, start, end, wall time and exit code), `runner.log`, `runner.out` |
| `calculations/btpcz_gas/`, `calculations/btpcz_pcm/` | Pipeline results of the second crystal form, in the layout of `calculations/btpc_*/` (geometry, wavefunction, frequencies, orbitals, Multiwfn, xz-compressed grids, run records) |
| `cases/btpcz_gas.json`, `cases/btpcz_pcm.json` | Case files of the second crystal form |
| `structure/btpcz/` | Source, checksum and specification of the second form's ion pair, and its starting geometry `btpcz_cluster.xyz` |
| `scripts/revision1/` | The revision scripts, the job lists (`jobs_*.json`), tests, the environment installer `setup_env.sh`, the VM creation script `vm/create_vm.sh`, and the environment locks of two local machines (`locks/`); `CHANGES.md` lists the edits made for deposit |
| `environment/revision1/vm1/`, `vm2/` | Environment of each compute VM: explicit conda lock (`env.lock.txt`), package versions and CPU (`env.versions.txt`), installer log (`setup.log`) |

## What was computed

All calculations use the level of theory of the published calculations unless a job says otherwise
(`scripts/revision1/rev_common.py`): B3LYP/6-311++G(d,p) (PySCF `b3lyp`, spherical d functions), density fitting
with PySCF's default auxiliary basis, integration grid level 4, SCF convergence 10⁻¹⁰ Ha for single points, and for
water IEF-PCM with ε = 78.3553 and Lebedev order 29. Optimisations use geomeTRIC 1.1.1 with Gaussian's default
criteria (all five required unless a job records an acceptance on negligible forces). Molecules: the published ion
pair (`btpc`), its second crystal form (`btpcz`), and for the excited states the free picrate anion and neutral
1H-benzotriazole.

1. **Second crystal form** (P2₁/n, COD 2230422; `structure/btpcz/README.md`). The same protocol as the published
   pair, in the gas phase and in IEF-PCM water: B3LYP/6-31G(d) pre-optimisation from the X-ray pair,
   B3LYP/6-311++G(d,p) optimisation, proton check, frontier orbitals and cubes (`c0_btpcz_gas_core`,
   `c0_btpcz_pcm_core`); finite-difference Hessians with a 0.005 Bohr step (gas phase in `c0_btpcz_gas_core`;
   IEF-PCM in `c0_btpcz_pcm_fdhess` with the helper workers `fdh1`–`fdh5`) and, before the full IEF-PCM Hessian,
   a check along the ten softest modes (`c0_btpcz_pcm_softmode`); Multiwfn analyses (`c0_btpcz_*_analysis`);
   Hirshfeld and CM5 charges (`charges_btpcz`); the single points of point 2 below (`headline_btpcz`).
2. **Where the LUMO lies.** Single points at both published minima in the gas phase and in IEF-PCM with
   ε = 2, 4, 8, 20 and 78.3553, with the cation and picrate shares of the frontier orbitals and the permittivity at
   which the LUMO moves from the cation to the picrate (`headline_gasgeom`, `headline_pcmgeom`); vertical ΔSCF
   states N+1 and N−1 (UKS doublets) in both models at both minima (`deltascf_gasgeom`, `deltascf_pcmgeom`);
   Hirshfeld and CM5 charges of the published wavefunctions (`charges_stored`); PBE0 and CAM-B3LYP at both minima
   (`sp8_pbe0`, `sp8_camb3lyp`).
3. **Dispersion.** B3LYP-D3(BJ) re-optimisations in the gas phase and in IEF-PCM, seeded with the published
   Hessians (`opt_d3_gas`, `opt_d3_pcm`), and Hessian checks at the new minima along the ten softest modes
   (`softmode_d3_gas`, `softmode_d3_pcm`).
4. **Proton transfer along N4–H4N···O1.** Relaxed scans of r(H4N···O1) in the gas phase, in three chains from the
   minimum and from the neutral form, plus two points on the oxygen side (`scan_gas_up`, `scan_gas_down`,
   `scan_gas_neutral`, `scan_gas_oside`); unconstrained optimisations started with the proton on O1
   (`neutral_start_gas`, `neutral_start_pcm`); single points with ωB97X-D3(BJ), M06-2X, PBE0, CAM-B3LYP and
   DF-MP2/SCS-MP2 at the scan geometries (`sp_scan_*`, `sp_oside_all`) and at the IEF-PCM salt and neutral minima
   (`sp_min_pcm`, `sp_neutral_pcm`, `sp_gas_pcmgeoms`); one-dimensional proton vibrational levels on every profile
   (sinc-DVR, `dvr_gas_b3lyp`, `dvr_gas_all`); and a coupled-cluster correction to the DF-MP2 salt–neutral energy
   difference, E[CCSD(T)] ≈ E[DF-MP2/6-311++G(d,p)] + (E_corr[CCSD(T)] − E_corr[MP2]) in a smaller basis (6-31G\*
   or 6-31+G\*\*) and in frozen natural orbitals keeping 50, 65 or 80 % of the virtual orbitals, in IEF-PCM water
   with the reaction field frozen at the HF density and in the gas phase (`cc_*`).
5. **Inter-ion torsion.** Relaxed scans of the C1–O1···N4–C12 torsion from −180° to −60° in 15° steps in the gas
   phase, with B3LYP (`scan_tors_b3lyp`) and B3LYP-D3(BJ) (`scan_tors_d3`, `scan_tors_d3_up`,
   `scan_tors_d3_up60`, three parts of one scan).
6. **Excited states.** TD-DFT in the Tamm–Dancoff approximation, 20 singlets, with B3LYP and CAM-B3LYP: the ion pair
   at both published minima, in IEF-PCM with non-equilibrium linear-response solvation (`td_pair_b3lyp`,
   `td_pair_camb3lyp`), and the free picrate anion and neutral 1H-benzotriazole optimised in IEF-PCM water
   (`td_free`); hole and electron fractions on the two ions.

Each script's docstring gives the method and its options in full. The job lists (`scripts/revision1/jobs_vm.json`,
`jobs_vm2.json`) give each job's command and, in the `item`, `note` and `note_*` fields, its purpose and every change
made while the queue ran. Comments that mention plan items,
review findings and internal notes refer to working files of the project that are not part of this deposit.

## Files of a job

Every job folder holds `queue.log` (the queue's log: the expanded command and the job's output). The others hold
`result.json` and `log.txt` as well, except `c0_prep`, `c0_btpcz_gas_analysis`, `c0_btpcz_pcm_analysis` and
`fdh1`–`fdh5`, whose output went into the case folders or into `c0_btpcz_pcm_fdhess/fd/`:

- `result.json`: the job's record. Command line, host (CPU, cores, memory), threads, memory limit, software
  versions and conda builds, the SHA-256 of the scripts that ran, the input files with their SHA-256, the settings
  (basis, functional, grid, Lebedev order, convergence), every item computed (energies in Ha, orbital energies in
  eV, fragment shares and charges, geometries, timings in s) and a summary.
- `log.txt`: the script's own log, with times in IST (UTC+5:30).
- PySCF output: `pyscf/<item>.log` or `pyscf_*.log` in the job folder.
- Geometries (Å, XYZ): `points/` (final geometry of each scan point or optimisation), `traj/` (every optimisation
  step), `geoms/` (the geometries a single-point job used), `free/` (free species of `td_free`), `scan.xyz`,
  `start_*.xyz`, `opt.xyz`, `opt_traj.xyz`.
- Tables and plots made by the scripts: `profile.csv` (scans), `table.csv` and `table.md` (single points, relative
  energies in kcal/mol), `levels.csv` and `*.png` (DVR); `constraints_p<angle>.txt` (the geomeTRIC constraint of
  each torsion-scan point).
- `final.molden` (`opt_d3_gas`, `opt_d3_pcm`): the wavefunction at the new minimum, padded for Multiwfn.
- Second crystal form: `preopt/` and `opt/` (the optimisation records), `fd/` (the finite-difference Hessian in
  chunks of three Cartesian columns, `cols_XXX_YYY.npz`, with a manifest; `*.claim` files name the helper worker
  that computed a chunk), `softmode/`, `logs/` (logs of the pipeline stages), `mwfn_gas/`, `mwfn_pcm/` (Multiwfn
  output of the charge jobs). The case results themselves are in `calculations/btpcz_gas/` and
  `calculations/btpcz_pcm/`.

Not deposited: the `chk/` folders of the revision jobs (222 PySCF checkpoints and 20 geomeTRIC constraint files,
746 MiB; the results the jobs took from them are in `result.json`, and they can be recomputed), scratch files, the
per-job exit-code files (the exit codes are in `runner/<vm>/status.json`), a copy of one input geometry that was
placed on the second VM (identical to `sp_neutral_pcm/geoms/neutral_pcm.xyz`), and, for the second crystal form,
`opt.chk` and the pipeline pictures, as for the published cases.

## Paths in the files

The jobs ran in a working folder on each compute VM that held the project copy (`repo/`), the environment
(`env/`), the queue (`queue/`) and scratch space (`scratch/`). In the deposited text files that folder is written
`<run-dir>`, the VM host names `<vm-1>` and `<vm-2>`, and the VM user `<user>`; nothing else was changed
(`docs/provenance.md` lists the rules and counts for the data, `scripts/revision1/CHANGES.md` for the scripts).
Paths relative to the project copy map to this repository as in `docs/provenance.md`, and further:

| Path in the files | Here |
|---|---|
| `dft/revision/<file>` | `scripts/revision1/<file>` |
| `<run-dir>/queue/out/<job>/` | `calculations/revision1/<job>/` (both VMs) |
| `<run-dir>/queue/logs/<job>.log` | `calculations/revision1/<job>/queue.log` |
| `<run-dir>/queue/status.json`, `runner.log`, `runner.out` | `calculations/revision1/runner/<vm>/` |
| `molecules/btpcz/cases/<case>.json` | `cases/<case>.json` |
| `molecules/btpcz/input/btpcz_cluster.xyz` | `structure/btpcz/btpcz_cluster.xyz` |
| `molecules/btpcz/input/2230422.cif` | not deposited (`structure/btpcz/README.md`) |
| `molecules/btpcz/results/<case>/<file>` | `calculations/<case>/<group>/<file>`, grouped as for `btpc_*` (`docs/provenance.md`); `softmode_check.json` is in `frequencies/` |

The PySCF checkpoints `calculations/btpcz_*/wavefunction/final.chk` were treated as those of the published cases
(rule L3 in `docs/provenance.md`): the path of the `opt.xyz` stored in the `mol` record was replaced by the
project-relative path `molecules/btpcz/results/<case>/opt.xyz`, padded to the same length.

## How the jobs were run

- **Environment.** Python 3.11.16, PySCF 2.14.0 (conda-forge build `py311h990ede3_1302`, the x86-64 level-3 build of
  the release used for the published calculations), geomeTRIC 1.1.1, NumPy 2.4.6, SciPy 1.17.1, h5py 3.16.0,
  pyscf-dispersion 1.5.0 (D3(BJ)), libxc 7.1.2, OpenBLAS 0.3.34, installed with micromamba 2.10.0 by
  `scripts/revision1/setup_env.sh`. `environment/revision1/vm1/env.lock.txt` is the explicit lock of what was
  installed; the second VM was installed from that lock and its lock is identical. matplotlib 3.11.2 was added to
  the first VM's environment on 22 Nov 2026, after the lock was written, for the pipeline plots of the second form
  and the DVR figures (its version is in `result.json` of `dvr_gas_b3lyp` and `dvr_gas_all`); no number depends on
  it. The published calculations used the level-4 build of the same PySCF release (`py311hc37a911_1402`);
  `headline_gasgeom` and `headline_pcmgeom` reproduce the published total energies to 4 × 10⁻¹⁰ Ha (gas phase) and 1.1 × 10⁻⁸ Ha (IEF-PCM).
- **Multiwfn** 2026.11.13 (official Linux noGUI binary, `Multiwfn_2026.11.13_bin_Linux_noGUI.zip`, SHA-256
  `1a7e610374d305c5470f8e5d11a3017a81a8c4fac356c4208c5da65fa4c98b20`; `environment/revision1/vm1/setup.log`) was
  used for the second crystal form and the charge jobs, because version 2026.9.13 of the published calculations was
  no longer offered. Settings as for the published runs (`aug2D` 6.0, `RDG_maxrho` 0.1), `nthreads` 4, and
  `iESPcode` 2 (LIBRETA) for both models of the second form (`calculations/btpcz_*/multiwfn/settings.ini`).
- **Hosts.** Two cloud VMs of the same type: 32 vCPU (AMD EPYC 7B13), 128 GB memory, Ubuntu 24.04. The first VM
  ran 55 jobs from 21 Nov 2026, 03:10 am, to 23 Nov 2026, 06:07 am (IST), up to eight at a time, mostly with 4
  threads each; the second ran the two `*_f080` coupled-cluster jobs on 23 Nov 2026, 08:10 am to 01:18 pm, with 16
  threads each. The environment was first installed and checked on two local machines (a 4-core x86-64 laptop and
  a 4-core ARM64 board); their locks are in `scripts/revision1/locks/`. ParaView was not installed on the VMs, so the
  3-D renders of the pipeline were not made for the second form.
- **Job queue.** `scripts/revision1/jobqueue.py` (one machine, no daemon) ran the job lists
  `scripts/revision1/jobs_vm.json` (first VM) and `jobs_vm2.json` (second VM) from the project copy:

  ```sh
  <run-dir>/env/bin/python dft/revision/jobqueue.py run --jobs dft/revision/jobs_vm.json --root <run-dir>/queue \
      --slots 8 --threads 36
  ```

  (the slots and thread budget recorded under `runner` in `runner/vm1/status.json`). It starts a job when its
  dependencies have finished and its threads fit the budget, exports `THREADS` and pins every BLAS and OpenMP pool
  to it, gives it a scratch folder, and records it in `status.json`. Each script writes `result.json` as it goes and
  skips finished items when restarted. All 57 jobs are recorded as finished with exit code 0; the `note_fix`,
  `note_split` and `note_threads` fields of the job lists record the jobs that were changed and restarted while the
  queue ran. `jobs_main.json`, `jobs_proton.json` and `jobs_piserver.json` are the earlier plans for the two
  local machines; `jobs_vm.json` merged them (the `source_file` field of each job) and is the list that ran.
- **Script versions.** Every `result.json` records the SHA-256 of the scripts that produced it. They match the
  deposited scripts, with two kinds of exception. (1) `rev_common.py`, `setup_env.sh` and a few other files were
  edited for deposit in comments and defaults only; `scripts/revision1/CHANGES.md` gives each file's SHA-256 before
  and after the edit (for `rev_common.py`, the value before is the one in `result.json`). (2) The four jobs `sp_scan_wb97xd3bj`, `sp_scan_m062x`,
  `sp_scan_pbe0` and `sp_scan_camb3lyp` ran with an earlier revision of `singlepoints.py` (SHA-256
  `320ef254e22886e2ba446ece761f953cc1a75bb8e2748c45f929a2ccc526d114`). The deposited version differs in the DF-MP2
  branch, which was made energy-only on 22 Nov 2026 (`note_fix` of `sp_scan_dfmp2` and `sp_min_pcm` in
  `jobs_vm.json`); those four jobs do not run DF-MP2. The earlier revision is not deposited.
- **Second crystal form, run records.** `calculations/btpcz_*/run/RESULTS.md` is the pipeline's generated report.
  It was written before the IEF-PCM Hessian had finished, so the IEF-PCM report has no frequency table, and its
  methods paragraph repeats the template's "grid level 3" and "computed analytically", which do not describe these
  runs (grid level 4, finite-difference Hessians; `result.json` of the `c0_*` jobs and `run/STAGES.md` record what
  was done). `frequencies/freq.json` and `hessian.npy` of `btpcz_pcm` come from `c0_btpcz_pcm_fdhess`.
  For `btpcz_gas`, the harmonic analysis with translations and rotations projected out (`scripts/dft/vibrations.py`)
  gives one imaginary mode, −59.1 cm⁻¹, while PySCF's own analysis of the same Hessian gives none
  (`freq.json`: `frequencies_cm-1` and `frequencies_pyscf_cm-1`; `run/STAGES.md`); the pipeline keeps the projected
  set. The `btpcz_pcm` Hessian has no imaginary mode (lowest 5.5 cm⁻¹).

## Facts, tables and figures

The results of the revision jobs enter the facts register through section 10 of `scripts/analysis/build_facts.py`:
the modules `scripts/analysis/rev/proton.py`, `structure.py`, `electronic.py`, `experimental.py` and `model.py`
read the job records (`calculations/revision1/<job>/result.json` and the files next to them, under their project
paths `data/revision_runs/<queue>/out/<job>/`), check them against one another and against the published register,
and write the `rev.*` keys of `analysis/facts.json` (1483 keys) and the `analysis/derived/rev_*.csv` tables. The
proton DVR levels are recomputed with `scripts/revision1/dvr.py`. The eleven revision tables
(`tables/tbl_rev_*.md`, `tables/tblS_rev_*.md`) are made from the register by `scripts/tables/rev_*_tables.py` and
`fix_widths.py`, and the revision figures by `scripts/figures/fig_rev_eps.py`, `fig_rev_proton.py`,
`fig_rev_torsion.py` and `graphical_abstract.py`. `scripts/reproduce.sh` rebuilds the register and the tables and
checks that nothing changes.

Some job records hold input checksums of files that were changed for deposit (paths written as `<run-dir>`);
`docs/scrubbed_files.tsv` gives the checksum of each such file as the run wrote it, and the fact modules use it when
they compare a file with a record.

## Running one job again

The scripts expect the layout of the project copy: they take the folder two levels above
`dft/revision/<script>.py` as the project root and read case files from `molecules/<id>/cases/`, the published
minima from `molecules/btpc/results/` and the published wavefunctions and Hessians from
`data/r2/molecules/btpc/results/`. A small tree of links in that layout is enough. For example, for
`headline_gasgeom` (Linux x86-64; about 1.2 h with 4 threads on the VM):

```sh
# 1. the environment of record (x86-64; on ARM64 use scripts/revision1/locks/linux-aarch64_piserver.lock.txt)
bash scripts/revision1/setup_env.sh "$HOME/btpc-env" --lock environment/revision1/vm1/env.lock.txt

# 2. a run tree in the project layout, outside the repository
REPO=$PWD; W=$(mktemp -d)
mkdir -p "$W/dft/revision" "$W/molecules/btpc/cases" "$W/molecules/btpc/results/btpc_gas" \
         "$W/data/r2/molecules/btpc/results/btpc_gas"
ln -s "$REPO"/scripts/dft/*.py "$W/dft/"                  # the published pipeline modules (dft/common.py: make_mf)
cp "$REPO"/scripts/revision1/*.py "$W/dft/revision/"      # copies: the scripts locate the project root from their own path
ln -s "$REPO/cases/btpc_gas.json" "$W/molecules/btpc/cases/"
ln -s "$REPO/calculations/btpc_gas/geometry/opt.xyz" "$W/molecules/btpc/results/btpc_gas/"
ln -s "$REPO/calculations/btpc_gas/wavefunction/final.chk" "$W/data/r2/molecules/btpc/results/btpc_gas/"

# 3. the job's command from scripts/revision1/jobs_vm.json, with {python}, {root} and {id} filled in
cd "$W"
THREADS=4 MAX_MEMORY_MB=8000 "$HOME/btpc-env/env/bin/python" dft/revision/headline.py \
    --job headline_gasgeom --out "$W/out/headline_gasgeom" --geom gas --refine 1
```

The job checks itself: its gas-phase point must reproduce the published energy (−1317.20823094 Ha) to 10⁻⁶ Ha.
Compare `out/headline_gasgeom/result.json` with `calculations/revision1/headline_gasgeom/result.json`. Other jobs
need the inputs named in their command and listed under `inputs` in their `result.json` (with SHA-256); a job that
reads another job's output (`{root}/out/<job>`) can be pointed at `calculations/revision1/<job>/`. The jobs of the
second crystal form (`c0_*`, `fdh*`) also run the published pipeline through `run_case.sh` in the project layout
(`molecules/btpcz/results/<case>/`), which the pipeline deposited in `scripts/` (adapted to this repository's
layout, `docs/provenance.md`) does not provide; recompute them in a copy of the project layout. This recipe was
written from the scripts and the records; it was not run while the deposit was being built.
