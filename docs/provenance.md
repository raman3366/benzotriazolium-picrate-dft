# Provenance

## Where the files come from

The calculations ran in September 2026 on a 48-vCPU Intel Sapphire Rapids VM (IBM Cloud profile mx3d-48x480), in a
project workspace that also held other studies. Gas phase: optimisation, Hessian and analyses 16–17 Sep; IEF-PCM:
optimisation 16–17 Sep, Hessian with 0.005 Bohr 17 Sep, Hessian with 0.01 Bohr 18 Sep (UTC; `run/STAGES.md` and
`run/btpc_*.out` have the times). The run folders were archived unchanged, and this repository was built from that
archive on 16 Nov 2026 by copying an explicit list of files (no version history of the workspace is included). Where
the workspace kept a second, tracked copy of a result file, the two copies were compared and are byte-identical.

The analysis register `analysis/facts.json`, its dossier `facts.md`, the derived CSV files and the tables were then
rebuilt inside this repository by the deposited scripts (docs/REPRODUCE.md).

## What is not deposited, and why

| Item | Why |
|---|---|
| The CIF of CCDC 835950 and the X-ray ion pair `btpc_cluster.xyz` | CCDC licence; obtain the CIF from the CCDC and rebuild the pair (`structure/README.md`, checksums given) |
| ParaView VTK copies of the cubes (`vtk/`, 135 MiB) | format copies of the deposited cubes, rebuilt by `scripts/dft/render.py` |
| `opt.chk` (last optimisation SCF) | superseded by `final.chk` |
| Pipeline pictures (`fig/*.png`) | superseded by the figures of the article (added at v1.0.0) |
| The September 2026 hand-over package and its unpacked copy | contains the CIF, the Multiwfn binary and internal notes |
| Multiwfn and ParaView binaries, the conda environments, fonts, pandoc | third-party software; versions, download names and checksums in `environment/` |
| `scripts/dft` modules for fleet management, progress, packaging, reports and charts (`fleet.py`, `status.py`, `progress.py`, `package.py`, `make_docx.py`, `report.py`, `summary.py`, `charts.py`) and the SMILES builder of other studies (`build_from_smiles.py`) | not needed to compute or check any result; the pipeline skips the `report` stage when `report.py` is absent |
| `check_captions.py` and the figkit demos and change log | the first checks captions in the manuscript, which is not deposited; the demos use the data of another study |
| The full literature list of the project | internal working file; the five records the scripts parse are deposited verbatim as `analysis/literature_claims.jsonl` |
| Manuscript sources, internal review notes, plans and audits | internal working files |
| Bioassay data | not part of this deposit (see README) |

## Workspace paths and this layout

The scripts were written for the project workspace and quote its paths in their outputs (the `source` of each entry
of `facts.json`, the inputs in `tables/INDEX.md`, code comments). `scripts/paths.py` maps them; in short:

| Workspace path | Here |
|---|---|
| `molecules/btpc/cases/<case>.json` | `cases/<case>.json` |
| `molecules/btpc/input/BTPC.cif` | `structure/BTPC.cif` (not deposited; or `$BTPC_CIF`) |
| `molecules/btpc/input/btpc_cluster.xyz` | `structure/btpc_cluster.xyz` (rebuilt, not deposited) |
| `molecules/btpc/results/<case>/<file>` (tracked copy) and `data/r2/molecules/btpc/results/<case>/<file>` (raw archive) | `calculations/<case>/<group>/<file>`, group by file type: `geometry/`, `wavefunction/`, `frequencies/`, `orbitals/`, `multiwfn/`, `grids/` (as `<file>.xz`), `run/` |
| `…/results/<case>/logs/<log>` | `calculations/<case>/run/logs/<log>` |
| `…/results/<case>/mwfn_home/settings.ini` | `calculations/<case>/multiwfn/settings.ini` |
| `data/r2/logs/btpc_gas.out`, `btpc_mep_retry.out` | `calculations/btpc_gas/run/` |
| `data/r2/logs/btpc_pcm.out`, `btpc_pcm_freq_step0.01.out` | `calculations/btpc_pcm/run/` |
| `papers/btpc/data/{facts.json,facts.md,derived/}` | `analysis/` |
| `papers/btpc/data/*.py` | `scripts/analysis/` |
| `papers/btpc/lit/candidates.jsonl` | `analysis/literature_claims.jsonl` (the five records used) |
| `papers/btpc/tables/` | `tables/` (script: `scripts/tables/make_tables.py`) |
| `papers/btpc/figures/` and `…/figures/scripts/` | `figures/` and `scripts/figures/` |
| `papers/_shared/figkit/` | `scripts/figures/figkit/` |
| `papers/_shared/tools/{include,crossref}.lua` | `scripts/tables/pandoc/` |
| `dft/`, `run_case.sh`, `setup_linux.sh` | `scripts/dft/`, `scripts/run_case.sh`, `scripts/setup_linux.sh` |

`paths.legacy_root()` builds `.layout/`, a git-ignored tree of relative symbolic links that shows the deposited files
under their workspace paths; the analysis, table and figure scripts run against it, so their code is unchanged apart
from the lines listed below. The DFT pipeline runs from the repository root instead and writes a recomputation to
`results/<case>/`, never into `calculations/`.

## Changes made to the deposited copies

Results were not edited. The changes below remove absolute paths, host names and user names of the compute
environment and citations of internal files, and adapt the scripts to this layout. Every rule was applied by a script
that checked the expected number of replacements.

### Data files

| Rule | Files | What changed | Count |
|---|---|---|---:|
| L1 | calculations/btpc_gas/run/btpc_mep_retry.out, calculations/btpc_gas/run/logs/optimize.log, calculations/btpc_gas/run/logs/plots.log, calculations/btpc_gas/run/logs/render.log, calculations/btpc_gas/run/logs/report.log, calculations/btpc_pcm/run/logs/frequencies.log, calculations/btpc_pcm/run/logs/optimize.log, calculations/btpc_pcm/run/logs/plots.log, calculations/btpc_pcm/run/logs/render.log, calculations/btpc_pcm/run/logs/report.log | Absolute path of the run directory on the compute VM (home folder of the VM user, project folder) replaced by `<run-dir>/` | 23 |
| L2 | calculations/btpc_gas/run/logs/optimize.log, calculations/btpc_pcm/run/logs/optimize.log | Absolute path of the local scratch volume of the compute VM replaced by `<scratch-dir>/` | 2 |
| L3 | calculations/btpc_gas/wavefunction/final.chk, calculations/btpc_pcm/wavefunction/final.chk | PySCF checkpoint: the `atom` field of the `mol` record (the absolute path of the opt.xyz the molecule was read from) replaced by the workspace-relative path `molecules/btpc/results/<case>/opt.xyz`, padded with JSON whitespace to the same length. Every byte offset and HDF5 object length is unchanged; all other datasets and all other `mol` fields were compared and are identical | 9 |

Not changed: the comment line of `calculations/btpc_pcm/geometry/preopt.xyz` ("harvested from the grid-3 IEF-PCM
optimisation trajectory on OVH (frame 16)") is a primary record and is kept; OVH is the cloud provider of the earlier
grid-3 runs, which are not archived. `run/status.json` of the IEF-PCM case is kept as written although its
`optimize` stage still says "running" (README).

### Scripts and the facts register

Rules P change paths only (behaviour unchanged); rules S, F and M remove internal references. Rules F change strings
that `scripts/analysis/build_facts.py` writes into `analysis/facts.json`; M those that `build_facts_md.py` writes
into `facts.md`.

| Rule | Files | What changed | Count |
|---|---|---|---:|
| P1 | scripts/dft/common.py, scripts/dft/render.py | Pipeline root: the directory above the script folder replaced by the repository root (scripts/paths.py `REPO`), so case files are found in `cases/` and a recomputation writes `results/<case>/` | 2 |
| P2 | scripts/dft/finalize_opt.py, scripts/dft/optimize.py, scripts/dft/pipeline.py | Starting geometry: `ROOT / case["xyz"]` replaced by `paths.resolve(case["xyz"])` (structure/btpc_cluster.xyz) | 4 |
| P3 | scripts/dft/frequencies.py, scripts/dft/pipeline.py | Script locations `ROOT/dft/...` replaced by `paths.SCRIPTS/dft/...` | 2 |
| P4 | scripts/dft/finalize_opt.py | The scratch mount searched for geomeTRIC trajectories made configurable (`SCRATCH_ROOT`, same default); `import os` added | 2 |
| P5 | scripts/dft/verify_grid.py | Grid check reads the archived geometry and energy from calculations/ (`paths.calc_file`) | 3 |
| P6 | scripts/dft/build_cluster.py | Cluster builder writes structure/<id>_cluster.xyz and structure/<id>_gas.json; the case file keeps the workspace-relative paths it always recorded; usage lines | 7 |
| P7 | scripts/dft/curate_cases.py | Case curation reads and writes structure/<id>_gas.json; usage line | 3 |
| P8 | scripts/analysis/build_facts.py, scripts/analysis/build_facts_md.py, scripts/analysis/wfn_topology.py | Analysis scripts run against the workspace-layout view (`paths.legacy_root()`); the figure-script constants are read through `paths.resolve`; usage lines | 6 |
| P9 | scripts/tables/make_tables.py | Table script runs against the view; pandoc found through `PANDOC` or the PATH; `use_file` keeps link paths (no `resolve()`); usage lines and the rebuild line of INDEX.md (which also held an absolute workspace path) | 8 |
| P10 | scripts/figures/figkit/style.py | figkit style: repository root replaced by the view | 1 |
| P11 | scripts/figures/*.py | Usage lines of the figure scripts (each held an absolute workspace path) rewritten as `python scripts/figures/<name>.py` | 23 |
| P12 | scripts/run_case.sh | run_case.sh: runs from the repository root, calls scripts/dft/pipeline.py, scratch mount configurable (`SCRATCH_ROOT`, same default); the comment naming the compute host rewritten | 4 |
| P13 | scripts/setup_linux.sh | setup_linux.sh: runs from the repository root; Python packages pinned to the recorded versions; micromamba version selectable (`MICROMAMBA_VERSION`) | 3 |
| P14 | scripts/dft/*.py (10 files) | Usage lines of the pipeline scripts: `python dft/...` written as `python scripts/dft/...` | 10 |
| S3 | scripts/analysis/build_facts.py, scripts/figures/figkit/render3d.py, scripts/figures/figkit/style.py, scripts/tables/make_tables.py | Comments naming the storage service of the project archive written as "raw results archive" | 5 |
| S4 | scripts/figures/figkit/geometry.py | Docstring wording "other modules and agents" written as "other modules and scripts" | 1 |
| S5 | scripts/dft/curate_cases.py | `META` of the case curation reduced to this compound (the entries of two other studies removed) | 2 |
| F1 | scripts/analysis/build_facts.py | Value of `protocol.host`: the host name of the compute VM replaced by "IBM Cloud mx3d-48x480 VM" (hardware description unchanged) | 1 |
| F2 | scripts/analysis/build_facts.py | Source of `protocol.host`: run paths written as `<run-dir>` and `<scratch-dir>`, the host name removed, the install log marked as not deposited | 1 |
| F3 | scripts/analysis/build_facts.py | Note of `gas.opt_wall`: host name replaced by "the 48-vCPU compute VM" | 1 |
| F4 | scripts/analysis/build_facts.py | Value of `inconsistency.pcm_hessian_step_text`: citation of an internal notes file removed; the evidence is the dft/frequencies.py comment | 1 |
| F5 | scripts/analysis/build_facts.py | Value of `inconsistency.lebedev_text`: citations of internal notes files removed; the evidence is the dft/common.py default | 1 |
| F6 | scripts/analysis/build_facts.py | Value of `inconsistency.gas_freq_policy`: citation of internal notes files written as "The project notes (not deposited)" | 1 |
| F7 | scripts/analysis/build_facts.py | Note of `protocol.pcm_lebedev_order`: citation of an internal progress file removed (the dft/common.py default remains) | 1 |
| F8 | scripts/analysis/build_facts.py | Source of `protocol.preopt`: citation of an internal progress file removed | 1 |
| F9 | scripts/analysis/build_facts.py | Note of `protocol.cubes`: citation of an internal progress file removed | 1 |
| F10 | scripts/analysis/build_facts.py | Source of `model.preopt_basis_functions`: internal progress file replaced by the basis-set count it rests on (25 heavy atoms × 14 + 8 H × 2 = 366) | 1 |
| F11 | scripts/analysis/build_facts.py | Sources of `software.python`, `software.numpy_scipy_matplotlib`, `software.libxc`: the workspace package list replaced by environment/compute-env.yml | 3 |
| F12 | scripts/analysis/build_facts.py | Source and note of `software.paraview`: install log marked as not deposited; workspace tool list replaced by environment/software.md | 2 |
| F13 | scripts/analysis/build_facts.py | Sources of `software.pyscf`, `software.geometric`, `software.multiwfn`: install log marked as not deposited | 3 |
| F14 | scripts/analysis/build_facts.py | Sources citing the grid-verification log: marked as not deposited | 3 |
| F15 | scripts/analysis/build_facts.py | Source of `inconsistency.pcm_convergence_text`: hand-over package README marked as not deposited | 1 |
| F16 | scripts/analysis/build_facts.py | Source of `inconsistency.chart_title`: pipeline chart summary marked as not deposited | 1 |
| M1 | scripts/analysis/build_facts_md.py | facts.md: citation of internal notes files about Lebedev order 17 replaced by the dft/common.py default | 1 |
| M2 | scripts/analysis/build_facts_md.py | facts.md: the storage service of the project archive written as "raw results archive" / "archive copies" | 2 |
| M3 | scripts/analysis/build_facts_md.py | facts.md: grid-verification log marked as not deposited | 1 |
| M4 | scripts/analysis/build_facts_md.py | facts.md: rebuild commands written for this repository | 1 |
| M5 | analysis/facts.md (fact-check log) | facts.md fact-check log (hand-written, kept across rebuilds): citations of an internal progress file written as "the project notes (not deposited)"; grid-verification log marked as not deposited | 3 |

Effect on `analysis/facts.json`: all 1007 keys are present and 1003 values are identical to the workspace register.
The four text values changed by F1, F4, F5 and F6 are `protocol.host`, `inconsistency.pcm_hessian_step_text`,
`inconsistency.lebedev_text` and `inconsistency.gas_freq_policy`; no number changed. Before this deposit the workspace
script was also changed in one place that affects the register: names in the `source` of computed entries are now
sorted (longest first, then alphabetically) instead of following Python's set order, so the register is the same on
every run.

The scripts keep `/mnt/scratch` only as the default value of the scratch mount (`SCRATCH_ROOT`), a generic path.
Comments that mention internal review rounds (for example "QA r4-m3", "audit/response_r1.md", "PLAN.md",
"OUTLINE.md") refer to notes that are not deposited and were left as they are.

## Checksums of the uncompressed grids

`scripts/unpack_grids.sh` checks the unpacked grids against `docs/grids.sha256`:

| File | SHA-256 |
|---|---|
| `calculations/btpc_gas/grids/density.cube` | `ad3c31819032f4d65e98849aff91cd6e6044c0de76e3af2fe625f133fae6a547` |
| `calculations/btpc_gas/grids/elf_benzotriazolium_ring_plane.txt` | `74364530c2c745b54e3f85b6e9fb69ee6a91f84de38a63e33495a651d342bf58` |
| `calculations/btpc_gas/grids/elf_picrate_ring_plane.txt` | `549b708245added588b614fb5d9e33234a0c86cb629ca10e9660fa2c2c3032f0` |
| `calculations/btpc_gas/grids/esp.cube` | `a7f60eeef6bf82b9551f2d0822272cd896211a023907100518c03e6a982063ce` |
| `calculations/btpc_gas/grids/homo.cube` | `12c8a6663714f9209c15723fed2f9789e46ce3b38ce92c477a56e1e400e03818` |
| `calculations/btpc_gas/grids/lol_benzotriazolium_ring_plane.txt` | `64f0daa83dec6544e7aa4291956cc53e18002eea961f10a256f91fa9830607ee` |
| `calculations/btpc_gas/grids/lol_picrate_ring_plane.txt` | `0302a1eb3d77622aec56876bb3724d4b7566135812c45ace6815c6b5a518ca8b` |
| `calculations/btpc_gas/grids/lumo.cube` | `48aa8c1d406140cda322af5a5d5a32df830d5b242896823ab0f3cb01eea90913` |
| `calculations/btpc_gas/grids/rdg.cube` | `edc3eddd44c0fbadc6224f1b88b0a52c59ed7487fd46c2973202d50ff085c828` |
| `calculations/btpc_gas/grids/rdg_scatter.txt` | `63a7f65170b28c29c24657541798c7ed60cc91f06ab9b15e19a546018913813f` |
| `calculations/btpc_gas/grids/sl2r.cube` | `43d5eaa012bd88cd7b348262d9e4bde1307a56a280c2a716970dd0280f0b4928` |
| `calculations/btpc_pcm/grids/density.cube` | `211a34606c9f0a27709f91ec5a04b1235669b88bc9743161f859e409ede00915` |
| `calculations/btpc_pcm/grids/elf_benzotriazolium_ring_plane.txt` | `64e2bbf22a89eeb7b091edb711529eafb21dab91d0fe5a40a38c9f993ab039b4` |
| `calculations/btpc_pcm/grids/elf_picrate_ring_plane.txt` | `bb4e505704ce4e7150c2e10ba58427b54cc9a085f4354561c0886836460ecc16` |
| `calculations/btpc_pcm/grids/esp.cube` | `4da660e76dced3597c779472c31dfcde37ee325acc3807c4a1a2d6a707798319` |
| `calculations/btpc_pcm/grids/homo.cube` | `8a7cb76219f55ba9bec257c04ffb3ec19e9858e40edf200b71de1ab47d53f57e` |
| `calculations/btpc_pcm/grids/lol_benzotriazolium_ring_plane.txt` | `8e70453bd9457399cd779dc36a88456f44f64df13f21cb007bb3b3fc07c190a9` |
| `calculations/btpc_pcm/grids/lol_picrate_ring_plane.txt` | `76f3e796ce880239282eb8921e4e1da3c39bfc9829a90f238c0f47ce31463ffa` |
| `calculations/btpc_pcm/grids/lumo.cube` | `9020d3a61504c4491f8b1bfad00db87c72b7d6ee119a6cfc6f3d5910201c6971` |
| `calculations/btpc_pcm/grids/rdg.cube` | `8151552393fbf00babc6c41458e80cf49bde1e00437c58a7f378045557afe553` |
| `calculations/btpc_pcm/grids/rdg_scatter.txt` | `61870a4ede9565cf8f7cf2794c6c3afd13715a7a83cf4b8a1890a544dbe9d202` |
| `calculations/btpc_pcm/grids/sl2r.cube` | `375e980097c0cdb82a07a15475de827193bcdf88c275b9c58aa9305d83cbb473` |

## Revision calculations (added 23 Nov 2026)

The calculations of the revision of the article (`docs/revision.md`) were copied on 23 Nov 2026 from the project's
archive of the two compute VMs: the queue folders (job outputs, per-job logs, runner records), the pipeline results
of the second crystal form, the environment records of each VM, and the project folder of the revision scripts.
They were copied by an explicit allowlist; results were not edited.

| Project path | Here |
|---|---|
| `<run-dir>/queue/out/<job>/` of either VM | `calculations/revision1/<job>/` |
| `<run-dir>/queue/logs/<job>.log` | `calculations/revision1/<job>/queue.log` |
| `<run-dir>/queue/status.json`, `runner.log`, `runner.out` | `calculations/revision1/runner/vm1/`, `vm2/` |
| `<run-dir>/env.lock.txt`, `env.versions.txt`, the installer log | `environment/revision1/vm1/`, `vm2/` (`setup.log`) |
| `molecules/btpcz/results/<case>/<file>` | `calculations/<case>/<group>/<file>`, grouped as above; `softmode_check.json` in `frequencies/` |
| `molecules/btpcz/cases/<case>.json` | `cases/<case>.json` |
| `molecules/btpcz/input/btpcz_cluster.xyz` | `structure/btpcz/btpcz_cluster.xyz` |
| `molecules/btpcz/input/2230422.cif`, `README.md` | not deposited; `structure/btpcz/README.md` gives the source, checksum and specification |
| `dft/revision/` | `scripts/revision1/` (edits in `scripts/revision1/CHANGES.md`) |

Not deposited: the `chk/` folders of the revision jobs (222 PySCF checkpoints and 20 geomeTRIC constraint files,
746 MiB), the per-job exit-code files (the codes are in `status.json`), empty lock files, a copy
of `sp_neutral_pcm/geoms/neutral_pcm.xyz` placed on the second VM (byte-identical), and for the second crystal form
`opt.chk` (2 files, 31 MiB) and the pipeline pictures (15 files, 7.7 MiB), as for the published cases. The grids of
the second form were xz-compressed like the others; their SHA-256 before compression are in `docs/grids.sha256`.

Changes to the deposited copies (text files; counted replacements):

| Rule | Files | What changed | Count |
|---|---|---|---:|
| D1 | 239 files in `calculations/revision1/`, 2 in `environment/revision1/` | Host name of the first compute VM (its domain names the cloud project) replaced by `<vm-1>` | 445 |
| D2 | 9 files in `calculations/revision1/`, 2 in `environment/revision1/` | Host name of the second compute VM replaced by `<vm-2>` | 13 |
| D3 | 344 files in `calculations/revision1/`, 9 in `calculations/btpcz_*/`, 4 in `environment/revision1/` | Absolute path of the working folder on the compute VMs (home folder of the VM user) replaced by `<run-dir>` | 4170 |
| D4 | `calculations/revision1/c0_prep/queue.log` | User and group of the VM user in `ls -l` output replaced by `<user>` | 2 |
| L3 | `calculations/btpcz_gas/wavefunction/final.chk`, `calculations/btpcz_pcm/wavefunction/final.chk` | As rule L3 above: the absolute path of the `opt.xyz` in the `mol` record replaced by `molecules/btpcz/results/<case>/opt.xyz`, padded with JSON whitespace to the same length; offsets and object lengths unchanged | 4 |

The command `$HOME/r007/...` recorded for the host-preparation job `c0_prep` is kept as written (`r007` is the name
of the working folder, not a path of any machine).
