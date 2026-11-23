# Reproducing the data

Two levels. Level A needs no quantum chemistry and takes a few minutes; it rebuilds every number of the article and
every table from the deposited files and checks that nothing changes. Level B recomputes the calculations.

Run everything from the repository root.

## Level A: check, rebuild the numbers and the tables (minutes)

Requirements: Linux or macOS, `xz`, `sha256sum`, Python 3.12 with the packages of
`environment/paper-requirements.txt` (no PySCF needed), and the CIF of CCDC 835950. pandoc 3.11 or newer is
optional (for the pandoc part of `make_tables.py --check`). Older versions, such as 3.1.3 (the Ubuntu 24.04 package),
do not write the table ids the check looks for, so `reproduce.sh` runs `--check` only with pandoc 3.11 or newer;
with an older pandoc, or none, it builds the tables without the pandoc checks and prints a note.

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r environment/paper-requirements.txt
# Obtain CCDC 835950 free of charge from https://www.ccdc.cam.ac.uk/structures/ (doi:10.5517/ccx1w4y)
# and save it as structure/BTPC.cif (or set BTPC_CIF=/path/to/BTPC.cif).
PYTHON=.venv/bin/python scripts/reproduce.sh
```

**Your CIF will not have the SHA-256 of the copy used here.** The CCDC writes the download date into every CIF it
serves (`_audit_update_record`; the copy used here reads "2011-07-21 deposited with the CCDC." and "2026-06-25
downloaded from the CCDC." on one line), and that copy has CRLF line ends. `reproduce.sh` therefore only reports how
the CIF compares: identical to the copy used (SHA-256 `17065c7a…5349b5`), the same apart from the download line and
the line ends (SHA-256 `9d93fe1c…4c67d9` after removing the line with "downloaded from the CCDC" and the carriage
returns), or different (a warning). It does not stop there. The checks that count come next: the rebuilt ion pair
`structure/btpc_cluster.xyz` must have SHA-256 `dfa71f6f…f95476`, and the rebuilt case file must be identical to
`cases/btpc_gas.json`.

`scripts/reproduce.sh` runs these steps; each can also be run by hand:

```sh
sha256sum -c SHA256SUMS                                   # every line "OK"
scripts/unpack_grids.sh                                   # xz -dk calculations/*/grids/*.xz, checked against docs/grids.sha256
tr -d '\r' < structure/BTPC.cif | grep -v 'downloaded from the CCDC' | sha256sum
#   a guide only (see above): 9d93fe1c427457780670219100894a7ff2b48af9050038d552bc9e82104c67d9
python scripts/dft/build_cluster.py structure/BTPC.cif btpc --title Benzotriazolium-picrate
echo "dfa71f6f1a4b15cf9b9f8f4cf4240f88ef284b073fef5971cfe9c3e5a0f95476  structure/btpc_cluster.xyz" | sha256sum -c
python scripts/dft/curate_cases.py                        # then: cmp structure/btpc_gas.json cases/btpc_gas.json
ref=$(mktemp) && cp analysis/facts.json "$ref"
python scripts/analysis/build_facts.py                    # rewrites analysis/facts.json and analysis/derived/
python scripts/tools/compare_facts.py "$ref" analysis/facts.json
python scripts/analysis/build_facts_md.py                 # rewrites analysis/facts.md
python scripts/tables/make_tables.py --check              # rewrites tables/*.md; --check needs pandoc >= 3.11
                                                          # (PANDOC=/path/to/pandoc); with an older one, omit --check
python scripts/tables/tables_to_csv.py                    # rewrites tables/*.csv
sha256sum -c SHA256SUMS                                   # the rebuilt files must still be "OK"
```

The analysis and table scripts read the deposited files through `.layout/`, a git-ignored tree of symbolic links in
the layout of the project workspace (`scripts/paths.py`; docs/provenance.md). The rebuilt ion pair, the draft case
file, the unpacked grids and `.layout/` are git-ignored.

### Expected outputs

Recorded on 16 Nov 2026 (IST) on a 4-core ARM64 machine (Raspberry Pi 5) with Python 3.12.3 and the packages of
`environment/paper-requirements.txt`, pandoc 3.11, and the CIF the calculations used (downloaded on 25 Jun 2026,
SHA-256 `17065c7a…5349b5`):

1. **Ion pair from the CIF** (`build_cluster.py`):

   ```
   [cluster] assembled 33 atoms, seed fragment 0, H-bond score 1.37, all fragments hydrogen-bonded: True
             fragment 0 (18 atoms) via symmetry image 1_555
             fragment 1 (15 atoms) via symmetry image 4_666
               N4#4_666-H4N#4_666...O1         H...A 1.775 A  D...A 2.623 A  angle 154.5 deg
   dfa71f6f1a4b15cf9b9f8f4cf4240f88ef284b073fef5971cfe9c3e5a0f95476  structure/btpc_cluster.xyz
   ```

   `curate_cases.py` then writes `structure/btpc_gas.json` byte-identical to `cases/btpc_gas.json`. Copies of the
   CIF with another download date, with and without the carriage returns, gave the same ion pair, case file, facts
   register, tables and checksums (16 Nov 2026); `reproduce.sh` reports such a copy as "same content as the copy the
   calculations used (download date or line ends differ)".
2. **Facts register** (`build_facts.py`): "wrote 1007 facts" after about 28 s (16 Nov 2026; since 23 Nov 2026,
   with the revision results: "revision additions: 1483", "wrote 2490 facts", about 30 s, and the fourteen
   `analysis/derived/rev_*.csv` files, byte-identical to the deposited ones; `compare_facts.py` reports "2490 / 2490
   keys; 0 values differ"). `analysis/facts.json` and the three CSV
   files in `analysis/derived/` come out byte-identical to the deposited ones (also with a different
   `PYTHONHASHSEED`); `compare_facts.py` reports "1007 / 1007 keys; 0 values differ". `build_facts_md.py` rewrites
   `analysis/facts.md` unchanged (940 lines).
3. **Tables** (`make_tables.py --check`): 17 tables (main Tables 1–7, SI Tables S1–S10) and `INDEX.md`, all
   byte-identical to the deposited files; the pandoc checks pass (HTML for every table, DOCX for `tbl_geometry`,
   cross-references numbered 1–7 and S1–S10). The check of the main-table order against the built manuscript is
   skipped, because the manuscript build is not deposited. With pandoc 3.1.3 on the PATH, `reproduce.sh` builds
   the same tables without the pandoc checks and says so.
4. **Checksums**: `sha256sum -c SHA256SUMS` reports OK for all 282 files, before and after the rebuild (1048 files
   once the revision calculations were added, 1150 with the paper outputs of the revised article, both on
   23 Nov 2026; the gate gave the same results then).

The register in the project workspace held the same 1007 values; four text values differ only by the substitutions
listed in docs/provenance.md (rules F1, F4, F5, F6). With the revision (23 Nov 2026): the same 2490 values, the same
four text values changed, and the notes of the ten `rev.exp.*` keys (rule F17); the `rev_*.csv` files and all 28
tables are byte-identical to the project's. The eleven revision tables are rebuilt by `scripts/tables/rev_*_tables.py`
and `fix_widths.py` (`reproduce.sh` step 5, "11 revision tables rebuilt").

### Figures

Each figure of the article is made by one script in `scripts/figures/` (for example
`python scripts/figures/fig_mep.py`) and written to `figures/`. The scripts check their own layout strictly
(`figkit.qa`) and need the Liberation Sans font (package `fonts-liberation`) where the toolkit looks for it:

```sh
mkdir -p tools/fonts/extracted/usr/share/fonts/truetype
ln -s /usr/share/fonts/truetype/liberation tools/fonts/extracted/usr/share/fonts/truetype/liberation
```

Without it they fall back to DejaVu Sans, and some figures stop with "text outside the canvas" or "text overlap".
Checked on 16 Nov 2026 with the font in place: `fig_descriptors.py` (PNG pixel-identical to the figure of the
article), `fig_scheme.py`, `fig_hbond_trajectory.py`, `figS_trajectory_full.py` and `figS_mo_composition.py` run
without error. The 3-D figures need PyVista with off-screen rendering. On 23 Nov 2026 `fig_rev_eps.py`,
`fig_rev_proton.py`, `fig_rev_torsion.py` and `graphical_abstract.py` were run inside this repository: the PNG files
are byte-identical to the deposited figures of the article, and the PDF files differ only in their creation date. The
graphical abstract script also writes a TIFF copy (git-ignored).

## Revision calculations

`sha256sum -c SHA256SUMS` covers the files of the revision calculations (`docs/revision.md`) as well, and
`scripts/unpack_grids.sh` unpacks and checks the grids of the second crystal form with the others (44 in all).
`structure/btpcz/README.md` shows how to rebuild the second form's starting geometry and case files from the COD
CIF; `docs/revision.md` shows how to run one revision job again. The facts register, the derived CSV files and the
tables above include the revision results (`rev.*` keys).

## Level B: full recomputation (Linux x86-64)

Recorded cost on the 48-vCPU VM: gas optimisation 14 min on 4 threads; IEF-PCM optimisation about 13.4 h on
4 threads; gas Hessian 149 min (16 workers × 2 threads); IEF-PCM Hessian 346 min (10 workers × 2 threads). In all,
about 250–300 vCPU-hours (an estimate from these wall times). Memory: `DFT_MAX_MEMORY_MB=40000` (gas) and `12000`
(IEF-PCM) were used; PySCF writes about 12 GB of scratch per gradient, so put scratch on a large local disk
(`SCRATCH_ROOT`, default `/mnt/scratch` when it exists, else `scratch/`).

```sh
./scripts/setup_linux.sh          # pinned conda environment, Multiwfn 2026.9.13, ParaView 6.1.1 -> tools/
# the final runs started from the archived preopt.xyz (the earlier grid-3 runs are not archived); the pipeline
# reuses a preopt.xyz found in the results folder and then skips the 6-31G(d) pre-optimisation
for c in btpc_gas btpc_pcm; do mkdir -p results/$c && cp calculations/$c/geometry/preopt.xyz results/$c/; done
NO_AUTO_PCM=1 DFT_MAX_MEMORY_MB=40000 ./scripts/run_case.sh btpc_gas 4
NO_AUTO_PCM=1 DFT_MAX_MEMORY_MB=12000 ./scripts/run_case.sh btpc_pcm 4
# IEF-PCM: if geomeTRIC stalls on negligible forces, accept the last step and continue:
#   tools/env/bin/python scripts/dft/finalize_opt.py btpc_pcm && ./scripts/run_case.sh btpc_pcm 4 -- --from proton_check
# the reported IEF-PCM Hessian used a 0.01 Bohr step:
NO_AUTO_PCM=1 DFT_MAX_MEMORY_MB=12000 FD_STEP=0.01 ./scripts/run_case.sh btpc_pcm 10 -- --from frequencies
```

- `run_case.sh` sets `GRID_LEVEL=4` and `PCM_LEBEDEV_ORDER=29`. Calling `scripts/dft/pipeline.py` directly falls back
  to the code defaults (grid level 3, Lebedev order 17), a different level of theory.
- `NO_AUTO_PCM=1` stops the gas-phase run from writing a new IEF-PCM case file over `cases/btpc_pcm.json`.
- Results go to `results/<case>/` in the flat layout of the workspace; compare them with `calculations/<case>/`.
- The pipeline's `report` stage is skipped (its script is not deposited).

Expected values (from `opt_summary.json`, `fmo.json`, `freq.json`):

| Quantity | Gas phase | IEF-PCM (water) |
|---|---|---|
| E | −1317.20823094 Ha (3 steps) | −1317.25368856 Ha (70 steps, accepted on negligible forces) |
| HOMO / LUMO / gap | −6.818 / −3.342 / 3.476 eV | −6.903 / −3.167 / 3.735 eV |
| Imaginary modes | none; lowest 12.19 cm⁻¹ | one, at −3.22 cm⁻¹ |
| N4···O1 | 2.478 Å | 2.672 Å |

Tolerance: SCF energies to about 10⁻⁶ Ha, and wavenumbers below 50 cm⁻¹ to about 2 cm⁻¹ (changing the
finite-difference step from 0.005 to 0.01 Bohr moved no mode by more than 1.89 cm⁻¹). The IEF-PCM optimisation ends
on a very flat surface, so a rerun may stop at a slightly different point; the geometry and energy should agree
within these tolerances.

`scripts/dft/verify_grid.py <case>` recomputes the SCF at the archived geometry with grid levels 3 and 4 (and
Lebedev orders 17 and 29 for IEF-PCM) and reports which setting reproduces the archived energy.
