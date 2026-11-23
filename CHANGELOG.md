# Changelog

## 0.1.0 (16 Nov 2026), private build, not released

- First build from the project archive: case files, both calculations (geometry, wavefunction, frequencies, orbitals,
  Multiwfn, grids, run records), the facts register with its scripts, the tables with CSV copies, the DFT pipeline,
  the analysis, table and figure scripts, pinned environments and documentation.
- Checked: the ion pair rebuilds byte for byte from the CIF; `scripts/analysis/build_facts.py` rebuilds
  `analysis/facts.json` byte for byte; `scripts/tables/make_tables.py --check` rebuilds every table unchanged and
  passes its pandoc checks; `SHA256SUMS` verifies (docs/REPRODUCE.md).
- Not yet included: the rendered figures (`figures/`) and the DFT supplementary information (`si/`), which come with
  the revised article at v1.0.0.

## 1.0.0 (23 Nov 2026)

- Added on 23 Nov 2026 (private): the calculations of the revision of the article (`docs/revision.md`).
  - `calculations/revision1/`: the 57 revision jobs (proton transfer along N4–H4N···O1 with scans, other
    functionals, DF-MP2, a CCSD(T) correction and proton vibrational levels; the location of the LUMO against the
    permittivity, ΔSCF and Hirshfeld/CM5 charges; B3LYP-D3(BJ) minima with soft-mode Hessian checks; inter-ion
    torsion scans; TD-DFT), each with its record (`result.json`), logs and output files, and the job-queue records
    of the two compute VMs. PySCF checkpoints of these jobs are not included.
  - The second crystal form (P2₁/n, COD 2230422): `calculations/btpcz_gas/`, `calculations/btpcz_pcm/` in the
    layout of the published cases, `cases/btpcz_*.json`, `structure/btpcz/` (source, checksum, specification and
    starting geometry; the CIF is not included). The starting geometry and case files rebuild byte for byte from
    the COD CIF.
  - `scripts/revision1/` (scripts, job lists, tests, environment installer, VM script, local environment locks;
    edits for deposit in `CHANGES.md`) and `environment/revision1/` (explicit conda locks of the two VMs).
  - `docs/revision.md`; README section "Revision calculations"; `docs/provenance.md` (mapping and rules D1–D4, L3);
    `docs/grids.sha256` and `scripts/unpack_grids.sh` cover the 22 new grids (44 in all).
  - `scripts/sync_paper_outputs.sh` copies the final facts register, tables, figures, DFT supplementary
    information and analysis scripts from the project workspace; to be run once when they are final (not run yet).
- Added on 23 Nov 2026 (private): the paper outputs of the revised article.
  - `analysis/facts.json` with the 1483 `rev.*` keys of the revision calculations (2490 in all), `facts.md`, the
    fourteen `analysis/derived/rev_*.csv` files and `analysis/experimental.json`, rebuilt inside this repository by
    `scripts/analysis/build_facts.py` with the revision fact modules (`scripts/analysis/rev/`); every value equals
    the project register apart from the text values of rules F1, F4–F6.
  - The eleven revision tables (`scripts/tables/rev_*_tables.py`, `fix_widths.py`) and their CSV copies; all 28
    tables byte-identical to the project's.
  - The rendered figures (PNG and PDF, 23 figures) and the graphical abstract, with the scripts of the revision
    figures; the DFT supplementary information (`si/`).
  - `scripts/paths.py` maps the revision paths; `docs/scrubbed_files.tsv` gives the as-run checksums of the files
    changed for deposit; `scripts/reproduce.sh` rebuilds the revision tables too; docs updated (`docs/provenance.md`
    rules P15–P22, F17, F18, S10).
- Still to come before the release: the DOI.
