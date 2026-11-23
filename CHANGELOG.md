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

## 1.0.0 (unreleased)

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
- Still to come before the release: the revised facts register, tables, figures and DFT supplementary information,
  and the DOI.
