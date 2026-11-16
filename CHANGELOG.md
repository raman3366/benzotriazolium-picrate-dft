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

## 1.0.0 (planned, about 1 Dec 2026)

- Public release with the Zenodo DOI; adds the final figures and the DFT supplementary information.
