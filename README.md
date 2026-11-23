# Benzotriazolium picrate (BTPC): DFT data on the N4–H4N···O1 ion pair in the gas phase and in IEF-PCM water

Data and scripts for the density functional study of the benzotriazolium picrate ion pair (BTPC;
C₆H₆N₃⁺·C₆H₂N₃O₇⁻; CCDC 835950). The ion pair held together by the charge-assisted N4–H4N···O1 hydrogen bond was
rebuilt by symmetry from the crystal structure and optimised at B3LYP/6-311++G(d,p) twice: in the gas phase and in
an IEF-PCM continuum for water. The deposit holds both optimisations, both finite-difference Hessians, the
wavefunctions, the frontier-orbital, electrostatic-potential, NCI, ELF and LOL analyses, the run logs, and the
register of every number in the associated article (`analysis/facts.json`) with the scripts that rebuild it, the
tables, the figures and the DFT supplementary information of the revised article. It also holds the calculations added during the revision of the article: a second crystal form of the salt
(P2₁/n, COD 2230422) treated in the same way, and 57 further jobs on proton transfer, the location of the LUMO,
dispersion, the inter-ion torsion and excited states ("Revision calculations" below; `docs/revision.md`).

> **Status: private build, not a release.** The facts register, tables, figures and DFT supplementary information
> are those of the revised article (23 Nov 2026). The DOI is added at release v1.0.0; nothing here is final before
> then.

## Contents

| Folder | Contents |
|---|---|
| `structure/` | Citation and checksum of the crystal structure (not redistributed), the ion-pair specification, and the checksum of the rebuilt starting geometry; `structure/btpcz/`: the same for the second form, with its starting geometry |
| `cases/` | One JSON file per calculation: atom labels and order, fragments, the watched hydrogen bond, ELF/LOL planes, solvent (`btpc_*`: the published form; `btpcz_*`: the second form) |
| `calculations/btpc_gas/`, `calculations/btpc_pcm/` | Per model: `geometry/` (start, optimised, trajectory), `wavefunction/` (Molden, PySCF checkpoint), `frequencies/` (Hessian, harmonic analysis), `orbitals/`, `multiwfn/` (inputs, settings and outputs), `grids/` (cubes and grids, xz-compressed), `run/` (logs and stage records) |
| `calculations/btpcz_gas/`, `calculations/btpcz_pcm/` | The second crystal form (revision), in the same layout |
| `calculations/revision1/` | The 57 revision jobs: one folder per job with its record (`result.json`), logs and output files, and the job-queue records (`README.md` there lists every job) |
| `analysis/` | `facts.json` (2490 numbers, each with unit, source file and note: 1007 of the published calculations and 1483 `rev.*` keys of the revision calculations), `facts.md` (a readable dossier of the published part), `derived/` (MO fragment shares, crystal hydrogen bonds, and the `rev_*.csv` tables of the revision results), `experimental.json` (the co-authors' experimental values quoted in the article), `literature_claims.jsonl` (the five literature values the scripts parse) |
| `tables/` | Every table of the article and its SI as Markdown (`tbl*.md`: 17 made by `scripts/tables/make_tables.py`, 11 of the revision by `scripts/tables/rev_*_tables.py`) with CSV copies; `INDEX.md` lists the facts keys and files behind the first 17 |
| `figures/` | Rendered figures of the article and SI (PNG and PDF) and the graphical abstract (`graphical_abstract/`); every figure is made by a script in `scripts/figures/` |
| `si/` | The DFT supplementary information of the revised article as PDF |
| `scripts/` | DFT pipeline (`dft/`, `run_case.sh`, `setup_linux.sh`), analysis, tables, figures, `paths.py`, `reproduce.sh`, `unpack_grids.sh`; `revision1/`: the scripts, job lists and environment installer of the revision calculations |
| `environment/` | Pinned software environments and the Multiwfn and ParaView details; `revision1/`: the environment locks of the revision calculations |
| `docs/` | `protocol.md`, `data_dictionary.md` (formats, units, atom indexing), `provenance.md` (what was changed for deposit), `REPRODUCE.md`, `revision.md` (the revision calculations) |

`SHA256SUMS` lists the SHA-256 of every file; `CHANGELOG.md` lists the versions.

## The model

- Crystal structure: CCDC 835950, monoclinic P2₁/c, Z = 4, 298 K (Dhanabal et al., *Optik* 124 (2013) 1966–1970).
- The two ions of the asymmetric unit are not hydrogen-bonded to each other, so the pair is assembled from symmetry
  images: picrate at x, y, z (1_555) and benzotriazolium at 1 + x, ½ − y, ½ + z (4_666, symmetry operation 4 of the
  CIF, x, −y − ½, z − ½, translated by one cell along each axis). The pair keeps the
  N4–H4N···O1 bond to the phenolate oxygen (X-ray: H···O 1.775 Å, N···O 2.623 Å, 154.5°). The atom order is the CIF
  label order stored in `cases/*.json` (`cif_labels`): picrate atoms 0–17, cation atoms 18–32.
- Level of theory: B3LYP (PySCF `b3lyp`, VWN-RPA local correlation as in Gaussian)/6-311++G(d,p) with spherical d
  functions, density fitting, PySCF integration grid level 4.
- Models: gas phase (`btpc_gas`) and IEF-PCM, water, ε = 78.3553, Lebedev order 29 (`btpc_pcm`).
- Optimiser geomeTRIC 1.1.1; Hessians by central finite differences of analytic gradients; analyses with Multiwfn
  2026.9.13.

| Quantity | Gas phase | IEF-PCM (water) |
|---|---|---|
| Total energy | −1317.20823094 Ha | −1317.25368856 Ha |
| HOMO / LUMO / gap | −6.818 / −3.342 / 3.476 eV | −6.903 / −3.167 / 3.735 eV |
| N4···O1 (X-ray 2.623 Å) | 2.478 Å | 2.672 Å |
| Imaginary modes | none (lowest 12.19 cm⁻¹) | one, 3.22i cm⁻¹ (see below) |

`docs/protocol.md` gives every setting; `docs/data_dictionary.md` gives the formats and units (energies in Ha,
coordinates in Å, cube grids in Bohr, wavenumbers in cm⁻¹, orbital energies in eV, atom indices 0-based in the JSON
files and 1-based in the Multiwfn inputs).

## Revision calculations (November 2026)

During the revision of the article further calculations were run with PySCF 2.14.0 at the level of theory above
unless stated otherwise (`docs/revision.md` gives the details, `calculations/revision1/README.md` the job list):

- **The second crystal form**, P2₁/n (Zeng, Li and Wang, *Acta Cryst.* E67 (2011) o1464; COD 2230422), whose
  asymmetric unit is itself an N6–H6···O1 hydrogen-bonded ion pair: optimisations in the gas phase and in IEF-PCM
  water, finite-difference Hessians, frontier orbitals, Multiwfn analyses and charges (`calculations/btpcz_gas/`,
  `calculations/btpcz_pcm/`, `structure/btpcz/`).
- **The location of the LUMO** at both published minima: permittivity series (ε = 1 to 78.3553), ΔSCF N+1/N−1
  states, Hirshfeld and CM5 charges, PBE0 and CAM-B3LYP.
- **Dispersion:** B3LYP-D3(BJ) re-optimisations in both models with soft-mode Hessian checks.
- **Proton transfer** along N4–H4N···O1: relaxed scans, optimisations started from the neutral form, single points
  with ωB97X-D3(BJ), M06-2X, PBE0, CAM-B3LYP and DF-MP2/SCS-MP2, one-dimensional proton vibrational levels, and a
  CCSD(T) correction to the DF-MP2 salt–neutral energy difference with frozen natural orbitals.
- **The inter-ion torsion** C1–O1···N4–C12: relaxed scans with B3LYP and B3LYP-D3(BJ).
- **Excited states:** TD-DFT (Tamm–Dancoff) with B3LYP and CAM-B3LYP for the ion pair and the free ions.

The jobs ran through a small job queue on two 32-vCPU cloud VMs (55 and 2 jobs, 21–23 Nov 2026) in an environment
installed from the explicit conda lock in `environment/revision1/`. Every job folder holds the job's record
(`result.json`: inputs with SHA-256, settings, software versions, script SHA-256, results) and its logs. The PySCF
checkpoints of these jobs are not deposited. The scripts are in `scripts/revision1/`; `docs/revision.md` shows
how to run one job again. Their results enter the facts register as the `rev.*` keys (`scripts/analysis/rev/`, read
from the job records), and from there the revision tables and figures, like every other number of the article.

## Provenance and known limitations

- **The final optimisations did not start from the crystal geometry.** Both started from structures of earlier runs
  at integration grid level 3: the gas phase from the converged grid-3 geometry, IEF-PCM from frame 16 of the grid-3
  IEF-PCM trajectory (comment lines of `geometry/preopt.xyz`). Those earlier runs are not archived, so the archived
  chain is CIF → ion pair → [not archived] → `preopt.xyz` → `opt.xyz`.
- **The IEF-PCM optimisation was accepted on negligible forces** (RMS 1.40 × 10⁻⁴, maximum 3.85 × 10⁻⁴ Ha Bohr⁻¹)
  after 70 steps without meeting the displacement criteria (`geometry/opt_summary.json`, `converged_on`;
  `scripts/dft/finalize_opt.py`). `run/status.json` of this case still shows the `optimize` stage as "running"; that
  record is stale and is kept as written; `run/STAGES.md` has the outcome.
- **The IEF-PCM Hessian has one imaginary mode**: −3.22 cm⁻¹ with a 0.01 Bohr step (`freq.json`, `hessian.npy`) and
  −2.48 cm⁻¹ with 0.005 Bohr (`freq_step0.005.json`; only the wavenumbers of this check are archived, not its
  Hessian). The gas-phase Hessian (0.005 Bohr) has no imaginary mode.
- **Gas-phase MEP.** Multiwfn's LIBRETA ESP code crashed on the gas-phase wavefunction, so that case used the
  original ESP code (`iESPcode = 1`); IEF-PCM used LIBRETA (`iESPcode = 2`). See `run/btpc_mep_retry.out`.
- **Not recorded:** the libxc version inside the PySCF build.
- **Paths and host names were removed from the deposited copies** of the logs, scripts and `facts.json`; the
  substitutions are listed in `docs/provenance.md`. The `source` entries of `facts.json` and the inputs listed in
  `tables/INDEX.md` quote paths of the project workspace (for example `molecules/btpc/results/btpc_pcm/opt.xyz`);
  `scripts/paths.py` and `docs/provenance.md` map them to this layout.
- Comments in the scripts that mention review rounds (for example "QA r4-m3", "audit/response_r1.md", "PLAN.md")
  refer to internal checking rounds of the manuscript; those notes are not part of the deposit.
- No bioassay data are included. The biological measurements reported in the article are in its supplementary
  material and are available from the corresponding author.

## How to reproduce

- **Level A** (minutes, no quantum chemistry): check the files, rebuild the ion pair from the CIF, rebuild every
  number of `analysis/facts.json` and every table, and compare. `scripts/reproduce.sh` runs it.
- **Level B** (full recomputation, Linux x86-64, about 250–300 vCPU-hours): install the pinned environment, Multiwfn
  and ParaView with `scripts/setup_linux.sh` and run `scripts/run_case.sh` for both cases.

Commands, expected outputs and tolerances are in `docs/REPRODUCE.md`.

## Software and citations

Please cite the programs that produced the data:

- PySCF: Q. Sun et al., *WIREs Comput. Mol. Sci.* 8 (2018) e1340 (doi:10.1002/wcms.1340); Q. Sun et al.,
  *J. Chem. Phys.* 153 (2020) 024109 (doi:10.1063/5.0006074); Q. Sun et al., *J. Chem. Phys.* 165 (2026)
  (doi:10.1063/5.0337441).
- geomeTRIC: L.-P. Wang, C. Song, *J. Chem. Phys.* 144 (2016) 214108 (doi:10.1063/1.4952956).
- Multiwfn (required by its licence): T. Lu, F. Chen, *J. Comput. Chem.* 33 (2012) 580–592; T. Lu, *J. Chem. Phys.*
  161 (2024) 082503.
- libxc (through PySCF), ParaView (3-D renders of the pipeline), matplotlib and PyVista (figures).

## Licences

- Code in `scripts/`: MIT License (`LICENSE`).
- Data and documentation: Creative Commons Attribution 4.0 International (`LICENSE-DATA`).
- The crystal structure is not redistributed and is not covered by these licences; material derived from the CCDC
  deposition is not covered either.

## The crystal structure

The structure was determined by Dhanabal, Amirthaganesan, Dhandapani and Das, *Optik* 124 (2013) 1966–1970,
[doi:10.1016/j.ijleo.2012.06.018](https://doi.org/10.1016/j.ijleo.2012.06.018), and deposited as CCDC 835950
([doi:10.5517/ccx1w4y](https://doi.org/10.5517/ccx1w4y)). The CIF is not redistributed here. Obtain it free of charge
from the CCDC (https://www.ccdc.cam.ac.uk/structures/), save it as `structure/BTPC.cif`, and let
`scripts/dft/build_cluster.py` rebuild the ion pair from it. `structure/README.md` gives the checksums of the CIF the
calculations used and of the rebuilt geometry. A fresh download will not match the CIF checksum byte for byte: the
CCDC writes the download date into each CIF it serves, and the copy used here has CRLF line ends.
`structure/README.md` shows how to compare the content; the checksum of the rebuilt ion pair is the check that
counts.

## How to cite

Cite the dataset as given in `CITATION.cff`. The DOI and the article reference are added at release.

## Contact

- DFT data: Anantha Raman Lakshmipathi, Rajalakshmi Institute of Technology, Chennai, India (raman3366@gmail.com;
  ORCID [0000-0003-3465-286X](https://orcid.org/0000-0003-3465-286X)).
- The article: B. Sahaya Infant Lasalle (corresponding author), SSN Research Centre, Sri Sivasubramaniya Nadar
  College of Engineering, Kalavakkam, India.
