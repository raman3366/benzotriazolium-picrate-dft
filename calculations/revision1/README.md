# Revision calculations: job outputs

The outputs of the 57 jobs run during the revision of the article (November 2026): 606 files,
47.5 MiB. `docs/revision.md` describes the methods, the environment, the hosts, the path placeholders
and how to run a job again.

## Layout

```
calculations/revision1/
  README.md                this file
  runner/vm1/, runner/vm2/ job-queue records of the two compute VMs: status.json (per job: expanded command, host,
                           start and end in IST, wall hours, threads, exit code), runner.log, runner.out
  <job>/                   one folder per job (57)
    queue.log              the queue's log of the job (command, then the job's output)
    result.json, log.txt   the job's record and log (not for c0_prep, c0_btpcz_*_analysis and fdh1-fdh5)
    ...                    every other file the job wrote, under the same names, except the PySCF checkpoints
                           (chk/), scratch and the geomeTRIC constraint files kept with the checkpoints
```

Jobs of the first VM came from the queue folder `<run-dir>/queue/out/<job>/` of that VM, the two `*_f080` jobs from
the second VM's; both are here under the job name. The job lists that define every job are
`scripts/revision1/jobs_vm.json` (first VM, 55 jobs) and `scripts/revision1/jobs_vm2.json` (second VM, 2 jobs).
Threads are those of the job list; wall hours are from `status.json`. Times in the logs are IST (UTC+5:30).

Paths in the text files: `<run-dir>` is the working folder on the compute VM, `<vm-1>` and `<vm-2>` the VM host
names, `<user>` the VM user (`docs/revision.md`, "Paths in the files").

## Jobs

### Second crystal form (P2₁/n, COD 2230422)

| Job | What it computes | Script | VM | Threads | Wall h | Files |
|---|---|---|---|---:|---:|---:|
| `c0_prep` | Host preparation: links the environment and Multiwfn for `run_case.sh`; prints the PySCF and geomeTRIC versions | `(shell)` | vm1 | 1 | 0.00 | 1 |
| `c0_btpcz_gas_core` | Gas phase: B3LYP/6-31G(d) pre-optimisation, optimisation, proton check, finite-difference Hessian (0.005 Bohr), frontier orbitals and cubes; case results in `calculations/btpcz_gas/` | `c0_case.py` | vm1 | 4 | 27.88 | 54 |
| `c0_btpcz_pcm_core` | IEF-PCM water: pre-optimisation, optimisation, proton check, frontier orbitals and cubes; case results in `calculations/btpcz_pcm/` | `c0_case.py` | vm1 | 4 | 14.86 | 17 |
| `c0_btpcz_pcm_softmode` | IEF-PCM minimum: Hessian check along the 10 softest modes of the gas-phase Hessian of the same form (`softmode_check.json`) | `c0_case.py` | vm1 | 4 | 4.62 | 7 |
| `c0_btpcz_pcm_fdhess` | IEF-PCM minimum: full finite-difference Hessian (0.005 Bohr), with `fdh1`–`fdh5`; `hessian.npy` and `freq.json` in `calculations/btpcz_pcm/frequencies/` | `c0_case.py` | vm1 | 4 | 16.82 | 59 |
| `fdh1` | Helper worker of `c0_btpcz_pcm_fdhess`: column chunks from the end of the list (`fd_helper.py`) | `fd_helper.py` | vm1 | 4 | 6.21 | 1 |
| `fdh2` | As `fdh1` | `fd_helper.py` | vm1 | 4 | 5.01 | 1 |
| `fdh3` | As `fdh1` | `fd_helper.py` | vm1 | 4 | 4.92 | 1 |
| `fdh4` | As `fdh1` | `fd_helper.py` | vm1 | 4 | 4.92 | 1 |
| `fdh5` | As `fdh1` | `fd_helper.py` | vm1 | 4 | 2.38 | 1 |
| `c0_btpcz_gas_analysis` | Gas phase: Multiwfn MEP, NCI, ELF and LOL, plots and `RESULTS.md` through `run_case.sh` (no ParaView renders) | `run_case.sh` | vm1 | 4 | 0.05 | 1 |
| `c0_btpcz_pcm_analysis` | IEF-PCM: as `c0_btpcz_gas_analysis` | `run_case.sh` | vm1 | 4 | 0.04 | 1 |
| `charges_btpcz` | Hirshfeld and CM5 charges (PySCF and Multiwfn) of both minima, from their `final.chk` | `charges.py` | vm1 | 4 | 0.01 | 5 |
| `headline_btpcz` | Cross single points and permittivity series at both minima: where the LUMO lies | `headline.py` | vm1 | 4 | 2.06 | 5 |

### Where the LUMO lies (published ion pair)

| Job | What it computes | Script | VM | Threads | Wall h | Files |
|---|---|---|---|---:|---:|---:|
| `headline_gasgeom` | Gas-phase minimum: single points at ε = 1, 2, 4, 8, 20, 78.3553 and one refinement point; frontier-orbital cation/picrate shares; ε at which the LUMO moves to the picrate | `headline.py` | vm1 | 4 | 1.17 | 4 |
| `headline_pcmgeom` | IEF-PCM minimum: as `headline_gasgeom` | `headline.py` | vm1 | 4 | 1.23 | 4 |
| `deltascf_gasgeom` | ΔSCF at the gas-phase minimum, gas phase and IEF-PCM: vertical N+1 and N−1 (UKS) states, where the electron and the hole go | `deltascf.py` | vm1 | 4 | 3.11 | 5 |
| `deltascf_pcmgeom` | ΔSCF at the IEF-PCM minimum, IEF-PCM and gas phase | `deltascf.py` | vm1 | 4 | 3.23 | 5 |
| `charges_stored` | Hirshfeld and CM5 charges (PySCF and Multiwfn) from the published wavefunctions, with Mulliken and Löwdin | `charges.py` | vm1 | 4 | 0.01 | 5 |
| `sp8_pbe0` | PBE0 at both minima, gas phase and IEF-PCM: energies and LUMO location | `headline.py` | vm1 | 4 | 0.65 | 5 |
| `sp8_camb3lyp` | CAM-B3LYP at both minima, gas phase and IEF-PCM: energies and LUMO location | `headline.py` | vm1 | 4 | 0.71 | 5 |

### Dispersion (B3LYP-D3(BJ))

| Job | What it computes | Script | VM | Threads | Wall h | Files |
|---|---|---|---|---:|---:|---:|
| `opt_d3_gas` | Gas-phase re-optimisation seeded with the published Hessian; analysis and `final.molden` at the new minimum | `opt_d3.py` | vm1 | 4 | 0.94 | 8 |
| `opt_d3_pcm` | IEF-PCM re-optimisation to all five criteria, seeded with the published Hessian; analysis and `final.molden` | `opt_d3.py` | vm1 | 4 | 23.34 | 8 |
| `softmode_d3_gas` | Hessian check at the D3(BJ) gas-phase minimum along the 10 softest modes (21 gradients, 0.02 Bohr) | `softmode_check.py` | vm1 | 4 | 2.46 | 4 |
| `softmode_d3_pcm` | Hessian check at the D3(BJ) IEF-PCM minimum along the 10 softest modes | `softmode_check.py` | vm1 | 4 | 5.23 | 4 |

### Proton transfer along N4–H4N···O1

| Job | What it computes | Script | VM | Threads | Wall h | Files |
|---|---|---|---|---:|---:|---:|
| `scan_gas_up` | Relaxed scan of r(H4N···O1), gas phase, B3LYP: from the minimum outward to 1.70 Å | `proton_scan.py` | vm1 | 4 | 10.36 | 19 |
| `scan_gas_down` | The same scan from the minimum inward to 1.08 Å | `proton_scan.py` | vm1 | 4 | 18.77 | 28 |
| `scan_gas_neutral` | The same scan from the neutral form (proton on O1) up to 1.28 Å; overlaps `scan_gas_down` | `proton_scan.py` | vm1 | 4 | 22.83 | 26 |
| `scan_gas_oside` | Two further points on the oxygen side, r(H···O) 0.90 and 0.94 Å | `proton_scan.py` | vm1 | 4 | 6.19 | 11 |
| `neutral_start_gas` | Unconstrained optimisation from the neutral form (O1–H 0.98 Å), gas phase | `neutral_start.py` | vm1 | 4 | 4.65 | 7 |
| `neutral_start_pcm` | Unconstrained optimisation from the neutral form, IEF-PCM water | `neutral_start.py` | vm1 | 4 | 18.05 | 7 |
| `sp_scan_wb97xd3bj` | ωB97X-D3(BJ) single points at the geometries of the three scan chains | `singlepoints.py` | vm1 | 4 | 1.64 | 20 |
| `sp_scan_m062x` | M06-2X single points at the scan geometries, with a grid-level-5 check | `singlepoints.py` | vm1 | 4 | 2.49 | 21 |
| `sp_scan_pbe0` | PBE0 single points at the scan geometries | `singlepoints.py` | vm1 | 4 | 1.35 | 20 |
| `sp_scan_camb3lyp` | CAM-B3LYP single points at the scan geometries | `singlepoints.py` | vm1 | 4 | 1.64 | 20 |
| `sp_scan_dfmp2` | DF-MP2 and SCS-MP2 (RI basis aug-cc-pVTZ-RI) at the scan geometries | `singlepoints.py` | vm1 | 4 | 0.59 | 20 |
| `sp_oside_all` | The five methods above at the two oxygen-side points | `singlepoints.py` | vm1 | 4 | 1.10 | 13 |
| `sp_min_pcm` | The five methods in IEF-PCM water at the published IEF-PCM minimum (salt) | `singlepoints.py` | vm1 | 4 | 0.05 | 11 |
| `sp_neutral_pcm` | The five methods in IEF-PCM water at the IEF-PCM neutral minimum of `neutral_start_pcm` | `singlepoints.py` | vm1 | 4 | 1.19 | 11 |
| `sp_gas_pcmgeoms` | Gas-phase DF-MP2 at the IEF-PCM salt and neutral geometries | `singlepoints.py` | vm1 | 4 | 0.12 | 8 |
| `dvr_gas_b3lyp` | Proton vibrational levels (sinc-DVR) on the B3LYP gas-phase profile, with analytic self-tests | `dvr.py` | vm1 | 1 | 0.00 | 6 |
| `dvr_gas_all` | Proton vibrational levels on the profiles of every method | `dvr.py` | vm1 | 1 | 0.01 | 12 |
| `cc_salt_631gd` | CCSD(T) − MP2 correction, 6-31G*, frozen natural orbitals 50 and 65 %, IEF-PCM water: salt minimum | `cc_correction.py` | vm1 | 10 | 11.64 | 4 |
| `cc_neutral_631gd` | The same, 50 %: neutral minimum | `cc_correction.py` | vm1 | 4 | 1.93 | 4 |
| `cc_neutral_631gd_f065` | The same, 65 %: neutral minimum | `cc_correction.py` | vm1 | 4 | 4.19 | 4 |
| `cc_salt_631gd_f080` | The same, 80 %: salt minimum | `cc_correction.py` | vm2 | 16 | 5.07 | 4 |
| `cc_neutral_631gd_f080` | The same, 80 %: neutral minimum | `cc_correction.py` | vm2 | 16 | 5.14 | 4 |
| `cc_salt_631ppgd` | CCSD(T) − MP2 correction, 6-31+G**, 50 %, IEF-PCM water: salt minimum | `cc_correction.py` | vm1 | 8 | 12.51 | 4 |
| `cc_neutral_631ppgd` | The same: neutral minimum | `cc_correction.py` | vm1 | 8 | 12.08 | 4 |
| `cc_gas_631gd` | CCSD(T) − MP2 correction, 6-31G*, 50 %, gas phase, at both IEF-PCM geometries | `cc_correction.py` | vm1 | 6 | 9.20 | 5 |

### Inter-ion torsion C1–O1···N4–C12 (gas phase)

| Job | What it computes | Script | VM | Threads | Wall h | Files |
|---|---|---|---|---:|---:|---:|
| `scan_tors_b3lyp` | Relaxed scan, B3LYP, −180° to −60° in 15° steps | `conf_scan.py` | vm1 | 4 | 25.55 | 32 |
| `scan_tors_d3` | Relaxed scan, B3LYP-D3(BJ) from the D3(BJ) minimum, −120° to −180° | `conf_scan.py` | vm1 | 4 | 7.81 | 20 |
| `scan_tors_d3_up` | The same scan, −105° to −75° | `conf_scan.py` | vm1 | 12 | 2.54 | 14 |
| `scan_tors_d3_up60` | The same scan, −60°, started from the −90° point | `conf_scan.py` | vm1 | 4 | 3.43 | 8 |

### Excited states (TD-DFT, Tamm–Dancoff, 20 singlets)

| Job | What it computes | Script | VM | Threads | Wall h | Files |
|---|---|---|---|---:|---:|---:|
| `td_free` | Free picrate anion and neutral 1H-benzotriazole optimised in IEF-PCM water; B3LYP and CAM-B3LYP | `tddft.py` | vm1 | 4 | 3.78 | 11 |
| `td_pair_b3lyp` | Ion pair at both minima (gas phase; IEF-PCM with non-equilibrium response), B3LYP; hole and electron fractions | `tddft.py` | vm1 | 4 | 8.15 | 5 |
| `td_pair_camb3lyp` | As `td_pair_b3lyp`, CAM-B3LYP | `tddft.py` | vm1 | 4 | 13.09 | 5 |
