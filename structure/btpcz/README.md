# Second crystal form: P2₁/n (COD 2230422)

Benzotriazolium picrate crystallises in two forms. The calculations of the article use the P2₁/c form (CCDC 835950,
`structure/README.md`, cases `btpc_gas` and `btpc_pcm`). The P2₁/n form below was added during the revision of the
article as cases `btpcz_gas` and `btpcz_pcm` (`cases/`, `calculations/btpcz_gas/`, `calculations/btpcz_pcm/`;
`docs/revision.md`).

| Item | Value |
|---|---|
| Structure | benzotriazolium picrate, C₆H₆N₃⁺·C₆H₂N₃O₇⁻ (C₁₂H₈N₆O₇), monoclinic P2₁/n (No. 14), a 14.4113(15), b 3.7608(4), c 24.941(3) Å, β 90.598(2)°, V 1351.7(3) Å³, Z = 4, 298 K, R1 0.0616 (values from the CIF) |
| Article | B. Zeng, J. Li, G.-d. Wang, *Benzotriazolium picrate*, Acta Cryst. E67 (2011) o1464, [doi:10.1107/S1600536811018307](https://doi.org/10.1107/S1600536811018307) |
| Source of the CIF | Crystallography Open Database entry 2230422, <https://www.crystallography.net/cod/2230422.cif> |
| CIF used | downloaded 20 Nov 2026, unmodified, 14,909 bytes, COD file revision `$Revision: 301806 $` (`$Date: 2025-08-19 $` in the header) |
| SHA-256 of the CIF used | `5a662db05f379cb539f14026df38c99fca4be7ea5b03d1e7efb60c05216a804d` |
| Licence | The COD dedicates its data to the public domain (CC0) and asks users to acknowledge the original authors; cite Zeng et al. (2011) wherever this model is used |

The CIF is not included in this repository. Download it from the address above and save it as
`structure/2230422.cif` (CIF files are git-ignored here). A later COD file revision may differ from the copy used;
the checksum of the rebuilt ion pair below is the check that counts.

## Hydrogen-bonded ion pair (33 atoms)

- The asymmetric unit is itself the hydrogen-bonded ion pair: both moieties are taken at x, y, z (1_555), and no
  symmetry operation is needed (the P2₁/c pair of CCDC 835950 needs the 4_666 image of the cation).
- The pair is defined by the N6–H6···O1 hydrogen bond to the phenolate oxygen: H···O 1.870 Å, N···O 2.770 Å,
  N–H···O 160.0° in the X-ray geometry, the `_geom_hbond` entry of the CIF (1.87(3) Å, 2.770(3) Å, 160(3)°,
  symmetry code `.`). For comparison, the P2₁/c pair: N4–H4N···O1, 1.775 Å, 2.623 Å, 154.5°.
- The cation's other N–H (N4–H4) donates a bifurcated bond to O1 and to the ortho-nitro oxygen O2 of a neighbouring
  picrate (H···O 2.035 and 2.161 Å, N···O 2.799 and 2.791 Å). Without `--require-hbond`, `build_cluster.py`
  assembles that pair instead (the cation at the 2_646 image), because two contacts outscore one. The N6–H6···O1 pair
  was chosen because its bond is shorter and more linear, it is the single charge-assisted bond to the phenolate
  oxygen that also defines the P2₁/c pair, and the CIF lists it within the asymmetric unit.
- Hydrogen positions are the refined X-ray values (N4–H4 0.85 Å, N6–H6 0.94 Å; C–H riding at 0.93 Å).
- Atom order and labels: `cif_labels` in `cases/btpcz_gas.json` and `cases/btpcz_pcm.json`; picrate atoms 0–17,
  cation atoms 18–32. The watched proton is H6 on N6, acceptor O1 (`proton` in the case files).

`btpcz_cluster.xyz` in this folder is the starting geometry of both cases (Cartesian X-ray coordinates, Å), named
in the case files as `molecules/btpcz/input/btpcz_cluster.xyz`, its path in the project workspace. It is derived from
the COD entry and is therefore included, unlike the P2₁/c pair.

## Rebuilding the starting geometry and the case files

```sh
python scripts/dft/build_cluster.py structure/2230422.cif btpcz --require-hbond N6:O1 \
    --title "Benzotriazolium picrate, P2_1/n polymorph (Zeng et al. 2011, COD 2230422)"
sha256sum structure/btpcz_cluster.xyz structure/btpcz/btpcz_cluster.xyz
# expected for both: e96dc32b31a9f8d19f9225c4040a1015b2b2ecf0436d00b5b68e45b073a389e9
```

The first command writes `structure/btpcz_cluster.xyz` and a draft case file `structure/btpcz_gas.json` (both
git-ignored). `scripts/revision1/make_btpcz_cases.py` then adds the names, ring planes, the watched proton and the
notes, and writes the IEF-PCM (water) twin. It expects the project layout: the draft as
`molecules/btpcz/cases/btpcz_gas.json` and the geometry as `molecules/btpcz/input/btpcz_cluster.xyz` below the
directory that holds `dft/revision/make_btpcz_cases.py`.

Checked on 23 Nov 2026 (Python 3.12.3 with `environment/paper-requirements.txt`): the rebuilt
`btpcz_cluster.xyz` is byte-identical to the one here, and `make_btpcz_cases.py` run on the rebuilt draft writes
case files byte-identical to `cases/btpcz_gas.json` and `cases/btpcz_pcm.json`.
