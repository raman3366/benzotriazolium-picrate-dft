# Crystal structure and ion-pair specification

- Structure: benzotriazolium picrate, C₆H₆N₃⁺·C₆H₂N₃O₇⁻ (C₁₂H₈N₆O₇), monoclinic P2₁/c, Z = 4, 298 K; the CIF calls the
  compound BTZP.
- Source: T. Dhanabal, G. Amirthaganesan, M. Dhandapani and S. K. Das, *Optik* 124 (2013) 1966–1970,
  doi:10.1016/j.ijleo.2012.06.018; deposition CCDC 835950, doi:10.5517/ccx1w4y.
- The CIF is not included (the CCDC licence does not allow redistribution of CSD data). Download it free of charge
  from the CCDC (https://www.ccdc.cam.ac.uk/structures/) and save it here as `BTPC.cif` (the file name is written
  into the rebuilt geometry, so keep it). The scripts look for `structure/BTPC.cif`; the environment variable
  `BTPC_CIF` can name another location.
- Your copy will not have the SHA-256 of the copy the calculations used. The CCDC writes the download date into
  every CIF it serves (`_audit_update_record`; the copy used here reads "2011-07-21 deposited with the CCDC." and
  "2026-06-25 downloaded from the CCDC." on one line), and that copy has CRLF line ends. Compare the content
  without the download line and the carriage returns instead:

  ```
  tr -d '\r' < BTPC.cif | grep -v 'downloaded from the CCDC' | sha256sum
  9d93fe1c427457780670219100894a7ff2b48af9050038d552bc9e82104c67d9   (expected)
  17065c7ac118c838b386de8448cb7d631073c64e02cd2cd4d7e929e34c5349b5   (sha256 of the copy used, as downloaded)
  ```

  This comparison is a guide only: `scripts/reproduce.sh` reports it but does not stop on it. The check that counts
  is the rebuilt ion pair below, whose SHA-256 does not depend on the download date or the line ends.

## Hydrogen-bonded ion pair (33 atoms)

The two ions of the asymmetric unit are not hydrogen-bonded to each other (their closest approach is O4···N5
2.99 Å, and N5 carries no hydrogen), so the pair is assembled from symmetry images:

| Moiety | Atoms (0-based) | Symmetry operation (CIF code) |
|---|---|---|
| picrate anion | 0–17 | x, y, z (1_555) |
| benzotriazolium cation | 18–32 | 1 + x, ½ − y, ½ + z (4_666: operation 4 of the CIF, x, −y − ½, z − ½, plus one cell along each axis) |

The pair keeps the bond the cation donates to the phenolate oxygen, N4–H4N···O1: H···O 1.775 Å, N···O 2.623 Å,
N–H···O 154.5° in the X-ray geometry. The cation's second N–H (N6–H6N) donates to a further picrate and is not part of
the pair. The atom order and labels are `cif_labels` in `cases/btpc_gas.json` and `cases/btpc_pcm.json`; labels with
the suffix `#4_666` belong to the symmetry-generated cation.

## Rebuilding the starting geometry

```
python scripts/dft/build_cluster.py structure/BTPC.cif btpc --title Benzotriazolium-picrate
python scripts/dft/curate_cases.py
```

The first command writes `structure/btpc_cluster.xyz` (Cartesian X-ray coordinates of the pair, Å) and a draft case
file `structure/btpc_gas.json`; the second adds the names, planes and notes. Expected results (the same with a CIF
of another download date or with LF line ends):

```
sha256  dfa71f6f1a4b15cf9b9f8f4cf4240f88ef284b073fef5971cfe9c3e5a0f95476  structure/btpc_cluster.xyz
structure/btpc_gas.json is byte-identical to cases/btpc_gas.json
```

`btpc_cluster.xyz` is the starting geometry named in the case files (`"xyz": "molecules/btpc/input/btpc_cluster.xyz"`,
the path in the project workspace; `scripts/paths.py` maps it here). Both rebuilt files are git-ignored. The
final optimisations started from `calculations/<case>/geometry/preopt.xyz`, not from this geometry (see README).
