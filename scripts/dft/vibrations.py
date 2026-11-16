"""
Harmonic frequencies from a Cartesian Hessian, with an explicit translation/rotation projection.

Why this exists: PySCF's `thermo.harmonic_analysis` selects which 6 modes to discard in a way that fails for
floppy clusters.  For PZC's gas-phase structure it discarded three genuinely imaginary modes (-93.7, -82.9,
-42.4 cm-1) and kept three spurious zero modes in their place, then reported "0 imaginary" -- certifying a
saddle point as a minimum, which is exactly the claim these deliverables rest on.

This routine instead builds the 6 translation/rotation vectors in mass-weighted coordinates, projects them
out with a QR-orthonormalised projector, and diagonalises what remains.  Validated against PySCF on 3ATB,
where the two agree to 0.1 cm-1 across all 69 modes.
"""
import numpy as np
from pyscf.data import elements

AU2CM = 219474.6313702
AMU2ME = 1822.888486209


def frequencies(mol, hess):
    """(freqs_cm1, n_imaginary, tr_residual_cm1) from a (natm, natm, 3, 3) Cartesian Hessian."""
    natm = mol.natm
    mass = np.array([elements.MASSES[mol.atom_charge(i)] for i in range(natm)]) * AMU2ME
    x = mol.atom_coords()

    Hm = hess.transpose(0, 2, 1, 3).reshape(3 * natm, 3 * natm)
    w = np.repeat(mass ** -0.5, 3)
    Hm = Hm * w[:, None] * w[None, :]
    Hm = 0.5 * (Hm + Hm.T)

    com = (mass[:, None] * x).sum(0) / mass.sum()
    r = x - com
    TR = np.zeros((6, 3 * natm))
    for i in range(natm):
        s = np.sqrt(mass[i])
        TR[0, 3 * i + 0] = TR[1, 3 * i + 1] = TR[2, 3 * i + 2] = s
        for k, axis in enumerate(np.eye(3)):
            TR[3 + k, 3 * i:3 * i + 3] = s * np.cross(axis, r[i])
    q, _ = np.linalg.qr(TR.T)
    P = np.eye(3 * natm) - q @ q.T

    ev = np.linalg.eigvalsh(P @ Hm @ P)
    f = np.sort(np.sign(ev) * np.sqrt(np.abs(ev)) * AU2CM)
    # the 6 modes closest to zero are the projected translations/rotations
    order = np.argsort(np.abs(f))
    tr_idx = set(order[:6].tolist())
    vib = np.array([v for i, v in enumerate(f) if i not in tr_idx])
    tr_res = float(np.abs(f[order[:6]]).max())
    return vib, int((vib < 0).sum()), tr_res
