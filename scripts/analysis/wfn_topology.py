"""
Density-based descriptors of the N4-H4N...O1 contact from the stored wavefunctions (final.molden), in numpy.

Local post-processing only: the Kohn-Sham orbitals written by PySCF at the end of each optimisation are
re-evaluated on points with papers/_shared/figkit/molden.py (validated against the PySCF cubes), and
everything else is finite differences of those orbital values.  No quantum-chemistry calculation is run.

    rho, grad rho, Laplacian of rho       sums over the occupied orbitals (n_i = 2)
    G  = 1/2 sum_i n_i |grad phi_i|^2     positive-definite kinetic energy density
    V  = 1/4 lap(rho) - 2 G               local virial theorem
    H  = G + V                            total energy density
    ELF (Becke) = 1 / (1 + (D/D0)^2), D = G - |grad rho|^2 / (8 rho) (+1e-5, Multiwfn ELF_addminimal),
                  D0 = 3/10 (3 pi^2)^(2/3) rho^(5/3)

Validation (run this file directly): rho against density.cube (PySCF cubegen), ELF against the Multiwfn plane
file elf_benzotriazolium_ring_plane.txt, and sign(lambda2)rho at the bond critical point against the lowest-RDG
point of the NCI region between H4N and O1.  build_facts.py imports this module.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # scripts/ (paths.py)
import paths  # noqa: E402
ROOT = paths.legacy_root(quiet=True)   # the dataset files under their workspace paths
sys.path.insert(0, str(ROOT / 'papers/_shared'))
from figkit import molden  # noqa: E402

BOHR = 0.529177210903
H_FD = 5e-4          # Bohr, finite-difference step for orbital derivatives
CF = 0.3 * (3 * math.pi ** 2) ** (2.0 / 3.0)


class Wfn:
    def __init__(self, path):
        self.mf = molden.read_molden(path)
        self.occ_idx = np.flatnonzero(self.mf.mo_occ > 1e-8)
        self.n = self.mf.mo_occ[self.occ_idx]
        self.C = self.mf.mo_coeff[:, self.occ_idx]

    def phi(self, P):
        return self.mf.ao_on_points(np.atleast_2d(P)) @ self.C          # (npts, nocc)

    def local(self, P):
        """rho, grad rho (3), Laplacian, G at points P (Bohr); orbital derivatives by central differences."""
        P = np.atleast_2d(np.asarray(P, float))
        f0 = self.phi(P)
        grad = np.zeros((3,) + f0.shape)
        lap = np.zeros(f0.shape)
        for k in range(3):
            e = np.zeros(3); e[k] = H_FD
            fp, fm = self.phi(P + e), self.phi(P - e)
            grad[k] = (fp - fm) / (2 * H_FD)
            lap += (fp - 2 * f0 + fm) / H_FD ** 2
        rho = (self.n * f0 ** 2).sum(1)
        grho = 2 * np.einsum('i,kpi,pi->pk', self.n, grad, f0)
        g2 = (grad ** 2).sum(0)                                             # |grad phi_i|^2
        G = 0.5 * (self.n * g2).sum(1)
        laprho = 2 * (self.n * (g2 + f0 * lap)).sum(1)
        return rho, grho, laprho, G

    def elf(self, P):
        rho, grho, _, G = self.local(P)
        D = G - (grho ** 2).sum(1) / (8 * rho) + 1e-5
        D0 = CF * rho ** (5.0 / 3.0)
        return 1.0 / (1.0 + (D / D0) ** 2)

    def hessian_rho(self, x, h=2e-3):
        Hm = np.zeros((3, 3))
        for k in range(3):
            e = np.zeros(3); e[k] = h
            gp = self.local(x + e)[1][0]; gm = self.local(x - e)[1][0]
            Hm[:, k] = (gp - gm) / (2 * h)
        return 0.5 * (Hm + Hm.T)

    def bcp(self, a, b, n_line=201, tol=1e-10, maxit=50):
        """(3,-1) critical point of rho between nuclei a and b (Bohr): start at the density minimum on the line."""
        t = np.linspace(0.15, 0.85, n_line)
        line = a[None, :] + t[:, None] * (b - a)[None, :]
        x = line[int(np.argmin(self.local(line)[0]))].copy()
        for _ in range(maxit):
            g = self.local(x)[1][0]
            Hm = self.hessian_rho(x)
            step = -np.linalg.solve(Hm, g)
            x = x + step
            if np.linalg.norm(step) < 1e-7 and np.linalg.norm(g) < tol ** 0.5:
                break
        rho, g, lap, G = (v[0] for v in self.local(x))
        lam = np.linalg.eigvalsh(self.hessian_rho(x))
        V = lap / 4 - 2 * G
        return {'x': x, 'rho': float(rho), 'grad': float(np.linalg.norm(g)), 'lap': float(lap), 'G': float(G),
                'V': float(V), 'H': float(G + V), 'lambda': lam, 'sign_l2_rho': float(np.sign(lam[1]) * rho),
                'VG': float(abs(V) / G), 'elf': float(self.elf(x)[0])}


if __name__ == '__main__':
    import json
    RAW = ROOT / 'data/r2/molecules/btpc/results'
    case = json.load(open(ROOT / 'molecules/btpc/cases/btpc_pcm.json'))
    LAB = [l.split('#')[0] for l in case['cif_labels']]; L = {l: i for i, l in enumerate(LAB)}
    for c in ('btpc_gas', 'btpc_pcm'):
        w = Wfn(RAW / c / 'final.molden')
        X = w.mf.coords
        # 1. density against the PySCF cube
        with open(RAW / c / 'density.cube') as fh:
            fh.readline(); fh.readline()
            nat, *org = fh.readline().split(); org = np.array([float(v) for v in org[:3]])
            ax = [fh.readline().split() for _ in range(3)]; n = [int(a[0]) for a in ax]
            sp = np.array([[float(v) for v in a[1:4]] for a in ax])
            for _ in range(abs(int(nat))): fh.readline()
            cube = np.fromstring(fh.read(), sep=' ').reshape(n)
        rng = np.random.default_rng(1)
        idx = np.argwhere(cube > 1e-3); idx = idx[rng.choice(len(idx), 400, replace=False)]
        P = org + idx @ sp
        r = w.local(P)[0]; ref = cube[tuple(idx.T)]
        print(c, 'rho vs cube: max rel dev', float(np.abs(r / ref - 1).max()), 'max abs', float(np.abs(r - ref).max()))
        # 2. ELF against Multiwfn plane
        pl = np.loadtxt(RAW / c / 'elf_benzotriazolium_ring_plane.txt')
        sel = pl[rng.choice(len(pl), 400, replace=False)]
        sel = sel[sel[:, 5] > 0.02]
        e = w.elf(sel[:, :3] / BOHR)
        print(c, 'ELF vs Multiwfn plane: max abs dev', float(np.abs(e - sel[:, 5]).max()), 'n', len(sel))
        # 3. bond critical point H4N...O1
        b = w.bcp(X[L['H4N']], X[L['O1']])
        print(c, {k: (np.round(v, 5) if not isinstance(v, float) else round(v, 5)) for k, v in b.items()})
        print('   d(H4N-BCP) %.3f A, d(BCP-O1) %.3f A' % (np.linalg.norm(b['x'] - X[L['H4N']]) * BOHR, np.linalg.norm(b['x'] - X[L['O1']]) * BOHR))
