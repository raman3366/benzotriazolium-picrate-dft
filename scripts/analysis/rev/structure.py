"""Revision facts: structure (keys rev.struct.*).

What is extracted (all inputs read-only):
  A. B3LYP-D3(BJ) re-optimisations of the P2_1/c ion pair, gas phase and IEF-PCM (water): opt_d3_gas, opt_d3_pcm.
  B. Soft-mode Hessian checks at those minima (softmode_d3_gas, softmode_d3_pcm) and at the P2_1/n IEF-PCM minimum
     (c0_btpcz_pcm_softmode).
  C. Relaxed gas-phase scans of the pseudo-torsion C1-O1...N4-C12, B3LYP (scan_tors_b3lyp) and B3LYP-D3(BJ)
     (scan_tors_d3 = -180..-120, scan_tors_d3_up = -105..-75, scan_tors_d3_up60 = -60).
  D. Torsion angles of the crystal pair, the published B3LYP minima and the D3(BJ) minima (Reviewer 1, point 5),
     and C1-O1...N4-C12 along the published IEF-PCM optimisation.
  E. The ion pair of the second crystal form, P2_1/n (COD 2230422): c0_btpcz_gas_core, c0_btpcz_pcm_core,
     c0_btpcz_pcm_fdhess, c0_btpcz_pcm_softmode and data/revision_runs/repo_results/btpcz/.

Every value is read from a job output or a structure file, or derived here by the formula stated in its note.
Geometry is recomputed from the xyz files with the definitions of build_facts.py (ring planes: six picrate ring
carbons and the nine ring atoms of the cation; nitro twist: O-N-O plane against the picrate ring plane; X-ray
connectivity: covalent radii x 1.25) and checked against the values the job scripts stored and against the published
facts. Cases: rev_gas / rev_pcm (P2_1/c pair, new calculations), rev_xray (P2_1/c crystal pair), btpc_gas / btpc_pcm
(published B3LYP minima), btpcz_gas / btpcz_pcm / btpcz_xray (P2_1/n pair).
Derived tables papers/btpc/data/derived/rev_struct_*.csv are written only in the full build_facts.py run.
"""
import csv
import hashlib
import math
import re

import numpy as np

from .common import *

R_KCAL = 1.987204259e-3          # molar gas constant, kcal mol-1 K-1 (CODATA 2018; as revision/round1/drafts/torsions.py)
T_ROOM = 298.15
AMU2ME = 1822.888486209
MASS = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999}   # as build_facts.py (the masses of freq.json)
COV = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66}            # X-ray connectivity rule of build_facts.py
BOND_FACTOR = 1.25

PIC_RING = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6']
CAT9 = ['N4', 'N5', 'N6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12']
NITRO = {'N1': ('N1', 'O2', 'O3'), 'N2': ('N2', 'O4', 'O5'), 'N3': ('N3', 'O6', 'O7')}
# torsions of Table S4b (revision/round1/drafts/torsions.py, same atoms and order)
TORSIONS = [('C1', 'O1', 'N4', 'C12'), ('C6', 'C1', 'O1', 'H4N'),
            ('O2', 'N1', 'C2', 'C1'), ('O3', 'N1', 'C2', 'C3'), ('O4', 'N2', 'C4', 'C3'), ('O5', 'N2', 'C4', 'C5'),
            ('O7', 'N3', 'C6', 'C1'), ('O6', 'N3', 'C6', 'C5'),
            ('O1', 'C1', 'C2', 'N1'), ('O1', 'C1', 'C6', 'N3'), ('C6', 'C1', 'C2', 'C3'),
            ('C12', 'N4', 'N5', 'N6'), ('N4', 'N5', 'N6', 'C7'), ('H4N', 'N4', 'N5', 'N6'), ('H6N', 'N6', 'N5', 'N4')]
# agreement with the values the job scripts stored (xyz files carry 10 decimals): distances in A, angles in deg
TOL = {'da': 1e-6, 'dh': 1e-6, 'ha': 1e-6, 'angle': 1e-3, 'tau': 1e-3, 'interplanar': 1e-3}
SCAN_NITRO = [('O2', 'N1', 'C2', 'C1'), ('O4', 'N2', 'C4', 'C3'), ('O7', 'N3', 'C6', 'C1')]
MODEL_TAG = {'xray': 'X-ray', 'gas': 'published B3LYP gas-phase minimum', 'pcm': 'published B3LYP IEF-PCM (water) structure',
             'd3gas': 'B3LYP-D3(BJ) gas-phase minimum', 'd3pcm': 'B3LYP-D3(BJ) IEF-PCM (water) minimum'}
MODEL_CASE = {'xray': 'rev_xray', 'gas': 'btpc_gas', 'pcm': 'btpc_pcm', 'd3gas': 'rev_gas', 'd3pcm': 'rev_pcm'}

# ---------------------------------------------------------------------------------------------------- file helpers
CASE_BTPC = ROOT / 'molecules/btpc/cases/btpc_pcm.json'
CASE_BTPCZ = ROOT / 'molecules/btpcz/cases/btpcz_pcm.json'
XRAY_BTPC = ROOT / 'molecules/btpc/input/btpc_cluster.xyz'
XRAY_BTPCZ = ROOT / 'molecules/btpcz/input/btpcz_cluster.xyz'
STORED = {m: ROOT / f'molecules/btpc/results/btpc_{m}' for m in ('gas', 'pcm')}
RAW_PCM_TRAJ = ROOT / 'data/r2/molecules/btpc/results/btpc_pcm/opt_traj.xyz'
BTPCZ = {m: RUNS / f'repo_results/btpcz/btpcz_{m}' for m in ('gas', 'pcm')}
FACTS_DISK = ROOT / 'papers/btpc/data/facts.json'
DERIVED = ROOT / 'papers/btpc/data/derived'
_disk = {}


def jload(p):
    return json.load(open(p))


def sha256(p):
    return sha256_as_run(p)   # dataset: checksum as the run wrote the file (rev/common.py)


def read_xyz(p):
    L = open(p).read().split('\n')
    n = int(L[0].split()[0])
    return [l.split()[0] for l in L[2:2 + n]], np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])


def read_traj(p):
    L = open(p).read().splitlines()
    out, i = [], 0
    while i < len(L):
        if not L[i].strip():
            i += 1
            continue
        n = int(L[i])
        out.append((L[i + 1], np.array([[float(v) for v in l.split()[1:4]] for l in L[i + 2:i + 2 + n]])))
        i += 2 + n
    return out


def labels_of(case_path):
    c = jload(case_path)
    lab = [l.split('#')[0] for l in c['cif_labels']]
    return c, lab, {l: i for i, l in enumerate(lab)}


def check_symbols(syms, lab, what):
    for s, l in zip(syms, lab):
        if not (l.startswith(s) and not l[len(s)].isalpha()):
            raise AssertionError(f'{what}: element {s} does not match label {l}')
    if len(syms) != len(lab):
        raise AssertionError(f'{what}: {len(syms)} atoms, {len(lab)} labels')


def pub(facts, key):
    """Published value: from the facts being built (full run) or from facts.json on disk (standalone run)."""
    if key in facts:
        return facts[key]['value']
    if not _disk:
        _disk.update(jload(FACTS_DISK))
    return _disk[key]['value']


def close(a, b, tol, what):
    if not abs(a - b) <= tol:
        raise AssertionError(f'cross-check failed: {what}: {a!r} vs {b!r} (tol {tol:g})')


def fl(v, nd=None):
    """Plain-Python list of floats (15 significant digits, as build_facts.add)."""
    return [float(f'{float(x):.15g}') if nd is None else round(float(x), nd) for x in v]


# ---------------------------------------------------------------------------------------------------- geometry
def dist(X, i, j):
    return float(np.linalg.norm(X[i] - X[j]))


def angle(X, i, j, k):
    v1, v2 = X[i] - X[j], X[k] - X[j]
    return math.degrees(math.acos(float(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2))))


def dih(X, a, b, c, d):
    p0, p1, p2, p3 = X[a], X[b], X[c], X[d]
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(float(np.dot(np.cross(b1n, v), w)), float(np.dot(v, w))))


def plane(X, idx):
    c = X[idx].mean(0)
    return c, np.linalg.svd(X[idx] - c)[2][-1]


def plane_angle(n1, n2):
    return math.degrees(math.acos(min(1.0, abs(float(np.dot(n1, n2))))))


def elevation(X, a, b, idx):
    """Angle between the vector a->b and the least-squares plane of atoms idx."""
    _, n = plane(X, idx)
    v = X[b] - X[a]
    return math.degrees(math.asin(abs(float(np.dot(v, n))) / float(np.linalg.norm(v))))


def kabsch(A, B):
    a, b = A - A.mean(0), B - B.mean(0)
    U, _, Vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Rm = Vt.T @ np.diag([1, 1, d]) @ U.T
    return (Rm @ a.T).T, b


def kabsch_rmsd(A, B):
    a, b = kabsch(A, B)
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))


def aligned_step(Xnew, Xold):
    """RMS and max atomic displacement (A) after best-fit superposition (geomeTRIC calc_drms_dmax, align=True)."""
    a, b = kabsch(Xnew, Xold)
    d = np.linalg.norm(a - b, axis=1)
    return float(np.sqrt((d ** 2).mean())), float(d.max())


def tau_norm(t):
    """Torsions of the scan in (-270, 90]: the -180 point stored as +179.999 becomes -180.001."""
    return t - 360.0 if t > 90.0 else t


def pair_geometry(X, I, donor, H, acc='O1', fused='C12'):
    """Hydrogen bond, inter-ion pseudo-torsion C1-O1...donor-fused, ring planes and nitro twists."""
    dh, ha, da = dist(X, I[donor], I[H]), dist(X, I[H], I[acc]), dist(X, I[donor], I[acc])
    c1, n1 = plane(X, [I[a] for a in PIC_RING])
    c2, n2 = plane(X, [I[a] for a in CAT9])
    g = {'dh': dh, 'ha': ha, 'da': da, 'angle': angle(X, I[donor], I[H], I[acc]),
         'tau': dih(X, I['C1'], I['O1'], I[donor], I[fused]), 'interplanar': plane_angle(n1, n2),
         'centroid': float(np.linalg.norm(c1 - c2)), 'proton_on': donor if dh < ha else acc,
         'asym': ha - dh}
    for nm, (N, Oa, Ob) in NITRO.items():
        g[f'twist_{nm}'] = plane_angle(plane(X, [I[N], I[Oa], I[Ob]])[1], n1)
    return g


class XrayRef:
    """X-ray connectivity of one crystal pair (bonds and angles), for agreement statistics as build_facts.py."""

    def __init__(self, syms, X):
        self.syms, self.X = syms, X
        n = len(syms)
        self.bonds = [(i, j) for i in range(n) for j in range(i + 1, n)
                      if np.linalg.norm(X[i] - X[j]) < BOND_FACTOR * (COV[syms[i]] + COV[syms[j]])]
        nb = {i: set() for i in range(n)}
        for i, j in self.bonds:
            nb[i].add(j)
            nb[j].add(i)
        self.angles = [(i, j, k) for j in range(n) for i in nb[j] for k in nb[j] if i < k]
        self.heavy = [i for i in range(n) if syms[i] != 'H']

    def stats(self, Y):
        hv = lambda t: all(self.syms[a] != 'H' for a in t)
        db = [dist(Y, i, j) - dist(self.X, i, j) for i, j in self.bonds if hv((i, j))]
        da = [angle(Y, *t) - angle(self.X, *t) for t in self.angles if hv(t)]
        rms = lambda v: float(np.sqrt(np.mean(np.square(v))))
        return {'bond_rmsd_heavy': rms(db), 'n_bonds_heavy': len(db), 'angle_rmsd_heavy': rms(da),
                'n_angles_heavy': len(da), 'kabsch_rmsd_heavy': kabsch_rmsd(Y[self.heavy], self.X[self.heavy]),
                'n_heavy': len(self.heavy)}


def normal_modes(syms, X_A, H4):
    """Mass-weighted, TR-projected harmonic analysis of a (n, n, 3, 3) Cartesian Hessian (build_facts.modes_for)."""
    n = len(syms)
    x = X_A / BOHR
    m = np.array([MASS[e] for e in syms]) * AMU2ME
    Hm = H4.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)
    w = np.repeat(m ** -0.5, 3)
    Hm = Hm * w[:, None] * w[None, :]
    Hm = 0.5 * (Hm + Hm.T)
    r = x - (m[:, None] * x).sum(0) / m.sum()
    TR = np.zeros((6, 3 * n))
    for i in range(n):
        sq = np.sqrt(m[i])
        TR[0, 3 * i] = TR[1, 3 * i + 1] = TR[2, 3 * i + 2] = sq
        for k, ax in enumerate(np.eye(3)):
            TR[3 + k, 3 * i:3 * i + 3] = sq * np.cross(ax, r[i])
    q, _ = np.linalg.qr(TR.T)
    P = np.eye(3 * n) - q @ q.T
    ev, V = np.linalg.eigh(P @ Hm @ P)
    f = np.sign(ev) * np.sqrt(np.abs(ev)) * AU2CM
    order = np.argsort(np.abs(f))
    tr = set(order[:6].tolist())
    keep = [i for i in np.argsort(f) if i not in tr]
    return x, m, f[keep], V[:, keep]


def rigid_basis(x, m, groups):
    B = []
    for g in groups:
        gm = m[g]
        com = (gm[:, None] * x[g]).sum(0) / gm.sum()
        for k in range(3):
            v = np.zeros((len(x), 3))
            v[g, k] = np.sqrt(gm)
            B.append(v.ravel())
        for ax in np.eye(3):
            v = np.zeros((len(x), 3))
            for i in g:
                v[i] = np.sqrt(m[i]) * np.cross(ax, x[i] - com)
            B.append(v.ravel())
    Q, _ = np.linalg.qr(np.array(B).T)
    return Q


def nitro_torsion_vector(x, m, I, N, C, Os):
    """Mass-weighted unit vector of a rigid turn of the O-N-O group about the C-N bond."""
    ax = x[I[N]] - x[I[C]]
    ax /= np.linalg.norm(ax)
    v = np.zeros_like(x)
    for o in Os:
        v[I[o]] = np.sqrt(m[I[o]]) * np.cross(ax, x[I[o]] - x[I[N]])
    v = v.ravel()
    return v / np.linalg.norm(v)


# ====================================================================================================== build
def build(add, facts):
    full_run = 'gas.e_total' in facts            # build_facts.py run (published facts already added) vs standalone
    case_b, LB, IB = labels_of(CASE_BTPC)
    PIC = sorted(case_b['fragments']['picrate anion'])
    CAT = sorted(case_b['fragments']['benzotriazolium cation'])
    sx, XR = read_xyz(XRAY_BTPC)
    check_symbols(sx, LB, rel(XRAY_BTPC))
    G = {'xray': XR}
    for m in ('gas', 'pcm'):
        s, X = read_xyz(STORED[m] / 'opt.xyz')
        check_symbols(s, LB, rel(STORED[m] / 'opt.xyz'))
        G[m] = X
    for m in ('gas', 'pcm'):
        s, X = read_xyz(job_dir(f'opt_d3_{m}') / 'opt.xyz')
        check_symbols(s, LB, f'opt_d3_{m}/opt.xyz')
        G['d3' + m] = X
    GEO = {k: pair_geometry(X, IB, 'N4', 'H4N') for k, X in G.items()}
    xref = XrayRef(sx, XR)
    # the definitions here reproduce the published facts (Table 2, Table 3)
    for k in ('xray', 'gas', 'pcm'):
        g = GEO[k]
        for mine, fk in (('dh', 'hb_N4_H4N_O1_dh'), ('ha', 'hb_N4_H4N_O1_ha'), ('da', 'hb_N4_H4N_O1_da'),
                         ('angle', 'hb_N4_H4N_O1_angle'), ('interplanar', 'interplanar_picrate_cation'),
                         ('centroid', 'centroid_distance_picrate_cation'), ('twist_N1', 'nitro_twist_N1'),
                         ('twist_N2', 'nitro_twist_N2'), ('twist_N3', 'nitro_twist_N3')):
            close(g[mine], pub(facts, f'{k}.{fk}'), 1e-6, f'{k}.{fk}')
    for k in ('gas', 'pcm'):
        st = xref.stats(G[k])
        for q in ('bond_rmsd_heavy', 'angle_rmsd_heavy', 'kabsch_rmsd_heavy'):
            close(st[q], pub(facts, f'{k}.{q}'), 1e-9, f'{k}.{q}')
    ctx = dict(facts=facts, LB=LB, IB=IB, PIC=PIC, CAT=CAT, G=G, GEO=GEO, xref=xref, rows={})
    _d3_minima(add, ctx)
    _softmode(add, ctx)
    _scans(add, ctx)
    _torsions(add, ctx)
    _btpcz(add, ctx)
    if full_run:
        _write_csvs(ctx)


# ---------------------------------------------------------------------------------------------- A. D3(BJ) minima
def _d3_minima(add, ctx):
    facts, IB, G, GEO = ctx['facts'], ctx['IB'], ctx['G'], ctx['GEO']
    thr = {'energy': pub(facts, 'protocol.opt_threshold_energy'), 'grms': pub(facts, 'protocol.opt_threshold_grad_rms'),
           'gmax': pub(facts, 'protocol.opt_threshold_grad_max'), 'drms': pub(facts, 'protocol.opt_threshold_disp_rms'),
           'dmax': pub(facts, 'protocol.opt_threshold_disp_max')}
    for m in ('gas', 'pcm'):
        job = f'opt_d3_{m}'
        r, src = result(job)
        jd = job_dir(job)
        o, fin = r['items']['opt'], r['items']['final']
        cs, P, k = f'rev_{m}', f'rev.struct.d3_{m}.', 'd3' + m
        model = 'Gas phase' if m == 'gas' else 'IEF-PCM (water)'
        g = GEO[k]
        xyz = rel(jd / 'opt.xyz')
        # integrity: the final single point belongs to this opt.xyz, and the job's own geometry agrees with ours
        if fin['opt_xyz_sha256'] != sha256(jd / 'opt.xyz'):
            raise AssertionError(f'{job}: final single point was not run on opt.xyz')
        for mine, jk in (('da', 'd_N...O_A'), ('dh', 'd_H-N_A'), ('ha', 'd_H...O_A'), ('angle', 'angle_N-H...O_deg'),
                         ('tau', 'torsion_deg'), ('interplanar', 'inter_ring_angle_deg')):
            close(g[mine], o['geometry'][jk], TOL[mine], f'{job} {jk}')
        sc = r['settings']['geometric_criteria']
        for a, b in (('energy', 'convergence_energy'), ('grms', 'convergence_grms'), ('gmax', 'convergence_gmax'),
                     ('drms', 'convergence_drms'), ('dmax', 'convergence_dmax')):
            close(sc[b], thr[a], 0, f'{job} threshold {b} = protocol.opt_threshold_*')
        lvl = (f"B3LYP-D3(BJ)/6-311++G(d,p), {model}, grid {r['settings']['grid_level']}"
               + (f", Lebedev {r['settings']['lebedev_order']}" if m == 'pcm' else '')
               + f", PySCF {r['software']['pyscf']}, geomeTRIC {r['software']['geometric']}")
        add(P + 'level', lvl, '', cs, f'{src}:settings, software', 'level of theory and code of the re-optimisation')
        add(P + 'opt_converged_on', o['converged_on'], '', cs, f'{src}:items.opt.converged_on',
            f"converged = {o['converged']}; settings.all_five_required = {r['settings']['all_five_required']}; "
            f"criteria {r['settings']['criteria_name']}; started from the published B3LYP structure "
            f"{rel(STORED[m] / 'opt.xyz')} with the stored finite-difference Hessian as the initial Hessian "
            f"({r['items']['opt']['seed_hessian']['n_internal_eigenvalues_floored']} internal eigenvalues floored at "
            f"{r['items']['opt']['seed_hessian']['floor_Ha_bohr2']:g} Ha Bohr-2)")
        add(P + 'opt_steps', o['steps'], '', cs, f'{src}:items.opt.steps', 'energy and gradient evaluations, step 0 included')
        hist = o['history']
        last = o['last_step']
        add(P + 'opt_last_dE', last['dE_Ha'], 'Ha', cs, f'{src}:items.opt.last_step.dE_Ha',
            f"energy change of the last step; threshold {thr['energy']:g} Ha (protocol.opt_threshold_energy)")
        add(P + 'opt_final_grad_rms', last['grms'], 'Ha/Bohr', cs, f'{src}:items.opt.last_step.grms',
            f"geomeTRIC norm (RMS of the per-atom gradient norms); threshold {thr['grms']:g}")
        add(P + 'opt_final_grad_max', last['gmax'], 'Ha/Bohr', cs, f'{src}:items.opt.last_step.gmax',
            f"largest per-atom gradient norm; threshold {thr['gmax']:g}")
        traj = read_traj(jd / 'opt_traj.xyz')
        if len(traj) != o['steps'] or float(np.abs(traj[-1][1] - G[k]).max()) > 1e-6:
            raise AssertionError(f'{job}: opt_traj.xyz does not end at opt.xyz or has {len(traj)} frames')
        drms, dmax = aligned_step(traj[-1][1], traj[-2][1])
        note_d = (f"last step, after best-fit superposition of the last two frames of {rel(jd / 'opt_traj.xyz')} "
                  f"(geomeTRIC measures the step this way, calc_drms_dmax align=True); the unaligned Cartesian values "
                  f"stored in result.json, {last['drms_A']:.2e} / {last['dmax_A']:.2e} A, include the overall turn of "
                  f"the pair")
        add(P + 'opt_last_disp_rms', drms, 'Å', cs, f"{rel(jd / 'opt_traj.xyz')} (last two frames)",
            note_d + f"; threshold {thr['drms']:g} A")
        add(P + 'opt_last_disp_max', dmax, 'Å', cs, f"{rel(jd / 'opt_traj.xyz')} (last two frames)",
            note_d + f"; threshold {thr['dmax']:g} A")
        all5 = (abs(last['dE_Ha']) < thr['energy'] and last['grms'] < thr['grms'] and last['gmax'] < thr['gmax']
                and drms < thr['drms'] and dmax < thr['dmax'])
        add(P + 'opt_all_five_met', 'yes' if all5 else 'no', '', cs, f'{src}:items.opt.last_step; opt_traj.xyz',
            'derived: |dE|, gradient RMS/max (geomeTRIC norms) and aligned step RMS/max of the last step against the '
            'protocol.opt_threshold_* values')
        add(P + 'e_total', fin['E_Ha'], 'Ha', cs, f'{src}:items.final.E_Ha',
            f"single point at opt.xyz with conv_tol {r['settings']['conv_tol_final']:g}; the optimisation's own last "
            f"energy {o['E_final_Ha']:.10f} Ha; includes the D3(BJ) term")
        add(P + 'e_dispersion', fin['E_dispersion_Ha'], 'Ha', cs, f'{src}:items.final.E_dispersion_Ha',
            'D3(BJ) dispersion energy at the minimum (pyscf-dispersion)')
        add(P + 'e_start', hist[0]['E_Ha'], 'Ha', cs, f'{src}:items.opt.history[0].E_Ha',
            'B3LYP-D3(BJ) energy of the published B3LYP structure (step 0 = molecules/btpc/results/btpc_'
            f"{m}/opt.xyz, checked to 1e-6 A)")
        close(float(np.abs(traj[0][1] - G[m]).max()), 0.0, 1e-6, f'{job} step 0 = published structure')
        add(P + 'opt_energy_change', (fin['E_Ha'] - hist[0]['E_Ha']) * HA_KCAL, 'kcal mol-1', cs,
            f'{src}:items.final.E_Ha, items.opt.history[0].E_Ha',
            'derived: E(D3(BJ) minimum) - E(published B3LYP structure), both with B3LYP-D3(BJ), x 627.509474')
        heavy = [i for i, l in enumerate(ctx['LB']) if not l.startswith('H')]
        add(P + 'rmsd_vs_b3lyp', kabsch_rmsd(G[k], G[m]), 'Å', cs, f"{xyz}; {rel(STORED[m] / 'opt.xyz')}",
            f"derived: best-fit RMSD over all 33 atoms against the published B3LYP {model} structure (result.json "
            f"rmsd_vs_stored_B3LYP_min_A {o['rmsd_vs_stored_B3LYP_min_A']:.4f} A)")
        add(P + 'rmsd_heavy_vs_b3lyp', kabsch_rmsd(G[k][heavy], G[m][heavy]), 'Å', cs,
            f"{xyz}; {rel(STORED[m] / 'opt.xyz')}", f'derived: best-fit RMSD over the {len(heavy)} heavy atoms')
        # geometry
        nb = f'derived from {xyz}'
        add(P + 'proton_on', g['proton_on'], '', cs, xyz, nb + '; atom nearer to H4N (N4 = ion pair, O1 = neutral pair)')
        add(P + 'hb_N4_H4N_O1_dh', g['dh'], 'Å', cs, xyz, nb)
        add(P + 'hb_N4_H4N_O1_ha', g['ha'], 'Å', cs, xyz, nb)
        add(P + 'hb_N4_H4N_O1_da', g['da'], 'Å', cs, xyz, nb)
        add(P + 'hb_N4_H4N_O1_angle', g['angle'], 'deg', cs, xyz, nb)
        add(P + 'proton_asymmetry', g['asym'], 'Å', cs, xyz, nb + '; d(H4N...O1) - d(N4-H4N); 0 = centred proton')
        add(P + 'tors_C1_O1_N4_C12', g['tau'], 'deg', cs, xyz,
            nb + '; pseudo-torsion through the non-bonded O1...N4 contact: the turn of the two ions about N4-H4N...O1')
        add(P + 'interplanar_picrate_cation', g['interplanar'], 'deg', cs, xyz,
            nb + '; least-squares planes of the six picrate ring C and the nine cation ring atoms')
        add(P + 'centroid_distance_picrate_cation', g['centroid'], 'Å', cs, xyz, nb)
        for nm in NITRO:
            add(P + f'nitro_twist_{nm}', g[f'twist_{nm}'], 'deg', cs, xyz, nb + f'; O-{nm}-O plane against the picrate ring plane')
        st = ctx['xref'].stats(G[k])
        add(P + 'bond_rmsd_heavy', st['bond_rmsd_heavy'], 'Å', cs, f'{xyz}; {rel(XRAY_BTPC)}',
            f"derived: {st['n_bonds_heavy']} bonds between non-H atoms, DFT - X-ray (as gas.bond_rmsd_heavy)")
        add(P + 'angle_rmsd_heavy', st['angle_rmsd_heavy'], 'deg', cs, f'{xyz}; {rel(XRAY_BTPC)}',
            f"derived: {st['n_angles_heavy']} angles with no H atom (as gas.angle_rmsd_heavy)")
        add(P + 'kabsch_rmsd_heavy', st['kabsch_rmsd_heavy'], 'Å', cs, f'{xyz}; {rel(XRAY_BTPC)}',
            f"derived: best-fit superposition RMSD of the {st['n_heavy']} heavy atoms against the X-ray pair")
        # differences against the published B3LYP structure of the same model and against X-ray
        for q, fk, u in (('da', 'hb_N4_H4N_O1_da', 'Å'), ('ha', 'hb_N4_H4N_O1_ha', 'Å'), ('dh', 'hb_N4_H4N_O1_dh', 'Å'),
                         ('angle', 'hb_N4_H4N_O1_angle', 'deg'), ('interplanar', 'interplanar_picrate_cation', 'deg')):
            add(P + f'diff_vs_b3lyp_{fk}', g[q] - pub(facts, f'{m}.{fk}'), u, cs, f'{xyz}; facts {m}.{fk}',
                f'derived: rev.struct.d3_{m}.{fk} - {m}.{fk} (B3LYP-D3(BJ) minus B3LYP, same model)')
        add(P + 'diff_vs_b3lyp_tors_C1_O1_N4_C12', g['tau'] - GEO[m]['tau'], 'deg', cs,
            f"{xyz}; {rel(STORED[m] / 'opt.xyz')}", 'derived: D3(BJ) minus published B3LYP pseudo-torsion')
        add(P + 'diff_vs_xray_hb_N4_H4N_O1_da', g['da'] - GEO['xray']['da'], 'Å', cs, f'{xyz}; {rel(XRAY_BTPC)}',
            'derived: N4...O1 of the D3(BJ) minimum minus X-ray')
        # frontier orbitals and charges at the minimum
        fr = fin['frontier']
        orb = fr['orbitals']
        fsrc = f'{src}:items.final.frontier'
        add(P + 'homo_ev', fr['HOMO_eV'], 'eV', cs, fsrc + '.HOMO_eV', 'final single point at the D3(BJ) minimum')
        add(P + 'lumo_ev', fr['LUMO_eV'], 'eV', cs, fsrc + '.LUMO_eV', 'final single point at the D3(BJ) minimum')
        add(P + 'gap_ev', fr['gap_eV'], 'eV', cs, fsrc + '.gap_eV', 'E(LUMO) - E(HOMO)')
        for o_, key in (('HOMO', 'homo'), ('LUMO', 'lumo'), ('LUMO+1', 'lumo_p1')):
            add(P + f'{key}_picrate_pct', orb[o_]['mulliken']['anion'], '%', cs, fsrc + f'.orbitals.{o_}.mulliken.anion',
                f"Mulliken share of the picrate; cation {orb[o_]['mulliken']['cation']} %; Lowdin picrate "
                f"{orb[o_]['lowdin']['anion']} %")
        add(P + 'lumo_on', 'picrate' if fr['LUMO_on'] == 'anion' else fr['LUMO_on'], '', cs, fsrc + '.LUMO_on',
            'ion carrying the LUMO (Mulliken share > 50 %)')
        cv = fr['lowest_virtual_on']['cation']
        add(P + 'cation_lowest_virtual', cv['label'], '', cs, fsrc + '.lowest_virtual_on.cation',
            f"lowest virtual orbital with a cation share > 50 % (Mulliken {cv['share_pct']} %)")
        add(P + 'cation_lowest_virtual_ev', cv['eV'], 'eV', cs, fsrc + '.lowest_virtual_on.cation.eV', '')
        add(P + 'ct_gap_picrate_to_cation_ev', cv['eV'] - fr['HOMO_eV'], 'eV', cs, fsrc,
            'derived: E(lowest cation virtual) - E(HOMO), HOMO on the picrate (as gas.ct_gap_picrate_to_cation_ev)')
        for q, fk in (('gap_ev', 'gap_ev'), ('homo_ev', 'homo_ev'), ('lumo_ev', 'lumo_ev')):
            v = {'gap_ev': fr['gap_eV'], 'homo_ev': fr['HOMO_eV'], 'lumo_ev': fr['LUMO_eV']}[q]
            add(P + f'diff_vs_b3lyp_{fk}', v - pub(facts, f'{m}.{fk}'), 'eV', cs, f'{fsrc}; facts {m}.{fk}',
                f'derived: rev.struct.d3_{m}.{fk} - {m}.{fk}')
        add(P + 'dipole_debye', fin['dipole_debye']['total'], 'D', cs, f'{src}:items.final.dipole_debye.total', '')
        for sch in ('hirshfeld', 'cm5'):
            add(P + f'{sch}_charge_picrate', fin['charges'][sch]['anion'], 'e', cs, f'{src}:items.final.charges.{sch}.anion',
                f"{sch} charge summed over the 18 picrate atoms; cation {fin['charges'][sch]['cation']:+.4f} e")
    # the two D3(BJ) minima against each other
    d3 = {m: result(f'opt_d3_{m}')[0]['items']['final'] for m in ('gas', 'pcm')}
    add('rev.struct.d3_pcm_minus_gas_gap_ev', d3['pcm']['frontier']['gap_eV'] - d3['gas']['frontier']['gap_eV'], 'eV',
        'rev_pcm', 'opt_d3_pcm, opt_d3_gas result.json:items.final.frontier.gap_eV', 'derived: IEF-PCM minus gas phase, B3LYP-D3(BJ)')


# ---------------------------------------------------------------------------------------------- B. soft-mode checks
def _softmode_assess(ref_rmsd, ref_max, proj_max):
    """The check tests the curvature along the softest modes of a reference Hessian. It speaks for the target
    structure only if the reference was computed near it and those modes are still soft there."""
    if ref_rmsd <= 0.1 and proj_max <= 1.5 * ref_max:
        return 'informative'
    return 'not informative'


def _softmode(add, ctx):
    G, LB = ctx['G'], ctx['LB']
    rows = ctx['rows'].setdefault('softmode', [])
    specs = [('d3_gas', 'softmode_d3_gas', None, 'rev_gas', G['gas'], G['d3gas'], 'opt_d3_gas'),
             ('d3_pcm', 'softmode_d3_pcm', None, 'rev_pcm', G['pcm'], G['d3pcm'], 'opt_d3_pcm')]
    zg = read_xyz(BTPCZ['gas'] / 'opt.xyz')[1]
    zp = read_xyz(BTPCZ['pcm'] / 'opt.xyz')[1]
    specs.append(('btpcz_pcm', 'c0_btpcz_pcm_softmode', 'softmode/result.json', 'btpcz_pcm', zg, zp, None))
    for tag, job, sub, cs, Xref, Xtgt, optjob in specs:
        jd = job_dir(job)
        p = jd / (sub or 'result.json')
        r = jload(p)
        src = rel(p)
        s, ref = r['summary'], r['reference']
        P = f'rev.struct.softmode_{tag}.'
        if optjob is not None:
            tgt = job_dir(optjob) / 'opt.xyz'
            if r['inputs']['geometry']['sha256'] != sha256(tgt):
                raise AssertionError(f'{job}: target geometry is not {rel(tgt)}')
        else:
            if r['inputs']['geometry']['sha256'] != sha256(BTPCZ['pcm'] / 'opt.xyz'):
                raise AssertionError(f'{job}: target geometry is not the btpcz_pcm opt.xyz')
        nm = int(r['settings']['nmodes'])
        ref_f = ref['lowest_freqs_at_target_geometry_cm-1'][:nm]
        proj = s['projected_freqs_cm-1']
        ref_rmsd = kabsch_rmsd(Xtgt, Xref)
        verdict = _softmode_assess(ref_rmsd, max(ref_f), max(proj))
        refname = {'d3_gas': 'the published B3LYP gas-phase Hessian (molecules/btpc/results/btpc_gas, data/r2 hessian.npy)',
                   'd3_pcm': 'the published B3LYP IEF-PCM Hessian (data/r2/molecules/btpc/results/btpc_pcm/hessian.npy)',
                   'btpcz_pcm': 'the btpcz gas-phase finite-difference Hessian (molecules/btpcz/results/btpcz_gas/hessian.npy)'}[tag]
        add(P + 'reference', refname, '', cs, f'{src}:inputs.ref_hessian', f"sha256 {r['inputs']['ref_hessian']['sha256'][:12]}")
        add(P + 'n_modes', nm, '', cs, f'{src}:settings.nmodes', 'softest vibrational modes of the reference Hessian, taken at the target geometry')
        add(P + 'step_bohr', r['settings']['step_max_atom_bohr'], 'bohr', cs, f'{src}:settings.step_max_atom_bohr',
            f"largest atomic displacement of each +/- step; {r['settings']['differences']}; grid response "
            f"{r['settings']['grid_response']}")
        ng = sum(1 for kk in r.get('timings', {}) if kk == 'x0' or kk.startswith('mode'))
        if ng:
            add(P + 'n_gradients', ng, '', cs, f'{src}:timings', '2 x n_modes displaced gradients + 1 at the structure itself')
        add(P + 'ref_geometry_rmsd', ref_rmsd, 'Å', cs, f'{src}; target and reference opt.xyz',
            'derived: best-fit RMSD over all atoms between the structure checked and the structure of the reference Hessian'
            + (f" (the job's own unaligned value: {ref['rms_geometry_difference_A']:.3g} A)" if tag != 'btpcz_pcm' else
               ' (the job could not compare: it had no reference geometry and stored 0)'))
        add(P + 'ref_freqs', fl(ref_f), 'cm-1', cs, f'{src}:reference.lowest_freqs_at_target_geometry_cm-1',
            'lowest wavenumbers of the reference Hessian, TR projected at the target geometry (negative = imaginary)')
        add(P + 'projected_freqs', fl(proj), 'cm-1', cs, f'{src}:summary.projected_freqs_cm-1',
            'eigenvalues of the Hessian projected on the soft subspace (Rayleigh-Ritz values), from central differences of analytic gradients')
        add(P + 'diagonal_freqs', fl(s['diagonal_freqs_cm-1']), 'cm-1', cs, f'{src}:summary.diagonal_freqs_cm-1',
            'curvature along each reference mode separately, as a wavenumber')
        add(P + 'lowest_projected', min(proj), 'cm-1', cs, f'{src}:summary.projected_freqs_cm-1', 'derived: lowest value')
        add(P + 'highest_projected', max(proj), 'cm-1', cs, f'{src}:summary.projected_freqs_cm-1', 'derived: highest value')
        add(P + 'n_imag_projected', s['n_imaginary_projected'], '', cs, f'{src}:summary.n_imaginary_projected', '')
        add(P + 'updated_lowest', s['updated_full_lowest_freqs_cm-1'][0], 'cm-1', cs,
            f'{src}:summary.updated_full_lowest_freqs_cm-1[0]',
            'lowest wavenumber of the reference Hessian updated with the new columns (estimates the coupling to the rest of the spectrum)')
        add(P + 'n_imag_updated', s['n_imaginary_updated_full'], '', cs, f'{src}:summary.n_imaginary_updated_full', '')
        add(P + 'asymmetry', s['asymmetry_rel'], '', cs, f'{src}:summary.asymmetry_rel',
            '||K - K^T|| / ||K|| of the projected Hessian: numerical noise of the differences')
        add(P + 'asymmetry_pct', 100 * s['asymmetry_rel'], '%', cs, f'{src}:summary.asymmetry_rel', 'derived: x 100')
        rr = s['residual_relative']
        add(P + 'residual_relative_range', fl([min(rr), max(rr)]), '', cs, f'{src}:summary.residual_relative',
            'min and max of ||(1 - QQ^T) H q_k|| / ||H q_k||: share of each column outside the soft subspace')
        add(P + 'residual_grad_rms', r['residual_gradient']['rms_atom_norm'], 'Ha/Bohr', cs,
            f'{src}:residual_gradient.rms_atom_norm', 'gradient at the structure itself (RMS of per-atom norms)')
        add(P + 'verdict_script', s['verdict'], '', cs, f'{src}:summary.verdict', 'verdict printed by softmode_check.py')
        add(P + 'assessment', verdict, '', cs, f'{src}; derived',
            'derived rule: informative only if the reference Hessian belongs to a structure within 0.1 A (best-fit RMSD) of '
            'the one checked and the highest projected wavenumber stays below 1.5 x the highest reference wavenumber, '
            'i.e. the modes probed are still the soft modes of the structure checked')
        for i in range(nm):
            rows.append({'check': tag, 'job': job, 'mode': i + 1, 'ref_cm-1': ref_f[i], 'diagonal_cm-1': s['diagonal_freqs_cm-1'][i],
                         'projected_cm-1': proj[i], 'residual_relative': rr[i], 'assessment': verdict})


# ---------------------------------------------------------------------------------------------- C. torsion scans
SCAN_JOBS = {'b3lyp': ['scan_tors_b3lyp'], 'd3': ['scan_tors_d3', 'scan_tors_d3_up', 'scan_tors_d3_up60']}


def _scan_points(method, LB, IB):
    """Points of one method: list of dicts sorted by the target torsion; reference = unconstrained minimum."""
    pts, ref = [], None
    for job in SCAN_JOBS[method]:
        r, src = result(job)
        jd = job_dir(job)
        want = 'none' if method == 'b3lyp' else 'd3bj'
        if (r['settings']['disp'] or 'none') != want or r['settings']['model'] != 'gas':
            raise AssertionError(f'{job}: settings {r["settings"]["disp"]}/{r["settings"]["model"]}')
        st = r['items']['start']
        if job != 'scan_tors_d3_up60':
            e = (st['E_Ha'], st['torsion_deg'], src)
            if ref is None:
                ref = e
            close(e[0], ref[0], 1e-8, f'{job}: start energy = minimum of {SCAN_JOBS[method][0]}')
        for key, v in r['items'].items():
            if key == 'start' or not isinstance(v, dict) or 'target_deg' not in v:
                continue
            p = jd / v['xyz']
            s, X = read_xyz(p)
            check_symbols(s, LB, rel(p))
            g = pair_geometry(X, IB, 'N4', 'H4N')
            for mine, jk in (('da', 'd_N...O_A'), ('dh', 'd_H-N_A'), ('ha', 'd_H...O_A'), ('tau', 'torsion_deg'),
                             ('interplanar', 'inter_ring_angle_deg')):
                close(g[mine], v['geometry'][jk], TOL[mine], f'{job} {key} {jk}')
            close(tau_norm(g['tau']), v['target_deg'], 0.01, f'{job} {key}: constrained torsion reached')
            pts.append({'method': method, 'job': job, 'src': src, 'xyz': rel(p), 'target': v['target_deg'],
                        'tau': tau_norm(g['tau']), 'E_Ha': v['E_Ha'], 'converged': bool(v['converged']),
                        'steps': int(v['steps']), 'start_from': v['start_from'].split('/')[-1],
                        'proton_transferred': bool(v['geometry']['proton_transferred']), 'g': g,
                        'nitro': [dih(X, *(IB[a] for a in t)) for t in SCAN_NITRO]})
    if method == 'd3':
        # the -60 point started from the converged -90 point; its job recomputed that point's energy
        up = result('scan_tors_d3_up')[0]['items']['p-90.0']['E_Ha']
        st60 = result('scan_tors_d3_up60')[0]['items']['start']['E_Ha']
        close(st60, up, 1e-7, 'scan_tors_d3_up60 start = scan_tors_d3_up p-90')
    pts.sort(key=lambda d: d['target'])
    grid = [d['target'] for d in pts]
    if len(set(grid)) != len(grid):
        raise AssertionError(f'{method}: duplicate scan points {grid}')
    return pts, ref


def _scans(add, ctx):
    facts, LB, IB, GEO = ctx['facts'], ctx['LB'], ctx['IB'], ctx['GEO']
    rt = R_KCAL * T_ROOM
    add('rev.struct.rt_298_kcal', rt, 'kcal mol-1', None, 'derived', f'RT at {T_ROOM} K, R = {R_KCAL} kcal mol-1 K-1 (CODATA 2018)')
    tau_x = GEO['xray']['tau']
    prof = {}
    rows = ctx['rows'].setdefault('torscan', [])
    for method in ('b3lyp', 'd3'):
        pts, ref = _scan_points(method, LB, IB)
        P, cs = f'rev.struct.scan_{method}.', 'rev_gas'
        lvl = 'B3LYP' if method == 'b3lyp' else 'B3LYP-D3(BJ)'
        jobs = ', '.join(SCAN_JOBS[method])
        srcs = '; '.join(sorted({d['src'] for d in pts}))
        if method == 'b3lyp':
            close(ref[0], pub(facts, 'gas.e_total'), 1e-8, 'B3LYP scan start = published gas minimum')
            mintag = 'published B3LYP gas-phase minimum (gas.e_total)'
        else:
            close(ref[0], result('opt_d3_gas')[0]['items']['final']['E_Ha'], 1e-8, 'D3 scan start = opt_d3_gas minimum')
            mintag = 'B3LYP-D3(BJ) gas-phase minimum (rev.struct.d3_gas.e_total)'
        dE = [(d['E_Ha'] - ref[0]) * HA_KCAL for d in pts]
        prof[method] = (pts, dE, ref)
        add(P + 'jobs', jobs, '', cs, srcs, f'relaxed {lvl} gas-phase scan of C1-O1...N4-C12; constraint by geomeTRIC '
            f"$set dihedral, all other coordinates optimised to the GAU criteria at each point; chains run outward from "
            f"the minimum, each point started from its neighbour")
        add(P + 'e_ref', ref[0], 'Ha', cs, ref[2] + ':items.start.E_Ha', f'reference energy: {mintag}, recomputed by the scan job')
        add(P + 'tau_ref', ref[1], 'deg', cs, ref[2] + ':items.start.torsion_deg', 'C1-O1...N4-C12 at that unconstrained minimum')
        add(P + 'grid', fl([d['target'] for d in pts]), 'deg', cs, srcs, 'constrained values of C1-O1...N4-C12, ascending')
        add(P + 'tau', fl([d['tau'] for d in pts]), 'deg', cs, srcs + ' (points/*.xyz)',
            'C1-O1...N4-C12 recomputed from each point (the -180 point is reported as about +180 by atan2 and shifted by -360)')
        add(P + 'dE', fl(dE), 'kcal mol-1', cs, srcs, f'derived: (E(point) - E(unconstrained minimum)) x {HA_KCAL}, same grid order')
        for q, nm, u in (('da', 'd_N4_O1', 'Å'), ('dh', 'd_N4_H4N', 'Å'), ('ha', 'd_H4N_O1', 'Å'),
                         ('angle', 'a_N4_H4N_O1', 'deg'), ('interplanar', 'interplanar', 'deg')):
            add(P + nm, fl([d['g'][q] for d in pts]), u, cs, srcs + ' (points/*.xyz)', 'derived from each scan point, grid order')
        for t, nm in zip(SCAN_NITRO, ('tors_O2_N1_C2_C1', 'tors_O4_N2_C4_C3', 'tors_O7_N3_C6_C1')):
            j = SCAN_NITRO.index(t)
            add(P + nm, fl([d['nitro'][j] for d in pts]), 'deg', cs, srcs + ' (points/*.xyz)',
                'derived: nitro torsion at each scan point, grid order')
        add(P + 'steps', [d['steps'] for d in pts], '', cs, srcs, 'optimisation evaluations per point, grid order')
        add(P + 'n_points', len(pts), '', cs, srcs, '')
        add(P + 'n_converged', sum(d['converged'] for d in pts), '', cs, srcs, 'points that met all five criteria')
        on_n = sum(1 for d in pts if d['g']['proton_on'] == 'N4' and not d['proton_transferred'])
        add(P + 'n_proton_on_N4', on_n, '', cs, srcs, 'points with H4N nearer to N4 than to O1 (and proton_transferred false in the job)')
        lo = int(np.argmin(dE))
        hi = int(np.argmax(dE))
        add(P + 'dE_lowest_point', dE[lo], 'kcal mol-1', cs, srcs, f"derived: lowest scan point, at {pts[lo]['target']:g} deg")
        add(P + 'tau_lowest_point', pts[lo]['target'], 'deg', cs, srcs, 'derived')
        add(P + 'dE_highest_point', dE[hi], 'kcal mol-1', cs, srcs, f"derived: highest scan point, at {pts[hi]['target']:g} deg")
        add(P + 'tau_highest_point', pts[hi]['target'], 'deg', cs, srcs, 'derived')
        add(P + 'span', max(dE) - min(min(dE), 0.0), 'kcal mol-1', cs, srcs,
            f"derived: highest minus lowest energy over {pts[0]['target']:g} to {pts[-1]['target']:g} deg, unconstrained minimum included")
        within = [d['target'] for d, e in zip(pts, dE) if e < rt]
        add(P + 'n_within_rt', len(within), '', cs, srcs, f'derived: scan points less than RT ({rt:.3f} kcal mol-1) above the minimum')
        add(P + 'within_rt_range', fl([min(within), max(within)]), 'deg', cs, srcs, 'derived: first and last such grid value')
        steps_ok = all(abs((b_ - a_) - 15.0) < 1e-9 for a_, b_ in zip(within, within[1:]))
        add(P + 'contiguous_within_rt', 'yes' if steps_ok else 'no', '', cs, srcs, 'derived: the points below RT form one unbroken run of the 15 deg grid')
        i165 = [d['target'] for d in pts].index(-165.0)
        add(P + 'dE_m165', dE[i165], 'kcal mol-1', cs, pts[i165]['src'], f'the point nearest the X-ray torsion ({tau_x:.1f} deg)')
        add(P + 'interplanar_m165', pts[i165]['g']['interplanar'], 'deg', cs, pts[i165]['xyz'], 'derived')
        add(P + 'd_N4_O1_m165', pts[i165]['g']['da'], 'Å', cs, pts[i165]['xyz'], 'derived')
        i150 = [d['target'] for d in pts].index(-150.0)
        w = (tau_x - (-165.0)) / 15.0
        add(P + 'dE_at_xray_tau', dE[i165] + w * (dE[i150] - dE[i165]), 'kcal mol-1', cs, srcs,
            f'derived: linear interpolation between the -165 and -150 deg points at the X-ray value {tau_x:.2f} deg')
        for q, nm, u in (('da', 'd_N4_O1', 'Å'), ('dh', 'd_N4_H4N', 'Å'), ('ha', 'd_H4N_O1', 'Å'), ('interplanar', 'interplanar', 'deg')):
            v = [d['g'][q] for d in pts]
            add(P + nm + '_range', fl([min(v), max(v)]), u, cs, srcs, 'derived: min and max over the scan points')
        for j, nm in enumerate(('tors_O2_N1_C2_C1', 'tors_O4_N2_C4_C3', 'tors_O7_N3_C6_C1')):
            v = [d['nitro'][j] for d in pts]
            add(P + nm + '_span', max(v) - min(v), 'deg', cs, srcs, 'derived: max - min over the scan points')
        for d, e in zip(pts, dE):
            rows.append({'method': lvl, 'job': d['job'], 'target_deg': d['target'], 'tau_deg': d['tau'], 'E_Ha': d['E_Ha'],
                         'dE_kcal': e, 'd_N4_O1_A': d['g']['da'], 'd_N4_H4N_A': d['g']['dh'], 'd_H4N_O1_A': d['g']['ha'],
                         'a_N4_H4N_O1_deg': d['g']['angle'], 'interplanar_deg': d['g']['interplanar'],
                         'tors_O2_N1_C2_C1': d['nitro'][0], 'tors_O4_N2_C4_C3': d['nitro'][1], 'tors_O7_N3_C6_C1': d['nitro'][2],
                         'converged': d['converged'], 'steps': d['steps'], 'proton_on': d['g']['proton_on'],
                         'start_from': d['start_from']})
    b, d3 = prof['b3lyp'], prof['d3']
    if [p['target'] for p in b[0]] != [p['target'] for p in d3[0]]:
        raise AssertionError('B3LYP and D3(BJ) scans have different grids')
    diff = [e2 - e1 for e1, e2 in zip(b[1], d3[1])]
    add('rev.struct.scan_d3_minus_b3lyp_dE', fl(diff), 'kcal mol-1', 'rev_gas', 'rev.struct.scan_d3.dE; rev.struct.scan_b3lyp.dE',
        'derived: D3(BJ) relative energy minus B3LYP relative energy at each grid point (each relative to its own minimum)')
    add('rev.struct.scan_d3_minus_b3lyp_dE_m180', diff[0], 'kcal mol-1', 'rev_gas', 'rev.struct.scan_*.dE[0]', 'derived: at -180 deg')
    add('rev.struct.scan_d3_minus_b3lyp_dE_m60', diff[-1], 'kcal mol-1', 'rev_gas', 'rev.struct.scan_*.dE[-1]', 'derived: at -60 deg')


# ---------------------------------------------------------------------------------------------- D. torsion angles
def tname(t):
    nonb = {('O1', 'N4'), ('N4', 'O1'), ('O1', 'H4N'), ('H4N', 'O1')}
    s = t[0]
    for a, b in zip(t, t[1:]):
        s += ('...' if (a, b) in nonb else '-') + b
    return s


def _torsions(add, ctx):
    G, IB, GEO = ctx['G'], ctx['IB'], ctx['GEO']
    rows = ctx['rows'].setdefault('torsions', [])
    srcs = {'xray': rel(XRAY_BTPC), 'gas': rel(STORED['gas'] / 'opt.xyz'), 'pcm': rel(STORED['pcm'] / 'opt.xyz'),
            'd3gas': rel(job_dir('opt_d3_gas') / 'opt.xyz'), 'd3pcm': rel(job_dir('opt_d3_pcm') / 'opt.xyz')}
    pic = [IB[a] for a in PIC_RING]
    cat = [IB[a] for a in CAT9]
    for t in TORSIONS:
        row = {'torsion': tname(t)}
        for mdl, X in G.items():
            v = dih(X, *(IB[a] for a in t))
            row[mdl] = v
            note = f'derived from {srcs[mdl]} ({MODEL_TAG[mdl]})'
            if 'O1' in t and ('N4' in t or 'H4N' in t) and t[0] in ('C1', 'C6'):
                note += '; pseudo-torsion through a non-bonded contact'
            if mdl == 'xray':
                note += '; the CIF has no torsion loop, so no standard uncertainty'
            add(f"rev.struct.tors_{mdl}_{'_'.join(t)}", v, 'deg', MODEL_CASE[mdl], srcs[mdl], note)
        rows.append(row)
    for mdl, X in G.items():
        close(dih(X, IB['C1'], IB['O1'], IB['N4'], IB['C12']), GEO[mdl]['tau'], 1e-9, 'tau definition')
        add(f'rev.struct.line_O1_N4_vs_picrate_{mdl}', elevation(X, IB['O1'], IB['N4'], pic), 'deg', MODEL_CASE[mdl], srcs[mdl],
            'derived: angle between the O1->N4 line and the least-squares plane of the six picrate ring carbons')
        add(f'rev.struct.line_O1_N4_vs_cation_{mdl}', elevation(X, IB['N4'], IB['O1'], cat), 'deg', MODEL_CASE[mdl], srcs[mdl],
            'derived: angle between the N4->O1 line and the least-squares plane of the nine cation ring atoms')
    for nm in ('line_O1_N4_vs_picrate', 'line_O1_N4_vs_cation'):
        rows.append({'torsion': nm, **{mdl: (elevation(X, IB['O1'], IB['N4'], pic) if 'picrate' in nm
                                            else elevation(X, IB['N4'], IB['O1'], cat)) for mdl, X in G.items()}})
    # X-ray ring orientation for the torsion table (published as xray.interplanar_picrate_cation etc.)
    # C1-O1...N4-C12 along the published IEF-PCM optimisation (Table S8 steps)
    traj = read_traj(RAW_PCM_TRAJ)
    if float(np.abs(traj[-1][1] - G['pcm']).max()) > 1e-6:
        raise AssertionError('published IEF-PCM trajectory does not end at opt.xyz')
    tau = [dih(X, IB['C1'], IB['O1'], IB['N4'], IB['C12']) for _, X in traj]
    add('rev.struct.pcm_traj_tors_C1_O1_N4_C12', fl(tau, 3), 'deg', 'btpc_pcm', rel(RAW_PCM_TRAJ),
        f'derived: C1-O1...N4-C12 at each of the {len(traj)} frames (steps 0-{len(traj) - 1}) of the published B3LYP IEF-PCM optimisation')
    add('rev.struct.pcm_traj_tors_step0', tau[0], 'deg', 'btpc_pcm', rel(RAW_PCM_TRAJ), 'derived: step 0')
    add('rev.struct.pcm_traj_tors_change', tau[-1] - tau[0], 'deg', 'btpc_pcm', rel(RAW_PCM_TRAJ), 'derived: last minus first frame')


# ---------------------------------------------------------------------------------------------- E. P2_1/n pair
def _btpcz(add, ctx):
    facts = ctx['facts']
    cz, LZ, IZ = labels_of(CASE_BTPCZ)
    czg = jload(ROOT / 'molecules/btpcz/cases/btpcz_gas.json')
    if [l.split('#')[0] for l in czg['cif_labels']] != LZ:
        raise AssertionError('btpcz gas and pcm case files list different atoms')
    PIC = sorted(cz['fragments']['picrate anion'])
    CAT = sorted(cz['fragments']['benzotriazolium cation'])
    pr = cz['proton']
    donor, Hh, acc = LZ[pr['donor_N']], LZ[pr['H']], LZ[pr['acceptor_O']]
    sx, XZ = read_xyz(XRAY_BTPCZ)
    check_symbols(sx, LZ, rel(XRAY_BTPCZ))
    xref = XrayRef(sx, XZ)
    # the fused carbon next to the donor nitrogen, from the X-ray connectivity (btpc: N4-C12)
    fused = [LZ[j] for i, j in xref.bonds + [(b, a) for a, b in xref.bonds] if LZ[i] == donor and LZ[j].startswith('C')]
    if len(fused) != 1:
        raise AssertionError(f'btpcz: donor {donor} has carbon neighbours {fused}')
    fused = fused[0]
    cpath = rel(CASE_BTPCZ)
    sg = re.search(r'space group (\S+?),', cz['notes'])
    add('rev.struct.btpcz.source', f"COD {cz['cod']}, doi:{cz['source_doi']}", '', 'btpcz_xray', cpath + ':cod, source_doi',
        'second crystal form of BTPC (Zeng, Li and Wang, Acta Cryst. E67, o1464, 2011); molecules/btpcz/input/README.md')
    add('rev.struct.btpcz.space_group', sg.group(1) if sg else '?', '', 'btpcz_xray', cpath + ':notes', '')
    for k, nm, u in zip(range(6), ('cell_a', 'cell_b', 'cell_c', 'cell_alpha', 'cell_beta', 'cell_gamma'),
                        ('Å', 'Å', 'Å', 'deg', 'deg', 'deg')):
        add(f'rev.struct.btpcz.{nm}', cz['cell'][k], u, 'btpcz_xray', cpath + f':cell[{k}]', 'from the COD CIF')
    add('rev.struct.btpcz.defining_hbond', f'{donor}-{Hh}...{acc}', '', 'btpcz_xray', cpath + ':proton, hbonds',
        f'the 33-atom asymmetric unit is itself the hydrogen-bonded pair; donor {donor} is bonded to {fused}, so '
        f'{donor}/{Hh}/{fused} of this form correspond to N4/H4N/C12 of the P2_1/c pair; picrate labels coincide')
    add('rev.struct.btpcz.label_map', f'{donor}=N4, {Hh}=H4N, {fused}=C12', '', 'btpcz_xray', cpath,
        'derived from the X-ray connectivity: donor N and the fused ring carbon bonded to it')
    geo = {'xray': pair_geometry(XZ, IZ, donor, Hh, acc, fused)}
    Xz = {'xray': XZ}
    cores = {'gas': 'c0_btpcz_gas_core', 'pcm': 'c0_btpcz_pcm_core'}
    rows = ctx['rows'].setdefault('minima', [])
    for m in ('xray', 'gas', 'pcm'):
        cs = f'btpcz_{m}'
        P = f'rev.struct.btpcz_{m}.'
        if m == 'xray':
            src = rel(XRAY_BTPCZ)
            X = XZ
        else:
            p = BTPCZ[m] / 'opt.xyz'
            s, X = read_xyz(p)
            check_symbols(s, LZ, rel(p))
            src = rel(p)
            Xz[m] = X
            geo[m] = pair_geometry(X, IZ, donor, Hh, acc, fused)
        g = geo[m]
        nb = f'derived from {src}'
        add(P + 'proton_on', g['proton_on'], '', cs, src, nb + f'; atom nearer to {Hh} ({donor} = ion pair, {acc} = neutral pair)')
        add(P + f'hb_{donor}_{Hh}_{acc}_dh', g['dh'], 'Å', cs, src, nb + f'; d({donor}-{Hh})')
        add(P + f'hb_{donor}_{Hh}_{acc}_ha', g['ha'], 'Å', cs, src, nb + f'; d({Hh}...{acc})')
        add(P + f'hb_{donor}_{Hh}_{acc}_da', g['da'], 'Å', cs, src, nb)
        add(P + f'hb_{donor}_{Hh}_{acc}_angle', g['angle'], 'deg', cs, src, nb)
        add(P + f'tors_C1_O1_{donor}_{fused}', g['tau'], 'deg', cs, src,
            nb + f'; pseudo-torsion C1-O1...{donor}-{fused}, the counterpart of C1-O1...N4-C12 of the P2_1/c pair; '
            'the mirror image of a structure has the opposite sign, so compare magnitudes')
        add(P + 'interplanar_picrate_cation', g['interplanar'], 'deg', cs, src, nb + '; six picrate ring C, nine cation ring atoms')
        add(P + 'centroid_distance_picrate_cation', g['centroid'], 'Å', cs, src, nb)
        for nm in NITRO:
            add(P + f'nitro_twist_{nm}', g[f'twist_{nm}'], 'deg', cs, src, nb)
        if m == 'xray':
            continue
        # optimisation, energies, frequencies, frontier orbitals
        core, csrc = result(cores[m])
        osum = jload(BTPCZ[m] / 'opt_summary.json')
        ojob = jload(job_dir(cores[m]) / 'opt/result.json')
        osrc = rel(job_dir(cores[m]) / 'opt/result.json')
        if core['items']['optimize']['opt_xyz_sha256'] != sha256(BTPCZ[m] / 'opt.xyz'):
            raise AssertionError(f'{cores[m]}: opt.xyz differs from repo_results')
        for mine, jk in (('da', 'd_N...O_A'), ('dh', 'd_H-N_A'), ('ha', 'd_H...O_A'), ('tau', 'torsion_deg'),
                         ('interplanar', 'inter_ring_angle_deg')):
            close(g[mine], ojob['items']['opt']['geometry'][jk], TOL[mine], f'{cores[m]} {jk}')
        st = ojob['settings']
        lvl = (f"B3LYP/6-311++G(d,p), {'Gas phase' if m == 'gas' else 'IEF-PCM (water)'}, grid {st['grid_level']}"
               + (f", Lebedev {st['lebedev_order']}" if m == 'pcm' else '') + ' (no dispersion correction, as the published pair)')
        add(P + 'level', lvl, '', cs, f'{osrc}:settings', 'pre-optimised at B3LYP/6-31G(d) (GAU_LOOSE) from the crystal pair')
        add(P + 'preopt_steps', core['items']['preopt']['steps'], '', cs, f'{csrc}:items.preopt.steps', core['items']['preopt']['how'])
        add(P + 'opt_converged_on', core['items']['optimize']['converged_on'], '', cs, f'{csrc}:items.optimize.converged_on',
            f"cap {st['maxsteps_total']} evaluations")
        add(P + 'opt_steps', core['items']['optimize']['steps'], '', cs, f'{csrc}:items.optimize.steps', 'evaluations, step 0 included')
        hist = ojob['items']['opt']['history']
        E = [h['E_Ha'] for h in hist]
        tail = 8
        alt = [abs(E[i] - E[i - 1]) for i in range(len(E) - tail + 1, len(E))]
        add(P + 'opt_tail_mean_abs_dE', float(np.mean(alt)) * HA_KCAL, 'kcal mol-1', cs, f'{osrc}:items.opt.history',
            f'derived: mean |E(i) - E(i-1)| over the last {tail - 1} steps, x {HA_KCAL}; a value far above what the step '
            'sizes and gradients allow signals energies alternating between two SCF solutions')
        add(P + 'opt_last_dE', ojob['items']['opt']['last_step']['dE_Ha'], 'Ha', cs, f'{osrc}:items.opt.last_step.dE_Ha',
            'energy change from the previous evaluation')
        add(P + 'opt_dE_two_steps', E[-1] - E[-3], 'Ha', cs, f'{osrc}:items.opt.history',
            'derived: E(last) - E(last but two), the change between evaluations on the same branch of an alternation')
        fj = jload(BTPCZ[m] / 'fmo.json')
        add(P + 'e_total', fj['E_total_Ha'], 'Ha', cs, rel(BTPCZ[m] / 'fmo.json') + ':E_total_Ha',
            f"optimisation's last energy {osum['E_final_Ha']:.10f} Ha (opt_summary.json)")
        rmsd = xref.stats(X)
        add(P + 'bond_rmsd_heavy', rmsd['bond_rmsd_heavy'], 'Å', cs, f'{src}; {rel(XRAY_BTPCZ)}',
            f"derived: {rmsd['n_bonds_heavy']} bonds between non-H atoms, DFT - X-ray of this crystal form")
        add(P + 'angle_rmsd_heavy', rmsd['angle_rmsd_heavy'], 'deg', cs, f'{src}; {rel(XRAY_BTPCZ)}',
            f"derived: {rmsd['n_angles_heavy']} angles with no H atom")
        add(P + 'kabsch_rmsd_heavy', rmsd['kabsch_rmsd_heavy'], 'Å', cs, f'{src}; {rel(XRAY_BTPCZ)}',
            f"derived: best-fit RMSD of the {rmsd['n_heavy']} heavy atoms against the X-ray pair of this crystal form")
        fq = jload(BTPCZ[m] / 'freq.json')
        fsrc = rel(BTPCZ[m] / 'freq.json')
        add(P + 'freq_method', fq['hessian_method'], '', cs, fsrc + ':hessian_method',
            f"job {'c0_btpcz_gas_core' if m == 'gas' else 'c0_btpcz_pcm_fdhess'}")
        add(P + 'n_imag', fq['n_imaginary'], '', cs, fsrc + ':n_imaginary', '')
        add(P + 'freq_lowest', fq['frequencies_cm-1'][0], 'cm-1', cs, fsrc + ':frequencies_cm-1[0]', 'negative = imaginary')
        add(P + 'freq_lowest6', fl(fq['frequencies_cm-1'][:6]), 'cm-1', cs, fsrc + ':frequencies_cm-1[:6]', 'unscaled harmonic')
        add(P + 'zpe', fq['ZPE_Ha'], 'Ha', cs, fsrc + ':ZPE_Ha', 'real modes only')
        add(P + 'g298', fq['G298_Ha'], 'Ha', cs, fsrc + ':G298_Ha', '298.15 K, 1 atm, RRHO')
        # raw finite-difference Hessian: asymmetry before symmetrisation (numerical quality)
        fdjob = 'c0_btpcz_gas_core' if m == 'gas' else 'c0_btpcz_pcm_fdhess'
        fdd = job_dir(fdjob) / 'fd'
        n3 = 3 * len(LZ)
        C = np.zeros((n3, n3))
        seen = set()
        for f in sorted(fdd.glob('cols_*.npz')):
            z = np.load(f)
            for kk in z.files:
                C[:, int(kk)] = z[kk].ravel()
                seen.add(int(kk))
        if len(seen) != n3:
            raise AssertionError(f'{fdjob}: {len(seen)} of {n3} Hessian columns')
        H4 = np.load(BTPCZ[m] / 'hessian.npy')
        close(float(np.abs(H4.transpose(0, 2, 1, 3).reshape(n3, n3) - 0.5 * (C + C.T)).max()), 0.0, 1e-10,
              f'{m}: hessian.npy = symmetrised raw columns')
        add(P + 'fd_hessian_asymmetry', float(np.linalg.norm(C - C.T) / np.linalg.norm(C)), '', cs, rel(fdd) + '/cols_*.npz',
            'derived: ||C - C^T|| / ||C|| of the raw finite-difference Hessian columns before symmetrisation (numerical noise '
            'and any column computed on a different SCF solution)')
        add(P + 'fd_hessian_asymmetry_pct', 100 * float(np.linalg.norm(C - C.T) / np.linalg.norm(C)), '%', cs,
            rel(fdd) + '/cols_*.npz', 'derived: x 100')
        x, mm, f, V = normal_modes(s, X, H4)
        close(f[0], fq['frequencies_cm-1'][0], 1e-6, f'{m}: lowest wavenumber recomputed')
        q1 = V[:, 0]
        v1 = q1.reshape(-1, 3)
        add(P + 'mode1_rigid_fraction', float(np.linalg.norm(rigid_basis(x, mm, [PIC, CAT]).T @ q1) ** 2), '', cs, rel(BTPCZ[m] / 'hessian.npy'),
            'derived: share of the lowest mode in the 12 rigid-body motions of the two ions (as saddle.imag_mode_rigid_fraction)')
        add(P + 'mode1_cation_weight', float((v1[CAT] ** 2).sum()), '', cs, rel(BTPCZ[m] / 'hessian.npy'),
            'derived: share of the mass-weighted amplitude of the lowest mode on the cation')
        n1o = float((v1[[IZ['O2'], IZ['O3']]] ** 2).sum())
        add(P + 'mode1_O2_O3_weight', n1o, '', cs, rel(BTPCZ[m] / 'hessian.npy'),
            'derived: share of the mass-weighted amplitude on O2 and O3 (the ortho-nitro group on C2)')
        tv = nitro_torsion_vector(x, mm, IZ, 'N1', 'C2', ('O2', 'O3'))
        add(P + 'mode1_nitro_N1_torsion_overlap', float(np.dot(q1, tv) ** 2), '', cs, rel(BTPCZ[m] / 'hessian.npy'),
            'derived: squared overlap of the lowest mode with a rigid turn of the O2-N1-O3 group about C2-N1')
        # frontier orbitals (Mulliken fragment shares from fmo.json)
        fc = fj['fragment_composition']
        fmsrc = rel(BTPCZ[m] / 'fmo.json')
        add(P + 'homo_ev', fj['HOMO_eV'], 'eV', cs, fmsrc + ':HOMO_eV', '')
        add(P + 'lumo_ev', fj['LUMO_eV'], 'eV', cs, fmsrc + ':LUMO_eV', '')
        add(P + 'gap_ev', fj['gap_eV'], 'eV', cs, fmsrc + ':gap_eV', '')
        for o_, key in (('HOMO', 'homo'), ('LUMO', 'lumo'), ('HOMO-1', 'homo_m1'), ('LUMO+1', 'lumo_p1')):
            add(P + f'{key}_picrate_pct', fc[o_]['picrate anion'], '%', cs, fmsrc + f':fragment_composition.{o_}',
                f"Mulliken; benzotriazolium/benzotriazole fragment {fc[o_]['benzotriazolium cation']} %"
                + ('; the fragments are neutral here (proton on O1)' if g['proton_on'] == acc else ''))
        add(P + 'dipole_debye', fj['dipole_Debye'], 'D', cs, fmsrc + ':dipole_Debye', '')
        rows.append({'structure': f'P2_1/n {m}', **{k: g[k] for k in ('dh', 'ha', 'da', 'angle', 'tau', 'interplanar')},
                     'proton_on': g['proton_on'], 'n_imag': fq['n_imaginary'], 'freq_lowest': fq['frequencies_cm-1'][0],
                     'E_Ha': fj['E_total_Ha'], 'gap_eV': fj['gap_eV']})
    rows.append({'structure': 'P2_1/n xray', **{k: geo['xray'][k] for k in ('dh', 'ha', 'da', 'angle', 'tau', 'interplanar')},
                 'proton_on': geo['xray']['proton_on']})
    # energies against the P2_1/c pair at the same level and composition
    e = {m: jload(BTPCZ[m] / 'fmo.json')['E_total_Ha'] for m in ('gas', 'pcm')}
    for m in ('gas', 'pcm'):
        cs = f'btpcz_{m}'
        P = f'rev.struct.btpcz_{m}.'
        add(P + 'e_minus_btpc', (e[m] - pub(facts, f'{m}.e_total')) * HA_KCAL, 'kcal mol-1', cs,
            f"{rel(BTPCZ[m] / 'fmo.json')}:E_total_Ha; facts {m}.e_total",
            f'derived: E(P2_1/n pair) - E(published P2_1/c pair), same functional, basis, grid and model, x {HA_KCAL}'
            + ('; the P2_1/n structure is the neutral O1-H6...N6 pair, the P2_1/c one the ion pair' if geo[m]['proton_on'] == acc else ''))
        fq = jload(BTPCZ[m] / 'freq.json')
        add(P + 'e_zpe_minus_btpc', ((e[m] + fq['ZPE_Ha']) - (pub(facts, f'{m}.e_total') + pub(facts, f'{m}.zpe'))) * HA_KCAL,
            'kcal mol-1', cs, f"{rel(BTPCZ[m] / 'freq.json')}:ZPE_Ha; facts {m}.zpe",
            'derived: difference of E + ZPE (harmonic, unscaled)')
        add(P + 'g298_minus_btpc', (fq['G298_Ha'] - pub(facts, f'{m}.g298')) * HA_KCAL, 'kcal mol-1', cs,
            f"{rel(BTPCZ[m] / 'freq.json')}:G298_Ha; facts {m}.g298",
            'derived: difference of G(298.15 K), RRHO; the many modes below 100 cm-1 are treated as harmonic oscillators, '
            'so this difference carries several tenths of a kcal mol-1 or more of model error')
        add(P + 'gap_minus_btpc', jload(BTPCZ[m] / 'fmo.json')['gap_eV'] - pub(facts, f'{m}.gap_ev'), 'eV', cs,
            f"{rel(BTPCZ[m] / 'fmo.json')}:gap_eV; facts {m}.gap_ev", 'derived')
        add(P + 'da_minus_btpc', geo[m]['da'] - pub(facts, f'{m}.hb_N4_H4N_O1_da'), 'Å', cs, f'{rel(BTPCZ[m] / "opt.xyz")}; facts {m}.hb_N4_H4N_O1_da',
            f'derived: {donor}...{acc} of the P2_1/n pair minus N4...O1 of the P2_1/c pair')
        add(P + 'abs_tors_minus_btpc', abs(geo[m]['tau']) - abs(ctx['GEO'][m]['tau']), 'deg', cs, f'{rel(BTPCZ[m] / "opt.xyz")}; {rel(STORED[m] / "opt.xyz")}',
            f'derived: |C1-O1...{donor}-{fused}| of the P2_1/n pair minus |C1-O1...N4-C12| of the published P2_1/c pair')
    add('rev.struct.btpcz_xray.abs_tors_minus_btpc', abs(geo['xray']['tau']) - abs(ctx['GEO']['xray']['tau']), 'deg', 'btpcz_xray',
        f'{rel(XRAY_BTPCZ)}; {rel(XRAY_BTPC)}', 'derived: |pseudo-torsion| of the P2_1/n crystal pair minus that of the P2_1/c crystal pair')
    add('rev.struct.btpcz_xray.da_minus_btpc', geo['xray']['da'] - ctx['GEO']['xray']['da'], 'Å', 'btpcz_xray',
        f'{rel(XRAY_BTPCZ)}; {rel(XRAY_BTPC)}', f'derived: {donor}...{acc} (P2_1/n) minus N4...O1 (P2_1/c), crystal pairs')
    # the gas-phase P2_1/n structure against the neutral pair reached from the P2_1/c pair (proton-topic job, read only)
    try:
        ns, nsrc = result('neutral_start_gas')
        en = ns['summary']['gas']['E_Ha']
        add('rev.struct.btpcz_gas.e_minus_neutral_gas', (e['gas'] - en) * HA_KCAL, 'kcal mol-1', 'btpcz_gas',
            f'{rel(BTPCZ["gas"] / "fmo.json")}:E_total_Ha; {nsrc}:summary.gas.E_Ha',
            f"derived: E(P2_1/n gas structure) - E(neutral pair optimised from the P2_1/c pair, job neutral_start_gas, "
            f"proton on {ns['summary']['gas']['proton_on']}), x {HA_KCAL}")
    except FileNotFoundError:
        pass
    # IEF-PCM optimisation path of the P2_1/n pair: the turn about the hydrogen bond
    for m in ('gas', 'pcm'):
        tj = job_dir(cores[m])
        taus = []
        for sub in ('preopt', 'opt'):
            for _, X in read_traj(tj / sub / 'opt_traj.xyz'):
                taus.append(dih(X, IZ['C1'], IZ['O1'], IZ[donor], IZ[fused]))
        add(f'rev.struct.btpcz_{m}.tors_path_range', fl([min(taus), max(taus)]), 'deg', f'btpcz_{m}',
            f'{rel(tj)}/{{preopt,opt}}/opt_traj.xyz',
            f'derived: smallest and largest C1-O1...{donor}-{fused} met along the pre-optimisation and optimisation '
            f'({len(taus)} frames); start {taus[0]:.1f} deg (crystal), end {taus[-1]:.1f} deg')
    ctx['btpcz'] = {'geo': geo, 'donor': donor, 'H': Hh, 'acc': acc, 'fused': fused}


# ---------------------------------------------------------------------------------------------- CSV output
def _write_csvs(ctx):
    DERIVED.mkdir(parents=True, exist_ok=True)

    def dump(name, rows, header=None):
        if not rows:
            return
        cols = header or list(rows[0].keys())
        with open(DERIVED / name, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore', lineterminator='\n')
            w.writeheader()
            for r in rows:
                w.writerow({k: (f'{v:.10g}' if isinstance(v, float) else v) for k, v in r.items()})

    dump('rev_struct_torscan.csv', ctx['rows']['torscan'])
    dump('rev_struct_torsions.csv', ctx['rows']['torsions'], ['torsion', 'xray', 'gas', 'pcm', 'd3gas', 'd3pcm'])
    dump('rev_struct_softmode.csv', ctx['rows']['softmode'])
    # key geometry of every structure (P2_1/c: X-ray, B3LYP and D3(BJ) minima; P2_1/n: X-ray and B3LYP minima)
    rows = []
    for k in ('xray', 'gas', 'd3gas', 'pcm', 'd3pcm'):
        g = ctx['GEO'][k]
        rows.append({'structure': f'P2_1/c {k}', **{q: g[q] for q in ('dh', 'ha', 'da', 'angle', 'tau', 'interplanar')},
                     'proton_on': g['proton_on']})
    rows += ctx['rows']['minima']
    dump('rev_struct_minima.csv', rows, ['structure', 'proton_on', 'dh', 'ha', 'da', 'angle', 'tau', 'interplanar',
                                         'n_imag', 'freq_lowest', 'E_Ha', 'gap_eV'])
