"""Harmonic normal modes of the BTPC ion pair from the stored finite-difference Hessians (post-processing only).

    modes(model) -> dict(freq, vecs_mw, disp, masses, symbols, xyz, freq_json, max_dev)

Recipe (PLAN.md Fig. 6): hessian.npy (33, 33, 3, 3) in Ha/Bohr^2 -> transpose(0, 2, 1, 3).reshape(99, 99),
symmetrised, mass-weighted with the masses PySCF itself used (H 1.008, C 12.011, N 14.007, O 15.999; facts.md
fact-check log, correction 3: with these the recomputation from opt.xyz reproduces freq.json to ~1e-7 cm-1; the
IUPAC mean masses leave a 0.10 cm-1 residual), the 6 translations and rotations about the centre of mass
projected out (QR-orthonormalised projector), np.linalg.eigh, nu = sign(lambda) sqrt(|lambda|) x
sqrt(E_h / (a0^2 u)) / (2 pi c), and the 6 eigenvalues nearest zero dropped. The result is checked against
freq.json:frequencies_cm-1 (every mode within 0.5 cm-1, as required before drawing).
Cartesian displacement of mode k: eigenvector / sqrt(m), normalised so the largest atomic |d| is 1.
"""
from __future__ import annotations

import json

import numpy as np

import _struct_common as C
from figkit import geometry as geo, style

MASS = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999}             # the masses PySCF used (fact-check log, 3)
E_H = 4.3597447222071e-18       # J
A0 = 5.29177210903e-11          # m
AMU = 1.66053906660e-27         # kg
C_CM = 2.99792458e10            # cm s-1
CONV = np.sqrt(E_H / (A0 ** 2 * AMU)) / (2 * np.pi * C_CM)            # 5140.487 cm-1 per sqrt(Ha/(Bohr^2 u))
TOL_CM = 0.5                    # task requirement; the actual deviation is ~1e-7 cm-1 (printed by the scripts)

_CACHE = {}


def freq_json(model, name='freq.json'):
    p = C.REPO / 'molecules' / 'btpc' / 'results' / C.CASES[model] / name
    d = json.loads(p.read_text())
    f = d['frequencies_cm-1']
    if isinstance(f, str):
        f = json.loads(f)
    return np.array(f, float), d


def modes(model):
    if model in _CACHE:
        return _CACHE[model]
    g = C.geometries()
    sym = g['symbols']
    X = g[model] / geo.BOHR                                                  # Bohr
    H4 = np.load(C.REPO / 'data' / 'r2' / 'molecules' / 'btpc' / 'results' / C.CASES[model] / 'hessian.npy')
    n = len(sym)
    assert H4.shape == (n, n, 3, 3), H4.shape
    Hc = H4.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)
    Hc = 0.5 * (Hc + Hc.T)
    m = np.array([MASS[s] for s in sym])
    m3 = np.repeat(m, 3)
    Hmw = Hc / np.sqrt(np.outer(m3, m3))
    com = (m[:, None] * X).sum(0) / m.sum()
    R = X - com
    T = np.zeros((3 * n, 6))
    sq = np.sqrt(m)
    for k in range(3):
        T[k::3, k] = sq
        e = np.zeros(3)
        e[k] = 1.0
        T[:, 3 + k] = (np.cross(e, R) * sq[:, None]).ravel()
    Q, _ = np.linalg.qr(T)
    P = np.eye(3 * n) - Q @ Q.T
    w, V = np.linalg.eigh(P @ Hmw @ P)
    drop = np.argsort(np.abs(w))[:6]
    keep = np.setdiff1d(np.arange(3 * n), drop)
    w, V = w[keep], V[:, keep]
    order = np.argsort(w)
    w, V = w[order], V[:, order]
    freq = np.sign(w) * np.sqrt(np.abs(w)) * CONV
    fj, meta = freq_json(model)
    assert len(fj) == len(freq) == 3 * n - 6
    dev = np.abs(np.sort(freq) - np.sort(fj))
    if dev.max() > TOL_CM:
        raise AssertionError(f'{model}: recomputed wavenumbers differ from freq.json by up to {dev.max():.3f} cm-1')
    disp = (V / np.sqrt(m3)[:, None]).T.reshape(-1, n, 3)                    # (modes, atoms, 3)
    amp = np.linalg.norm(disp, axis=2).max(1)
    disp = disp / amp[:, None, None]
    out = dict(freq=freq, vecs_mw=V.T.reshape(-1, n, 3), disp=disp, masses=m, symbols=sym, freq_json=fj,
               max_dev=float(dev.max()), meta=meta)
    _CACHE[model] = out
    return out


def h_weight(md, k, h):
    """Fraction of the mass-weighted eigenvector norm on atom h."""
    v = md['vecs_mw'][k]
    return float((v[h] ** 2).sum() / (v ** 2).sum())


def stretch_fraction(md, k, a, h, model):
    """Fraction of the H displacement (relative to its partner a) that lies along the a-h bond."""
    X = C.geometries()[model]
    b = X[h] - X[a]
    b /= np.linalg.norm(b)
    d = md['disp'][k][h] - md['disp'][k][a]
    return float((d @ b) ** 2 / (d @ d))


def find_mode(md, target, tol=0.5):
    k = int(np.argmin(np.abs(md['freq'] - target)))
    if abs(md['freq'][k] - target) > tol:
        raise AssertionError(f'no mode within {tol} cm-1 of {target} (nearest {md["freq"][k]:.2f})')
    return k


def xh_stretches(model):
    """Local X-H stretch modes: for each H atom the mode with the largest H weight among the modes above
    2000 cm-1 (all X-H stretches of this pair). Returns dict label -> (mode index, wavenumber, H weight)."""
    md = modes(model)
    g = C.geometries()
    labels = g['labels']
    hi = np.flatnonzero(md['freq'] > 2000)
    out = {}
    for h in [i for i, s in enumerate(g['symbols']) if s == 'H']:
        ws = [h_weight(md, k, h) for k in hi]
        k = int(hi[int(np.argmax(ws))])
        out[labels[h]] = (k, float(md['freq'][k]), float(max(ws)))
    return out


def benzo_ch_modes(model, n=4):
    """The n modes above 2000 cm-1 with the largest summed weight on the benzo hydrogens H8-H11."""
    md = modes(model)
    hs = [C.I[h] for h in ('H8', 'H9', 'H10', 'H11')]
    hi = np.flatnonzero(md['freq'] > 2000)
    w = np.array([sum(h_weight(md, k, h) for h in hs) for k in hi])
    return sorted(int(k) for k in hi[np.argsort(w)[::-1][:n]])


# ---------------------------------------------------------------------------------------------------------------
# displacement arrows on the 2D ball-and-stick drawing
# ---------------------------------------------------------------------------------------------------------------
def draw_mode(ax, model, k, atoms=None, label_atoms=(), hbonds=(), scale_A=1.0, min_len=0.1, stubs=True,
              colour=None, arrow_lw=0.9, head=(3.2, 1.7), strict=True, phase=None, arrow_atoms=None,
              extra_segments=(), min_head_ratio=None, head_rule_outside=False):
    """Ball-and-stick drawing of `model` (common orientation) with the Cartesian displacement of mode k drawn as
    arrows (projected (x, y) of the displacement x scale_A; the largest atomic displacement is 1 before scaling;
    arrows shorter than min_len Å in projection are skipped). atoms: subset to draw (bonds leaving the subset
    are drawn as short grey stubs); arrow_atoms: atoms that get arrows (default all drawn atoms).
    Returns (drawing, n_arrows)."""
    from matplotlib.patches import FancyArrowPatch
    from matplotlib.lines import Line2D
    import matplotlib.patheffects as pe
    from figkit import draw2d
    colour = C.MODE_COLOUR if colour is None else colour
    g = C.geometries()
    md = modes(model)
    R, _ = C.views()[model]
    Xv = C.oriented(model)
    D = md['disp'][k]
    if phase is not None:                     # (atom, direction) the chosen phase should move that atom along
        i, v = phase
        if D[i] @ np.asarray(v, float) < 0:
            D = -D
    Dv = D @ np.asarray(R).T
    sym, labels = g['symbols'], g['labels']
    sub = list(range(len(sym))) if atoms is None else sorted(atoms)
    pos = {a: n for n, a in enumerate(sub)}
    all_bonds = geo.detect_bonds(sym, g[model])
    bonds = [(pos[i], pos[j]) for i, j in all_bonds if i in pos and j in pos]
    hb = [(pos[d], pos[h], pos[a]) for d, h, a in hbonds if d in pos and h in pos and a in pos]
    arrows, segs = [], list(extra_segments)       # extra_segments: label obstacles drawn by the caller
    # min_head_ratio: also skip arrows shorter (as printed) than this many head lengths, i.e. arrows that would be
    # mostly arrowhead (min_head_ratio=2: the shaft is at least as long as the head).  head_rule_outside=True
    # applies the rule to the part of the arrow outside the atom's ball (the arrow starts at the atom centre, so
    # a short arrow is half hidden in the ball and shows as a lone head; QA round 2, r2-m4)
    T = ax.transData
    pt_per_A = abs(T.transform((1.0, 0.0))[0] - T.transform((0.0, 0.0))[0]) * 72.0 / ax.figure.dpi
    head_pt = head[0]
    if head_rule_outside:          # the head as printed: the mitred tip and the stroke add to head_length
        head_pt = head[0] + arrow_lw / 2 / np.sin(np.arctan2(head[1], head[0])) + arrow_lw / 2
    for i in (sub if arrow_atoms is None else [a for a in sub if a in set(arrow_atoms)]):
        v = Dv[i, :2] * scale_A
        if np.hypot(*v) < min_len:
            continue
        vis_len = np.hypot(*v) - (style.BALL_RADII.get(sym[i], 0.35) if head_rule_outside else 0.0)
        if min_head_ratio is not None and vis_len * pt_per_A < min_head_ratio * head_pt:
            continue
        p0 = Xv[i, :2]
        arrows.append((p0, p0 + v))
        segs.append((p0, p0 + v, arrow_lw / 2 + 1.2))
    if stubs:
        for i, j in all_bonds:
            if (i in pos) != (j in pos):
                a, b = (i, j) if i in pos else (j, i)
                p0 = Xv[a, :2]
                p1 = Xv[a, :2] + 0.42 * (Xv[b, :2] - Xv[a, :2])
                ax.add_line(Line2D([p0[0], p1[0]], [p0[1], p1[1]], color='0.72', lw=1.6, zorder=1.5,
                                   solid_capstyle='round'))
                segs.append((p0, p1, 1.0))
    d = draw2d.draw_molecule(ax, [sym[i] for i in sub], Xv[sub], [labels[i] for i in sub], bonds=bonds,
                             hbonds=hb, label_atoms=list(label_atoms), label_h=True, hbond_distances=False,
                             fit=False, strict=strict, extra_segments=segs)
    if d.dropped:
        raise AssertionError(f'labels dropped: {d.dropped_labels()}')
    for p0, p1 in arrows:
        fa = FancyArrowPatch(p0, p1, arrowstyle=f'-|>,head_length={head[0]},head_width={head[1]}',
                             mutation_scale=1.0, color=colour, lw=arrow_lw, zorder=100, shrinkA=0, shrinkB=0,
                             joinstyle='miter', capstyle='butt')
        fa.set_path_effects([pe.Stroke(linewidth=arrow_lw + 1.2, foreground='white'), pe.Normal()])
        ax.add_patch(fa)
    return d, len(arrows)


# ---------------------------------------------------------------------------------------------------------------
# rigid-body decomposition of a mode (Fig. 6d, e; QA round 1, M3)
# ---------------------------------------------------------------------------------------------------------------
def rigid_decomposition(model, k, ions=(('picrate', C.PIC), ('cation', C.CAT))):
    """Least-squares fit of the Cartesian displacement of mode k, ion by ion, by a rigid motion d_a = t + w x r_a
    (r_a from the ion's centre of mass), in the mass-weighted metric; i.e. the projection of the mass-weighted
    eigenvector on the 12-dimensional space of rigid translations and rotations of the two ions (the space of
    `*.mode1_rigid_fraction`).  Returns dict(ions={name: dict(idx, com (Å), t, w (per unit of the normalised
    displacement `disp`), trans_share, rot_share)}, rigid_fraction)."""
    md = modes(model)
    X = C.geometries()[model]
    m = md['masses']
    q = md['vecs_mw'][k]
    amp = np.linalg.norm(q / np.sqrt(m)[:, None], axis=1).max()     # disp = q / sqrt(m) / amp
    out, rig = {}, 0.0
    for name, idx in ions:
        idx = list(idx)
        mm, x = m[idx], X[idx]
        com = (mm[:, None] * x).sum(0) / mm.sum()
        r = x - com
        A = np.zeros((3 * len(idx), 6))
        for a in range(len(idx)):
            sq = np.sqrt(mm[a])
            A[3 * a:3 * a + 3, :3] = sq * np.eye(3)
            rx = r[a]
            A[3 * a:3 * a + 3, 3:] = sq * np.array([[0, rx[2], -rx[1]], [-rx[2], 0, rx[0]], [rx[1], -rx[0], 0]])
        sol, *_ = np.linalg.lstsq(A, q[idx].ravel(), rcond=None)
        fit_t, fit_r = A[:, :3] @ sol[:3], A[:, 3:] @ sol[3:]
        out[name] = dict(idx=idx, com=com, t=sol[:3] / amp, w=sol[3:] / amp,
                         trans_share=float(fit_t @ fit_t), rot_share=float(fit_r @ fit_r))
        rig += float((A @ sol) @ (A @ sol))
    return dict(ions=out, rigid_fraction=rig / float((q ** 2).sum()))
