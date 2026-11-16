"""ELF / LOL 2 x 2 figure builder shared by fig_elf_lol (cation plane) and figS_elf_lol_picrate (picrate plane).

Rows gas phase / IEF-PCM (water), columns ELF / LOL, one colour bar 0-1. All four panels share one rotation, one
centre (the centroid of a chosen ring of each model) and one window, so the reference ion is drawn identically in
every panel. The rotation is the grid-aligned one (0/90/180/270 deg; an oblique window would be cut by the
Multiwfn grid) nearest to the in-plane angle of the Fig. 3 camera (build/views_struct.json); the mirror follows
that camera too (plane normal towards the camera <=> no mirror).
"""
from __future__ import annotations

import json
import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from figkit import style, maps
import _elec_common as C

LEVELS = np.linspace(0, 1, 21)                       # filled contours every 0.05
LINE_LEVELS = np.round(np.linspace(0.1, 0.9, 9), 2)  # lines every 0.1
ATOM_TOL = 0.5


def camera_angle(plane, X, model, idx):
    """In-plane rotation (deg) and mirror that best reproduce the Fig. 3 camera for the atoms idx."""
    v = json.loads((C.BUILD / 'views_struct.json').read_text())[model]
    d, up, f = (np.array(v[k]) for k in ('direction', 'up', 'focal'))
    right = np.cross(-d, up)
    right /= np.linalg.norm(right)
    S = np.c_[(X - f) @ right, (X - f) @ up][idx]
    uv, _ = plane.project(X)
    A = uv[idx] - uv[idx].mean(0)
    B = S - S.mean(0)
    mirror = bool(plane.normal @ d < 0)
    if mirror:
        A = A @ np.diag([1.0, -1.0])
    U, _, Vt = np.linalg.svd(A.T @ B)
    R = (U @ Vt).T
    assert np.linalg.det(R) > 0
    return float(np.degrees(np.arctan2(R[1, 0], R[0, 0]))), mirror


def load(plane_name, defining, ring, ring_tol):
    """Plane grids + geometries, with checks: ELF/LOL share one frame, the three defining atoms lie in the plane,
    the ring atoms within ring_tol, map ranges equal facts `*.{elf,lol}_<plane>_range`."""
    G = C.geometries()
    out = {}
    for m in C.MODELS:
        f = maps.plane_files(G['case'][m]['paths'].raw, plane_name)
        with warnings.catch_warnings():
            warnings.simplefilter('error')              # a non-orthonormal frame or an npz fallback would warn
            elf, lol = maps.read_plane(f['elf']), maps.read_plane(f['lol'])
        assert np.abs(elf.origin - lol.origin).max() < 1e-6 and np.abs(elf.normal - lol.normal).max() < 1e-9
        assert np.array_equal(elf.u, lol.u) and np.array_equal(elf.v, lol.v)
        assert elf.fit_residual < 1e-3
        X = G[m]
        _, dist = elf.project(X)
        assert np.abs(dist[defining]).max() < 2e-3, (m, dist[defining])
        assert np.abs(dist[ring]).max() < ring_tol, (m, np.abs(dist[ring]).max())
        for kind, P in (('elf', elf), ('lol', lol)):
            lo, hi = (float(s) for s in C.fact(f'{m}.{kind}_{plane_name}_range').split('..'))
            assert f'{np.nanmin(P.F):.4f}' == f'{lo:.4f}' and f'{np.nanmax(P.F):.4f}' == f'{hi:.4f}', (m, kind)
        out[m] = dict(elf=elf, lol=lol, X=X, dist=dist)
    return G, out


def common_geometry(G, data, centre_idx, cam_idx, margin=1.0):
    """One rotation, centre (centroid of centre_idx) and window for all four panels."""
    angles = {m: camera_angle(data[m]['elf'], data[m]['X'], m, cam_idx) for m in C.MODELS}
    for m, (a, mi) in angles.items():
        print(f'{m}: Fig. 3 camera = in-plane rotation {a:.1f} deg, mirror {mi}')
    mirror = angles['gas'][1]
    assert angles['pcm'][1] == mirror
    mean = np.degrees(np.angle(sum(np.exp(1j * np.radians(a)) for a, _ in angles.values())))
    rot = float(90 * np.round(mean / 90) % 360)
    print(f'mean camera angle {mean:.1f} deg -> grid-aligned rotation {rot:.0f} deg')

    def geom(m, crop, **kw):
        X = data[m]['X']
        return maps.map_geometry(data[m]['elf'], G['symbols'], X, atom_tol=ATOM_TOL, rotate=rot, crop=crop,
                                 centre=X[centre_idx].mean(0), mirror=mirror, **kw)

    wins = [np.array(geom(m, 'atoms', crop_margin=margin)['window']) for m in C.MODELS]
    U = (min(w[0] for w in wins), max(w[1] for w in wins), min(w[2] for w in wins), max(w[3] for w in wins))
    fitted = [np.array(geom(m, U)['window']) for m in C.MODELS]           # each cut to its own grid
    win = (max(f[0] for f in fitted), min(f[1] for f in fitted), max(f[2] for f in fitted), min(f[3] for f in fitted))
    geoms = {}
    for m in C.MODELS:
        geoms[m] = geom(m, win)
        xy, ins = geoms[m]['xy'], geoms[m]['inplane']
        inside = (xy[:, 0] > win[0] + 0.15) & (xy[:, 0] < win[1] - 0.15) & (xy[:, 1] > win[2] + 0.15) & \
                 (xy[:, 1] < win[3] - 0.15)
        lost = [G['labels'][i] for i in np.nonzero(ins & ~inside)[0]]
        assert not lost, f'{m}: in-plane atoms outside the common window: {lost}'
    return geoms, win, rot, mirror


def figure(G, data, geoms, win, stem, label_h=(), hbonds=(C.HB,), width_mm=175.0,
           left=10.3, gapx=3.0, cb_gap=2.2, cb_w=3.0, cb_lab=9.9, bottom=10.2, title=5.0, gapy=3.2, extra=None):
    """Draw and save the 2 x 2 figure. label_h: hydrogen atoms (indices) that get a label; extra: optional
    callback extra(axes) (axes = {(model, 'ELF'|'LOL'): Axes}) run before saving."""
    sym, labels = G['symbols'], G['labels']
    Wmm = width_mm
    pw = (Wmm - left - gapx - cb_gap - cb_w - cb_lab) / 2
    ph = pw * (win[3] - win[2]) / (win[1] - win[0])
    Hmm = bottom + 2 * ph + gapy + 2 * title
    fig = plt.figure(figsize=(C.mm(Wmm), C.mm(Hmm)))
    axes = {}
    for r, m in enumerate(C.MODELS):
        y = bottom + (1 - r) * (ph + gapy + title)
        for c, kind in enumerate(('ELF', 'LOL')):
            axes[(m, kind)] = fig.add_axes([(left + c * (pw + gapx)) / Wmm, y / Hmm, pw / Wmm, ph / Hmm])
    cax = fig.add_axes([(left + 2 * pw + gapx + cb_gap) / Wmm, bottom / Hmm, cb_w / Wmm,
                        (2 * ph + gapy + title) / Hmm])
    out = {}
    letters = iter('abcd')
    for m in C.MODELS:
        ins = geoms[m]['inplane']
        lab_atoms = [i for i in np.nonzero(ins)[0] if sym[i] != 'H'] + [h for h in label_h if ins[h]]
        for kind in ('ELF', 'LOL'):
            ax = axes[(m, kind)]
            bottom_row, left_col = m == C.MODELS[-1], kind == 'ELF'
            pm = maps.plot_plane_map(ax, data[m][kind.lower()], sym, data[m]['X'], labels, kind=kind,
                                     geom=geoms[m], levels=LEVELS, line_levels=LINE_LEVELS, atom_tol=ATOM_TOL,
                                     hbonds=list(hbonds), colorbar=False, label_atoms=lab_atoms,
                                     xlabel='x (Å)' if bottom_row else None, ylabel='y (Å)' if left_col else None,
                                     strict=True)
            out[(m, kind)] = pm
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.tick_params(labelbottom=bottom_row, labelleft=left_col, length=2.5)
            C.panel_title(ax, next(letters), f'{style.MODEL_LABELS[m]}, {kind}', pad_pt=4.8)   # QA r4-m3: about 1 mm between title descenders and frame
    cb = fig.colorbar(out[('gas', 'ELF')].mappable, cax=cax, ticks=np.linspace(0, 1, 11))
    cb.set_label('ELF, LOL')
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=style.TICK_PT, length=2.5)
    if extra is not None:
        extra(axes)
    C.save(fig, stem, width_mm)
    return out
