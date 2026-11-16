"""Fig. 7 fig_nci: non-covalent interactions of the benzotriazolium picrate ion pair.

(a, b) RDG = 0.5 isosurfaces of the B3LYP/6-311++G(d,p) density (Multiwfn rdg.cube), coloured by sign(lambda2)rho
(sl2r.cube) over -0.035...+0.02 a.u., gas phase and IEF-PCM (water), camera of the structure figures, one print
scale; numbered callouts on the patches.  (c, d) RDG against sign(lambda2)rho for all grid points
(rdg_scatter.npz), same colour map and range, the pipeline spikes (nci_spikes.csv) numbered like the patches.
Numbers: 1 N4-H4N...O1, 2 H4N...O7 and O1...O7 (one surface), 3 C11-H11...O2 (gas), 4 H11...O1 (gas),
5 O1...O2 (within the picrate), 6 triazole ring.  Patch and spike assignments are computed here and asserted
against facts.json.  A spike whose lowest-RDG points lie between two numbered pairs gets both numbers (IEF-PCM
-0.0145 a.u.: O1/O2 = 5 and O1/O7 = 2, as in Table 7); bins are never merged across sign(lambda2)rho = 0.

Callouts (QA round 2, R2-M1, r2-m3): O1 is labelled below its own sphere, every badge circle keeps >= 1 mm from
every atom label, and no leader crosses a bond (bonds projected and cut back to the ball rims); badge 6 (triazole
ring centre) is placed first, its leader leaving the ring over the N5 ball, which has no exocyclic bond.  All
three rules are asserted (check_callouts).  QA round 3 (r3-m2): that leader is drawn as passing BEHIND the N5 ball
(the part over the projected ball is left out, clip_leaders), and check_callouts also asserts that no drawn leader
passes over a ball other than the one it starts from.

Render field (QA round 1, M6): Multiwfn sets RDG = 100 wherever rho >= 0.1 a.u. (RDG_maxrho).  In the gas phase the
N4-H4N...O1 contact reaches that density, so the RDG = 0.5 contour snapped to the unscreened voxels and the
surface of interaction 1 was faceted and dented.  For the grid points between H4N and O1 (cylinder of radius
0.8 Å around the H4N->O1 axis, 10-90 % of the H...O distance) that were screened and have rho < 0.15 a.u., the
RDG is recomputed from final.molden (density and finite-difference gradient) before contouring; the recomputation
is validated on the unscreened points of the same cylinder (max |dRDG| < 1e-4) and it must not create an
RDG < 0.5 voxel next to a voxel that stays screened.  The scatter plots (c, d) keep the Multiwfn points.

    python scripts/figures/fig_nci.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from dataclasses import replace

from figkit import style, compose, callouts, rdg, render3d as r3, molden
from figkit.labels import LabelPlacer
import _render_common as C
import _struct_common as S

STEM = 'fig_nci'
ISO = 0.5
VRANGE = rdg.NCI_RANGE                          # (-0.035, 0.02) a.u., shared with the scatter colours
XLIM, YLIM = (-0.10, 0.06), (0.0, 1.5)
RAIL_BOTTOM, RAIL_Y = -0.25, -0.125          # strip below RDG = 0 for the spike numbers
BAND_MM = 47.0                                  # height of the scatter row
I = C.I
TRIAZOLE = [I['N4'], I['N5'], I['N6'], I['C7'], I['C12']]
# numbered interactions: atom pair (patch/spike assignment), text for the caption.  Patch 2 is ONE surface that
# covers both H4N...O7 and the O1...O7 contact inside the picrate (checked vertex by vertex: nearest atom pairs
# H4N/O7 45 % gas / 43 % PCM, O1/O7 11 % / 26 %), so O1...O7 belongs to 2, not to 5.
NUMBERS = {
    1: (('H4N', 'O1'), 'N4–H4N···O1'),
    2: (('H4N', 'O7'), 'H4N···O7 and O1···O7'),
    3: (('H11', 'O2'), 'C11–H11···O2'),
    4: (('H11', 'O1'), 'H11···O1'),
    5: (('O1', 'O2'), 'O1···O2'),
    6: (('N4', 'N6'), 'triazole ring'),
}
PAIR_NUMBER = {frozenset(('H4N', 'O1')): 1, frozenset(('H4N', 'O7')): 2, frozenset(('H11', 'O2')): 3,
               frozenset(('H11', 'O1')): 4, frozenset(('O1', 'O2')): 5, frozenset(('O1', 'O7')): 2,
               frozenset(('N4', 'N6')): 6, frozenset(('N4', 'C12')): 6}
# interaction numbers are drawn as small circled numbers (in the renders and in the scatter plots), so they can
# never be read as part of an atom label such as O1
BADGE = dict(boxstyle='circle,pad=0.14', fc='white', ec='#1A1A1A', lw=0.55)
BADGE_KW = dict(fontsize=style.TICK_PT, fontweight='bold', bbox=BADGE)
GAS_ONLY = (3, 4)
# atom labels: O1 below (left) of its own sphere, never between the sphere and surface 1 (QA round 2, R2-M1)
LABEL_PREFER = {'O1': dict(prefer=(-0.25, -1.0), prefer_weight=4.0)}
# gas phase: H4N above its ball (placed first), then N4 above its ball next to it, outside the triazole ring whose
# top is crossed by the leader of badge 6; badge 1 keeps its place under surface 1
LABEL_PREFER_GAS = {'H4N': dict(prefer=(0.0, 1.0), prefer_weight=2.0),
                    'N4': dict(prefer=(0.2, 1.0), prefer_weight=3.0)}
LABEL_BADGE_MM = 1.0               # least gap between an atom label and a badge circle (QA round 2, R2-M1)
RING_LEADER_MM = 15.0              # badge 6 (triazole ring centre): leader out of the ring without crossing a bond


def mm_in_px(ax):
    """Image pixels per mm at print size on a render axes."""
    T = ax.transData
    return ax.figure.dpi / 25.4 / float(abs(T.transform((1, 0))[0] - T.transform((0, 0))[0]))


def bond_segments(rd, into=0.0):
    """Projected bonds of the render (image px), each cut back to the rims of its two balls (into=0: the part of
    a bond that is seen) or reaching the fraction `into` of the ball radius into each ball (so that a leader
    must also keep clear of the point where a bond meets a ball).  A leader crossing one of these reads as an
    extra bond (QA round 2, r2-m3)."""
    mol = rd.mol
    xy, _ = rd.project(mol.coords)
    out = []
    for i, j in mol.bonds:
        a, b = xy[i], xy[j]
        L = float(np.hypot(*(b - a)))
        ri = (1 - into) * style.BALL_RADII[mol.symbols[i]] * rd.px_per_A
        rj = (1 - into) * style.BALL_RADII[mol.symbols[j]] * rd.px_per_A
        if L <= ri + rj:
            continue                                   # the bond is hidden between the two balls
        u = (b - a) / L
        out.append(np.r_[a + ri * u, b - rj * u])
    return np.array(out)


def ring_exit(rd, atom):
    """Screen direction (x right, y UP) from the triazole ring centre through `atom`."""
    xy, _ = rd.project(rd.mol.coords)
    d = xy[rd.mol.index(atom)] - xy[TRIAZOLE].mean(0)
    return np.array([d[0], -d[1]]) / np.hypot(*d)


def ball_discs(rd):
    """Projected atom balls of the render: centres (image px) and radii (image px)."""
    xy, _ = rd.project(rd.mol.coords)
    r = np.array([style.BALL_RADII.get(s, 0.35) for s in rd.mol.symbols]) * rd.px_per_A
    return np.asarray(xy, float), r


def clip_leaders(ax, outs):
    """QA round 3, r3-m2: a leader that crosses an atom ball (badge 6 over N5) is drawn as passing BEHIND it: the
    part over the projected ball is left out (butt ends at the rim), so it can be read neither as a line struck
    through the atom nor as an extra bond.  The placed leader (full segment) stays in o['leader'], invisible, for
    the placement checks; the drawn pieces are o['pieces']."""
    from matplotlib.lines import Line2D
    rd = ax._figkit_render
    xy, rad = ball_discs(rd)
    n_cut = 0
    for o in outs:
        ln = o['leader']
        o['pieces'] = [ln] if ln is not None else []
        if ln is None:
            continue
        xd, yd = ln.get_data()
        p0, p1 = np.array([xd[0], yd[0]]), np.array([xd[1], yd[1]])
        pieces = S.segment_outside_discs(p0, p1, xy, rad, own_tol=0.5)
        if len(pieces) == 1 and np.allclose(pieces[0][0], p0) and np.allclose(pieces[0][1], p1):
            continue
        ln.set_visible(False)
        o['pieces'] = []
        for a, b in pieces:
            pc = Line2D([a[0], b[0]], [a[1], b[1]], lw=ln.get_linewidth(), color=ln.get_color(),
                        zorder=ln.get_zorder(), solid_capstyle='butt', path_effects=ln.get_path_effects())
            ax.add_line(pc)
            o['pieces'].append(pc)
        n_cut += 1
    return n_cut


def check_callouts(ax, labs, badges, bonds, m):
    """QA round 2 checks: no leader (atom labels or badges) crosses a bond; every badge circle is at least
    LABEL_BADGE_MM from every atom label.  QA round 3 (r3-m2): no drawn leader piece passes over an atom ball
    (other than the ball the leader starts from)."""
    from figkit.callouts import _segs_cross
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    xy, rad = ball_discs(ax._figkit_render)
    for o in labs + badges:
        if o['leader'] is None:
            continue
        xd, yd = o['leader'].get_data()
        start = np.array([xd[0], yd[0]])
        assert not _segs_cross(start, np.array([xd[1], yd[1]]), bonds), \
            f'{m}: the leader of "{o["text"].get_text()}" crosses a bond'
        for pc in o.get('pieces', [o['leader']]):
            px, py = pc.get_data()
            hit = S.over_discs((px[0], py[0]), (px[1], py[1]), xy, rad, start, tol=0.5)
            assert not hit, (f'{m}: the leader of "{o["text"].get_text()}" is drawn over the ball of '
                             f'{[ax._figkit_render.mol.labels[k] for k in hit]}')
    px_mm = fig.dpi / 25.4
    worst = np.inf
    for b in badges:
        bb = b['text'].get_bbox_patch().get_window_extent(r)
        for a in labs:
            ab = a['text'].get_window_extent(r)
            dx = max(ab.x0 - bb.x1, bb.x0 - ab.x1, 0.0)
            dy = max(ab.y0 - bb.y1, bb.y0 - ab.y1, 0.0)
            gap = float(np.hypot(dx, dy)) / px_mm
            worst = min(worst, gap)
            assert gap >= LABEL_BADGE_MM, (f'{m}: badge {b["text"].get_text()} is {gap:.2f} mm from atom label '
                                           f'{a["text"].get_text()} (< {LABEL_BADGE_MM} mm)')
    print(f'{m}: no leader crosses a bond; closest atom label to a badge {worst:.2f} mm')


# ---------------------------------------------------------------------------------------------------------------
# render field: the rho >= 0.1 a.u. screening undone locally in the N4-H4N...O1 contact (QA round 1, M6)
# ---------------------------------------------------------------------------------------------------------------
RDG_CONST = 1.0 / (2.0 * (3.0 * np.pi ** 2) ** (1.0 / 3.0))
UNSCREEN_RHO_MAX = 0.15
UNSCREEN_RADIUS = 0.8                             # Å around the H4N -> O1 axis


def _density_and_gradient(mf, P, h=1e-3):
    """rho and grad rho (a.u.) at points P (Bohr) from the occupied MOs of a Molden file (central differences)."""
    occ = np.flatnonzero(mf.mo_occ > 0)
    w = mf.mo_occ[occ]

    def rho(Q):
        return (w[:, None] * mf.mo_on_points(occ, Q) ** 2).sum(0)
    g = np.zeros((len(P), 3))
    for k in range(3):
        e = np.zeros(3)
        e[k] = h
        g[:, k] = (rho(P + e) - rho(P - e)) / (2 * h)
    return rho(P), g


def unscreened_cubes(model, mol):
    """(rdg Cube, sl2r Cube, info): Multiwfn's cubes with the RDG of the screened grid points between H4N and O1
    (rho < UNSCREEN_RHO_MAX) recomputed from the wavefunction.  sl2r.cube is not screened and is used as it is."""
    d = r3.case_dir(C.CASES[model])
    rc, sc = r3.read_cube(d / 'rdg.cube'), r3.read_cube(d / 'sl2r.cube')
    mf = molden.read_molden(d / 'final.molden')
    assert np.abs(mf.coords * molden.BOHR - mol.coords).max() < 1e-6, 'final.molden is not in the opt.xyz frame'
    H, O = mol.coords[I['H4N']], mol.coords[I['O1']]
    L = float(np.linalg.norm(O - H))
    u = (O - H) / L
    X = rc.origin + np.indices(rc.data.shape).reshape(3, -1).T @ rc.axes
    t = (X - H) @ u
    rad = np.linalg.norm((X - H) - np.outer(t, u), axis=1)
    sel = np.flatnonzero((rad < UNSCREEN_RADIUS) & (t > 0.1 * L) & (t < 0.9 * L))
    rho, g = _density_and_gradient(mf, X[sel] / molden.BOHR)
    rdg_new = RDG_CONST * np.linalg.norm(g, axis=1) / rho ** (4.0 / 3.0)
    R1, S1 = rc.data.ravel(), sc.data.ravel()
    screened = R1[sel] >= 99.99
    ok = ~screened & (R1[sel] < 3.0)
    dev = float(np.abs(rdg_new[ok] - R1[sel][ok]).max())
    assert dev < 1e-4, f'{model}: recomputed RDG differs from rdg.cube by {dev}'
    assert np.abs(rho - np.abs(S1[sel])).max() < 1e-3, 'density differs from |sl2r|'
    assert np.all(rho[screened] >= 0.1 - 1e-3), 'screened point below the screening density'
    fill = screened & (rho < UNSCREEN_RHO_MAX)
    new = R1.copy()
    new[sel[fill]] = rdg_new[fill]
    G = new.reshape(rc.data.shape)
    still = G >= 99.99
    for v in sel[fill]:                            # no new surface may close against a still-screened voxel
        i, j, k = np.unravel_index(v, G.shape)
        if G[i, j, k] < ISO:
            nb = [still[i + a, j + b, k + c] for a, b, c in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0),
                                                              (0, 0, 1), (0, 0, -1))]
            assert not any(nb), f'{model}: unscreened voxel with RDG < {ISO} next to a screened one'
    info = dict(n_filled=int(fill.sum()), n_below_iso=int((rdg_new[fill] < ISO).sum()),
                min_rdg=float(rdg_new[fill].min()) if fill.any() else None, check=dev)
    return replace(rc, data=G), sc, info


# ---------------------------------------------------------------------------------------------------------------
# 3D patches
# ---------------------------------------------------------------------------------------------------------------
def patches(scene, mol):
    """{number: vertex array (A)} for the numbered RDG = 0.5 patches (connected components of the isosurface)."""
    surf = next(ms for ms, k in scene.items if 'sl2r' in ms.point_data)
    conn = surf.connectivity('all')
    ids = np.asarray(conn.point_data['RegionId'])
    P = np.asarray(conn.points)
    ring_c = mol.coords[TRIAZOLE].mean(0)
    found = {}
    for k in range(ids.max() + 1):
        pts = P[ids == k]
        cen = pts.mean(0)
        if np.linalg.norm(cen - ring_c) < 0.5:
            n = 6
        else:
            near = [mol.labels[i] for i in np.argsort(np.linalg.norm(mol.coords - cen, axis=1))[:2]]
            n = PAIR_NUMBER.get(frozenset(near))
            if n == 6:
                n = None
        if n is None:
            continue
        assert n not in found, f'two patches for interaction {n}'
        found[n] = pts
    return found


# ---------------------------------------------------------------------------------------------------------------
# RDG spikes
# ---------------------------------------------------------------------------------------------------------------
def spike_pairs(note):
    """Atom pairs of a facts spike note: 'lowest-RDG points lie between O2/O1, O7/O1 (derived ...)'."""
    m = re.search(r'lie between (.*?) \(derived', note)
    return [tuple(t.strip().split(' ')[0].split('/')) for t in m.group(1).split(',')]


def spikes(model, mol):
    """Pipeline spikes with the facts assignment ('lowest-RDG points lie between A/B, C/D ...') turned into
    interaction numbers (every numbered pair counts, so one spike can carry two numbers); neighbouring bins
    (< 0.01 a.u. apart) with the same numbers ON THE SAME SIDE of zero are merged into the lowest-RDG one."""
    d = r3.case_dir(C.CASES[model])
    raw = rdg.read_spikes(d)
    out = []
    for k, s in enumerate(raw, start=1):
        f = C.facts()[f'{model}.nci_spike{k}']
        assert abs(f['value'] - s['value']) < 1e-9, (model, k)
        nums = []
        for a, b in spike_pairs(f['note']):
            n = PAIR_NUMBER.get(frozenset((a, b)))
            if n is not None and n not in nums:
                nums.append(n)
        assert nums, (model, k, f['note'])
        out.append(dict(s, numbers=nums, label=' + '.join(NUMBERS[n][1] for n in nums),
                        key=f'{model}.nci_spike{k}'))
    merged = []
    for s in sorted(out, key=lambda s: s['value']):
        if merged and merged[-1]['numbers'] == s['numbers'] and abs(merged[-1]['value'] - s['value']) < 0.0101 \
                and np.sign(merged[-1]['value']) == np.sign(s['value']):
            if s['rdg'] < merged[-1]['rdg']:
                merged[-1] = s
            continue
        merged.append(s)
    return merged


def label_spikes(ax, sp, name):
    """Mark each spike (lowest-RDG point of its sign(l2)rho bin) with a small dot and put its circled interaction
    number on a rail in the empty strip below RDG = 0, straight under the spike where possible (neighbours
    closer than a badge are spread symmetrically), joined to the dot by a thin leader.  Returns the badge texts."""
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    T = ax.transData
    y_rail = RAIL_Y
    # badge diameter in data x units, measured on a drawn badge
    probe = ax.text(0, y_rail, '8', ha='center', va='center', zorder=6, **BADGE_KW)
    fig.canvas.draw()
    bb = probe.get_bbox_patch().get_window_extent(rend)
    probe.remove()
    inv = T.inverted()
    dx = abs(inv.transform((bb.x1, 0))[0] - inv.transform((bb.x0, 0))[0])
    dy = abs(inv.transform((0, bb.y1))[1] - inv.transform((0, bb.y0))[1])
    sep = 1.25 * dx
    # one badge per (spike, number); a spike with two numbers gets two badges side by side
    entries = [(s['value'], j, s, n) for s in sp for j, n in enumerate(s['numbers'])]
    entries.sort(key=lambda e: (e[0], e[1]))
    xs = np.array([e[0] for e in entries], float)
    for _ in range(200):                                   # 1D spreading: push neighbours apart
        moved = False
        for i in range(len(xs) - 1):
            gap = xs[i + 1] - xs[i]
            if gap < sep - 1e-12:
                d = (sep - gap) / 2
                xs[i] -= d
                xs[i + 1] += d
                moved = True
        lo, hi = XLIM[0] + 0.6 * dx, XLIM[1] - 0.6 * dx
        xs = np.clip(xs, lo, hi)
        if not moved:
            break
    assert np.all(np.diff(xs) >= sep - 1e-9), 'spike numbers do not fit on the rail'
    out = []
    for s in sp:
        ax.plot(s['value'], s['rdg'], 'o', ms=2.6, mfc='#1A1A1A', mec='white', mew=0.4, zorder=5)
    for xb, (_, _, s, n) in zip(xs, entries):
        ax.plot([xb, s['value']], [y_rail + 0.5 * dy, s['rdg']], color='0.2', lw=0.5, zorder=4,
                solid_capstyle='butt')
        out.append(ax.text(xb, y_rail, str(n), ha='center', va='center', zorder=6, **BADGE_KW))
    ax.axhline(0, color='0.75', lw=0.4, zorder=1)
    return out


def annotate_render(ax, m, found_m, strict=True):
    """Atom labels and circled interaction numbers around one render.  Order: badge 6 first (the centre of the
    triazole ring; to cross no bond its leader must leave the ring over the N5 ball, which has no exocyclic
    bond, so it needs a longer leader than the others), then the atom labels (clear of that leader and at least
    LABEL_BADGE_MM from its badge), then the other badges.  QA round 2: no leader crosses a bond (r2-m3: they
    read as extra bonds) and every badge circle keeps LABEL_BADGE_MM from every atom label (R2-M1: 'O1' next
    to badge 1 read 'O1(1)'); O1 is labelled below its own sphere."""
    rd = ax._figkit_render
    items = {}
    for n in sorted(found_m):
        pts = found_m[n]
        tgt = callouts.patch_target(rd, pts)
        assert tgt is not None, (m, n)
        vis = rd.visible(pts, towards_camera=0.02)
        fp, _ = rd.project(pts[vis])
        items[n] = dict(target=tgt, text=str(n), footprint=fp, dot=False, kw=dict(BADGE_KW))
    bonds = bond_segments(rd)
    mm_px = mm_in_px(ax)
    grow = lambda boxes: [b + np.array([-1, -1, 1, 1]) * (LABEL_BADGE_MM + 0.05) * mm_px for b in boxes]

    def patch_box(o):
        fig = ax.figure
        r = fig.canvas.get_renderer()
        bb = o['text'].get_bbox_patch().get_window_extent(r)
        (x0, y0), (x1, y1) = ax.transData.inverted().transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        return np.array([min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)])

    kw = dict(min_leader_mm=0.8, pad_mm=0.6, clearance_mm=0.6, patch_extent=True, strict=strict)
    # (bonds reaching half-way into the balls, so the leader cannot slip out past a ring-bond junction)
    six = callouts.place(ax, [dict(items[6], prefer=tuple(ring_exit(rd, 'N5')), prefer_weight=20.0)],
                         max_leader_mm=RING_LEADER_MM, n_dirs=72, name=f'NCI number 6 {m}',
                         avoid_segments=bond_segments(rd, into=0.5), **kw)
    lead6 = [np.r_[o['leader'].get_data()[0][0], o['leader'].get_data()[1][0], o['leader'].get_data()[0][1],
                   o['leader'].get_data()[1][1]] for o in six if o['leader'] is not None]
    obst = np.vstack([bonds] + [np.atleast_2d(l) for l in lead6])
    ax.figure.canvas.draw()
    labs = []
    order = ('O1', 'O7', 'O2', 'N4', 'H4N') if m == 'pcm' else ('O1', 'O7', 'O2', 'H4N', 'N4', 'H11')
    for lab in order:
        k = rd.mol.index(lab)
        it = dict(pos=rd.mol.coords[k], text=lab, dot=False,
                  radius=style.BALL_RADII[rd.mol.symbols[k]] * rd.px_per_A)
        it.update(LABEL_PREFER.get(lab, {}))
        if m == 'gas':
            it.update(LABEL_PREFER_GAS.get(lab, {}))
        labs.append(it)
    # the leader of badge 6 as a chain of small boxes, so that no atom label sits on it
    beads = []
    for x0, y0, x1, y1 in lead6:
        n = max(2, int(np.hypot(x1 - x0, y1 - y0) / (0.4 * mm_px)) + 1)
        for t in np.linspace(0, 1, n):
            c = np.array([x0 + t * (x1 - x0), y0 + t * (y1 - y0)])
            beads.append(np.r_[c - 0.5 * mm_px, c + 0.5 * mm_px])
    placed_labs = callouts.place(ax, labs, min_leader_mm=0.3, max_leader_mm=5.0, no_leader_below_mm=0.9,
                                 name=f'atom labels {m}', avoid_segments=obst,
                                 avoid_boxes=grow([patch_box(o) for o in six]) + beads, strict=strict)
    rest = callouts.place(ax, [items[n] for n in sorted(items) if n != 6], max_leader_mm=9.0,
                          name=f'NCI numbers {m}', avoid_segments=obst,
                          avoid_boxes=grow([o['box'] for o in placed_labs]), **kw)
    n_cut = clip_leaders(ax, six + rest + placed_labs)
    print(f'{m}: {n_cut} leader(s) drawn as passing behind an atom ball')
    check_callouts(ax, placed_labs, six + rest, bonds, m)


def main():
    t0 = time.time()
    style.apply()
    mols = C.molecules()
    V = C.views()
    scenes, found = {}, {}
    for m in C.MODELS:
        rc, sc, info = unscreened_cubes(m, mols[m])
        print(f'{m}: RDG recomputed for {info["n_filled"]} screened points between H4N and O1 (check '
              f'{info["check"]:.1e}); {info["n_below_iso"]} below RDG = {ISO}, lowest {info["min_rdg"]:.3f}')
        scenes[m] = r3.nci_scene(rc, sc, mols[m], iso=ISO, vrange=VRANGE)
        found[m] = patches(scenes[m], mols[m])
        expect = set(NUMBERS) - (set(GAS_ONLY) if m == 'pcm' else set())
        assert set(found[m]) == expect, (m, sorted(found[m]))
        # the H-bond patch holds the strongest attraction; its value must match the real-space region in facts
        sl = scenes[m].meta['sl2r_min']
        print(f'{m}: patches {sorted(found[m])}, sign(l2)rho on the surface down to {sl:.4f}')
    vw = {m: V[m].fit(scenes[m], mols[m], pad=0.3) for m in C.MODELS}
    pxA = r3.px_per_angstrom(list(vw.values()), C.PANEL2_MM)
    rr = {}
    for m in C.MODELS:
        full = r3.render(scenes[m], vw[m], px_per_A=pxA, crop=False)
        vis = C.check_view(vw[m], mols[m], scenes[m], full, [(I['H4N'], I['O1']), (I['H4N'], I['O7'])])
        print(m, 'visible fraction of the H-bond patches', {f'{mols[m].labels[a]}-{mols[m].labels[b]}': round(v, 2)
                                                           for (a, b), v in vis.items()})
        rr[m] = r3.crop_together(full, pad=int(C.mm_to_px(1.6)))
    fig, axes = compose.grid([rr['gas'], rr['pcm']], ncols=2, width=style.DOUBLE,
                             titles=[style.MODEL_LABELS['gas'], style.MODEL_LABELS['pcm']],
                             colourbar_scale=scenes['pcm'].colour, colourbar_len=0.55, extra_bottom_mm=BAND_MM)
    for m, ax in zip(C.MODELS, axes):
        annotate_render(ax, m, found[m])

    # ---------------- (c, d) RDG scatter ----------------
    x0, y0, w, h = fig._figkit_grid['bottom_band']
    W, H = fig.get_size_inches()
    left_in, gap_in = 10.5 * style.MM, 3.0 * style.MM
    head_in = compose._header_mm() * style.MM
    bottom_in = 9.6 * style.MM
    aw = (w - left_in - gap_in) / 2
    ah = h - head_in - bottom_in - 1.0 * style.MM
    numbering = {NUMBERS[n][1]: n for n in NUMBERS}
    for j, m in enumerate(C.MODELS):
        ax_x = x0 + left_in + j * (aw + gap_in)
        ax = fig.add_axes([ax_x / W, (y0 + bottom_in) / H, aw / W, ah / H])
        data = rdg.load_rdg(r3.case_dir(C.CASES[m]), xlim=XLIM, ylim=YLIM)
        sp = spikes(m, mols[m])
        rdg.plot_rdg(ax, data['x'], data['y'], xlim=XLIM, ylim=YLIM, vrange=VRANGE, cmap=r3.NCI_CMAP,
                     spikes=None, key=False, ylabel='RDG' if j == 0 else None)
        ax.set_ylim(RAIL_BOTTOM, YLIM[1])                  # strip below RDG = 0 holds the spike numbers
        for ln in ax.lines:                                # the dashed zero line stops at the strip
            if list(ln.get_xdata()) == [0, 0]:
                ln.set_ydata([(0 - RAIL_BOTTOM) / (YLIM[1] - RAIL_BOTTOM), 1])
        label_spikes(ax, sp, name=f'spike numbers {m}')
        if j == 1:
            ax.set_yticklabels([])
        ax.set_yticks([0, 0.5, 1.0, 1.5])
        ax.set_xticks([-0.10, -0.05, 0.0, 0.05])
        ax.xaxis.set_major_formatter(compose._formatter('{:.2f}'))
        compose._header(fig, (x0 if j == 0 else ax_x) / W, (y0 + h) / H, f'({"cd"[j]})',
                        style.MODEL_LABELS[m])
        print(m, 'spikes', [(s['numbers'], s['value'], s['rdg']) for s in sp])
        if m == 'gas':                      # QA round 1, M5: the +0.0045 spike (Table 7) is numbered 4
            assert any(s['value'] > 0 and 4 in s['numbers'] for s in sp), sp
        if m == 'pcm':                      # the -0.0145 spike carries 5 (O1...O2) and 2 (O1...O7), as Table 7
            assert any(abs(s['value'] + 0.0145) < 1e-6 and s['numbers'] == [5, 2] for s in sp), sp
    C.save(fig, STEM)
    print(f'total {time.time() - t0:.1f} s')


if __name__ == '__main__':
    main()
