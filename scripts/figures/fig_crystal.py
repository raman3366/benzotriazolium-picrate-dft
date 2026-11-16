"""Fig. 2 fig_crystal: the modelled ion pair in its crystal context (P2_1/c, CCDC 835950, 298 K).

(a) asymmetric unit as deposited: the two ions are not hydrogen bonded; closest contact O4...N5 (dotted);
(b) the only ring stacking, triazole...triazole between two inversion-related cations, tilted 38 deg from edge-on;
    dotted centroid line, its length 4.106 Å given in the key (QA round 4, r4-m4: a label placed near the rings read
    as a label of N5), interplanar separation 3.341 Å as an arrow
    beside the rings with ticks on the two ring planes; the key for (b) and (c) to the right of (b);
(c) five ions of the N-H...O chain along b (P <- C -> P <- C -> P), drawn with b vertical as a column right of (a)
    and (b), the modelled pair in colour, the other ions in paler element colours, with an arrow for b.
QA round 1 (M1, m2-m4): the unit-cell projection down a (19 wireframe ions) was dropped; it read as a tangle at
print size and the chain panel carries the packing message.  QA round 2 (r2-m1): two columns instead of two rows
(no empty block; 94 instead of 131 mm high), scales solved so that both columns fill the width and end together.

Everything is generated from BTPC.cif with gemmi (CIF symmetry operators, symmetry codes n_klm = operator n plus
translation (k-5, l-5, m-5)); the chain is found by a breadth-first search over N-H...O bonds (H...A < 2.5 Å,
D-H...A > 100 deg). The script asserts that the highlighted pair equals btpc_cluster.xyz (RMSD < 0.01 Å) and that
every printed distance equals facts.json at the printed precision.

    python scripts/figures/fig_crystal.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gemmi
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from figkit import style, draw2d, geometry as geo
import _struct_common as C

STEM = 'fig_crystal'
WIDTH_MM = 175.0
GAP_MM = 4.0
HEAD_MM = 5.0
CIF = C.REPO / 'molecules' / 'btpc' / 'input' / 'BTPC.cif'
TEAL = draw2d.HBOND_COLOUR
GREY = C.HBOND_GREY
LIGHT = 0.62                    # lightening of context ions
CELL_COLOUR = '0.25'
S1_MAX = 6.2                   # mm per Å: cap for (a) and (b)
ROW_GAP_MM = 2.5               # between (a) and the header of (b)
STACK_TILT = 38.0               # (b): tilt from edge-on (QA round 1, m2: 22 deg overlapped the rings)


# ---------------------------------------------------------------------------------------------------------------
# crystal
# ---------------------------------------------------------------------------------------------------------------
class Crystal:
    def __init__(self):
        b = gemmi.cif.read(str(CIF)).sole_block()
        ss = gemmi.make_small_structure_from_block(b)
        self.cell = ss.cell
        self.ops = [gemmi.Op(t.strip("'")) for t in b.find_values('_symmetry_equiv_pos_as_xyz')]
        want = [t.strip().replace(' ', '') for t in str(C.fact('xray.symmetry_operations')).split(';')]
        assert [o.triplet() for o in self.ops] == want, ([o.triplet() for o in self.ops], want)
        self.labels = [s.label for s in ss.sites]
        self.symbols = [s.type_symbol for s in ss.sites]
        self.F = np.array([[s.fract.x, s.fract.y, s.fract.z] for s in ss.sites])
        self.M = np.array(ss.cell.orth.mat.tolist())
        X = self.F @ self.M.T
        bonds = geo.detect_bonds(self.symbols, X)
        comp = _components(len(self.symbols), bonds)
        assert len(comp) == 2, 'expected two ions in the asymmetric unit'
        self.sites = {('P' if 'O1' in [self.labels[k] for k in c] else 'C'): sorted(c) for c in comp}
        self.bonds = {k: [(i, j) for i, j in bonds if i in v and j in v] for k, v in self.sites.items()}

    def ion(self, kind, op, shift):
        """Ion `kind` ('P' picrate, 'C' cation) transformed by CIF operator number `op` (1-based) and translated
        by `shift` (fractional); its symmetry code is op_(5+k)(5+l)(5+m)."""
        o = self.ops[op - 1]
        R = np.array(o.rot, float) / o.DEN
        t = np.array(o.tran, float) / o.DEN
        idx = self.sites[kind]
        Fi = self.F[idx] @ R.T + t + np.asarray(shift, float)
        code = f"{op}_{''.join(str(5 + int(round(s))) for s in shift)}"
        return Ion(kind, code, idx, Fi, Fi @ self.M.T, self)

    def images(self, lo=-1.6, hi=2.6):
        out = []
        for op in range(1, len(self.ops) + 1):
            for kind in ('P', 'C'):
                base = self.ion(kind, op, (0, 0, 0))
                c0 = base.F.mean(0)
                for sh in itertools.product(range(-3, 4), repeat=3):
                    c = c0 + np.array(sh)
                    if np.all(c > lo) and np.all(c < hi):
                        out.append(self.ion(kind, op, sh))
        return out


class Ion:
    def __init__(self, kind, code, idx, F, X, cr):
        self.kind, self.code, self.idx, self.F, self.X, self.cr = kind, code, idx, F, X, cr
        self.lab = [cr.labels[k] for k in idx]
        self.sym = [cr.symbols[k] for k in idx]

    def pos(self, label):
        return self.X[self.lab.index(label)]

    def same(self, other):
        return self.kind == other.kind and np.allclose(self.F.mean(0), other.F.mean(0), atol=1e-4)

    def __repr__(self):
        return f'{self.kind}({self.code})'


def _components(n, bonds):
    adj = {i: [] for i in range(n)}
    for i, j in bonds:
        adj[i].append(j)
        adj[j].append(i)
    seen, comps = set(), []
    for i in range(n):
        if i in seen:
            continue
        st, c = [i], []
        while st:
            k = st.pop()
            if k in seen:
                continue
            seen.add(k)
            c.append(k)
            st += adj[k]
        comps.append(c)
    return comps


def nh_o(ions, acceptors=('O1', 'O6'), max_ha=2.5, min_angle=100.0):
    """N-H...O bonds between the ions (indices into `ions`): dicts don, d, h, acc, a, ha, da, angle."""
    out = []
    for i, ci in enumerate(ions):
        if ci.kind != 'C':
            continue
        for d, h in (('N4', 'H4N'), ('N6', 'H6N')):
            D, H = ci.pos(d), ci.pos(h)
            for j, pj in enumerate(ions):
                if pj.kind != 'P':
                    continue
                for a in acceptors:
                    A = pj.pos(a)
                    ha = float(np.linalg.norm(A - H))
                    if ha > max_ha:
                        continue
                    ang = float(np.degrees(np.arccos(np.dot(D - H, A - H) / (np.linalg.norm(D - H) * ha))))
                    if ang > min_angle:
                        out.append(dict(don=i, d=d, h=h, acc=j, a=a, ha=ha, da=float(np.linalg.norm(A - D)),
                                        angle=ang))
    return out


def triazole_centroid(ion):
    return np.mean([ion.pos(n) for n in ('N4', 'N5', 'N6', 'C7', 'C12')], 0)


def triazole_normal(ion):
    P = np.array([ion.pos(n) for n in ('N4', 'N5', 'N6', 'C7', 'C12')])
    return np.linalg.svd(P - P.mean(0))[2][2]


def find(ions, ion):
    for k, o in enumerate(ions):
        if o.same(ion):
            return k
    return None


# ---------------------------------------------------------------------------------------------------------------
# drawing several ions in one depth-sorted ball-and-stick picture
# ---------------------------------------------------------------------------------------------------------------
class Scene:
    """Ions concatenated for one draw2d.draw_molecule call, in a viewing frame (R rows: right, up, to viewer)."""

    def __init__(self, ions, model, R, centre):
        self.ions, self.model = ions, model
        self.R, self.centre = np.asarray(R, float), np.asarray(centre, float)
        self.sym, self.lab, self.bonds, self.colours, self.X3, self.off = [], [], [], {}, [], []
        n = 0
        for k, ion in enumerate(ions):
            self.off.append(n)
            self.sym += ion.sym
            self.lab += ion.lab
            loc = {s: m for m, s in enumerate(ion.idx)}
            self.bonds += [(n + loc[i], n + loc[j]) for i, j in ion.cr.bonds[ion.kind]]
            if k not in model:
                for m, s in enumerate(ion.sym):
                    self.colours[n + m] = draw2d.lighten(draw2d.element_colour(s), LIGHT)
            self.X3.append(ion.X)
            n += len(ion.sym)
        self.X3 = np.vstack(self.X3)
        self.V = self.view(self.X3)

    def view(self, P):
        return (np.atleast_2d(P) - self.centre) @ self.R.T

    def i(self, k, label):
        """Atom index in the scene of `label` in ion k."""
        return self.off[k] + self.ions[k].lab.index(label)

    def box(self, pad=0.0):
        r = np.array([style.BALL_RADII.get(s, 0.35) for s in self.sym])
        lo = (self.V[:, :2] - r[:, None]).min(0) - pad
        hi = (self.V[:, :2] + r[:, None]).max(0) + pad
        return np.array([lo[0], hi[0], lo[1], hi[1]])


def draw_scene(ax, sc, hbonds, hb_colours, label_atoms=(), hb_labels=(), extra_items=(), extra_segments=(),
               lines=()):
    """hbonds (d, h, a) scene indices; hb_labels: list of (triple, text, colour) placed next to the line;
    lines: extra (p0, p1, kw) drawn in view coordinates (also label obstacles)."""
    segs = list(extra_segments)
    for p0, p1, kw in lines:
        ax.add_line(Line2D([p0[0], p1[0]], [p0[1], p1[1]], **kw))
        segs.append((p0, p1, kw.get('lw', 0.8) / 2 + 0.4))
    items = list(extra_items)
    for (d, h, a), text, col in hb_labels:
        p0, p1 = sc.V[h, :2], sc.V[a, :2]
        mid = (p0 + p1) / 2
        dv = p1 - p0
        perp = np.array([-dv[1], dv[0]])
        items.append(dict(text=text, anchor=mid, radius_pt=0.8, prefer=perp if perp[1] >= 0 else -perp,
                          priority=0, kw=dict(fontsize=style.SMALL_PT, color=col), key=('hbl', d, h, a),
                          max_ring=6))
    d = draw2d.draw_molecule(ax, sc.sym, sc.V, sc.lab, bonds=sc.bonds, hbonds=hbonds, label_atoms=list(label_atoms),
                             label_h=True, hbond_distances=False, hbond_colours=hb_colours, fit=False,
                             strict=True, atom_colours=sc.colours, extra_items=items, extra_segments=segs,
                             hbond_lw=1.0)
    assert not d.dropped, d.dropped_labels()
    return d


# ---------------------------------------------------------------------------------------------------------------
# the three panels
# ---------------------------------------------------------------------------------------------------------------
def panel_a(cr):
    """Asymmetric unit as deposited."""
    P, Cc = cr.ion('P', 1, (0, 0, 0)), cr.ion('C', 1, (0, 0, 0))
    ions = [P, Cc]
    # closest inter-ion contact (heavy atoms)
    best = None
    for a, pa in zip(P.lab, P.X):
        for b, pb in zip(Cc.lab, Cc.X):
            if a[0] == 'H' or b[0] == 'H':
                continue
            dd = float(np.linalg.norm(pa - pb))
            if best is None or dd < best[0]:
                best = (dd, a, b)
    dist, a_lab, b_lab = best
    assert {a_lab, b_lab} == {'O4', 'N5'} and C.fact('xray.asu_closest_interion_contact') == 'O4...N5'
    C.check_equal(dist, 'xray.asu_closest_interion_distance', 3, 'ASU closest contact')
    assert not nh_o(ions, acceptors=('O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7'), max_ha=2.9), 'ASU is H-bonded?'
    X = np.vstack([P.X, Cc.X])
    sym = P.sym + Cc.sym
    _, R, c = geo.best_view(sym, X, hbonds=[])
    # picrate on the left, cation on the right (as in every other figure)
    v = (np.vstack([P.X.mean(0), Cc.X.mean(0)]) - c) @ R.T
    th = np.arctan2(v[1, 1] - v[0, 1], v[1, 0] - v[0, 0])
    Rz = np.array([[np.cos(-th), -np.sin(-th), 0], [np.sin(-th), np.cos(-th), 0], [0, 0, 1]])
    sc = Scene(ions, model={0, 1}, R=Rz @ R, centre=c)
    trip = (sc.i(1, 'N5'), sc.i(1, 'N5'), sc.i(0, 'O4'))
    lab = [sc.i(0, n) for n in ('O1', 'O4')] + [sc.i(1, n) for n in ('N4', 'N5', 'N6', 'H4N', 'H6N')]

    def draw(ax):
        return draw_scene(ax, sc, [trip], {trip: GREY}, label_atoms=lab,
                          hb_labels=[(trip, f'{dist:.3f} Å', '0.35')])
    return sc, draw, dict(dist=dist)


def model_pair(cr):
    P0 = cr.ion('P', 1, (0, 0, 0))
    sym_c = str(C.fact('model.symmetry_cation')).split()[0]           # '4_666'
    op, t = int(sym_c.split('_')[0]), [int(ch) - 5 for ch in sym_c.split('_')[1]]
    C0 = cr.ion('C', op, t)
    assert str(C.fact('model.symmetry_picrate')).split()[0] == P0.code
    # the highlighted pair is the calculated cluster (atom order of the case file)
    g = C.geometries()
    lab = g['labels']
    X = np.array([P0.pos(l) if i < 18 else C0.pos(l) for i, l in enumerate(lab)])
    _, rmsd = geo.align(X, g['xray'])
    assert rmsd < 0.01, rmsd
    return P0, C0, rmsd


def build_chain(cr, P0, C0):
    """P-1 <- C0 -> P0 <- C1 -> P1 (arrows: N-H donation)."""
    pool = cr.images(-2.6, 3.6)
    cats = [o for o in pool if o.kind == 'C']
    pics = [o for o in pool if o.kind == 'P']

    def accepts(cat, which):          # picrate receiving cat's N4-H4N...O1 ('N4') or N6-H6N...O6 ('N6')
        for p in pics:
            hb = nh_o([cat, p])
            for x in hb:
                if x['d'] == which and x['a'] == ('O1' if which == 'N4' else 'O6'):
                    return p, x
        raise AssertionError(f'no acceptor for {cat} {which}')

    def donor_to(pic, which):         # cation donating to pic's O1 ('N4') or O6 ('N6')
        for c in cats:
            hb = nh_o([c, pic])
            for x in hb:
                if x['d'] == which and x['a'] == ('O1' if which == 'N4' else 'O6'):
                    return c, x
        raise AssertionError(f'no donor for {pic} {which}')

    p0, hb1 = accepts(C0, 'N4')
    assert p0.same(P0)
    pm1, hb2 = accepts(C0, 'N6')
    c1, hb3 = donor_to(P0, 'N6')
    p1, hb4 = accepts(c1, 'N4')
    C.check_equal(hb1['da'], 'xray.hb01_N4_H4N_O1_da', 3)
    C.check_equal(hb2['da'], 'xray.hb02_N6_H6N_O6_da', 3)
    C.check_equal(hb3['da'], 'xray.hb02_N6_H6N_O6_da', 3)
    # the second picrate of the modelled cation is the 2_1-screw image 2_746 of the modelled picrate
    assert pm1.code == '2_746', pm1.code
    return [p1, c1, P0, C0, pm1]


def stacking_partner(cr, C0):
    """The cation whose triazole stacks on that of C0 (across an inversion centre); checked against facts."""
    best = None
    for o in cr.images(-2.6, 3.6):
        if o.kind != 'C' or o.same(C0):
            continue
        dd = float(np.linalg.norm(triazole_centroid(o) - triazole_centroid(C0)))
        if best is None or dd < best[0]:
            best = (dd, o)
    cg, partner = best
    C.check_equal(cg, 'xray.stacking_triazole_cg_cg', 3, 'triazole stacking centroid distance')
    perp = abs(float((triazole_centroid(partner) - triazole_centroid(C0)) @ triazole_normal(C0)))
    C.check_equal(perp, 'xray.stacking_triazole_perp', 3, 'triazole stacking interplanar separation')
    # inversion-related: the partner's triazole normal is parallel to C0's and the centroid of the pair is a
    # centre of symmetry (operator 3 = -x,-y,-z plus a lattice translation)
    mid = (partner.F.mean(0) + C0.F.mean(0)) / 2
    assert np.allclose(2 * mid, np.round(2 * mid), atol=1e-6), 'stacked pair is not centrosymmetric'
    return partner, cg, perp


def panel_d(cr, C0, partner):
    """The stacked cation pair seen nearly edge-on: dotted centroid line and the interplanar separation."""
    g0, g1 = triazole_centroid(C0), triazole_centroid(partner)
    n = triazole_normal(C0)
    if (g1 - g0) @ n < 0:
        n = -n
    slip = (g1 - g0) - ((g1 - g0) @ n) * n
    ex = slip / np.linalg.norm(slip)                       # screen x: the slip between the two rings
    ey = n                                                 # screen y: the ring normal
    ez = np.cross(ex, ey)
    tilt = np.radians(STACK_TILT)                          # rotate about screen x so the rings read as rings
    ey2 = np.cos(tilt) * ey + np.sin(tilt) * ez
    ez2 = np.cross(ex, ey2)
    R = np.array([ex, ey2, ez2])
    ions = [C0, partner]
    sc = Scene(ions, {0}, R, (g0 + g1) / 2)
    G0, G1 = sc.view(g0)[0, :2], sc.view(g1)[0, :2]
    perp = float((g1 - g0) @ n)
    cg = float(np.linalg.norm(g1 - g0))
    # interplanar separation: an arrow along the ring normal BESIDE the rings (right of every atom), from ring
    # plane 0 to ring plane 1, with short ticks lying in the two planes (along the slip direction = screen x)
    x_side = float(sc.V[:, 0].max()) + 0.9 - float(G0[0])
    Q0 = sc.view(g0 + x_side * ex)[0, :2]
    Q1 = sc.view(g0 + x_side * ex + perp * n)[0, :2]
    assert abs(Q0[0] - Q1[0]) < 1e-6                       # the normal projects onto screen y
    tick = 0.45
    lines = [(G0, G1, dict(color='black', lw=1.1, ls=(0, (0.1, 1.9)), dash_capstyle='round', zorder=80))]
    for Q in (Q0, Q1):
        lines.append((Q - [tick, 0], Q + [tick, 0], dict(color='0.2', lw=0.6, zorder=80)))
    lines.append((Q0, Q1, dict(color='0.2', lw=0.6, zorder=80)))          # the arrow shaft (label obstacle)
    mid = (G0 + G1) / 2
    dv = G1 - G0
    nrm = np.array([-dv[1], dv[0]]) / np.linalg.norm(dv)
    items = [dict(text=f'{perp:.3f} Å', anchor=(Q0 + Q1) / 2, radius_pt=1.2, prefer=(1, 0), priority=0,
                  key='perp', max_ring=3, leader=False, kw=dict(fontsize=style.SMALL_PT, color='0.2'))]
    radii = np.array([style.BALL_RADII.get(x, 0.35) for x in sc.sym])

    def draw(ax):
        from matplotlib.patches import FancyArrowPatch
        import matplotlib.patheffects as pe
        ax.add_patch(FancyArrowPatch(tuple(Q0), tuple(Q1), arrowstyle='<|-|>,head_length=2.4,head_width=1.1',
                                     mutation_scale=1.0, lw=0.6, color='0.2', shrinkA=0, shrinkB=0, zorder=81))
        for P_ in (G0, G1):
            ax.plot([P_[0]], [P_[1]], 'o', ms=2.2, color='black', zorder=92)     # above the leader's halo
        d = draw_scene(ax, sc, [], {}, extra_items=items, lines=lines)
        # QA round 4, r4-m4 / figures_r2 #3: the centroid distance is given in the key, not as a label; placed
        # anywhere near the rings it read as a label of N5 or as a bond to the centroid dot
        return d
    return sc, draw, dict(cg=cg, perp=perp)


def chain_view(ions, P0, C0):
    """Screen x along the chain (b, sign so the modelled picrate is left of its cation); rotation about b chosen
    to minimise projected overlaps plus the picture height (1 overlap unit per Å)."""
    X = np.vstack([o.X for o in ions])
    sym = sum((o.sym for o in ions), [])
    ex = np.array([0.0, 1.0, 0.0])
    if (C0.X.mean(0) - P0.X.mean(0)) @ ex < 0:
        ex = -ex
    heavy = np.array([s_ != 'H' for s_ in sym])
    Xh = X[heavy] - X.mean(0)
    u = np.array([1.0, 0.0, 0.0])
    w = np.cross(ex, u)
    w /= np.linalg.norm(w)
    u = np.cross(w, ex)
    best = None
    for phi in np.radians(np.arange(0, 360, 2)):
        vd = np.cos(phi) * u + np.sin(phi) * w          # towards the viewer
        up = np.cross(vd, ex)
        R = np.array([ex, up, vd])
        P = Xh @ R.T
        d = np.linalg.norm(P[:, None, :2] - P[None, :, :2], axis=2)
        np.fill_diagonal(d, 9)
        pen = float(np.maximum(0.95 - d, 0).sum())
        pz = ((P0.pos('N3') - P0.pos('N1')) @ up < 0) * 0.5     # N3 nitro above N1, as in Fig. 3
        h = P[:, 1].max() - P[:, 1].min()
        score = pen + pz + 1.0 * h
        if best is None or score < best[0]:
            best = (score, R, phi)
    return best[1], X.mean(0), np.degrees(best[2])


def panel_c(cr, P0, C0):
    ions = build_chain(cr, P0, C0)
    R, c, phi = chain_view(ions, P0, C0)
    # QA round 2, r2-m1: the chain is drawn with b vertical (up), as a column beside (a) and (b), so that the
    # figure has no empty block (the horizontal chain used 60 % of the width and left an 82 x 50 mm hole)
    R = np.array([-R[1], R[0], R[2]])
    k_p0, k_c0 = 2, 3
    sc = Scene(ions, {k_p0, k_c0}, R, c)
    hbs, cols, labs = [], {}, []
    for x in nh_o(ions):
        t = (sc.i(x['don'], x['d']), sc.i(x['don'], x['h']), sc.i(x['acc'], x['a']))
        hbs.append(t)
        cols[t] = TEAL if x['d'] == 'N4' else GREY
        if x['don'] == k_c0:                    # the two bonds donated by the modelled cation
            labs.append((x['d'], x['da']))
    assert len(hbs) == 4, hbs
    lab = [sc.i(k_p0, n) for n in ('O1', 'O6', 'O7')] + \
          [sc.i(k_c0, n) for n in ('N4', 'N5', 'N6', 'H4N', 'H6N')]

    da = dict(labs)
    C.check_equal(da['N4'], 'xray.hb01_N4_H4N_O1_da', 3)
    C.check_equal(da['N6'], 'xray.hb02_N6_H6N_O6_da', 3)

    # axis arrow for b (QA round 1, m4): 2 Å long, along b (vertical), in the free corner of the chain picture
    bdir = sc.R @ (cr.M[:, 1] / np.linalg.norm(cr.M[:, 1]))
    assert abs(abs(bdir[1]) - 1) < 1e-9, bdir                 # screen y is along b by construction
    box = sc.box()
    best = None
    for cx in (box[0] + 0.5, box[1] - 0.5):
        for cy in (box[2] + 0.4, box[3] - 2.4):
            q0 = np.array([cx, cy + (2.0 if bdir[1] < 0 else 0.0)])
            q1 = q0 + np.array([0.0, 2.0 * np.sign(bdir[1])])
            dmin = min(float(np.min(np.linalg.norm(sc.V[:, :2] - (q0 + t * (q1 - q0)), axis=1)))
                       for t in np.linspace(0, 1, 11))
            if best is None or dmin > best[0]:
                best = (dmin, q0, q1)
    dmin, p0, p1 = best
    assert dmin > 1.5, f'b arrow too close to an atom ({dmin:.2f} Å)'
    b_item = dict(text='$b$', anchor=p1, radius_pt=1.5, prefer=(0, np.sign(bdir[1])), priority=0, key='b',
                  max_ring=2, leader=False, kw=dict(fontsize=style.BASE_PT, color=CELL_COLOUR))

    def draw(ax):
        from matplotlib.patches import FancyArrowPatch
        ax.add_patch(FancyArrowPatch(tuple(p0), tuple(p1), arrowstyle='-|>,head_length=3.0,head_width=1.5',
                                     mutation_scale=1.0, lw=0.8, color=CELL_COLOUR, shrinkA=0, shrinkB=0,
                                     zorder=81))
        return draw_scene(ax, sc, hbs, cols, label_atoms=lab, extra_items=[b_item],
                          extra_segments=[(p0, p1, 1.5)])
    return sc, draw, dict(phi=phi, ions=ions, k=(k_p0, k_c0), da=da)


# ---------------------------------------------------------------------------------------------------------------
def key_handles(da, cg):
    """Key for (b) and (c): the balls show C, N and O in full (modelled pair) and in paler (other ions) colours,
    as drawn (QA round 1, m3: the grey 'other ions' swatch did not match the pale element colours of (c))."""
    from figkit.draw2d import element_colour, lighten, darken

    def balls(light):
        out = []
        for e in ('C', 'N', 'O'):
            c = element_colour(e)
            c = lighten(c, LIGHT) if light else c
            out.append(Line2D([], [], ls='none', marker='o', ms=4.6, mfc=c, mec=darken(c, 0.55), mew=0.45))
        return tuple(out)
    return [(balls(False), 'Modelled ion pair'), (balls(True), 'Other ions'),
            (Line2D([], [], color=TEAL, lw=1.0, ls=(0, (2.2, 1.4))), f'N4–H4N···O1 (N···O {da["N4"]:.3f} Å)'),
            (Line2D([], [], color=GREY, lw=1.0, ls=(0, (2.2, 1.4))), f'N6–H6N···O6 (N···O {da["N6"]:.3f} Å)'),
            (Line2D([], [], color='black', lw=1.1, ls=(0, (0.1, 1.9)), dash_capstyle='round'),
             f'Triazole stacking, centroids {cg:.3f} Å apart')]


def make_key(axk, hk, loc='upper left'):
    from matplotlib.legend_handler import HandlerTuple
    leg = axk.legend(handles=[h for h, t in hk], labels=[t for h, t in hk], loc=loc, frameon=False,
                     fontsize=style.TICK_PT, handlelength=2.6, handletextpad=0.6, borderaxespad=0.0,
                     labelspacing=0.55, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)},
                     title='Key to (b) and (c)', title_fontsize=style.TICK_PT, alignment='left')
    leg.get_title().set_fontweight('bold')
    return leg


def key_size_mm(hk):
    """(width, height) in mm of the key as make_key draws it."""
    f = plt.figure(figsize=(C.mm(120.0), C.mm(80.0)))
    axk = f.add_axes([0, 0, 1, 1])
    C.hide_frame(axk)
    leg = make_key(axk, hk)
    f.canvas.draw()
    bb = leg.get_window_extent(f.canvas.get_renderer())
    plt.close(f)
    return bb.width / f.dpi / style.MM + 0.3, bb.height / f.dpi / style.MM + 0.3


def main():
    style.apply()
    cr = Crystal()
    P0, C0, rmsd = model_pair(cr)
    print(f'model pair {P0.code} + {C0.code}: RMSD vs btpc_cluster.xyz {rmsd:.2e} Å')
    partner, cg, perp = stacking_partner(cr, C0)
    sa, fa, ia = panel_a(cr)
    sc_, fc, ic = panel_c(cr, P0, C0)
    sd, fd, idd = panel_d(cr, C0, partner)
    print(f'(a) O4...N5 {ia["dist"]:.3f} Å; (c) chain {ic["ions"]}, rotation {ic["phi"]:.0f} deg; '
          f'(b) stacking partner {partner.code}: {cg:.3f} / {perp:.3f} Å')
    # layout (QA round 2, r2-m1): left column (a) over [(b) | key], right column the chain (c) with b vertical;
    # scales chosen so that the two columns fill the width and end at the same height
    hk = key_handles(ic['da'], cg)
    KEY_W, KEY_H = key_size_mm(hk)
    W_FIT = WIDTH_MM - 2 * 25.4 * float(plt.rcParams['savefig.pad_inches']) - 0.1   # savefig adds the pad
    s1, s2 = 5.8, 4.5
    for _ in range(40):
        ba, bd = C.measure_drawing(fa, sa.box(), s1), C.measure_drawing(fd, sd.box(), s1)
        bc = C.measure_drawing(fc, sc_.box(), s2)
        wa, ha = (ba[1] - ba[0]) * s1 + 1.0, (ba[3] - ba[2]) * s1 + 1.0
        wd, hd = (bd[1] - bd[0]) * s1 + 1.0, (bd[3] - bd[2]) * s1 + 1.0
        wc, hc = (bc[1] - bc[0]) * s2 + 1.0, (bc[3] - bc[2]) * s2 + 1.0
        L = W_FIT - GAP_MM - wc                                       # width of the left column
        h_left = HEAD_MM + ha + ROW_GAP_MM + HEAD_MM + max(hd, KEY_H)
        h_right = HEAD_MM + hc
        ok = wa <= L + 1e-6 and wd + GAP_MM + KEY_W <= L + 1e-6 and abs(h_left - h_right) < 0.6
        if ok:
            break
        s1 = min(S1_MAX, s1 * min(L / wa, (L - GAP_MM - KEY_W) / wd) * 0.998)
        s2 = min(s1, s2 * (h_left - HEAD_MM - 1.0) / (hc - 1.0))
    assert ok, (s1, s2, wa, wd, L, h_left, h_right)
    H = max(h_left, h_right)
    print(f'scales (a, b) {s1:.2f}, (c) {s2:.2f} mm per Å; left column {L:.1f} mm; height {H:.1f} mm '
          f'(left {h_left:.1f}, right {h_right:.1f})')
    fig = plt.figure(figsize=(C.mm(WIDTH_MM), C.mm(H)))
    # (a) top left
    ax, w, h = C.drawing_axes(fig, 0.0, H - HEAD_MM - ha, ba, s1)
    fa(ax)
    C.header(fig, 0.0, H - HEAD_MM + 1.4, 'a', 'Asymmetric unit as deposited')
    # (b) below (a), key to its right, top-aligned with the drawing of (b)
    y_b = H - HEAD_MM - ha - ROW_GAP_MM
    ax, w, h = C.drawing_axes(fig, 0.0, y_b - HEAD_MM - hd, bd, s1)
    fd(ax)
    C.header(fig, 0.0, y_b - HEAD_MM + 1.4, 'b', 'Triazole stacking')
    x_key = L - KEY_W
    axk = C.mm_box(fig, x_key, y_b - HEAD_MM - KEY_H, KEY_W, KEY_H)
    C.hide_frame(axk)
    make_key(axk, hk, loc='upper left')
    # (c) right column
    x_c = W_FIT - wc
    ax, w, h = C.drawing_axes(fig, x_c, H - HEAD_MM - hc, bc, s2)
    fc(ax)
    C.header(fig, x_c, H - HEAD_MM + 1.4, 'c', 'N–H···O chain along $b$')
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
