"""Fig. 4 fig_xray_comparison: optimised versus X-ray geometry.

(a, b) stick overlays of the crystal ion pair (black) with the gas-phase / IEF-PCM (water) structure after a
Kabsch superposition on the 16 picrate heavy atoms, drawn in the paper's common orientation (X-ray camera);
(c) parity plot of the 26 heavy-atom bond lengths; (d) inter-ring angle and nitro twists.
Every printed number is computed from the XYZ files and asserted against facts.json.

    python scripts/figures/fig_xray_comparison.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from figkit import style, draw2d, geometry as geo, plots
from figkit.labels import LabelPlacer, _segments_cross
import _struct_common as C

STEM = 'fig_xray_comparison'
WIDTH_MM = 175.0
GAP_MM = 4.0
HEAD_MM = 5.5            # header strip above each row
LABEL_ATOMS = ['O1', 'N1', 'N2', 'N3']          # picrate: the superposed geometries coincide, one label each
CATION_LABELS = ['N4', 'N5', 'N6']              # cation: labelled on both geometries, in each geometry's colour
# QA round 1, M2: the X-ray structure is drawn ON TOP with a thin line and its own atom dots, the model underneath
# with a wide line, so an X-ray bond lying under a model bond can no longer peek out as a tapered sliver and every
# X-ray atom stays visible; cation labels sit outside both rings (ring interiors are label obstacles)
LW_MODEL, LW_XRAY = 2.0, 0.8
HALO_MIN_ANGLE = 20.0    # QA round 4, r4-m2: the X-ray halo is kept only where it crosses a model bond steeply
LEADER_END_PT = 2.5      # a leader may not end closer than this to a bond it does not belong to
LEADER_CLEAR_PT = 1.5    # a leader keeps this far from every bond it does not cross
CROSS_MIN_ANGLE = 50.0   # the one crossing allowed for an atom inside the other cation's ring: steep ...
CROSS_ATOM_CLEAR_A = 0.25  # ... and at least this far (Å) from every atom


def _seg_point_dist(S, P):
    """Distance from points P (N, 2) to segments S (M, 4) -> (N, M)."""
    a, b = S[None, :, 0:2], S[None, :, 2:4]
    d = b - a
    t = np.clip(((P[:, None, :] - a) * d).sum(-1) / np.maximum((d ** 2).sum(-1), 1e-12), 0, 1)
    q = a + t[..., None] * d
    return np.hypot(*(P[:, None, :] - q).transpose(2, 0, 1))


class _StrictLeaders(LabelPlacer):
    """LabelPlacer for the cation labels (QA round 4, r4-m2): a leader may not end on a bond of either structure and
    may not run along one. Where an atom of one cation lies inside the triazole ring of the other (N6 of the gas-phase
    cation in (a), N6 of the X-ray cation in (b)), no leader can leave it without crossing one bond; such a leader
    may cross exactly one bond of the other structure, steeply (>= CROSS_MIN_ANGLE) and away from every atom, and is
    drawn with a white halo. `avoid` holds the bonds of both structures in data coordinates; bonds that meet at the
    label's own atom are exempt at that atom."""
    avoid = ()
    atoms = ()

    def _evaluate(self, it, w, h, placed_boxes, placed_leaders, discs, segs, bnd, static_boxes):
        ok, cost, Cc, boxes, ring_id, a, r = super()._evaluate(it, w, h, placed_boxes, placed_leaders, discs, segs,
                                                                bnd, static_boxes)
        lead = ring_id >= self.leader_from
        if lead.any() and len(self.avoid):
            T = self.ax.transData.transform
            S = np.array([np.r_[T(p0), T(p1)] for p0, p1 in self.avoid])
            own = (np.hypot(S[:, 0] - a[0], S[:, 1] - a[1]) < 2.0) | (np.hypot(S[:, 2] - a[0], S[:, 3] - a[1]) < 2.0)
            S = S[~own]
            L = self._leaders(a, r, Cc, boxes)
            cr = _segments_cross(L, S)                                   # (Ncand, Nseg)
            ncross = cr.sum(1)
            near_end = (_seg_point_dist(S, L[:, 2:4]) < LEADER_END_PT * self._ppt()).any(1)
            # a leader that does not cross a bond must also keep clear of it along its whole length (the N5 leader of
            # the old (a) started on the X-ray N4-N5 line): segment-segment distance from the four endpoint distances
            dLS = np.minimum(np.minimum(_seg_point_dist(S, L[:, 0:2]), _seg_point_dist(S, L[:, 2:4])),
                             np.minimum(_seg_point_dist(L, S[:, 0:2]).T, _seg_point_dist(L, S[:, 2:4]).T))
            graze = ((dLS < LEADER_CLEAR_PT * self._ppt()) & ~cr).any(1)
            bad = near_end | graze | (ncross >= 2)
            one = np.nonzero(ncross == 1)[0]
            if len(one):
                A = np.array([T(x) for x in self.atoms])
                for k in one:
                    j = int(np.nonzero(cr[k])[0][0])
                    u = L[k, 2:4] - L[k, 0:2]; v = S[j, 2:4] - S[j, 0:2]
                    ang = np.degrees(np.arccos(min(1.0, abs(float(u @ v)) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-12))))
                    # intersection point of the two lines
                    M = np.array([u, -v]).T
                    t = np.linalg.solve(M, S[j, 0:2] - L[k, 0:2])[0]
                    X = L[k, 0:2] + t * u
                    dmin = np.hypot(A[:, 0] - X[0], A[:, 1] - X[1]).min()
                    if ang < CROSS_MIN_ANGLE or dmin < CROSS_ATOM_CLEAR_A * self._scale():
                        bad[k] = True
            ok &= ~(lead & bad)
            cost = cost + 3.0 * (lead & (ncross == 1))
        return ok, cost, Cc, boxes, ring_id, a, r


def _halo_only_at_crossings(ax, geoms, m, symbols):
    """Remove the white halo of the X-ray layer where it runs nearly parallel over a model bond: there it cut the
    wide model stick into a tapered sliver that read as a stereo wedge (QA round 4, r4-m2)."""
    Xm = np.asarray(geoms[style.MODEL_LABELS[m]])
    sym_bonds = [(Xm[i, :2], Xm[j, :2]) for i, j in geo.detect_bonds(symbols, Xm)
                 if symbols[i] != 'H' and symbols[j] != 'H']
    removed = 0
    for ln in list(ax.lines):
        if ln.get_color() != 'white' or ln.get_zorder() < 12:          # the X-ray layer is drawn second (zorder 12+)
            continue
        x, y = ln.get_data()
        p0, p1 = np.array([x[0], y[0]]), np.array([x[1], y[1]])
        u = (p1 - p0) / np.linalg.norm(p1 - p0)
        for q0, q1 in sym_bonds:
            v = (q1 - q0) / np.linalg.norm(q1 - q0)
            ang = np.degrees(np.arccos(min(1.0, abs(float(u @ v)))))
            mid = (p0 + p1) / 2
            dist = _seg_point_dist(np.r_[q0, q1][None], mid[None])[0, 0]
            if ang < HALO_MIN_ANGLE and dist < 0.5 * np.linalg.norm(p1 - p0):
                ln.remove()
                removed += 1
                break
    return removed


def overlay_geoms(g):
    out, rmsd = {}, {}
    Xx = C.oriented('xray')
    for m in C.MODELS:
        al, r = geo.align(g[m], g['xray'], idx=C.PIC_HEAVY)
        C.check_equal(r, f'{m}.kabsch_rmsd_picrate_heavy', 3, f'picrate superposition RMSD ({m})')
        rmsd[m] = r
        out[m] = {style.MODEL_LABELS[m]: C.oriented('xray', al), 'X-ray': Xx}      # X-ray layer on top
    return out, rmsd


def draw_overlay_panel(ax, g, geoms, m):
    cols = {'X-ray': style.MODEL_COLOURS['xray'], style.MODEL_LABELS[m]: style.MODEL_COLOURS[m]}
    d = draw2d.draw_overlay(ax, g['symbols'], geoms, cols, include_h=False, labels=g['labels'],
                            label_atoms=LABEL_ATOMS, legend=False, fit=False, strict=True,
                            lw={style.MODEL_LABELS[m]: LW_MODEL, 'X-ray': LW_XRAY})
    assert not d.dropped
    nh = _halo_only_at_crossings(ax, geoms, m, g['symbols'])
    print(f'{m}: X-ray halo removed on {nh} bonds lying nearly parallel over a model bond')
    # X-ray atoms as slightly larger dots on top of everything, so an atom under a model bond stays visible
    Xx = np.asarray(geoms['X-ray'])
    from matplotlib.patches import Circle
    for i, x in enumerate(g['symbols']):
        if x != 'H':
            ax.add_patch(Circle(Xx[i, :2], 0.085, facecolor=cols['X-ray'], edgecolor='white', lw=0.35,
                                zorder=40))
    # the cations do not coincide: label N4, N5, N6 on each geometry, coloured like its sticks, so every vertex
    # can be assigned to its structure
    sym = g['symbols']
    heavy = [i for i, x in enumerate(sym) if x != 'H']
    pl = _StrictLeaders(ax, name='cation labels', n_rings=13, leader_from=1, n_dirs=48)
    pl.avoid = [(np.asarray(X)[i, :2], np.asarray(X)[j, :2]) for X in geoms.values()
                for i, j in geo.detect_bonds(g['symbols'], np.asarray(X)) if i in heavy and j in heavy]
    pl.atoms = [np.asarray(X)[i, :2] for X in geoms.values() for i in heavy]
    rings = [[C.I[a] for a in ('N4', 'N5', 'N6', 'C7', 'C12')],
             [C.I[a] for a in ('C7', 'C8', 'C9', 'C10', 'C11', 'C12')]]
    for name, X in geoms.items():
        X = np.asarray(X)
        for i, j in geo.detect_bonds(sym, X):
            if i in heavy and j in heavy:
                pl.add_segment(X[i, :2], X[j, :2], hw_pt=LW_MODEL / 2 + 0.8)
        for i in heavy:
            pl.add_disc(X[i, :2], r_data=0.09, r_pt=0.3)
        for rg in rings:                      # ring interiors (projected): no label inside either ring.  Spokes
            P = X[rg, :2]                     # from the centre to every atom and bond midpoint block label boxes
            c = P.mean(0)                     # but, unlike discs, still let a leader leave an atom inside a ring
            for q in list(P) + [(P[k] + P[(k + 1) % len(P)]) / 2 for k in range(len(P))]:
                pl.add_segment(c, q, hw_pt=1.0)
    cen_all = np.mean([np.asarray(X)[C.CAT9].mean(0) for X in geoms.values()], axis=0)
    for it in d.placed:
        pl.add_text_obstacle(it.artist)
    for name, X in geoms.items():
        X = np.asarray(X)
        col = draw2d.darken(cols[name], 0.12)
        for lab in CATION_LABELS:
            i = C.I[lab]
            pl.add_item(lab, X[i, :2], radius_data=0.12, radius_pt=0.3, prefer=X[i, :2] - cen_all[:2], priority=1,
                        kw=dict(fontsize=style.SMALL_PT, color=col, leader_color=col, leader_lw=0.5),
                        key=(name, lab), max_ring=13)
    placed, dropped = pl.place()
    assert not dropped, [it.text for it in dropped]
    # leaders that cross a bond of the other structure (see _StrictLeaders) get a white halo, so the crossing reads
    # as a pointer passing over the stick, not as a bond; no leader crosses more than one bond or ends on one
    import matplotlib.patheffects as pe
    segs = np.array([np.r_[p0, p1] for p0, p1 in pl.avoid])
    n_cross = 0
    for it in placed:
        ln = it.leader_artist
        if ln is None:
            continue
        x, y = ln.get_data()
        Lk = np.array([[x[0], y[0], x[-1], y[-1]]])
        a_ = np.array(it.anchor)
        own = (np.hypot(segs[:, 0] - a_[0], segs[:, 1] - a_[1]) < 1e-6) | (np.hypot(segs[:, 2] - a_[0], segs[:, 3] - a_[1]) < 1e-6)
        c_ = int(_segments_cross(Lk, segs[~own]).sum())
        assert c_ <= 1, (it.text, it.key, c_)
        if c_:
            n_cross += 1
            ln.set_path_effects([pe.Stroke(linewidth=ln.get_linewidth() + 1.4, foreground='white'), pe.Normal()])
            ln.set_zorder(45)
    print(f'{m}: cation labels placed; {n_cross} leader(s) cross one bond of the other structure (haloed)')
    return d


def bars_panel(ax, g):
    # nitro groups named by the ring carbon they sit on, as in Table 2 (QA round 1, m6: N1 and N3 are both ortho
    # to O1, so 'O1 side' did not distinguish them)
    cats = ['Inter-ring', 'N1 nitro\n(ortho, on C2)', 'N2 nitro\n(para, on C4)', 'N3 nitro\n(ortho, on C6)']
    assert [C.NITRO[n][1] for n in ('N1', 'N2', 'N3')] == [C.I['C2'], C.I['C4'], C.I['C6']]
    series = {}
    for m in ('xray', 'gas', 'pcm'):
        X = g[m]
        vals = [C.interplanar(X)] + [C.nitro_twist(X, n) for n in ('N1', 'N2', 'N3')]
        C.check_equal(vals[0], f'{m}.interplanar_picrate_cation', 1)
        for n, v in zip(('N1', 'N2', 'N3'), vals[1:]):
            C.check_equal(v, f'{m}.nitro_twist_{n}', 1)
        series[m] = vals
    plots.grouped_bars(ax, cats, series, ylabel='Angle (°)', fmt='{:.1f}',
                       labels={'xray': 'X-ray', 'gas': style.MODEL_LABELS['gas'], 'pcm': style.MODEL_LABELS['pcm']},
                       legend=False, width=0.88)
    for t in ax.texts:                               # horizontal value labels (6 pt fits the widened bars)
        t.set_rotation(0)
    ax.set_ylim(0, 66)
    ax.tick_params(axis='x', length=0)
    return series


def parity_panel(ax, g):
    sym, labels = g['symbols'], g['labels']
    rows = plots.bond_table(sym, g['xray'], {'gas': g['gas'], 'pcm': g['pcm']}, labels, exclude_h=True)
    assert len(rows) == 26, len(rows)
    ref = np.array([r['ref'] for r in rows])
    mod = {m: np.array([r[m] for r in rows]) for m in C.MODELS}
    for m in C.MODELS:
        rmsd = float(np.sqrt(np.mean((mod[m] - ref) ** 2)))
        C.check_equal(rmsd, f'{m}.bond_rmsd_heavy', 4, f'heavy-atom bond RMSD ({m})')
    lim = (1.17, 1.48)
    assert ref.min() > lim[0] and max(v.max() for v in mod.values()) < lim[1]
    # one marker per model across the paper (QA round 2, r2-m6): gas filled triangle, IEF-PCM open square
    hs = plots.parity_plot(ax, ref, mod, band=0.02, lim=lim, rmsd_fmt='{:.4f}', legend_loc='upper left',
                           marker_kw={m: C.model_marker(m, 3.3) for m in C.MODELS})
    # legend with the two models only; the ±0.02 Å band is labelled directly along its lower edge
    ax.get_legend().remove()
    leg = ax.legend(handles=hs, loc='upper left', handletextpad=0.3, borderaxespad=0.4, frameon=False)
    fig = ax.figure
    fig.canvas.draw()
    xb = 1.285
    p0 = ax.transData.transform((xb, xb + 0.02))
    p1 = ax.transData.transform((xb + 0.01, xb + 0.01 + 0.02))
    ang = float(np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0])))
    ax.annotate('±0.02 Å', (xb, xb + 0.02), xytext=(-0.8, 0.8), textcoords='offset points', rotation=ang,
                rotation_mode='anchor', ha='center', va='bottom', fontsize=style.SMALL_PT, color='0.35')
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05))
    # outliers: the three largest |DFT - X-ray| bonds of each model (union), cf. *.bond_largest_dev_heavy
    top = set()
    for m in C.MODELS:
        dev = np.abs(mod[m] - ref)
        top |= set(np.argsort(dev)[::-1][:3].tolist())
    for m in C.MODELS:
        listed = [t.split()[0].replace('-', '–') for t in str(C.fact(f'{m}.bond_largest_dev_heavy')).split(';')][:3]
        mine = [rows[i]['name'] for i in np.argsort(np.abs(mod[m] - ref))[::-1][:3]]
        norm = lambda s: '–'.join(sorted(s.split('–')))
        assert sorted(map(norm, listed)) == sorted(map(norm, mine)), (listed, mine)
    fig.canvas.draw()
    pl = LabelPlacer(ax, name='parity labels', n_rings=5, clearance_pt=1.0)
    for m in C.MODELS:
        for x, y in zip(ref, mod[m]):
            pl.add_disc((x, y), r_pt=2.2)
    pl.add_segment((lim[0], lim[0]), (lim[1], lim[1]), hw_pt=0.5)
    for t in leg.get_texts():
        pl.add_text_obstacle(t)
    for i in sorted(top):
        r = rows[i]
        y = max(r['gas'], r['pcm']) if (r['gas'] - r['ref']) > 0 else min(r['gas'], r['pcm'])
        above = (r['gas'] - r['ref']) > 0
        name = '–'.join(sorted(r['name'].split('–'), key=lambda s: (s[0], int(''.join(c for c in s if c.isdigit())))))
        pl.add_item(name, (r['ref'], y), radius_pt=2.5, prefer=(-0.4, 1) if above else (0.4, -1), priority=1,
                    kw=dict(fontsize=style.SMALL_PT, color='0.15'), key=('bond', i), max_ring=4)
    placed, dropped = pl.place()
    assert not dropped, [d.text for d in dropped]
    return rows


def main():
    style.apply()
    g = C.geometries()
    sym = g['symbols']
    geoms, rmsd = overlay_geoms(g)

    # ---- row 1: overlays at one scale, cropped to what is drawn ----
    def overlay_boxes(s):
        boxes = []
        for m in C.MODELS:
            allxy = np.vstack([X[[i for i, x in enumerate(sym) if x != 'H']] for X in geoms[m].values()])
            b0 = np.array([allxy[:, 0].min(), allxy[:, 0].max(), allxy[:, 1].min(), allxy[:, 1].max()])
            boxes.append(C.measure_drawing(lambda ax, m=m: draw_overlay_panel(ax, g, geoms[m], m), b0, s))
        y0, y1 = min(b[2] for b in boxes), max(b[3] for b in boxes)
        return [np.array([b[0], b[1], y0, y1]) for b in boxes]

    # content width = WIDTH_MM minus the two savefig pads, so the saved figure is exactly the planned width
    # (QA round 1, m5: it came out 175.5 mm)
    W_FIT = WIDTH_MM - 2 * 25.4 * float(plt.rcParams['savefig.pad_inches']) - 0.15
    s = 6.0
    for _ in range(6):
        boxes = overlay_boxes(s)
        w = sum(b[1] - b[0] for b in boxes) * s + GAP_MM + 2 * 0.5 * 2
        if W_FIT - 0.25 <= w <= W_FIT:
            break
        s *= (W_FIT - 0.12) / w
    h1 = (boxes[0][3] - boxes[0][2]) * s + 1.0
    h2 = 74.0                       # row 2 plots incl. tick labels
    H = HEAD_MM + h1 + 3.0 + HEAD_MM + h2
    fig = plt.figure(figsize=(C.mm(WIDTH_MM), C.mm(H)))
    y_row1 = H - HEAD_MM - h1
    x = 0.0
    for k, (m, b) in enumerate(zip(C.MODELS, boxes)):
        ax, w, h = C.drawing_axes(fig, x, y_row1, b, s)
        draw_overlay_panel(ax, g, geoms[m], m)
        hl = [Line2D([], [], color=style.MODEL_COLOURS['xray'], lw=LW_XRAY, marker='o', ms=2.4, mec='white',
                     mew=0.3, label='X-ray'),
              Line2D([], [], color=style.MODEL_COLOURS[m], lw=LW_MODEL,
                     label=f'{style.MODEL_LABELS[m]}, picrate RMSD {rmsd[m]:.3f} Å')]
        C.header(fig, x, H - HEAD_MM + 1.6, 'ab'[k], handles=hl)
        x += w + GAP_MM

    # ---- row 2: parity plot and angle bars ----
    y2 = 9.5                                        # room for the x tick labels and axis labels
    hp = h2 - y2 - 1.0
    axc = C.mm_box(fig, 13.0, y2, hp, hp)
    rows = parity_panel(axc, g)
    xd = 13.0 + hp + 17.0
    axd = C.mm_box(fig, xd, y2, W_FIT - xd - 0.3, hp)
    bars_panel(axd, g)
    yh2 = y2 + hp + 2.2
    C.header(fig, 0.0, yh2, 'c')
    hd = [plt.Rectangle((0, 0), 1, 1, color=style.MODEL_COLOURS[m], label=lab)
          for m, lab in (('xray', 'X-ray'), ('gas', style.MODEL_LABELS['gas']), ('pcm', style.MODEL_LABELS['pcm']))]
    C.header(fig, xd - 12.0, yh2, 'd', handles=hd, handle_len=1.0)
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
