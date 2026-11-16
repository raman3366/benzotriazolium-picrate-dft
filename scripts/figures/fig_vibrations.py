"""Fig. 6 fig_vibrations: harmonic vibrational evidence (no intensities were computed).

(a) X-H stretching wavenumbers, gas phase and IEF-PCM (water), with the free-minus-bonded N-H red shifts;
(b, c) the N4-H4N stretch drawn as Cartesian displacement arrows on the structure in the paper's common
orientation; (d, e) the lowest mode of each model (gas 12 cm-1, IEF-PCM 3.2i cm-1), which is a libration of the
two ions: each ion's displacement is fitted by a rigid rotation + translation about its own centre of mass
(_vib.rigid_decomposition, the 12-dimensional space of `*.mode1_rigid_fraction`), and the rotation is drawn as a
hoop round the ion's rotation axis (true 3D projection; far half behind the molecule, near half in front with the
arrowhead giving the sense) with the axis as dash-dotted stubs.  The translations are under 5 % of each mode
(asserted) and are not drawn (QA round 1, M3: per-atom arrows of these modes were mostly arrowheads).
Normal modes are recomputed from hessian.npy (mass-weighted, translations/rotations projected) and checked
against freq.json (every wavenumber within 0.5 cm-1) before anything is drawn; printed wavenumbers are the
freq.json values and are asserted against facts.json.

    python scripts/figures/fig_vibrations.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from figkit import style
import _struct_common as C
import _vib as VB

STEM = 'fig_vibrations'
WIDTH_MM = 175.0
A_W = 55.0                 # column (a), incl. its tick labels
GAP_MM = 3.0
HEAD_MM = 5.0
HOOP_LW = 1.0              # rotation hoops of (d, e)
NEAR = ['N4', 'N5', 'N6', 'C7', 'C12', 'H4N', 'H6N', 'O1', 'C1', 'C2', 'C6', 'N3', 'O7', 'O6']
HB_ALL = {'gas': [C.HB, (C.I['C11'], C.I['H11'], C.I['O2'])], 'pcm': [C.HB]}


def wn(v, imag_ok=True):
    """Wavenumber text: nearest integer, imaginary as '3.2i'."""
    if v < 0 and imag_ok:
        return f'{abs(v):.1f}i'
    return f'{v:.0f}'


def mode_table():
    """Indices and freq.json wavenumbers of the modes shown, checked against facts.json."""
    T = {}
    for m in C.MODELS:
        md = VB.modes(m)
        fj = md['freq_json']
        xs = VB.xh_stretches(m)
        k_n4 = VB.find_mode(md, float(C.fact(f'{m}.nu_N4_H4N_stretch')))
        k_n6 = VB.find_mode(md, float(C.fact(f'{m}.nu_N6_H6N_stretch')))
        assert xs['H4N'][0] == k_n4 and xs['H6N'][0] == k_n6
        assert VB.stretch_fraction(md, k_n4, C.I['N4'], C.I['H4N'], m) > 0.5
        C.check_close(fj[k_n4], f'{m}.nu_N4_H4N_stretch', 1e-6)
        C.check_close(fj[k_n6], f'{m}.nu_N6_H6N_stretch', 1e-6)
        C.check_close(fj[k_n6] - fj[k_n4], f'diff.nh_redshift_{m}', 1e-6)
        k_low = 0
        C.check_close(fj[k_low], f'{m}.freq_lowest', 1e-6)
        T[m] = dict(n4=k_n4, n6=k_n6, low=k_low, h3=xs['H3'][0], h5=xs['H5'][0], benzo=VB.benzo_ch_modes(m),
                    fj=fj, maxdev=md['max_dev'])
        assert xs['H3'][2] > 0.8 and xs['H5'][2] > 0.8           # local picrate C-H stretches
        # the C-H wavenumbers against their facts.json keys (QA round 3, r3-m1; they replace the PLAN section 7
        # planning values of rounds 1-2): picrate C3-H3 and C5-H5 exactly, the four benzo modes at the
        # 1-decimal precision of the key, and their spread
        C.check_close(fj[xs['H3'][0]], f'{m}.nu_C3_H3_stretch', 1e-6)
        C.check_close(fj[xs['H5'][0]], f'{m}.nu_C5_H5_stretch', 1e-6)
        bz = [fj[k] for k in VB.benzo_ch_modes(m)]
        want = [float(v) for v in C.fact(f'{m}.nu_benzo_CH_stretches')]
        assert [f'{v:.1f}' for v in bz] == [f'{v:.1f}' for v in want], (m, bz, want)
        C.check_close(max(bz) - min(bz), f'{m}.benzo_ch_spread', 1e-6)
        print(f'{m}: max |recomputed - freq.json| = {md["max_dev"]:.3f} cm-1; N4-H4N {fj[k_n4]:.1f}, '
              f'N6-H6N {fj[k_n6]:.1f}, C3-H3 {fj[T[m]["h3"]]:.1f}, C5-H5 {fj[T[m]["h5"]]:.1f}, benzo '
              f'{", ".join(f"{fj[k]:.1f}" for k in T[m]["benzo"])}, lowest {fj[k_low]:.2f}')
    # largest IEF-PCM - gas shift of the six C-H stretches (benzo modes paired in ascending order)
    ch = {m: [T[m]['fj'][T[m]['h3']], T[m]['fj'][T[m]['h5']]] + [T[m]['fj'][k] for k in T[m]['benzo']]
          for m in C.MODELS}
    C.check_close(max(abs(p - g) for p, g in zip(ch['pcm'], ch['gas'])), 'diff.ch_stretch_max_shift', 1e-6)
    return T


def stretch_chart(ax, T):
    rows = [('N6–H6N\n(free)', 'n6'), ('N4–H4N\n(bonded)', 'n4'), ('C3–H3', 'h3'), ('C5–H5', 'h5'),
            ('Benzo C–H\n(4 modes)', 'benzo')]
    ny = len(rows)
    off = {'gas': 0.16, 'pcm': -0.16}
    mk = {m: C.model_marker(m, 3.9) for m in C.MODELS}          # one marker per model in the paper (r2-m6)
    for r, (name, key) in enumerate(rows):
        y = ny - 1 - r
        ax.axhline(y, color='0.93', lw=0.5, zorder=0)
        ks = {m: (T[m][key] if isinstance(T[m][key], list) else [T[m][key]]) for m in C.MODELS}
        for m in C.MODELS:
            x = [T[m]['fj'][k] for k in ks[m]]
            if len(x) > 1:          # the four coupled benzo C-H stretches: one tick per mode
                ax.plot(x, [y + off[m]] * len(x), ls='none', marker='|', ms=5.5, mew=0.45,
                        color=style.MODEL_COLOURS[m], zorder=3)          # thin: 4 modes within 24-30 cm-1
            else:
                ax.plot(x, [y + off[m]] * len(x), ls='none', zorder=3, **mk[m])
    ax.set_yticks(range(ny))
    ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=style.TICK_PT)
    ax.tick_params(axis='y', length=0)
    ax.set_ylim(-0.6, ny - 0.4)
    ax.set_xlim(2050, 3750)
    ax.xaxis.set_major_locator(plt.MultipleLocator(500))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(100))
    ax.set_xlabel('Wavenumber (cm$^{-1}$)')
    # free minus bonded N-H: horizontal double arrows between the two N-H rows
    y_n6, y_n4 = ny - 1, ny - 2
    for m, yy in (('gas', y_n4 + 0.62), ('pcm', y_n4 + 0.38)):
        x0, x1 = T[m]['fj'][T[m]['n4']], T[m]['fj'][T[m]['n6']]
        ax.annotate('', (x0, yy), (x1, yy), arrowprops=dict(arrowstyle='<|-|>,head_length=2.2,head_width=1.0',
                                                            color=style.MODEL_COLOURS[m], lw=0.7, shrinkA=0,
                                                            shrinkB=0, mutation_scale=1))
        ax.plot([x0, x0], [y_n4 + off[m] + 0.09, yy], color=style.MODEL_COLOURS[m], lw=0.5, ls=(0, (1.5, 1.2)),
                zorder=1)
        ax.plot([x1, x1], [y_n6 + off[m] - 0.09, yy], color=style.MODEL_COLOURS[m], lw=0.5, ls=(0, (1.5, 1.2)),
                zorder=1)
        dv = x1 - x0
        va = 'bottom' if m == 'gas' else 'top'
        # labels lifted clear of their arrows; the IEF-PCM one shifted left, away from the dotted line at the
        # N6-H6N wavenumber (QA round 1, m8)
        xc = (x0 + x1) / 2 if m == 'gas' else x0 + 30.0          # IEF-PCM: starts right of the dotted line at x0
        ax.annotate(f'Δν = {dv:.0f} cm$^{{-1}}$', (xc, yy), xytext=(0, 2.2 if m == 'gas' else -2.4),
                    textcoords='offset points', ha='center' if m == 'gas' else 'left', va=va,
                    fontsize=style.SMALL_PT, color=style.MODEL_COLOURS[m])


def rotation_glyph(ax, model, ion, sign, colour=C.MODE_COLOUR, n=240, span=(-200.0, 55.0)):
    """Hoop round the rotation axis of one ion (dict from _vib.rigid_decomposition) in the drawing coordinates of
    `model` (C.oriented).  Radius = RMS distance of the ion's heavy atoms from the axis.  The part of the hoop on
    the far side of the ion's centre is drawn behind the molecule (lighter), the near part on top of it and ends in
    the arrowhead, so the arrow gives the sense of rotation of the side facing the reader."""
    from matplotlib.patches import FancyArrowPatch
    from matplotlib.path import Path as MPath
    from matplotlib.lines import Line2D
    import matplotlib.patheffects as pe
    R, focal = C.views()[model]
    R = np.asarray(R)
    P = (ion['com'] - focal) @ R.T
    w = sign * (R @ ion['w'])
    a = w / np.linalg.norm(w)
    X = C.geometries()[model]
    heavy = [i for i in ion['idx'] if C.geometries()['symbols'][i] != 'H']
    V = (X[heavy] - focal) @ R.T - P
    rho = float(np.sqrt(np.mean(np.sum((V - np.outer(V @ a, a)) ** 2, axis=1))))
    ref = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = np.cross(a, ref)
    u /= np.linalg.norm(u)
    v = np.cross(a, u)                                  # (u, v, a) right-handed: +theta = rotation about +a
    th_f = np.arctan2(v[2], u[2])                       # nearest point of the hoop to the reader
    th = th_f + np.radians(np.linspace(span[0], span[1], n))
    pts = P + rho * (np.outer(np.cos(th), u) + np.outer(np.sin(th), v))
    front = pts[:, 2] >= P[2]
    k0 = int(np.argmax(front))                          # first point of the near half (the arc ends in it)
    halo = [pe.Stroke(linewidth=HOOP_LW + 1.3, foreground='white'), pe.Normal()]
    back = pts[:k0 + 1]
    ax.add_line(Line2D(back[:, 0], back[:, 1], color=colour, lw=HOOP_LW * 0.8, alpha=0.55, zorder=1.9,
                       solid_capstyle='round'))
    fr = pts[k0:]
    path = MPath(fr[:, :2])
    fa = FancyArrowPatch(path=path, arrowstyle='-|>,head_length=3.4,head_width=1.8', mutation_scale=1.0,
                         color=colour, lw=HOOP_LW, zorder=60, shrinkA=0, shrinkB=0, joinstyle='miter',
                         capstyle='butt')
    fa.set_path_effects(halo)
    ax.add_patch(fa)
    # axis: short dash-dotted stubs just outside the hoop on both sides
    for sg in (1, -1):
        q0, q1 = P + sg * a * (rho * 0.55), P + sg * a * (rho + 1.0)
        ln = Line2D([q0[0], q1[0]], [q0[1], q1[1]], color='0.3', lw=0.6, ls=(0, (3.0, 1.2, 0.8, 1.2)),
                    zorder=61 if (q1[2] >= P[2]) else 1.9)
        ln.set_path_effects([pe.Stroke(linewidth=1.8, foreground='white'), pe.Normal()])
        ax.add_line(ln)
    return dict(rho=rho, axis_view=a, segments=[(fr[i, :2], fr[i + 1, :2], 1.2) for i in range(0, len(fr) - 1, 8)])


def check_offaxis_caption():
    """Caption of (b) (no facts keys; QA round 2, r2-m10): the gas-phase N4-H4N stretch moves H4N 19 deg off the
    N4-H4N axis, and the bond points 37 deg out of the page in the common camera.  Recomputed and asserted."""
    md = VB.modes('gas')
    k = VB.find_mode(md, float(C.fact('gas.nu_N4_H4N_stretch')), tol=0.06)
    X = C.geometries()['gas']
    n, h = C.I['N4'], C.I['H4N']
    b = (X[h] - X[n]) / np.linalg.norm(X[h] - X[n])
    d = md['disp'][k][h]
    off = float(np.degrees(np.arccos(abs(d @ b) / np.linalg.norm(d))))
    tilt = float(np.degrees(np.arcsin(abs((np.asarray(C.views()['gas'][0]) @ b)[2]))))
    print(f'gas N4-H4N stretch: H4N moves {off:.1f} deg off the bond axis; bond {tilt:.1f} deg out of the page')
    # QA round 4, r4-m7: asserted against the facts key the caption quotes, not a hard-coded 19; the 37 deg tilt is no
    # longer quoted, only that the bond points out of the page
    C.check_close(off, 'gas.nu_N4_H4N_stretch_offaxis_deg', 1e-3, 'gas N4-H4N stretch: H4N off the bond axis')
    assert tilt > 20.0, tilt


def main():
    check_offaxis_caption()
    style.apply()
    g = C.geometries()
    T = mode_table()
    near = [C.I[a] for a in NEAR]
    # phases: stretch drawn with H4N moving towards O1; lowest modes with the cation centroid moving up-screen
    Xg = {m: g[m] for m in C.MODELS}
    ph_n4 = {m: (C.I['H4N'], Xg[m][0] - Xg[m][31]) for m in C.MODELS}
    ph_low = {m: (C.I['N6'], C.views()[m][0][1]) for m in C.MODELS}

    def f_top(m):
        return lambda ax: VB.draw_mode(ax, m, T[m]['n4'], atoms=near, label_atoms=['N4', 'H4N', 'O1', 'N5'],
                                       hbonds=HB_ALL[m], phase=ph_n4[m])

    heavy = [i for i, x in enumerate(g['symbols']) if x != 'H']

    # (d, e): rigid-body decomposition of the lowest mode; the phase drawn moves N6 up-screen
    dec, sgn = {}, {}
    for m in C.MODELS:
        dec[m] = VB.rigid_decomposition(m, T[m]['low'])
        key = 'gas.mode1_rigid_fraction' if m == 'gas' else 'saddle.imag_mode_rigid_fraction'
        C.check_close(dec[m]['rigid_fraction'], key, 1e-6, f'{m} rigid-body fraction of the lowest mode')
        cat = dec[m]['ions']['cation']
        d_n6 = cat['t'] + np.cross(cat['w'], g[m][C.I['N6']] - cat['com'])
        sgn[m] = 1.0 if d_n6 @ C.views()[m][0][1] >= 0 else -1.0
        for name, ion in dec[m]['ions'].items():
            assert ion['trans_share'] < 0.05 and ion['rot_share'] > 0.3, (m, name, ion['trans_share'],
                                                                          ion['rot_share'])
            print(f'{m} lowest mode, {name}: rotation {ion["rot_share"]:.3f}, translation {ion["trans_share"]:.3f} '
                  f'of the mode')
        # caption statements: the two ions turn in opposite senses (axes more than 150 deg apart); in IEF-PCM the
        # cation axis lies within 11 deg of the line joining the two ring centroids
        wp, wc = dec[m]['ions']['picrate']['w'], dec[m]['ions']['cation']['w']
        ang = lambda a, b: float(np.degrees(np.arccos(np.clip(a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1, 1))))
        link = g[m][C.CAT9].mean(0) - g[m][C.PIC_RING].mean(0)
        a_pc, a_cl = ang(wp, wc), min(ang(wc, link), 180 - ang(wc, link))
        print(f'{m}: angle between the ion rotation axes {a_pc:.1f} deg; cation axis to ring-centroid line {a_cl:.1f}')
        pre = 'saddle.imag_mode' if m == 'pcm' else 'gas.mode1'                      # QA round 4, r4-m7: facts keys
        C.check_close(a_pc, f'{pre}_axes_angle_deg', 1e-3, f'{m}: angle between the ion rotation axes')
        C.check_close(a_cl, f'{pre}_cation_axis_to_centroid_line_deg', 1e-3, f'{m}: cation axis to centroid line')
        assert a_pc > 150.0 and (m == 'gas' or a_cl < 11.0), (m, a_pc, a_cl)       # caption: 'within 11 deg'

    def f_bot(m):
        def f(ax):
            segs = []
            for ion in dec[m]['ions'].values():
                gl = rotation_glyph(ax, m, ion, sgn[m])
                segs += gl['segments']
            return VB.draw_mode(ax, m, T[m]['low'], atoms=None, label_atoms=['O1', 'N4', 'N6'],
                                hbonds=HB_ALL[m], arrow_atoms=[], extra_segments=segs)
        return f

    def box(m, idx):
        X = C.oriented(m)[idx]
        return np.array([X[:, 0].min(), X[:, 0].max(), X[:, 1].min(), X[:, 1].max()])

    W_FIT = WIDTH_MM - 2 * 25.4 * float(plt.rcParams['savefig.pad_inches'])      # savefig adds the pad
    right_w = W_FIT - A_W - GAP_MM
    s1, b1 = C.fit_row_scale([f_top(m) for m in C.MODELS], right_w, GAP_MM, [box(m, near) for m in C.MODELS], 8.0)
    s2, b2 = C.fit_row_scale([f_bot(m) for m in C.MODELS], right_w, GAP_MM, [box(m, list(range(33)))
                                                                           for m in C.MODELS], 4.5)
    h1 = (b1[0][3] - b1[0][2]) * s1 + 1.0
    h2 = (b2[0][3] - b2[0][2]) * s2 + 1.0
    H = HEAD_MM + h1 + 2.0 + HEAD_MM + h2
    print(f'scales {s1:.2f} / {s2:.2f} mm per Å; height {H:.1f} mm')
    fig = plt.figure(figsize=(C.mm(WIDTH_MM), C.mm(H)))
    x0 = A_W + GAP_MM
    titles = {}
    for m in C.MODELS:
        fj = T[m]['fj']
        titles[(m, 'n4')] = f'{style.MODEL_LABELS[m]}, {wn(fj[T[m]["n4"]])} cm$^{{-1}}$'
        titles[(m, 'low')] = f'{style.MODEL_LABELS[m]}, {wn(fj[T[m]["low"]])} cm$^{{-1}}$'
    # printed wavenumbers equal facts.json at the printed precision
    for m in C.MODELS:
        C.check_equal(T[m]['fj'][T[m]['n4']], f'{m}.nu_N4_H4N_stretch', 0)
    C.check_equal(T['gas']['fj'][0], 'gas.freq_lowest', 0)
    C.check_equal(T['pcm']['fj'][0], 'pcm.freq_lowest', 1)
    for row, (key, boxes, s, fns) in enumerate((('n4', b1, s1, f_top), ('low', b2, s2, f_bot))):
        y_top = H if row == 0 else H - HEAD_MM - h1 - 2.0
        hh = h1 if row == 0 else h2
        x = x0
        for j, (m, b) in enumerate(zip(C.MODELS, boxes)):
            ax, w, h = C.drawing_axes(fig, x, y_top - HEAD_MM - hh, b, s)
            fns(m)(ax)
            C.header(fig, x, y_top - HEAD_MM + 1.4, 'bcde'[2 * row + j], titles[(m, key)])
            x += w + GAP_MM
    # (a)
    y_bot = 10.5
    axa = C.mm_box(fig, 15.5, y_bot, A_W - 15.5 - 1.5, H - HEAD_MM - 1.0 - y_bot)
    stretch_chart(axa, T)
    hl = [Line2D([], [], ls='none', label=style.MODEL_LABELS[m], **C.model_marker(m, 3.9)) for m in C.MODELS]
    C.header(fig, 0.0, H - HEAD_MM + 1.4, 'a', handles=hl, handle_len=1.0)
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
