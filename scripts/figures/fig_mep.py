"""Fig. 11 fig_mep: molecular electrostatic potential of the benzotriazolium picrate ion pair.

(a, b) ESP (esp.cube, kcal/mol) mapped on the rho = 0.001 a.u. isodensity surface (density.cube), gas phase and
IEF-PCM (water), camera of the structure figures, ONE colour range for both models (outer multiples of 5 of the
surface extremes, green at zero).  No extremum spheres: the global minimum and maximum of each surface (Multiwfn,
mep_extrema.csv) are named by callouts outside the surface, with short leaders to the point (or, when it lies on
the far side, to the atom it belongs to).  (c) Surface area per 2.5 kcal/mol interval of the ESP (Multiwfn
vertices, vtx.txt).  (d) Every surface minimum and maximum found by Multiwfn; |V| > 30 kcal/mol (and the
gas-phase O1 minimum) labelled by the nearest atom, from a label column beside each strip with a leader to every
marker.  All printed numbers are read from the data files and asserted against facts.json.

    python scripts/figures/fig_mep.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy import ndimage

from figkit import style, compose, callouts, plots, render3d as r3
import _render_common as C

STEM = 'fig_mep'
ISO = 0.001
BAND_MM = 58.0                  # (d) needs the height for its labelled extrema (QA round 1, M7)
LABEL_ABOVE = 30.0              # kcal/mol: extrema labelled by site in (d)
# white room round each render for the two callouts (global minimum on the left rim, maximum at the upper right)
MARGINS_MM = dict(left=8.5, top=1.0, right=1.0, bottom=0.5)
MAX_BEHIND_RIM_A = 0.5          # an extremum hidden behind the outline by less than this is pointed at on the rim
CALLOUT_PT = style.TICK_PT


def extrema(model):
    """Rows of the curated mep_extrema.csv with clean site labels; the global pair checked against facts."""
    rows = r3.read_mep_extrema(r3.case_dir(C.CASES[model], raw=False) / 'mep_extrema.csv')
    raw = r3.read_mep_extrema(r3.case_dir(C.CASES[model]) / 'mep_extrema.csv')
    assert [(r['type'], r['kcal']) for r in rows] == [(r['type'], r['kcal']) for r in raw]
    for r in rows:
        r['site'] = r3.clean_label(r['near'].split('(')[0].strip())
        assert '#' not in r['site']
    lo = min(rows, key=lambda r: r['kcal'])
    hi = max(rows, key=lambda r: r['kcal'])
    assert lo['is_global'] and lo['type'] == 'minimum'
    # the CSV rounds to 3 decimals and facts now hold Multiwfn's unrounded values (round 2, QA r4-m1), so compare at
    # the 1 decimal the figure prints
    C.check_equal(lo['kcal'], f'{model}.mep_min_kcal', 1)
    C.check_equal(hi['kcal'], f'{model}.mep_max_kcal', 1)
    assert lo['site'] == C.fact(f'{model}.mep_min_site') and hi['site'] == C.fact(f'{model}.mep_max_site')
    n_min = sum(r['type'] == 'minimum' for r in rows)
    assert n_min == C.fact(f'{model}.mep_n_minima') and len(rows) - n_min == C.fact(f'{model}.mep_n_maxima')
    return rows, lo, hi


def callout_target(rd, pos):
    """Where the leader of a global extremum ends: the projected surface point when it faces the camera;
    otherwise (an in-plane lone-pair or N-H extremum just behind the outline) the nearest point of the outline,
    provided the extremum is less than MAX_BEHIND_RIM_A from it.  Returns ((x, y) image px, how, distance A)."""
    P = np.atleast_2d(pos)
    xy = rd.project(P)[0][0]
    if bool(rd.visible(P, towards_camera=0.1, window=2)[0]):
        return xy, 'on the visible surface', 0.0
    content = callouts.content_mask(rd.image)
    d, idx = ndimage.distance_transform_edt(content, return_indices=True)
    i, j = int(xy[1]), int(xy[0])
    dist_A = float(d[i, j] / rd.px_per_A)
    if dist_A > MAX_BEHIND_RIM_A:
        raise AssertionError(f'extremum at {pos} is {dist_A:.2f} A inside the outline and hidden: no honest target')
    bi, bj = idx[:, i, j]
    v = np.array([j - bj, i - bi], float)
    v /= max(np.hypot(*v), 1e-9)
    return np.array([bj + 0.5, bi + 0.5]) + 2.0 * v, 'on the outline', dist_A


def strip_positions(values, x0, dy, dx):
    """Beeswarm offsets: points closer than dy in value are spread sideways by dx steps."""
    order = np.argsort(values)
    xs = np.full(len(values), float(x0))
    placed = []
    for k in order:
        v = values[k]
        for step in (0, 1, -1, 2, -2, 3, -3):
            x = x0 + step * dx
            if all(abs(v - pv) >= dy or abs(x - px) >= 0.99 * dx for pv, px in placed):
                break
        xs[k] = x
        placed.append((v, x))
    return xs


def blockers(ax, segments, others, r_pt=2.4):
    """Markers of `others` ((x, y) data) that lie within r_pt points of any segment ((x0, y0), (x1, y1))."""
    T = ax.transData
    k = 72.0 / ax.figure.dpi
    out = []
    for p0, p1 in segments:
        a, b = T.transform(p0) * k, T.transform(p1) * k
        d = b - a
        for o_ in others:
            o = T.transform(o_) * k
            t = np.clip(np.dot(o - a, d) / max(np.dot(d, d), 1e-9), 0, 1)
            if np.linalg.norm(a + t * d - o) < r_pt and o_ not in out:
                out.append(o_)
    return out


def column_positions(y0, sep, lim):
    """Label heights for a column: as close as possible to the wanted heights y0, at least `sep` apart and inside
    `lim` (1D spreading of the ordered labels)."""
    order = np.argsort(y0)
    y = np.array(y0, float)[order]
    for _ in range(5000):
        moved = False
        for i in range(len(y) - 1):
            if y[i + 1] - y[i] < sep - 1e-7:
                d = (sep - (y[i + 1] - y[i])) / 2
                y[i] -= d
                y[i + 1] += d
                moved = True
        clamped = np.clip(y, lim[0], lim[1])
        if np.any(np.abs(clamped - y) > 1e-9):
            y = clamped
            moved = True
        if not moved:
            break
    assert np.all(np.diff(y) >= sep - 1e-6) and y[0] >= lim[0] - 1e-6 and y[-1] <= lim[1] + 1e-6, y
    out = np.empty_like(y)
    out[order] = y
    return out


def main():
    t0 = time.time()
    style.apply()
    mols = C.molecules()
    V = C.views()
    scenes, ext = {}, {}
    for m in C.MODELS:
        d = r3.case_dir(C.CASES[m])
        scenes[m] = r3.mep_scene(d / 'density.cube', d / 'esp.cube', mols[m], iso=ISO)
        ext[m] = extrema(m)
    # one colour range for both models: outer multiples of 5 of every surface value (grid and Multiwfn)
    lo = min(min(s.meta['vmin_surface'] for s in scenes.values()), min(e[1]['kcal'] for e in ext.values()))
    hi = max(max(s.meta['vmax_surface'] for s in scenes.values()), max(e[2]['kcal'] for e in ext.values()))
    vr = (5 * math.floor(lo / 5), 5 * math.ceil(hi / 5))
    for m in C.MODELS:
        r3.recolour_mep(scenes[m], vrange=vr, centre_zero=True)
        scenes[m].colour.ticks = [t for t in range(-50, 101, 25) if vr[0] <= t <= vr[1]]
    print('colour range', vr, 'surface extremes', {m: (round(s.meta['vmin_surface'], 2), round(s.meta['vmax_surface'], 2))
                                                  for m, s in scenes.items()})
    vw = {m: V[m].fit(scenes[m], pad=0.15) for m in C.MODELS}
    pxA = r3.px_per_angstrom(list(vw.values()), C.PANEL2_MM)
    rr = {}
    for m in C.MODELS:
        full = r3.render(scenes[m], vw[m], px_per_A=pxA, crop=False)
        rc = r3.crop_together(full, pad=4)
        # margins in image px: the render is scaled to the panel width later, so convert mm at the final scale
        rr[m] = callouts.expand_canvas(rc, **{k: C.mm_to_px(v) for k, v in MARGINS_MM.items()})
    fig, axes = compose.grid([rr['gas'], rr['pcm']], ncols=2, width=style.DOUBLE,
                             titles=[style.MODEL_LABELS['gas'], style.MODEL_LABELS['pcm']],
                             colourbar_scale=scenes['pcm'].colour, colourbar_len=0.55, extra_bottom_mm=BAND_MM,
                             colourbar_linear=True)
    notes = {}
    for m, ax in zip(C.MODELS, axes):
        rd = ax._figkit_render
        items = []
        for r, kind in ((ext[m][1], 'min'), (ext[m][2], 'max')):
            tgt, how, dA = callout_target(rd, r['pos'])
            notes[(m, kind)] = f'{how} ({dA:.2f} A)'
            items.append(dict(target=tgt, text=f'{C.num(r["kcal"], 1, sign=kind == "max")}\n{r["site"]}',
                              prefer=(-1, 0) if kind == 'min' else (1, 0.6),
                              kw=dict(fontsize=CALLOUT_PT, linespacing=1.05)))
        callouts.place(ax, items, min_leader_mm=1.0, max_leader_mm=4.0, name=f'MEP extrema {m}')
    print('callout targets', notes)

    # ---------------- (c) ESP histogram, (d) extrema strip ----------------
    x0, y0, w, h = fig._figkit_grid['bottom_band']
    W, H = fig.get_size_inches()
    head_in = compose._header_mm() * style.MM
    bottom_in = 9.6 * style.MM
    ah = h - head_in - bottom_in - 1.0 * style.MM
    left_c, w_c = x0 + 10.5 * style.MM, 92.0 * style.MM
    left_d = left_c + w_c + 14.0 * style.MM
    w_d = x0 + w - left_d - 1.0 * style.MM
    axc = fig.add_axes([left_c / W, (y0 + bottom_in) / H, w_c / W, ah / H])
    data = {}
    for m in C.MODELS:
        a = plots.esp_surface_areas(r3.case_dir(C.CASES[m]))
        ck = a['check']
        assert abs(ck['total_A2'] / C.fact(f'{m}.mep_surface_area') - 1) < 0.02, ck
        assert abs(ck['positive_A2'] / C.fact(f'{m}.mep_positive_area') - 1) < 0.03, ck
        assert abs(ck['negative_A2'] / C.fact(f'{m}.mep_negative_area') - 1) < 0.03, ck
        data[m] = a
    bins = np.arange(vr[0], vr[1] + 1e-9, 2.5)
    plots.esp_histogram(axc, data, bins=bins, xlabel='Surface potential $V$ (kcal mol$^{-1}$)',
                        ylabel='Area (Å$^2$ per 2.5 kcal mol$^{-1}$)')
    axc.set_xticks([t for t in range(-50, 101, 25)])
    axc.xaxis.set_major_formatter(compose._formatter('{:.0f}'))
    compose._header(fig, x0 / W, (y0 + h) / H, '(c)', None)

    axd = fig.add_axes([left_d / W, (y0 + bottom_in) / H, w_d / W, ah / H])
    axd.set_xlim(-0.8, 1.8)
    axd.set_xticks([0, 1])
    axd.set_xticklabels(['Gas\nphase', 'IEF-PCM\n(water)'])
    axd.set_ylim(vr[0] - 7, vr[1] + 6)                   # room for the label columns at both ends
    axd.set_yticks([t for t in range(-50, 101, 25)])
    axd.yaxis.set_major_formatter(compose._formatter('{:.0f}'))
    axd.set_ylabel('Extremum (kcal mol$^{-1}$)')
    axd.axhline(0, color='0.6', lw=0.5, ls=(0, (3, 2)), zorder=1)
    marker = {'minimum': 'v', 'maximum': '^'}
    pts = []
    for j, m in enumerate(C.MODELS):
        rows = ext[m][0]
        vals = np.array([r['kcal'] for r in rows])
        xs = strip_positions(vals, j, dy=3.5, dx=0.08)
        for kind in ('minimum', 'maximum'):
            sel = np.array([r['type'] == kind for r in rows])
            axd.scatter(xs[sel], vals[sel], marker=marker[kind], s=11, c=style.MODEL_COLOURS[m], lw=0.3,
                        edgecolors='white', zorder=3)
        pts += [(j, m, r, x) for r, x in zip(rows, xs)]
    from matplotlib.lines import Line2D
    # labels (QA round 1, M7): every labelled extremum gets a leader from a label column beside its strip (gas
    # phase on the left, IEF-PCM on the right); extrema of one type nearest the same atom share one label with
    # one leader each (IEF-PCM: two maxima at N6); the gas-phase O1 minimum is labelled although |V| < 30
    fig.canvas.draw()
    yl0, yl1 = axd.get_ylim()
    kcal_per_pt = (yl1 - yl0) / (axd.bbox.height * 72 / fig.dpi)
    min_sep = (style.SMALL_PT * 0.95 + 0.8) * kcal_per_pt
    lab_rows = []
    for j, m in enumerate(C.MODELS):
        groups = {}
        for jj, mm_, r, x in pts:
            if mm_ != m:
                continue
            if abs(r['kcal']) > LABEL_ABOVE or (r['type'] == 'minimum' and r['site'] == 'O1'):
                groups.setdefault((r['type'], r['site']), []).append((x, r))
        items = []
        for (kind, site), mem in groups.items():
            glob = any(r is ext[m][1] or r is ext[m][2] for x, r in mem)
            xm = float(np.mean([x for x, r in mem]))
            # label on the side of the strip where the marker sits; a marker on the centre line gets the outer side
            # unless a leader to that side would run across a neighbouring marker
            side = -1 if xm < j - 1e-6 else 1 if xm > j + 1e-6 else (-1 if j == 0 else 1)
            own = [(x, r['kcal']) for x, r in mem]
            others = [(x, r['kcal']) for jj, mm_, r, x in pts if mm_ == m and (x, r['kcal']) not in own]
            yw = float(np.mean([r['kcal'] for x, r in mem]))

            def blocked(sd, yy=yw):
                return blockers(axd, [((j + sd * 0.27, yy), t) for t in own], others)
            if abs(xm - j) < 1e-6 and blocked(side) and not blocked(-side):
                side = -side
            bl = blocked(side)
            if bl:        # blocked on both sides: aim the label away from the blocking marker, so the leader slants
                yw += (-0.8 if np.mean([b[1] for b in bl]) > yw else 0.8) * min_sep
            items.append(dict(text=site, bold=glob, targets=own, side=side, y=yw, others=others))
        lab_rows.append((j, m, items))
        for side in (-1, 1):
            col = [it for it in items if it['side'] == side]
            if not col:
                continue
            x_text = j + side * 0.30
            ys = column_positions([it['y'] for it in col], min_sep, (yl0 + 0.62 * min_sep, yl1 - 0.62 * min_sep))
            for it, yl in zip(col, ys):
                axd.text(x_text, yl, it['text'], ha='right' if side < 0 else 'left', va='center',
                         fontsize=style.SMALL_PT, fontweight='bold' if it['bold'] else 'normal', zorder=5)
                x_lead = x_text - side * 0.03
                for tx, ty in it['targets']:
                    axd.annotate('', xy=(tx, ty), xytext=(x_lead, yl),
                                 arrowprops=dict(arrowstyle='-', lw=0.45, color='0.4', shrinkA=0.0, shrinkB=2.6),
                                 zorder=2)
                it['y_label'] = float(yl)
                bad = blockers(axd, [((x_lead, yl), t) for t in it['targets']], it['others'], r_pt=2.0)
                assert not bad, (m, it['text'], bad)          # no leader runs across another marker
        print(m, 'labels', [(it['text'], it['side'], round(it['y'], 2), round(it['y_label'], 1), len(it['targets']))
                            for it in items])
    assert any(it['text'] == 'O1' for j, m, items in lab_rows if m == 'gas' for it in items)
    compose._header(fig, (left_d - 11.5 * style.MM) / W, (y0 + h) / H, '(d)', None)
    # key in the panel header, away from the data (QA round 1, M7)
    hs = [Line2D([], [], ls='none', marker=marker[k], ms=3.6, mfc='0.25', mec='white', mew=0.3, label=k)
          for k in ('maximum', 'minimum')]
    fig.legend(handles=hs, loc='upper left', ncol=2, handletextpad=0.1, columnspacing=0.8, borderaxespad=0.0,
               borderpad=0.0, fontsize=style.SMALL_PT, frameon=False,
               bbox_to_anchor=((left_d - 11.5 * style.MM + 7.5 * style.MM) / W, (y0 + h - 0.35 * style.MM) / H))
    C.save(fig, STEM)
    print(f'total {time.time() - t0:.1f} s')


if __name__ == '__main__':
    main()
