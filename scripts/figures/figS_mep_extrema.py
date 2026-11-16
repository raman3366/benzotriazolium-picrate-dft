"""SI figure figS_mep_extrema: where every MEP surface extremum of the benzotriazolium picrate ion pair lies.

ESP (esp.cube, kcal/mol) on the rho = 0.001 a.u. surface (density.cube), gas phase (top) and IEF-PCM (water)
(bottom), each seen from the front (the camera of Figs. 3, 7, 8 and 11) and from the back (the same camera turned
by 180 deg about the screen vertical), same colour range as main-text Fig. 10.  Every extremum located by Multiwfn
(mep_extrema.csv; 27 per model) is drawn ONCE, in the view where it faces the camera, as its Multiwfn index
(key to Table S5): minima in black squares, maxima in white circles.  Badges sit on the extremum; where two would
collide they are pushed apart and joined to the extremum (small dot) by a short leader.  An extremum hidden in
both views (in the cleft between the ions) is drawn at its projection with a dashed outline.  Every value and
index is checked against facts.json.

    python scripts/figures/figS_mep_extrema.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from matplotlib.lines import Line2D

from figkit import style, compose, qa, render3d as r3
import _render_common as C

STEM = 'figS_mep_extrema'
ISO = 0.001
BADGE_PT = style.SMALL_PT        # 6.5 pt (QA round 1, m17: 6 pt digits printed ~5.5 pt high)
BADGE_GAP_MM = 0.8              # clear gap between badges (was 0.35 mm: badges touched)
SIDES = ('front', 'back')
STYLE = {'minimum': dict(box='square,pad=0.22', fc='#1A1A1A', ec='#1A1A1A', tc='white'),
         'maximum': dict(box='circle,pad=0.14', fc='white', ec='#1A1A1A', tc='black')}


def extrema(model):
    """All Multiwfn extrema (curated csv == raw csv), each checked against its facts key <model>.mep_minNN /
    mep_maxNN; the counts against mep_n_minima / mep_n_maxima."""
    rows = r3.read_mep_extrema(r3.case_dir(C.CASES[model], raw=False) / 'mep_extrema.csv')
    raw = r3.read_mep_extrema(r3.case_dir(C.CASES[model]) / 'mep_extrema.csv')
    assert [(r['type'], r['index'], r['kcal']) for r in rows] == [(r['type'], r['index'], r['kcal']) for r in raw]
    for r in rows:
        key = f"{model}.mep_{'min' if r['type'] == 'minimum' else 'max'}{r['index']:02d}"
        C.check_equal(r['kcal'], key, 3)
    n_min = sum(r['type'] == 'minimum' for r in rows)
    assert n_min == C.fact(f'{model}.mep_n_minima') and len(rows) - n_min == C.fact(f'{model}.mep_n_maxima')
    return rows


def hidden_depth(rd, pos, step=0.1, limit=8.0):
    """How far (A) a surface point must move towards the camera before it is in front of the rendered surface
    (0 = visible; limit = hidden by the whole molecule)."""
    for k in range(int(limit / step) + 1):
        if bool(rd.visible(np.atleast_2d(pos), towards_camera=0.1 + k * step, window=2)[0]):
            return k * step
    return limit


def spread(anchors, radii, lo, hi, gap, iters=400):
    """Badge centres: start on the anchors, push overlapping badges apart (pairwise, symmetric), pull each back
    towards its anchor, keep inside the box lo..hi.  Returns (N, 2) positions."""
    P = np.array(anchors, float)
    A = P.copy()
    R = np.asarray(radii, float)
    for _ in range(iters):
        moved = False
        for i in range(len(P)):
            for j in range(i + 1, len(P)):
                d = P[j] - P[i]
                dist = float(np.hypot(*d))
                need = R[i] + R[j] + gap
                if dist < need:
                    u = d / dist if dist > 1e-9 else np.array([1.0, 0.0])
                    P[i] -= 0.5 * (need - dist) * u
                    P[j] += 0.5 * (need - dist) * u
                    moved = True
        P += 0.02 * (A - P)
        P = np.clip(P, lo + R[:, None], hi - R[:, None])
        if not moved:
            break
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            assert np.hypot(*(P[j] - P[i])) >= R[i] + R[j] + 0.5 * gap, 'badges still overlap'
    return P


def draw_badges(ax, items, mm):
    """items: dicts with xy (image px), text, kind, hidden.  One badge per item, spread apart, leaders where a
    badge had to move off its extremum."""
    fig = ax.figure
    arts = []
    for it in items:
        st = STYLE[it['kind']]
        bb = dict(boxstyle=st['box'], fc=st['fc'], ec=st['ec'], lw=0.5)
        if it['hidden']:
            bb.update(ls=(0, (1.6, 1.0)), fc='white', ec='#1A1A1A')
        arts.append(ax.text(*it['xy'], it['text'], ha='center', va='center', fontsize=BADGE_PT,
                            fontweight='bold', color='black' if it['hidden'] else st['tc'], bbox=bb, zorder=60,
                            clip_on=False))
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    radii = []
    for a in arts:
        e = a.get_bbox_patch().get_window_extent(rend)
        (x0, y0), (x1, y1) = inv.transform([(e.x0, e.y0), (e.x1, e.y1)])
        radii.append(0.5 * max(abs(x1 - x0), abs(y1 - y0)))
    H, W = ax._figkit_render.image.shape[:2]
    P = spread([it['xy'] for it in items], radii, np.array([0.0, 0.0]), np.array([W, H], float), gap=BADGE_GAP_MM * mm)
    moved = 0
    for a, it, p, r in zip(arts, items, P, radii):
        a.set_position(p)
        d = np.asarray(p) - np.asarray(it['xy'])
        if np.hypot(*d) > 1.0 * r:                       # moved off its point: leader + dot at the true position
            moved += 1
            end = np.asarray(p) - d / np.hypot(*d) * r * 0.95
            ax.add_line(Line2D([it['xy'][0], end[0]], [it['xy'][1], end[1]], lw=0.5, color='#1A1A1A', zorder=58))
            ax.plot(*it['xy'], 'o', ms=1.8, mfc='#1A1A1A', mec='white', mew=0.4, zorder=59)
    return arts, moved


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
    # the colour range of main-text Fig. 10: outer multiples of 5 of every surface value (grid and Multiwfn), both models
    lo = min(min(s.meta['vmin_surface'] for s in scenes.values()), min(r['kcal'] for e in ext.values() for r in e))
    hi = max(max(s.meta['vmax_surface'] for s in scenes.values()), max(r['kcal'] for e in ext.values() for r in e))
    vr = (5 * math.floor(lo / 5), 5 * math.ceil(hi / 5))
    for m in C.MODELS:
        r3.recolour_mep(scenes[m], vrange=vr, centre_zero=True)
        scenes[m].colour.ticks = [t for t in range(-50, 101, 25) if vr[0] <= t <= vr[1]]
    base = {(m, s): (V[m] if s == 'front' else V[m].flipped()) for m in C.MODELS for s in SIDES}
    vw = {k: v.fit(scenes[k[0]], pad=0.3) for k, v in base.items()}
    pxA = r3.px_per_angstrom(list(vw.values()), C.PANEL2_MM)
    rr = {k: r3.render(scenes[k[0]], v, px_per_A=pxA, crop=False) for k, v in vw.items()}
    # every extremum goes to the view where it faces the camera (front first); hidden in both -> least hidden
    place = {k: [] for k in vw}
    hidden_note = []
    for m in C.MODELS:
        for r in ext[m]:
            depth = {s: hidden_depth(rr[(m, s)], r['pos']) for s in SIDES}
            side = 'front' if depth['front'] == 0 else ('back' if depth['back'] == 0 else min(SIDES, key=depth.get))
            hidden = depth[side] > 0
            if hidden:
                hidden_note.append((m, r['type'], r['index'], round(depth[side], 1)))
            place[(m, side)].append(dict(r, hidden=hidden))
    assert sum(len(v) for v in place.values()) == sum(len(e) for e in ext.values())
    print('hidden in both views:', hidden_note)
    print('extrema per view:', {f'{m}/{sd}': len(v) for (m, sd), v in place.items()})
    crops = {}
    for m in C.MODELS:
        pair = r3.crop_together([rr[(m, s)] for s in SIDES], pad=int(C.mm_to_px(1.2)))
        crops.update({(m, s): p for s, p in zip(SIDES, pair)})
    order = [(m, s) for m in C.MODELS for s in SIDES]
    # main-text Fig. 10 since round 1 (the descriptor chart was dropped); 'main-text' because inside the SI a bare
    # 'Fig. 10' would read as an SI figure (figures_r1 #3, #18)
    titles = [f'{style.MODEL_LABELS[m]}, {"front (view of main-text Fig. 10)" if s == "front" else "back"}'
              for m, s in order]
    fig, axes = compose.grid([crops[k] for k in order], ncols=2, width=style.DOUBLE, titles=titles,
                             gap_mm=(2.5, 2.0), colourbar_scale=scenes['pcm'].colour, colourbar_len=0.5,
                             colourbar_linear=True)
    n_moved = 0
    for k, ax in zip(order, axes):
        rd = ax._figkit_render
        T = ax.transData
        mm = fig.dpi / 25.4 / float(abs(T.transform((1, 0))[0] - T.transform((0, 0))[0]))
        items = []
        for r in place[k]:
            xy = rd.project(np.atleast_2d(r['pos']))[0][0]
            items.append(dict(xy=xy, text=str(r['index']), kind=r['type'], hidden=r['hidden']))
        _, mv = draw_badges(ax, items, mm)
        n_moved += mv
    print('badges moved off their extremum:', n_moved)
    # key to the symbols, left of the colour bar
    hs = [Line2D([], [], ls='none', marker='s', ms=5.0, mfc='#1A1A1A', mec='#1A1A1A', label='minimum'),
          Line2D([], [], ls='none', marker='o', ms=5.4, mfc='white', mec='#1A1A1A', mew=0.6, label='maximum')]
    cb = fig._figkit_grid['colourbar']
    W, H = fig.get_size_inches()
    bx = cb.ax.get_position()
    fig.legend(handles=hs, loc='center right', bbox_to_anchor=(bx.x0 - 4 * style.MM / W, bx.y0 + bx.height / 2),
               ncol=1, handletextpad=0.3, labelspacing=0.35, borderaxespad=0.0, fontsize=style.TICK_PT)
    C.save(fig, STEM)
    print(f'total {time.time() - t0:.1f} s')


if __name__ == '__main__':
    main()
