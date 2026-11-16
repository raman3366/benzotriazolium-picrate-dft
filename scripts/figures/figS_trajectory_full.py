"""SI figure figS_trajectory_full: the full course of the 70-step IEF-PCM (water) optimisation.

(a) energy relative to the final structure; (b) N4-H4N; (c) H4N...O1, N4...O1 and H4N...O7; (d) angle between
the ring planes and the three nitro-group twists; (e) H11...O2 (the C-H...O contact that exists only in the
gas-phase minimum). From data/r2/.../btpc_pcm/opt_traj.xyz via _struct_common.pcm_trajectory, which asserts the
eight saddle.traj_stepNN checkpoints; the final values are asserted against facts.json. Steps 40-70 are shaded as
in Fig. 5. The gas-phase optimisation took three energy evaluations (gas.opt_steps) and is not plotted.

    python scripts/figures/figS_trajectory_full.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from figkit import style, plots
from figkit.draw2d import HBOND_COLOUR
import _struct_common as C

STEM = 'figS_trajectory_full'
WIDTH_MM = 85.0
SPAN = (40, 70)
SPAN_COLOUR = '#EDEDED'
# no fragment hues (picrate purple, cation green) for series (QA round 1, m16); N4...O1 brown as in Fig. 5; the
# single-series panels (a, b, e) and the inter-ring angle in IEF-PCM blue, as in Fig. 5b, d (QA round 2, r2-m7)
TEAL = HBOND_COLOUR
PCM_BLUE = style.MODEL_COLOURS['pcm']
# QA round 3, r3-m5: in Fig. 5c mid grey means N4-H4N; here H4N...O7 is dark grey with OPEN markers, a style used
# nowhere else in Fig. 5 or this figure, so grey does not change meaning between the two (N4-H4N is the
# single-series blue of (b), as in Fig. 5)
O7_GREY = '0.25'
OPEN = {'H4N···O7'}


def main():
    style.apply()
    T = C.pcm_trajectory()
    s = T['step']
    # final values against facts.json at table precision
    for k, key in (('N4–H4N', 'pcm.hb_N4_H4N_O1_dh'), ('H4N···O1', 'pcm.hb_N4_H4N_O1_ha'),
                   ('N4···O1', 'pcm.hb_N4_H4N_O1_da'), ('H4N···O7', 'pcm.hb_N4_H4N_O7_ha'),
                   ('H11···O2', 'pcm.hb_C11_H11_O2_ha')):
        C.check_equal(T[k][-1], key, 3, k)
    C.check_equal(T['interplanar'][-1], 'pcm.interplanar_picrate_cation', 1)
    for n in ('N1', 'N2', 'N3'):
        C.check_equal(T[f'twist_{n}'][-1], f'pcm.nitro_twist_{n}', 1)
    panels = [
        ('a', 'Energy relative to the final structure (kcal mol$^{-1}$)',
         [('', T['dE_kcal'], PCM_BLUE)], (-0.05, 0.5), 0.2),
        ('b', 'N4–H4N (Å)', [('', T['N4–H4N'], PCM_BLUE)], (1.046, 1.055), 0.002),
        ('c', 'Hydrogen-bond contacts (Å)',
         [('N4···O1', T['N4···O1'], '#8C510A'), ('H4N···O7', T['H4N···O7'], O7_GREY),
          ('H4N···O1', T['H4N···O1'], TEAL)], (1.5, 2.9), 0.4),
        ('d', 'Ring and nitro-group angles (°)',
         [('Inter-ring', T['interplanar'], PCM_BLUE), ('NO$_2$ N1', T['twist_N1'], '#B2182B'),
          ('NO$_2$ N3', T['twist_N3'], '0.42'), ('NO$_2$ N2', T['twist_N2'], '0.6')], (-3, 65), 20),
        ('e', 'H11···O2 (Å)', [('', T['H11···O2'], PCM_BLUE)], (2.9, 4.4), 0.5),
    ]
    ph, hd, gap, bottom, top = 17.0, 4.6, 1.2, 10.0, 1.0
    H = top + len(panels) * (hd + ph) + (len(panels) - 1) * gap + bottom
    fig = plt.figure(figsize=(C.mm(WIDTH_MM), C.mm(H)))
    xl, wr = 11.0, 61.2
    axes = []
    for k, (letter, title, series, ylim, major) in enumerate(panels):
        y_head = H - top - 3.2 - k * (hd + ph + gap)
        ax = C.mm_box(fig, xl, y_head - 1.4 - ph, wr, ph)
        C.header(fig, 0.0, y_head, letter, title)
        ax.axvspan(*SPAN, color=SPAN_COLOUR, lw=0, zorder=0)
        ax.set_xlim(-1.5, 71.5)
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        items = []
        for lab, y, col in series:
            if lab in OPEN:          # open markers (r3-m5)
                ax.plot(s, y, '-o', color=col, ms=1.5, lw=0.7, mfc='white', mec=col, mew=0.5, zorder=3)
            else:
                ax.plot(s, y, '-o', color=col, ms=1.4, lw=0.7, zorder=3)
            if lab:
                items.append((lab, s[-1], y[-1], col))
        lo, hi = ylim
        assert all(lo <= np.min(y) and np.max(y) <= hi for _, y, _ in series), (title, ylim)
        ax.set_ylim(lo, hi)
        ax.yaxis.set_major_locator(MultipleLocator(major))
        if letter == 'a':
            ax.axhline(0, color='0.6', lw=0.4, ls=(0, (3, 2)), zorder=1)
        if k < len(panels) - 1:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel('Optimisation step (IEF-PCM)')
        axes.append((ax, items))
    for ax, items in axes:
        if items:
            plots.end_labels(ax, items, fontsize=style.SMALL_PT, pad_pt=5.5)
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
