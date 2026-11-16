"""Fig. 5 fig_hbond_trajectory: where the proton sits, and what moved during the IEF-PCM optimisation.

(a) N4-H4N, H4N...O1 and N4...O1 in the P2_1/c crystal (X-ray), IEF-PCM (water) and the gas phase (facts.json),
    with the two N...O1 distances of the P2_1/n polymorph (zeng2011benzotriazolium) and the N...O range below
    which periodic GGA / meta-GGA misplaced protons (abramov2024it), both parsed from the abstract-checked claims
    in papers/btpc/lit/candidates.jsonl;
(b-d) the 70-step IEF-PCM optimisation from opt_traj.xyz: (b) energy relative to the final structure,
    (c) change of the three hydrogen-bond distances from step 0, (d) angle between the picrate and cation ring
    planes; steps 40-70 shaded, where 94% of the displacement lies along the 3.2i cm-1 mode
    (saddle.drift_40_to_70_overlap_imag). The eight saddle.traj_stepNN checkpoints are asserted.

    python scripts/figures/fig_hbond_trajectory.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

from figkit import style, plots
from figkit.draw2d import HBOND_COLOUR
import _struct_common as C

STEM = 'fig_hbond_trajectory'
WIDTH_MM = 175.0
TEAL = HBOND_COLOUR
# one marker per model across the paper (QA round 2, r2-m6): X-ray filled circle, gas filled triangle, IEF-PCM open
# square (_struct_common.model_marker)
MK = {m: C.model_marker(m, 4.4) for m in ('xray', 'pcm', 'gas')}
MK['lit'] = dict(marker='D', ms=3.8)
# N4...O1 in brown, not the picrate purple of Figs. 9 and 10 (QA round 1, m7); N4-H4N in neutral grey and the
# single-series curves of (b, d) in IEF-PCM blue, not the near-black of the X-ray points in (a) (QA round 2, r2-m7)
DIST_COLOURS = {'N4–H4N': '0.5', 'H4N···O1': TEAL, 'N4···O1': '#8C510A'}
PCM_BLUE = style.MODEL_COLOURS['pcm']
SPAN = (40, 70)
SPAN_COLOUR = '#EDEDED'


def panel_a(ax):
    rows = [('N4–H4N', 'dh'), ('H4N···O1', 'ha'), ('N4···O1', 'da')]
    off = {'xray': 0.24, 'pcm': 0.0, 'gas': -0.24}
    lit = C.literature_claim('zeng2011benzotriazolium', r'N4-H4\.\.\.O1 (\d\.\d+) A, N6-H6\.\.\.O1 (\d\.\d+) A')
    limit = C.literature_claim('abramov2024it', r'shorter than about (\d\.\d+) A')[0]
    ny = len(rows)
    y_no = 0
    # periodic-DFT strip in the N...O row
    x0 = 0.8
    ax.add_patch(plt.Rectangle((x0, y_no - 0.5), limit - x0, 1.0, facecolor=SPAN_COLOUR, edgecolor='none',
                               zorder=0))
    ax.plot([limit, limit], [y_no - 0.5, y_no + 0.5], color='0.55', lw=0.6, ls=(0, (2, 1.5)), zorder=0.5)
    # stated as the benchmark result it is (QA round 1, m7), not as a rule, and not as a failure zone of this method
    # (referee_r1 m11: the benchmark tested periodic GGA/meta-GGA, not B3LYP)
    ax.text((x0 + limit) / 2 - 0.1, y_no - 0.12, f'N···O < {limit:.2f} Å: periodic GGA/\nmeta-GGA benchmark (not B3LYP)',
            ha='center', va='center', fontsize=style.SMALL_PT, color='0.35', linespacing=1.15)
    vals = {}
    for r, (name, key) in enumerate(rows):
        y = ny - 1 - r
        if r:
            ax.axhline(y + 0.5, color='0.85', lw=0.5, zorder=0)
        for m in ('xray', 'pcm', 'gas'):
            v = float(C.fact(f'{m}.hb_N4_H4N_O1_{key}'))
            vals[(m, key)] = v
            ax.plot([v], [y + off[m]], ls='none', zorder=3, **MK[m])
    # the two polymorph values are 0.03 Å apart: stack them so the markers do not overlap
    for v, dy in zip(sorted(lit), (off['xray'], off['xray'] + 0.17)):
        ax.plot([v], [y_no + dy], ls='none', mfc='white', mec=C.LIT_GREY, mew=0.9, zorder=3, **MK['lit'])
    ax.set_yticks(range(ny))
    # the N...O row also holds the polymorph's N6...O1 distance, hence 'N···O1' (QA round 1, m7)
    ax.set_yticklabels([r[0] if r[1] != 'da' else 'N···O1' for r in rows][::-1])
    ax.tick_params(axis='y', length=0)
    ax.set_ylim(-0.5, ny - 0.5)
    ax.set_xlim(x0, 2.9)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.set_xlabel('Distance (Å)')
    # the gas-phase N...O lies inside the strip, the IEF-PCM one outside (checked, it is the claim of the panel)
    assert vals[('gas', 'da')] < limit < vals[('pcm', 'da')]
    return vals, lit, limit


def main():
    style.apply()
    T = C.pcm_trajectory()
    s = T['step']
    overlap = float(C.fact('saddle.drift_40_to_70_overlap_imag'))
    nu_imag = float(C.fact('pcm.freq_lowest'))
    assert nu_imag < 0
    fig = plt.figure(figsize=(C.mm(WIDTH_MM), C.mm(100.0)))
    W, H = WIDTH_MM, 100.0
    # ---- (a) ----
    ax_a = C.mm_box(fig, 17.0, 27.0, 50.0, 64.0)
    vals, lit, limit = panel_a(ax_a)
    C.header(fig, 0.0, 95.0, 'a', 'Proton position')
    hl = [Line2D([], [], ls='none', label='X-ray ($P2_1/c$)', **MK['xray']),
          Line2D([], [], ls='none', label=style.MODEL_LABELS['pcm'], **MK['pcm']),
          Line2D([], [], ls='none', label=style.MODEL_LABELS['gas'], **MK['gas']),
          Line2D([], [], ls='none', mfc='white', mec=C.LIT_GREY, mew=0.9, label='X-ray ($P2_1/n$ polymorph)',
                 **MK['lit'])]
    fig.legend(handles=hl, loc='upper left', bbox_to_anchor=(0.0, 17.0 / H), ncol=2, frameon=False,
               fontsize=style.TICK_PT, handletextpad=0.3, columnspacing=1.0, borderaxespad=0.0, labelspacing=0.5)
    # ---- (b-d) ----
    xl, wr = 90.0, 73.2            # left edge and width of the right-hand axes (end labels to their right)
    ph, hd, gap = 21.0, 5.0, 1.8   # panel height, header strip, gap
    y_top = 95.0
    axes = []
    titles = [('b', 'Energy relative to the final structure (kcal mol$^{-1}$)'),
              ('c', 'Change from step 0 (Å)'),
              ('d', 'Angle between the ring planes (°)')]
    for k, (letter, title) in enumerate(titles):
        y_head = y_top - k * (hd + ph + gap)
        ax = C.mm_box(fig, xl, y_head - 1.4 - ph, wr, ph)
        C.header(fig, xl - 12.0, y_head, letter, title)
        axes.append(ax)
    ax_b, ax_c, ax_d = axes
    for ax in axes:
        ax.axvspan(*SPAN, color=SPAN_COLOUR, lw=0, zorder=0)
        ax.set_xlim(-1.5, 71.5)
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(5))
    for ax in (ax_b, ax_c):
        ax.tick_params(labelbottom=False)
    ax_d.set_xlabel('Optimisation step (IEF-PCM)')
    # (b) energy
    ax_b.plot(s, T['dE_kcal'], '-o', color=PCM_BLUE, ms=1.8, lw=0.8, zorder=3)
    ax_b.axhline(0, color='0.6', lw=0.4, ls=(0, (3, 2)), zorder=1)
    ax_b.set_ylim(-0.05, 0.5)
    ax_b.yaxis.set_major_locator(MultipleLocator(0.2))
    ax_b.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax_b.text(55.0, 0.44, f'{100 * overlap:.0f}% along the\n{abs(nu_imag):.1f}i cm$^{{-1}}$ mode', ha='center',
              va='top', fontsize=style.SMALL_PT, color='0.3', linespacing=1.15)
    # (c) distances, change from step 0
    items = []
    for lab in ('N4–H4N', 'H4N···O1', 'N4···O1'):
        d = T[lab] - T[lab][0]
        ax_c.plot(s, d, '-o', color=DIST_COLOURS[lab], ms=1.6, lw=0.8, zorder=3)
        items.append((lab, s[-1], d[-1], DIST_COLOURS[lab]))
    ax_c.axhline(0, color='0.6', lw=0.4, ls=(0, (3, 2)), zorder=1)
    ax_c.set_ylim(-0.055, 0.025)
    ax_c.yaxis.set_major_locator(MultipleLocator(0.02))
    ax_c.yaxis.set_minor_locator(MultipleLocator(0.01))
    # caption: no distance moves more than 0.05 Å from its starting value (max |d - d0| = 0.048 Å, H4N...O1);
    # the full min-max range of H4N...O1 is 0.051 Å, so the wording must be 'from its starting value'
    dev = max(np.abs(T[lab] - T[lab][0]).max() for lab in ('N4–H4N', 'H4N···O1', 'N4···O1'))
    assert dev < 0.05, dev
    # (d) inter-ring angle
    ax_d.plot(s, T['interplanar'], '-o', color=PCM_BLUE, ms=1.8, lw=0.8, zorder=3)
    ax_d.set_ylim(0, 65)
    ax_d.yaxis.set_major_locator(MultipleLocator(20))
    ax_d.yaxis.set_minor_locator(MultipleLocator(10))
    C.check_equal(T['interplanar'][-1], 'pcm.interplanar_picrate_cation', 1)
    plots.end_labels(ax_c, items, fontsize=style.SMALL_PT, pad_pt=5.5)       # start beyond the right spine
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
