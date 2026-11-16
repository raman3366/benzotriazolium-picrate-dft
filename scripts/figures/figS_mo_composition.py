"""SI figure figS_mo_composition: Mulliken composition of HOMO-7 ... LUMO+7 over the picrate and the
benzotriazolium cation, gas phase (a) and IEF-PCM (water) (b).

One row per orbital (names once, on the left); each panel has the orbital energies (eV, 3 decimals) as tick labels,
a 0-100 % bar split into picrate (purple, from the left) and cation (green), and the picrate share (%, 1 decimal)
at the right, printed from the facts string. Data: papers/btpc/data/derived/mo_fragment_composition_btpc_{gas,pcm}.csv; energies and shares
asserted against facts `*.orbitals_near_gap` (in _elec_common.mo_table). Mulliken shares slightly outside
0-100 % (gas LUMO+6/LUMO+7, diffuse virtuals) are printed as computed and drawn clipped to 0-100 %.

    python scripts/figures/figS_mo_composition.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator

from figkit import style
import _elec_common as C

STEM = 'figS_mo_composition'
WIDTH_MM = 120.0
N_BELOW, N_ABOVE = 7, 7
ROW_MM = 3.7
BAR = 0.72


def facts_shares(model):
    """{orbital name (true minus sign): picrate share as printed in facts `*.orbitals_near_gap` (string)}."""
    import re
    items = re.findall(r'(HOMO(?:-\d+)?|LUMO(?:\+\d+)?) (-?\d+\.\d+) eV \((-?\d+\.\d+)% picrate\)',
                       C.fact(f'{model}.orbitals_near_gap'))
    assert len(items) == N_BELOW + N_ABOVE + 2, (model, len(items))
    return {name.replace('-', C.MINUS): pct for name, ev, pct in items}


def main():
    style.apply()
    tabs = {m: C.mo_table(m) for m in C.MODELS}
    n = N_BELOW + N_ABOVE + 2
    rows = {}
    for m in C.MODELS:
        t = tabs[m]
        sel = np.arange(t['n_occ'] - 1 - N_BELOW, t['n_occ'] + 1 + N_ABOVE)[::-1]     # LUMO+7 at the top
        rows[m] = sel
    names = [C.orbital_name(k, tabs['gas']['n_occ']) for k in rows['gas']]
    assert names == [C.orbital_name(k, tabs['pcm']['n_occ']) for k in rows['pcm']]

    # layout (mm)
    W = WIDTH_MM
    # measured at 6.5 pt: 'HOMO−7' 10.2 mm, '−8.612' 7.5 mm, '100.0' 6.1 mm, '% picrate' 10.0 mm
    names_w, e_w, bar_w, pct_w, gap, right = 11.4, 9.0, 34.4, 7.2, 3.5, 10.6
    top, legend_h, bottom = 7.4, 4.4, 9.5
    ph = n * ROW_MM
    H = bottom + ph + top + legend_h
    fig = plt.figure(figsize=(C.mm(W), C.mm(H)))
    x0 = {'gas': names_w + e_w}
    x0['pcm'] = x0['gas'] + bar_w + pct_w + gap + e_w
    assert abs(x0["pcm"] + bar_w + right - W) < 1.2, x0['pcm'] + bar_w + right
    axes = {}
    for m in C.MODELS:
        ax = fig.add_axes([x0[m] / W, bottom / H, bar_w / W, ph / H])
        axes[m] = ax
        t = tabs[m]
        y = np.arange(n)
        pic = np.clip(t['pic'][rows[m]], 0, 100)
        ax.barh(y, pic, height=BAR, color=C.PIC_COLOUR, lw=0, zorder=2)
        ax.barh(y, 100 - pic, left=pic, height=BAR, color=C.CAT_COLOUR, lw=0, zorder=2)
        ax.set_ylim(n - 0.5, -0.5)
        ax.set_xlim(0, 100)
        ax.set_yticks(y)
        ax.set_yticklabels([C.fmt(t['energy'][k], 3) for k in rows[m]], fontsize=style.SMALL_PT)
        ax.tick_params(axis='y', length=0, pad=2)
        ax.xaxis.set_major_locator(MultipleLocator(50))
        ax.xaxis.set_minor_locator(MultipleLocator(25))
        ax.set_xlabel('Mulliken share (%)')
        # occupied / virtual boundary
        b = N_ABOVE + 0.5
        ax.axhline(b, color='0.15', lw=0.8, zorder=3, clip_on=False, xmin=-0.02, xmax=1.02)
        for sp in ('top', 'right', 'left'):
            ax.spines[sp].set_visible(False)
        # picrate share printed from the facts string itself (1 decimal, `*.orbitals_near_gap`), not re-rounded
        # from the CSV: 0.35 % would print as 0.3 through float rounding while facts and Table S7 give 0.4
        # (QA round 1, M9); every printed string must equal the facts string exactly
        shares = facts_shares(m)
        for yy, k in zip(y, rows[m]):
            nm = C.orbital_name(k, t['n_occ'])
            txt = C.fmt(float(shares[nm]), 1)
            assert txt == shares[nm].replace('-', C.MINUS), (m, nm, txt, shares[nm])
            assert abs(t['pic'][k] - float(shares[nm])) <= 0.05 + 1e-9, (m, nm, t['pic'][k], shares[nm])
            ax.annotate(txt, (100, yy), xytext=(2.0, 0), textcoords='offset points', ha='left',
                        va='center', fontsize=style.SMALL_PT, color=C.PIC_COLOUR, annotation_clip=False)
        # column headers
        ax.annotate('$E$ (eV)', (0, -0.5), xytext=(-2.0, 2.0), textcoords='offset points', ha='right', va='bottom',
                    fontsize=style.SMALL_PT, color='0.3', annotation_clip=False)
        ax.annotate('% picrate', (100, -0.5), xytext=(1.0, 2.0), textcoords='offset points', ha='left',
                    va='bottom', fontsize=style.SMALL_PT, color='0.3', annotation_clip=False)
        C.panel_title(ax, 'ab'[C.MODELS.index(m)], style.MODEL_LABELS[m], y=1.0, pad_pt=9.0)
    # orbital names once, on the far left; HOMO and LUMO in bold
    ax = axes['gas']
    for yy, nm in enumerate(names):
        ax.annotate(nm, (0, yy), xycoords=('axes fraction', 'data'),
                    xytext=(-(e_w + 0.8) / 25.4 * 72, 0), textcoords='offset points', ha='right', va='center',
                    fontsize=style.SMALL_PT, fontweight='bold' if nm in ('HOMO', 'LUMO') else 'normal',
                    annotation_clip=False)
    handles = [Patch(facecolor=C.PIC_COLOUR, edgecolor='none', label='Picrate'),
               Patch(facecolor=C.CAT_COLOUR, edgecolor='none', label='Benzotriazolium cation')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False,
               handlelength=1.6, columnspacing=2.0, fontsize=style.TICK_PT, borderaxespad=0.1)
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
