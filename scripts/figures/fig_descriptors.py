"""Fig. 10 fig_descriptors: gas phase versus IEF-PCM (water), frontier-orbital descriptors and charge separation.

(a) E_HOMO, E_LUMO, gap, chi, eta, omega (eV); (b) global softness S = 1/(2 eta) (eV^-1); (c) dipole moment (D);
(d) Mulliken charge summed over the picrate (e). Horizontal bars so that every value label reads horizontally.
Values from facts.json (`*.homo_ev`, `*.lumo_ev`, `*.gap_ev`, `*.chi_ev`, `*.eta_ev`, `*.omega_ev`,
`*.softness_ev_inv`, `*.dipole_debye`, `*.mulliken_charge_picrate`), cross-checked against
molecules/btpc/results/<case>/fmo.json; the gap character (charge transfer in the gas phase, picrate-local in
IEF-PCM) is asserted from `*.ct_gap_picrate_to_cation_ev` / `*.local_gap_picrate_ev` and the LUMO composition.

    python scripts/figures/fig_descriptors.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator

from figkit import style
import _elec_common as C

STEM = 'fig_descriptors'
WIDTH_MM = 175.0
BAR = 0.36                      # bar thickness (row pitch 1)
OFF = {'gas': -0.2, 'pcm': 0.2}  # gas bar above, IEF-PCM bar below (y axis inverted)

ENERGY = [('homo_ev', r'$E_\mathrm{HOMO}$', 'HOMO_eV'),
          ('lumo_ev', r'$E_\mathrm{LUMO}$', 'LUMO_eV'),
          ('gap_ev', r'$\Delta E_\mathrm{gap}$', 'gap_eV'),
          ('chi_ev', r'$\chi$', 'electronegativity_chi_eV'),
          ('eta_ev', r'$\eta$', 'hardness_eta_eV'),
          ('omega_ev', r'$\omega$', 'electrophilicity_omega_eV')]
# key, title, decimals, fmo.json key, (tick step, limit step)
SMALL = [('softness_ev_inv', 'Global softness $S$ (eV$^{-1}$)', 3, 'softness_S_eV-1', (0.1, 0.05)),
         ('dipole_debye', 'Dipole moment (D)', 2, 'dipole_Debye', (10, 5)),
         ('mulliken_charge_picrate', 'Mulliken charge on the picrate (e)', 3, None, (0.5, 0.25))]
TITLE_MM, TICKS_MM, LEGEND_MM, XLAB_MM = 5.0, 4.6, 5.2, 10.0


def values():
    fj = {m: json.loads((C.REPO / 'molecules' / 'btpc' / 'results' / C.CASES[m] / 'fmo.json').read_text())
          for m in C.MODELS}
    v = {}
    for m in C.MODELS:
        for key, _, jk in ENERGY:
            v[(m, key)] = C.check_close(fj[m][jk], f'{m}.{key}', 1e-9)
        for key, _, _, jk, _ in SMALL:
            f = C.ffact(f'{m}.{key}')
            if jk:
                C.check_close(fj[m][jk], f'{m}.{key}', 1e-9)
            v[(m, key)] = f
    # gap character, as annotated
    assert abs(v[('gas', 'gap_ev')] - C.ffact('gas.ct_gap_picrate_to_cation_ev')) < 1e-9
    assert abs(v[('pcm', 'gap_ev')] - C.ffact('pcm.local_gap_picrate_ev')) < 1e-9
    assert C.ffact('gas.lumo_picrate_pct') < 1 and C.ffact('pcm.lumo_picrate_pct') > 99
    assert C.ffact('gas.homo_picrate_pct') > 99 and C.ffact('pcm.homo_picrate_pct') > 99
    return v


def bars(ax, rows, v, nd, notes=None):
    """rows: list of facts keys (one row each); returns the value-label artists."""
    arts = []
    for r, key in enumerate(rows):
        for m in C.MODELS:
            x = v[(m, key)]
            y = r + OFF[m]
            ax.barh(y, x, height=BAR, color=style.MODEL_COLOURS[m], lw=0, zorder=2)
            n = nd[key] if isinstance(nd, dict) else nd
            s = C.fmt(x, n)
            if notes and (m, key) in notes:
                s += f'  {notes[(m, key)]}'
            arts.append(ax.annotate(s, (x, y), xytext=(2.0 if x >= 0 else -2.0, 0), textcoords='offset points',
                                    ha='left' if x >= 0 else 'right', va='center', fontsize=style.SMALL_PT,
                                    color='0.1', zorder=3))
    ax.axvline(0, color='0.3', lw=0.6, zorder=3)
    ax.set_ylim(len(rows) - 0.45, -0.55)
    ax.tick_params(axis='y', length=0)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    return arts


def main():
    style.apply()
    v = values()
    # layout (mm): panel (a) on the left; (b)-(d) stacked on the right with their bars as thick as those of (a)
    # and the column spanning the height of (a), the title of (b) level with that of (a):
    #   3 ph + 2 TITLE + 2 TICKS = ha,  ph = 1.1 ha / 6.1
    na = len(ENERGY)
    span_a, span_k = na + 0.1, 1.1                      # y-axis spans (row pitch 1)
    ha = (2 * TITLE_MM + 2 * TICKS_MM) / (1 - 3 * span_k / span_a)
    ph = span_k * ha / span_a
    W = WIDTH_MM
    H = LEGEND_MM + TITLE_MM + ha + XLAB_MM
    fig = plt.figure(figsize=(C.mm(W), C.mm(H)))
    la, wa, ba = 10.2, 96.4, XLAB_MM
    ax = fig.add_axes([la / W, ba / H, wa / W, ha / H])
    notes = {('gas', 'gap_ev'): '(picrate → cation)', ('pcm', 'gap_ev'): '(within picrate)'}
    bars(ax, [k for k, _, _ in ENERGY], v, 2, notes)
    ax.set_yticks(range(na))
    ax.set_yticklabels([lab for _, lab, _ in ENERGY], fontsize=style.BASE_PT)
    ax.set_xlim(-8.8, 9.4)
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_xlabel('Energy (eV)')
    for r in range(na - 1):
        ax.axhline(r + 0.5, color='0.9', lw=0.5, zorder=0)
    C.panel_title(ax, 'a', 'Frontier-orbital energies and descriptors')

    lx = la + wa + 9.0
    wx = W - lx - 2.3                                   # saved width <= 175 mm (QA round 1, m5)
    for k, (key, title, nd, _, (tick, lim)) in enumerate(SMALL):
        yb = ba + (2 - k) * (ph + TITLE_MM + TICKS_MM)
        axk = fig.add_axes([lx / W, yb / H, wx / W, ph / H])
        bars(axk, [key], v, nd)
        axk.set_yticks([])
        axk.spines['left'].set_visible(False)
        vals = [v[(m, key)] for m in C.MODELS]
        if min(vals) < 0:
            axk.set_xlim(lim * np.floor(min(vals) * 1.17 / lim), 0)
        else:
            axk.set_xlim(0, lim * np.ceil(max(vals) * 1.17 / lim))
        axk.xaxis.set_major_locator(MultipleLocator(tick))
        axk.xaxis.set_minor_locator(MultipleLocator(tick / 2))
        C.panel_title(axk, 'bcd'[k], title)

    handles = [Patch(facecolor=style.MODEL_COLOURS[m], edgecolor='none', label=style.MODEL_LABELS[m])
               for m in C.MODELS]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False,
               handlelength=1.6, columnspacing=2.0, fontsize=style.TICK_PT, borderaxespad=0.1)
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
