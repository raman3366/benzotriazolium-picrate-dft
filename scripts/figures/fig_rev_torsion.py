"""fig_rev_torsion (revision round 1): relaxed gas-phase scans of the pseudo-torsion C1-O1...N4-C12.

(a) energy relative to the unconstrained minimum of each functional, B3LYP and B3LYP-D3(BJ), with RT at 298.15 K;
(b) N4...O1 along the scans; (c) angle between the picrate and cation ring planes along the scans.
The X-ray torsion is a dotted vertical line; the X-ray N4...O1 and ring-plane angle are drawn on it, the two
unconstrained minima as open stars. Every value comes from papers/btpc/data/facts.json (keys rev.struct.scan_*,
rev.struct.tors_xray_*, rev.struct.d3_gas.*, gas.*, xray.*); every number printed on the figure is asserted
against facts.json.

    python scripts/figures/fig_rev_torsion.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

from figkit import style
import _struct_common as C

STEM = 'fig_rev_torsion'
WIDTH_MM = 85.0
HEIGHT_MM = 128.0
S = 'rev.struct.'
# both scans are gas phase: the functional is told apart by colour and line, not by the model colours
SER = {'b3lyp': dict(label='B3LYP', color=style.MODEL_COLOURS['gas'], ls='-', marker='^', ms=4.2, mfc=style.MODEL_COLOURS['gas'],
                     mec='white', mew=0.3),
       'd3': dict(label='B3LYP-D3(BJ)', color='#7F2704', ls=(0, (3.2, 1.6)), marker='o', ms=3.4, mfc='white',
                  mec='#7F2704', mew=0.9)}
XRAY = C.model_marker('xray', 4.4)
STAR = dict(marker='*', ms=7.5, mfc='white', mew=0.9, ls='none', zorder=5)
VLINE = dict(color='0.45', lw=0.6, ls=(0, (1.2, 1.6)), zorder=0.5)


def series(method):
    P = S + f'scan_{method}.'
    g = np.array(C.fact(P + 'grid'), float)
    tau = np.array(C.fact(P + 'tau'), float)
    assert np.abs(tau - g).max() < 0.01, 'constrained torsions reached'
    assert C.fact(P + 'n_converged') == C.fact(P + 'n_points') == len(g)
    assert C.fact(P + 'n_proton_on_N4') == len(g)
    return dict(tau=tau, dE=np.array(C.fact(P + 'dE'), float), no=np.array(C.fact(P + 'd_N4_O1'), float),
                ring=np.array(C.fact(P + 'interplanar'), float), tau_ref=float(C.fact(P + 'tau_ref')))


def main():
    style.apply()
    D = {m: series(m) for m in ('b3lyp', 'd3')}
    # the unconstrained minima: B3LYP = published gas-phase minimum, D3(BJ) = rev.struct.d3_gas
    MIN = {'b3lyp': dict(no=float(C.fact('gas.hb_N4_H4N_O1_da')), ring=float(C.fact('gas.interplanar_picrate_cation'))),
           'd3': dict(no=float(C.fact(S + 'd3_gas.hb_N4_H4N_O1_da')), ring=float(C.fact(S + 'd3_gas.interplanar_picrate_cation')))}
    C.check_equal(D['b3lyp']['tau_ref'], S + 'tors_gas_C1_O1_N4_C12', 3)
    C.check_equal(D['d3']['tau_ref'], S + 'd3_gas.tors_C1_O1_N4_C12', 3)
    tau_x = float(C.fact(S + 'tors_xray_C1_O1_N4_C12'))
    no_x = float(C.fact('xray.hb_N4_H4N_O1_da'))
    ring_x = float(C.fact('xray.interplanar_picrate_cation'))
    rt = float(C.fact(S + 'rt_298_kcal'))

    fig = plt.figure(figsize=(C.mm(WIDTH_MM), C.mm(HEIGHT_MM)))
    xl, wa = 14.5, 67.0
    boxes = {'c': (10.5, 27.0), 'b': (46.0, 27.0), 'a': (81.5, 34.0)}
    titles = {'a': 'Energy relative to the minimum (kcal mol$^{-1}$)', 'b': 'N4···O1 (Å)',
              'c': 'Angle between the ring planes (°)'}
    ax = {}
    for k, (y0, h) in boxes.items():
        ax[k] = C.mm_box(fig, xl, y0, wa, h)
        C.header(fig, 1.0, y0 + h + 1.6, k, titles[k])
    for k, a in ax.items():
        a.set_xlim(-186, -54)
        a.xaxis.set_major_locator(MultipleLocator(30))
        a.xaxis.set_minor_locator(MultipleLocator(15))
        a.axvline(tau_x, **VLINE)
        if k != 'c':
            a.tick_params(labelbottom=False)
    ax['c'].set_xlabel('C1–O1···N4–C12 (°)')

    # (a) energy
    a = ax['a']
    a.axhline(0, color='0.75', lw=0.4, zorder=0.4)
    a.axhline(rt, color='0.55', lw=0.6, ls=(0, (4, 2)), zorder=0.6)
    for m in ('b3lyp', 'd3'):
        d = D[m]
        a.plot(d['tau'], d['dE'], lw=0.9, zorder=3, **SER[m])
        a.plot([d['tau_ref']], [0.0], mec=SER[m]['color'], **STAR)
    a.set_ylim(-0.08, 1.0)
    a.yaxis.set_major_locator(MultipleLocator(0.2))
    a.yaxis.set_minor_locator(MultipleLocator(0.1))
    a.text(-118.0, rt + 0.025, f'$RT$ (298 K) = {C.fmt_num(rt, 2)}', ha='center', va='bottom', fontsize=style.SMALL_PT,
           color='0.35')
    C.check_equal(rt, S + 'rt_298_kcal', 2)
    a.text(tau_x + 2.0, 0.95, f'X-ray\n{C.fmt_num(tau_x, 1)}°', ha='left', va='top', fontsize=style.SMALL_PT,
           color='0.3', linespacing=1.1)
    C.check_equal(tau_x, S + 'tors_xray_C1_O1_N4_C12', 1)

    # (b) N4...O1
    a = ax['b']
    for m in ('b3lyp', 'd3'):
        d = D[m]
        a.plot(d['tau'], d['no'], lw=0.9, zorder=3, **SER[m])
        a.plot([d['tau_ref']], [MIN[m]['no']], mec=SER[m]['color'], **STAR)
    a.plot([tau_x], [no_x], ls='none', zorder=4, **XRAY)
    a.set_ylim(2.44, 2.65)
    a.yaxis.set_major_locator(MultipleLocator(0.05))
    a.yaxis.set_minor_locator(MultipleLocator(0.025))
    a.text(tau_x + 4.0, no_x, f'X-ray {C.fmt_num(no_x, 3)} Å', ha='left', va='center', fontsize=style.SMALL_PT, color='0.3')
    C.check_equal(no_x, 'xray.hb_N4_H4N_O1_da', 3)

    # (c) ring planes
    a = ax['c']
    for m in ('b3lyp', 'd3'):
        d = D[m]
        a.plot(d['tau'], d['ring'], lw=0.9, zorder=3, **SER[m])
        a.plot([d['tau_ref']], [MIN[m]['ring']], mec=SER[m]['color'], **STAR)
    a.plot([tau_x], [ring_x], ls='none', zorder=4, **XRAY)
    a.set_ylim(0, 80)
    a.yaxis.set_major_locator(MultipleLocator(20))
    a.yaxis.set_minor_locator(MultipleLocator(10))
    a.annotate(f'X-ray {C.fmt_num(ring_x, 1)}°', xy=(tau_x, ring_x), xytext=(-160.5, 42.0), ha='left', va='center',
               fontsize=style.SMALL_PT, color='0.3',
               arrowprops=dict(arrowstyle='-', color='0.45', lw=0.5, shrinkA=1.0, shrinkB=2.5))
    C.check_equal(ring_x, 'xray.interplanar_picrate_cation', 1)

    handles = [Line2D([], [], lw=0.9, **SER['b3lyp']), Line2D([], [], lw=0.9, **SER['d3']),
               Line2D([], [], color='0.3', label='Unconstrained minimum', **{k: v for k, v in STAR.items() if k != 'zorder'}),
               Line2D([], [], ls='none', label='X-ray ($P2_1/c$)', **XRAY)]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(xl / WIDTH_MM - 0.01, 1.0), ncol=2, frameon=False,
               fontsize=style.TICK_PT, handlelength=2.4, handletextpad=0.4, columnspacing=1.2, borderaxespad=0.0,
               labelspacing=0.35)
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
