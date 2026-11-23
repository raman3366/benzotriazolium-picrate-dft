"""fig_rev_eps (revision round 1): where the LUMO of the ion pair lives as the permittivity of the continuum grows.

B3LYP/6-311++G(d,p) single points at the two published P2_1/c minima (gas-phase minimum, IEF-PCM minimum), the
permittivity eps of the IEF-PCM continuum varied at fixed geometry (eps = 1: no continuum; 78.3553: water).
(a) E_HOMO and E_LUMO versus eps (log axis); open markers: LUMO on the cation, filled: LUMO on the picrate.
(b) Mulliken share of the LUMO on the picrate atoms (%).
(c) The two acceptor levels: lowest cation-localised and lowest picrate-localised virtual orbital; the LUMO is the
    lower of the two, so the LUMO changes ion where they cross.
(d) Delta = E(cation level) - E(picrate level) against f = (eps - 1)/eps with the least-squares lines.
Vertical dotted lines: the interpolated switch eps of each geometry (facts rev.elec.eps.<geom>.switch_eps).
If facts.json holds a key starting 'lit.eps_crystal' (a measured crystal permittivity, e.g. from ref. [12]), it is
marked in (a)-(c); no such value is available yet, and the script says so.

Data: papers/btpc/data/derived/rev_elec_eps_series.csv (written by papers/btpc/data/rev/electronic.py); every plotted
point is checked against facts.json (rev.elec.eps.*) and every printed number is asserted against facts.json.

    python scripts/figures/fig_rev_eps.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullLocator, MultipleLocator

from figkit import style
import _elec_common as C

STEM = 'fig_rev_eps'
WIDTH_MM = 175.0
GEOMS = ('gasgeom', 'pcmgeom')
GEO_LABEL = {'gasgeom': 'Gas-phase minimum', 'pcmgeom': 'IEF-PCM minimum'}
GEO_COLOUR = {'gasgeom': style.MODEL_COLOURS['gas'], 'pcmgeom': style.MODEL_COLOURS['pcm']}
GEO_MARKER = {'gasgeom': 'o', 'pcmgeom': 's'}
CSV_PATH = C.DERIVED / 'rev_elec_eps_series.csv'


def etag(eps):
    return 'e' + f'{eps:g}'.replace('.', 'p')


def load():
    rows = [r for r in csv.DictReader(CSV_PATH.open()) if r['series'] in GEOMS]
    out = {}
    for g in GEOMS:
        rs = sorted((r for r in rows if r['series'] == g), key=lambda r: float(r['eps']))
        d = {k: np.array([float(r[c]) for r in rs]) for k, c in
             (('eps', 'eps'), ('f', 'f_(eps-1)/eps'), ('homo', 'HOMO_eV'), ('lumo', 'LUMO_eV'),
              ('pic', 'LUMO_picrate_pct'), ('cat_virt', 'lowest_cation_virtual_eV'),
              ('pic_virt', 'lowest_picrate_virtual_eV'), ('delta', 'delta_cat_minus_pic_eV'))}
        d['on'] = [r['LUMO_on'] for r in rs]
        assert len(rs) == C.fact(f'rev.elec.eps.{g}.n_points'), g
        for k, e in enumerate(d['eps']):
            key = f'rev.elec.eps.{g}.{etag(e)}'
            for q, fk in (('homo', 'homo_ev'), ('lumo', 'lumo_ev'), ('pic', 'lumo_picrate_pct'),
                          ('delta', 'delta_cat_pic_ev')):
                C.check_close(d[q][k], f'{key}.{fk}', 1e-8)
            assert d['on'][k] == C.fact(f'{key}.lumo_on'), (g, e)
            assert abs(d['f'][k] - (e - 1) / e) < 1e-12
            assert abs(min(d['cat_virt'][k], d['pic_virt'][k]) - d['lumo'][k]) < 1e-9, (g, e)
        out[g] = d
    return out


def switch_lines(ax, sw, ymax_frac=1.0, labels=False, y_text=None):
    arts = []
    for g in GEOMS:
        ax.axvline(sw[g], color=GEO_COLOUR[g], lw=0.8, ls=':', zorder=1, ymax=ymax_frac)
        if labels:
            arts.append(ax.annotate(f'{C.fmt(sw[g], 2)}', (sw[g], y_text), xytext=(-2, 0), textcoords='offset points',
                                    ha='right', va='center', fontsize=style.SMALL_PT, color='0.15'))
    return arts


def eps_axis(ax, wat, xlabel=True):
    ax.set_xscale('log')
    ticks = [1, 2, 4, 8, 20, wat]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels([f'{t:g}' if t != wat else f'{wat:.1f}' for t in ticks])
    ax.set_xlim(0.88, wat * 1.25)
    if xlabel:
        ax.set_xlabel('Permittivity of the continuum, $\\varepsilon$')
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


def markers(ax, g, x, y, on, size=4.2, z=4):
    for xi, yi, o in zip(x, y, on):
        filled = o == 'picrate'
        ax.plot([xi], [yi], marker=GEO_MARKER[g], ms=size, mec=GEO_COLOUR[g], mew=0.9,
                mfc=GEO_COLOUR[g] if filled else 'white', ls='none', zorder=z)


def crystal_marks(axes):
    keys = [k for k in C.facts() if k.startswith('lit.eps_crystal')]
    if not keys:
        print('NOTE: no measured crystal permittivity in facts.json (lit.eps_crystal*); none is marked')
        return
    for k in keys:
        v = C.ffact(k)
        for ax in axes:
            ax.axvline(v, color='0.35', lw=1.2, ls='-', alpha=0.35, zorder=0)


def main():
    style.apply()
    D = load()
    wat = C.ffact('rev.elec.eps.water')
    sw = {g: C.ffact(f'rev.elec.eps.{g}.switch_eps') for g in GEOMS}
    for g in GEOMS:
        assert C.fact(f'rev.elec.eps.{g}.switch') == 'cation->picrate'
        assert C.ffact(f'rev.elec.eps.{g}.switch_eps_lo') < sw[g] < C.ffact(f'rev.elec.eps.{g}.switch_eps_hi')
        assert D[g]['eps'][-1] == wat and D[g]['eps'][0] == 1.0

    W = WIDTH_MM
    LEG, TITLE, PH, XLAB, L, PW, GAP, R = 7.0, 5.0, 47.0, 11.0, 11.5, 69.5, 21.5, 3.0
    assert abs(L + 2 * PW + GAP + R - W) < 1e-9
    H = LEG + 2 * (TITLE + PH + XLAB)
    fig = plt.figure(figsize=(C.mm(W), C.mm(H)))

    def ax_at(col, row):
        x = L + col * (PW + GAP)
        y = XLAB + (1 - row) * (TITLE + PH + XLAB)
        return fig.add_axes([x / W, y / H, PW / W, PH / H])

    # ---------------------------------------------------------------- (a) HOMO and LUMO
    ax = ax_at(0, 0)
    for g in GEOMS:
        d = D[g]
        ax.plot(d['eps'], d['lumo'], color=GEO_COLOUR[g], lw=1.1, zorder=3)
        markers(ax, g, d['eps'], d['lumo'], d['on'])
        ax.plot(d['eps'], d['homo'], color=GEO_COLOUR[g], lw=1.1, ls=(0, (4, 2)), zorder=3)
        markers(ax, g, d['eps'], d['homo'], ['picrate'] * len(d['eps']), size=3.2)
    eps_axis(ax, wat)
    ax.set_ylim(-7.3, -2.7)
    ax.yaxis.set_major_locator(MultipleLocator(1.0))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.set_ylabel('Orbital energy (eV)')
    switch_lines(ax, sw)
    ax.text(wat * 1.12, D['gasgeom']['lumo'][-1] + 0.32, 'LUMO', ha='right', va='bottom', fontsize=style.SMALL_PT)
    ax.text(wat * 1.12, D['gasgeom']['homo'][-1] + 0.25, 'HOMO', ha='right', va='bottom', fontsize=style.SMALL_PT)
    C.panel_title(ax, 'a', 'Frontier orbitals')
    axa = ax

    # ---------------------------------------------------------------- (b) picrate share of the LUMO
    ax = ax_at(1, 0)
    for g in GEOMS:
        d = D[g]
        ax.plot(d['eps'], d['pic'], color=GEO_COLOUR[g], lw=1.1, zorder=3)
        markers(ax, g, d['eps'], d['pic'], d['on'])
    eps_axis(ax, wat)
    ax.set_ylim(-5, 112)
    ax.yaxis.set_major_locator(MultipleLocator(25))
    ax.set_ylabel('LUMO on the picrate (%)')
    switch_lines(ax, sw, labels=True, y_text=75)
    for g in GEOMS:
        C.check_equal(sw[g], f'rev.elec.eps.{g}.switch_eps', 2, 'printed switch eps')
    C.panel_title(ax, 'b', 'Picrate share of the LUMO')
    axb = ax

    # ---------------------------------------------------------------- (c) the two acceptor levels
    ax = ax_at(0, 1)
    for g in GEOMS:
        d = D[g]
        ax.plot(d['eps'], d['cat_virt'], color=GEO_COLOUR[g], lw=1.1, ls=(0, (1.2, 1.4)), marker=GEO_MARKER[g],
                ms=3.4, mfc='white', mec=GEO_COLOUR[g], mew=0.8, zorder=3)
        ax.plot(d['eps'], d['pic_virt'], color=GEO_COLOUR[g], lw=1.1, marker=GEO_MARKER[g], ms=3.4,
                mfc=GEO_COLOUR[g], mec=GEO_COLOUR[g], zorder=3)
    eps_axis(ax, wat)
    ax.set_ylim(-4.0, -2.2)
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.25))
    ax.set_ylabel('Orbital energy (eV)')
    switch_lines(ax, sw)
    ax.text(wat * 1.12, D['gasgeom']['cat_virt'][-1] + 0.07, 'cation level', ha='right', va='bottom',
            fontsize=style.SMALL_PT)
    ax.text(wat * 1.12, D['gasgeom']['pic_virt'][-1] - 0.12, 'picrate level', ha='right', va='top',
            fontsize=style.SMALL_PT)
    C.panel_title(ax, 'c', 'Lowest virtual orbital on each ion')
    axc = ax

    # ---------------------------------------------------------------- (d) Delta against f
    ax = ax_at(1, 1)
    fgrid = np.linspace(0, 1, 50)
    for g in GEOMS:
        d = D[g]
        b, a = C.ffact(f'rev.elec.eps.{g}.fit_slope_ev'), C.ffact(f'rev.elec.eps.{g}.fit_intercept_ev')
        bb, aa = np.polyfit(d['f'], d['delta'], 1)
        assert abs(bb - b) < 1e-9 and abs(aa - a) < 1e-9
        ax.plot(fgrid, a + b * fgrid, color=GEO_COLOUR[g], lw=0.8, zorder=2)
        markers(ax, g, d['f'], d['delta'], d['on'])
    ax.axhline(0, color='0.45', lw=0.6, zorder=1)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-1.45, 1.15)
    ax.xaxis.set_major_locator(MultipleLocator(0.25))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.set_xlabel('$f = (\\varepsilon - 1)/\\varepsilon$')
    ax.set_ylabel('$\\Delta$ = $E_\\mathrm{cation}$ $-$ $E_\\mathrm{picrate}$ (eV)')
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ytxt = {'gasgeom': 0.95, 'pcmgeom': 0.80}
    for g in GEOMS:
        r2 = C.ffact(f'rev.elec.eps.{g}.fit_r2')
        C.check_equal(r2, f'rev.elec.eps.{g}.fit_r2', 3, 'printed R2')
        ax.text(0.03, ytxt[g], f'$R^2$ = {r2:.3f}', transform=ax.transAxes, ha='left', va='top',
                fontsize=style.SMALL_PT, color=GEO_COLOUR[g])
    ax.text(0.97, 0.03, 'LUMO on the cation below the line,\non the picrate above it', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=style.SMALL_PT, color='0.3')
    C.panel_title(ax, 'd', 'Level ordering against the Pekar factor')

    crystal_marks([axa, axb, axc])

    handles = [Line2D([], [], color=GEO_COLOUR[g], marker=GEO_MARKER[g], mfc=GEO_COLOUR[g], mec=GEO_COLOUR[g],
                      lw=1.1, ms=4, label=GEO_LABEL[g]) for g in GEOMS]
    handles += [Line2D([], [], color='0.35', marker='o', mfc='white', mec='0.35', lw=0, ms=4,
                       label='open: on the cation'),
                Line2D([], [], color='0.35', marker='o', mfc='0.35', mec='0.35', lw=0, ms=4,
                       label='filled: on the picrate'),
                Line2D([], [], color='0.35', lw=0.8, ls=':', label='switch $\\varepsilon$')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=5, frameon=False,
               handlelength=1.8, columnspacing=1.4, fontsize=style.TICK_PT, borderaxespad=0.1)
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
