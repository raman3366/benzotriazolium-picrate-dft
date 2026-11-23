"""fig_rev_proton: where the proton sits, by level of theory (revision, COMPTC-D-26-02710 round 1).

(a) Gas phase, one small panel per level (B3LYP, wB97X-D3(BJ), M06-2X, PBE0, CAM-B3LYP, DF-MP2, SCS-MP2): the
    energy along the relaxed B3LYP scan of r(H4N...O1) (16 points; the other levels are single points at those
    geometries), the interpolating spline used by the one-dimensional DVR, the proton levels v = 0 and v = 1 and the
    ground-state density |psi_0|^2 (scaled to a common height). Above each panel: depth of the ionic well (or none)
    and the O-side probability of v = 0.
(b) IEF-PCM (water): E(neutral pair) - E(salt) at the two B3LYP IEF-PCM minima for every level, and the CCSD(T)
    estimate (DF-MP2 + FNO-CCSD(T) correction at 80 % of the virtual space) with the range of its extrapolations.

Inputs: papers/btpc/data/facts.json (rev.proton.*) and papers/btpc/data/derived/rev_proton_dvr_{grid,levels}.csv,
both written by papers/btpc/data/build_facts.py (module rev/proton.py). Every printed number is asserted against
facts.json; the curves are checked against the profile points of facts.json.

    python scripts/figures/fig_rev_proton.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator

from figkit import style
import _struct_common as C

STEM = 'fig_rev_proton'
WIDTH_MM, HEIGHT_MM = 175.0, 150.0
DERIVED = C.PAPER / 'data' / 'derived'
LEVELS = [('b3lyp', 'B3LYP'), ('wb97xd3bj', 'ωB97X-D3(BJ)'), ('m062x', 'M06-2X'), ('pbe0', 'PBE0'),
          ('camb3lyp', 'CAM-B3LYP'), ('dfmp2', 'DF-MP2'), ('scsmp2', 'SCS-MP2')]
GAS = style.MODEL_COLOURS['gas']
PCM = style.MODEL_COLOURS['pcm']
INK2 = '0.30'
XLIM = (0.86, 1.74)
YLIM = (-0.3, 8.4)
PSI_H = 1.5                     # common display height of |psi_0|^2 (kcal/mol units of the axis)
P = 'rev.proton.'


def F(key):
    return C.fact(P + key)


def fmt(v, nd, sign=False):
    s = f'{v:+.{nd}f}' if sign else f'{v:.{nd}f}'
    return s.replace('-', C.MINUS)


def read_csv(name):
    with open(DERIVED / name) as f:
        return list(csv.DictReader(f))


def load():
    grid = read_csv('rev_proton_dvr_grid.csv')
    g = {k: np.array([float(r[k]) for r in grid]) for k in grid[0]}
    lev = {}
    for r in read_csv('rev_proton_dvr_levels.csv'):
        lev.setdefault(r['method'], {})[int(r['n'])] = float(r['E_kcal'])
    return g, lev


def segments(x, y, mask):
    """Line segments of y where mask holds (NaN elsewhere), for levels drawn only in the allowed region."""
    out = np.where(mask, y, np.nan)
    return x, out


def panel_a(fig, ax, k, name, g, lev, r_pts, il, first_col):
    V = g[f'V_kcal_{k}']
    psi = g[f'psi2_v0_{k}']
    r = g['r_HO_A']
    E0, E1 = lev[k][0], lev[k][1]
    C.check_equal(E0, f'{P}dvr.{k}.zpe_kcal', 4, f'{k}: v = 0 level in the CSV')
    C.check_equal(E1, f'{P}dvr.{k}.v1_kcal', 4, f'{k}: v = 1 level in the CSV')
    # computed points, on the DVR zero (minimum of the spline): dE relative to the lowest point + spline undershoot
    dE = np.array(F(f'gas.{k}.dE_kcal'))
    vp = dE - dE[il] + float(F(f'dvr.{k}.spline_undershoot_kcal'))
    x_grid, x_pts = g['x_A'], np.array(F('scan.x_A'))
    v_at = np.interp(x_pts, x_grid, V)
    assert np.abs(v_at - vp).max() < 0.02, (k, np.abs(v_at - vp).max())     # the spline passes through the points
    ax.axvline(float(F('scan.transfer_r_HO_A')), color='0.6', lw=0.5, ls=(0, (3, 2)), zorder=1)
    ax.plot(r, V, color=GAS, lw=1.0, zorder=3)
    ax.plot(r_pts, vp, ls='none', zorder=4, **C.model_marker('gas', 3.0))
    m0 = V <= E0
    ax.plot(*segments(r, np.full_like(r, E0), m0), color='0.15', lw=0.8, zorder=5)
    ax.plot(*segments(r, np.full_like(r, E1), V <= E1), color='0.35', lw=0.6, ls=(0, (2.5, 1.5)), zorder=2)
    sel = psi > 4e-3 * psi.max()          # no hairline tail along the level
    ax.fill_between(r[sel], E0, E0 + PSI_H * psi[sel] / psi.max(), color=GAS, alpha=0.30, lw=0, zorder=2)
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    if not first_col:
        ax.tick_params(labelleft=False)
    # title block above the panel: level, ionic well, P(O side)
    pO = C.check_equal(float(F(f'dvr.{k}.p_O_v0_pct')), f'{P}dvr.{k}.p_O_v0_pct', 1)
    if F(f'gas.{k}.ionic_well') == 'yes':
        well = C.check_equal(float(F(f'gas.{k}.well_depth_kcal')), f'{P}gas.{k}.well_depth_kcal', 2)
        wtxt = f'ionic well {fmt(well, 2)}'
    else:
        assert F(f'gas.{k}.ionic_well') == 'no'
        wtxt = 'no ionic well'
    bb = ax.get_position()
    fig.text(bb.x0, bb.y1 + C.mm(3.6) / fig.get_size_inches()[1], name, fontsize=style.TICK_PT,
             ha='left', va='bottom')
    fig.text(bb.x0, bb.y1 + C.mm(0.8) / fig.get_size_inches()[1],
             f'{wtxt}; $P$(O side) {fmt(pO, 1)} %', fontsize=style.SMALL_PT, color=INK2, ha='left', va='bottom')
    return E0


def panel_b(ax):
    rows = [(name, float(F(f'pcm.{k}.dE_neutral_minus_salt_kcal'))) for k, name in LEVELS]
    best = float(F('cc.best_kcal'))
    lo, hi = [float(v) for v in F('cc.best_range_kcal')]
    n = len(rows) + 1
    ys = np.arange(n)[::-1]
    ax.axvline(0, color='0.35', lw=0.6, zorder=1)
    for (name, v), y in zip(rows, ys[:-1]):
        key = next(k for k, nm in LEVELS if nm == name)
        C.check_equal(v, f'{P}pcm.{key}.dE_neutral_minus_salt_kcal', 2)
        ax.plot([0, v], [y, y], color=PCM, lw=0.6, alpha=0.5, zorder=2)
        ax.plot([v], [y], ls='none', zorder=4, **C.model_marker('pcm', 4.2))
        ax.annotate(fmt(v, 2, sign=True), (v, y), xytext=(5 if v > 0 else -5, 0), textcoords='offset points',
                    ha='left' if v > 0 else 'right', va='center', fontsize=style.SMALL_PT)
    y = ys[-1]
    C.check_equal(best, f'{P}cc.best_kcal', 2)
    assert lo <= best <= hi
    ax.plot([lo, hi], [y, y], color='0.1', lw=1.6, solid_capstyle='butt', zorder=3)
    for xe in (lo, hi):
        ax.plot([xe, xe], [y - 0.18, y + 0.18], color='0.1', lw=0.8, zorder=3)
    ax.plot([0, lo], [y, y], color='0.1', lw=0.6, alpha=0.5, zorder=2)
    ax.plot([best], [y], ls='none', marker='s', ms=4.2, mfc='0.1', mec='white', mew=0.4, zorder=4)
    ax.annotate(f'{fmt(best, 2, sign=True)} ({fmt(lo, 2)} to {fmt(hi, 2)})', (lo, y), xytext=(-5, 0),
                textcoords='offset points', ha='right', va='center', fontsize=style.SMALL_PT)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows] + ['DF-MP2 + Δ[CCSD(T)]'])
    ax.tick_params(axis='y', length=0)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(-4.8, 5.4)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.set_xlabel('$E$(neutral pair) − $E$(salt) (kcal mol$^{-1}$)')
    ax.text(-0.15, n - 0.55, '← neutral pair lower', ha='right', va='center', fontsize=style.SMALL_PT, color=INK2)
    ax.text(0.15, n - 0.55, 'salt lower →', ha='left', va='center', fontsize=style.SMALL_PT, color=INK2)
    ax.grid(axis='x', color='0.9', lw=0.4, zorder=0)
    ax.set_axisbelow(True)


def main():
    style.apply()
    g, lev = load()
    n = int(F('scan.n_points'))
    r_pts = np.array(F('scan.r_HO_A'))
    assert len(r_pts) == n
    il = int(np.argmin(np.abs(r_pts - float(F('gas.lowest_point_r_HO_A')))))
    fig = plt.figure(figsize=(C.mm(WIDTH_MM), C.mm(HEIGHT_MM)))
    H = HEIGHT_MM
    # ---- (a) small multiples ----
    C.header(fig, 0.0, H - 4.0, 'a', 'Gas phase: energy along the proton scan, proton levels and ground-state density')
    x0, w, gap = 15.0, 36.6, 3.6
    ah, top1, top2 = 25.0, H - 14.0, H - 52.0
    axes = []
    for i, (k, name) in enumerate(LEVELS):
        row, col = divmod(i, 4)
        top = top1 if row == 0 else top2
        ax = C.mm_box(fig, x0 + col * (w + gap), top - ah, w, ah)
        panel_a(fig, ax, k, name, g, lev, r_pts, il, col == 0)
        axes.append(ax)
    fig.text(C.mm(3.0) / C.mm(WIDTH_MM), (top1 - ah / 2 + top2 - ah / 2) / 2 / H, 'Energy (kcal mol$^{-1}$)',
             rotation=90, ha='center', va='center', fontsize=style.BASE_PT)
    xl = x0 + 1.5 * (w + gap) + w / 2
    fig.text(xl / WIDTH_MM, (top2 - ah - 9.0) / H, '$r$(H4N···O1) (Å)', ha='center', va='baseline',
             fontsize=style.BASE_PT)
    # key in the eighth cell
    kx = x0 + 3 * (w + gap)
    handles = [Line2D([], [], color=GAS, lw=1.0, label='Spline through the points',
                      **{k_: v_ for k_, v_ in C.model_marker('gas', 3.0).items() if k_ != 'color'}),
               Line2D([], [], color='0.15', lw=0.8, label='Proton level $v$ = 0'),
               Line2D([], [], color='0.35', lw=0.6, ls=(0, (2.5, 1.5)), label='Proton level $v$ = 1'),
               Patch(facecolor=GAS, alpha=0.30, lw=0, label='|$\\psi_0$|$^2$ (scaled)'),
               Line2D([], [], color='0.6', lw=0.5, ls=(0, (3, 2)), label='$r$(N4–H4N) = $r$(H4N···O1)')]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(kx / WIDTH_MM, (top2 + 4.0) / H), frameon=False,
               fontsize=style.SMALL_PT, handlelength=2.0, handletextpad=0.5, labelspacing=0.55, borderaxespad=0.0)
    fig.text((kx + 1.0) / WIDTH_MM, (top2 - ah - 1.0) / H, 'Left: proton on O1 (neutral)\nRight: proton on N4 (salt)\n'
             'Ionic well in kcal mol$^{-1}$', ha='left', va='bottom', fontsize=style.SMALL_PT, color=INK2,
             linespacing=1.3)
    # ---- (b) ----
    yb = top2 - ah - 18.0
    C.header(fig, 0.0, yb, 'b', 'IEF-PCM (water): neutral pair relative to the salt, at their B3LYP minima')
    bh = 40.0
    ax_b = C.mm_box(fig, 38.0, yb - 4.0 - bh, WIDTH_MM - 38.0 - 4.0, bh)
    panel_b(ax_b)
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
