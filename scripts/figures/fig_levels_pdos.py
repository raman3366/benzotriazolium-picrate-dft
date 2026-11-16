"""Fig. 9 fig_levels_pdos: Kohn-Sham orbital ladder with fragment character, and fragment-resolved DOS.

One energy axis for the three panels: (a) gas-phase DOS (drawn leftwards), (b) the HOMO-7 ... LUMO+7 ladder of
both models, each level split by its Mulliken share on the picrate (purple) and the benzotriazolium cation
(green), (c) IEF-PCM (water) DOS. Dotted connectors follow the highest occupied and lowest virtual orbital of
each ion from the gas phase to IEF-PCM; the cation pi* (gas LUMO) and the picrate pi* (gas LUMO+1) cross.

Data: papers/btpc/data/derived/mo_fragment_composition_btpc_{gas,pcm}.csv (every MO; energies = fmo.json).
Asserted: HOMO/LUMO energies (fmo.json, facts), HOMO-7...LUMO+7 energies and picrate shares
(`*.orbitals_near_gap`), the ion-resolved levels (`*.cation_highest_occ_ev`, `*.cation_lowest_virt_ev`,
`*.picrate_highest_occ_ev`, `*.picrate_lowest_virt_ev`), their shifts (`diff.*`) and both gaps.

    python scripts/figures/fig_levels_pdos.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

from figkit import style, plots
import _elec_common as C

STEM = 'fig_levels_pdos'
WIDTH_MM = 175.0
YLIM = (-9.8, 0.0)              # eV; down to where both DOS have decayed (QA round 1, m10: -9.35 cut the IEF-PCM
                                # DOS at ~15 states/eV; next levels -9.40/-9.44 PCM, -10.07 gas)
FWHM = 0.3                      # eV, Gaussian broadening of the DOS
N_BELOW, N_ABOVE = 7, 7
XC = {'gas': 0.0, 'pcm': 1.25}  # ladder column centres (data units)
HW = 0.25                       # half length of a level line
SIDE_GAP = 0.09                 # gap between near-degenerate levels drawn side by side (visible; QA round 1, m10)
DEGEN = 0.07                    # eV: levels closer than this are drawn side by side (gas HOMO-5/-4 0.066 eV)
LW, LW_FRONT = 0.9, 1.7         # pt
DASH = (0, (2.2, 1.1))
GREY = '0.25'
MARK_GREY = '0.6'                # HOMO / LUMO marker lines of the DOS panels


def ion_levels(t, model):
    """Indices of the highest occupied / lowest virtual orbital of each ion (larger Mulliken share)."""
    no = t['n_occ']
    cat_major = t['cat'] > 50
    occ = np.arange(no)
    vir = np.arange(no, len(t['energy']))
    out = dict(cat_occ=int(occ[cat_major[occ]][-1]), cat_vir=int(vir[cat_major[vir]][0]),
               pic_occ=int(occ[~cat_major[occ]][-1]), pic_vir=int(vir[~cat_major[vir]][0]))
    e = t['energy']
    C.check_close(e[out['cat_occ']], f'{model}.cation_highest_occ_ev', 1e-4)
    C.check_close(e[out['cat_vir']], f'{model}.cation_lowest_virt_ev', 1e-4)
    C.check_close(e[out['pic_occ']], f'{model}.picrate_highest_occ_ev', 1e-4)
    C.check_close(e[out['pic_vir']], f'{model}.picrate_lowest_virt_ev', 1e-4)
    return out


def slots(e_sel, keys, inner_first):
    """x offsets (relative to the column centre) and half lengths of the levels; near-degenerate levels
    (< DEGEN eV apart, single linkage) share the column side by side."""
    order = np.argsort(e_sel)
    groups, cur = [], [order[0]]
    for a, b in zip(order[:-1], order[1:]):
        if e_sel[b] - e_sel[a] < DEGEN:
            cur.append(b)
        else:
            groups.append(cur)
            cur = [b]
    groups.append(cur)
    off, half = np.zeros(len(e_sel)), np.full(len(e_sel), HW)
    for g in groups:
        k = len(g)
        if k == 1:
            continue
        # connected levels (keys) go to the side facing the other column
        g = sorted(g, key=lambda i: (keys[i] in inner_first['keys']) == inner_first['right'])
        L = (2 * HW - (k - 1) * SIDE_GAP) / k
        for j, i in enumerate(g):
            off[i] = -HW + j * (L + SIDE_GAP) + L / 2
            half[i] = L / 2
    return off, half


def draw_ladder(ax, tabs, ions):
    ends = {}                                  # (model, index) -> (x_left, x_right, y)
    for m in C.MODELS:
        t = tabs[m]
        no, e = t['n_occ'], t['energy']
        sel = np.arange(no - 1 - N_BELOW, no + 1 + N_ABOVE)
        keys = list(sel)
        conn = {ions[m][k] for k in ('cat_occ', 'cat_vir', 'pic_occ', 'pic_vir')}
        off, half = slots(e[sel], keys, dict(keys=conn, right=(m == 'gas')))
        for j, k in enumerate(sel):
            pic = float(np.clip(t['pic'][k], 0, 100)) / 100
            x0 = XC[m] + off[j] - half[j]
            x1 = XC[m] + off[j] + half[j]
            xm = x0 + pic * (x1 - x0)
            y = e[k]
            front = k in (t['homo'], t['lumo'])
            kw = dict(lw=LW_FRONT if front else LW, solid_capstyle='butt', dash_capstyle='butt', zorder=3)
            if k >= no:
                kw['ls'] = DASH
            if xm - x0 > 1e-3:
                ax.plot([x0, xm], [y, y], color=C.PIC_COLOUR, **kw)
            if x1 - xm > 1e-3:
                ax.plot([xm, x1], [y, y], color=C.CAT_COLOUR, **kw)
            ends[(m, int(k))] = (x0, x1, y)
    return ends


def draw_connectors(ax, ends, ions, tabs, starts=None):
    """Dotted gas -> IEF-PCM connectors for the frontier level of each ion, with the shift in eV.  starts:
    optional {(model, index): x} where a connector starts further right (after a level name drawn between the
    level and the connector)."""
    spec = [('cat_occ', 'cation_highest_occ_ev', C.CAT_COLOUR), ('pic_occ', 'picrate_highest_occ_ev', C.PIC_COLOUR),
            ('cat_vir', 'cation_lowest_virt_ev', C.CAT_COLOUR), ('pic_vir', 'picrate_lowest_virt_ev', C.PIC_COLOUR)]
    out = {}
    starts = starts or {}
    for key, fkey, col in spec:
        g = ends[('gas', ions['gas'][key])]
        p = ends[('pcm', ions['pcm'][key])]
        xa0 = starts.get(('gas', ions['gas'][key]), g[1])
        (xa, ya), (xb, yb) = (xa0 + 0.015, g[2]), (p[0] - 0.015, p[2])
        ax.plot([xa, xb], [ya, yb], color=col, lw=0.8, ls=(0, (1.0, 1.4)), dash_capstyle='round', zorder=2)
        d = yb - ya
        C.check_equal(d, f'diff.{fkey}', 2)
        out[key] = dict(a=(xa, ya), b=(xb, yb), d=d, colour=col)
    return out


def label_connectors(ax, con):
    """Shift labels placed where the inter-column space is empty (checked by the QA overlap test)."""
    txt = lambda c: C.fmt(c['d'], 2, sign=True) + ' eV'
    kw = dict(fontsize=style.SMALL_PT, zorder=6)

    def at(c, t):
        (xa, ya), (xb, yb) = c['a'], c['b']
        return xa + t * (xb - xa), ya + t * (yb - ya)

    # cation pi: steep rise; label on the left (upper) side of the line
    x, y = at(con['cat_occ'], 0.66)
    ax.text(x - 0.03, y + 0.06, txt(con['cat_occ']), color=con['cat_occ']['colour'], ha='right', va='bottom', **kw)
    # picrate HOMO: nearly flat; label above its middle
    x, y = at(con['pic_occ'], 0.5)
    ax.text(x, y + 0.12, txt(con['pic_occ']), color=con['pic_occ']['colour'], ha='center', va='bottom', **kw)
    # the two virtual connectors cross: the cation pi* label sits above its rising line near the IEF-PCM end
    # (nothing else is there), the picrate pi* label below its falling line near the same end (empty gap region)
    x, y = at(con['cat_vir'], 0.84)
    ax.text(x - 0.04, y + 0.09, txt(con['cat_vir']), color=con['cat_vir']['colour'], ha='right', va='bottom', **kw)
    x, y = at(con['pic_vir'], 0.9)
    ax.text(x - 0.02, y - 0.14, txt(con['pic_vir']), color=con['pic_vir']['colour'], ha='right', va='top', **kw)


def level_labels(ax, ends, tabs, ions):
    """Orbital names on the outer side of each column, in the colour of the ion that carries the orbital.
    Returns {(model, index): x data} for names drawn between a level and its connector (the connector starts
    there)."""
    starts = {}
    fig = ax.figure
    for m in C.MODELS:
        t = tabs[m]
        ks = {t['homo'], t['lumo'], *ions[m].values()}
        for k in sorted(ks):
            x0, x1, y = ends[(m, k)]
            col = C.CAT_COLOUR if t['cat'][k] > 50 else C.PIC_COLOUR
            name = C.orbital_name(k, t['n_occ'])
            front = k in (t['homo'], t['lumo'])
            if m == 'gas' and x1 < XC[m] + HW - 1e-6 or m == 'gas' and x0 <= XC[m] - HW + 1e-6:
                ax.text(XC[m] - HW - 0.035, y, name, ha='right', va='center', color=col, fontsize=style.SMALL_PT,
                        fontweight='bold' if front else 'normal', zorder=6)
            elif m == 'gas':
                # a level drawn on the inner half of a side-by-side pair (gas HOMO-5 next to HOMO-4): its name
                # sits level with its own line, right of its end, and the connector to IEF-PCM starts after the
                # name (QA round 2, r2-m2: the name used to hang under the line, level with HOMO-6)
                txt = ax.text(x1 + 0.025, y, name, ha='left', va='center', color=col, fontsize=style.SMALL_PT,
                              fontweight='bold' if front else 'normal', zorder=6)
                fig.canvas.draw()
                bb = txt.get_window_extent(fig.canvas.get_renderer())
                starts[(m, k)] = float(ax.transData.inverted().transform((bb.x1, bb.y0))[0]) + 0.01
            else:
                ax.text(XC[m] + HW + 0.035, y, name, ha='left', va='center', color=col, fontsize=style.SMALL_PT,
                        fontweight='bold' if front else 'normal', zorder=6)
    return starts


def gap_arrows(ax, tabs):
    notes = {'gas': 'picrate → cation', 'pcm': 'within picrate'}
    for m in C.MODELS:
        t = tabs[m]
        eh, el = t['energy'][t['homo']], t['energy'][t['lumo']]
        gap = el - eh
        C.check_equal(gap, f'{m}.gap_ev', 2)
        # the HOMO-LUMO gap is the charge-transfer gap in the gas phase and the picrate-local gap in IEF-PCM
        C.check_equal(gap, f'{m}.ct_gap_picrate_to_cation_ev' if m == 'gas' else f'{m}.local_gap_picrate_ev', 3)
        x = XC[m] - 0.12
        ax.add_patch(FancyArrowPatch((x, eh), (x, el), arrowstyle='<|-|>', mutation_scale=6, lw=0.7, color=GREY,
                                     shrinkA=2.0, shrinkB=2.0, zorder=4))
        ym = 0.5 * (eh + el)
        ax.text(x + 0.05, ym + 0.08, f'{C.fmt(gap, 2)} eV', ha='left', va='bottom', fontsize=style.TICK_PT,
                color='0.1', zorder=6)
        ax.text(x + 0.05, ym - 0.08, notes[m], ha='left', va='top', fontsize=style.SMALL_PT, color='0.3',
                zorder=6)


def pdos(ax, t, g, mirror):
    """Stacked picrate / cation DOS (states per eV; both spins) against energy on the vertical axis."""
    e = t['energy']
    pic = np.clip(t['pic'], 0, 100)
    cat = np.clip(t['cat'], 0, 100)
    s = pic + cat
    wp, wc = pic / s, cat / s                  # clipped to 0-100 % and renormalised (diffuse virtuals only)
    s_ = FWHM / (2 * np.sqrt(2 * np.log(2)))
    G = np.exp(-0.5 * ((g[None] - e[:, None]) / s_) ** 2) / (s_ * np.sqrt(2 * np.pi))
    dp = 2 * (wp[:, None] * G).sum(0)
    dc = 2 * (wc[:, None] * G).sum(0)
    tot = 2 * plots.broaden(e, g, FWHM)
    assert np.allclose(dp + dc, tot, atol=1e-9)
    mid = 0.5 * (e[t['homo']] + e[t['lumo']])
    occ = g <= mid
    for part, alpha in ((occ, 0.55), (~occ, 0.2)):
        # include the boundary point so the fills meet
        idx = np.nonzero(part)[0]
        idx = np.unique(np.r_[idx, idx.max() + 1 if idx.max() + 1 < len(g) else idx.max()])
        ax.fill_betweenx(g[idx], 0, dp[idx], color=C.PIC_COLOUR, alpha=alpha, lw=0, zorder=2)
        ax.fill_betweenx(g[idx], dp[idx], dp[idx] + dc[idx], color=C.CAT_COLOUR, alpha=alpha, lw=0, zorder=2)
    ax.plot(dp, g, color=C.PIC_COLOUR, lw=0.6, zorder=3)
    ax.plot(dp + dc, g, color=C.CAT_COLOUR, lw=0.6, zorder=3)
    # HOMO and LUMO: thin solid grey lines named at the outer end (QA round 2, r2-m2: dotted grey lines here had
    # the style of the connector swatch in the key)
    for v, name in ((e[t['homo']], 'HOMO'), (e[t['lumo']], 'LUMO')):
        ax.axhline(v, color=MARK_GREY, lw=0.45, zorder=1)
        ax.annotate(name, (0.0 if mirror else 1.0, v), xycoords=('axes fraction', 'data'),
                    xytext=(2.0 if mirror else -2.0, 1.2), textcoords='offset points',
                    ha='left' if mirror else 'right', va='bottom', fontsize=style.SMALL_PT, color='0.35',
                    zorder=6)
    if mirror:
        ax.invert_xaxis()
    return tot.max()


def main():
    style.apply()
    tabs = {m: C.mo_table(m) for m in C.MODELS}
    ions = {m: ion_levels(tabs[m], m) for m in C.MODELS}
    # the storyline: gas LUMO is the cation pi*, IEF-PCM cation pi* is LUMO+3 and cation pi is HOMO-1
    assert ions['gas']['cat_vir'] == tabs['gas']['lumo'] and ions['pcm']['cat_vir'] == tabs['pcm']['lumo'] + 3
    assert ions['pcm']['cat_occ'] == tabs['pcm']['homo'] - 1 and ions['gas']['cat_occ'] == tabs['gas']['homo'] - 5

    # clipped Mulliken shares: report where they occur (caption)
    for m in C.MODELS:
        t = tabs[m]
        bad = np.nonzero((t['pic'] < -0.05) | (t['pic'] > 100.05))[0]
        inwin = [k for k in bad if YLIM[0] - 1.0 < t['energy'][k] < YLIM[1] + 1.0]
        print(f'{m}: shares outside 0-100 % near the window: ' +
              ', '.join(f"{C.orbital_name(k, t['n_occ'])} {t['energy'][k]:.3f} eV ({t['pic'][k]:.1f} %)"
                        for k in inwin))
        assert all(t['energy'][k] > -1.4 for k in inwin)

    W = WIDTH_MM
    H = 98.0
    fig = plt.figure(figsize=(C.mm(W), C.mm(H)))
    # millimetre layout
    left, bottom, top_strip = 9.8, 10.5, 9.5
    wd, wl, gap = 38.0, 79.6, 2.0
    right_lab = 7.0
    h = H - bottom - top_strip
    x_a = left
    x_b = x_a + wd + gap
    x_c = x_b + wl + gap
    assert abs(x_c + wd + right_lab - W) < 4.5, x_c + wd + right_lab
    rect = lambda x, w: [x / W, bottom / H, w / W, h / H]
    ax_a = fig.add_axes(rect(x_a, wd))
    ax_b = fig.add_axes(rect(x_b, wl), sharey=ax_a)
    ax_c = fig.add_axes(rect(x_c, wd), sharey=ax_a)

    g = np.linspace(YLIM[0] - 1.0, YLIM[1] + 1.0, 4000)
    mx = max(pdos(ax_a, tabs['gas'], g, mirror=True), pdos(ax_c, tabs['pcm'], g, mirror=False))
    xmax = 5 * np.ceil(mx * 1.03 / 5)
    ax_a.set_xlim(xmax, 0)
    ax_c.set_xlim(0, xmax)
    for ax in (ax_a, ax_c):
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.set_xlabel('DOS (states eV$^{-1}$)')
    ax_a.set_ylim(*YLIM)
    ax_a.yaxis.set_major_locator(MultipleLocator(1))
    ax_a.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax_a.set_ylabel('Orbital energy (eV)')
    ax_c.yaxis.tick_right()
    ax_c.yaxis.set_label_position('right')
    ax_c.tick_params(axis='y', which='both', left=False, labelleft=False, right=True, labelright=True)

    # ladder
    ax_b.set_xlim(XC['gas'] - HW - 0.47, XC['pcm'] + HW + 0.47)
    for sp in ax_b.spines.values():
        sp.set_visible(False)
    ax_b.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax_b.set_xticks([XC['gas'], XC['pcm']])
    ax_b.set_xticklabels([style.MODEL_LABELS['gas'], style.MODEL_LABELS['pcm']])
    ax_b.tick_params(axis='x', length=0, pad=4)
    ends = draw_ladder(ax_b, tabs, ions)
    starts = level_labels(ax_b, ends, tabs, ions)
    con = draw_connectors(ax_b, ends, ions, tabs, starts=starts)
    label_connectors(ax_b, con)
    gap_arrows(ax_b, tabs)

    # panel titles
    C.panel_title(ax_a, 'a', style.MODEL_LABELS['gas'], pad_pt=4.0)      # 1 mm above the spine (m10)
    C.panel_title(ax_b, 'b', 'Orbital levels, HOMO−7 to LUMO+7', pad_pt=4.0)
    C.panel_title(ax_c, 'c', style.MODEL_LABELS['pcm'], pad_pt=4.0)
    # key strip across the top
    # occupied / virtual: level line style AND the DOS shade (QA round 1, m10); the connector swatch is neutral
    # and says that its colour is that of the ion
    from matplotlib.legend_handler import HandlerTuple
    handles = [Patch(facecolor=C.PIC_COLOUR, alpha=0.7, edgecolor='none'),
               Patch(facecolor=C.CAT_COLOUR, alpha=0.7, edgecolor='none'),
               (Line2D([], [], color='0.3', lw=LW), Patch(facecolor='0.3', alpha=0.55, edgecolor='none')),
               (Line2D([], [], color='0.3', lw=LW, ls=DASH), Patch(facecolor='0.3', alpha=0.2, edgecolor='none')),
               Line2D([], [], color='0.3', lw=0.8, ls=(0, (1.0, 1.4)))]
    names = ['Picrate', 'Benzotriazolium cation', 'occupied', 'virtual',
             'same ion, gas phase → IEF-PCM (colour = ion)']
    fig.legend(handles=handles, labels=names, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=5,
               frameon=False, handlelength=2.2, columnspacing=1.3, handletextpad=0.5, fontsize=style.TICK_PT,
               borderaxespad=0.1, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
