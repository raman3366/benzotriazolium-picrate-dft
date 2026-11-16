"""SI figure figS_bond_angle_deviations: DFT minus X-ray for every heavy-atom bond length and bond angle.

(a) the 26 heavy-atom bonds, (b) the 36 heavy-atom bond angles of the X-ray ion pair (btpc_cluster.xyz), gas
phase and IEF-PCM (water) (opt.xyz), sorted by the IEF-PCM deviation; grey bands +-0.02 Å and +-2 deg.
The RMSDs are asserted against facts.json (*.bond_rmsd_heavy, *.angle_rmsd_heavy) before drawing.

    python scripts/figures/figS_bond_angle_deviations.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator

from figkit import style, plots, geometry as geo
import _struct_common as C

STEM = 'figS_bond_angle_deviations'
WIDTH_MM = 175.0
ROW_MM = 2.9
# one marker per model across the paper (QA round 2, r2-m6): gas filled triangle, IEF-PCM open square
MKW = {m: C.model_marker(m, 3.3) for m in C.MODELS}


def _nat(label):
    """Natural sort key of a CIF label: element letters, then number (C1 < C2 < C10 < N1 < O1)."""
    import re
    m = re.match(r'([A-Za-z]+?)(\d+)(.*)$', label)
    return (m.group(1), int(m.group(2)), m.group(3)) if m else (label, 0, '')


def bond_name(a, b):
    return '–'.join(sorted((a, b), key=_nat))


def angle_rows(sym, labels, X, models):
    heavy = [i for i, s in enumerate(sym) if s != 'H']
    bonds = [(i, j) for i, j in geo.detect_bonds(sym, X) if i in heavy and j in heavy]
    nb = {i: [] for i in heavy}
    for i, j in bonds:
        nb[i].append(j)
        nb[j].append(i)
    rows = []
    for j in heavy:
        for a in range(len(nb[j])):
            for b in range(a + 1, len(nb[j])):
                i, k = sorted((nb[j][a], nb[j][b]), key=lambda t: _nat(labels[t]))
                r = dict(name=f'{labels[i]}–{labels[j]}–{labels[k]}', ref=geo.angle(X, i, j, k))
                for m, Y in models.items():
                    r[m] = geo.angle(Y, i, j, k)
                    r['d_' + m] = r[m] - r['ref']
                rows.append(r)
    return rows


def main():
    style.apply()
    g = C.geometries()
    sym, labels = g['symbols'], g['labels']
    models = {'gas': g['gas'], 'pcm': g['pcm']}
    brows = plots.bond_table(sym, g['xray'], models, labels, exclude_h=True)
    for r in brows:
        r['name'] = bond_name(labels[r['i']], labels[r['j']])
    arows = angle_rows(sym, labels, g['xray'], models)
    assert len(brows) == 26 and len(arows) == 36, (len(brows), len(arows))
    for m in C.MODELS:
        rb = float(np.sqrt(np.mean([r['d_' + m] ** 2 for r in brows])))
        ra = float(np.sqrt(np.mean([r['d_' + m] ** 2 for r in arows])))
        C.check_equal(rb, f'{m}.bond_rmsd_heavy', 4, f'bond RMSD ({m})')
        C.check_equal(ra, f'{m}.angle_rmsd_heavy', 3, f'angle RMSD ({m})')
    n = max(len(brows), len(arows))
    top, bottom = 14.5, 10.0
    H = top + n * ROW_MM + bottom
    fig = plt.figure(figsize=(C.mm(WIDTH_MM), C.mm(H)))
    specs = [('a', 'Bond lengths (26)', brows, 0.02, 'Å', 13.5, 68.0),
             ('b', 'Bond angles (36)', arows, 2.0, '°', 104.5, 69.5)]
    for letter, title, rows, band, unit, x0, w in specs:
        rows = sorted(rows, key=lambda r: r['d_pcm'])
        names = [r['name'] for r in rows]
        devs = {m: [r['d_' + m] for r in rows] for m in C.MODELS}
        h = n * ROW_MM              # both panels as tall as the longer list (QA round 1, m13: no empty corner)
        ax = C.mm_box(fig, x0, H - top - h, w, h)
        plots.deviation_dots(ax, names, devs, band=band, colours=style.MODEL_COLOURS,
                             labels=style.MODEL_LABELS, legend=False, marker_kw=MKW,
                             xlabel=f'DFT {C.MINUS} X-ray ({unit})')
        ax.set_yticklabels(names, fontsize=style.SMALL_PT)
        if unit == 'Å':
            ax.xaxis.set_major_locator(MultipleLocator(0.02))
        else:
            ax.xaxis.set_major_locator(MultipleLocator(2))
        rms = {m: float(np.sqrt(np.mean(np.square(devs[m])))) for m in C.MODELS}
        fmt = '{:.4f}' if unit == 'Å' else '{:.2f}'
        sp = ' ' if unit == 'Å' else ''                    # 0.0195 Å but 1.43° (no space before a degree sign)
        hs = [Line2D([], [], ls='none', label=f'{style.MODEL_LABELS[m]}, RMSD {fmt.format(rms[m])}{sp}{unit}',
                     **MKW[m]) for m in C.MODELS]
        C.header(fig, x0 - 13.5 if unit == 'Å' else x0 - 17.5, H - 3.6, letter, title)
        leg_items = hs + [Patch(color='0.9', label=f'±{band:g}{sp}{unit}')]
        ax.legend(handles=leg_items, loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=2, frameon=False,
                  fontsize=style.SMALL_PT, handletextpad=0.3, columnspacing=0.9, borderaxespad=0.2)
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
