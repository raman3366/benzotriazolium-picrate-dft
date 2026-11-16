"""Fig. 3 fig_structures: optimised ion pair in the gas phase (a) and in IEF-PCM (water) (b).

Same camera for both panels (the paper's common orientation, `_struct_common.views`: picrate on the left in a
common position, cation on the right), one print scale (mm per Å) for both, CIF labels on every heavy atom plus
H4N, H6N and H11, hydrogen bonds with H···A distances (3 decimals) computed from opt.xyz and checked against
facts.json before saving.

    python scripts/figures/fig_structures.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt

from figkit import style, draw2d, geometry as geo
import _struct_common as C

STEM = 'fig_structures'
WIDTH_MM = 175.0
MARGIN_A = 0.8           # room for labels around the atoms (Å)
GAP_MM = 3.0


def build(s, sym, labels, V, hbonds):
    boxes = [C.extent(sym, V[m], MARGIN_A) for m in C.MODELS]
    y0 = min(b[2] for b in boxes)
    y1 = max(b[3] for b in boxes)
    boxes = [np.array([b[0], b[1], y0, y1]) for b in boxes]
    wA = sum(b[1] - b[0] for b in boxes)
    W = wA * s + GAP_MM
    H = (y1 - y0) * s + 8.0
    fig = plt.figure(figsize=(C.mm(W), C.mm(H)))
    axes, s, h = C.molecule_row(fig, boxes, 0.0, 0.0, W, gap_mm=GAP_MM, scale_mm_per_A=s)
    label_atoms = [i for i, x in enumerate(sym) if x != 'H'] + [C.I['H4N'], C.I['H6N'], C.I['H11']]
    for m, ax in zip(C.MODELS, axes):
        ax.patch.set_visible(False)
        d = draw2d.draw_molecule(ax, sym, V[m], labels, hbonds=hbonds[m], label_atoms=label_atoms, label_h=True,
                                 hbond_distances=True, hbond_fmt='{:.3f} Å', fit=False, strict=True)
        assert not d.dropped, d.dropped_labels()
    C.crop_row(fig, axes, gap_mm=GAP_MM, top_mm=5.0)
    r = fig.canvas.get_renderer()
    bbs = [C.content_bbox(ax, r) for ax in axes]
    top = max(b.y1 for b in bbs) + 2.0 * fig.dpi / 72
    for k, (m, bb) in enumerate(zip(C.MODELS, bbs)):
        C.panel_title(fig, bb.x0, top, 'ab'[k], style.MODEL_LABELS[m])
    return fig


def main():
    style.apply()
    g = C.geometries()
    sym, labels = g['symbols'], g['labels']
    hbonds = {'gas': [C.HB, (C.I['C11'], C.I['H11'], C.I['O2'])], 'pcm': [C.HB]}
    # the numbers printed on the figure, checked against facts.json at the printed precision
    for m in C.MODELS:
        C.check_equal(geo.distance(g[m], 31, 0), f'{m}.hb_N4_H4N_O1_ha', 3)
    C.check_equal(geo.distance(g['gas'], 27, 15), 'gas.hb_C11_H11_O2_ha', 3)
    V = {m: C.oriented(m) for m in C.MODELS}
    fig, s = C.fit_to_width(lambda s: build(s, sym, labels, V, hbonds), WIDTH_MM, 7.0)
    print(f'scale {s:.3f} mm per Å')
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
