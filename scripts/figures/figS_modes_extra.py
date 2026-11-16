"""SI figure figS_modes_extra: displacement vectors of further N-H modes (method of Fig. 6).

Columns: N4-H4N in-plane bend, N4-H4N out-of-plane wag, free N6-H6N stretch; rows: gas phase (top, as in every
other figure of the paper) and IEF-PCM (water). Modes recomputed from hessian.npy by _vib.modes (PySCF masses,
checked against freq.json); each mode is located from its facts.json wavenumber (*.nu_N4_H4N_inplane_bend,
*.nu_N4_H4N_oop_wag, *.nu_N6_H6N_stretch); of the two listed wag modes the one with the larger H4N weight is drawn.
Common orientation of Fig. 3; largest atomic displacement drawn 1.0 Å.

    python scripts/figures/figS_modes_extra.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt

from figkit import style
import _struct_common as C
import _vib as VB

STEM = 'figS_modes_extra'
WIDTH_MM = 175.0
GAP_MM = 3.0
HEAD_MM = 5.0
ROWLAB_MM = 5.0
NEAR = ['N4', 'N5', 'N6', 'C7', 'C12', 'H4N', 'H6N', 'O1', 'C1', 'C2', 'C6', 'N3', 'O7', 'O6']
CATION = list(range(18, 33))
# (c, f): the triazole side of the cation only, like the crops of (a, b, d, e) (QA round 2, r2-m4: the whole
# cation made rows (c, f) taller and left an empty band under (a, b))
TRIAZOLE_SIDE = ['N4', 'N5', 'N6', 'C7', 'C12', 'H4N', 'H6N']          # bonds to C8, C11 drawn as stubs
HB_ALL = {'gas': [C.HB, (C.I['C11'], C.I['H11'], C.I['O2'])], 'pcm': [C.HB]}
COLS = [('nu_N4_H4N_inplane_bend', 'N4–H4N bend', 'H4N'),         # in-plane bend
        ('nu_N4_H4N_oop_wag', 'N4–H4N wag', 'H4N'),                # out-of-plane wag
        ('nu_N6_H6N_stretch', 'N6–H6N stretch', 'H6N')]            # the free N-H
MIN_LEN = 0.2              # projected arrows shorter than this (Å) are omitted (heads only otherwise);
                           # arrows whose part OUTSIDE the atom ball is shorter than two head lengths as
                           # printed are omitted too (m15; QA round 2, r2-m4: the H6N stub of (a))


def pick(model, key, hlab):
    md = VB.modes(model)
    vals = [float(t) for t in str(C.fact(f'{model}.{key}')).split(',')]
    ks = [VB.find_mode(md, v, tol=0.06) for v in vals]
    k = max(ks, key=lambda kk: VB.h_weight(md, kk, C.I[hlab]))
    return k, float(md['freq_json'][k]), VB.h_weight(md, k, C.I[hlab])


def header_width_mm(letter, title):
    """Width (mm) of a panel header '(x) title' as C.header draws it."""
    f = plt.figure(figsize=(C.mm(WIDTH_MM), C.mm(20.0)))
    arts = C.header(f, 0.0, 10.0, letter, title)
    f.canvas.draw()
    r = f.canvas.get_renderer()
    x1 = max(a.get_window_extent(r).x1 for a in arts) / f.dpi / style.MM
    plt.close(f)
    return x1


def main():
    style.apply()
    g = C.geometries()
    near = [C.I[a] for a in NEAR]
    tri = [C.I[a] for a in TRIAZOLE_SIDE]
    modes = {}
    for m in C.MODELS:
        md = VB.modes(m)
        print(f'{m}: max |recomputed - freq.json| = {md["max_dev"]:.2e} cm-1')
        for key, name, hlab in COLS:
            modes[(m, key)] = pick(m, key, hlab)
            k, nu, w = modes[(m, key)]
            print(f'  {name}: mode {k}, {nu:.1f} cm-1, H weight {w:.2f}')

    def phase(m, key):
        X = g[m]
        R = C.views()[m][0]
        if key == 'nu_N6_H6N_stretch':                       # H6N moving away from N6
            return (C.I['H6N'], X[C.I['H6N']] - X[C.I['N6']])
        if key == 'nu_N4_H4N_oop_wag':                       # H4N moving away from N4 ON THE SCREEN, so the
            b = X[C.I['H4N']] - X[C.I['N4']]                  # arrow does not end on the N4 ball (m15)
            return (C.I['H4N'], R[0] * (R[0] @ b) + R[1] * (R[1] @ b))
        return (C.I['H4N'], R[1])                             # H4N moving up-screen

    def fn(m, key):
        k = modes[(m, key)][0]
        if key == 'nu_N6_H6N_stretch':
            return lambda ax: VB.draw_mode(ax, m, k, atoms=tri, label_atoms=['N5', 'N6', 'H6N'],
                                           hbonds=(), phase=phase(m, key), min_len=MIN_LEN,
                                           min_head_ratio=2.0, head_rule_outside=True)
        return lambda ax: VB.draw_mode(ax, m, k, atoms=near, label_atoms=['N4', 'H4N', 'O1'],
                                       hbonds=HB_ALL[m], phase=phase(m, key), min_len=MIN_LEN,
                                       min_head_ratio=2.0, head_rule_outside=True)

    def box(m, idx):
        X = C.oriented(m)[idx]
        return np.array([X[:, 0].min(), X[:, 0].max(), X[:, 1].min(), X[:, 1].max()])

    W_FIT = WIDTH_MM - 2 * 25.4 * float(plt.rcParams['savefig.pad_inches'])
    rows = {}
    for m in C.MODELS:
        fns = [fn(m, key) for key, _, _ in COLS]
        boxes = [box(m, tri if key == 'nu_N6_H6N_stretch' else near) for key, _, _ in COLS]
        rows[m] = C.fit_row_scale(fns, W_FIT - ROWLAB_MM, GAP_MM, boxes, s0=6.0)
    s = min(rows[m][0] for m in C.MODELS)                    # one scale for all six drawings
    # a column is as wide as its drawing or its header, whichever is wider (the triazole-only (c, f) are narrower
    # than their headers); shrink the scale until the widest row fits
    title = {(m, key): f'{name}, {modes[(m, key)][1]:.0f} cm$^{{-1}}$' for m in C.MODELS for key, name, _ in COLS}
    letter = {(m, key): 'abcdef'[3 * i + j] for i, m in enumerate(C.MODELS) for j, (key, _, _) in enumerate(COLS)}
    hw = {k: header_width_mm(letter[k], title[k]) for k in title}
    avail = W_FIT - ROWLAB_MM
    for _ in range(8):
        boxes = {m: [C.measure_drawing(fn(m, key), b, s) for (key, _, _), b in zip(COLS, rows[m][1])]
                 for m in C.MODELS}
        worst, fixed = 0.0, 0.0
        for m in C.MODELS:
            ws = [(b[1] - b[0]) * s + 1.0 for b in boxes[m]]
            tot = sum(max(w, hw[(m, key)]) for w, (key, _, _) in zip(ws, COLS)) + GAP_MM * (len(COLS) - 1)
            if tot > worst:
                worst = tot
                fixed = sum(hw[(m, key)] for w, (key, _, _) in zip(ws, COLS) if hw[(m, key)] >= w) + \
                    GAP_MM * (len(COLS) - 1)
        if worst <= avail:
            break
        s *= (avail - 0.3 - fixed) / (worst - fixed)
    for m in C.MODELS:
        y0, y1 = min(b[2] for b in boxes[m]), max(b[3] for b in boxes[m])
        rows[m] = [np.array([b[0], b[1], y0, y1]) for b in boxes[m]]
    hrow = {m: (rows[m][0][3] - rows[m][0][2]) * s + 1.0 for m in C.MODELS}
    H = sum(HEAD_MM + hrow[m] for m in C.MODELS) + 2.0
    print(f'scale {s:.2f} mm per Å; height {H:.1f} mm')
    fig = plt.figure(figsize=(C.mm(WIDTH_MM), C.mm(H)))
    y_top = H
    letters = iter('abcdef')
    for m in C.MODELS:
        x = ROWLAB_MM
        for (key, name, hlab), b in zip(COLS, rows[m]):
            ax, w, h = C.drawing_axes(fig, x, y_top - HEAD_MM - hrow[m], b, s)
            fn(m, key)(ax)
            nu = modes[(m, key)][1]
            C.header(fig, x, y_top - HEAD_MM + 1.4, next(letters), f'{name}, {nu:.0f} cm$^{{-1}}$')
            x += max(w, hw[(m, key)]) + GAP_MM
        # row labels in regular weight, as the compose.grid row titles of figS_orbitals_extra (r2-m8)
        fig.text(0.6 / WIDTH_MM, (y_top - HEAD_MM - hrow[m] / 2) / H, style.MODEL_LABELS[m], rotation=90,
                 ha='left', va='center', fontsize=style.BASE_PT)
        y_top -= HEAD_MM + hrow[m] + 2.0
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
