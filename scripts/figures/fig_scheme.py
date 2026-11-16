"""Fig. 1 fig_scheme: picric acid + 1H-benzotriazole -> benzotriazolium picrate.

RDKit structures (sanitised, formal charges asserted: picric acid 0, benzotriazole 0, picrate -1, cation +1;
CoordGen 2D coordinates, one bond length for all four species) drawn as vector matplotlib by figkit.scheme, so
qa.save checks every symbol. Nitro groups are condensed (NO2) so the only charges drawn are the ionic ones
(O1 minus, N4 plus). The proton that moves is teal: on O1 of picric acid, with a curved arrow to N4 of the base,
and as H4N of the salt, where it donates N4-H4N...O1 (teal dots) as in the crystal. CIF labels (O1, N4, N5, N6,
H4N, H6N) on the salt link the scheme to every later figure. The pKa values are parsed from the abstract-checked
claims in papers/btpc/lit/candidates.jsonl (never typed in); dpKa is computed from them.

    python scripts/figures/fig_scheme.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

from figkit import style, scheme
from figkit.draw2d import HBOND_COLOUR, crop_to_content
from figkit.labels import LabelPlacer
import _struct_common as C

STEM = 'fig_scheme'
WIDTH_MM = 175.0
SYM_PT = 9.0                 # atom symbols
NAME_PT = style.BASE_PT      # compound names
SUB_PT = style.TICK_PT       # pKa lines, systematic name
CIF_PT = style.SMALL_PT      # CIF atom labels
CIF_COLOUR = '0.38'
TEAL = HBOND_COLOUR
LW = 0.8

# SMILES; atom-map numbers tag the CIF atoms (picrate: O1 = 1, C1 = 11, C2 = 12, C6 = 16; cation: N4 = 4, N5 = 5,
# N6 = 6, C7 = 7, C12 = 12)
SMI = {
    'acid': ('[OH:1][c:11]1[c:12]([N+](=O)[O-])[cH][c]([N+](=O)[O-])[cH][c:16]1[N+](=O)[O-]', 0),
    'base': ('c1cc[c:7]2[nH:6][n:5][n:4][c:12]2c1', 0),
    'picrate': ('[O-:1][c:11]1[c:12]([N+](=O)[O-])[cH][c]([N+](=O)[O-])[cH][c:16]1[N+](=O)[O-]', -1),
    'cation': ('c1cc[c:7]2[nH:6][n:5][nH+:4][c:12]2c1', +1),
}


def literature_pka():
    """pKa values from the abstract-checked literature claims (papers/btpc/lit/candidates.jsonl)."""
    rows = {}
    for line in (C.PAPER / 'lit' / 'candidates.jsonl').read_text().splitlines():
        if line.strip():
            d = json.loads(line)
            rows[d['key']] = d
    out = {}
    for key, pat in (('davis1963spectrophotometric', r'pK of picric acid in water at 25 C was determined '
                                                      r'spectrophotometrically as (\d\.\d+)'),
                     ('wang2000photochemistry', r'benzotriazolium cation to pKa = (\d\.\d+)')):
        d = rows[key]
        assert d.get('abstract_checked') is True, key
        m = re.search(pat, ' '.join(d['claims']))
        assert m, (key, d['claims'])
        out[key] = float(m.group(1))
    return out['davis1963spectrophotometric'], out['wang2000photochemistry']


def literature_dpka_range():
    """Range of dpKa over the published values (figures_r2 #9): the flash-photolysis pKa of the cation with the
    spread of the picric acid values gives the lower end, the older value of 1.6 (quoted by Alkorta et al.) the upper."""
    rows = {}
    for line in (C.PAPER / 'lit' / 'candidates.jsonl').read_text().splitlines():
        if line.strip():
            d = json.loads(line)
            rows[d['key']] = d
    claims = lambda k: ' '.join(rows[k]['claims'])
    lo_acid, hi_acid = map(float, re.search(r'spanning roughly (\d\.\d+)-(\d\.\d+)', claims('davis1963spectrophotometric')).groups())
    old = float(re.search(r'Basicity pKa of benzotriazole (\d\.\d+)', claims('alkorta2012theoretical')).group(1))
    new = float(re.search(r'pKa = (\d\.\d+)', claims('wang2000photochemistry')).group(1))
    return new - hi_acid, old - lo_acid


def picrate_like(kind):
    s = scheme.Structure.from_smiles(*SMI[kind])
    O1, C1, C2, C6 = s.idx(1), s.idx(11), s.idx(12), s.idx(16)
    s.orient(C1, O1, 0.0)
    if s.xy[C6, 1] < s.xy[C2, 1]:                 # N3 nitro (on C6) at the top, N1 nitro (on C2) at the bottom,
        s.xy[:, 1] = 2 * s.xy[O1, 1] - s.xy[:, 1]  # as in the 3D figures
    if kind == 'acid':
        s.rotate_atom(s.h_of(O1), O1, 60.0)       # O-H up and to the right (C-O-H 120 deg)
    s.xy -= ring_centre(s, C1)                    # origin at the ring centre
    return s


def ring_centre(s, atom):
    rings = [r for r in s.mol.GetRingInfo().AtomRings() if atom in r]
    return s.xy[list(min(rings, key=len))].mean(0)


def triazole_like(kind):
    s = scheme.Structure.from_smiles(*SMI[kind])
    N4, N5, C12 = s.idx(4), s.idx(5), s.idx(12)
    c5 = ring_centre(s, N4)
    v = s.xy[N4] - c5
    a = np.pi - np.arctan2(v[1], v[0])            # ring centre -> N4 along -x (N4-H points to the picrate O1)
    R = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    s.xy = (s.xy - c5) @ R.T
    if s.xy[N5, 1] < s.xy[C12, 1]:                # N5 at the top, benzo ring below (as in the 3D figures)
        s.xy[:, 1] = -s.xy[:, 1]
    s.xy -= s.xy[N4]                              # origin at N4
    return s


def build(bond, pk_acid, pk_cat, measure=False):
    """Draw the scheme with `bond` mm per bond. Returns (fig, content width mm, info)."""
    acid, base = picrate_like('acid'), triazole_like('base')
    pic, cat = picrate_like('picrate'), triazole_like('cation')
    H_MM = 52.0
    fig = plt.figure(figsize=(C.mm(WIDTH_MM + 40), C.mm(H_MM)))
    ax = scheme.mm_axes(fig)
    y0 = 30.0                                     # the common axis: picrate ring centre, O1, N4 of the base
    gap_plus, gap_arrow, arrow_len = 4.5, 4.0, 17.0
    teal = {'color': TEAL}
    info = {}

    x = 2.0
    # ---- picric acid ----
    xa = x - (acid.extent()[0]) * bond + 3.5      # left O2N text reaches ~3.5 mm past the N2 atom
    hA = acid.h_of(acid.idx(1))
    dA = scheme.draw_structure(ax, acid, origin=(xa, y0), bond=bond, fontsize=SYM_PT, lw=LW,
                               atom_colours={hA: TEAL}, atom_weights={hA: 'bold'})
    info['acid'] = dA
    x_right_acid = max(t.get_window_extent(fig.canvas.get_renderer()).x1 for t in dA.texts) / fig.dpi * 25.4
    # ---- plus sign ----
    xp = x_right_acid + gap_plus
    plus = ax.text(xp, y0, '+', ha='left', va='center_baseline', fontsize=11)
    fig.canvas.draw()
    xp1 = plus.get_window_extent(fig.canvas.get_renderer()).x1 / fig.dpi * 25.4
    # ---- 1H-benzotriazole (origin at N4) ----
    xb = xp1 + gap_plus + 1.6                     # N4 symbol half width
    dB = scheme.draw_structure(ax, base, origin=(xb, y0), bond=bond, fontsize=SYM_PT, lw=LW)
    info['base'] = dB
    xb1 = max(max(p[0] for p in dB.pos), max(t.get_window_extent(fig.canvas.get_renderer()).x1 / fig.dpi * 25.4
                                             for t in dB.texts))
    # ---- reaction arrow ----
    xr0 = xb1 + gap_arrow
    xr1 = xr0 + arrow_len
    ax.add_patch(FancyArrowPatch((xr0, y0), (xr1, y0), arrowstyle='-|>,head_length=4.0,head_width=1.8',
                                 mutation_scale=1.0, lw=LW, color='black', shrinkA=0, shrinkB=0))
    dpka = pk_cat - pk_acid
    info['dpka'] = dpka
    lo, hi = literature_dpka_range()                  # figures_r2 #9: the range, not the lower estimate alone
    assert round(lo, 1) == -0.4 and round(hi, 1) == 1.4, (lo, hi)
    info['dpka_range'] = (lo, hi)
    ax.text((xr0 + xr1) / 2, y0 + 1.6, 'Δp$K_\\mathrm{a}$ ≈', ha='center', va='bottom', fontsize=SUB_PT)
    ax.text((xr0 + xr1) / 2, y0 - 1.6, f'{lo:.1f} to {hi:.1f}'.replace('-', '−'), ha='center', va='top', fontsize=SUB_PT)
    # ---- salt: picrate ... cation ----
    xs = xr1 + gap_arrow - pic.extent()[0] * bond + 3.5
    dP = scheme.draw_structure(ax, pic, origin=(xs, y0), bond=bond, fontsize=SYM_PT, lw=LW)
    O1 = dP.pos[pic.idx(1)]
    hbond_len = 1.55 * bond                       # O1 ... N4 in the drawing (O1...H4N ~ 0.55 bond)
    xc = O1[0] + hbond_len + bond                 # N4 of the cation; H4N sits one bond to its left
    h4 = cat.h_of(cat.idx(4))
    h6 = cat.h_of(cat.idx(6))
    dC = scheme.draw_structure(ax, cat, origin=(xc, y0), bond=bond, fontsize=SYM_PT, lw=LW,
                               atom_colours={h4: TEAL}, atom_weights={h4: 'bold'})
    info['pic'], info['cat'] = dP, dC
    # N4-H4N...O1: teal dots between the symbol boxes
    bo, bh = dP.boxes[pic.idx(1)], dC.boxes[h4]
    p0 = np.array([bo[0] + bo[2] + 0.5, y0])
    # the charge of O1 sits at its upper right; start the dots after the O symbol box
    p1 = np.array([bh[0] - bh[2] - 0.3, y0])
    hb = Line2D([p0[0], p1[0]], [p0[1], p1[1]], color=TEAL, lw=1.1, ls=(0, (0.1, 2.2)), dash_capstyle='round',
                solid_capstyle='round')
    ax.add_line(hb)
    # ---- curved arrow: the acid proton goes to N4 of the base ----
    pH = dA.pos[hA]
    pN = dB.pos[base.idx(4)]
    arc = FancyArrowPatch((pH[0] + 1.0, pH[1] + 1.9), (pN[0] - 0.4, pN[1] + 2.3),
                          connectionstyle='arc3,rad=-0.55', arrowstyle='-|>,head_length=3.2,head_width=1.5',
                          mutation_scale=1.0, lw=0.9, color=TEAL, shrinkA=0, shrinkB=0)
    ax.add_patch(arc)
    # apex of the arc for its label
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    ab = arc.get_window_extent(r)
    xm = (ab.x0 + ab.x1) / 2 / fig.dpi * 25.4
    ym = ab.y1 / fig.dpi * 25.4
    ax.text(xm, ym + 0.8, 'H$^+$', ha='center', va='bottom', fontsize=SUB_PT, color=TEAL)
    # ---- CIF labels on the salt ----
    pl = LabelPlacer(ax, name='scheme labels', n_rings=4, clearance_pt=1.2, bounds=None)
    for d in (dA, dB, dP, dC):
        scheme.add_obstacles(pl, d)
    pl.add_segment(p0, p1, hw_pt=1.0)
    lab = {pic.idx(1): 'O1', cat.idx(4): 'N4', cat.idx(5): 'N5', cat.idx(6): 'N6', h4: 'H4N', h6: 'H6N'}
    # QA round 1, m1: H4N up and to the left, clear of the N4 charge; H6N right next to its H (smaller exclusion
    # radius: an H symbol is narrower than an N); O1 labelled on the acid too, where the caption names it
    prefer = {'O1': (0.2, -1), 'N4': (0.2, -1), 'N5': (0, 1), 'N6': (1, 0.2), 'H4N': (-0.8, 1), 'H6N': (1, 0.3)}
    radius = {'H4N': 4.6, 'H6N': 3.3}
    for i, t in lab.items():
        d = dP if t == 'O1' else dC
        pl.add_item(t, d.pos[i], radius_pt=radius.get(t, 4.6), prefer=prefer[t], priority=1, key=t, max_ring=3,
                    leader=False, kw=dict(fontsize=CIF_PT, color=CIF_COLOUR))
    pl.add_item('O1', dA.pos[acid.idx(1)], radius_pt=4.6, prefer=(0.2, -1), priority=1, key='O1 acid',
                max_ring=3, leader=False, kw=dict(fontsize=CIF_PT, color=CIF_COLOUR))
    placed, dropped = pl.place()
    assert not dropped, [dd.text for dd in dropped]
    # ---- names and pKa lines under each block ----
    fig.canvas.draw()
    r = fig.canvas.get_renderer()

    def block_box(ds, extra=()):
        arts = [t for d in ds for t in d.texts] + [ln for d in ds for ln in d.lines] + list(extra)
        bbs = [a.get_window_extent(r) for a in arts]
        return (min(b.x0 for b in bbs) / fig.dpi * 25.4, max(b.x1 for b in bbs) / fig.dpi * 25.4,
                min(b.y0 for b in bbs) / fig.dpi * 25.4, max(b.y1 for b in bbs) / fig.dpi * 25.4)

    boxes = {'acid': block_box([dA]), 'base': block_box([dB]),
             'salt': block_box([dP, dC], [it.artist for it in placed if it.key != 'O1 acid'])}
    y_name = min(b[2] for b in boxes.values()) - 3.4
    names = {'acid': ['picric acid', f'p$K_\\mathrm{{a}}$ = {pk_acid:.2f}'],
             'base': ['1$\\mathit{H}$-benzotriazole', f'p$K_\\mathrm{{a}}$ = {pk_cat:.2f} (cation)'],
             'salt': ['benzotriazolium picrate',
                      '(1$\\mathit{H}$-benzotriazol-3-ium 2,4,6-trinitrophenolate)']}
    for k, lines in names.items():
        xc_ = (boxes[k][0] + boxes[k][1]) / 2
        if k == 'salt':                              # centre the name under the whole ion pair
            xc_ = (boxes[k][0] + boxes[k][1]) / 2
        ax.text(xc_, y_name, lines[0], ha='center', va='baseline', fontsize=NAME_PT)
        ax.text(xc_, y_name - 3.6, lines[1], ha='center', va='baseline', fontsize=SUB_PT)
    info['boxes'] = boxes
    # shrink the figure to what is drawn; 1 data unit stays 1 mm
    W, H = crop_to_content(fig, ax, pad_mm=0.3)
    return fig, W * 25.4, info


def main():
    style.apply()
    pk_acid, pk_cat = literature_pka()
    print(f'pKa picric acid {pk_acid}, benzotriazolium {pk_cat}; dpKa {pk_cat - pk_acid:.2f}')
    # formal charges of the four species are asserted in Structure.from_smiles; check the CIF-based protonation:
    cat = scheme.Structure.from_smiles(*SMI['cation'])
    n4, n5, n6 = cat.idx(4), cat.idx(5), cat.idx(6)
    a = cat.mol.GetAtomWithIdx
    assert a(n4).GetFormalCharge() == 1 and a(n4).GetTotalNumHs(True) == 1 and a(n6).GetTotalNumHs(True) == 1
    assert a(n5).GetTotalNumHs(True) == 0                    # N5 carries no hydrogen (CIF)
    target = WIDTH_MM - 2 * 25.4 * float(plt.rcParams['savefig.pad_inches'])
    bond, pts = 5.0, []
    w = None
    for _ in range(6):
        fig, w, info = build(bond, pk_acid, pk_cat)
        pts.append((bond, w))
        if target - 0.6 <= w <= target:
            break
        plt.close(fig)
        if len(pts) >= 2:
            (b1, w1), (b2, w2) = pts[-2], pts[-1]
            bond = b2 + (target - 0.3 - w2) * (b2 - b1) / (w2 - w1)
        else:
            bond = bond * (target - 0.3) / w * 1.02
    print(f'bond {bond:.3f} mm, content width {w:.1f} mm')
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
