"""Fig. 8 fig_fmo: HOMO and LUMO of the benzotriazolium picrate ion pair, gas phase (a, b) and IEF-PCM (c, d).

Isosurfaces +-0.02 a.u. of the PySCF orbital cubes (homo.cube, lumo.cube, B3LYP/6-311++G(d,p)), drawn with the
camera of the structure figures (one camera per model, picrate in a common position), one print scale for all
four panels, HOMO and LUMO of a model pixel-aligned.  Orbital energies (fmo.json) and Mulliken fragment shares
(fmo.json fragment_composition) are printed in the panel titles and asserted against facts.json.

    python scripts/figures/fig_fmo.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from figkit import style, compose, callouts, draw2d, render3d as r3
import _render_common as C

STEM = 'fig_fmo'
ISO = 0.02
LOBE_OPACITY = 0.82                 # QA round 1, m12: 0.6 let overlapping orange and violet lobes turn brown
# QA round 3, r3-m3: the N4-H4N...O1 dashes in the paper's hydrogen-bond teal (draw2d.HBOND_COLOUR, as in
# Figs. 3 and 6), not the render default dark grey
HB_TEAL = dict(hbond_colour=draw2d.HBOND_COLOUR)
LABELS = ('O1', 'N4')          # the N4-H4N...O1 partners, as orientation cues (only where seen)


def isosurface_inside(cube_path, iso=ISO):
    """The +-iso surfaces must not touch the cube boundary (a clipped lobe would be drawn cut open)."""
    cube = r3.read_cube(cube_path)
    a = np.abs(cube.data)
    edge = max(a[0].max(), a[-1].max(), a[:, 0].max(), a[:, -1].max(), a[:, :, 0].max(), a[:, :, -1].max())
    if edge >= iso:
        raise AssertionError(f'{cube_path}: |psi| = {edge:.4f} on the grid boundary >= isovalue {iso}')
    return edge


def main():
    t0 = time.time()
    style.apply()
    mols = C.molecules()
    V = C.views()
    orbs = ('homo', 'lumo')
    scenes, titles = {}, []
    for m in C.MODELS:
        d = r3.case_dir(C.CASES[m])
        fmo = json.loads((r3.case_dir(C.CASES[m], raw=False) / 'fmo.json').read_text())
        scenes[m] = {}
        for o in orbs:
            edge = isosurface_inside(d / f'{o}.cube')
            scenes[m][o] = r3.orbital_scene(d / f'{o}.cube', mols[m], iso=ISO, opacity=LOBE_OPACITY,
                                           molecule_kw=HB_TEAL)
            e = C.check_equal(fmo[f'{o.upper()}_eV'], f'{m}.{o}_ev', 2)
            comp = fmo['fragment_composition'][o.upper()]
            pic = C.check_equal(comp['picrate anion'], f'{m}.{o}_picrate_pct', 1)
            cat = comp['benzotriazolium cation']
            assert f'{cat:.1f}' == f'{100 - pic:.1f}', (m, o, cat, pic)
            share = f'{pic:.1f}% picrate' if pic >= 50 else f'{cat:.1f}% cation'
            titles.append(f'{style.MODEL_LABELS[m]}, {o.upper()} ({C.num(e, 2)} eV), {share}')
            print(f'{m} {o}: E = {e:.4f} eV, picrate {pic:.1f}%, max |psi| on the grid edge {edge:.4f}')
    # one camera per model, framed on both orbitals; one print scale for all four panels
    vw = {m: V[m].fit(scenes[m]['homo'], scenes[m]['lumo'], mols[m], pad=0.25) for m in C.MODELS}
    pxA = r3.px_per_angstrom(list(vw.values()), C.PANEL2_MM)
    rr = {}
    for m in C.MODELS:
        pair = [r3.render(scenes[m][o], vw[m], px_per_A=pxA, crop=False) for o in orbs]
        rr[m] = r3.crop_together(pair)
        print(m, 'faces: picrate %.2f, cation %.2f' % (C.ring_face(vw[m], mols[m], C.PIC_RING),
                                                        C.ring_face(vw[m], mols[m], C.CAT9)))
    fig, axes = compose.grid([rr['gas'][0], rr['gas'][1], rr['pcm'][0], rr['pcm'][1]], ncols=2,
                             width=style.DOUBLE, titles=titles, gap_mm=(2.5, 2.0))
    # orientation cues: the two hydrogen-bond partners, labelled on the white background (never on a lobe), and
    # only where the atom is seen: an atom under a lobe gets no cue (QA round 2, r2-m5: the N4 leader of (b)
    # ended inside the orange lobe that hides N4)
    for ax, (m, o) in zip(axes, [(m, o) for m in C.MODELS for o in orbs]):
        rd = ax._figkit_render
        cube = r3.read_cube(r3.case_dir(C.CASES[m]) / f'{o}.cube')
        items, hidden = [], []
        for lab in LABELS:
            k = rd.mol.index(lab)
            rad = style.BALL_RADII[rd.mol.symbols[k]]
            if not C.atom_visible(cube, rd.view, rd.mol.coords[k], rad, ISO):
                hidden.append(lab)
                continue
            items.append(dict(pos=rd.mol.coords[k], text=lab, dot=False, radius=rad * rd.px_per_A))
        print(f'{m} {o}: cues {[it["text"] for it in items]}, hidden under a lobe {hidden}')
        assert items, f'{m} {o}: no orientation cue left'
        callouts.place(ax, items, min_leader_mm=0.3, max_leader_mm=7.0, no_leader_below_mm=0.9,
                       name='atom labels')
    C.save(fig, STEM)
    print(f'total {time.time() - t0:.1f} s')


if __name__ == '__main__':
    main()
