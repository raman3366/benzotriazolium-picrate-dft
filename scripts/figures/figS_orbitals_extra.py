"""SI figure figS_orbitals_extra: the benzotriazolium pi orbital, LUMO+1 and the pi* orbital of both models.

Isosurfaces +-0.02 a.u. of molecular orbitals evaluated from final.molden (B3LYP/6-311++G(d,p)) with
figkit.molden on a box with a 7 Bohr margin (the production orbital cubes, 4.5 Bohr margin, exist only for the
HOMO and LUMO).  This is post-processing of the converged wavefunction, not a new calculation.  Gate (PLAN.md):
the evaluator must reproduce the production homo.cube and lumo.cube to max |dpsi| < 1e-4 a.u.; the full-grid check
is build/r3d/validate_molden.log and this script repeats it on 20 000 random grid points of both cubes before
drawing.  Columns: the highest occupied orbital located on the cation (gas phase HOMO-5, IEF-PCM HOMO-1; round 2,
figures_r2 #4: pairing HOMO-1 by label set a picrate orbital over a cation one), LUMO+1, and the lowest virtual
orbital located on the cation (gas phase: the LUMO;
IEF-PCM: LUMO+3); rows: gas phase, IEF-PCM (water).  Camera of the structure figures (one per model), one print
scale for all six panels.  Energies and Mulliken fragment shares come from
papers/btpc/data/derived/mo_fragment_composition_btpc_<model>.csv and are asserted against facts.json
(*.homo_m1_ev, *.lumo_p1_ev, *.cation_lowest_virt_ev, *.orbitals_near_gap) and against the molden energies.
Cubes are cached in build/cubes/ (derived, regenerable).

    python scripts/figures/figS_orbitals_extra.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import csv
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from figkit import style, compose, callouts, draw2d, molden, render3d as r3
import _render_common as C

STEM = 'figS_orbitals_extra'
ISO = 0.02
LOBE_OPACITY = 0.82                 # QA round 1, m12: 0.6 let overlapping orange and violet lobes turn brown
# QA round 3, r3-m3: the N4-H4N...O1 dashes in the paper's hydrogen-bond teal (draw2d.HBOND_COLOUR, as in
# Figs. 3 and 6), not the render default dark grey
HB_TEAL = dict(hbond_colour=draw2d.HBOND_COLOUR)
MARGIN_BOHR, SPACING_BOHR = 7.0, 0.3
CUBES = C.BUILD / 'cubes'
HARTREE_EV = 27.211386245988
LABELS = ('O1', 'N4')
PREFER = {'O1': (-0.3, -1), 'N4': (0.3, 1)}      # O1 label below, N4 label above: never side by side
COLUMNS = ('cation_occ', 'LUMO+1', 'cation')    # figures_r2 #4: first column pairs the cation pi orbital by character


def composition(model):
    """Rows of the derived Mulliken table keyed by label ('HOMO-1', 'LUMO', ...), with index0 = molden MO."""
    path = C.PAPER / 'data' / 'derived' / f'mo_fragment_composition_{C.CASES[model]}.csv'
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {r['label']: dict(index=int(r['index0']), e=float(r['energy_eV']), occ=float(r['occupation']),
                             pic=float(r['picrate_pct']), cat=float(r['cation_pct'])) for r in rows}


def near_gap(model):
    """facts <model>.orbitals_near_gap parsed into {label: (energy eV, picrate %)}."""
    txt = C.fact(f'{model}.orbitals_near_gap')
    return {m.group(1): (float(m.group(2)), float(m.group(3)))
            for m in re.finditer(r'(HOMO(?:-\d+)?|LUMO(?:\+\d+)?) (-?\d+\.\d+) eV \((-?\d+\.\d+)% picrate\)', txt)}


def cation_occupied(comp):
    """Label of the highest occupied orbital with a cation share above 50 % (the benzotriazolium pi orbital:
    HOMO-5 in the gas phase, HOMO-1 in IEF-PCM)."""
    occ = sorted((r['e'], lab) for lab, r in comp.items() if r['occ'] > 1.5 and r['cat'] > 50)
    return occ[-1][1]


def cation_virtual(comp):
    """Label of the lowest virtual orbital with a cation share above 50 %."""
    virt = sorted((r['e'], lab) for lab, r in comp.items() if r['occ'] < 0.5 and r['cat'] > 50)
    return virt[0][1]


def check_gate(mf, d):
    """Repeat the validation gate on 20 000 random points of the production HOMO and LUMO cubes."""
    rng = np.random.default_rng(0)
    worst = 0.0
    for name, mo in (('homo', mf.homo), ('lumo', mf.lumo)):
        cube = r3.read_cube(d / f'{name}.cube')
        idx = np.stack([rng.integers(0, n, 20000) for n in cube.shape], axis=1)
        P = (cube.origin + idx @ cube.axes) / molden.BOHR                     # Angstrom -> Bohr
        ref = cube.data[idx[:, 0], idx[:, 1], idx[:, 2]]
        v = mf.mo_on_points([mo], P)[0]
        s = np.sign((v * ref).sum())
        worst = max(worst, float(np.abs(s * v - ref).max()))
    if worst >= 1e-4:
        raise AssertionError(f'molden evaluator misses the production cubes by {worst:.2e} a.u. (gate 1e-4)')
    return worst


def mo_cube(model, mf, mo):
    """Cube of MO `mo` (0-based molden index) on the wide box; cached in build/cubes/."""
    CUBES.mkdir(parents=True, exist_ok=True)
    path = CUBES / f'{C.CASES[model]}_mo{mo:03d}_m{MARGIN_BOHR:g}_s{SPACING_BOHR:g}.cube'
    src = r3.case_dir(C.CASES[model]) / 'final.molden'
    if not path.exists() or path.stat().st_mtime < src.stat().st_mtime:
        grid = molden.box_grid(mf.coords, margin=MARGIN_BOHR, spacing=SPACING_BOHR)
        vals = mf.mo_on_grid([mo], grid)[0]
        molden.write_cube(path, mf, grid, vals, title=f'{C.CASES[model]} MO {mo}',
                          comment=f'figkit.molden from final.molden, margin {MARGIN_BOHR} Bohr')
        print('wrote', path.relative_to(C.REPO), grid.shape)
    return path


def boundary_max(cube):
    a = np.abs(cube.data)
    return max(a[0].max(), a[-1].max(), a[:, 0].max(), a[:, -1].max(), a[:, :, 0].max(), a[:, :, -1].max())


def main():
    t0 = time.time()
    style.apply()
    mols = C.molecules()
    V = C.views()
    scenes, titles, cubes = {}, {}, {}
    for m in C.MODELS:
        d = r3.case_dir(C.CASES[m])
        mf = molden.read_molden(d / 'final.molden')
        gate = check_gate(mf, d)
        comp, ng = composition(m), near_gap(m)
        assert comp['HOMO']['index'] == mf.homo, (m, comp['HOMO'], mf.homo)
        cat = cation_virtual(comp)
        cat_occ = cation_occupied(comp)
        labs = {'cation_occ': cat_occ, 'LUMO+1': 'LUMO+1', 'cation': cat}
        C.check_equal(comp[cat_occ]['e'], f'{m}.cation_highest_occ_ev', 3)
        if m == 'pcm':
            assert cat_occ == 'HOMO-1', cat_occ
            C.check_equal(comp['HOMO-1']['e'], f'{m}.homo_m1_ev', 3)
        C.check_equal(comp['LUMO+1']['e'], f'{m}.lumo_p1_ev', 3)
        C.check_equal(comp[cat]['e'], f'{m}.cation_lowest_virt_ev', 3)
        scenes[m] = {}
        for col in COLUMNS:
            lab = labs[col]
            r = comp[lab]
            e_molden = mf.mo_energy[r['index']] * HARTREE_EV
            assert abs(e_molden - r['e']) < 5e-4, (m, lab, e_molden, r['e'])
            fe, fp = ng[lab]
            assert f"{r['e']:.3f}" == f'{fe:.3f}' and f"{r['pic']:.1f}" == f'{fp:.1f}', (m, lab, r, ng[lab])
            cube = r3.read_cube(mo_cube(m, mf, r['index']))
            edge = boundary_max(cube)
            if edge >= ISO:
                raise AssertionError(f'{m} {lab}: |psi| = {edge:.4f} on the box boundary >= {ISO}')
            scenes[m][col] = r3.orbital_scene(cube, mols[m], iso=ISO, opacity=LOBE_OPACITY,
                                           molecule_kw=HB_TEAL)
            cubes[(m, col)] = cube
            share = f"{r['pic']:.1f}% picrate" if r['pic'] >= 50 else f"{r['cat']:.1f}% cation"
            titles[(m, col)] = f"{lab.replace('-', C.MINUS)} ({C.num(r['e'], 2)} eV), {share}"
            print(f'{m} {col:7s} {lab:7s} MO {r["index"]}: E = {r["e"]:.3f} eV, picrate {r["pic"]:.1f}%, '
                  f'max |psi| on the box edge {edge:.4f}')
        print(f'{m}: gate max |dpsi| on 2 x 20000 cube points {gate:.1e} a.u.')
    vw = {m: V[m].fit(*scenes[m].values(), mols[m], pad=0.25) for m in C.MODELS}
    pxA = r3.px_per_angstrom(list(vw.values()), (175.0 - 2 * C.MARGIN_MM - 4.2 - 2 * 2.0) / 3)
    rr = {}
    for m in C.MODELS:
        row = [r3.render(scenes[m][col], vw[m], px_per_A=pxA, crop=False) for col in COLUMNS]
        rr[m] = r3.crop_together(row)
    order = [(m, col) for m in C.MODELS for col in COLUMNS]
    fig, axes = compose.grid([rr[m][k] for m in C.MODELS for k in range(len(COLUMNS))], ncols=3,
                             width=style.DOUBLE, titles=[titles[k] for k in order], gap_mm=(2.0, 2.0),
                             row_titles=[style.MODEL_LABELS[m] for m in C.MODELS])
    # orientation cues only on atoms that are seen (QA round 2, r2-m5: in (a) the O1 and N4 leaders both ended
    # inside lobes); an atom under a lobe gets no cue
    for ax, (m, col) in zip(axes, order):
        rd = ax._figkit_render
        items, hidden = [], []
        for lab in LABELS:
            k = rd.mol.index(lab)
            rad = style.BALL_RADII[rd.mol.symbols[k]]
            if not C.atom_visible(cubes[(m, col)], rd.view, rd.mol.coords[k], rad, ISO):
                hidden.append(lab)
                continue
            items.append(dict(pos=rd.mol.coords[k], text=lab, dot=False, prefer=PREFER[lab],
                              radius=rad * rd.px_per_A))
        print(f'{m} {col}: cues {[it["text"] for it in items]}, hidden under a lobe {hidden}')
        if items:
            callouts.place(ax, items, min_leader_mm=0.3, max_leader_mm=7.0, no_leader_below_mm=0.9,
                           name='atom labels')
    C.save(fig, STEM)
    print(f'total {time.time() - t0:.1f} s')


if __name__ == '__main__':
    main()
