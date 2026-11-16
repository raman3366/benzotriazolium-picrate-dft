"""Shared helpers for the BTPC 3D render figures: fig_fmo (Fig. 8), fig_nci (Fig. 7), fig_mep (Fig. 11),
figS_mep_extrema and figS_orbitals_extra.

    python scripts/figures/<stem>.py      (from the repository root; writes into figures/)

Camera: ONE orientation per model, the one used by the structure drawings of the paper (Fig. 3,
scripts/_struct_common.py): face view of the picrate heavy atoms tilted by 35 deg towards the screen direction
210 deg, rolled so that picrate is on the left and the cation on the right; the gas-phase camera is the IEF-PCM one
carried over by a Kabsch fit on the picrate heavy atoms + N4 + H4N.  In this view both rings are seen at 34-44 deg
from face-on in both models and the N4-H4N...O1 NCI disc is 51 % (IEF-PCM) and 69 % (gas) face-on, with about
two thirds of its vertices unobstructed (checked with the z-buffer, see check_view()).  When the structure
scripts have written build/views_struct.json, that file is used (so the renders follow Fig. 3 exactly) and
compared with the recipe here.

Everything printed on a figure is computed from the data files and checked against facts.json (`fact`,
`check_equal`).
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

from figkit import style, compose, render3d as r3

REPO = style.REPO
PAPER = REPO / 'papers' / 'btpc'
FIG = PAPER / 'figures'
BUILD = FIG / 'build'
FACTS_PATH = PAPER / 'data' / 'facts.json'
CASES = {'gas': 'btpc_gas', 'pcm': 'btpc_pcm'}
MODELS = ('gas', 'pcm')
MINUS = '−'

PIC = list(range(18))
CAT = list(range(18, 33))
PIC_HEAVY = list(range(16))
PIC_RING = [4, 5, 6, 8, 9, 3]            # C1 C2 C3 C4 C5 C6
CAT9 = list(range(18, 27))               # N4 N6 C7 C12 N5 C11 C10 C9 C8
HB = (18, 31, 0)                         # N4, H4N, O1
I = dict(O1=0, N3=1, O7=2, C6=3, C1=4, C2=5, C3=6, O6=7, C4=8, C5=9, N1=10, N2=11, O4=12, O5=13, O3=14,
         O2=15, H5=16, H3=17, N4=18, N6=19, C7=20, C12=21, N5=22, C11=23, C10=24, C9=25, C8=26, H11=27,
         H8=28, H10=29, H9=30, H4N=31, H6N=32)

TILT, PHI = 35.0, 210.0                  # the camera recipe of _struct_common.py
MARGIN_MM = compose.OUTER_MM
PANEL2_MM = (175.0 - 2 * MARGIN_MM - 2.5) / 2      # two render panels across a double column


# ---------------------------------------------------------------------------------------------------------------
# facts
# ---------------------------------------------------------------------------------------------------------------
_FACTS = None


def facts():
    global _FACTS
    if _FACTS is None:
        _FACTS = json.loads(FACTS_PATH.read_text())
    return _FACTS


def fact(key):
    return facts()[key]['value']


def num(v, nd, sign=False):
    """Number with a true minus sign (and an explicit + when sign=True)."""
    s = f'{v:+.{nd}f}' if sign else f'{v:.{nd}f}'
    return s.replace('-', MINUS)


def check_equal(value, key, nd, what=''):
    """Assert `value` printed with nd decimals equals facts.json[key] at the same precision; returns the fact."""
    fv = float(fact(key))
    a, b = f'{value:.{nd}f}', f'{fv:.{nd}f}'
    if a != b:
        raise AssertionError(f'{what or key}: computed {a} but facts.json {key} = {b}')
    return fv


# ---------------------------------------------------------------------------------------------------------------
# molecules and the common camera
# ---------------------------------------------------------------------------------------------------------------
_MOLS = None
_VIEWS = None


def molecules():
    global _MOLS
    if _MOLS is None:
        _MOLS = {m: r3.load_molecule(CASES[m]) for m in MODELS}
        for m, mol in _MOLS.items():
            assert all('#' not in l for l in mol.labels), mol.labels
            assert [mol.labels[i] for i in HB] == ['N4', 'H4N', 'O1']
            assert all(mol.labels[i] == k for k, i in I.items())
    return _MOLS


def _recipe_views():
    mols = molecules()
    mp = mols['pcm']
    V = r3.auto_view(mp, mode='face', atoms=PIC_HEAVY)
    V = V.roll(PHI).elevation(TILT).roll(-PHI)
    for _ in range(2):
        s = V.screen(np.vstack([mp.coords[PIC_RING].mean(0), mp.coords[CAT9].mean(0)]))
        d = s[1, :2] - s[0, :2]
        theta = float(np.degrees(np.arctan2(d[1], d[0])))
        if abs(theta) < 1e-3:
            break
        V = V.roll(-theta)
    Vg = r3.matched_view(V, mp, mols['gas'], atoms=PIC_HEAVY + [18, 31])
    return {'pcm': V, 'gas': Vg}


def _angle(a, b):
    return math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(a, b))))))


def views():
    """{'gas': View, 'pcm': View} (unframed: call .fit(scene, ...) per figure)."""
    global _VIEWS
    if _VIEWS is not None:
        return _VIEWS
    rec = _recipe_views()
    out = dict(rec)
    js = BUILD / 'views_struct.json'
    if js.exists():
        d = json.loads(js.read_text())
        dev = max(max(_angle(rec[m].direction, d[m]['direction']), _angle(rec[m].up, d[m]['up'])) for m in MODELS)
        if dev > 0.05:
            print(f'NOTE: build/views_struct.json differs from the recipe here by {dev:.2f} deg; using the file so '
                  f'the renders keep the orientation of the structure figures')
            for m in MODELS:
                out[m] = r3.View(focal=np.array(d[m]['focal']), direction=np.array(d[m]['direction']),
                                 up=np.array(d[m]['up']))
    _VIEWS = out
    (BUILD / 'views_render3d.json').parent.mkdir(parents=True, exist_ok=True)
    (BUILD / 'views_render3d.json').write_text(json.dumps({m: v.as_dict() for m, v in out.items()}, indent=1))
    return out


def ring_face(view, mol, idx):
    """|cos| of the angle between a ring normal and the viewing direction (1 = face-on)."""
    P = mol.coords[idx]
    _, _, vt = np.linalg.svd(P - P.mean(0))
    return abs(float(vt[2] @ view.direction))


def check_view(view, mol, scene, render, pairs, min_visible=0.4, radius=0.8):
    """z-buffer check that the NCI patch between each atom pair (a, b) is largely unobstructed in the render.
    Returns {pair: visible fraction}; raises if a patch is hidden."""
    surf = next(ms for ms, k in scene.items if 'sl2r' in ms.point_data)
    P = np.asarray(surf.points)
    out = {}
    for a, b in pairs:
        c = 0.5 * (mol.coords[a] + mol.coords[b])
        pts = P[np.linalg.norm(P - c, axis=1) < radius]
        if not len(pts):
            continue
        f = float(render.visible(pts, towards_camera=0.02).mean())
        out[(a, b)] = f
        if f < min_visible:
            raise AssertionError(f'NCI patch {mol.labels[a]}...{mol.labels[b]} only {f:.0%} visible in this view')
    return out


# ---------------------------------------------------------------------------------------------------------------
# saving
# ---------------------------------------------------------------------------------------------------------------
def save(fig, stem):
    """compose.save in strict mode (600 dpi PNG + PDF, figkit.qa checks), then a 200 dpi print-size preview."""
    out = compose.save(fig, FIG / stem, strict=True)
    width_mm = fig._figkit_grid['width'] / style.MM
    preview(stem, width_mm)
    for p in out:
        print('wrote', p.relative_to(REPO))
    return out


def preview(stem, width_mm, dpi=200):
    from PIL import Image
    src = FIG / f'{stem}.png'
    BUILD.mkdir(parents=True, exist_ok=True)
    with Image.open(src) as im:
        w_mm = im.width / (im.info.get('dpi', (600, 600))[0] / 25.4)
        tw = int(round(w_mm / 25.4 * dpi))
        th = int(round(im.height * tw / im.width))
        im.convert('RGB').resize((tw, th), Image.LANCZOS).save(BUILD / f'{stem}_print200.png', dpi=(dpi, dpi))
        h_mm = im.height / (im.info.get('dpi', (600, 600))[1] / 25.4)
    print(f'preview build/{stem}_print200.png ({w_mm:.1f} x {h_mm:.1f} mm; target width {width_mm:.0f} mm)')


def mm_to_px(mm, px_per_A=None, dpi=style.DPI):
    """Image pixels for a length in mm at print size, assuming the render is printed at `dpi` (px_per_angstrom
    targets 600 dpi at the panel width)."""
    return mm / 25.4 * dpi


# ---------------------------------------------------------------------------------------------------------------
# orbital figures: is an atom seen, or hidden under a (translucent) lobe?  (QA round 2, r2-m5)
# ---------------------------------------------------------------------------------------------------------------
def _lobe_on_ray(cube, P0, direction, iso, length=15.0, step=0.04):
    """True where the ray from P0 towards the camera meets |psi| >= iso."""
    from scipy.ndimage import map_coordinates
    t = np.arange(0.0, length, step)
    P = P0[None] + t[:, None] * direction[None]
    idx = np.linalg.solve(np.asarray(cube.axes, float).T, (P - cube.origin).T)
    v = map_coordinates(cube.data, idx, order=1, mode='constant', cval=0.0)
    return bool((np.abs(v) >= iso).any())


def atom_visible(cube, view, pos, radius, iso, n_rim=4, rim=0.7):
    """True when an atom ball (centre pos, radius Å) is seen, not covered by an isosurface |psi| = iso of `cube`:
    the rays towards the camera from the front of the ball at its centre and at n_rim points `rim` x radius off
    the centre (in the screen plane) are tested, and at most one of them may meet a lobe.  The render's z-buffer
    cannot answer this, because translucent lobes do not write it."""
    d = np.asarray(view.direction, float)
    right = np.asarray(view.right, float)
    up = np.cross(d, right)
    starts = [np.asarray(pos, float) + radius * d]
    for k in range(n_rim):
        a = 2 * np.pi * k / n_rim
        starts.append(starts[0] + rim * radius * (np.cos(a) * right + np.sin(a) * up))
    hits = [_lobe_on_ray(cube, s, d, iso) for s in starts]
    return (not hits[0]) and sum(hits) <= 1
