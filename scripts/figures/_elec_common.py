"""Shared helpers for the BTPC 2D electronic-structure figures (this share only): orbital ladder + PDOS,
descriptor chart, ELF/LOL maps, MO composition.

    python scripts/figures/<stem>.py      (from the repository root; writes into figures/)

Every number printed on a figure is computed from the data files or read from facts.json; `fact`,
`check_equal` and `fmt` make that explicit. `save` runs figkit's strict QA and writes a 200-dpi print-size
preview into build/.
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np

from figkit import style, geometry as geo, qa

REPO = style.REPO
PAPER = REPO / 'papers' / 'btpc'
FIG = PAPER / 'figures'
BUILD = FIG / 'build'
FACTS_PATH = PAPER / 'data' / 'facts.json'
DERIVED = PAPER / 'data' / 'derived'
CASES = {'gas': 'btpc_gas', 'pcm': 'btpc_pcm'}
MODELS = ('gas', 'pcm')
MINUS = '−'

PIC_COLOUR = '#762A83'          # picrate (PRGn purple end)
CAT_COLOUR = '#1B7837'          # benzotriazolium cation (PRGn green end)
FRAG_COLOURS = {'picrate': PIC_COLOUR, 'cation': CAT_COLOUR}
FRAG_NAMES = {'picrate': 'Picrate', 'cation': 'Benzotriazolium'}

# atom indices (0-based; picrate 0-17 then cation 18-32, CIF order within each ion)
I = dict(O1=0, N3=1, O7=2, C6=3, C1=4, C2=5, C3=6, O6=7, C4=8, C5=9, N1=10, N2=11, O4=12, O5=13, O3=14,
         O2=15, H5=16, H3=17, N4=18, N6=19, C7=20, C12=21, N5=22, C11=23, C10=24, C9=25, C8=26, H11=27,
         H8=28, H10=29, H9=30, H4N=31, H6N=32)
HB = (18, 31, 0)                # N4, H4N, O1

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


def ffact(key):
    return float(fact(key))


def fmt(v, nd, sign=False):
    """Number with a true minus sign (and an explicit + when sign=True)."""
    s = f'{v:+.{nd}f}' if sign else f'{v:.{nd}f}'
    return s.replace('-', MINUS)


def check_equal(value, key, nd, what=''):
    """Assert `value` equals facts.json[key] when both are printed with `nd` decimals; returns the facts value."""
    fv = ffact(key)
    a, b = f'{value:.{nd}f}', f'{fv:.{nd}f}'
    if a != b:
        raise AssertionError(f'{what or key}: computed {a} but facts.json {key} = {b}')
    return fv


def check_close(value, key, tol, what=''):
    fv = ffact(key)
    if abs(value - fv) > tol:
        raise AssertionError(f'{what or key}: computed {value!r} vs facts.json {fv!r} (tol {tol})')
    return fv


# ---------------------------------------------------------------------------------------------------------------
# orbitals
# ---------------------------------------------------------------------------------------------------------------
def mo_table(model):
    """Every MO of a model from the derived Mulliken-composition CSV (validated against fmo.json).

    Returns dict(index0, label, energy (eV), occ, pic (%), cat (%), n_occ, homo (index), lumo (index))."""
    path = DERIVED / f'mo_fragment_composition_{CASES[model]}.csv'
    rows = list(csv.DictReader(path.open()))
    idx = np.array([int(r['index0']) for r in rows])
    assert np.all(idx == np.arange(len(rows))), path
    e = np.array([float(r['energy_eV']) for r in rows])
    occ = np.array([float(r['occupation']) for r in rows])
    pic = np.array([float(r['picrate_pct']) for r in rows])
    cat = np.array([float(r['cation_pct']) for r in rows])
    assert np.all(np.diff(e) >= -1e-9), 'energies not ascending'
    assert np.allclose(pic + cat, 100.0, atol=0.02), 'shares do not add to 100 %'
    n_occ = int((occ > 1.0).sum())
    assert np.all(occ[:n_occ] == 2.0) and np.all(occ[n_occ:] == 0.0)
    labels = [r['label'] for r in rows]
    assert labels[n_occ - 1] == 'HOMO' and labels[n_occ] == 'LUMO'
    # cross-checks against fmo.json and facts.json
    fj = json.loads((REPO / 'molecules' / 'btpc' / 'results' / CASES[model] / 'fmo.json').read_text())
    assert int(fj['n_occ']) == n_occ
    ej = np.sort(np.asarray(fj['mo_energies_eV'], float))
    assert np.allclose(ej[:len(e)], e, atol=1e-4), 'CSV energies differ from fmo.json'
    check_close(e[n_occ - 1], f'{model}.homo_ev', 1e-4)
    check_close(e[n_occ], f'{model}.lumo_ev', 1e-4)
    check_ladder(model, labels, e, pic, n_occ)
    return dict(index0=idx, label=labels, energy=e, occ=occ, pic=pic, cat=cat, n_occ=n_occ,
                homo=n_occ - 1, lumo=n_occ, source=str(path.relative_to(REPO)))


def check_ladder(model, labels, e, pic, n_occ):
    """HOMO-7 ... LUMO+7 energies (3 decimals) and picrate shares (1 decimal) equal facts `*.orbitals_near_gap`."""
    import re
    s = fact(f'{model}.orbitals_near_gap')
    items = re.findall(r'(HOMO(?:-\d+)?|LUMO(?:\+\d+)?) (-?\d+\.\d+) eV \((-?\d+\.\d+)% picrate\)', s)
    assert len(items) == 16, (model, len(items))
    pos = {l: k for k, l in enumerate(labels)}
    for name, ev, pct in items:
        k = pos[name]
        assert abs(e[k] - float(ev)) < 6e-4, (model, name, e[k], ev)
        assert f'{pic[k]:.1f}' == f'{float(pct):.1f}' or abs(pic[k] - float(pct)) < 0.051, (model, name, pic[k], pct)


def orbital_name(k, n_occ):
    """HOMO−n / LUMO+n name (true minus sign) of MO index k."""
    if k < n_occ:
        d = n_occ - 1 - k
        return 'HOMO' if d == 0 else f'HOMO{MINUS}{d}'
    d = k - n_occ
    return 'LUMO' if d == 0 else f'LUMO+{d}'


# ---------------------------------------------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------------------------------------------
_GEOM = None


def geometries():
    """dict(symbols, labels, gas, pcm, case) with (33, 3) arrays in Å and clean CIF labels."""
    global _GEOM
    if _GEOM is None:
        cp = geo.read_case(CASES['pcm'])
        cg = geo.read_case(CASES['gas'])
        pcm = geo.read_xyz(cp['paths'].opt_xyz)
        gas = geo.read_xyz(cg['paths'].opt_xyz)
        labels = list(cp['labels'])
        assert labels == list(cg['labels']) and list(pcm.symbols) == list(gas.symbols)
        assert all('#' not in l for l in labels)
        assert all(labels[i] == n for n, i in I.items())
        _GEOM = dict(symbols=list(pcm.symbols), labels=labels, gas=np.asarray(gas.xyz, float),
                     pcm=np.asarray(pcm.xyz, float), case={'pcm': cp, 'gas': cg})
    return _GEOM


# ---------------------------------------------------------------------------------------------------------------
# saving + print-size preview
# ---------------------------------------------------------------------------------------------------------------
def save(fig, stem, width_mm, preview=True, glyphs=True):
    """qa.save in strict mode (600 dpi PNG + PDF, missing-glyph check on), then a 200-dpi preview at print width.

    FIG_DRAFT=1 in the environment writes to build/draft/ and reports issues without failing (development)."""
    draft = bool(os.environ.get('FIG_DRAFT'))
    target = (BUILD / 'draft' / stem) if draft else (FIG / stem)
    out = qa.save(fig, target, strict=not draft, max_width_mm=width_mm, glyphs=glyphs)
    png = [p for p in out if p.suffix == '.png'][0]
    if preview:
        make_preview(png, stem + ('_draft' if draft else ''), width_mm)
    for p in out:
        print('wrote', p.relative_to(REPO))
    return out


def make_preview(src, stem, width_mm, dpi=200):
    from PIL import Image
    BUILD.mkdir(parents=True, exist_ok=True)
    with Image.open(src) as im:
        d = im.info.get('dpi', (600, 600))[0]
        w_mm = im.width / (d / 25.4)
        h_mm = im.height / (d / 25.4)
        tw = int(round(w_mm / 25.4 * dpi))
        th = int(round(im.height * tw / im.width))
        im.convert('RGB').resize((tw, th), Image.LANCZOS).save(BUILD / f'{stem}_print200.png', dpi=(dpi, dpi))
    print(f'preview build/{stem}_print200.png  ({w_mm:.1f} x {h_mm:.1f} mm; target width {width_mm:.0f} mm)')
    if abs(w_mm - width_mm) > 1.5:
        print(f'  NOTE: printed width {w_mm:.1f} mm differs from the target {width_mm:.0f} mm')
    return w_mm, h_mm


def mm(x):
    return x * style.MM


def panel_title(ax, letter, title='', gap_pt=2.5, y=1.0, pad_pt=1.5, **kw):
    """Bold '(a)' at the upper-left corner of the axes, followed by a normal-weight title at BASE_PT, placed after
    the measured width of the letter (no overlap). Returns (letter_text, title_text or None)."""
    from matplotlib.transforms import offset_copy
    fig = ax.figure
    tr = offset_copy(ax.transAxes, fig=fig, x=0, y=pad_pt, units='points')
    t1 = ax.text(0.0, y, f'({letter})', transform=tr, ha='left', va='baseline', fontsize=style.PANEL_PT,
                 fontweight='bold', **kw)
    if not title:
        return t1, None
    r = fig.canvas.get_renderer()
    w_pt = t1.get_window_extent(r).width * 72 / fig.dpi
    tr2 = offset_copy(ax.transAxes, fig=fig, x=w_pt + gap_pt, y=pad_pt, units='points')
    t2 = ax.text(0.0, y, title, transform=tr2, ha='left', va='baseline', fontsize=style.BASE_PT, **kw)
    return t1, t2
