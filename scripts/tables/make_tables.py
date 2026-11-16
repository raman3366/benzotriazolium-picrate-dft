#!/usr/bin/env python3
"""make_tables.py: every table of the benzotriazolium picrate (BTPC) paper, main text and SI.

Run (from the repository root; docs/REPRODUCE.md):
    python scripts/tables/make_tables.py [--check]

Paths below are workspace paths (papers/btpc/tables = tables/ here); scripts/paths.py maps them to this dataset.

Writes papers/btpc/tables/<stem>.md (one pandoc pipe table per file; caption line
'Table: <caption> <notes> {#tbl:<stem>}' after the table; notes and footnotes are sentences of the caption) and
papers/btpc/tables/INDEX.md (placement, label, facts keys and other inputs of every table).
--check additionally runs pandoc (environment variable PANDOC, else pandoc on the PATH) over the outputs (HTML for every table, DOCX and HTML for one table,
and a synthetic main/SI document through the manuscript filters include.lua + crossref.lua); test files go to
papers/btpc/tables/build/.

Numbers come from papers/btpc/data/facts.json (value, unit, source, note of each key). The SI coordinate,
frequency and all-bond tables also read the structure files (opt.xyz, the X-ray cluster xyz), freq.json and
freq_step0.005.json, and several tables cross-check facts against the files they were derived from
(nci_spikes.csv, mep_extrema.csv, crystal_hbonds.csv, mo_fragment_composition_*.csv).

The script stops with a message naming the table and the key on: a missing facts key, a unit other than the
expected one, a note or value that does not match the pattern the table relies on, a failed cross-check, or
markup that breaks the house style (hyphen used as minus, '#' symmetry suffixes, 'Schiff', ASCII ranges...).
"""
from __future__ import annotations

import argparse
import ast
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # scripts/ (paths.py)
import paths                                       # noqa: E402
REPO = paths.legacy_root()                         # the dataset files under their workspace paths (scripts/paths.py)
PAPER = REPO / 'papers' / 'btpc'
HERE = PAPER / 'tables'                            # = tables/ of the dataset
sys.path.insert(0, str(REPO / 'papers' / '_shared'))
from figkit import geometry as geo                 # noqa: E402  (clean_labels, read_xyz)

FACTS_PATH = PAPER / 'data' / 'facts.json'
BUILD_FACTS = PAPER / 'data' / 'build_facts.py'
DERIVED = PAPER / 'data' / 'derived'
MOL = REPO / 'molecules' / 'btpc'
RES = MOL / 'results'
RAW = REPO / 'data' / 'r2' / 'molecules' / 'btpc' / 'results'   # raw results archive (Multiwfn surfanalysis.txt)
PANDOC = Path(os.environ.get('PANDOC') or shutil.which('pandoc') or REPO / 'tools' / 'bin' / 'pandoc')
SHARED_TOOLS = REPO / 'papers' / '_shared' / 'tools'

MINUS = '−'
NDASH = '–'
NBSP = '\u00a0'
MODEL = {'xray': 'X-ray', 'pcm': 'IEF-PCM (water)', 'gas': 'Gas phase'}
CASE = {'pcm': 'btpc_pcm', 'gas': 'btpc_gas'}
WIDTH_CAP = 60                                     # max dashes per column (relative widths in pandoc)


# =================================================================================================
# fail-loud access to facts.json
# =================================================================================================
class TableError(RuntimeError):
    pass


FACTS = json.loads(FACTS_PATH.read_text(encoding='utf-8'))


class _Ctx:
    stem = None
    keys: list = []
    files: list = []
    checks: list = []


CTX = _Ctx()


def fail(msg):
    raise TableError(f'[{CTX.stem or "make_tables"}] {msg}')


def rec(key):
    if key not in FACTS:
        fail(f"facts.json has no key '{key}'")
    if key not in CTX.keys:
        CTX.keys.append(key)
    return FACTS[key]


def val(key, unit=None):
    r = rec(key)
    if unit is not None and r['unit'] != unit:
        fail(f"key '{key}' has unit {r['unit']!r}, expected {unit!r}")
    return r['value']


def fnum(key, unit=None):
    v = val(key, unit)
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        fail(f"key '{key}' is not a number: {v!r}")
    return float(v)


def fint(key, unit=''):
    v = val(key, unit)
    if isinstance(v, bool) or not isinstance(v, int):
        fail(f"key '{key}' is not an integer: {v!r}")
    return v


def fstr(key, unit=''):
    v = val(key, unit)
    if not isinstance(v, str):
        fail(f"key '{key}' is not a string: {v!r}")
    return v


def note(key):
    return rec(key)['note']


def _grab(text, pattern, key, what, flags=0):
    m = re.search(pattern, text, flags)
    if not m:
        fail(f"{what} of '{key}' does not match /{pattern}/: {text!r}")
    return m


def from_note(key, pattern, flags=0):
    return _grab(note(key), pattern, key, 'note', flags)


def from_value(key, pattern, flags=0):
    return _grab(fstr(key, rec(key)['unit']), pattern, key, 'value', flags)


def numbered_keys(prefix, pattern):
    """Existing keys '<prefix><n>' (n = pattern group), sorted by n, recorded, checked to run 1..N."""
    ks = [k for k in FACTS if re.fullmatch(re.escape(prefix) + pattern, k)]
    if not ks:
        fail(f"no facts keys match '{prefix}{pattern}'")
    ks.sort(key=lambda k: int(k[len(prefix):]))
    nums = [int(k[len(prefix):]) for k in ks]
    if nums != list(range(1, len(nums) + 1)):
        fail(f"keys '{prefix}NN' are not numbered 1..{len(nums)}: {nums}")
    for k in ks:
        rec(k)
    return ks


def use_file(path):
    p = Path(os.path.normpath(Path(path).absolute()))   # keep the links of the .layout view (no resolve())
    rel = str(p.relative_to(REPO))
    if rel not in CTX.files:
        CTX.files.append(rel)
    if not p.exists():
        fail(f'missing input file {rel}')
    return p


def check(ok, what):
    """Record a cross-check (listed in INDEX.md) and stop if it fails."""
    if not ok:
        fail(f'cross-check failed: {what}')
    if what not in CTX.checks:
        CTX.checks.append(what)


def close(a, b, tol, what):
    check(abs(a - b) <= tol, f'{what} (|{a!r} - {b!r}| <= {tol:g})')


# =================================================================================================
# number and text formatting
# =================================================================================================
def fmt(x, nd, sign=False):
    """Fixed decimals with a true minus sign; '+' only when sign=True; no signed zero."""
    s = f'{x:+.{nd}f}' if sign else f'{x:.{nd}f}'
    if float(s) == 0.0:
        s = s.lstrip('+-')
    return s.replace('-', MINUS)


_SUP = str.maketrans('0123456789-+', '⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺')


def sci(x, nd=1):
    """1.403e-4 -> '1.4 × 10⁻⁴' (true minus in the mantissa, Unicode superscript exponent)."""
    if x == 0:
        return '0'
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if round(abs(m), nd) >= 10:
        m, e = m / 10, e + 1
    return f'{fmt(m, nd)} × 10{str(e).translate(_SUP)}'


def sci_str(s, nd=None):
    """'1.0e-6' -> '1.0 × 10⁻⁶', keeping the decimals written in the source string."""
    m = re.fullmatch(r'(-?\d+(?:\.(\d+))?)[eE]([-+]?\d+)', s.strip())
    if not m:
        fail(f'not a number in e-notation: {s!r}')
    d = len(m.group(2) or '') if nd is None else nd
    return f'{fmt(float(m.group(1)), d)} × 10{str(int(m.group(3))).translate(_SUP)}'


def wn(x, nd=1):
    """Wavenumber: imaginary modes (stored negative) as '3.2i'."""
    return f'{abs(x):.{nd}f}i' if x < 0 else f'{x:.{nd}f}'


def pct(x, nd=1):
    return fmt(x, nd)


def words(n):
    """Small counts in words, as in running text ('nine ring atoms'); 10 and above as digits."""
    return ('zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine')[n] if 0 <= n < 10 else str(n)


def hb(s):
    """'N4-H4N...O1' -> 'N4–H4N···O1' (en dash, three middle dots)."""
    s = s.replace('...', '···')
    return re.sub(r'(?<=[A-Za-z0-9])-(?=[A-Z])', NDASH, s)


FRACTIONS = {'1/2': '½', '3/2': '³⁄₂', '1/4': '¼', '3/4': '¾'}    # one glyph style in every code (QA r3, r3-m6)


def symop(s):
    """'x-1, -y+1/2, z-1/2' -> '*x* − 1, −*y* + ½, *z* − ½'."""
    comps = [c.replace(' ', '') for c in s.split(',')]
    if len(comps) != 3:
        fail(f'not a symmetry operation: {s!r}')
    out = []
    for c in comps:
        terms = re.findall(r'([+-]?)(\d+(?:/\d+)?|[xyz])', c)
        if ''.join(sg + t for sg, t in terms) != c:
            fail(f'cannot parse symmetry component {c!r} of {s!r}')
        txt = ''
        for k, (sg, t) in enumerate(terms):
            if '/' in t and t not in FRACTIONS:                      # no mixed styles (QA r3, r3-m6)
                fail(f'no glyph for the fraction {t!r} in {s!r}')
            t = FRACTIONS.get(t, f'*{t}*' if t in ('x', 'y', 'z') else t)
            if k == 0:
                txt += (MINUS if sg == '-' else '') + t
            else:
                txt += f'{NBSP}{MINUS if sg == "-" else "+"}{NBSP}{t}'   # a component never breaks across lines
        out.append(txt)
    return ', '.join(out)


def formula(s):
    """'C12 H8 N6 O7' or 'C6H2N3O7' -> 'C~12~H~8~N~6~O~7~'."""
    parts = re.findall(r'([A-Z][a-z]?)(\d*)', s.replace(' ', ''))
    if ''.join(e + n for e, n in parts) != s.replace(' ', ''):
        fail(f'cannot parse formula {s!r}')
    return ''.join(e + (f'~{n}~' if n and n != '1' else '') for e, n in parts)


def heavy_count(f):
    """Number of non-H atoms in a compact formula such as 'C6H2N3O7'."""
    return sum(int(n or 1) for e, n in re.findall(r'([A-Z][a-z]?)(\d*)', f) if e != 'H')


def spacegroup(s):
    m = re.fullmatch(r'([A-Z])(\d)_(\d)/([a-z])', s)
    if not m:
        fail(f'unexpected space-group symbol {s!r}')
    return f'*{m.group(1)}*{m.group(2)}~{m.group(3)}~/*{m.group(4)}*'


def tidy(s):
    """Typography of free-text facts values: arrows, ±, ≥, superscripts, ×, true minus, e-notation, ranges."""
    s = s.replace('->', '→').replace('+-', '±').replace('>=', '≥').replace('<=', '≤')
    s = re.sub(r'(\d)\^2\b', r'\1²', s)
    s = re.sub(r'(\d)\^3\b', r'\1³', s)
    s = re.sub(r'(?<=[\w)]) x (?=\d)', ' × ', s)
    s = s.replace('Treutler-Ahlrichs', 'Treutler–Ahlrichs')
    s = re.sub(r'(-?\d+(?:\.\d+)?)[eE]([-+]?\d+)', lambda m: sci_str(m.group(0)), s)
    s = re.sub(r'(?<=\d)\.\.(?=[-+]?\d)', ' to ', s)
    s = re.sub(r'(^|[\s(=,/])-(?=\d)', lambda m: m.group(1) + MINUS, s)
    s = re.sub(r'(?<=\d) deg\b', '°', s)
    return s


def _nat(label):
    m = re.match(r'([A-Za-z]+?)(\d+)(.*)$', label)
    return (m.group(1), int(m.group(2)), m.group(3)) if m else (label, 0, '')


def _elem(label):
    return re.match(r'[A-Z][a-z]?', label).group(0)


def _order(label):
    """Heavy atoms before H; within each, element letters then number (C1 < C2 < C12 < N1 < O1; H last)."""
    return (_elem(label) == 'H', _nat(label))


def bond_name(a, b):
    return NDASH.join(sorted((a, b), key=_order))


def angle_name(a, b, c):
    x, y = sorted((a, c), key=_order)
    return f'{x}{NDASH}{b}{NDASH}{y}'


# =================================================================================================
# table model, pipe-table writer and style lint
# =================================================================================================
@dataclass
class Table:
    stem: str
    where: str                     # 'main' or 'SI'
    title: str                     # short description for INDEX.md
    header: list
    align: str                     # one of l / r / c per column
    rows: list
    caption: str
    notes: list = field(default_factory=list)
    widths: list | None = None     # optional relative column widths (dash counts); default: content width
    keys: list = field(default_factory=list)
    files: list = field(default_factory=list)
    checks: list = field(default_factory=list)
    number: str = ''


REGISTRY = []
MAIN_ORDER = ['tbl_crystal_model', 'tbl_geometry', 'tbl_hbonds', 'tbl_vib_thermo', 'tbl_nci', 'tbl_fmo_descriptors',
              'tbl_mep']                   # citation order in manuscript.md (QA r4-m6)
STATS = {}                         # counts found while building, quoted in INDEX.md


def table(stem, where, title):
    def deco(fn):
        REGISTRY.append((stem, where, title, fn))
        return fn
    return deco


def _vis(s):
    return len(re.sub(r'[*~^\\]', '', s))


def _esc(c):
    return str(c).replace('|', r'\|')


def render(t: Table) -> str:
    ncol = len(t.header)
    if len(t.align) != ncol:
        fail(f'{len(t.align)} alignments for {ncol} columns')
    for r in t.rows:
        if len(r) != ncol:
            fail(f'row with {len(r)} cells in a {ncol}-column table: {r}')
    cells = [[_esc(c) for c in t.header]] + [[_esc(c) for c in r] for r in t.rows]
    widths = [max(3, min(WIDTH_CAP, max(len(r[i]) for r in cells))) for i in range(ncol)]

    def line(r):
        out = []
        for i, c in enumerate(r):
            pad = ' ' * max(widths[i] - len(c), 0)
            out.append(pad + c if t.align[i] == 'r' else c + pad)
        return '| ' + ' | '.join(out) + ' |'

    if t.widths:
        if len(t.widths) != ncol or min(t.widths) < 3:
            fail(f'bad widths {t.widths} for {ncol} columns')
        dashes = [w - 2 for w in t.widths]
    else:
        # relative widths (pandoc reads them from the dash counts): the widest body cell of each column, at least
        # the longest word of its header (headers may wrap, values should not), capped for long text columns
        def longest_word(s):
            return max((len(w) for w in re.sub(r'[*~^\\]', '', s).split()), default=0)
        dashes = [max(3, min(WIDTH_CAP, max([_vis(r[i]) for r in cells[1:]] or [0])), longest_word(cells[0][i]))
                  for i in range(ncol)]
    sep = []
    for a, w in zip(t.align, dashes):
        d = '-' * (w + 2)
        sep.append({'l': ':' + d[1:], 'r': d[:-1] + ':', 'c': ':' + d[2:] + ':'}[a])
    lines = [line(cells[0]), '|' + '|'.join(sep) + '|'] + [line(r) for r in cells[1:]]
    # notes and footnotes are part of the caption (round 1, figures_r1 #1): set as body paragraphs they were
    # line-numbered in the PDF and Word output and could be separated from their table by a float
    cap = ' '.join([t.caption] + [n.strip() for n in t.notes])
    text = '\n'.join(lines) + '\n\n' + f'Table: {cap} {{#tbl:{t.stem}}}\n'
    return text


def lint(text, stem):
    bad = []
    body = text.replace(f'{{#tbl:{stem}}}', '')
    checks = [
        (r'(?i)schiff', 'the word "Schiff"'),
        (r'#', "a '#' (symmetry suffix or stray label)"),
        ('—', 'an em dash'),
        (r'\.\.\.', 'three periods instead of ···'),
        (r'(?<=\d)\.\.(?=[-+−]?\d)', "an ASCII '..' range"),
        (r'->|\+-|>=|<=', 'an ASCII arrow / ± / ≥ / ≤'),
        (r'(?m)(?:^|[\s(=,;/:])-(?=\d)', 'a hyphen used as a minus sign'),
        (r'\d[eE][-+]\d', 'e-notation'),
        (r'(?<=[\w)]) x (?=\d)', "an 'x' used as a multiplication sign"),
        (r'(?<=\d) deg\b', "'deg' instead of °"),
        (r'\bkcal/mol|cm-1\b|\bA\^', 'an ASCII unit'),
        (r'\{#', 'a second label'),
    ]
    for pat, what in checks:
        m = re.search(pat, body)
        if m:
            ctx = body[max(0, m.start() - 40):m.end() + 40].replace('\n', ' ')
            bad.append(f'{what}: ...{ctx}...')
    for ln in body.splitlines():
        for cell in (ln.split(' | ') if ln.startswith('|') else [ln]):
            c = cell.replace(r'\~', '')
            if c.count('~') % 2 or c.count('^') % 2:
                bad.append(f'unbalanced ~ or ^ markup: {cell!r}')
    if bad:
        fail('style lint:\n  ' + '\n  '.join(bad))


# =================================================================================================
# shared structure data
# =================================================================================================
def case_json(model):
    return json.loads(use_file(MOL / 'cases' / f'{CASE[model]}.json').read_text())


def labels_and_ions():
    """Clean CIF labels (atom order of every xyz) and label -> 'picrate' | 'cation'."""
    cj = case_json('pcm')
    labels = geo.clean_labels(cj['cif_labels'])
    check(geo.clean_labels(case_json('gas')['cif_labels']) == labels, 'gas and IEF-PCM case files list the same atom labels')
    ion = {}
    for name, idx in cj['fragments'].items():
        short = 'picrate' if 'picrate' in name else ('cation' if 'cation' in name else None)
        if short is None:
            fail(f'unexpected fragment name {name!r}')
        for i in idx:
            ion[labels[i]] = short
    check(len(ion) == len(labels), 'every atom belongs to one of the two ions')
    return labels, ion


def read_geom(path, labels):
    fr = geo.read_xyz(use_file(path))
    check(len(fr.symbols) == len(labels), f'{Path(path).name}: {len(labels)} atoms')
    for s, lab in zip(fr.symbols, labels):
        if _elem(lab) != s:
            fail(f'{path}: element {s} does not match label {lab}')
    return fr


def xray_xyz_path():
    return REPO / case_json('pcm')['xyz']


def bond_rule():
    """Covalent radii and factor of the X-ray connectivity used by build_facts.py (read from its source)."""
    src = use_file(BUILD_FACTS).read_text()
    m = re.search(r'^COV = (\{[^}]*\})', src, re.M)
    f = re.search(r'< ([\d.]+) \* \(COV\[SYM\[i\]\] \+ COV\[SYM\[j\]\]\)', src)
    if not m or not f:
        fail('cannot read the bond rule (COV, factor) from build_facts.py')
    return ast.literal_eval(m.group(1)), float(f.group(1))


def connectivity(sym, X):
    cov, factor = bond_rule()
    n = len(sym)
    bonds = [(i, j) for i in range(n) for j in range(i + 1, n)
             if np.linalg.norm(X[i] - X[j]) < factor * (cov[sym[i]] + cov[sym[j]])]
    nb = {i: set() for i in range(n)}
    for i, j in bonds:
        nb[i].add(j)
        nb[j].add(i)
    angles = [(i, j, k) for j in range(n) for i in sorted(nb[j]) for k in sorted(nb[j]) if i < k]
    return bonds, angles


def dist(X, i, j):
    return float(np.linalg.norm(X[i] - X[j]))


def ang(X, i, j, k):
    v1, v2 = X[i] - X[j], X[k] - X[j]
    return math.degrees(math.acos(float(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2))))


def cif_su(key, unit):
    """CIF value with its standard uncertainty, e.g. '1.233(3)', checked against the key's value."""
    v = fnum(key, unit)
    m = from_note(key, r'with s\.u\.:? (-?[\d.]+)\((\d+)\)')
    close(float(m.group(1)), v, 1e-9, f'{key}: s.u. string matches the value')
    return f'{fmt(float(m.group(1)), len(m.group(1).split(".")[1]) if "." in m.group(1) else 0)}({m.group(2)})'


def method():
    f = from_value('protocol.functional', r'^(B3LYP)\b').group(1)
    b = from_value('protocol.basis', r'^(6-311\+\+G\(d,p\))').group(1)
    return f'{f}/{b}'


# =================================================================================================
# MAIN-TEXT TABLES
# =================================================================================================
@table('tbl_crystal_model', 'main', 'Crystal data (CIF) and the ion-pair model')
def t_crystal_model():
    ccdc = fint('compound.ccdc')
    rows = [['*Crystal structure*', '']]
    rows.append(['Sum formula', formula(fstr('compound.formula_sum'))])
    m = from_note('compound.formula_moiety', r'picrate anion (\w+)\(-\), benzotriazolium cation (\w+)\(\+\)')
    moieties = set(fstr('compound.formula_moiety').replace(' ', '').split(','))
    check(moieties == {m.group(1), m.group(2)}, 'moiety formula of the CIF = the two ions named in its note')
    rows.append(['Ions', f'{formula(m.group(2))}^+^ (benzotriazolium), {formula(m.group(1))}^{MINUS}^ (picrate)'])
    # relative molar mass: dimensionless (QA round 2, r2-m11); the CIF formula weight in g mol-1 has the same value
    rows.append(['*M*~r~', fmt(fnum('xray.mr', 'g mol-1'), 2)])
    rows.append(['Crystal system', fstr('xray.crystal_system')])
    rows.append(['Space group', spacegroup(fstr('xray.space_group'))])
    rows.append(['*a*, *b*, *c* (Å)', ', '.join(cif_su(f'xray.cell_{k}', 'Å') for k in 'abc')])
    rows.append(['*β* (°)', cif_su('xray.cell_beta', 'deg')])
    rows.append(['*V* (Å³)', cif_su('xray.cell_volume', 'Å^3')])
    rows.append(['*Z*, *Z*′', f"{fint('xray.z')}, {fint('xray.z_prime')}"])
    t = from_note('xray.temperature', r'CIF: (\d+)\((\d+)\) K')
    close(float(t.group(1)), fnum('xray.temperature', 'K'), 1e-9, 'temperature s.u. string matches the value')
    rows.append(['*T* (K)', f'{t.group(1)}({t.group(2)})'])
    rows.append(['Radiation, *λ* (Å)', f"{fstr('xray.radiation')}, {fmt(fnum('xray.wavelength', 'Å'), 5)}"])
    rows.append(['*D*~calc~ (g cm⁻³)', fmt(fnum('xray.density_calc', 'g cm-3'), 3)])
    rows.append(['*μ* (mm⁻¹)', fmt(fnum('xray.mu', 'mm-1'), 3)])
    rows.append(['Reflections: measured, unique, *I* > 2*σ*(*I*)',
                 ', '.join(str(fint(f'xray.reflections_{k}')) for k in ('measured', 'unique', 'observed'))])
    rows.append(['*R*~int~', fmt(fnum('xray.r_int', ''), 4)])
    rows.append(['*R*~1~, *wR*~2~ [*I* > 2*σ*(*I*)]', f"{fmt(fnum('xray.r1_gt', ''), 4)}, {fmt(fnum('xray.wr2_gt', ''), 4)}"])
    rows.append(['Goodness of fit', fmt(fnum('xray.gof', ''), 3)])
    ht = fstr('xray.hydrogen_treatment')
    check('refined coordinates' in ht and 'refined Uiso' in ht, 'CIF: every H atom has refined coordinates and Uiso')
    rows.append(['H atoms', 'coordinates and *U*~iso~ refined'])
    rows.append(['Hydrogen-bond table in the CIF', fstr('xray.hbond_table_deposited')])
    rows.append(['*Ion-pair model*', ''])
    rows.append(['Picrate at', symop(from_value('model.symmetry_picrate', r'\((.+)\)$').group(1))])
    rows.append(['Cation at', symop(from_value('model.symmetry_cation', r'\((.+)\)$').group(1))])
    d = from_value('model.defining_hbond', r'^(\S+)\s').group(1)
    check(d == 'N4-H4N...O1', 'defining hydrogen bond of the model is N4-H4N...O1')
    rows.append(['Hydrogen bond between the ions',
                 f"{hb(d)} (N···O {fmt(fnum('xray.hb_N4_H4N_O1_da', 'Å'), 3)} Å in the crystal)"])
    mh = from_value('model.missing_hbond', r'^(\S+)\s').group(1)
    check(mh == 'N6-H6N...O6' and fstr('xray.hb02_N6_H6N_O6') == mh, 'contact missing from the model = crystal contact hb02')
    from_note('model.missing_hbond', r'2_1-screw image')
    rows.append(['Crystal hydrogen bond not in the model',
                 f"{hb(mh)} to a second picrate (N···O {fmt(fnum('xray.hb02_N6_H6N_O6_da', 'Å'), 3)} Å)"])
    rows.append(['Atoms, electrons', f"{fint('model.n_atoms')}, {fint('model.n_electrons')}"])
    rows.append(['Charge, spin multiplicity', f"{fint('model.charge', 'e')}, {fint('model.multiplicity')}"])
    rows.append(['Basis functions', str(fint('model.n_basis'))])
    rows.append(['Level of theory', method()])
    check(fstr('protocol.solvation_model') == 'IEF-PCM' and fstr('protocol.solvent') == 'water', 'solvation model is IEF-PCM, water')
    rows.append(['Models', f"gas phase and IEF-PCM (water, *ε* = {fmt(fnum('protocol.solvent_eps', ''), 4)})"])
    crit = from_value('xray.hb_search_criteria', r'H\.\.\.A < ([\d.]+) Å and D-H\.\.\.A > ([\d.]+) deg')
    notes = [
        f"Crystal data from the CIF deposited with the original structure report (<https://doi.org/{fstr('compound.source_doi')}>); "
        'standard uncertainties in the last digit are given in parentheses. The CIF contains no hydrogen-bond table, so '
        'the N···O distances were derived here from the deposited coordinates and symmetry operations '
        f'(H···A < {crit.group(1)} Å, D–H···A > {crit.group(2)}°). Positions of the two ions refer to the deposited '
        'asymmetric unit. In the model the second N–H group of the cation, N6–H6N, is free.'
    ]
    return dict(header=['Quantity', 'Value'], align='ll', rows=rows, notes=notes,
                caption=f'Crystal data of benzotriazolium picrate (CCDC {ccdc}) and the ion-pair model used in the calculations.')


# bonds of Table 2 in the order of the facts keys r_<a>_<b>; value = CIF key with s.u. where one exists
T2_BONDS = [('O1', 'C1', 'bond_O1_C1'), ('C1', 'C2', 'bond_C1_C2'), ('C1', 'C6', 'bond_C6_C1'), ('C2', 'C3', 'bond_C2_C3'),
            ('C5', 'C6', 'bond_C6_C5'), ('C2', 'N1', 'bond_C2_N1'), ('C4', 'N2', 'bond_C4_N2'), ('C6', 'N3', 'bond_N3_C6'),
            ('N1', 'O2', 'bond_N1_O2'), ('N1', 'O3', 'bond_N1_O3'), ('N2', 'O4', 'bond_N2_O4'), ('N2', 'O5', 'bond_N2_O5'),
            ('N3', 'O6', 'bond_N3_O6'), ('N3', 'O7', 'bond_N3_O7'), ('N4', 'N5', 'bond_N4_N5'), ('N5', 'N6', 'bond_N6_N5'),
            ('N4', 'C12', 'bond_N4_C12'), ('N6', 'C7', 'bond_N6_C7'), ('C7', 'C12', 'bond_C7_C12'), ('N4', 'H4N', 'bond_N4_H4N'),
            ('N6', 'H6N', 'bond_N6_H6N')]
T2_ANGLES = [('C2', 'C1', 'C6', 'angle_C2_C1_C6'), ('N4', 'N5', 'N6', 'angle_N4_N5_N6'),
             ('N5', 'N4', 'C12', 'angle_N5_N4_C12'), ('N5', 'N6', 'C7', 'angle_N5_N6_C7')]


@table('tbl_geometry', 'main', 'Selected bond lengths, angles, ring orientation and agreement statistics vs X-ray')
def t_geometry():
    rows = [['*Bond lengths (Å)*', '', '', '']]
    for a, b, ck in T2_BONDS:
        rc = fnum(f'xray.r_{a}_{b}', 'Å')
        if ck:
            x = cif_su(f'xray.{ck}', 'Å')
            close(fnum(f'xray.{ck}', 'Å'), rc, 0.006 if 'H' in (_elem(a), _elem(b)) else 0.002,
                  f'CIF bond {a}-{b} vs the X-ray cluster coordinates')
        else:
            x = fmt(rc, 3)
        rows.append([bond_name(a, b), x, fmt(fnum(f'gas.r_{a}_{b}', 'Å'), 3), fmt(fnum(f'pcm.r_{a}_{b}', 'Å'), 3)])
    rows.append(['*Bond angles (°)*', '', '', ''])
    for a, b, c, ck in T2_ANGLES:
        x = cif_su(f'xray.{ck}', 'deg')
        close(fnum(f'xray.{ck}', 'deg'), fnum(f'xray.a_{a}_{b}_{c}', 'deg'), 0.15, f'CIF angle {a}-{b}-{c} vs the X-ray cluster')
        rows.append([angle_name(a, b, c), x, fmt(fnum(f'gas.a_{a}_{b}_{c}', 'deg'), 1), fmt(fnum(f'pcm.a_{a}_{b}_{c}', 'deg'), 1)])
    rows.append(['*Ring orientation*', '', '', ''])
    ip = from_note('pcm.interplanar_picrate_cation', r'least-squares planes of the picrate C6 ring and the (\d+) ring atoms of the cation')
    rows.append(['Angle between the ring planes (°)'] +
                [fmt(fnum(f'{m}.interplanar_picrate_cation', 'deg'), 1) for m in ('xray', 'gas', 'pcm')])
    rows.append(['Ring centroid distance (Å)'] +
                [fmt(fnum(f'{m}.centroid_distance_picrate_cation', 'Å'), 3) for m in ('xray', 'gas', 'pcm')])
    for k in (1, 2, 3):
        pos = from_note(f'pcm.nitro_twist_N{k}', r'\((ortho|para), on (C\d)\b')
        from_note(f'gas.nitro_twist_N{k}', rf'\({pos.group(1)}, on {pos.group(2)}\b')
        rows.append([f'Nitro twist, N{k} (*{pos.group(1)}*, on {pos.group(2)}) (°)'] +
                    [fmt(fnum(f'{m}.nitro_twist_N{k}', 'deg'), 1) for m in ('xray', 'gas', 'pcm')])
    nb = int(from_note('pcm.bond_rmsd_heavy', r'(\d+) bonds between non-H atoms').group(1))
    na = int(from_note('pcm.angle_rmsd_heavy', r'(\d+) angles with no H atom').group(1))
    nall = int(from_note('pcm.kabsch_rmsd_heavy', r'all (\d+) heavy atoms').group(1))
    fr = from_value('model.fragments', r'picrate anion: \d+ atoms .*?, (C\d+H\d+N\d+O\d+)\); benzotriazolium cation: \d+ atoms .*?, (C\d+H\d+N\d+)\)')
    npic, ncat = heavy_count(fr.group(1)), heavy_count(fr.group(2))
    check(npic + ncat == nall, f'picrate ({npic}) + cation ({ncat}) heavy atoms = {nall}')
    for m in ('gas',):
        from_note(f'{m}.bond_rmsd_heavy', rf'\b{nb} bonds between non-H atoms')
        from_note(f'{m}.angle_rmsd_heavy', rf'\b{na} angles with no H atom')
    rows.append(['*Agreement with X-ray*', '', '', ''])
    stat = [(f'Bond RMSD, {nb} heavy-atom bonds (Å)', 'bond_rmsd_heavy', 'Å', 4),
            (f'Bond MAE, {nb} heavy-atom bonds (Å)', 'bond_mae_heavy', 'Å', 4),
            (f'Angle RMSD, {na} heavy-atom angles (°)', 'angle_rmsd_heavy', 'deg', 2),
            (f'Superposition RMSD, {nall} heavy atoms (Å)', 'kabsch_rmsd_heavy', 'Å', 3),
            (f'Superposition RMSD, {npic} picrate heavy atoms (Å)', 'kabsch_rmsd_picrate_heavy', 'Å', 3),
            (f'Superposition RMSD, {ncat} cation heavy atoms (Å)', 'kabsch_rmsd_cation_heavy', 'Å', 3)]
    for lab, k, u, nd in stat:
        rows.append([lab, NDASH, fmt(fnum(f'gas.{k}', u), nd), fmt(fnum(f'pcm.{k}', u), nd)])
    nall_b = int(from_note('pcm.bond_rmsd_all', r'(\d+) bonds incl').group(1))
    notes = [
        'X-ray bond lengths and angles are the CIF values, with standard uncertainties in parentheses; the '
        'ring-orientation values were computed from the deposited coordinates. X-ray N–H lengths are short, as usual '
        'for hydrogen positions from X-ray '
        f'diffraction. Ring planes are least-squares planes of the six picrate ring carbons and the {words(int(ip.group(1)))} ring atoms '
        'of the cation; a nitro twist is the angle between the O–N–O plane and the picrate ring plane. RMSD and MAE '
        'refer to DFT − X-ray over the bonds and angles of the X-ray connectivity that contain no hydrogen; over all '
        f"{nall_b} bonds, X–H included, the bond RMSD is {fmt(fnum('gas.bond_rmsd_all', 'Å'), 4)} Å (gas phase) and "
        f"{fmt(fnum('pcm.bond_rmsd_all', 'Å'), 4)} Å (IEF-PCM). Superposition RMSDs are for a best fit of the atoms named."
    ]
    T = fnum('xray.temperature', 'K')
    return dict(header=['Parameter', MODEL['xray'], MODEL['gas'], MODEL['pcm']], align='lrrr', rows=rows, notes=notes,
                caption=f'Selected bond lengths (Å), bond angles (°) and ring orientation of benzotriazolium picrate from X-ray '
                        f'diffraction ({T:.0f} K) and from {method()} optimisation in the gas phase and in IEF-PCM (water).')


@table('tbl_hbonds', 'main', 'Hydrogen-bond geometry: crystal and both models')
def t_hbonds():
    rows = []
    cross = {('N4', 'H4N', 'O1'): 'hb01', ('N4', 'H4N', 'O7'): 'hb03', ('C11', 'H11', 'O2'): 'hb09'}
    for d, h, a in [('N4', 'H4N', 'O1'), ('N4', 'H4N', 'O7'), ('C11', 'H11', 'O2'), ('C11', 'H11', 'O1')]:
        name = f'{d}{NDASH}{h}···{a}'
        if (d, h, a) in cross:
            c = cross[(d, h, a)]
            check(fstr(f'xray.{c}_{d}_{h}_{a}') == f'{d}-{h}...{a}', f'crystal contact {c} is {name}')
            for q in ('dh', 'ha', 'da', 'angle'):
                close(fnum(f'xray.hb_{d}_{h}_{a}_{q}'), fnum(f'xray.{c}_{d}_{h}_{a}_{q}'), 1e-6,
                      f'{name} {q}: X-ray cluster = crystal search')
        for k, m in enumerate(('xray', 'gas', 'pcm')):
            p = f'{m}.hb_{d}_{h}_{a}'
            rows.append([name if k == 0 else '', MODEL[m], fmt(fnum(p + '_dh', 'Å'), 3), fmt(fnum(p + '_ha', 'Å'), 3),
                         fmt(fnum(p + '_da', 'Å'), 3), fmt(fnum(p + '_angle', 'deg'), 1)])
    ops = {}
    for c, d, h, a in [('hb01', 'N4', 'H4N', 'O1'), ('hb03', 'N4', 'H4N', 'O7'), ('hb09', 'C11', 'H11', 'O2'),
                       ('hb02', 'N6', 'H6N', 'O6'), ('hb07', 'C9', 'H9', 'O5')]:
        k = f'xray.{c}_{d}_{h}_{a}'
        m = from_note(k, r'acceptor at symmetry \d_\d{3} = \(([^)]*)\); donor at x,y,z; intermolecular; (present in|not in) the DFT cluster')
        ops[a] = (m.group(1), m.group(2) == 'present in')
    check(ops['O1'][0] == ops['O7'][0] == ops['O2'][0] and ops['O1'][1], 'O1, O7 and O2 belong to the modelled picrate')
    check(not ops['O6'][1] and not ops['O5'][1], 'the O6 and O5 acceptors are not in the model')
    for c, d, h, a in [('hb02', 'N6', 'H6N', 'O6'), ('hb07', 'C9', 'H9', 'O5')]:
        p = f'xray.{c}_{d}_{h}_{a}'
        check(fstr(p) == f'{d}-{h}...{a}', f'crystal contact {c} is {d}-{h}...{a}')
        rows.append([f'{d}{NDASH}{h}···{a}', 'X-ray only', fmt(fnum(p + '_dh', 'Å'), 3), fmt(fnum(p + '_ha', 'Å'), 3),
                     fmt(fnum(p + '_da', 'Å'), 3), fmt(fnum(p + '_angle', 'deg'), 1)])
    asym_x = fnum('xray.nh_minus_ho_asym', 'Å')                     # QA round 3, r3-m1: the key, not _ha - _dh
    hcut = float(from_value('xray.hb_search_criteria', r'H\.\.\.A < ([\d.]+) Å').group(1))
    for m in ('xray', 'pcm', 'gas'):
        close(fnum(f'{m}.nh_minus_ho_asym', 'Å'), fnum(f'{m}.hb_N4_H4N_O1_ha', 'Å') - fnum(f'{m}.hb_N4_H4N_O1_dh', 'Å'),
              1e-9, f'{m}: asymmetry key = H...O1 - N4-H4N')
    notes = [
        'D = donor, A = acceptor. X-ray contacts were derived from the deposited coordinates. Relative to a cation at '
        f"*x*, *y*, *z*, the acceptors O1, O7 and O2 belong to the picrate at {symop(ops['O1'][0])}, the one in the model; "
        f"O6 belongs to the picrate at {symop(ops['O6'][0])} and O5 to the picrate at {symop(ops['O5'][0])}, neither of "
        'which is part of the model. X-ray N–H and C–H lengths are short, so *d*(D···A) is the reliable crystal quantity. '
        f"The proton asymmetry *d*(H4N···O1) − *d*(N4–H4N) is {fmt(asym_x, 3)} Å (X-ray), "
        f"{fmt(fnum('gas.nh_minus_ho_asym', 'Å'), 3)} Å (gas phase) and {fmt(fnum('pcm.nh_minus_ho_asym', 'Å'), 3)} Å "
        '(IEF-PCM); zero would correspond to a centred proton. '
        # figures_r2 #5: C11-H11...O1 is in the table only because the gas-phase NCI analysis finds a surface there
        f"C11–H11···O1 is listed for comparison with the NCI analysis (Section 3.5), which finds a weak surface there in "
        f"the gas phase; in the crystal "
        f"({fmt(fnum('xray.hb_C11_H11_O1_ha', 'Å'), 3)} Å) and in IEF-PCM ({fmt(fnum('pcm.hb_C11_H11_O1_ha', 'Å'), 3)} Å) "
        f'its H···O distance exceeds the {hcut:g} Å limit of the contact search.'
    ]
    check(fnum('xray.hb_C11_H11_O1_ha', 'Å') > hcut and fnum('pcm.hb_C11_H11_O1_ha', 'Å') > hcut,
          'C11-H11...O1 fails the 2.9 Å contact criterion in the crystal and in IEF-PCM')
    # figures_r2 #7: units in the column headers
    return dict(header=['D–H···A', 'Structure', '*d*(D–H) (Å)', '*d*(H···A) (Å)', '*d*(D···A) (Å)', '∠D–H···A (°)'],
                align='llrrrr', rows=rows, notes=notes,
                caption='Hydrogen-bond and contact geometry in the crystal and in the ion pairs '
                        f'optimised at the {method()} level in the gas phase and in IEF-PCM (water).')


def lowest_modes(m):
    lst = [float(x) for x in fstr(f'{m}.freq_lowest10', 'cm-1').split(',')]
    first = fnum(f'{m}.freq_lowest', 'cm-1')
    close(lst[0], first, 0.006, f'{m}: first of freq_lowest10 = freq_lowest')
    return [first] + lst[1:]


def wag_list(key):
    return [float(x) for x in fstr(key, 'cm-1').split(',')]


def step_005():
    """Displacement step of the first IEF-PCM Hessian, read from the note of saddle.freq_imag_step0p005."""
    return from_note('saddle.freq_imag_step0p005', r'step ([\d.]+) Bohr').group(1)


@table('tbl_vib_thermo', 'main', 'Harmonic vibrational data (N–H modes, lowest modes) and thermochemistry')
def t_vib_thermo():
    ha_kcal = float(from_note('pcm.zpe_kcal', r'x ([\d.]+)').group(1))
    TK = fnum('protocol.thermo_T', 'K')
    p = fnum('protocol.thermo_p', 'Pa')
    check(p == 101325, 'thermochemistry pressure is 101325 Pa (1 atm)')
    rows = []
    def row(label, f):
        rows.append([label, f('gas'), f('pcm')])
    row('*E* (Ha)', lambda m: fmt(fnum(f'{m}.e_total', 'Ha'), 8))
    row('Finite-difference step of the Hessian (Bohr)',
        lambda m: from_value(f'{m}.hessian_method', r'step ([\d.]+) Bohr').group(1))
    for m in ('pcm', 'gas'):
        lm = lowest_modes(m)
        check(sum(x < 0 for x in lm) == fint(f'{m}.n_imag'), f'{m}: n_imag = negative entries among the lowest modes')
    row('Imaginary modes', lambda m: str(fint(f'{m}.n_imag')))
    row('Three lowest wavenumbers', lambda m: ', '.join(wn(x) for x in lowest_modes(m)[:3]))
    row('Rigid-body share of these three modes',
        lambda m: ', '.join(fmt(fnum(f'{m}.mode{k}_rigid_fraction', ''), 3) for k in (1, 2, 3)))
    row('Modes below 100 cm⁻¹', lambda m: str(fint(f'{m}.n_modes_below_100')) + ('^a^' if m == 'pcm' else ''))
    w0 = lambda x: fmt(x, 0)                      # N-H rows to the nearest cm-1, as in the text (figures_r1 #15)
    row('N4–H4N stretch (hydrogen-bonded)', lambda m: w0(fnum(f'{m}.nu_N4_H4N_stretch', 'cm-1')))
    row('N4–H4N in-plane bend', lambda m: w0(fnum(f'{m}.nu_N4_H4N_inplane_bend', 'cm-1')))
    row('N4–H4N out-of-plane wag', lambda m: ', '.join(w0(x) for x in wag_list(f'{m}.nu_N4_H4N_oop_wag')))
    row('N6–H6N stretch (free)', lambda m: w0(fnum(f'{m}.nu_N6_H6N_stretch', 'cm-1')))
    row('N6–H6N in-plane bend', lambda m: w0(fnum(f'{m}.nu_N6_H6N_inplane_bend', 'cm-1')))
    row('N6–H6N out-of-plane wags', lambda m: ', '.join(w0(x) for x in wag_list(f'{m}.nu_N6_H6N_oop_wag')))
    for m in ('pcm', 'gas'):
        close(fnum(f'diff.nh_redshift_{m}', 'cm-1'),
              fnum(f'{m}.nu_N6_H6N_stretch', 'cm-1') - fnum(f'{m}.nu_N4_H4N_stretch', 'cm-1'), 1e-6,
              f'{m}: red shift = free N6-H6N minus bonded N4-H4N stretch')
    row('*ν*(N6–H6N) − *ν*(N4–H4N)', lambda m: fmt(fnum(f'diff.nh_redshift_{m}', 'cm-1'), 0))
    for m in ('pcm', 'gas'):
        close(fnum(f'{m}.zpe', 'Ha') * ha_kcal, fnum(f'{m}.zpe_kcal', 'kcal mol-1'), 1e-6, f'{m}: ZPE in kcal/mol')
        close(fnum(f'{m}.e_total', 'Ha') + fnum(f'{m}.thermal_corr_h', 'Ha'), fnum(f'{m}.h298', 'Ha'), 2e-6, f'{m}: H = E + (H - E)')
        close(fnum(f'{m}.e_total', 'Ha') + fnum(f'{m}.thermal_corr_g', 'Ha'), fnum(f'{m}.g298', 'Ha'), 2e-6, f'{m}: G = E + (G - E)')
        close((fnum(f'{m}.h298', 'Ha') - fnum(f'{m}.g298', 'Ha')) * ha_kcal * 1000 / TK,
              fnum(f'{m}.entropy', 'cal mol-1 K-1'), 1e-3, f'{m}: S = (H - G)/T')
    row('Zero-point energy (kcal mol⁻¹)', lambda m: fmt(fnum(f'{m}.zpe_kcal', 'kcal mol-1'), 2))
    row('*H* − *E* (kcal mol⁻¹)', lambda m: fmt(fnum(f'{m}.thermal_corr_h', 'Ha') * ha_kcal, 2))
    row('*G* − *E* (kcal mol⁻¹)', lambda m: fmt(fnum(f'{m}.thermal_corr_g', 'Ha') * ha_kcal, 2))
    for m in ('gas', 'pcm'):
        close(fnum(f'{m}.entropy_rrho_check', 'cal mol-1 K-1'), fnum(f'{m}.entropy', 'cal mol-1 K-1'), 0.05,
              f'{m}: RRHO entropy recomputed for the quasi-RRHO treatment = (H - G)/T')
    row('*S*, harmonic (cal mol⁻¹ K⁻¹)', lambda m: fmt(fnum(f'{m}.entropy', 'cal mol-1 K-1'), 1))
    row('*S*, quasi-RRHO (cal mol⁻¹ K⁻¹)^b^', lambda m: fmt(fnum(f'{m}.entropy_qrrho', 'cal mol-1 K-1'), 1))
    nreal = from_note('pcm.n_modes_below_100', r'sum over (\d+) real modes').group(1)
    rigid = from_note('pcm.mode1_rigid_fraction', r'(\d+)-dim space of rigid translations/rotations of the two ions').group(1)
    from_note('protocol.thermo_p', r'no low-frequency correction')
    check(fint('pcm.n_imag_pyscf') == 0 and fint('pcm.n_imag') == 1, 'PySCF counted 0 imaginary modes, the projection 1')
    notes = [
        'Imaginary modes are written with i. No infrared or Raman intensities were computed. Differences are taken '
        'before rounding. Rigid-body share: fraction of the mass-weighted normal mode lying in the '
        f'{rigid}-dimensional space of translations and rotations of the two ions. The thermal terms use the '
        'rigid-rotor, harmonic-oscillator model; because the soft modes below 100 cm⁻¹ enter it as harmonic '
        'oscillators, *S* and *G* are the least certain entries of the table.',
        # figures_r2 #6 / language_r2 #17: the step test and PySCF's count are in Section 2.4, not repeated here
        f'^a^ Includes the imaginary mode (Section 2.4); the IEF-PCM thermochemistry uses the {nreal} real modes.',
        '^b^ Each vibrational term interpolated between the harmonic-oscillator and a free-rotor entropy, with the '
        'weight 1/[1 + (100 cm⁻¹/*ν*)⁴], after Grimme; translation and rotation as in the harmonic value.'
    ]
    return dict(header=['Quantity', MODEL['gas'], MODEL['pcm']], align='lrr', rows=rows, notes=notes,
                caption='Harmonic vibrational data (unscaled wavenumbers in cm⁻¹) and thermochemistry at '
                        f'{TK:.2f} K and {p / 101325:g} atm of the benzotriazolium picrate ion pair, {method()}.')


@table('tbl_fmo_descriptors', 'main', 'Frontier orbitals, Koopmans-type descriptors, fragment character, dipole, charge')
def t_fmo():
    for m in ('pcm', 'gas'):
        g = lambda k, u='eV': fnum(f'{m}.{k}', u)
        close(g('gap_ev'), g('lumo_ev') - g('homo_ev'), 1e-9, f'{m}: gap = LUMO - HOMO')
        close(g('ip_ev'), -g('homo_ev'), 1e-9, f'{m}: IP = -E_HOMO')
        close(g('ea_ev'), -g('lumo_ev'), 1e-9, f'{m}: EA = -E_LUMO')
        close(g('chi_ev'), (g('ip_ev') + g('ea_ev')) / 2, 1e-9, f'{m}: chi = (IP + EA)/2')
        close(g('mu_ev'), -g('chi_ev'), 1e-9, f'{m}: mu = -chi')
        close(g('eta_ev'), (g('ip_ev') - g('ea_ev')) / 2, 1e-9, f'{m}: eta = (IP - EA)/2')
        close(g('softness_ev_inv', 'eV-1'), 1 / (2 * g('eta_ev')), 1e-9, f'{m}: S = 1/(2 eta)')
        close(g('omega_ev'), g('mu_ev') ** 2 / (2 * g('eta_ev')), 1e-9, f'{m}: omega = mu^2/(2 eta)')
        close(g('ct_gap_picrate_to_cation_ev'), g('cation_lowest_virt_ev') - g('picrate_highest_occ_ev'), 1e-9,
              f'{m}: CT gap = cation lowest virtual - picrate highest occupied')
    fm = fstr('protocol.fmo_descriptors')
    for part in ('IP = -E_HOMO', 'EA = -E_LUMO', 'chi = (IP+EA)/2', 'mu = -chi', 'eta = (IP-EA)/2', 'S = 1/(2 eta)',
                 'omega = mu^2/(2 eta)'):
        check(part in fm, f'protocol.fmo_descriptors contains "{part}"')
    rows = []
    def row(label, f):
        rows.append([label, f('gas'), f('pcm')])
    e3 = lambda k: (lambda m: fmt(fnum(f'{m}.{k}', 'eV'), 3))
    row('*E*~HOMO~ (eV)', e3('homo_ev'))
    row('*E*~LUMO~ (eV)', e3('lumo_ev'))
    row('Δ*E*~gap~ (eV)', e3('gap_ev'))

    def character(m):
        h, l = fnum(f'{m}.homo_picrate_pct', '%'), fnum(f'{m}.lumo_picrate_pct', '%')
        check(h > 90 and (l > 90 or l < 10), f'{m}: HOMO and LUMO are each localised on one ion')
        return 'within the picrate' if l > 50 else 'picrate → cation'
    row('Character of the HOMO–LUMO gap', character)
    row('IP (eV)', e3('ip_ev'))
    row('EA (eV)', e3('ea_ev'))
    row('*χ* (eV)', e3('chi_ev'))
    row('*μ* (eV)', e3('mu_ev'))
    row('*η* (eV)', e3('eta_ev'))
    row('*S* (eV⁻¹)', lambda m: fmt(fnum(f'{m}.softness_ev_inv', 'eV-1'), 3))
    row('*ω* (eV)', e3('omega_ev'))
    row('Picrate share of the HOMO (%)', lambda m: pct(fnum(f'{m}.homo_picrate_pct', '%')))
    row('Picrate share of the LUMO (%)', lambda m: pct(fnum(f'{m}.lumo_picrate_pct', '%')))
    row('Picrate → cation gap (eV)', e3('ct_gap_picrate_to_cation_ev'))
    row('Gap within the picrate (eV)', e3('local_gap_picrate_ev'))
    row('Gap within the cation (eV)', e3('local_gap_cation_ev'))
    row('Dipole moment (D)', lambda m: fmt(fnum(f'{m}.dipole_debye', 'D'), 2))
    row('Mulliken charge of the picrate (e)', lambda m: fmt(fnum(f'{m}.mulliken_charge_picrate', 'e'), 3))
    row('Löwdin charge of the picrate (e)', lambda m: fmt(fnum(f'{m}.lowdin_charge_picrate', 'e'), 3))
    from_value('protocol.fragment_composition_method', r'^Mulliken partition of each MO over the two fragments')
    from_note('pcm.ct_gap_picrate_to_cation_ev', r'E\(lowest cation-localised virtual\) - E\(highest picrate-localised occupied\)')
    notes = [
        'IP and EA are *E*~HOMO~ and *E*~LUMO~ with the sign reversed; *χ*, *μ*, *η*, *S* and *ω* follow from '
        '[@eq:chi_eta; @eq:s_omega], and a softness defined as 1/*η* would double the *S* row. Fragment shares come from a '
        'Mulliken partition of each orbital over the basis functions of the two ions. The gaps within each ion and '
        'the picrate → cation gap are taken between the highest occupied and the lowest unoccupied orbitals located '
        'on the ions named; Table S7 lists the orbitals from HOMO−7 to LUMO+7. In IEF-PCM the orbital energies include '
        'the reaction field of the continuum. Fragment charges obtained with a diffuse basis set are indicative only; '
        'the Löwdin charges are given as a second partition.'
    ]
    return dict(header=['Quantity', MODEL['gas'], MODEL['pcm']], align='lrr', rows=rows, notes=notes,
                caption='Frontier-orbital energies, Koopmans-type descriptors, fragment character of the orbitals, dipole '
                        f'moment and picrate charge of the ion pair, {method()}.')


def mep_extrema(m):
    """All surface extrema of one model from the facts notes: list of dicts (sorted minima, then maxima, by index)."""
    out = []
    for typ in ('min', 'max'):
        for k in numbered_keys(f'{m}.mep_{typ}', r'\d\d'):
            mm = from_note(k, r'^(-?[\d.]+) a\.u\.; at \(([^)]*)\) Å; nearest (.+?); (picrate|cation)(; GLOBAL)?'
                              r'(?:; kcal mol-1 value unrounded from Multiwfn surfanalysis\.txt .*)?$')
            kcal = fnum(k, 'kcal mol-1')
            au = float(mm.group(1))
            out.append(dict(key=k, type=typ, index=int(k[-2:]), kcal=kcal, au_str=mm.group(1), au=au,
                            nearest=mm.group(3), ion=mm.group(4), glob=bool(mm.group(5))))
    return out


@table('tbl_mep', 'main', 'Electrostatic potential on the 0.001 a.u. surface: extrema and surface statistics')
def t_mep():
    iso = from_value('protocol.mep_surface', r'electron-density isosurface ([\d.]+) a\.u\.').group(1)
    thr = from_note('pcm.mep_polar_area', r'\|ESP\| > (\d+) kcal mol-1').group(1)
    from_note('gas.mep_polar_area', rf'\|ESP\| > {thr} kcal mol-1')
    o1min = {}
    for m in ('pcm', 'gas'):
        ex = mep_extrema(m)
        mins = [e for e in ex if e['type'] == 'min']
        maxs = [e for e in ex if e['type'] == 'max']
        check(len(mins) == fint(f'{m}.mep_n_minima') and len(maxs) == fint(f'{m}.mep_n_maxima'),
              f'{m}: number of extremum keys = mep_n_minima / mep_n_maxima')
        gmin = [e for e in mins if e['glob']]
        gmax = [e for e in maxs if e['glob']]
        check(len(gmin) == 1 and len(gmax) == 1, f'{m}: one global minimum and one global maximum')
        close(gmin[0]['kcal'], fnum(f'{m}.mep_min_kcal', 'kcal mol-1'), 1e-9, f'{m}: global minimum = mep_min_kcal')
        close(gmax[0]['kcal'], fnum(f'{m}.mep_max_kcal', 'kcal mol-1'), 1e-9, f'{m}: global maximum = mep_max_kcal')
        check(gmin[0]['nearest'].split()[0] == fstr(f'{m}.mep_min_site'), f'{m}: global-minimum site = nearest atom')
        check(gmax[0]['nearest'].split()[0] == fstr(f'{m}.mep_max_site'), f'{m}: global-maximum site = nearest atom')
        at_o1 = [e for e in mins if e['nearest'].split()[0] == 'O1']
        check(len(at_o1) <= 1, f'{m}: at most one minimum nearest O1')
        o1min[m] = at_o1[0] if at_o1 else None
        s2p, s2n, s2t = (fnum(f'{m}.mep_sigma2_{k}', '(kcal mol-1)^2') for k in ('pos', 'neg', 'tot'))
        close(s2p + s2n, s2t, 0.01, f'{m}: sigma2_tot = sigma2_+ + sigma2_-')
        close(s2p * s2n / s2t ** 2, fnum(f'{m}.mep_nu', ''), 5e-5, f'{m}: nu = sigma2_+ sigma2_- / sigma2_tot^2')
        close(s2t * fnum(f'{m}.mep_nu', ''), fnum(f'{m}.mep_sigma2_nu', '(kcal mol-1)^2'), 0.01, f'{m}: sigma2_tot nu')
        close(fnum(f'{m}.mep_positive_area', 'Å^2') + fnum(f'{m}.mep_negative_area', 'Å^2'),
              fnum(f'{m}.mep_surface_area', 'Å^2'), 0.01, f'{m}: positive + negative area = total')
        close(100 * fnum(f'{m}.mep_polar_area', 'Å^2') / fnum(f'{m}.mep_surface_area', 'Å^2'),
              fnum(f'{m}.mep_polar_area_pct', '%'), 0.01, f'{m}: polar area fraction')
    rows = []
    def row(label, f):
        rows.append([label, f('gas'), f('pcm')])
    row('Global minimum (kcal mol⁻¹)', lambda m: f"{fmt(fnum(f'{m}.mep_min_kcal', 'kcal mol-1'), 2)} ({fstr(f'{m}.mep_min_site')})")
    row('Global maximum (kcal mol⁻¹)', lambda m: f"{fmt(fnum(f'{m}.mep_max_kcal', 'kcal mol-1'), 2)} ({fstr(f'{m}.mep_max_site')})")
    row('Minimum at the phenolate O1 (kcal mol⁻¹)', lambda m: fmt(o1min[m]['kcal'], 2) if o1min[m] else 'none')
    row('Surface minima, maxima', lambda m: f"{fint(f'{m}.mep_n_minima')}, {fint(f'{m}.mep_n_maxima')}")
    two = lambda k, u: (lambda m: fmt(fnum(f'{m}.{k}', u), 2))
    row('Surface area (Å²)', two('mep_surface_area', 'Å^2'))
    row('Positive area (Å²)', two('mep_positive_area', 'Å^2'))
    row('Negative area (Å²)', two('mep_negative_area', 'Å^2'))
    row(f'Polar area, |*V*| > {thr} kcal mol⁻¹ (%)', lambda m: fmt(fnum(f'{m}.mep_polar_area_pct', '%'), 1))
    row('Enclosed volume (Å³)', two('mep_volume', 'Å^3'))
    row('Mean positive *V* (kcal mol⁻¹)', two('mep_avg_pos', 'kcal mol-1'))
    row('Mean negative *V* (kcal mol⁻¹)', two('mep_avg_neg', 'kcal mol-1'))
    row('*σ*²~+~ (kcal² mol⁻²)', two('mep_sigma2_pos', '(kcal mol-1)^2'))
    row(f'*σ*²~{MINUS}~ (kcal² mol⁻²)', two('mep_sigma2_neg', '(kcal mol-1)^2'))
    row('*σ*²~tot~ (kcal² mol⁻²)', two('mep_sigma2_tot', '(kcal mol-1)^2'))
    row('*ν*', lambda m: fmt(fnum(f'{m}.mep_nu', ''), 4))
    row('*σ*²~tot~*ν* (kcal² mol⁻²)', two('mep_sigma2_nu', '(kcal mol-1)^2'))
    row('*Π* (kcal mol⁻¹)', two('mep_Pi', 'kcal mol-1'))
    row('MPI (kcal mol⁻¹)', two('mep_mpi_kcal', 'kcal mol-1'))
    row('MPI (eV)', lambda m: fmt(fnum(f'{m}.mep_mpi', 'eV'), 3))
    from_value('protocol.mep_note', r'PCM surface charges are not part of the mapped potential')
    tr = traj_checkpoints()
    a0, a1 = tr[0][1][5], tr[-1][1][5]
    notes = [
        '*V* is the potential of the molecular charge density alone; for IEF-PCM this is the density polarised by '
        'the continuum, and the surface charges of the continuum are not included. *σ*²~+~, *σ*²~−~ and *σ*²~tot~ '
        'are the variances of the positive, negative and total surface potential, *ν* = *σ*²~+~*σ*²~−~/(*σ*²~tot~)², '
        '*Π* the mean absolute deviation of *V* from its surface average and MPI the molecular polarity index. The '
        'relative orientation of the two ions in the IEF-PCM model is poorly defined (during the optimisation the '
        f'angle between the ring planes went from {a0}° to {a1}° while the energy fell by only '
        f"{fmt(abs(fnum('pcm.opt_energy_drop_total', 'kcal mol-1')), 2)} kcal mol⁻¹), which affects all surface "
        'statistics.'
    ]
    return dict(header=['Quantity', MODEL['gas'], MODEL['pcm']], align='lrr', rows=rows, notes=notes,
                caption=f'Electrostatic potential *V* on the *ρ* = {iso} a.u. isodensity surface of the ion pair '
                        f'(Multiwfn surface analysis), {method()}.')


# ---- NCI helpers ------------------------------------------------------------------------------
PAIR_READING = {
    frozenset({'H4N', 'O1'}): 'N4–H4N···O1 hydrogen bond',
    frozenset({'H4N', 'O7'}): 'H4N···O7 (bifurcated branch)',
    frozenset({'O1', 'O2'}): 'O1···O2 within the picrate',
    frozenset({'O1', 'O7'}): 'O1···O7 within the picrate',
    frozenset({'H11', 'O2'}): 'C11–H11···O2',
    frozenset({'H11', 'O1'}): 'C11–H11···O1',
    frozenset({'H11', 'H4N'}): 'H11···H4N',
    frozenset({'N4', 'N6'}): 'triazole ring centre',
    frozenset({'N4', 'C12'}): 'triazole ring centre',
}


def nci_spikes(m):
    out = []
    for k in numbered_keys(f'{m}.nci_spike', r'\d+'):
        mm = from_note(k, r"min RDG ([\d.]+); (\d+) points in bin; pipeline label '([^']+)'; "
                          r'lowest-RDG points lie between (.+?) \(derived')
        pairs = []
        for p in mm.group(4).split(', '):
            p = p.replace(' (inter-ion)', '')
            a, b = p.split('/')
            pairs.append((a, b))
        out.append(dict(key=k, x=fnum(k, 'a.u.'), rdg=float(mm.group(1)), rdg_str=mm.group(1), n=int(mm.group(2)),
                        label=mm.group(3), pairs=pairs))
    path = use_file(RES / CASE[m] / 'nci_spikes.csv')
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    check(len(rows) == len(out), f'{m}: nci_spikes.csv has as many rows as nci_spike keys')
    for s, r in zip(out, rows):
        close(float(r['sign_l2_rho_au']), s['x'], 1e-9, f"{m}: {s['key']} = nci_spikes.csv")
        close(float(r['min_RDG']), s['rdg'], 1e-9, f"{m}: {s['key']} min RDG = nci_spikes.csv")
        check(int(r['n_points']) == s['n'] and r['type'] == s['label'], f"{m}: {s['key']} points/label = nci_spikes.csv")
    return out


def nci_regions(m):
    out = []
    for k in numbered_keys(f'{m}.nci_region', r'\d\d'):
        mm = from_note(k, r'connected RDG<0\.5 region \((\d+) grid points\); value = sign\(l2\)rho at its lowest-RDG '
                          r'point \(RDG ([\d.]+)\); range ([+-]?[\d.]+)\.\.([+-]?[\d.]+); lowest-RDG point nearest '
                          r'(\S+) and (\S+) \((inter-ion|within one ion)\)')
        out.append(dict(key=k, num=k[-2:], x=fnum(k, 'a.u.'), n=int(mm.group(1)), rdg=float(mm.group(2)),
                        rdg_str=mm.group(2), lo=float(mm.group(3)), hi=float(mm.group(4)), lo_str=mm.group(3),
                        hi_str=mm.group(4), atoms=(mm.group(5), mm.group(6)), inter=mm.group(7) == 'inter-ion'))
    return out


def _sgn(s):
    """Signed string as recorded ('-0.0537', '+0.0591') with a true minus."""
    s = s if s.startswith(('+', '-')) else '+' + s
    return s.replace('-', MINUS)


@table('tbl_nci', 'main', 'RDG spikes, their atoms, interpretation and real-space region')
def t_nci():
    width = float(from_value('protocol.nci_spike_method', r'in ([\d.]+) bins').group(1))
    crit = from_value('protocol.nci_spike_method', r'holds >(\d+) points and its minimum RDG < ([\d.]+)')
    maxrho = from_value('protocol.nci_settings', r'RDG_maxrho = ([\d.]+) a\.u\.').group(1)
    rows, notes_a = [], []
    letters = iter('abcdefgh')
    for m in ('gas', 'pcm'):
        rows.append([f'*{MODEL[m]}*', '', '', ''])
        spikes, regions = nci_spikes(m), nci_regions(m)
        for s in spikes:
            atoms = {a for p in s['pairs'] for a in p}
            pairs = {frozenset(p) for p in s['pairs']}
            lo, hi = s['x'] - width / 2, s['x'] + width / 2
            inside = [r for r in regions if r['lo'] <= hi and r['hi'] >= lo]
            s['regions'] = [r for r in inside if frozenset(r['atoms']) in pairs] or \
                           [r for r in inside if set(r['atoms']) & atoms]
            check(len(s['regions']) >= 1, f"{m}: {s['key']} lies in a real-space region with the same atoms")
            for p in s['pairs']:
                if frozenset(p) not in PAIR_READING:
                    fail(f"{s['key']}: no reading for the atom pair {p}")
        groups = []
        for s in spikes:
            g = groups[-1] if groups else None
            if (g and g[0]['pairs'] == s['pairs'] and [r['key'] for r in g[0]['regions']] == [r['key'] for r in s['regions']]
                    and (g[0]['x'] < 0) == (s['x'] < 0)):
                g.append(s)
            else:
                groups.append([s])
        first = True
        for g in groups:
            s = g[0]
            xs = fmt(s['x'], 4, True)
            if len(g) > 1:
                L = next(letters)
                xs += f'^{L}^'
                others = g[1:]
                notes_a.append(
                    f'^{L}^ Our spike search finds {words(len(others))} further bin{"s" if len(others) > 1 else ""} of the '
                    'same region at ' + ' and '.join(fmt(o['x'], 4, True) for o in others) + ' a.u. (lowest RDG '
                    + ' and '.join(o['rdg_str'] for o in others) + ').')
            readings = []
            for p in s['pairs']:
                r = PAIR_READING[frozenset(p)]
                if r not in readings:
                    readings.append(r)
            tail = ' within the picrate'
            if len(readings) > 1 and all(r.endswith(tail) for r in readings):
                readings = [' and '.join(r[:-len(tail)] for r in readings) + tail]
            reg = '; '.join(f"{r['num']} ({_sgn(r['lo_str'])} to {_sgn(r['hi_str'])})" for r in s['regions'])
            rows.append(['; '.join(readings), xs, s['rdg_str'], reg])
            first = False
    gas_hb = [r for r in nci_regions('gas') if set(r['atoms']) == {'H4N', 'O1'}]
    check(len(gas_hb) == 1, 'one gas-phase region between H4N and O1')
    notes = [
        f"Spikes found by binning sign(*λ*~2~)*ρ* in intervals of {width:g} a.u. (bins with more than {crit.group(1)} "
        f"points and a lowest RDG below {crit.group(2)}); the value is the bin centre. The interaction is read from "
        'the two atoms nearest to the lowest-RDG points of the bin. Region: the real-space RDG < 0.5 region whose lowest-RDG point lies '
        'between the same two atoms (or, where there is none, next to one of them) and whose sign(*λ*~2~)*ρ* range, '
        'given in parentheses, contains the bin (region numbers as in Table S6). '
        f"Points with *ρ* ≥ {maxrho} a.u. were excluded, so the gas-phase N4–H4N···O1 region, which "
        f"reaches {_sgn(gas_hb[0]['lo_str'])} a.u., may be slightly truncated."
    ] + notes_a
    return dict(header=['Interaction (our reading)', 'sign(*λ*~2~)*ρ* (a.u.)', 'Min. RDG', 'Region (range, a.u.)'],
                align='lrrl', rows=rows, notes=notes, widths=[38, 14, 9, 39],
                caption='Spikes of the reduced density gradient (RDG) of the ion pair in the gas phase and in IEF-PCM '
                        f'(water), {method()}, with our reading of the interaction and the real-space region that '
                        'contains them.')


# =================================================================================================
# SUPPLEMENTARY TABLES
# =================================================================================================
@table('tblS_coords', 'SI', 'Cartesian coordinates of both optimised ion pairs')
def ts_coords():
    labels, _ = labels_and_ions()
    frames = {}
    for m in ('pcm', 'gas'):
        fr = read_geom(RES / CASE[m] / 'opt.xyz', labels)
        close(fr.energy, fnum(f'{m}.e_total', 'Ha'), 5e-9, f'{m}: energy in the opt.xyz comment = e_total')
        frames[m] = fr
    rows = []
    for i, lab in enumerate(labels):
        rows.append([lab] + [fmt(v, 6) for m in ('gas', 'pcm') for v in frames[m].xyz[i]])
    ep, eg = fmt(fnum('pcm.e_total', 'Ha'), 8), fmt(fnum('gas.e_total', 'Ha'), 8)
    hp, hg = fmt(fnum('pcm.h298', 'Ha'), 6), fmt(fnum('gas.h298', 'Ha'), 6)          # moved here from Table 4 (referee m7)
    gp, gg = fmt(fnum('pcm.g298', 'Ha'), 6), fmt(fnum('gas.g298', 'Ha'), 6)
    notes = ['Each model is given in the Cartesian frame of its own optimisation, so the two sets of coordinates are '
             'not superimposed. Picrate atoms are listed first, then the cation.']
    # figures_r2 #7: units in the column headers, so that a continuation page is self-explanatory
    return dict(header=['Atom', '*x* (Å), gas phase', '*y* (Å), gas phase', '*z* (Å), gas phase', '*x* (Å), IEF-PCM',
                        '*y* (Å), IEF-PCM', '*z* (Å), IEF-PCM'], align='lrrrrrr', rows=rows, notes=notes,
                caption=f'Cartesian coordinates of the ion pairs optimised at the {method()} level in the gas phase, '
                        f'*E* = {eg} Ha, *H* = {hg} Ha and *G* = {gg} Ha, and in IEF-PCM (water), *E* = {ep} Ha, '
                        f'*H* = {hp} Ha and *G* = {gp} Ha (298.15 K, 1 atm, harmonic). Atom labels follow the CIF.')


def mode_assignments(m, freqs):
    out = {}

    def put(i, text):
        if i in out:
            fail(f'{m}: mode {i + 1} assigned twice ({out[i]!r}, {text!r})')
        out[i] = text

    def exact(key, text):
        v = fnum(key, 'cm-1')
        idx = [i for i, f in enumerate(freqs) if abs(f - v) < 1e-6]
        check(len(idx) == 1, f'{key} matches one mode of freq.json')
        put(idx[0], text)

    def approx(key, text):
        for v in wag_list(key):
            idx = [i for i, f in enumerate(freqs) if abs(f - v) < 0.051]
            check(len(idx) == 1, f'{key}: {v} matches one mode of freq.json')
            put(idx[0], text)

    exact(f'{m}.nu_N4_H4N_stretch', '*ν* N4–H4N (hydrogen-bonded)')
    exact(f'{m}.nu_N4_H4N_inplane_bend', '*δ* N4–H4N')
    approx(f'{m}.nu_N4_H4N_oop_wag', '*γ* N4–H4N')
    exact(f'{m}.nu_N6_H6N_stretch', '*ν* N6–H6N (free)')
    exact(f'{m}.nu_N6_H6N_inplane_bend', '*δ* N6–H6N')
    approx(f'{m}.nu_N6_H6N_oop_wag', '*γ* N6–H6N')
    # C-H stretches (QA round 3, r3-m1; keys added to facts.json after round 2): picrate C3-H3 and C5-H5 (full
    # precision), the four coupled benzo stretches (key: list at 1 decimal)
    exact(f'{m}.nu_C3_H3_stretch', '*ν* C3–H3')
    exact(f'{m}.nu_C5_H5_stretch', '*ν* C5–H5')
    bz = rec(f'{m}.nu_benzo_CH_stretches')
    check(bz['unit'] == 'cm-1' and isinstance(bz['value'], list) and len(bz['value']) == 4,
          f'{m}.nu_benzo_CH_stretches: four wavenumbers in cm-1')
    for v in bz['value']:
        idx = [i for i, f in enumerate(freqs) if abs(f - float(v)) < 0.051]
        check(len(idx) == 1, f'{m}.nu_benzo_CH_stretches: {v} matches one mode of freq.json')
        put(idx[0], '*ν* benzo C–H (coupled)')
    close(max(bz['value']) - min(bz['value']), fnum(f'{m}.benzo_ch_spread', 'cm-1'), 0.06,
          f'{m}: benzo C-H spread = benzo_ch_spread')
    for k in (1, 2, 3):
        fr = fnum(f'{m}.mode{k}_rigid_fraction', '')
        f_note = float(from_note(f'{m}.mode{k}_rigid_fraction', r'normal mode \d \((-?[\d.]+) cm-1\)').group(1))
        close(f_note, freqs[k - 1], 0.006, f'{m}: rigid fraction {k} refers to mode {k}')
        put(k - 1, f'inter-ion rigid body ({pct(100 * fr)} %)')          # language_r2 #16: space before %
    return out


@table('tblS_frequencies', 'SI', 'All 93 harmonic wavenumbers of both models, with N–H, C–H and rigid-body assignments')
def ts_frequencies():
    fq = {m: json.loads(use_file(RES / CASE[m] / 'freq.json').read_text()) for m in ('pcm', 'gas')}
    f005 = json.loads(use_file(RES / 'btpc_pcm' / 'freq_step0.005.json').read_text())
    F = {m: fq[m]['frequencies_cm-1'] for m in ('pcm', 'gas')}
    F5 = f005['frequencies_cm-1']
    for m in ('pcm', 'gas'):
        check(len(F[m]) == fint(f'{m}.n_modes'), f'{m}: freq.json has n_modes wavenumbers')
        check(F[m] == sorted(F[m]), f'{m}: freq.json wavenumbers are sorted')
        close(F[m][0], fnum(f'{m}.freq_lowest', 'cm-1'), 1e-9, f'{m}: lowest mode = freq_lowest')
        close(F[m][-1], fnum(f'{m}.freq_highest', 'cm-1'), 1e-9, f'{m}: highest mode = freq_highest')
        check(fq[m]['hessian_method'] == fstr(f'{m}.hessian_method'), f'{m}: freq.json hessian_method = facts')
    check(len(F5) == len(F['pcm']) and F5 == sorted(F5), 'freq_step0.005.json: 93 sorted wavenumbers')
    close(F5[0], fnum('saddle.freq_imag_step0p005', 'cm-1'), 1e-9, 'imaginary mode with the 0.005 Bohr step')
    dmax = max(abs(a - b) for a, b in zip(F['pcm'], F5))
    close(dmax, fnum('saddle.step_change_max', 'cm-1'), 1e-9, 'largest change on halving the step = saddle.step_change_max')
    steps = {m: re.search(r'step ([\d.]+) Bohr', fq[m]['hessian_method']).group(1) for m in ('pcm', 'gas')}
    s5 = re.search(r'step ([\d.]+) Bohr', f005['hessian_method']).group(1)
    A = {m: mode_assignments(m, F[m]) for m in ('pcm', 'gas')}
    rows = []
    for i in range(len(F['pcm'])):
        rows.append([str(i + 1), wn(F['gas'][i]), A['gas'].get(i, ''), wn(F['pcm'][i]), wn(F5[i]), A['pcm'].get(i, '')])
    notes = [
        'Modes are numbered in order of increasing wavenumber within each model, so one row does not pair the same '
        'vibration of the two models. Assignments are given for the N–H modes, the C–H stretches and the three lowest '
        'modes only, from the Hessian eigenvectors: *ν* stretch, *δ* in-plane bend, *γ* out-of-plane wag; the four '
        'C–H stretches of the benzo ring (C8–H8 to C11–H11) are coupled and are assigned together; "inter-ion rigid '
        'body" gives the share of the mode that is translation and rotation of the two ions. Imaginary modes are '
        f'written with i. No infrared or Raman intensities were computed. The IEF-PCM thermochemistry uses the '
        f'{steps["pcm"]} Bohr Hessian; the {s5} Bohr column shows the effect of the step.'
    ]
    return dict(header=['Mode', f'Gas phase, {steps["gas"]} Bohr', 'Assignment (gas phase)', f'IEF-PCM, {steps["pcm"]} Bohr',
                        f'IEF-PCM, {s5} Bohr', 'Assignment (IEF-PCM)'], align='rrlrrl', rows=rows,
                notes=notes, widths=[6, 10, 26, 10, 10, 26],
                caption=f'Unscaled harmonic wavenumbers (cm⁻¹) of the ion pair, {method()}, from finite-difference '
                        'Hessians; column headings give the model and the displacement step.')


@table('tblS_crystal_contacts', 'SI', 'All X–H···N/O contacts of the P2₁/c crystal (derived from the CIF)')
def ts_crystal_contacts():
    keys = [k for k in FACTS if re.fullmatch(r'xray\.hb\d\d_[A-Z][A-Za-z0-9]*_[A-Z][A-Za-z0-9]*_[A-Z][A-Za-z0-9]*', k)]
    keys.sort(key=lambda k: int(k[7:9]))
    check([int(k[7:9]) for k in keys] == list(range(1, len(keys) + 1)), 'crystal contacts xray.hbNN numbered 1..N')
    for k in keys:
        rec(k)
    rows = []
    path = use_file(DERIVED / 'crystal_hbonds.csv')
    with open(path) as fh:
        csv_rows = list(csv.DictReader(fh))
    check(len(csv_rows) == len(keys), 'crystal_hbonds.csv has one row per xray.hbNN contact')
    intra = []
    for k, cr in zip(keys, csv_rows):
        name = fstr(k)
        d, h, a = re.fullmatch(r'(\w+)-(\w+)\.\.\.(\w+)', name).groups()
        check((cr['D'], cr['H'], cr['A']) == (d, h, a), f'{k}: same atoms as crystal_hbonds.csv')
        mm = from_note(k, r'acceptor at symmetry \d_\d{3} = \(([^)]*)\); donor at x,y,z; (intermolecular|intramolecular)[^;]*; '
                          r'(?:.*; )?(present in|not in) the DFT cluster')
        op, kind, inmodel = mm.group(1), mm.group(2), mm.group(3) == 'present in'
        check(cr['symmetry_operation_of_A'].replace(' ', '') == op.replace(' ', ''), f'{k}: symmetry = crystal_hbonds.csv')
        vals = [fnum(f'{k}_{q}', u) for q, u in (('dh', 'Å'), ('ha', 'Å'), ('da', 'Å'), ('angle', 'deg'))]
        for v, col, nd in zip(vals, ('D-H_A', 'H...A_A', 'D...A_A', 'D-H...A_deg'), (3, 3, 3, 1)):
            check(fmt(v, nd) == fmt(float(cr[col]), nd), f'{k} {col} = crystal_hbonds.csv')
        where = symop(op) + (' (same ion)' if kind == 'intramolecular' else '')
        if kind == 'intramolecular':
            intra.append(hb(name))
        rows.append([hb(name), where, fmt(vals[0], 3), fmt(vals[1], 3), fmt(vals[2], 3), fmt(vals[3], 1),
                     NDASH if kind == 'intramolecular' else ('yes' if inmodel else 'no')])
    crit = from_value('xray.hb_search_criteria', r'H\.\.\.A < ([\d.]+) Å and D-H\.\.\.A > ([\d.]+) deg, acceptors N/O, '
                                                  r'all images within \+-(\d) cells')
    ops = fstr('xray.symmetry_operations').split(';')
    sg = spacegroup(fstr('xray.space_group'))
    notes = [
        f'Derived from the deposited coordinates and the {words(len(ops))} symmetry operations of {sg} (search: H···A < '
        f'{crit.group(1)} Å, D–H···A > {crit.group(2)}°, N and O acceptors, all images within ±{crit.group(3)} cells); '
        'the donor is at *x*, *y*, *z* and the acceptor position is given. ' + ' and '.join(intra) +
        ' are contacts within one picrate with angles near 100°, not hydrogen bonds (both atoms of each belong to '
        'the modelled picrate). X-ray X–H lengths are short, so '
        '*d*(D···A) is the reliable quantity.'
    ]
    return dict(header=['D–H···A', 'Acceptor at', '*d*(D–H) (Å)', '*d*(H···A) (Å)', '*d*(D···A) (Å)', '∠D–H···A (°)',
                        'Between the modelled ions'],
                align='llrrrrl', rows=rows, notes=notes,
                caption=f'All contacts X–H···A (A = N, O) of the {sg} crystal structure of benzotriazolium picrate with '
                        f'H···A < {crit.group(1)} Å and D–H···A > {crit.group(2)}°.')


@table('tblS_bonds_angles', 'SI', 'All bond lengths and bond angles: X-ray and both models, with deviations')
def ts_bonds_angles():
    labels, ion = labels_and_ions()
    xr = read_geom(xray_xyz_path(), labels)
    G = {m: read_geom(RES / CASE[m] / 'opt.xyz', labels).xyz for m in ('pcm', 'gas')}
    X = xr.xyz
    bonds, angles = connectivity(xr.symbols, X)
    nb_all = int(from_note('pcm.bond_rmsd_all', r'(\d+) bonds incl').group(1))
    na_all = int(from_note('pcm.angle_rmsd_all', r'(\d+) bond angles').group(1))
    check(len(bonds) == nb_all and len(angles) == na_all, f'X-ray connectivity: {nb_all} bonds, {na_all} angles')
    heavy = lambda t: all(xr.symbols[i] != 'H' for i in t)
    for m in ('pcm', 'gas'):
        db = np.array([dist(G[m], i, j) - dist(X, i, j) for i, j in bonds])
        da = np.array([ang(G[m], *t) - ang(X, *t) for t in angles])
        dbh = np.array([dist(G[m], i, j) - dist(X, i, j) for i, j in bonds if heavy((i, j))])
        dah = np.array([ang(G[m], *t) - ang(X, *t) for t in angles if heavy(t)])
        rms = lambda v: float(np.sqrt(np.mean(v ** 2)))
        close(rms(db), fnum(f'{m}.bond_rmsd_all', 'Å'), 1e-9, f'{m}: bond RMSD (all) reproduced')
        close(rms(dbh), fnum(f'{m}.bond_rmsd_heavy', 'Å'), 1e-9, f'{m}: bond RMSD (heavy) reproduced')
        close(rms(da), fnum(f'{m}.angle_rmsd_all', 'deg'), 1e-9, f'{m}: angle RMSD (all) reproduced')
        close(rms(dah), fnum(f'{m}.angle_rmsd_heavy', 'deg'), 1e-9, f'{m}: angle RMSD (heavy) reproduced')
    n_key_b = n_key_a = 0
    brows = {'picrate': [], 'cation': []}
    for i, j in bonds:
        a, b = labels[i], labels[j]
        if ion[a] != ion[b]:
            fail(f'bond {a}-{b} joins the two ions')
        vals = {'xray': dist(X, i, j), 'pcm': dist(G['pcm'], i, j), 'gas': dist(G['gas'], i, j)}
        for p, q in ((a, b), (b, a)):
            if f'pcm.r_{p}_{q}' in FACTS:
                n_key_b += 1
                for m in ('xray', 'pcm', 'gas'):
                    close(vals[m], fnum(f'{m}.r_{p}_{q}', 'Å'), 1e-9, f'{m}.r_{p}_{q} = distance from the coordinates')
        brows[ion[a]].append((sorted((a, b), key=_order), [bond_name(a, b), fmt(vals['xray'], 3), fmt(vals['gas'], 3),
                              fmt(vals['pcm'], 3), fmt(vals['gas'] - vals['xray'], 3, True),
                              fmt(vals['pcm'] - vals['xray'], 3, True)]))
    arows = {'picrate': [], 'cation': []}
    for i, j, k in angles:
        a, b, c = labels[i], labels[j], labels[k]
        vals = {'xray': ang(X, i, j, k), 'pcm': ang(G['pcm'], i, j, k), 'gas': ang(G['gas'], i, j, k)}
        for p, q in ((a, c), (c, a)):
            if f'pcm.a_{p}_{b}_{q}' in FACTS:
                n_key_a += 1
                for m in ('xray', 'pcm', 'gas'):
                    close(vals[m], fnum(f'{m}.a_{p}_{b}_{q}', 'deg'), 1e-9, f'{m}.a_{p}_{b}_{q} = angle from the coordinates')
        outer = sorted((a, c), key=_order)
        arows[ion[b]].append(((_order(b), _order(outer[0]), _order(outer[1])),
                              [angle_name(a, b, c), fmt(vals['xray'], 1), fmt(vals['gas'], 1), fmt(vals['pcm'], 1),
                               fmt(vals['gas'] - vals['xray'], 1, True), fmt(vals['pcm'] - vals['xray'], 1, True)]))
    check(n_key_b > 0 and n_key_a > 0, 'some bonds and angles are also facts keys')
    STATS['S4_keyed'] = (n_key_b, n_key_a)
    rows = []
    for title, dct, key in (('Bond lengths (Å)', brows, lambda t: tuple(_order(x) for x in t[0])),
                            ('Bond angles (°)', arows, lambda t: t[0])):
        for part in ('picrate', 'cation'):
            rows.append([f'*{title}, {part}*', '', '', '', '', ''])
            rows.extend(r for _, r in sorted(dct[part], key=key))
    nbh = sum(heavy(b) for b in bonds)
    nah = sum(heavy(t) for t in angles)
    notes = [
        f'The {len(bonds)} bonds and {len(angles)} angles of the X-ray connectivity. X-ray values were computed from the '
        'deposited coordinates of the modelled ion pair (without standard uncertainties); X-ray X–H lengths are short '
        f"(see Table 2 of the main text). RMSD over the {nbh} bonds without hydrogen: "
        f"{fmt(fnum('gas.bond_rmsd_heavy', 'Å'), 4)} Å (gas phase) and {fmt(fnum('pcm.bond_rmsd_heavy', 'Å'), 4)} Å "
        f"(IEF-PCM); over the {nah} angles without hydrogen: {fmt(fnum('gas.angle_rmsd_heavy', 'deg'), 1)}° and "
        f"{fmt(fnum('pcm.angle_rmsd_heavy', 'deg'), 1)}°. Table 2 of the main text quotes the values listed in the CIF, with their "
        'standard uncertainties; values recomputed here from the deposited coordinates can differ from those in the '
        'last digit.'                                                   # QA round 1, m19
    ]
    return dict(header=['Parameter', MODEL['xray'], MODEL['gas'], MODEL['pcm'], 'Gas phase − X-ray', 'IEF-PCM − X-ray'],
                align='lrrrrr', rows=rows, notes=notes,
                caption=f'All bond lengths (Å) and bond angles (°) of the ion pair in the crystal '
                        f"({fnum('xray.temperature', 'K'):.0f} K) and after {method()} optimisation in the gas "
                        'phase and in IEF-PCM (water); the last two columns give each model minus X-ray.')


@table('tblS_mep_extrema', 'SI', 'Every surface extremum of the electrostatic potential (Multiwfn indices)')
def ts_mep_extrema():
    iso = from_value('protocol.mep_surface', r'electron-density isosurface ([\d.]+) a\.u\.').group(1)
    ha_kcal = float(from_note('pcm.zpe_kcal', r'x ([\d.]+)').group(1))
    rows = []
    for m in ('gas', 'pcm'):
        ex = mep_extrema(m)
        path = use_file(RES / CASE[m] / 'mep_extrema.csv')
        with open(path) as fh:
            csvr = {(r['type'][:3], int(r['index'])): r for r in csv.DictReader(fh)}
        check(len(csvr) == len(ex), f'{m}: mep_extrema.csv has one row per extremum key')
        surf, kind = {}, None
        for ln in use_file(RAW / CASE[m] / 'surfanalysis.txt').read_text().splitlines():
            if 'surface minima' in ln: kind = 'min'; continue
            if 'surface maxima' in ln: kind = 'max'; continue
            tk = ln.replace('*', ' ').split()
            if kind and len(tk) == 7 and tk[0].isdigit():
                surf[(kind, int(tk[0]))] = float(tk[3])
        check(len(surf) == len(ex), f'{m}: surfanalysis.txt lists every extremum')
        first = True
        for e in ex:
            r = csvr[(e['type'], e['index'])]
            # QA r4-m1 / figures_r2 #2: facts hold Multiwfn's unrounded value (surfanalysis.txt); the CSV rounds to 3 decimals
            close(float(r['ESP_kcal_mol']), e['kcal'], 6e-4, f"{e['key']} = mep_extrema.csv to its 3 decimals")
            close(surf[(e['type'], e['index'])], e['kcal'], 1e-9, f"{e['key']} = surfanalysis.txt (unrounded)")
            close(e['au'] * ha_kcal, e['kcal'], 0.002, f"{e['key']}: a.u. x 627.51 = kcal/mol")
            near = re.sub(r' A\)', ')', e['nearest'])
            typ = ('minimum' if e['type'] == 'min' else 'maximum') + (' (global)' if e['glob'] else '')
            rows.append([MODEL[m] if first else '', str(e['index']), typ, fmt(e['kcal'], 2), e['au_str'].replace('-', MINUS),
                         near, e['ion']])
            first = False
    notes = [
        'No. is the index of the extremum in the Multiwfn output, used to mark its position on the surface. Minima '
        'and maxima are topological extrema of *V* on the surface, so a maximum can be negative (between two '
        'negative regions) and a minimum positive. Nearest atoms are given with their distance (Å) from the '
        'extremum. *V* is the potential of the molecular charge density alone (for IEF-PCM, polarised by the '
        'continuum, without its surface charges).'
    ]
    return dict(header=['Model', 'No.', 'Type', '*V* (kcal mol⁻¹)', '*V* (a.u.)', 'Nearest atoms (Å)', 'Ion'],
                align='lrlrrll', rows=rows, notes=notes,
                caption=f'Every extremum of the electrostatic potential *V* on the *ρ* = {iso} a.u. isodensity surface of '
                        f'the ion pair located by Multiwfn, gas phase and IEF-PCM (water), {method()}.')


@table('tblS_nci_regions', 'SI', 'Real-space RDG < 0.5 regions of both models')
def ts_nci_regions():
    labels, ion = labels_and_ions()
    rows = []
    for m in ('gas', 'pcm'):
        first = True
        for r in nci_regions(m):
            a, b = r['atoms']
            ions = {ion[a], ion[b]}
            check((len(ions) == 2) == r['inter'], f"{r['key']}: inter-ion flag agrees with the atom labels")
            where = 'between the ions' if r['inter'] else f'within the {ions.pop()}'
            rows.append([MODEL[m] if first else '', r['num'], fmt(r['x'], 4, True), r['rdg_str'],
                         f"{_sgn(r['lo_str'])} to {_sgn(r['hi_str'])}", str(r['n']), f'{a}, {b} ({where})'])
            first = False
    maxrho = from_value('protocol.nci_settings', r'RDG_maxrho = ([\d.]+) a\.u\.').group(1)
    notes = [
        'Regions are connected sets of grid points of the RDG cube with RDG < 0.5; sign(*λ*~2~)*ρ* is given at the '
        'lowest-RDG point of each region and as its range over the region. Atoms: the two atoms nearest to the '
        f'lowest-RDG point. Points with *ρ* ≥ {maxrho} a.u. were excluded from the analysis.'
    ]
    return dict(header=['Model', 'Region', 'sign(*λ*~2~)*ρ* (a.u.)', 'Lowest RDG', 'Range (a.u.)', 'Grid points',
                        'Nearest atoms'], align='lrrrlrl', rows=rows, notes=notes, widths=[12, 9, 10, 8, 20, 7, 28],
                caption=f'Real-space regions of low reduced density gradient (RDG < 0.5) of the ion pair in the gas '
                        f'phase and in IEF-PCM (water), {method()}.')


def ladder(m):
    s = fstr(f'{m}.orbitals_near_gap')
    items = re.findall(r'(HOMO|LUMO)([+-]\d+)? (-?[\d.]+) eV \((-?[\d.]+)% picrate\)', s)
    check(len(items) == 16, f'{m}: orbitals_near_gap lists HOMO-7 to LUMO+7')
    out = [(f'{a}{b}', float(e), e, float(p), p) for a, b, e, p in items]
    path = use_file(DERIVED / f'mo_fragment_composition_{CASE[m]}.csv')
    with open(path) as fh:
        csvr = {r['label']: r for r in csv.DictReader(fh)}
    for lab, e, _, p, _ in out:
        r = csvr.get(lab) or fail(f'{lab} missing in {path.name}')
        close(float(r['energy_eV']), e, 0.0005, f'{m} {lab}: energy = derived CSV')
        close(float(r['picrate_pct']), p, 0.055, f'{m} {lab}: picrate share = derived CSV')
    close(dict((o[0], o[1]) for o in out)['HOMO'], fnum(f'{m}.homo_ev', 'eV'), 0.0005, f'{m}: HOMO in the ladder = homo_ev')
    close(dict((o[0], o[1]) for o in out)['LUMO'], fnum(f'{m}.lumo_ev', 'eV'), 0.0005, f'{m}: LUMO in the ladder = lumo_ev')
    return out


@table('tblS_mo_ladder', 'SI', 'Orbital energies and picrate share, HOMO−7 to LUMO+7')
def ts_mo_ladder():
    L = {m: ladder(m) for m in ('pcm', 'gas')}
    check([o[0] for o in L['pcm']] == [o[0] for o in L['gas']], 'same orbital labels in both models')
    rows = []
    for k in range(len(L['pcm']) - 1, -1, -1):
        lab = L['pcm'][k][0].replace('-', MINUS)
        rows.append([lab] + [x for m in ('gas', 'pcm') for x in (L[m][k][2].replace('-', MINUS),
                                                                  L[m][k][4].replace('-', MINUS))])
    notes = ['The cation share is 100 % minus the picrate share. Shares come from a Mulliken partition of each orbital over '
             'the basis functions of the two ions; for diffuse unoccupied orbitals they can fall slightly outside 0–100 %.']
    return dict(header=['Orbital', '*E* (eV), gas phase', 'Picrate (%), gas phase', '*E* (eV), IEF-PCM',
                        'Picrate (%), IEF-PCM'], align='lrrrr', rows=rows, notes=notes,
                caption=f'Kohn–Sham orbital energies and Mulliken share of the picrate from HOMO{MINUS}7 to LUMO+7 of the '
                        f'ion pair in the gas phase and in IEF-PCM (water), {method()}.')


def traj_checkpoints():
    """saddle.traj_stepNN -> [(step, (E, N4-H4N, H4N...O1, N4...O1, angle, inter-ring, H11...O2) as strings)]."""
    keys = [k for k in FACTS if re.fullmatch(r'saddle\.traj_step\d\d', k)]
    keys.sort(key=lambda k: int(k[-2:]))
    if not keys:
        fail('no saddle.traj_stepNN keys')
    pat = (r'^E (-?[\d.]+) Ha; N4-H4N ([\d.]+); H4N\.\.\.O1 ([\d.]+); N4\.\.\.O1 ([\d.]+) Å; angle ([\d.]+) deg; '
           r'interplanar ([\d.]+) deg; H11\.\.\.O2 ([\d.]+) Å$')
    data = [(int(k[-2:]), from_value(k, pat).groups()) for k in keys]
    check(data[0][0] == 0 and data[-1][0] == fint('pcm.opt_steps'), 'trajectory checkpoints run from step 0 to the last step')
    close(float(data[-1][1][5]), fnum('pcm.interplanar_picrate_cation', 'deg'), 0.05, 'final inter-ring angle = facts')
    return data


@table('tblS_pcm_trajectory', 'SI', 'IEF-PCM optimisation: checkpoints and the imaginary-mode checks')
def ts_pcm_trajectory():
    ha_kcal = float(from_note('pcm.zpe_kcal', r'x ([\d.]+)').group(1))
    data = traj_checkpoints()
    nsteps = fint('pcm.opt_steps')
    Efin = float(data[-1][1][0])
    close(Efin, fnum('pcm.e_total', 'Ha'), 5e-9, 'final checkpoint energy = pcm.e_total')
    close((Efin - float(data[0][1][0])) * ha_kcal, fnum('pcm.opt_energy_drop_total', 'kcal mol-1'), 1e-4,
          'E(final) - E(step 0) = pcm.opt_energy_drop_total')
    rows = []
    for step, g in data:
        E = float(g[0])
        rows.append([str(step), fmt((E - Efin) * ha_kcal, 3), g[1], g[2], g[3], g[4], g[5], g[6]])
    cr = from_value('protocol.opt_criteria', r'\|dE\| < ([\d.e+-]+) Ha, RMS grad < ([\d.e+-]+), max grad < ([\d.e+-]+) '
                                             r'Ha/Bohr, RMS disp < ([\d.e+-]+), max disp < ([\d.e+-]+) Å')
    check(fstr('pcm.opt_converged').startswith('no'), 'IEF-PCM optimisation did not meet all criteria')
    g_rms, g_max = fnum('pcm.opt_final_grad_rms', 'Ha/Bohr'), fnum('pcm.opt_final_grad_max', 'Ha/Bohr')
    check(g_rms < float(cr.group(2)) and g_max < float(cr.group(3)), 'final gradients inside the thresholds')
    dE = fnum('pcm.opt_last_step_dE', 'Ha')
    d_rms, d_max = fnum('pcm.opt_last_step_disp_rms', 'Å'), fnum('pcm.opt_last_step_disp_max', 'Å')
    check(abs(dE) > float(cr.group(1)) and d_rms > float(cr.group(4)) and d_max > float(cr.group(5)),
          'last-step energy change and displacements outside the thresholds')
    cation_w = 1 - fnum('saddle.imag_mode_picrate_weight', '')
    from_note('saddle.imag_mode_picrate_weight', r'\(cation 0\.66\)')
    drift = []
    for k in ('saddle.drift_0_to_70_overlap_imag', 'saddle.drift_40_to_70_overlap_imag'):
        a, b = re.fullmatch(r'saddle\.drift_(\d+)_to_(\d+)_overlap_imag', k).groups()
        check(int(b) == nsteps, f'{k} ends at the final step')
        drift.append((a, b, 100 * fnum(k, '')))
    close(fnum('saddle.freq_imag_step0p01', 'cm-1'), fnum('pcm.freq_lowest', 'cm-1'), 1e-9, 'imaginary mode = pcm.freq_lowest')
    st01 = from_value('pcm.hessian_method', r'step ([\d.]+) Bohr').group(1)
    fa = json.loads(use_file(RES / 'btpc_pcm' / 'freq.json').read_text())
    fb = json.loads(use_file(RES / 'btpc_pcm' / 'freq_step0.005.json').read_text())
    check(fa['hessian_method'] == fstr('pcm.hessian_method'), 'freq.json is the 0.01 Bohr Hessian')
    close((fa['ZPE_Ha'] - fb['ZPE_Ha']) * ha_kcal, fnum('saddle.step_change_zpe', 'kcal mol-1'), 1e-9,
          'step_change_zpe = ZPE(0.01 Bohr) - ZPE(0.005 Bohr)')
    close((fa['G298_Ha'] - fb['G298_Ha']) * ha_kcal, fnum('saddle.step_change_g', 'kcal mol-1'), 1e-9,
          'step_change_g = G(0.01 Bohr) - G(0.005 Bohr)')
    notes = [
        # figures_r2 #6: every header carries its unit, so the old units sentence (wrong for Step and ΔE) is gone
        f"Energies relative to the final structure (step {nsteps}, *E* = {fmt(Efin, 8)} Ha). The inter-ring angle is the "
        'angle between the least-squares planes of the picrate ring and of the cation. The optimisation '
        f"was accepted on the gradient criteria (final RMS and maximum gradient {sci(g_rms, 2)} and {sci(g_max, 2)} "
        f"Ha Bohr⁻¹, thresholds {sci_str(cr.group(2))} and {sci_str(cr.group(3))}), whereas the energy change of the "
        f"last step ({sci(dE)} Ha, threshold {sci_str(cr.group(1))}) and its RMS and maximum displacements "
        f"({fmt(d_rms, 3)} and {fmt(d_max, 3)} Å, thresholds {fmt(float(cr.group(4)), 4)} and "
        f"{fmt(float(cr.group(5)), 4)} Å) were not met.",
        'The Hessian at the final geometry has one imaginary mode, '
        f"{wn(fnum('saddle.freq_imag_step0p01', 'cm-1'))} cm⁻¹ with a displacement step of {st01} Bohr and "
        f"{wn(fnum('saddle.freq_imag_step0p005', 'cm-1'))} cm⁻¹ with {step_005()} Bohr. No other wavenumber differs "
        f"between the two Hessians by more than {fmt(fnum('saddle.step_change_max', 'cm-1'), 2)} cm⁻¹, and the "
        f"zero-point energy and *G* from the {st01} Bohr Hessian differ from the {step_005()} Bohr values by "
        f"{fmt(fnum('saddle.step_change_zpe', 'kcal mol-1'), 3)} and {fmt(fnum('saddle.step_change_g', 'kcal mol-1'), 2)} "
        # referee_r2 M1: both wavenumbers lie below the resolution of the Hessians (curvature of the overall rotations)
        "kcal mol⁻¹. Both values of the imaginary mode lie below the numerical resolution of the Hessians, whose "
        "curvature along the overall rotations of the pair, zero for an exact Hessian at a stationary point, reaches "
        f"{wn(-fnum('pcm.hessian_rot_curvature_max_abs', 'cm-1'))} cm⁻¹ (Section 2.4 of the main text). The mode is "
        f"{pct(100 * fnum('saddle.imag_mode_rigid_fraction', ''))} % rigid-body motion of the two ions, with "
        f"{100 * cation_w:.0f} % of its mass-weighted amplitude on the cation. The displacement from step "
        f"{drift[0][0]} to step {drift[0][1]} projects {drift[0][2]:.0f} % onto this mode; for steps {drift[1][0]} to "
        f"{drift[1][1]} the projection is {drift[1][2]:.0f} %."
    ]
    return dict(header=['Step', 'Δ*E* (kcal mol⁻¹)', 'N4–H4N (Å)', 'H4N···O1 (Å)', 'N4···O1 (Å)', '∠N4–H4N···O1 (°)',
                        'Inter-ring angle (°)', 'H11···O2 (Å)'], align='rrrrrrrr', rows=rows, notes=notes,
                caption=f'Selected steps of the geometry optimisation of the ion pair in IEF-PCM (water), {method()}.')


@table('tblS_protocol', 'SI', 'Computational settings and software versions')
def ts_protocol():
    rows = []
    add = lambda item, setting: rows.append([item, setting])
    fn = fstr('protocol.functional')
    xc = _grab(fn, r'libxc (HYB_GGA_XC_B3LYP), ID (\d+), (VWN-RPA) local correlation, same definition as Gaussian',
               'protocol.functional', 'value')
    add('Functional', f'B3LYP (libxc {xc.group(1)}, ID {xc.group(2)}; {xc.group(3)} local correlation, as in Gaussian)')
    bpe = from_value('model.basis_per_element', r'^C, N, O: (\d+) functions each \((\w+), spherical d\); H: (\d+) functions \((\w+)\)$')
    add('Basis set', f"{fstr('protocol.basis')}; {fint('model.n_basis')} basis functions: {bpe.group(1)} for each C, N "
                     f"and O atom ({bpe.group(2)}), {bpe.group(3)} for each H atom ({bpe.group(4)})")
    check('RI-JK' in fstr('protocol.density_fitting'), 'density fitting is RI-JK')
    check('not logged' in fstr('protocol.aux_basis'), 'auxiliary basis is the (unlogged) PySCF default')
    add('Density fitting', 'RI-JK; auxiliary basis: the PySCF default (not recorded in the logs)')
    add('Integration grid', f"level {fint('protocol.grid_level')}: {tidy(fstr('protocol.grid_points'))}")
    add('SCF convergence', f"{sci(fnum('protocol.scf_conv_tol_opt_gas', 'Ha'), 0)} Ha (gas-phase optimisation), "
                           f"{sci(fnum('protocol.scf_conv_tol_opt_pcm', 'Ha'), 0)} Ha (IEF-PCM optimisation), "
                           f"{sci(fnum('protocol.scf_conv_tol_final', 'Ha'), 0)} Ha (final and Hessian calculations); "
                           f"at most {fint('protocol.scf_max_cycle')} cycles")
    npts = from_note('protocol.pcm_lebedev_order', r'^(\d+) points per atomic sphere').group(1)
    add('Continuum model', f"{fstr('protocol.solvation_model')}, {fstr('protocol.solvent')} "
                           f"(*ε* = {fmt(fnum('protocol.solvent_eps', ''), 4)}); {tidy(fstr('protocol.pcm_discretisation'))}; "
                           f"Lebedev order {fint('protocol.pcm_lebedev_order')} ({npts} points per atomic sphere)")
    add('Cavity radii', tidy(fstr('protocol.pcm_radii')))
    opt = fstr('protocol.optimizer')
    check(' via pyscf.geomopt.geometric_solver' in opt and 'maxiter' in opt, 'optimizer text has the expected parts')
    opt = opt.replace(' via pyscf.geomopt.geometric_solver', '').replace('(= 3N)', '(= 3*N*)')
    opt = re.sub(r'maxiter (\d+)', r'at most \1 steps', opt)
    add('Geometry optimiser', tidy(opt))
    cr = from_value('protocol.opt_criteria', r'\|dE\| < ([\d.e+-]+) Ha, RMS grad < ([\d.e+-]+), max grad < ([\d.e+-]+) '
                                             r'Ha/Bohr, RMS disp < ([\d.e+-]+), max disp < ([\d.e+-]+) Å \(all five required\)')
    check(fstr('gas.opt_converged').startswith('yes') and fstr('pcm.opt_converged').startswith('no'),
          'gas optimisation met all criteria, IEF-PCM did not')
    add('Convergence criteria', f'|Δ*E*| < {sci_str(cr.group(1))} Ha; RMS and maximum gradient < {sci_str(cr.group(2))} '
                                f'and {sci_str(cr.group(3))} Ha Bohr⁻¹; RMS and maximum displacement < '
                                f'{sci_str(cr.group(4))} and {sci_str(cr.group(5))} Å (geomeTRIC\'s five criteria with '
                                'Gaussian\'s default thresholds, but with displacements in Å, where Gaussian uses Bohr and '
                                'has no energy criterion); all five met in the gas phase, the gradient criteria only in IEF-PCM')
    pre = from_value('protocol.preopt', r'^(B3LYP/6-31G\(d\)) with loose criteria')
    from_note('protocol.preopt', r'earlier grid-3 runs')
    frame = from_value('pcm.opt_start', r'grid-3 IEF-PCM optimisation trajectory .*\(frame (\d+)\)').group(1)
    from_value('gas.opt_start', r'converged grid-3 geometry')
    add('Starting geometries', f'{pre.group(1)} pre-optimisation with loose criteria (pipeline default), done in earlier '
                               'runs with integration grid level 3; the final gas-phase optimisation started from the '
                               f'converged structure of its earlier run, the IEF-PCM one from an intermediate structure '
                               f'(step {frame}) of its earlier run')
    from_note('gas.opt_steps', r'number of energy\+gradient evaluations')
    add('Optimisation length', f"gas phase {fint('gas.opt_steps')} energy and gradient evaluations; IEF-PCM "
                               f"{fint('pcm.opt_steps')} steps")
    hm = tidy(fstr('protocol.hessian_method').replace('(6N =', '(6*N* ='))
    st = {m: from_value(f'{m}.hessian_method', r'step ([\d.]+) Bohr').group(1) for m in ('pcm', 'gas')}
    from_note('saddle.freq_imag_step0p005', r'first successful Hessian')
    check(fstr('protocol.fd_gradient_grid_response').startswith('no'), 'FD gradients without grid response')
    add('Hessian', f"{hm}; analytic gradients without the grid-weight response (PySCF default); step {st['pcm']} Bohr "
                   f"for IEF-PCM (an earlier IEF-PCM Hessian with {step_005()} Bohr serves as a check) and {st['gas']} Bohr "
                   'for the gas phase')
    pr = fstr('protocol.hessian_projection')
    check(', own routine dft/vibrations.py' in pr, 'projection text has the expected parts')
    add('Normal modes', tidy(pr.replace(', own routine dft/vibrations.py', ', own code')))
    from_note('protocol.frequency_scaling', r'No IR/Raman intensities were computed')
    add('Frequency scaling', tidy(fstr('protocol.frequency_scaling')) + '; no infrared or Raman intensities')
    from_note('protocol.thermo_p', r'rigid-rotor/harmonic-oscillator thermochemistry, no low-frequency correction')
    add('Thermochemistry', f"{fnum('protocol.thermo_T', 'K'):.2f} K, {fint('protocol.thermo_p', 'Pa')} Pa; rigid rotor, "
                           'harmonic oscillator; in addition a quasi-RRHO entropy (free-rotor interpolation, '
                           '*ν*~0~ = 100 cm⁻¹, *B*~av~ = 10⁻⁴⁴ kg m²) computed from the same wavenumbers')
    add('Descriptors', 'IP = −*E*~HOMO~, EA = −*E*~LUMO~, *χ* = (IP + EA)/2, *μ* = −*χ*, *η* = (IP − EA)/2, '
                       '*S* = 1/(2*η*), *ω* = *μ*²/(2*η*)')
    from_value('protocol.fragment_composition_method', r'^Mulliken partition of each MO over the two fragments')
    add('Fragment shares', 'Mulliken partition of each molecular orbital over the basis functions of the two ions')
    add('Atomic and fragment charges', 'Mulliken and Löwdin populations, with the overlap matrix rebuilt from the '
        'basis set in the wavefunction file')
    add('Density topology', 'bond critical point between H4N and O1 located by a Newton search on the gradient of *ρ*; '
        '*ρ*, its Laplacian, *G*, *V* = ∇²*ρ*/4 − 2*G* and ELF evaluated from the Kohn–Sham orbitals of the '
        'wavefunction file, orbital derivatives by central differences')
    # figures_r2 #8: no undefined ESP / MEP abbreviations
    cub = fstr('protocol.cubes').replace(' (PySCF cubegen)', ' (PySCF)')
    check('ESP 64^3' in cub and 'MEP on 0.001 a.u. density' in cub, 'cube text has the expected parts')
    cub = cub.replace('ESP 64^3', 'electrostatic potential 64^3').replace('MEP on 0.001 a.u. density',
                                                                         'electrostatic potential mapped on the 0.001 a.u. density')
    add('Cube files', tidy(cub))
    ms = fstr('protocol.mep_surface')
    check(ms.startswith('Multiwfn main function 12: '), 'MEP surface text has the expected parts')
    esp = {m: from_value(f'protocol.mep_esp_code_{m}', r'\((.+)\)$').group(1) for m in ('pcm', 'gas')}
    add('Electrostatic potential', tidy(ms.replace('Multiwfn main function 12: ', 'Multiwfn: ')
                                        .replace('ESP evaluated', 'potential evaluated')) +
        f"; code {esp['pcm']} (IEF-PCM), {esp['gas'].replace('Multiwfn original ESP code', 'the original Multiwfn routine')} (gas phase)")
    nci = fstr('protocol.nci_settings')
    check(nci.startswith('Multiwfn 20 -> 1 (NCI), '), 'NCI text has the expected parts')
    nci = nci.replace('Multiwfn 20 -> 1 (NCI), ', 'Multiwfn NCI analysis, ').replace('rho >=', '*ρ* ≥')
    add('NCI analysis', tidy(nci))
    add('NCI grids', f"IEF-PCM {tidy(fstr('protocol.nci_grid_pcm'))}; gas phase {tidy(fstr('protocol.nci_grid_gas'))}")
    sm = fstr('protocol.nci_spike_method')
    mm = _grab(sm, r'^dft/mwfn\.py:nci_post: (.*?); type label:', 'protocol.nci_spike_method', 'value')
    add('RDG spikes', tidy(mm.group(1).replace('sign(l2)rho', 'sign(*λ*~2~)*ρ*')))
    el = from_value('protocol.elf_lol_settings', r'plane through three atoms, (\d+) x (\d+) points, aug2D = ([\d.]+) Bohr, '
                                                 r'Becke definition')
    pp = from_value('protocol.elf_plane_picrate', r'through (\w+), (\w+), (\w+)')
    pc = from_value('protocol.elf_plane_cation', r'through (\w+), (\w+), (\w+)')
    add('ELF and LOL', f'Multiwfn, Becke definition; planes through {", ".join(pp.groups())} (picrate ring) and '
                       f'{", ".join(pc.groups())} (benzotriazolium); {el.group(1)} × {el.group(2)} points, '
                       f'{el.group(3)} Bohr extension beyond the atoms')
    sw = [('PySCF', 'software.pyscf', r'^([\d.]+)$'), ('geomeTRIC', 'software.geometric', r'^([\d.]+)$'),
          ('Multiwfn', 'software.multiwfn', r'^([\d.]+) '), ('Python', 'software.python', r'^([\d.]+)$')]
    parts = [f'{name} {from_value(k, p).group(1)}' for name, k, p in sw]
    parts.append(fstr('software.numpy_scipy_matplotlib'))
    add('Software (calculations)', '; '.join(parts) + f"; libxc {fstr('software.libxc')}")
    add('Software (analysis and figures)', fstr('software.figures_env') + '; our own code for the fragment shares, '
        'the Löwdin charges, the density topology and the quasi-RRHO entropy')
    # QA round 1, m20: for the IEF-PCM run the grid level and the cavity Lebedev order rest on the launch-script
    # defaults only (facts notes, open question O2); say so instead of presenting them as read from its logs
    from_note('protocol.grid_level', r'btpc_pcm was NOT in the verify run')
    from_note('protocol.grid_level', r'matched grid 4 for btpc_gas \(stored E reproduced')
    from_note('protocol.pcm_lebedev_order', r'btpc_pcm itself was not verified')
    notes = ['Settings read from the input files, logs and program sources of the calculations; the libxc version and '
             'the auxiliary basis were not recorded. The integration grid level of the IEF-PCM calculations and the '
             'Lebedev order of the cavity are the defaults of our job script; they are not recorded in the logs of '
             'those calculations (the gas-phase grid level was confirmed by re-evaluating the stored energy).']
    return dict(header=['Item', 'Setting'], align='ll', rows=rows, notes=notes,
                caption=f'Computational settings and software used for the {method()} calculations on benzotriazolium '
                        'picrate.')


@table('tblS_elf_lol', 'SI', 'ELF and LOL plane maps: ranges, ring-centroid values, O1 offset')
def ts_elf_lol():
    pp = from_value('protocol.elf_plane_picrate', r'through (\w+), (\w+), (\w+)').groups()
    pc = from_value('protocol.elf_plane_cation', r'through (\w+), (\w+), (\w+)').groups()
    el = from_value('protocol.elf_lol_settings', r'(\d+) x (\d+) points')
    rows = []
    def row(label, f):
        rows.append([label, f('gas'), f('pcm')])
    rng = lambda k: (lambda m: fstr(f'{m}.{k}').replace('..', NDASH))
    for fn in ('elf', 'lol'):
        F = fn.upper()
        row(f'{F} range, picrate plane', rng(f'{fn}_picrate_ring_range'))
        row(f'{F} range, cation plane', rng(f'{fn}_benzotriazolium_ring_range'))
        from_note(f'pcm.{fn}_triazole_centroid', r'interpolated at the ring centroid')
        for ring, lab in (('picrate', 'picrate ring'), ('benzo', 'benzo ring'), ('triazole', 'triazole ring')):
            row(f'{F} at the {lab} centroid', lambda m, fn=fn, ring=ring: fmt(fnum(f'{m}.{fn}_{ring}_centroid', ''), 3))
    from_note('pcm.cation_plane_offset_O1', r'H4N lies within 0\.01 Å of it')
    row('Distance of O1 from the cation plane (Å)', lambda m: fmt(fnum(f'{m}.cation_plane_offset_O1', 'Å'), 3))
    for m in ('pcm', 'gas'):
        for k in ('elf_picrate_ring_range', 'elf_benzotriazolium_ring_range', 'lol_picrate_ring_range',
                  'lol_benzotriazolium_ring_range'):
            check(re.fullmatch(r'\d\.\d{4}\.\.\d\.\d{4}', fstr(f'{m}.{k}')) is not None, f'{m}.{k} is a range string')
    notes = [f'Maps of {el.group(1)} × {el.group(2)} points in the plane through {", ".join(pp)} (picrate ring) and in the '
             f'plane through {", ".join(pc)} (benzotriazolium, which contains the whole cation). Centroid values '
             'interpolated from the maps. The O1 distance is signed with respect to the plane normal; H4N lies within '
             '0.01 Å of the cation plane.']
    return dict(header=['Quantity', MODEL['gas'], MODEL['pcm']], align='lrr', rows=rows, notes=notes,
                caption=f'ELF and LOL values sampled in the planes of the two rings of the ion pair, gas phase and '
                        f'IEF-PCM (water), {method()}.')


# =================================================================================================
# INDEX.md
# =================================================================================================
def compress_keys(keys):
    groups, order = {}, []
    for k in keys:
        p, s = k.split('.', 1)
        if s not in groups:
            groups[s] = []
            order.append(s)
        if p not in groups[s]:
            groups[s].append(p)
    items = [(tuple(groups[s]), s) for s in order]
    pfx = lambda ps: ps[0] if len(ps) == 1 else '{' + ','.join(ps) + '}'
    out, i = [], 0
    while i < len(items):
        ps, s = items[i]
        m = re.fullmatch(r'(.*?)(\d+)', s)
        if m:
            stem, n0, w = m.group(1), int(m.group(2)), len(m.group(2))
            j, last = i, n0
            while j + 1 < len(items):
                ps2, s2 = items[j + 1]
                m2 = re.fullmatch(r'(.*?)(\d+)', s2)
                if ps2 == ps and m2 and m2.group(1) == stem and int(m2.group(2)) == last + 1 and len(m2.group(2)) == w:
                    j, last = j + 1, last + 1
                else:
                    break
            if j > i:
                out.append(f'`{pfx(ps)}.{stem}{n0:0{w}d}`…`{stem}{last:0{w}d}`')
                i = j + 1
                continue
        out.append(f'`{pfx(ps)}.{s}`')
        i += 1
    return out


DEVIATIONS = [
    'Location and naming follow the task, not PLAN.md §2: script `papers/btpc/tables/make_tables.py`, outputs '
    '`papers/btpc/tables/<stem>.md`, label `{#tbl:<stem>}` (so `@tbl:tbl_geometry`, `@tbl:tblS_coords`), caption line '
    'after the table; since the round-1 revision the notes and footnotes are sentences of the caption, so that they '
    'are not set as line-numbered body text (figures_r1 #1).',
    'Precision of agreement statistics: bond RMSD/MAE 4 decimals and angle RMSD 2 decimals (as in the figure legends), '
    'because 3 and 1 decimals would print the two models as identical. Everything else follows STYLE.md.',
    'Table 3: the proton asymmetry d(H···O1) − d(N4–H4N) is in the footnote, not an extra row (it does not fit the '
    'D–H / H···A / D···A columns); all three values come from `*.nh_minus_ho_asym` (X-ray key added to facts.json '
    'after QA round 2). Crystal-only contacts N6–H6N···O6 and C9–H9···O5 '
    'are single "X-ray only" rows instead of rows with "absent" model cells.',
    'Table 4: an electronic-energy row was added (so that H − E and G − E can be read); the three lowest wavenumbers '
    'carry 1 decimal (3.2i), the N–H wavenumbers and red shifts are integers, as quoted in the text (figures_r1 #15). '
    'Thermal terms in kcal mol⁻¹ with the quasi-RRHO entropy beside the harmonic one; the absolute *H* and *G* moved '
    'to the caption of Table S1 (referee m7). The imaginary-mode check with the 0.005 Bohr step and the PySCF count '
    'are in Section 2.4 of the paper and in Table S8; footnote a points to them (round 2, figures_r2 #6).',
    'Table 6: a derived row "Character of the HOMO–LUMO gap" (within the picrate / picrate → cation) was added, from '
    'the picrate shares of HOMO and LUMO.',
    'Table 5: no column with the type labels of the spike search (they are thresholds on sign(λ₂)ρ, PLAN footnote); gas '
    'spikes 1–3 (same atoms, same region, same sign) are merged into one row with footnote a; the matching real-space '
    'region is found automatically: its lowest-RDG point lies between the same atom pair as the spike (or, if no region '
    'does, next to one of its atoms) and its sign(λ₂)ρ range overlaps the bin.',
    'Table S1 (`tblS_coords`): one table with the IEF-PCM and gas-phase coordinates side by side instead of parts (a) '
    'and (b); no element column (the CIF label starts with the element symbol).',
    'Table S2 (`tblS_frequencies`): columns Mode | gas phase | assignment (gas phase) | IEF-PCM 0.01 Bohr | IEF-PCM '
    '0.005 Bohr | assignment (IEF-PCM) (gas phase first since the round-1 revision): two assignment columns because '
    'mode k is a different vibration in each model. The C–H stretch assignments come from `*.nu_C3_H3_stretch`, '
    '`*.nu_C5_H5_stretch` and `*.nu_benzo_CH_stretches` (the four coupled benzo stretches are one assignment, "benzo '
    'C–H", on four rows; QA round 3, r3-m1).',
    'Table S4 (`tblS_bonds_angles`): every value is computed from the coordinate files with the connectivity rule of '
    'build_facts.py; instead of flagging derived cells, the script checks every cell that also exists as a facts key '
    '({S4_bonds} bonds, {S4_angles} angles) and reproduces the four RMSD keys to 1e-9.',
    'Table S8 (`tblS_pcm_trajectory`): the absolute-energy column is replaced by ΔE (kcal mol⁻¹, 3 decimals, relative '
    'to the final step, whose absolute energy is given in the footnote), because nine columns do not fit the text '
    'width in the PDF build; every header carries its unit (round 2, figures_r2 #6). The imaginary-mode checks '
    'are in the footnotes rather than in rows below the table.',
    'Table S9 (`tblS_protocol`): settings are reworded without internal code paths; the host machine is omitted; every '
    'number still comes from `protocol.*` / `software.*` / `model.*` and is checked against them.',
    'Column widths: pandoc takes relative widths from the dash counts of the separator line. The script sizes each '
    'column by its widest body cell (at least the longest word of its header) and sets explicit widths for Tables 7, '
    'S2 and S6, so that the PDF build (lualatex, 16 cm text) has no overfull cells and the DOCX has the same '
    'proportions.',
]


def summarise_checks(checks):
    """Readable check list: numeric details dropped, checks that differ only in numbers collapsed."""
    groups, order = {}, []
    for c in checks:
        c = re.sub(r' \(\|.*\| <= [^)]*\)$', '', c)
        tpl = re.sub(r'\d+', '0', c)
        if tpl not in groups:
            groups[tpl] = []
            order.append(tpl)
        if c not in groups[tpl]:
            groups[tpl].append(c)
    return [g[0] if len(g) == 1 else f'{g[0]} (and {len(g) - 1} more like it)' for g in (groups[t] for t in order)]


def write_index(tables):
    L = ['# Tables of the BTPC paper (benzotriazolium picrate, paper id `btpc`)', '',
         'Generated by `scripts/tables/make_tables.py` together with the tables; do not edit by hand. Rebuild with',
         '`python scripts/tables/make_tables.py --check` from the repository root (docs/REPRODUCE.md).',
         '',
         'Each `<stem>.md` holds one pandoc pipe table, then the caption line `Table: … {#tbl:<stem>}`, whose last '
         'sentences are the notes and footnotes of the table (round 1, figures_r1 #1). Insert a table into a manuscript with the include line given below; the path '
         'is relative to the manuscript: `papers/btpc/manuscript.md` uses `tables/<stem>.md`, the SI '
         '(`papers/btpc/si/manuscript.md`, built with the `S` number prefix) uses `../tables/<stem>.md`. Cite a table '
         'as `@tbl:<stem>` (for example `@tbl:tbl_geometry` → "Table 2"). The only cross-references inside the tables are '
         'the equation references `@eq:chi_eta` and `@eq:s_omega` in the note of `tbl_fmo_descriptors`, which resolve in '
         'the main text where that table sits; every other file builds in either document.', '',
         '**SI order matters.** The SI numbers below assume the SI includes the tables in this order. The draft '
         'caption of `figS_mep_extrema` (PLAN.md §4) calls the extremum table "Table S5"; keep `tblS_mep_extrema` '
         'fifth or change that caption.', '',
         '| No. | Stem (file) | Label | Where | Contents | Include line |', '|---|---|---|---|---|---|']
    for t in tables:
        inc = f'tables/{t.stem}.md' if t.where == 'main' else f'../tables/{t.stem}.md'
        L.append(f'| {t.number} | `{t.stem}.md` | `@tbl:{t.stem}` | {"main text" if t.where == "main" else "SI"} | '
                 f'{t.title} | `<!-- include: {inc} -->` |')
    L += ['', '## Facts keys and other inputs per table', '',
          'Keys are listed in the order the script reads them; `{a,b}.k` means `a.k` and `b.k`, `p.k01`…`k09` a run of '
          'numbered keys. "Other inputs" are the files read besides facts.json; "checks" are the cross-checks the script '
          'makes before writing (it stops if any fails).', '']
    for t in tables:
        L.append(f'### Table {t.number}: `{t.stem}` ({"main text" if t.where == "main" else "SI"})')
        L.append('')
        L.append(f'Caption: {t.caption}')
        L.append('')
        L.append(f'- Facts keys ({len(t.keys)}): ' + ', '.join(compress_keys(t.keys)))
        others = [f for f in t.files if f != str(FACTS_PATH.relative_to(REPO))]
        L.append('- Other inputs: ' + (', '.join(f'`{f}`' for f in others) if others else 'none'))
        if t.checks:
            L.append(f'- Checks ({len(t.checks)}): ' + '; '.join(summarise_checks(t.checks)))
        L.append('')
    L += ['## Conventions', '',
          '- Model names "IEF-PCM (water)" and "Gas phase"; column order X-ray, gas phase, IEF-PCM (round-1 revision, '
          'figures_r1 #13: the order of the figures and of the text), not the IEF-PCM-first order of PLAN §5.',
          '- Precision (STYLE.md): distances 3 decimals (Å), angles 1 decimal, orbital energies and descriptors 3 '
          'decimals (eV), softness 3 decimals, dipole 2 decimals, charge 3 decimals, shares 1 decimal (%), '
          'wavenumbers 1 decimal in the frequency tables with imaginary modes written 3.2i, energies 8 decimals (Ha), '
          'thermal terms 6 decimals (Ha), kcal mol⁻¹ 2 decimals, sign(λ₂)ρ 4 decimals and RDG 3 decimals as recorded, '
          'coordinates 6 decimals; RMSD/MAE statistics one decimal more than the quantity.',
          '- Typography: true minus (U+2212), en dash in bond names, ··· in hydrogen bonds, Unicode superscripts in '
          'units (crossref.lua turns them into Word superscripts), no `#` symmetry suffixes (labels cleaned with '
          '`figkit.geometry.clean_labels`); a lint pass in the script enforces these.',
          '- Bond names order atoms alphabetically by element and then by number (C12–N4, C7–N6), H last, as in '
          '`figS_bond_angle_deviations`; angle names put the central atom in the middle.', '',
          '## Deviations from PLAN.md (§5, §6) and why', '']
    nb, na = STATS.get('S4_keyed', ('?', '?'))
    L += [f'{k}. ' + d.replace('{S4_bonds}', str(nb)).replace('{S4_angles}', str(na)) for k, d in enumerate(DEVIATIONS, 1)]
    (HERE / 'INDEX.md').write_text('\n'.join(L) + '\n', encoding='utf-8')


# =================================================================================================
# pandoc checks
# =================================================================================================
def _run(cmd, cwd=None):
    r = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        raise TableError(f'pandoc failed ({" ".join(map(str, cmd))}):\n{r.stderr}')
    return r


def pandoc_check(tables):
    build = HERE / 'build'
    build.mkdir(exist_ok=True)
    report = []
    for t in tables:
        html = _run([PANDOC, '-f', 'markdown', '-t', 'html', HERE / f'{t.stem}.md']).stdout
        n_tab = html.count('<table')
        n_tr = len(re.findall(r'<tr[ >]', html))
        n_th = len(re.findall(r'<th[ >]', html))
        has_id = f'id="tbl:{t.stem}"' in html
        ok = n_tab == 1 and has_id and '<caption>' in html and n_tr == len(t.rows) + 1 and n_th == len(t.header)
        if not ok:
            raise TableError(f'{t.stem}: pandoc HTML has {n_tab} tables, {n_tr} rows (expected {len(t.rows) + 1}), '
                             f'{n_th} header cells (expected {len(t.header)}), table id present: {has_id}')
        report.append(f'{t.stem}: HTML ok ({len(t.rows)} rows x {len(t.header)} columns, caption + id)')
    one = next(t for t in tables if t.stem == 'tbl_geometry')
    for ext in ('docx', 'html'):
        out = build / f'{one.stem}.{ext}'
        _run([PANDOC, '-f', 'markdown', '-s', '--metadata', f'title={one.stem}', '-o', out, HERE / f'{one.stem}.md'])
    with zipfile.ZipFile(build / f'{one.stem}.docx') as z:
        doc = z.read('word/document.xml').decode('utf-8')
    if doc.count('<w:tbl>') != 1 or 'Selected bond lengths' not in doc:
        raise TableError('tbl_geometry.docx: expected one Word table and its caption')
    report.append(f'{one.stem}: DOCX ok (1 Word table, caption) and standalone HTML written')
    for where, prefix in (('main', None), ('SI', 'S')):
        ts = [t for t in tables if t.where == where]
        name = 'test_main' if where == 'main' else 'test_si'
        md = ['---', f'title: Table include test ({where})'] + ([f'number-prefix: {prefix}'] if prefix else []) + ['---', '']
        for t in ts:
            md += [f'See @tbl:{t.stem}.', '', f'<!-- include: ../{t.stem}.md -->', '']
        # equations cited in table notes live in the manuscript; stand-ins keep crossref.lua from failing here
        eqs = sorted({e for t in ts for e in re.findall(r'@(eq:[A-Za-z0-9_-]+)',
                                                           (HERE / f'{t.stem}.md').read_text(encoding='utf-8'))})
        md += [f'$$x$$ {{#{e}}}\n' for e in eqs]
        (build / f'{name}.md').write_text('\n'.join(md), encoding='utf-8')
        for ext in ('docx', 'html'):
            _run([PANDOC, f'{name}.md', '-f', 'markdown', '--lua-filter', SHARED_TOOLS / 'include.lua',
                  '--lua-filter', SHARED_TOOLS / 'crossref.lua', '--metadata', f'crossref-map={name}.map.json',
                  '-s', '-o', f'{name}.{ext}'], cwd=build)
        cmap = json.loads((build / f'{name}.map.json').read_text())
        got = [(k, v['num']) for k, v in cmap.items() if v['kind'] == 'tbl']
        want = [(f'tbl:{t.stem}', str(k)) for k, t in enumerate(ts, 1)]
        if got != want:
            raise TableError(f'{name}: crossref numbering {got} != {want}')
        html = (build / f'{name}.html').read_text()
        pre = prefix or ''
        for k, t in enumerate(ts, 1):
            if f'Table {pre}{k}' not in html.replace(' ', ' '):
                raise TableError(f'{name}.html: reference/caption "Table {pre}{k}" not found')
        report.append(f'{name}: include.lua + crossref.lua ok, DOCX and HTML, tables numbered '
                      f'{pre}1–{pre}{len(ts)} in order and every @tbl reference resolved')
    return report


# =================================================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--check', action='store_true', help='also run the pandoc checks (build/)')
    ap.add_argument('--only', nargs='*', help='build only these stems')
    args = ap.parse_args()
    tables, nmain, nsi = [], 0, 0
    # QA r4-m6: number the main tables in the order the manuscript cites them (crossref.json of the built paper),
    # not in the order their functions appear in this file
    rank = {s: k for k, s in enumerate(MAIN_ORDER)}
    check(sorted(s for s, w, _, _ in REGISTRY if w == 'main') == sorted(MAIN_ORDER), 'MAIN_ORDER lists every main table')
    reg = sorted([r for r in REGISTRY if r[1] == 'main'], key=lambda r: rank[r[0]]) + [r for r in REGISTRY if r[1] != 'main']
    cmap_path = PAPER / 'build' / 'crossref.json'
    if cmap_path.exists():
        cm = json.loads(cmap_path.read_text())
        built = [k[4:] for k, v in sorted(((k, v) for k, v in cm.items() if v.get('kind') == 'tbl'),
                                          key=lambda kv: int(kv[1]['num']))]
        check(built == MAIN_ORDER, f'main-table order = order in build/crossref.json ({built})')
    for stem, where, title, fn in reg:
        if where == 'main':
            nmain += 1
            number = str(nmain)
        else:
            nsi += 1
            number = f'S{nsi}'
        if args.only and stem not in args.only:
            continue
        CTX.stem, CTX.keys, CTX.files, CTX.checks = stem, [], [str(FACTS_PATH.relative_to(REPO))], []
        spec = fn()
        t = Table(stem=stem, where=where, title=title, number=number, keys=list(CTX.keys), files=list(CTX.files),
                  checks=list(CTX.checks), **spec)
        text = render(t)
        lint(text, stem)
        (HERE / f'{stem}.md').write_text(text, encoding='utf-8')
        tables.append(t)
        print(f'Table {number:>3}  {stem:<22} {len(t.rows):3d} rows  {len(t.keys):3d} keys  {len(t.checks):3d} checks')
    if not args.only:
        write_index(tables)
        print(f'INDEX.md written ({len(tables)} tables)')
    if args.check:
        for line in pandoc_check(tables):
            print('pandoc:', line)


if __name__ == '__main__':
    try:
        main()
    except TableError as e:
        print(f'ERROR {e}', file=sys.stderr)
        sys.exit(1)
