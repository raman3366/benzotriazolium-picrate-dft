"""Shared helpers for the BTPC structure, crystal, geometry-comparison, trajectory and vibration figures.

Imported by the scripts in this directory (fig_scheme, fig_crystal, fig_structures, fig_xray_comparison,
fig_hbond_trajectory, fig_vibrations, figS_*). Run any of them as

    python scripts/figures/<stem>.py      (from the repository root; writes into figures/)

Everything printed on a figure is either computed here from the data files or read from facts.json; the
helpers `fact`, `check_equal` and `fmt_num` make that explicit.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from figkit import style, geometry as geo, qa

REPO = style.REPO
PAPER = REPO / 'papers' / 'btpc'
FIG = PAPER / 'figures'
BUILD = FIG / 'build'
FACTS_PATH = PAPER / 'data' / 'facts.json'
CASES = {'gas': 'btpc_gas', 'pcm': 'btpc_pcm'}
MODELS = ('gas', 'pcm')
MINUS = '−'

# ---------------------------------------------------------------------------------------------------------------
# atom indices (0-based, CIF order; identical in opt.xyz, btpc_cluster.xyz, cubes and hessian.npy)
# ---------------------------------------------------------------------------------------------------------------
PIC = list(range(18))
CAT = list(range(18, 33))
PIC_HEAVY = list(range(16))
PIC_RING = [4, 5, 6, 8, 9, 3]            # C1 C2 C3 C4 C5 C6
CAT9 = list(range(18, 27))               # N4 N6 C7 C12 N5 C11 C10 C9 C8
HB = (18, 31, 0)                         # N4, H4N, O1
I = dict(O1=0, N3=1, O7=2, C6=3, C1=4, C2=5, C3=6, O6=7, C4=8, C5=9, N1=10, N2=11, O4=12, O5=13, O3=14,
         O2=15, H5=16, H3=17, N4=18, N6=19, C7=20, C12=21, N5=22, C11=23, C10=24, C9=25, C8=26, H11=27,
         H8=28, H10=29, H9=30, H4N=31, H6N=32)
NITRO = {'N1': (10, 5, (15, 14)), 'N2': (11, 8, (12, 13)), 'N3': (1, 3, (7, 2))}   # N, ring C, (O, O)

HBOND_GREY = '#9A9A9A'
MODE_COLOUR = '#C51B7D'
LIT_GREY = '#8C8C8C'
PIC_COLOUR = '#762A83'
CAT_COLOUR = '#1B7837'


def model_marker(model, ms):
    """One marker per model in every figure of the paper (QA round 2, r2-m6): X-ray filled circle, gas phase
    filled triangle, IEF-PCM open square (open, so that a gas-phase point under it stays visible).  Returns
    Line2D keyword arguments (marker, ms, colours); `ms` is the size of the circle, the others are matched."""
    c = style.MODEL_COLOURS[model]
    if model == 'xray':
        return dict(marker='o', ms=ms, color=c, mfc=c, mec='white', mew=0.4)
    if model == 'gas':
        return dict(marker='^', ms=1.2 * ms, color=c, mfc=c, mec='white', mew=0.3)
    if model == 'pcm':
        return dict(marker='s', ms=0.92 * ms, color=c, mfc='none', mec=c, mew=0.85)
    raise KeyError(model)


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
    """Value of a facts.json key (raises KeyError when missing)."""
    return facts()[key]['value']


def fmt_num(v, nd):
    s = f'{v:.{nd}f}'
    return s.replace('-', MINUS)


def check_equal(value, key, nd, what=''):
    """Assert that `value` printed with `nd` decimals equals the facts.json value at the same precision.
    Returns the facts value (float)."""
    fv = float(fact(key))
    a, b = f'{value:.{nd}f}', f'{fv:.{nd}f}'
    if a != b:
        raise AssertionError(f'{what or key}: computed {a} but facts.json {key} = {b}')
    return fv


def check_close(value, key, tol, what=''):
    fv = float(fact(key))
    if abs(value - fv) > tol:
        raise AssertionError(f'{what or key}: computed {value} vs facts.json {fv} (tol {tol})')
    return fv


# ---------------------------------------------------------------------------------------------------------------
# geometries
# ---------------------------------------------------------------------------------------------------------------
_GEOM = None


def geometries():
    """dict(symbols, labels, xray, gas, pcm, case) with (33, 3) arrays in Å and clean CIF labels."""
    global _GEOM
    if _GEOM is None:
        cp = geo.read_case(CASES['pcm'])
        cg = geo.read_case(CASES['gas'])
        xr = geo.read_xyz(cp['paths'].xray_xyz)
        pcm = geo.read_xyz(cp['paths'].opt_xyz)
        gas = geo.read_xyz(cg['paths'].opt_xyz)
        labels = cp['labels']
        assert labels == cg['labels'] and pcm.symbols == gas.symbols == xr.symbols
        assert all('#' not in l for l in labels)
        _GEOM = dict(symbols=list(pcm.symbols), labels=list(labels), xray=np.asarray(xr.xyz, float),
                     gas=np.asarray(gas.xyz, float), pcm=np.asarray(pcm.xyz, float), case=cp, case_gas=cg)
    return _GEOM


def plane_normal(X, idx):
    P = np.asarray(X)[list(idx)]
    c = P.mean(0)
    _, _, vt = np.linalg.svd(P - c)
    return vt[2], c


def interplanar(X):
    """Angle (deg) between the least-squares planes of the picrate ring and the 9 cation ring atoms."""
    n1, _ = plane_normal(X, PIC_RING)
    n2, _ = plane_normal(X, CAT9)
    c = abs(float(n1 @ n2))
    return float(np.degrees(np.arccos(min(1.0, c))))


def nitro_twist(X, name):
    """Twist (deg, 0-90) of a nitro group: angle between the O-N-O plane and the picrate ring plane."""
    n_idx, c_idx, (oa, ob) = NITRO[name]
    X = np.asarray(X)
    n_ring, _ = plane_normal(X, PIC_RING)
    v = np.cross(X[oa] - X[n_idx], X[ob] - X[n_idx])
    v /= np.linalg.norm(v)
    a = float(np.degrees(np.arccos(min(1.0, abs(float(v @ n_ring))))))
    return a


# ---------------------------------------------------------------------------------------------------------------
# common camera (PLAN.md section 2): picrate left, cation right, tilt 35 deg about the N4-H4N...O1 bond
# ---------------------------------------------------------------------------------------------------------------
_VIEWS = None
TILT = 35.0
PHI = 210.0


def views():
    """dict model -> (R, focal) for draw2d orient=(R, focal), models 'pcm', 'gas', 'xray'; plus the View objects.

    V_pcm = auto_view(mol_pcm, mode='hbond', atoms=PIC_HEAVY, hbonds=[HB], tilt=35), rolled so the vector from
    the picrate ring centroid to the cation 9-atom centroid points to screen right; V_gas / V_xray are the same
    camera carried over by a Kabsch fit (gas: picrate heavy atoms + N4 + H4N; X-ray: picrate heavy atoms)."""
    global _VIEWS
    if _VIEWS is not None:
        return _VIEWS
    from figkit import render3d as r3
    g = geometries()
    mol_pcm = r3.load_molecule(CASES['pcm'])
    mol_gas = r3.load_molecule(CASES['gas'])
    mol_x = r3.load_molecule(CASES['pcm'], xyz=g['case']['paths'].xray_xyz)
    # face view of the picrate, tilted by TILT towards the screen direction PHI (the auto_view 'hbond' tilt
    # direction, 30 deg, shows the IEF-PCM cation edge-on; PHI = 210 deg keeps both rings and the N4-H4N...O1
    # NCI disc visible in all three structures, see PLAN.md section 2)
    V = r3.auto_view(mol_pcm, mode='face', atoms=PIC_HEAVY)
    V = V.roll(PHI).elevation(TILT).roll(-PHI)
    s = V.screen(np.vstack([g['pcm'][PIC_RING].mean(0), g['pcm'][CAT9].mean(0)]))
    d = s[1, :2] - s[0, :2]
    theta = float(np.degrees(np.arctan2(d[1], d[0])))
    # roll(deg) rotates the picture counter-clockwise: rotate by -theta to bring the vector onto +x
    V = V.roll(-theta)
    s = V.screen(np.vstack([g['pcm'][PIC_RING].mean(0), g['pcm'][CAT9].mean(0)]))
    d = s[1, :2] - s[0, :2]
    if abs(np.degrees(np.arctan2(d[1], d[0]))) > 0.01:
        V = V.roll(2 * theta)                     # sign convention the other way round
        s = V.screen(np.vstack([g['pcm'][PIC_RING].mean(0), g['pcm'][CAT9].mean(0)]))
        d = s[1, :2] - s[0, :2]
        assert abs(np.degrees(np.arctan2(d[1], d[0]))) < 0.01
    Vg = r3.matched_view(V, mol_pcm, mol_gas, atoms=PIC_HEAVY + [18, 31])
    Vx = r3.matched_view(V, mol_pcm, mol_x, atoms=PIC_HEAVY)
    out = {}
    for k, v in (('pcm', V), ('gas', Vg), ('xray', Vx)):
        R = np.vstack([v.right, v.up, v.direction])
        out[k] = (R, np.asarray(v.focal, float))
    out['views'] = {'pcm': V, 'gas': Vg, 'xray': Vx}
    BUILD.mkdir(parents=True, exist_ok=True)
    (BUILD / 'views_struct.json').write_text(json.dumps({k: v.as_dict() for k, v in out['views'].items()},
                                                         indent=1))
    _VIEWS = out
    return out


def oriented(model, X=None):
    """Coordinates of a model (or of X, given in that model's frame) in the common viewing frame."""
    R, c = views()[model]
    X = geometries()[model] if X is None else X
    return geo.apply_orientation(X, R, c)


# ---------------------------------------------------------------------------------------------------------------
# saving + print-size preview
# ---------------------------------------------------------------------------------------------------------------
def save(fig, stem, width_mm, preview=True):
    """qa.save in strict mode (600 dpi PNG + PDF), then a 200-dpi preview at print width into build/."""
    import os
    draft = bool(os.environ.get('FIG_DRAFT'))           # development only: report QA issues without failing
    out = qa.save(fig, (BUILD / 'draft' / stem) if draft else (FIG / stem), strict=not draft,
                  max_width_mm=width_mm)
    if draft:
        for p in out:
            print('DRAFT', p.relative_to(REPO))
        return out
    if preview:
        make_preview(stem, width_mm)
    for p in out:
        print('wrote', p.relative_to(REPO))
    return out


def make_preview(stem, width_mm, dpi=200):
    from PIL import Image
    src = FIG / f'{stem}.png'
    BUILD.mkdir(parents=True, exist_ok=True)
    with Image.open(src) as im:
        w_mm = im.width / (im.info.get('dpi', (600, 600))[0] / 25.4)
        tw = int(round(w_mm / 25.4 * dpi))
        th = int(round(im.height * tw / im.width))
        im.convert('RGB').resize((tw, th), Image.LANCZOS).save(BUILD / f'{stem}_print200.png', dpi=(dpi, dpi))
    print(f'preview build/{stem}_print200.png  ({w_mm:.1f} mm wide; target {width_mm:.0f} mm)')


def mm(x):
    return x * style.MM


# ---------------------------------------------------------------------------------------------------------------
# side-by-side molecule panels at one common print scale
# ---------------------------------------------------------------------------------------------------------------
def extent(symbols, xy, margin=0.0, atom_scale=1.0):
    """(xmin, xmax, ymin, ymax) of the ball-and-stick discs of `xy` (Å) plus `margin` (Å)."""
    r = np.array([style.BALL_RADII.get(s, 0.35) * atom_scale for s in symbols])
    xy = np.asarray(xy)[:, :2]
    lo = (xy - r[:, None]).min(0) - margin
    hi = (xy + r[:, None]).max(0) + margin
    return np.array([lo[0], hi[0], lo[1], hi[1]])


def hide_frame(ax):
    """axis('off') plus invisible patch and spines, so a tight bounding box hugs what is drawn."""
    ax.axis('off')
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def molecule_row(fig, boxes, x0_mm, y0_mm, width_mm, gap_mm=3.0, scale_mm_per_A=None, valign='center'):
    """Axes side by side for drawings whose data extents are `boxes` [(xmin, xmax, ymin, ymax), ...] (Å),
    all at one scale (mm per Å). The row starts at (x0_mm, y0_mm) from the figure's lower-left corner and is
    `width_mm` wide (the scale is chosen to fill it unless given). Returns (axes, scale_mm_per_A, height_mm)."""
    boxes = [np.asarray(b, float) for b in boxes]
    wA = np.array([b[1] - b[0] for b in boxes])
    hA = np.array([b[3] - b[2] for b in boxes])
    s = scale_mm_per_A or (width_mm - gap_mm * (len(boxes) - 1)) / wA.sum()
    H = hA.max() * s
    W_fig, H_fig = fig.get_size_inches() / style.MM
    axes = []
    x = x0_mm
    for b, w, h in zip(boxes, wA, hA):
        ax = fig.add_axes([x / W_fig, y0_mm / H_fig, w * s / W_fig, H / H_fig])
        hide_frame(ax)
        cy = (b[2] + b[3]) / 2
        half = H / s / 2
        if valign == 'top':
            ax.set_ylim(b[3] - 2 * half, b[3])
        else:
            ax.set_ylim(cy - half, cy + half)
        ax.set_xlim(b[0], b[1])
        ax.set_aspect('equal', adjustable='box')
        ax.set_autoscale_on(False)
        axes.append(ax)
        x += w * s + gap_mm
    return axes, s, H


def content_bbox(ax, renderer=None):
    """Display-space bbox of everything drawn in `ax` (patches, lines, texts, collections, legend); the axes
    background and spines are ignored."""
    from matplotlib.transforms import Bbox
    from matplotlib.text import Text
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    from matplotlib.collections import Collection
    from matplotlib.legend import Legend
    from matplotlib.spines import Spine
    fig = ax.figure
    if renderer is None:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
    bbs = []
    for a in ax.get_children():
        if not a.get_visible() or a is ax.patch or isinstance(a, Spine):
            continue
        if isinstance(a, Text) and not a.get_text().strip():
            continue
        if isinstance(a, (Patch, Line2D, Text, Legend, Collection)):
            try:
                bb = a.get_window_extent(renderer)
            except Exception:
                continue
            if np.all(np.isfinite([bb.x0, bb.y0, bb.x1, bb.y1])) and (bb.width > 0 or bb.height > 0):
                bbs.append(bb)
    return Bbox.union(bbs)


def tight_width_mm(fig):
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    return bb.width / style.MM, bb.height / style.MM


def fit_to_width(build, width_mm, s0, tol_mm=0.6, max_iter=6, max_height_mm=None):
    """Call build(scale) -> fig repeatedly, adjusting the drawing scale (mm per Å) until the tight bounding box
    of the figure is `width_mm` wide (within tol, never wider). Returns (fig, scale)."""
    import matplotlib.pyplot as plt
    # savefig adds rcParams['savefig.pad_inches'] on every side of the tight box
    width_mm = width_mm - 2 * float(plt.rcParams.get('savefig.pad_inches', 0.0)) * 25.4
    pts = []
    s = s0
    fig = None
    for _ in range(max_iter):
        if fig is not None:
            plt.close(fig)
        fig = build(s)
        w, h = tight_width_mm(fig)
        if max_height_mm is not None and h > max_height_mm:
            w = max(w, width_mm * h / max_height_mm)
        pts.append((s, w))
        if width_mm - tol_mm <= w <= width_mm:
            return fig, s
        if len(pts) >= 2:
            (s1, w1), (s2, w2) = pts[-2], pts[-1]
            a = (w2 - w1) / (s2 - s1) if abs(s2 - s1) > 1e-9 else w2 / s2
            s = s2 + (width_mm - tol_mm / 2 - w2) / a
        else:
            s = s * (width_mm - tol_mm / 2) / w
    if w > width_mm:
        plt.close(fig)
        s *= 0.98 * width_mm / w
        fig = build(s)
    return fig, s


def panel_title(fig, x_px, y_px, letter, title=None, fontsize=None):
    """Bold panel letter and an optional plain title, with the text baseline at display point (x_px, y_px)."""
    fontsize = style.BASE_PT if fontsize is None else fontsize
    t = fig.text(x_px / fig.bbox.width, y_px / fig.bbox.height, f'({letter})', fontsize=style.PANEL_PT,
                 fontweight='bold', ha='left', va='baseline')
    arts = [t]
    if title:
        fig.canvas.draw()
        bb = t.get_window_extent(fig.canvas.get_renderer())
        arts.append(fig.text((bb.x1 + 3.0 * fig.dpi / 72) / fig.bbox.width, y_px / fig.bbox.height, title,
                             fontsize=fontsize, ha='left', va='baseline'))
    return arts


def crop_row(fig, axes, gap_mm=3.0, pad_mm=0.6, top_mm=0.0, common_y=True, x_offsets_mm=None):
    """Shrink a row of equal-aspect drawing axes to what is drawn in them, keeping the print scale: each axes gets
    the x range of its content, all share the y range of the union (common_y) so the drawings stay aligned, and
    the figure is resized to the row (plus `top_mm` for titles). Returns the figure size in mm."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    ext, ipd = [], []
    for ax in axes:
        bb = content_bbox(ax, r)
        inv = ax.transData.inverted()
        (x0, y0), (x1, y1) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
        ext.append([x0, x1, y0, y1])
        xl = ax.get_xlim()
        ipd.append(ax.bbox.width / fig.dpi / (xl[1] - xl[0]))
    if common_y:
        Y0, Y1 = min(e[2] for e in ext), max(e[3] for e in ext)
        for e in ext:
            e[2], e[3] = Y0, Y1
    sizes = []
    for e, k in zip(ext, ipd):
        pad = pad_mm * style.MM / k
        e[0] -= pad; e[1] += pad; e[2] -= pad; e[3] += pad
        sizes.append(((e[1] - e[0]) * k, (e[3] - e[2]) * k))
    W = sum(w for w, h in sizes) + gap_mm * style.MM * (len(axes) - 1)
    Hc = max(h for w, h in sizes)
    H = Hc + top_mm * style.MM
    fig.set_size_inches(W, H)
    x = 0.0
    for ax, e, (w, h) in zip(axes, ext, sizes):
        ax.set_position([x / W, (Hc - h) / H, w / W, h / H])
        ax.set_xlim(e[0], e[1])
        ax.set_ylim(e[2], e[3])
        x += w + gap_mm * style.MM
    fig.canvas.draw()
    return W / style.MM, H / style.MM


# ---------------------------------------------------------------------------------------------------------------
# two-pass drawings for multi-row figures: measure what a drawing occupies, then place it exactly
# ---------------------------------------------------------------------------------------------------------------
def measure_drawing(draw, box, s_mm_per_A, pad_A=2.0):
    """Draw with draw(ax) into a scratch axes covering `box` (xmin, xmax, ymin, ymax; Å) plus pad_A at
    s mm per Å, and return the data-space box of everything drawn (atoms, bonds, labels, arrows)."""
    import matplotlib.pyplot as plt
    b = np.asarray(box, float) + np.array([-pad_A, pad_A, -pad_A, pad_A])
    w, h = (b[1] - b[0]) * s_mm_per_A, (b[3] - b[2]) * s_mm_per_A
    fig = plt.figure(figsize=(w * style.MM, h * style.MM))
    ax = fig.add_axes([0, 0, 1, 1])
    hide_frame(ax)
    ax.set_xlim(b[0], b[1])
    ax.set_ylim(b[2], b[3])
    ax.set_aspect('equal', adjustable='box')
    ax.set_autoscale_on(False)
    draw(ax)
    fig.canvas.draw()
    bb = content_bbox(ax, fig.canvas.get_renderer())
    (x0, y0), (x1, y1) = ax.transData.inverted().transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
    plt.close(fig)
    return np.array([x0, x1, y0, y1])


def drawing_axes(fig, x_mm, y_mm, box, s_mm_per_A, pad_mm=0.5):
    """Axes at (x_mm, y_mm) (lower-left, from the figure's lower-left corner) showing `box` (Å) at s mm per Å
    plus pad_mm on every side. Returns (ax, width_mm, height_mm)."""
    pad = pad_mm / s_mm_per_A
    b = np.asarray(box, float) + np.array([-pad, pad, -pad, pad])
    w, h = (b[1] - b[0]) * s_mm_per_A, (b[3] - b[2]) * s_mm_per_A
    W, H = fig.get_size_inches() / style.MM
    ax = fig.add_axes([x_mm / W, y_mm / H, w / W, h / H])
    hide_frame(ax)
    ax.set_xlim(b[0], b[1])
    ax.set_ylim(b[2], b[3])
    ax.set_aspect('equal', adjustable='box')
    ax.set_autoscale_on(False)
    return ax, w, h


def mm_box(fig, x_mm, y_mm, w_mm, h_mm):
    """Axes rectangle in mm from the lower-left corner of the figure."""
    W, H = fig.get_size_inches() / style.MM
    return fig.add_axes([x_mm / W, y_mm / H, w_mm / W, h_mm / H])


def header(fig, x_mm, y_mm, letter, title=None, handles=None, fontsize=None, handle_len=1.6):
    """Panel header at (x_mm, y_mm) = baseline start: bold '(a)', optional title, optional inline legend
    (handles) after it. Returns the list of artists."""
    fontsize = style.TICK_PT if fontsize is None else fontsize
    W, H = fig.get_size_inches() / style.MM
    arts = [fig.text(x_mm / W, y_mm / H, f'({letter})', fontsize=style.PANEL_PT, fontweight='bold',
                     ha='left', va='baseline')]
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    x = arts[-1].get_window_extent(r).x1 / fig.dpi / style.MM + 1.2
    if title:
        arts.append(fig.text(x / W, y_mm / H, title, fontsize=style.BASE_PT, ha='left', va='baseline'))
        fig.canvas.draw()
        x = arts[-1].get_window_extent(r).x1 / fig.dpi / style.MM + 1.5
    if handles:
        leg = fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(x / W, (y_mm - 1.1) / H),
                         ncol=len(handles), fontsize=fontsize, frameon=False, handlelength=handle_len,
                         columnspacing=1.0, handletextpad=0.4, borderaxespad=0.0, borderpad=0.0)
        arts.append(leg)
    return arts


def fit_row_scale(draw_fns, width_mm, gap_mm, box_fns, s0=5.0, pad_mm=0.5, iters=5, common_y=True):
    """Scale (mm per Å) at which the measured drawings draw_fns[i](ax) (initial data boxes box_fns[i]) sit side
    by side in width_mm with gap_mm between them. Returns (scale, boxes) with boxes sharing one y range when
    common_y."""
    s = s0
    boxes = None
    for _ in range(iters):
        boxes = [measure_drawing(f, b, s) for f, b in zip(draw_fns, box_fns)]
        if common_y:
            y0, y1 = min(b[2] for b in boxes), max(b[3] for b in boxes)
            boxes = [np.array([b[0], b[1], y0, y1]) for b in boxes]
        w = sum(b[1] - b[0] for b in boxes) * s + gap_mm * (len(boxes) - 1) + 2 * pad_mm * len(boxes)
        if width_mm - 0.8 <= w <= width_mm:
            break
        s *= (width_mm - 0.4 - gap_mm * (len(boxes) - 1) - 2 * pad_mm * len(boxes)) / \
             (w - gap_mm * (len(boxes) - 1) - 2 * pad_mm * len(boxes))
    return s, boxes


# ---------------------------------------------------------------------------------------------------------------
# the IEF-PCM optimisation trajectory (Fig. 5, figS_trajectory_full)
# ---------------------------------------------------------------------------------------------------------------
TRAJ_CHECK_STEPS = (0, 25, 30, 40, 45, 50, 60, 70)


def pcm_trajectory():
    """Series along data/r2/.../btpc_pcm/opt_traj.xyz (71 frames): step, E (Ha), dE_kcal (relative to the final
    frame, geometry.HARTREE_KCAL), N4-H4N, H4N...O1, N4...O1, H4N...O7, H11...O2 (Å), angle N4-H4N...O1 and
    the inter-ring angle and nitro twists (deg). Asserts the eight saddle.traj_stepNN checkpoints of facts.json
    (distances to 0.001 Å, angles to 0.1 deg, E exact to 1e-8 Ha) and the total energy drop."""
    import re
    path = REPO / 'data' / 'r2' / 'molecules' / 'btpc' / 'results' / CASES['pcm'] / 'opt_traj.xyz'
    frames = geo.read_xyz_frames(path)
    assert len(frames) == int(fact('pcm.opt_steps')) + 1, len(frames)
    sym = geometries()['symbols']
    assert all(f.symbols == sym for f in frames)
    steps = []
    for f in frames:
        m = re.match(r'\s*step\s+(\d+)', f.comment)
        steps.append(int(m.group(1)))
    assert steps == list(range(len(frames)))
    E = np.array([f.energy for f in frames])
    out = dict(step=np.array(steps), E=E, dE_kcal=(E - E[-1]) * geo.HARTREE_KCAL)
    pairs = {'N4–H4N': (18, 31), 'H4N···O1': (31, 0), 'N4···O1': (18, 0), 'H4N···O7': (31, 2), 'H11···O2': (27, 15)}
    for k, (i, j) in pairs.items():
        out[k] = np.array([geo.distance(f.xyz, i, j) for f in frames])
    out['angle'] = np.array([geo.angle(f.xyz, 18, 31, 0) for f in frames])
    out['interplanar'] = np.array([interplanar(f.xyz) for f in frames])
    for n in ('N1', 'N2', 'N3'):
        out[f'twist_{n}'] = np.array([nitro_twist(f.xyz, n) for f in frames])
    num = r'(-?\d+\.\d+)'
    for s in TRAJ_CHECK_STEPS:
        v = str(fact(f'saddle.traj_step{s:02d}'))
        e = float(re.search(r'E ' + num, v).group(1))
        want = {'N4–H4N': r'N4-H4N ' + num, 'H4N···O1': r'H4N\.\.\.O1 ' + num, 'N4···O1': r'N4\.\.\.O1 ' + num,
                'H11···O2': r'H11\.\.\.O2 ' + num}
        assert abs(E[s] - e) < 1e-8, (s, E[s], e)
        for k, pat in want.items():
            w = float(re.search(pat, v).group(1))
            assert abs(out[k][s] - w) <= 0.0005 + 1e-9, (s, k, out[k][s], w)
        a = float(re.search(r'angle ' + num, v).group(1))
        ip = float(re.search(r'interplanar ' + num, v).group(1))
        assert abs(out['angle'][s] - a) <= 0.05 + 1e-9 and abs(out['interplanar'][s] - ip) <= 0.05 + 1e-9, \
            (s, out['angle'][s], a, out['interplanar'][s], ip)
    check_close(out['dE_kcal'][0] * -1, 'pcm.opt_energy_drop_total', 1e-6, 'energy drop over the optimisation')
    return out


def literature_claim(key, pattern):
    """Numbers parsed (regex groups) from an abstract-checked claim in papers/btpc/lit/candidates.jsonl."""
    import re
    for line in (PAPER / 'lit' / 'candidates.jsonl').read_text().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        if d['key'] == key:
            assert d.get('abstract_checked') is True, key
            m = re.search(pattern, ' '.join(d['claims']))
            assert m, (key, pattern)
            return [float(g) for g in m.groups()]
    raise KeyError(key)


# ---------------------------------------------------------------------------------------------------------------
# leaders that pass an atom ball (QA round 3, r3-m2)
# ---------------------------------------------------------------------------------------------------------------
def segment_outside_discs(p0, p1, centres, radii, own_tol=0.0):
    """Pieces [(a, b), ...] of the 2D segment p0 -> p1 that lie outside every disc (centres (N, 2), radii (N,)),
    in order from p0.  A disc that holds p0 (the leader's own atom, whose ball it starts from, or a surface point
    seen in front of a ball; distance <= radius + own_tol) is not cut.  Used to draw a leader that crosses an
    atom ball as passing BEHIND the ball (the part over the disc is left out) instead of across it."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    c = np.atleast_2d(np.asarray(centres, float))
    r = np.atleast_1d(np.asarray(radii, float))
    d = p1 - p0
    a = float(d @ d)
    if a == 0.0:
        return [(p0, p1)]
    own = np.hypot(*(c - p0).T) <= r + own_tol
    cuts = []
    for ci, ri in zip(c[~own], r[~own]):
        f = p0 - ci
        b = 2.0 * float(f @ d)
        q = float(f @ f) - ri * ri
        disc = b * b - 4.0 * a * q
        if disc <= 0.0:
            continue
        s = np.sqrt(disc)
        lo, hi = max((-b - s) / (2 * a), 0.0), min((-b + s) / (2 * a), 1.0)
        if hi > lo:
            cuts.append((lo, hi))
    pieces, t = [], 0.0
    for lo, hi in sorted(cuts):
        if lo > t:
            pieces.append((t, lo))
        t = max(t, hi)
    if t < 1.0:
        pieces.append((t, 1.0))
    return [(p0 + t0 * d, p0 + t1 * d) for t0, t1 in pieces if t1 > t0]


def over_discs(p0, p1, centres, radii, start, tol=0.0, n=400):
    """Indices of the discs (other than those holding `start`) that the drawn piece p0 -> p1 enters by more than
    `tol` (the check 'no leader drawn over a ball')."""
    c = np.atleast_2d(np.asarray(centres, float))
    r = np.atleast_1d(np.asarray(radii, float))
    own = np.hypot(*(c - np.asarray(start, float)).T) <= r + tol
    t = np.linspace(0.0, 1.0, n)[:, None]
    P = np.asarray(p0, float) * (1 - t) + np.asarray(p1, float) * t
    dmin = (np.linalg.norm(P[:, None, :] - c[None, :, :], axis=2) - r[None, :]).min(0)
    return [int(k) for k in np.nonzero((dmin < -tol) & ~own)[0]]
