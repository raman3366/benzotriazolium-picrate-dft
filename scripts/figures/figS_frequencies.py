"""SI figure figS_frequencies: positions of all harmonic wavenumbers, and the finite-difference step check.

(a) the 93 unscaled harmonic wavenumbers of the gas-phase (up) and IEF-PCM (water) (down) ion pairs as
    unit-height sticks (no intensities were computed); the N-H modes assigned in facts.json (*.nu_N4_H4N_*,
    *.nu_N6_H6N_*) are coloured (N4-H4N teal, N6-H6N magenta) and tagged nu (stretch), delta (in-plane bend) or
    gamma (out-of-plane wag);
(b) IEF-PCM: wavenumber from the 0.01 Bohr Hessian (freq.json, authoritative) minus that from the 0.005 Bohr
    Hessian (freq_step0.005.json), mode by mode; the largest change and its mode are asserted against
    saddle.step_change_max, _max_signed and _max_mode, the mean |change| against saddle.step_change_mean.

    python scripts/figures/figS_frequencies.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

from figkit import style
from figkit.draw2d import HBOND_COLOUR
from figkit.labels import LabelPlacer
import _struct_common as C
import _vib as VB

STEM = 'figS_frequencies'
WIDTH_MM = 175.0
GROUP_COLOURS = {'N4': HBOND_COLOUR, 'N6': C.MODE_COLOUR}
TYPES = (('stretch', 'ν'), ('inplane_bend', 'δ'), ('oop_wag', 'γ'))
XLIM = (-80.0, 3760.0)


def assigned_modes(model, f):
    """dict mode index -> (group, symbol) for the N-H modes that facts.json assigns."""
    out = {}
    for grp in ('N4', 'N6'):
        for key, sym in TYPES:
            v = C.fact(f'{model}.nu_{grp}_H{grp[1:]}N_{key}')
            vals = [float(t) for t in str(v).split(',')]
            for x in vals:
                k = int(np.argmin(np.abs(f - x)))
                assert abs(f[k] - x) < 0.051, (model, grp, key, x, f[k])
                out[k] = (grp, sym)
    return out


def main():
    style.apply()
    f = {m: VB.freq_json(m)[0] for m in C.MODELS}
    f005 = VB.freq_json('pcm', 'freq_step0.005.json')[0]
    assert all(len(v) == 93 for v in (f['gas'], f['pcm'], f005))
    C.check_close(f['pcm'][0], 'saddle.freq_imag_step0p01', 1e-9)
    C.check_close(f005[0], 'saddle.freq_imag_step0p005', 1e-9)
    dnu = f['pcm'] - f005
    kmax = int(np.argmax(np.abs(dnu)))
    C.check_equal(abs(dnu[kmax]), 'saddle.step_change_max', 2, 'largest step change')
    # the printed change and its mode against their facts.json keys (QA round 3, r3-m1): the signed change, the
    # mode's wavenumber in the authoritative 0.01 Bohr Hessian and the mean |change|
    C.check_close(dnu[kmax], 'saddle.step_change_max_signed', 1e-9)
    C.check_close(f['pcm'][kmax], 'saddle.step_change_max_mode', 1e-9)
    C.check_close(float(np.abs(dnu).mean()), 'saddle.step_change_mean', 1e-9)
    # the 0.005 Bohr wavenumber of the same mode (1522.6 cm-1) has no key of its own; the notes of both keys
    # quote it and it is asserted against them
    import re
    for key, pat in (('saddle.step_change_max', r'mode at (\d+\.\d) cm-1'),
                     ('saddle.step_change_max_mode', r'lies at (\d+\.\d) cm-1 in the 0\.005 Bohr Hessian')):
        nu_note = float(re.search(pat, C.facts()[key]['note']).group(1))
        assert f'{f005[kmax]:.1f}' == f'{nu_note:.1f}', (key, f005[kmax], nu_note)
    fig = plt.figure(figsize=(C.mm(WIDTH_MM), C.mm(96.0)))
    H = 96.0
    xl, wr = 14.0, WIDTH_MM - 14.0 - 1.2
    # ---- (a) sticks ----
    ax = C.mm_box(fig, xl, 47.0, wr, 41.0)
    C.header(fig, 0.0, H - 3.4, 'a', 'Harmonic wavenumbers (unscaled; stick height carries no intensity)')
    ax.axhline(0, color='0.3', lw=0.6, zorder=2)
    pl_items = []
    tips = []
    for m, sgn in (('gas', 1.0), ('pcm', -1.0)):
        asg = assigned_modes(m, f[m])
        for k, x in enumerate(f[m]):
            if k in asg:
                continue
            ax.plot([x, x], [0, sgn], color='0.62', lw=0.6, zorder=1, solid_capstyle='butt')
        for k, (grp, sym) in asg.items():
            x = f[m][k]
            ax.plot([x, x], [0, sgn], color=GROUP_COLOURS[grp], lw=1.1, zorder=3, solid_capstyle='butt')
            tips.append((x, sgn, grp, sym))
    ax.set_xlim(*XLIM)
    ax.set_ylim(-1.55, 1.55)
    ax.set_yticks([])
    ax.tick_params(labelbottom=True, labelsize=style.SMALL_PT, pad=1.5)       # QA r1 m14: tick labels on (a)
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.text(1920, 0.55, style.MODEL_LABELS['gas'], ha='center', va='center', fontsize=style.TICK_PT)
    ax.text(1920, -0.55, style.MODEL_LABELS['pcm'], ha='center', va='center', fontsize=style.TICK_PT)
    # tags at the stick tips, placed without overlaps
    fig.canvas.draw()
    pl = LabelPlacer(ax, name='mode tags', n_rings=4, clearance_pt=0.6)
    for ln in ax.lines:
        xy = ln.get_xydata()
        if len(xy) == 2:
            pl.add_segment(xy[0], xy[1], hw_pt=0.5)
    for t in ax.texts:
        pl.add_text_obstacle(t)
    # tags of two sticks of the same group and type closer than 10 cm-1 (IEF-PCM gamma(N4-H4N) at 959.8 and
    # 963.0 cm-1) are stacked vertically over the pair, not scattered by the placer (QA r1 m14)
    tips.sort()
    fixed = set()
    for a, b in zip(tips, tips[1:]):
        if a[1:] == b[1:] and abs(b[0] - a[0]) < 10.0:
            xm = (a[0] + b[0]) / 2
            for k, t in enumerate((a, b)):
                ax.annotate(t[3], (xm, t[1]), xytext=(0, t[1] * (3.2 + 7.0 * k)), textcoords='offset points',
                            ha='center', va='bottom' if t[1] > 0 else 'top', fontsize=style.SMALL_PT,
                            color=GROUP_COLOURS[t[2]])
            fixed.update({a, b})
    fig.canvas.draw()
    for t in ax.texts[-2 * (len(fixed) // 2):] if fixed else []:
        pl.add_text_obstacle(t)
    for x, sgn, grp, sym in tips:
        if (x, sgn, grp, sym) in fixed:
            continue
        pl.add_item(sym, (x, sgn), radius_pt=0.8, prefer=(0, sgn), priority=0, key=(x, sgn),
                    kw=dict(fontsize=style.SMALL_PT, color=GROUP_COLOURS[grp]), max_ring=4)
    placed, dropped = pl.place()
    assert not dropped, [d.text for d in dropped]
    hs = [Line2D([], [], color=GROUP_COLOURS['N4'], lw=1.1, label='N4–H4N (hydrogen-bonded)'),
          Line2D([], [], color=GROUP_COLOURS['N6'], lw=1.1, label='N6–H6N (free)'),
          Line2D([], [], color='0.62', lw=0.6, label='other modes')]
    # key on its own line under (a): colours on the left, the meaning of the tags on the right
    W = WIDTH_MM
    y_key = 38.6
    fig.legend(handles=hs, loc='center left', bbox_to_anchor=(xl / W, y_key / H), ncol=3, frameon=False,
               fontsize=style.SMALL_PT, handlelength=1.6, columnspacing=1.2, borderaxespad=0.0, borderpad=0.0)
    fig.text((xl + wr) / W, y_key / H, 'ν stretch, δ in-plane bend, γ out-of-plane wag', ha='right',
             va='center', fontsize=style.SMALL_PT)
    # ---- (b) step check ----
    axb = C.mm_box(fig, xl, 10.5, wr, 18.0)
    C.header(fig, 0.0, 31.4, 'b', 'IEF-PCM (water): ν(0.01 Bohr step) − ν(0.005 Bohr step)')
    axb.axhline(0, color='0.6', lw=0.5, ls=(0, (3, 2)), zorder=1)
    axb.plot(f['pcm'], dnu, 'o', ms=2.2, color='0.25', mec='none', zorder=3)
    # highlight in black: the 1520.8 cm-1 mode is not an N6-H6N mode, so the N6-H6N magenta would mislead (m14)
    axb.plot([f['pcm'][kmax]], [dnu[kmax]], 'o', ms=3.6, mfc='none', mec='black', mew=0.8, zorder=4)
    axb.annotate(f'{C.fmt_num(dnu[kmax], 2)} cm$^{{-1}}$ (mode at {f["pcm"][kmax]:.1f} cm$^{{-1}}$; '
                 f'{f005[kmax]:.1f} cm$^{{-1}}$ with the 0.005 Bohr step)',
                 (f['pcm'][kmax], dnu[kmax]), xytext=(6, -1),
                 textcoords='offset points', ha='left', va='center', fontsize=style.SMALL_PT, color='black')
    axb.set_xlim(*XLIM)
    lim = np.ceil(np.abs(dnu).max() * 2) / 2 + 0.3
    axb.set_ylim(-lim, lim)
    axb.yaxis.set_major_locator(MultipleLocator(1))
    axb.xaxis.set_major_locator(MultipleLocator(500))
    axb.xaxis.set_minor_locator(MultipleLocator(100))
    axb.set_ylabel('Δν (cm$^{-1}$)')
    axb.set_xlabel('Wavenumber (cm$^{-1}$)')
    print(f'largest |dnu| {abs(dnu[kmax]):.3f} cm-1 at {f["pcm"][kmax]:.1f} cm-1; mean |dnu| {np.abs(dnu).mean():.3f}')
    C.save(fig, STEM, WIDTH_MM)


if __name__ == '__main__':
    main()
