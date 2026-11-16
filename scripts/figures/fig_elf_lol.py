"""Fig. 12 fig_elf_lol: ELF and LOL in the plane of the benzotriazolium cation (through N4, N5 and C7),
gas phase (top row) and IEF-PCM (water) (bottom row), one colour scale 0-1.

Data: data/r2/molecules/btpc/results/btpc_{gas,pcm}/{elf,lol}_benzotriazolium_ring_plane.txt (Multiwfn, 300 x 300);
atoms from molecules/btpc/results/<case>/opt.xyz. All four panels share one rotation, centre (the centroid of the
nine cation ring atoms) and window, so the cation is drawn identically in every panel and the different position
of O1 in the two models reads directly. The rotation is the grid-aligned one closest to the camera of Fig. 3
(computed from build/views_struct.json and printed; see _elf_lol.py). Asserted: ELF and LOL share the plane frame;
the cation ring atoms lie in the plane; |O1 offset| equals `*.cation_plane_offset_O1` and |O7 offset| (printed in
(a, b)) `*.cation_plane_offset_O7`; the map ranges equal
`*.elf/lol_benzotriazolium_ring_range`.

    python scripts/figures/fig_elf_lol.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from figkit import style
import _elec_common as C
import _elf_lol as E

STEM = 'fig_elf_lol'
PLANE = 'benzotriazolium_ring'
CAT9 = list(range(18, 27))                           # N4 N6 C7 C12 N5 C11 C10 C9 C8
DEFINING = [C.I['N4'], C.I['N5'], C.I['C7']]         # protocol.elf_plane_cation


def main():
    style.apply()
    assert all(n in C.fact('protocol.elf_plane_cation') for n in ('N4', 'N5', 'C7'))
    G, data = E.load(PLANE, DEFINING, CAT9, ring_tol=0.01)
    for m in C.MODELS:
        # |offset| as in facts (its sign depends on the direction of the normal)
        d = abs(data[m]['dist'][C.I['O1']])
        C.check_equal(d if C.ffact(f'{m}.cation_plane_offset_O1') > 0 else -d, f'{m}.cation_plane_offset_O1', 3)
    geoms, win, rot, mirror = E.common_geometry(G, data, centre_idx=CAT9, cam_idx=list(range(18, 33)))
    # caption (QA round 2, r2-m9): the maps are the view of Fig. 3 turned anticlockwise in the page by about 40 deg
    # (gas phase) and 25 deg (IEF-PCM), not mirrored
    for m, want in (('gas', 40), ('pcm', 25)):
        a, mi = E.camera_angle(data[m]['elf'], data[m]['X'], m, list(range(18, 33)))
        turn = (rot - a + 180) % 360 - 180
        print(f'{m}: maps turned {turn:+.1f} deg (anticlockwise > 0) from the Fig. 3 view')
        assert not mi and abs(turn - want) < 2.5, (m, turn)
    for m in C.MODELS:
        assert all(geoms[m]['inplane'][i] for i in (C.I['H4N'], C.I['H6N'], C.I['O1'], C.I['N4']))
    # QA round 1, m11: the feature at the left edge of the gas-phase maps (green in ELF, a blue disc in LOL) is
    # the density of O7, 1.26 Å out of the plane; its projection gets a grey cross and tag in (a, b)
    o7 = C.I['O7']
    d7 = abs(float(data['gas']['dist'][o7]))
    xy7 = geoms['gas']['xy'][o7]
    assert not geoms['gas']['inplane'][o7] and win[0] < xy7[0] < win[0] + 0.5, (d7, xy7)
    # the printed offset (2 decimals) against its facts.json key (QA round 3, r3-m1); IEF-PCM: O7 far from the
    # plane (key asserted at 3 decimals), no feature
    C.check_equal(d7, 'gas.cation_plane_offset_O7', 2, 'gas-phase O7 offset from the cation plane')
    C.check_equal(abs(float(data['pcm']['dist'][o7])), 'pcm.cation_plane_offset_O7', 3, 'IEF-PCM O7 offset')
    assert C.ffact('pcm.cation_plane_offset_O7') > 2.0

    def tag_o7(axes):
        import matplotlib.patheffects as pe
        for kind in ('ELF', 'LOL'):
            ax = axes[('gas', kind)]
            ax.plot([xy7[0]], [xy7[1]], marker='x', ms=3.2, mew=0.8, color='0.25', zorder=30)
            ax.annotate(f'O7 ({d7:.2f} Å off plane)', xy7, xytext=(2.5, 4.0), textcoords='offset points',
                        ha='left', va='bottom', fontsize=style.SMALL_PT, color='0.2', zorder=31,
                        path_effects=[pe.withStroke(linewidth=1.6, foreground='white')])
    print(f'gas O7: {d7:.2f} Å from the plane, projected at ({xy7[0]:.2f}, {xy7[1]:.2f})')
    E.figure(G, data, geoms, win, STEM, label_h=(C.I['H4N'], C.I['H6N']), extra=tag_o7)


if __name__ == '__main__':
    main()
