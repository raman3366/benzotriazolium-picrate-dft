"""SI figure figS_elf_lol_picrate: ELF and LOL in the picrate ring plane (through C6, C4 and C2), gas phase (top
row) and IEF-PCM (water) (bottom row), one colour scale 0-1; the companion of Fig. 12 (fig_elf_lol).

Data: data/r2/molecules/btpc/results/btpc_{gas,pcm}/{elf,lol}_picrate_ring_plane.txt (Multiwfn, 300 x 300); atoms
from molecules/btpc/results/<case>/opt.xyz. Shared rotation, centre (picrate ring centroid) and window for the four
panels; rotation = grid-aligned angle nearest to the Fig. 3 camera for the picrate atoms (see _elf_lol.py).
Asserted: ELF/LOL share the plane frame; C6, C4, C2 lie in the plane and the other ring atoms within 0.05 Å; the
map ranges equal `*.elf/lol_picrate_ring_range`; the O7 offsets equal `*.picrate_plane_offset_O7` (O7 drawn in
the gas phase only).

    python scripts/figures/figS_elf_lol_picrate.py      (from the repository root; writes into figures/)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from figkit import style
import _elec_common as C
import _elf_lol as E

STEM = 'figS_elf_lol_picrate'
PLANE = 'picrate_ring'
PIC_RING = [C.I[n] for n in ('C1', 'C2', 'C3', 'C4', 'C5', 'C6')]
DEFINING = [C.I['C6'], C.I['C4'], C.I['C2']]          # protocol.elf_plane_picrate


def main():
    style.apply()
    assert all(n in C.fact('protocol.elf_plane_picrate') for n in ('C6', 'C4', 'C2'))
    G, data = E.load(PLANE, DEFINING, PIC_RING, ring_tol=0.05)
    for m in C.MODELS:
        d = data[m]['dist']
        print(f'{m}: C1 {d[C.I["C1"]]:+.3f} Å, O1 {d[C.I["O1"]]:+.3f} Å, H4N {d[C.I["H4N"]]:+.3f} Å from the plane')
    geoms, win, rot, mirror = E.common_geometry(G, data, centre_idx=PIC_RING, cam_idx=list(range(0, 18)))
    # caption (QA round 2, r2-m9): the maps are the view of Fig. 3 turned clockwise in the page by about 25 deg; the
    # plane is seen from the side of the Fig. 3 camera (the map is reflected to match it), so the picrate is not
    # mirrored with respect to Fig. 3
    for m in C.MODELS:
        a, mi = E.camera_angle(data[m]['elf'], data[m]['X'], m, list(range(0, 18)))
        turn = (rot - a + 180) % 360 - 180
        print(f'{m}: maps turned {turn:+.1f} deg (anticlockwise > 0) from the Fig. 3 view')
        assert mi == mirror and abs(turn + 25) < 2.5, (m, turn)
    for m in C.MODELS:
        ins = geoms[m]['inplane']
        assert all(ins[i] for i in PIC_RING + [C.I['O1'], C.I['H4N']])
        # O7 against its facts.json key (QA round 3, r3-m1): drawn in the gas phase (0.43 Å), not in IEF-PCM
        # (0.50 Å, just beyond the 0.5 Å cut; caption)
        C.check_equal(abs(float(data[m]['dist'][C.I['O7']])), f'{m}.picrate_plane_offset_O7', 3, f'{m} O7 offset')
        assert bool(ins[C.I['O7']]) == (m == 'gas'), (m, ins[C.I['O7']])
        print(m, 'in plane:', ', '.join(G['labels'][i] for i in np.nonzero(ins)[0]))
    E.figure(G, data, geoms, win, STEM, label_h=(C.I['H4N'],))


if __name__ == '__main__':
    main()
