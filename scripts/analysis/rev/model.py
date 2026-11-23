"""Revision facts: numbers derived from the model definition for the method justification (R1.4)."""
import json
from .common import *

def build(add, facts):
    case = json.load(open(ROOT / 'molecules/btpc/cases/btpc_pcm.json'))
    xyz = (ROOT / case['xyz']).read_text().split('\n')
    n = int(xyz[0].split()[0]); el = [l.split()[0] for l in xyz[2:2 + n]]
    nh = el.count('H'); nheavy = n - nh
    src = case['xyz'] + ' (element count)'
    # 6-311++G(d,p): one diffuse sp shell (s + 3 p) on C, N, O and one diffuse s on H
    add('rev.model.n_diffuse_functions', nheavy * 4 + nh * 1, '', None, src,
        f'derived: {nheavy} heavy atoms x 4 (diffuse sp shell) + {nh} H x 1 (diffuse s)')
    # aug-cc-pVTZ with spherical functions: C, N, O 46 functions (5s4p3d2f with aug), H 23 (4s3p2d)
    naug = nheavy * 46 + nh * 23
    add('rev.model.n_basis_augccpvtz', naug, '', None, src, f'derived: {nheavy} x 46 + {nh} x 23 spherical functions')
    nb = facts['model.n_basis']['value']
    r = naug / nb
    add('rev.model.augccpvtz_ratio', r, '', None, 'rev.model.n_basis_augccpvtz / model.n_basis', 'derived')
    add('rev.model.augccpvtz_cost_ratio_range', [r ** 3, r ** 4], '', None, 'rev.model.augccpvtz_ratio',
        'derived: ratio cubed and to the fourth power (cost scaling of a density-fitted hybrid between N^3 and N^4)')
    add('rev.model.gas_min_per_gradient', facts['gas.opt_wall']['value'] / facts['gas.opt_steps']['value'], 'min',
        'btpc_gas', 'gas.opt_wall / gas.opt_steps', 'derived: wall time per energy + gradient, 4 threads')
    add('rev.model.pcm_min_per_gradient', float(__import__('re').search(r'about ([\d.]+) h', facts['pcm.opt_wall']['value']).group(1)) * 60 / (facts['pcm.opt_steps']['value'] + 1), 'min', 'btpc_pcm',
        'pcm.opt_wall ("about 13.4 h") / (pcm.opt_steps + 1)', 'derived: about 11 min per energy + gradient; thread count not logged')

def _versions():
    import re
    p = RUNS / 'vm_extras' / 'env.versions.txt'
    out = {}
    for line in p.read_text().splitlines():
        m = re.match(r'^(\S+)\s+(\d[\w.]*)\s+\S+$', line)
        if m:
            out[m.group(1)] = m.group(2)
    return out, rel(p)

_build0 = build
def build(add, facts):
    _build0(add, facts)
    v, src = _versions()
    for pkg in ('pyscf', 'geometric', 'libxc-c', 'pyscf-dispersion', 'libopenblas'):
        add(f'rev.software.{pkg.replace("-", "_")}', v[pkg], '', None, f'{src}:{pkg}',
            'version in the conda environment of the revision calculations (same lock file on both VMs)')
