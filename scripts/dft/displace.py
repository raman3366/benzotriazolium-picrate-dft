"""
Escape a saddle point: displace the optimised geometry along each imaginary normal mode and write the
result as the starting geometry for a fresh optimisation.

Needed because a floppy hydrogen-bonded cluster can converge on gradient while sitting on a ridge -- PZC's
gas-phase structure came out with three imaginary modes (-93.7, -82.9, -42.4 cm-1).  Standard practice is to
step along the imaginary modes and re-optimise.

Usage: python scripts/dft/displace.py <case> [--scale 0.3]
Writes results/<case>/preopt.xyz (picked up automatically by optimize.py) and backs up the invalid results.
"""
import sys, os, shutil, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import *
from vibrations import AU2CM, AMU2ME
from pyscf.data import elements

case = load_case(sys.argv[1]); d = case['dir']
scale = float(sys.argv[sys.argv.index('--scale') + 1]) if '--scale' in sys.argv else 0.3

mol = mol_for_case(case, d / 'opt.xyz', verbose=0)
hess = np.load(d / 'hessian.npy')
natm = mol.natm
mass = np.array([elements.MASSES[mol.atom_charge(i)] for i in range(natm)]) * AMU2ME
x = mol.atom_coords()

Hm = hess.transpose(0, 2, 1, 3).reshape(3 * natm, 3 * natm)
w = np.repeat(mass ** -0.5, 3)
Hm = 0.5 * ((Hm * w[:, None] * w[None, :]) + (Hm * w[:, None] * w[None, :]).T)

com = (mass[:, None] * x).sum(0) / mass.sum(); r = x - com
TR = np.zeros((6, 3 * natm))
for i in range(natm):
    s = np.sqrt(mass[i])
    TR[0, 3*i], TR[1, 3*i+1], TR[2, 3*i+2] = s, s, s
    for k, ax in enumerate(np.eye(3)):
        TR[3 + k, 3*i:3*i+3] = s * np.cross(ax, r[i])
q, _ = np.linalg.qr(TR.T)
P = np.eye(3 * natm) - q @ q.T
ev, vec = np.linalg.eigh(P @ Hm @ P)
freq = np.sign(ev) * np.sqrt(np.abs(ev)) * AU2CM

# imaginary modes, excluding the six projected translations/rotations
order = np.argsort(np.abs(freq)); tr = set(order[:6].tolist())
imag = [i for i in range(len(freq)) if freq[i] < -1.0 and i not in tr]
if not imag:
    print(f'[displace] {case["name"]}: no imaginary modes -- nothing to do'); sys.exit(0)

print(f'[displace] {case["name"]}: {len(imag)} imaginary modes '
      f'{[round(float(freq[i]), 1) for i in imag]} cm-1')
step = np.zeros(3 * natm)
for i in imag:
    v = vec[:, i] * w                      # mass-weighted -> Cartesian
    step += v / np.linalg.norm(v)
step = step / np.linalg.norm(step) * scale * np.sqrt(natm)   # Bohr

new = (x + step.reshape(natm, 3)) * BOHR_A
syms = [mol.atom_symbol(i) for i in range(natm)]
stamp = datetime_stamp = __import__('datetime').datetime.now().strftime('%H%M')
bak = d.parent / f'{case["name"]}_saddle_{stamp}'
shutil.copytree(d, bak, ignore=shutil.ignore_patterns('*.cube', '*.chk', 'vtk'))
print(f'[displace] invalid saddle-point results archived to {bak.name}')
for f in ('freq.json', 'RESULTS.md', 'opt.xyz', 'opt_summary.json', 'status.json', 'hessian.npy'):
    (d / f).unlink(missing_ok=True)
write_xyz(d / 'preopt.xyz', syms, new,
          f'displaced {scale} along {len(imag)} imaginary modes of the saddle point; re-optimise from here')
stage_note(case, f'saddle point: displaced along {len(imag)} imaginary modes '
                 f'{[round(float(freq[i]),1) for i in imag]} cm-1, re-optimising')
print(f'[displace] wrote {d / "preopt.xyz"} -- relaunch the pipeline for this case')
