"""
Relaxed scan of the inter-ion torsion of the BTPC ion pair (plan section 5, item 6; reviewer 1 point 3).

The torsion is C1-O1...N4-C12 for btpc (atom indices 4, 0, 18, 21; the one tabulated by
revision/round1/audit/evidence/dft/dihedrals.py: X-ray -164.2, gas minimum -117.8, PCM minimum -121.6 deg).  By
default it is found from the case: the acceptor oxygen's carbon, the acceptor O, the donor N and the donor N's ring
carbon (case['proton']).  The scan runs from -180 to -60 deg in 15 deg steps (9 points), gas phase, with B3LYP or
B3LYP-D3(BJ) (--disp), as two chains outward from the grid point nearest to the start minimum: down to -180 and up to
-60, each point started from the optimised previous point.  Every other coordinate is relaxed with geomeTRIC (the
torsion held by a $set dihedral constraint) to the five criteria of the published optimisations.  For each point:
energy (also relative to the unconstrained start minimum, kcal/mol), achieved torsion, N...O and H positions (proton
transfer?), inter-ring angle, steps, convergence, and the optimised geometry (points/*.xyz).

Usage:
  python dft/revision/conf_scan.py --job scan_b3lyp --out <dir> [--disp none|d3bj] [--start stored|PATH]
         [--grid -180:-60:15] [--torsion auto|C1,O1,N4,C12] [--maxsteps 150] [--retry-unconverged] [--case btpc]
  For --disp d3bj, --start should be the B3LYP-D3(BJ) gas minimum (an opt_d3.py output directory or its opt.xyz).
--smoke: STO-3G, grid 2, three points (the start grid point and one each side), loose criteria, 3 steps per point.
"""
import sys, pathlib, argparse, time, shutil
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import rev_common as rc          # first: pins the thread pools
import numpy as np
import headline as hl
from opt_d3 import CONV, CONV_LOOSE, read_traj_last, geomeric_norms


def parser():
    p = rc.add_common_args(argparse.ArgumentParser(description=__doc__.split('\n')[1]))
    p.add_argument('--disp', default='none', help='none (B3LYP) or d3bj')
    p.add_argument('--xc', default=rc.XC)
    p.add_argument('--case', default='btpc')
    p.add_argument('--model', choices=('gas', 'pcm'), default='gas')
    p.add_argument('--start', default='stored')
    p.add_argument('--grid', default='-180:-60:15', help='start:stop:step in degrees (inclusive)')
    p.add_argument('--torsion', default='auto', help='auto, or four labels / 0-based indices, comma separated')
    p.add_argument('--maxsteps', type=int, default=None)
    p.add_argument('--retry-unconverged', action='store_true')
    return p


def torsion_atoms(spec, case, syms, X):
    if spec == 'auto':
        t = hl.hbond_torsion_atoms(case, syms, X)
        if not t:
            raise ValueError('no proton/H-bond definition in the case; give --torsion')
        return list(t)
    parts = [s.strip() for s in spec.split(',')]
    if all(p.isdigit() for p in parts):
        return [int(p) for p in parts]
    lab = hl.labels_of(case, len(syms))
    return [lab.index(p) for p in parts]


def grid_values(spec):
    a, b, s = (float(x) for x in spec.split(':'))
    n = int(round((b - a) / s))
    return [round(a + i * s, 6) for i in range(n + 1)]


def angdiff(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


def run(job, args):
    from pyscf.geomopt import geometric_solver
    lv = hl.level(args.smoke)
    disp = None if args.disp in (None, '', 'none') else args.disp
    solvent = None if args.model == 'gas' else 'water'
    conv = CONV_LOOSE if args.smoke else CONV
    maxsteps = args.maxsteps or (3 if args.smoke else 150)
    scf_tol = 1e-6 if args.smoke else (rc.CONV_TOL_OPT_PCM if solvent else rc.CONV_TOL_OPT_GAS)
    case = rc.load_case(f'{args.case}_{args.model}')
    start = hl.resolve_geom(args.model, args.case)[1] if args.start == 'stored' else pathlib.Path(args.start)
    if not start.is_absolute():
        start = rc.ROOT / start
    if start.is_dir():
        start = start / 'opt.xyz'
    job.input_file('start_geometry', start)
    syms0, X0 = rc.read_xyz(start)
    tor = torsion_atoms(args.torsion, case, syms0, X0)
    lab = hl.labels_of(case, len(syms0))
    grid = grid_values(args.grid)
    job.r['settings'] = {**lv, 'xc': args.xc, 'disp': disp, 'model': args.model, 'scf_conv_tol_opt': scf_tol,
                         'geometric_criteria': conv, 'maxsteps_per_point': maxsteps, 'grid_deg': grid,
                         'torsion_atoms_0based': tor, 'torsion_labels': '-'.join(lab[i] for i in tor),
                         'constraint': 'geomeTRIC $set dihedral (1-based indices in the constraints file)'}
    job.save()
    pdir = job.out / 'points'; pdir.mkdir(exist_ok=True)

    # ---- the unconstrained start minimum (energy reference)
    if not job.done('start'):
        mol = rc.mol_for_case(case, start, basis=lv['basis'], verbose=3, output=job.out / 'pyscf_start.log')
        mf = rc.make_mf(mol, xc=args.xc, disp=disp, solvent=solvent, conv_tol=lv['conv_tol'],
                        grid_level=lv['grid_level'], lebedev_order=lv['lebedev_order'],
                        chkfile=job.chk_dir / 'start.chk')
        wall = rc.run_scf(mf, log=job.log)
        g = mf.nuc_grad_method().kernel()
        rec = {'E_Ha': float(mf.e_tot), 'torsion_deg': hl.dihedral(X0, tor), 'scf_wall_s': round(wall, 1),
               **geomeric_norms(g), 'geometry': hl.geometry_metrics(case, syms0, X0), 'xyz': str(start)}
        if disp:
            rec['E_dispersion_Ha'] = float(np.asarray(mf.scf_summary.get('dispersion', np.nan)))
        job.put('start', rec)
        job.log(f"start: E = {mf.e_tot:.10f} Ha, torsion {rec['torsion_deg']:.1f} deg, grad rms/max "
                f"{rec['grms']:.1e}/{rec['gmax']:.1e}")
        mol.stdout.close()
    st = job.get('start')
    phi0 = st['torsion_deg']
    p0 = min(grid, key=lambda v: abs(angdiff(v, phi0)))
    down = [v for v in grid if v <= p0][::-1]
    up = [v for v in grid if v > p0]
    if args.smoke:
        down, up = down[:2], up[:1]
        if not up:
            down = down[:3]
    chains = [('down', down), ('up', up)]
    job.r['chains'] = {'start_torsion_deg': phi0, 'first_point_deg': p0, 'down': down, 'up': up}
    job.save()
    ok = True

    def optimise(val, guess_xyz, guess_chk):
        key = f'p{val:+.1f}'
        it = job.get(key)
        if it and it.get('status') == 'done' and (it.get('converged') or not args.retry_unconverged):
            return it
        ptraj = pdir / f'{key}_traj.xyz'
        geom = guess_xyz
        if ptraj.exists():
            s, Xl, nfr = read_traj_last(ptraj)
            if nfr:
                geom = [(a, tuple(x)) for a, x in zip(s, Xl)]
                job.log(f'{key}: continuing from frame {nfr} of {ptraj.name}')
        mol = rc.mol_for_case(case, geom, basis=lv['basis'], verbose=3, output=job.out / 'pyscf_scan.log')
        chk = job.chk_dir / f'{key}.chk'
        if not chk.exists() and guess_chk and pathlib.Path(guess_chk).exists():
            shutil.copy(guess_chk, chk)
        mf = rc.make_mf(mol, xc=args.xc, disp=disp, solvent=solvent, conv_tol=scf_tol,
                        grid_level=lv['grid_level'], lebedev_order=lv['lebedev_order'], chkfile=chk)
        if chk.exists():
            mf.init_guess = 'chkfile'
        cons = job.out / f'constraints_{key}.txt'
        cons.write_text('$set\ndihedral {} {} {} {} {:.4f}\n'.format(*[i + 1 for i in tor], val))
        hist = []
        rec = {'target_deg': val, 'start_from': str(guess_xyz) if not isinstance(guess_xyz, list) else 'previous point',
               'history': hist}
        job.put(key, rec, status='running')
        tf = open(ptraj, 'a')
        state = {'X_prev': None, 't': time.time(), 'edisp': None}

        def cb(env):
            m = env['mol']; e = float(env['energy']); g = np.asarray(env['gradients'])
            s = [m.atom_pure_symbol(i) for i in range(m.natm)]; X = m.atom_coords() * rc.BOHR_A
            nrm = geomeric_norms(g, X, state['X_prev']); state['X_prev'] = X
            try:
                state['edisp'] = float(np.asarray(env['g_scanner'].base.scf_summary.get('dispersion', np.nan)))
            except Exception:
                pass
            hist.append({'E_Ha': e, 'torsion_deg': hl.dihedral(X, tor), **nrm, 'wall_s': round(time.time() - state['t'], 1)})
            state['t'] = time.time()
            tf.write(f'{len(s)}\n{key} step {len(hist) - 1} E = {e:.10f}\n' +
                     ''.join(f'{a:2s} {x:16.10f} {y:16.10f} {z:16.10f}\n' for a, (x, y, z) in zip(s, X)))
            tf.flush()
            job.put(key, rec, status='running')
            job.log(f"{key} step {len(hist) - 1}: E = {e:.10f}  torsion {hist[-1]['torsion_deg']:.2f}  grad rms/max "
                    f"{nrm['grms']:.1e}/{nrm['gmax']:.1e}")

        t0 = time.time()
        converged, mol_opt = geometric_solver.kernel(mf, maxsteps=maxsteps, callback=cb, constraints=str(cons), **conv)
        tf.close()
        s, X = hl.mol_geometry(mol_opt)
        rc.write_xyz(pdir / f'{key}.xyz', s, X, f'{key} torsion {hl.dihedral(X, tor):.2f} E = {hist[-1]["E_Ha"]:.10f}')
        rec.update({'converged': bool(converged), 'steps': len(hist), 'E_Ha': hist[-1]['E_Ha'],
                    'torsion_deg': hl.dihedral(X, tor), 'geometry': hl.geometry_metrics(case, s, X),
                    'xyz': f'points/{key}.xyz', 'chk': str(chk), 'wall_min': round((time.time() - t0) / 60, 2),
                    'dE_vs_start_kcal': (hist[-1]['E_Ha'] - st['E_Ha']) * rc.HARTREE_KCAL})
        if disp and state['edisp'] is not None:
            rec['E_dispersion_Ha'] = state['edisp']
        job.put(key, rec)
        job.r['timings'][key] = round(time.time() - t0, 1)
        job.log(f"{key}: {'converged' if converged else 'NOT converged'} in {len(hist)} steps, E = {rec['E_Ha']:.10f} Ha "
                f"({rec['dE_vs_start_kcal']:+.2f} kcal/mol vs start), torsion {rec['torsion_deg']:.2f}, N...O "
                f"{rec['geometry'].get('d_N...O_A', float('nan')):.3f} A, proton transferred "
                f"{rec['geometry'].get('proton_transferred')}")
        mol.stdout.close()
        return job.get(key)

    first_rec = None
    for name, vals in chains:
        prev_xyz, prev_chk = start, job.chk_dir / 'start.chk'
        if name == 'up':
            if first_rec is None:
                first_rec = job.get(f'p{p0:+.1f}')
            if first_rec and first_rec.get('status') == 'done':
                prev_xyz, prev_chk = job.out / first_rec['xyz'], pathlib.Path(first_rec['chk'])
        for v in vals:
            r = optimise(v, prev_xyz, prev_chk)
            if name == 'down' and v == p0:
                first_rec = r
            ok &= bool(r.get('converged')) or args.smoke
            prev_xyz, prev_chk = job.out / r['xyz'], pathlib.Path(r['chk'])
    pts = sorted([v for k, v in job.r['items'].items() if k.startswith('p') and v.get('status') == 'done'],
                 key=lambda r: r['target_deg'])
    emin = min([r['E_Ha'] for r in pts] + [st['E_Ha']])
    job.r['summary'] = {'start': {'E_Ha': st['E_Ha'], 'torsion_deg': phi0},
                        'points': [{'target_deg': r['target_deg'], 'torsion_deg': r['torsion_deg'], 'E_Ha': r['E_Ha'],
                                    'dE_vs_start_kcal': r['dE_vs_start_kcal'],
                                    'dE_vs_lowest_kcal': (r['E_Ha'] - emin) * rc.HARTREE_KCAL,
                                    'converged': r['converged'], 'steps': r['steps'],
                                    'N...O_A': r['geometry'].get('d_N...O_A'), 'H-N_A': r['geometry'].get('d_H-N_A'),
                                    'proton_transferred': r['geometry'].get('proton_transferred'),
                                    'inter_ring_angle_deg': r['geometry'].get('inter_ring_angle_deg')} for r in pts]}
    job.save()
    return ok


if __name__ == '__main__':
    rc.install_sigterm_handler()
    rc.main_wrapper(parser, run, __file__)
