"""
B3LYP-D3(BJ) re-optimisation of the BTPC ion pair in the gas phase or in IEF-PCM water (plan section 5, items 4-5),
followed by the frontier-orbital analysis at the new minimum.

Optimiser: geomeTRIC 1.1.1 through PySCF (TRIC coordinates), with the five convergence criteria of the published
optimisations (dft/optimize.py: geomeTRIC's defaults = Gaussian's), all required:
  |dE| < 1.0e-6 Ha, RMS gradient < 3.0e-4, max gradient < 4.5e-4 Ha/Bohr, RMS step < 1.2e-3, max step < 1.8e-3 A;
SCF conv_tol 1e-9 (gas) / 1e-8 (PCM) during the optimisation, as dft/optimize.py; maxsteps 300.
Start: the stored B3LYP minimum of the same model (molecules/btpc/results/btpc_<model>/opt.xyz) or --start.
Initial Hessian (--seed-hessian auto): the stored finite-difference Hessian of that model
(data/r2/molecules/btpc/results/btpc_<model>/hessian.npy, (natm, natm, 3, 3) in Ha/Bohr^2, the layout of
dft/frequencies.py), reshaped to the (3N, 3N) Cartesian matrix geomeTRIC expects, with overall translations and
rotations given a curvature of 0.05 Ha/Bohr^2 (geomeTRIC's guess value) and internal eigenvalues floored at 5e-5
Ha/Bohr^2, so that the seed is positive definite and geomeTRIC does not discard it at the first update.  It is passed
as geomeTRIC's `hess_data`.  (Passing hessian='file:...' through PySCF's geometric_solver would instead compute an
analytic Hessian and overwrite that file -- see pyscf/geomopt/geometric_solver.py _make_hessian.)
At the minimum: a conv_tol 1e-10 single point with the HOMO/LUMO energies, fmo.py descriptors, fragment shares of
HOMO-1..LUMO+3, Mulliken/Loewdin/Hirshfeld/CM5 fragment charges, dipole, inter-ion geometry and its change from the
B3LYP minimum, the dispersion energy, a padded molden file (for Multiwfn) and chk/final.chk.

Restart: the trajectory (opt_traj.xyz) is written every step; a restarted job continues from its lowest-energy frame.
An optimisation that ends without meeting the criteria is kept as 'stopped' (not done), so rerunning the same command
continues it from the trajectory.

Caps and stalls (21 Nov 2026; needed because PySCF 2.14's geometric_solver.optimize() returns the last geometry without
any flag when maxsteps is reached, and btpc_pcm crawled for 70 steps on a flat inter-ion mode):
  --maxsteps N     total evaluations over all restarts (not per run); default 300, smoke 4.
  --max-hours H    total optimisation wall time over all restarts (sum of the per-step times in result.json).
  --negligible-stop W  (0 = off, the default) stop early on a stall: the gradient criteria of the active set (geomeTRIC
                   per-atom norms, RMS and max) met at each of the last W evaluations AND the energy spread over them
                   below --negligible-de (default 1e-6 Ha, geomeTRIC's per-step energy criterion, i.e. ten steps
                   together move less than one step may).  Only the displacement criteria are then unmet; the result
                   is accepted as "negligible forces (stall)".  With it on, a run that reaches a cap is also accepted
                   if the rule used for btpc_pcm and in proton_scan.py holds (gradient criteria met at the last step,
                   energy spread < 1e-5 Ha over the last 5); otherwise the job fails and can be resumed with a larger cap.
                   Calibration on the btpc_pcm log (data/r2/.../btpc_pcm/logs/optimize.log): the forces were inside the
                   thresholds from step 6, yet the energy fell by another 0.66 mHa (0.4 kcal/mol) by step 70 with
                   steps of up to 0.3 A, so a rule on forces alone would stop far too early; W = 10 with 1e-6 Ha never
                   triggers on that descent and does trigger on a true stall (trust radius collapsed, |dE| ~ 1e-7).
result.json item 'opt': converged (all five geomeTRIC criteria), accepted (converged or accepted on negligible forces),
converged_on (which rule), stop_reason.
Start: --start stored (the stored minimum of this case), xray (the case's crystal pair, case['xyz']), or an xyz path /
opt_d3 output directory.  A case without a stored minimum (e.g. btpcz, the P2_1/n polymorph) needs xray or a path; the
RMSD to the stored minimum is then omitted and the RMSD to the start geometry is given instead.
--basis B, --criteria gau|loose, --scf-tol T: for the 6-31G(d) pre-optimisation of dft/optimize.py (loose criteria,
SCF 1e-8), driven by c0_case.py.
Usage:
  python dft/revision/opt_d3.py --job opt_d3_pcm --out <dir> --model pcm [--case btpc] [--start stored|xray|PATH]
         [--seed-hessian auto|none|PATH.npy] [--maxsteps 300] [--max-hours H] [--negligible-stop W]
         [--disp d3bj|none] [--xc b3lyp] [--basis B] [--criteria gau|loose] [--scf-tol T] [--no-fmo] [--smoke]
--smoke: STO-3G, grid 2, Lebedev 17, geomeTRIC's loose criteria (as dft/optimize.py's pre-optimisation), 4 steps;
         the step limit is expected to be reached and is not counted as a failure.
"""
import sys, pathlib, argparse, time, json
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import rev_common as rc          # first: pins the thread pools
import numpy as np
import headline as hl
import charges as ch

CONV = dict(convergence_energy=1e-6, convergence_grms=3e-4, convergence_gmax=4.5e-4,
            convergence_drms=1.2e-3, convergence_dmax=1.8e-3)          # geomeTRIC defaults (Gaussian)
CONV_LOOSE = dict(convergence_energy=1e-6, convergence_grms=1.7e-3, convergence_gmax=2.5e-3,
                  convergence_drms=6.7e-3, convergence_dmax=1.0e-2)     # dft/optimize.py PRE_CONV (GAU_LOOSE)
SEED_FLOOR = 5e-5      # Ha/Bohr^2
SEED_EXT = 0.05        # Ha/Bohr^2
CAP_RULE = dict(window=5, dE=1e-5)   # at a cap: the btpc_pcm / proton_scan.py negligible-forces rule


class OptStop(Exception):
    """Raised from the optimisation callback to end geomeTRIC early (stall detected or a cap reached)."""


# ---------------------------------------------------------------- Hessian helpers (also used by softmode_check.py)
def tr_basis(X, masses=None):
    """Orthonormal basis (3N x 6, or 5 for a linear molecule) of overall translations and rotations at geometry X
    (any length unit), mass-weighted when masses are given."""
    X = np.asarray(X, float); n = len(X)
    w = np.ones(n) if masses is None else np.asarray(masses, float)
    sw = np.sqrt(w)
    c = (w[:, None] * X).sum(0) / w.sum()
    R = X - c
    vecs = []
    for k in range(3):
        v = np.zeros((n, 3)); v[:, k] = sw; vecs.append(v.ravel())
    for k in range(3):
        e = np.zeros(3); e[k] = 1.0
        vecs.append((np.cross(e, R) * sw[:, None]).ravel())
    U, s, _ = np.linalg.svd(np.array(vecs).T, full_matrices=False)
    return U[:, s > 1e-6 * s.max()]


def complement(T):
    """Orthonormal basis of the orthogonal complement of the columns of T."""
    U, _, _ = np.linalg.svd(T, full_matrices=True)
    return U[:, T.shape[1]:]


def hessian_matrix(path):
    """Stored Hessian (natm, natm, 3, 3) or (3N, 3N) -> symmetric (3N, 3N) in Ha/Bohr^2."""
    H = np.load(path)
    if H.ndim == 4:
        n = H.shape[0]
        H = H.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)
    return 0.5 * (H + H.T)


def regularised_seed(H, X_bohr, floor=SEED_FLOOR, ext=SEED_EXT):
    """Positive-definite seed: internal eigenvalues floored, overall translation/rotation set to `ext`."""
    T = tr_basis(X_bohr)
    N = complement(T)
    w, V = np.linalg.eigh(N.T @ H @ N)
    nfl = int((w < floor).sum())
    Hs = (N @ V) @ np.diag(np.maximum(w, floor)) @ (N @ V).T + ext * (T @ T.T)
    return 0.5 * (Hs + Hs.T), {'n_internal_eigenvalues_floored': nfl, 'lowest_internal_eigenvalues': w[:8].tolist(),
                               'floor_Ha_bohr2': floor, 'external_curvature_Ha_bohr2': ext,
                               'n_external': int(T.shape[1])}


def kabsch_rmsd(A, B):
    A = np.asarray(A, float) - np.mean(A, 0); B = np.asarray(B, float) - np.mean(B, 0)
    U, S, Vt = np.linalg.svd(A.T @ B)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return float(np.sqrt(((A @ R - B) ** 2).sum(1).mean()))


def read_traj_best(path):
    """Lowest-energy frame of a trajectory written by these scripts (comment line '... E = <Ha> ...')
    -> (symbols, xyz, number of frames, energy, frame index); falls back to the last frame without energies."""
    lines = open(path).read().splitlines()
    best, i, nfr = None, 0, 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1; continue
        n = int(lines[i]); fr = lines[i + 2:i + 2 + n]
        if len(fr) < n:
            break
        try:
            e = float(lines[i + 1].split('E =')[1].split()[0])
        except (IndexError, ValueError):
            e = float('inf')
        if best is None or e <= best[0] or e == float('inf'):
            best = (e, nfr, fr)
        nfr += 1; i += 2 + n
    if best is None:
        return None, None, 0, None, None
    syms = [l.split()[0] for l in best[2]]
    X = np.array([[float(v) for v in l.split()[1:4]] for l in best[2]])
    return syms, X, nfr, best[0], best[1]


def read_traj_last(path):
    """Last frame of an xyz trajectory -> (symbols, xyz, number of frames)."""
    lines = open(path).read().splitlines()
    frames, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1; continue
        n = int(lines[i]); fr = lines[i + 2:i + 2 + n]
        if len(fr) < n:
            break
        frames.append(fr); i += 2 + n
    if not frames:
        return None, None, 0
    syms = [l.split()[0] for l in frames[-1]]
    X = np.array([[float(v) for v in l.split()[1:4]] for l in frames[-1]])
    return syms, X, len(frames)


def geomeric_norms(g, X_new=None, X_old=None):
    """geomeTRIC-style gradient norms (per-atom vector norms) and step norms (A)."""
    ga = np.linalg.norm(np.asarray(g).reshape(-1, 3), axis=1)
    out = {'grms': float(np.sqrt((ga ** 2).mean())), 'gmax': float(ga.max())}
    if X_new is not None and X_old is not None:
        d = np.linalg.norm(np.asarray(X_new) - np.asarray(X_old), axis=1)
        out.update({'drms_A': float(np.sqrt((d ** 2).mean())), 'dmax_A': float(d.max())})
    return out


def stalled(history, window, grms, gmax, de):
    """Gradient criteria met at each of the last `window` evaluations and their energy spread below `de` Ha."""
    if window <= 0 or len(history) < window:
        return False
    tail = history[-window:]
    if any(h['grms'] > grms or h['gmax'] > gmax for h in tail):
        return False
    e = [h['E_Ha'] for h in tail]
    return max(e) - min(e) < de


def cap_rule(history, grms, gmax, window=CAP_RULE['window'], de=CAP_RULE['dE']):
    """The rule btpc_pcm was accepted on (and proton_scan.py applies at maxsteps): gradient criteria met at the last
    evaluation and an energy spread below `de` over the last `window`."""
    if len(history) < window:
        return False
    tail = history[-window:]
    e = [h['E_Ha'] for h in tail]
    return tail[-1]['grms'] <= grms and tail[-1]['gmax'] <= gmax and max(e) - min(e) < de


# ---------------------------------------------------------------- final analysis at a minimum
def analyse_minimum(job, case, mol, args, lv, disp, solvent, tag='final'):
    """conv_tol 1e-10 single point + FMO, charges, dipole, molden.  Returns the record."""
    frags = rc.fragments(case)
    chk = job.chk_dir / f'{tag}.chk'
    mf = rc.make_mf(mol, xc=args.xc, disp=disp, solvent=solvent, conv_tol=lv['conv_tol'],
                    grid_level=lv['grid_level'], lebedev_order=lv['lebedev_order'], chkfile=chk)
    guess = chk if chk.exists() else (job.chk_dir / 'opt.chk')
    wall = rc.run_scf(mf, chk_guess=guess, log=job.log)
    fr = hl.frontier(mf, frags)
    rec = {'E_Ha': float(mf.e_tot), 'scf_cycles': mf._cycles, 'scf_wall_s': round(wall, 1),
           'settings': rc.settings_of(mf), **rc.runtime_sizes(mf),
           'E_solvation_Ha': hl.solvation_energy(mf), 'frontier': fr,
           'descriptors_fmo_py': hl.fmo_descriptors(fr['HOMO_eV'], fr['LUMO_eV']), 'dipole_debye': hl.dipole_debye(mf)}
    if disp:
        rec['E_dispersion_Ha'] = float(np.asarray(mf.scf_summary.get('dispersion', np.nan)))
        rec['E_without_dispersion_Ha'] = rec['E_Ha'] - rec['E_dispersion_Ha']
    free = ch.free_atom_dms(mol, xc=args.xc)
    pops = ch.full_analysis(mol, mf.make_rdm1(), mf.grids, free, frags, case)
    rec['charges'] = {m: pops[m]['fragments'] for m in ('mulliken', 'lowdin', 'hirshfeld', 'cm5')}
    rec['charges_atoms'] = {m: pops[m]['atoms'] for m in ('mulliken', 'lowdin', 'hirshfeld', 'cm5')}
    if 'proton_H' in pops:
        rec['charge_proton_H'] = pops['proton_H']
    ch.write_molden(mol, mf.mo_coeff, mf.mo_energy, mf.mo_occ, job.out / f'{tag}.molden')
    rec['molden'] = f'{tag}.molden'
    return rec, mf


# ---------------------------------------------------------------- driver
def parser():
    p = rc.add_common_args(argparse.ArgumentParser(description=__doc__.split('\n')[1]))
    p.add_argument('--model', choices=('gas', 'pcm'), required=True)
    p.add_argument('--case', default='btpc')
    p.add_argument('--start', default='stored',
                   help='stored (the B3LYP minimum of this model), xray (the case\'s crystal pair) or an xyz path')
    p.add_argument('--seed-hessian', default='auto', help='auto (stored hessian.npy of this model), none, or a .npy')
    p.add_argument('--seed-floor', type=float, default=SEED_FLOOR)
    p.add_argument('--maxsteps', type=int, default=None, help='total evaluations over all restarts (default 300)')
    p.add_argument('--max-hours', type=float, default=None, help='total optimisation wall hours over all restarts')
    p.add_argument('--negligible-stop', type=int, default=0, metavar='W',
                   help='stop on a stall: gradient criteria met for W evaluations, energy spread < --negligible-de')
    p.add_argument('--negligible-de', type=float, default=1e-6, help='Ha (default 1e-6)')
    p.add_argument('--basis', default=None, help='override the basis (default 6-311++G(d,p); smoke STO-3G)')
    p.add_argument('--criteria', choices=('gau', 'loose'), default=None,
                   help='gau (Gaussian defaults, default) or loose (GAU_LOOSE, dft/optimize.py pre-optimisation)')
    p.add_argument('--scf-tol', type=float, default=None, help='SCF conv_tol during the optimisation')
    p.add_argument('--xc', default=rc.XC)
    p.add_argument('--disp', default='d3bj', help='d3bj (default) or none')
    p.add_argument('--no-fmo', action='store_true', help='skip the analysis at the minimum')
    return p


def run(job, args):
    from pyscf.geomopt import geometric_solver
    lv = hl.level(args.smoke)
    if args.basis:
        lv['basis'] = args.basis
    disp = None if args.disp in (None, '', 'none') else args.disp
    solvent = 'water' if args.model == 'pcm' else None
    loose = args.criteria == 'loose' or (args.criteria is None and args.smoke)
    conv = CONV_LOOSE if loose else CONV
    maxsteps = args.maxsteps or (4 if args.smoke else 300)
    scf_tol = args.scf_tol or (1e-6 if args.smoke else (rc.CONV_TOL_OPT_PCM if solvent else rc.CONV_TOL_OPT_GAS))
    case = rc.load_case(f'{args.case}_{args.model}')
    job.input_file('case', rc.find_case_file(case['name']))
    stored_xyz = hl.resolve_geom(args.model, args.case)[1]
    has_stored = stored_xyz.exists()
    if args.start == 'stored':
        if not has_stored:
            raise SystemExit(f"{case['name']} has no stored minimum ({stored_xyz}); pass --start xray or a path")
        start = stored_xyz
    elif args.start == 'xray':
        start = rc.ROOT / case['xyz']
    else:
        start = pathlib.Path(args.start)
    if not start.is_absolute():
        start = rc.ROOT / start
    if start.is_dir():
        start = start / 'opt.xyz'
    job.input_file('start_geometry', start)
    job.r['settings'] = {**lv, 'xc': args.xc, 'disp': disp, 'model': args.model, 'solvent': solvent,
                         'scf_conv_tol_opt': scf_tol, 'conv_tol_final': lv['conv_tol'], 'geometric_criteria': conv,
                         'criteria_name': 'GAU_LOOSE' if loose else 'GAU (geomeTRIC defaults = Gaussian)',
                         'all_five_required': True, 'maxsteps_total': maxsteps, 'max_hours_total': args.max_hours,
                         'negligible_stop': ({'window': args.negligible_stop, 'dE_Ha': args.negligible_de,
                                              'at_cap': CAP_RULE, 'gradient_norms': 'geomeTRIC per-atom'}
                                             if args.negligible_stop > 0 else None),
                         'optimizer': 'geomeTRIC (TRIC) via pyscf.geomopt.geometric_solver.kernel',
                         'pcm_method': rc.PCM_METHOD if solvent else None}
    job.save()
    traj = job.out / 'opt_traj.xyz'
    ok = True

    # ---------------- optimisation
    if not job.done('opt'):
        geom = start
        prev = job.get('opt') or {}
        history = prev.get('history', [])
        if traj.exists():
            syms, X, nfr, e_b, k_b = read_traj_best(traj)
            if nfr:
                geom = [(s, tuple(x)) for s, x in zip(syms, X)]
                job.log(f'restart: continuing from the lowest-energy frame {k_b} (E = {e_b:.10f}) of {nfr} in {traj.name}')
        mol = rc.mol_for_case(case, geom, basis=lv['basis'], verbose=3 if args.smoke else 4,
                              output=job.out / 'pyscf_opt.log')
        mf = rc.make_mf(mol, xc=args.xc, disp=disp, solvent=solvent, conv_tol=scf_tol,
                        grid_level=lv['grid_level'], lebedev_order=lv['lebedev_order'],
                        chkfile=job.chk_dir / 'opt.chk')
        if (job.chk_dir / 'opt.chk').exists():
            mf.init_guess = 'chkfile'
        kwargs = dict(conv)
        seed_info = None
        if args.seed_hessian != 'none':
            sd = hl.stored_dir(f'{args.case}_{args.model}')
            hp = ((sd / 'hessian.npy' if sd else None) if args.seed_hessian == 'auto'
                  else pathlib.Path(args.seed_hessian))
            if hp and hp.exists():
                H = hessian_matrix(hp)
                if H.shape != (3 * mol.natm, 3 * mol.natm):
                    raise ValueError(f'seed Hessian {hp} has shape {H.shape}, expected {(3 * mol.natm,) * 2}')
                Hs, seed_info = regularised_seed(H, mol.atom_coords(), floor=args.seed_floor)
                seed_info['path'] = str(hp)
                job.input_file('seed_hessian', hp)
                kwargs['hess_data'] = Hs.tolist()       # a list: geomeTRIC tests `if kwargs.get('hess_data')`
                kwargs['frequency'] = False             # no vibrational analysis of the seed
                job.log(f"seed Hessian {hp}: {seed_info['n_internal_eigenvalues_floored']} internal eigenvalues "
                        f"floored at {args.seed_floor:g} Ha/Bohr^2; lowest {np.round(seed_info['lowest_internal_eigenvalues'][:4], 6)}")
            elif args.seed_hessian != 'auto':
                raise FileNotFoundError(hp)
            else:
                job.log('no stored Hessian found; geomeTRIC guess Hessian')
        rec = {'start': str(start), 'seed_hessian': seed_info, 'history': history, 'n_atoms': mol.natm,
               'nao': int(mol.nao), 'settings': rc.settings_of(mf)}
        job.put('opt', rec, status='running')
        tf = open(traj, 'a')
        state = {'X_prev': None, 't': time.time(), 'last': None}
        thr_rms, thr_max = conv['convergence_grms'], conv['convergence_gmax']
        wall_cap = args.max_hours * 3600 if args.max_hours else None

        def wall_used():
            return sum(h.get('wall_s') or 0 for h in history)

        def cb(env):
            m = env['mol']; e = float(env['energy']); g = np.asarray(env['gradients'])
            syms = [m.atom_pure_symbol(i) for i in range(m.natm)]; X = m.atom_coords() * rc.BOHR_A
            nrm = geomeric_norms(g, X, state['X_prev'])
            gc_ = np.asarray(g).ravel()
            dE = e - history[-1]['E_Ha'] if history else None
            history.append({'step': len(history), 'E_Ha': e, 'dE_Ha': dE, **nrm,
                            'gmax_cart': float(np.abs(gc_).max()), 'grms_cart': float(np.sqrt((gc_ ** 2).mean())),
                            'wall_s': round(time.time() - state['t'], 1)})
            state['X_prev'] = X; state['t'] = time.time(); state['last'] = (syms, X)
            tf.write(f'{len(syms)}\nstep {len(history) - 1}  E = {e:.10f} Ha\n' +
                     ''.join(f'{s:2s} {x:16.10f} {y:16.10f} {z:16.10f}\n' for s, (x, y, z) in zip(syms, X)))
            tf.flush()
            job.put('opt', rec, status='running')
            job.log(f"opt step {len(history) - 1}: E = {e:.10f} Ha  dE {dE if dE is None else f'{dE:+.2e}'}  "
                    f"grad rms/max {nrm['grms']:.2e}/{nrm['gmax']:.2e}" +
                    (f"  step rms/max {nrm['drms_A']:.2e}/{nrm['dmax_A']:.2e} A" if 'drms_A' in nrm else ''))
            if args.negligible_stop > 0 and stalled(history, args.negligible_stop, thr_rms, thr_max, args.negligible_de):
                raise OptStop('stall')
            if wall_cap and wall_used() >= wall_cap:
                raise OptStop('max_hours')

        t0 = time.time()
        remaining = maxsteps - len(history)
        stop_reason = None
        if remaining <= 0:
            stop_reason = 'maxsteps'
        elif wall_cap and wall_used() >= wall_cap:
            stop_reason = 'max_hours'
        if stop_reason:
            converged, mol_opt = False, None
            job.log(f'{stop_reason} cap already reached ({len(history)} evaluations, {wall_used() / 3600:.2f} h); '
                    f'no further steps (raise --maxsteps / --max-hours to continue)')
        else:
            try:
                converged, mol_opt = geometric_solver.kernel(mf, maxsteps=remaining, callback=cb, **kwargs)
                if not converged:
                    stop_reason = 'maxsteps'
            except OptStop as ex:
                converged, mol_opt, stop_reason = False, None, str(ex)
        tf.close()
        if mol_opt is not None:
            syms, X = hl.mol_geometry(mol_opt)
        elif state['last'] is not None:
            syms, X = state['last']
        else:
            syms, X, _ = read_traj_last(traj)
        e_written = history[-1]['E_Ha']
        accepted, converged_on = bool(converged), ('geomeTRIC criteria (all five)' if converged else None)
        if not converged and args.negligible_stop > 0:
            last = history[-1]
            if stop_reason == 'stall':
                accepted = True
                converged_on = (f"negligible forces (stall): gradient criteria met at each of the last "
                                f"{args.negligible_stop} evaluations, energy spread < {args.negligible_de:g} Ha; "
                                f"last grad rms/max {last['grms']:.2e}/{last['gmax']:.2e}")
            elif stop_reason in ('maxsteps', 'max_hours') and cap_rule(history, thr_rms, thr_max):
                accepted = True
                converged_on = (f"negligible forces at the {stop_reason} cap: last grad rms/max {last['grms']:.2e}/"
                                f"{last['gmax']:.2e}, energy spread < {CAP_RULE['dE']:g} Ha over the last "
                                f"{CAP_RULE['window']} (the btpc_pcm rule)")
        if accepted and not converged:
            job.log(f'accepted: {converged_on}')
        if not accepted:
            # geomeTRIC returns the last evaluated geometry, which may be a rejected step: keep the lowest-energy one
            s_b, X_b, _, e_b, k_b = read_traj_best(traj)
            if s_b is not None and e_b < history[-1]['E_Ha'] - 1e-8:
                syms, X, e_written = s_b, X_b, e_b
                rec['written_geometry'] = f'lowest-energy evaluation (trajectory frame {k_b}, E = {e_b:.10f} Ha)'
                job.log(f"not converged: opt.xyz is the lowest-energy frame {k_b} (E = {e_b:.10f} Ha)")
        state_word = ('optimised' if converged else 'accepted on negligible forces' if accepted else 'NOT converged')
        rc.write_xyz(job.out / 'opt.xyz', syms, X,
                     f"{case['title']} -- {args.xc}{'-' + disp if disp else ''}/{lv['basis']} "
                     f"{state_word}, E = {e_written:.10f} Ha")
        s_start, X_start = rc.read_xyz(start)
        rec.update({'converged': bool(converged), 'accepted': bool(accepted), 'converged_on': converged_on,
                    'stop_reason': stop_reason if not converged else None,
                    'steps': len(history), 'wall_min': round((time.time() - t0) / 60, 2),
                    'wall_min_total': round(wall_used() / 60, 2),
                    'E_final_Ha': e_written, 'last_step': history[-1],
                    'rmsd_vs_start_A': kabsch_rmsd(X, X_start),
                    'geometry': hl.geometry_metrics(case, syms, X), 'xyz': 'opt.xyz'})
        if has_stored:
            s0, X0 = rc.read_xyz(stored_xyz)
            rec.update({'rmsd_vs_stored_B3LYP_min_A': kabsch_rmsd(X, X0),
                        'geometry_stored_B3LYP_min': hl.geometry_metrics(case, s0, X0)})
        else:
            rec.update({'rmsd_vs_stored_B3LYP_min_A': None, 'geometry_stored_B3LYP_min': None})
        # 'stopped' (not done) when not accepted in a production run: rerunning the same command continues it
        job.put('opt', rec, status='done' if (accepted or args.smoke) else 'stopped')
        job.r['timings']['opt'] = round(time.time() - t0, 1)
        rmsd_txt = (f"RMSD vs stored B3LYP minimum {rec['rmsd_vs_stored_B3LYP_min_A']:.3f} A" if has_stored
                    else f"RMSD vs start {rec['rmsd_vs_start_A']:.3f} A")
        job.log(f"optimisation {'CONVERGED' if converged else 'ACCEPTED (' + str(converged_on) + ')' if accepted else 'NOT converged (' + str(stop_reason) + ')'} "
                f"after {len(history)} evaluations, E(opt.xyz) = {e_written:.10f} Ha, {rmsd_txt}")
        mol.stdout.close()
    rec = job.get('opt')
    if not rec.get('accepted', rec.get('converged')) and not args.smoke:
        ok = False
        job.log(f"optimisation not accepted ({rec.get('stop_reason')}); rerun the same command (with a larger "
                f"--maxsteps / --max-hours if a cap was reached) to continue from {traj.name}")

    # ---------------- analysis at the minimum (only at an accepted geometry; redone if opt.xyz changed since)
    opt_sha = rc.sha256_file(job.out / 'opt.xyz') if (job.out / 'opt.xyz').exists() else None
    if (not args.no_fmo and (ok or args.smoke)
            and (not job.done('final') or (job.get('final') or {}).get('opt_xyz_sha256', opt_sha) != opt_sha)):
        syms, X = rc.read_xyz(job.out / 'opt.xyz')
        mol = rc.mol_for_case(case, job.out / 'opt.xyz', basis=lv['basis'], verbose=3 if args.smoke else 4,
                              output=job.out / 'pyscf_final.log')
        fin, mf = analyse_minimum(job, case, mol, args, lv, disp, solvent)
        fin['geometry'] = hl.geometry_metrics(case, syms, X)
        fin['from_converged_optimisation'] = bool(rec.get('converged'))
        fin['from_accepted_optimisation'] = bool(rec.get('accepted', rec.get('converged')))
        fin['opt_xyz_sha256'] = opt_sha
        job.put('final', fin)
        fr = fin['frontier']
        job.log(f"final: E = {fin['E_Ha']:.10f} Ha; HOMO {fr['HOMO_eV']:.4f} eV ({fr['HOMO_on']}), LUMO "
                f"{fr['LUMO_eV']:.4f} eV ({fr['LUMO_on']}), gap {fr['gap_eV']:.4f} eV; picrate charge Hirshfeld "
                f"{fin['charges']['hirshfeld']['anion']:+.3f}, CM5 {fin['charges']['cm5']['anion']:+.3f} e")
        mol.stdout.close()
    f = job.get('final') or {}
    job.r['summary'] = {'converged': rec.get('converged'), 'accepted': rec.get('accepted', rec.get('converged')),
                        'converged_on': rec.get('converged_on'), 'stop_reason': rec.get('stop_reason'),
                        'steps': rec.get('steps'), 'E_opt_Ha': rec.get('E_final_Ha'),
                        'E_final_sp_Ha': f.get('E_Ha'), 'E_dispersion_Ha': f.get('E_dispersion_Ha'),
                        'HOMO_eV': (f.get('frontier') or {}).get('HOMO_eV'),
                        'LUMO_eV': (f.get('frontier') or {}).get('LUMO_eV'),
                        'LUMO_on': (f.get('frontier') or {}).get('LUMO_on'),
                        'geometry': rec.get('geometry'), 'rmsd_vs_stored_B3LYP_min_A': rec.get('rmsd_vs_stored_B3LYP_min_A'),
                        'rmsd_vs_start_A': rec.get('rmsd_vs_start_A')}
    job.save()
    return ok


if __name__ == '__main__':
    rc.install_sigterm_handler()
    rc.main_wrapper(parser, run, __file__)
