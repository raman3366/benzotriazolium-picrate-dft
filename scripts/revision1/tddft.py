"""
TD-DFT excited states of the BTPC ion pair and of its free components (plan section 5, item 7).

Systems (--systems, comma list):
  pair-gas      ion pair at the stored gas-phase B3LYP minimum, gas phase
  pair-pcm      ion pair at the stored IEF-PCM minimum, IEF-PCM water
  picrate-pcm   free picrate anion (the ion pair's anion atoms, charge -1), optimised in IEF-PCM water
  bta-pcm       neutral 1H-benzotriazole (the cation's atoms without the H-bonded proton, case['proton']['H']),
                optimised in IEF-PCM water
  (picrate-gas / bta-gas: the same in the gas phase.)
Functionals (--xc): b3lyp and camb3lyp (CAM-B3LYP) by default; the free species are optimised once with B3LYP
(geomeTRIC, the five default criteria) and both functionals are applied at that geometry, as for the pair.
Response: the Tamm-Dancoff approximation (TDA) by default, --rpa for full TD-DFT (RPA); --nstates singlets (default 20).
In PySCF the response is dominated by the XC-kernel contraction on the DFT grid for every trial vector; full RPA needed
3.3x the trial vectors of TDA for the same 3 roots of the ion pair (msi, 20 Nov 2026: 407 vs 124), so TDA is the
production default and RPA the option (e.g. for oscillator strengths of the bright states).
Solvent: PySCF 2.14 linear-response IEF-PCM for TD-DFT, non-equilibrium by default (the fast part of the continuum
responds with the optical dielectric constant eps_opt = 1.78, fixed in pyscf/solvent/_attach_solvent.py; the ground
state is fully equilibrated with eps 78.3553); --equilibrium uses the static eps for the response instead.
Per state: excitation energy (eV, nm), oscillator strength (length gauge), the leading orbital pairs, and for the ion
pair the hole and electron fractions on cation and picrate and the charge-transfer number from the Omega matrix of the
transition density X+Y (X for TDA) with a Loewdin partition (Plasser & Lischka, J. Chem. Theory Comput. 8, 2777 (2012)):
  T = C_occ (X+Y) C_vir^T,  T_L = S^1/2 T S^1/2,  Omega_AB = sum_{mu in A, nu in B} T_L,mu nu^2,
  hole on A = sum_B Omega_AB / Omega,  electron on B = sum_A Omega_AB / Omega,  CT = sum_{A != B} Omega_AB / Omega.
Ground state per functional: the frontier-orbital analysis of headline.py (LUMO location with that functional).

Usage:
  python dft/revision/tddft.py --job td_pair --out <dir> [--systems pair-gas,pair-pcm,picrate-pcm,bta-pcm]
         [--xc b3lyp,camb3lyp] [--nstates 20] [--rpa] [--equilibrium] [--case btpc] [--smoke]
--smoke: STO-3G, grid 2, Lebedev 17, 3 states, free-species optimisation 3 steps with loose criteria.
"""
import sys, pathlib, argparse, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import rev_common as rc          # first: pins the thread pools
import numpy as np
import headline as hl
from opt_d3 import CONV, CONV_LOOSE, geomeric_norms

EV_NM = 1239.841984


def parser():
    p = rc.add_common_args(argparse.ArgumentParser(description=__doc__.split('\n')[1]))
    p.add_argument('--systems', default='pair-gas,pair-pcm,picrate-pcm,bta-pcm')
    p.add_argument('--xc', default='b3lyp,camb3lyp')
    p.add_argument('--nstates', type=int, default=None)
    p.add_argument('--rpa', action='store_true', help='full TD-DFT (RPA) instead of the default TDA')
    p.add_argument('--tda', action='store_true', help='TDA (the default; accepted for explicitness)')
    p.add_argument('--equilibrium', action='store_true', help='equilibrium LR-PCM response (default non-equilibrium)')
    p.add_argument('--case', default='btpc')
    p.add_argument('--td-conv', type=float, default=1e-6, help='TD-DFT eigenvalue convergence, Ha (default 1e-6)')
    return p


def orbital_label(k, nocc):
    return ('HOMO' + (f'-{nocc - 1 - k}' if k < nocc - 1 else '')) if k < nocc else ('LUMO' + (f'+{k - nocc}' if k > nocc else ''))


def analyse_states(td, mf, frags, rpa):
    mol = mf.mol
    occ = np.asarray(mf.mo_occ); C = np.asarray(mf.mo_coeff)
    oi = np.where(occ > 0)[0]; vi = np.where(occ == 0)[0]; nocc = len(oi)
    Co, Cv = C[:, oi], C[:, vi]
    S = mol.intor_symmetric('int1e_ovlp'); Sh = rc.lowdin_half(S)
    sl = mol.aoslice_by_atom()
    fidx = {f: np.concatenate([np.arange(sl[a][2], sl[a][3]) for a in atoms]) for f, atoms in frags.items()}
    e_eV = np.asarray(td.e) * rc.HARTREE_EV
    fosc = np.asarray(td.oscillator_strength(gauge='length'))
    conv = np.atleast_1d(getattr(td, 'converged', [True] * len(e_eV)))
    states = []
    for n, (x, y) in enumerate(td.xy):
        x = np.asarray(x).reshape(nocc, -1)
        y = np.zeros_like(x) if np.isscalar(y) or np.ndim(y) == 0 else np.asarray(y).reshape(nocc, -1)
        w = x ** 2 - y ** 2
        w = w / w.sum()
        top = np.argsort(w.ravel())[::-1][:4]
        pairs = [{'from': orbital_label(int(oi[i // w.shape[1]]), nocc),
                  'to': orbital_label(int(vi[i % w.shape[1]]), nocc), 'weight_pct': round(100 * float(w.ravel()[i]), 1)}
                 for i in top if w.ravel()[i] > 0.02]
        st = {'state': n + 1, 'E_eV': float(e_eV[n]), 'lambda_nm': float(EV_NM / e_eV[n]), 'f': float(fosc[n]),
              'converged': bool(conv[n]) if n < len(conv) else None, 'transitions': pairs}
        if len(frags) > 1:
            T = Co @ (x + y if rpa else x) @ Cv.T
            TL = Sh @ T @ Sh
            Om = {(a, b): float((TL[np.ix_(fidx[a], fidx[b])] ** 2).sum()) for a in frags for b in frags}
            tot = sum(Om.values())
            st['omega'] = {f'{a}->{b}': v / tot for (a, b), v in Om.items()}
            st['hole'] = {a: sum(Om[(a, b)] for b in frags) / tot for a in frags}
            st['electron'] = {b: sum(Om[(a, b)] for a in frags) / tot for b in frags}
            st['CT_number'] = sum(v for (a, b), v in Om.items() if a != b) / tot
        states.append(st)
    return states


def free_species(case, X, syms, which):
    """Atom list, charge, symbols and coordinates of the free picrate anion or neutral benzotriazole."""
    fr = rc.fragments(case)
    if which == 'picrate':
        idx = list(fr['anion']); charge = -1
    elif which == 'bta':
        H = case['proton']['H']
        idx = [i for i in fr['cation'] if i != H]; charge = 0
    else:
        raise ValueError(which)
    return idx, charge, [syms[i] for i in idx], X[idx]


def run(job, args):
    from pyscf.geomopt import geometric_solver
    lv = hl.level(args.smoke)
    nstates = args.nstates or (3 if args.smoke else 20)
    xcs = [x.strip() for x in args.xc.split(',') if x.strip()]
    systems = [s.strip() for s in args.systems.split(',') if s.strip()]
    if args.rpa and args.tda:
        raise SystemExit('--rpa and --tda are exclusive')
    rpa = bool(args.rpa)
    job.r['settings'] = {**lv, 'xc': xcs, 'nstates': nstates, 'response': 'full TD-DFT (RPA)' if rpa else 'TDA',
                         'td_conv_tol_Ha': args.td_conv, 'pcm_method': rc.PCM_METHOD,
                         'pcm_ground_state_eps': rc.SOLVENT_EPS['water'],
                         'pcm_response': 'equilibrium (static eps)' if args.equilibrium else
                                         'non-equilibrium, optical eps 1.78 (PySCF 2.14 default)',
                         'free_species_optimisation': 'B3LYP, geomeTRIC default criteria' if not args.smoke else 'smoke',
                         'ct_analysis': 'Omega matrix of X+Y (X for TDA), Loewdin partition'}
    job.save()
    ok = True
    case_gas = rc.load_case(f'{args.case}_gas')
    frags_pair = rc.fragments(case_gas)
    fdir = job.out / 'free'; fdir.mkdir(exist_ok=True)

    for sysname in systems:
        kind, model = sysname.rsplit('-', 1)
        solvent = 'water' if model == 'pcm' else None
        case = rc.load_case(f'{args.case}_{model}')
        gpath = hl.resolve_geom(model, args.case)[1]
        job.input_file(f'geometry_{model}', gpath)
        if kind == 'pair':
            geom, charge, frags = gpath, 0, frags_pair
        else:
            okey = f'opt/{sysname}'
            if not job.done(okey):
                syms, X = rc.read_xyz(gpath)
                idx, charge, s2, X2 = free_species(case, X, syms, kind)
                mol = rc.build_mol([(s, tuple(x)) for s, x in zip(s2, X2)], charge=charge, spin=0, basis=lv['basis'],
                                   verbose=3, output=job.out / f'pyscf_opt_{sysname}.log')
                mf = rc.make_mf(mol, xc=rc.XC, solvent=solvent, conv_tol=1e-6 if args.smoke else
                                (rc.CONV_TOL_OPT_PCM if solvent else rc.CONV_TOL_OPT_GAS),
                                grid_level=lv['grid_level'], lebedev_order=lv['lebedev_order'],
                                chkfile=job.chk_dir / f'opt_{sysname}.chk')
                hist = []

                def cb(env):
                    hist.append({'E_Ha': float(env['energy']), **geomeric_norms(env['gradients'])})
                    job.log(f"{okey} step {len(hist) - 1}: E = {hist[-1]['E_Ha']:.10f}  grad rms/max "
                            f"{hist[-1]['grms']:.1e}/{hist[-1]['gmax']:.1e}")
                t0 = time.time()
                job.put(okey, {'atoms_from_pair': idx, 'charge': charge, 'history': hist}, status='running')
                conv, mo = geometric_solver.kernel(mf, maxsteps=3 if args.smoke else 100, callback=cb,
                                                   **(CONV_LOOSE if args.smoke else CONV))
                s3, X3 = hl.mol_geometry(mo)
                rc.write_xyz(fdir / f'{sysname}.xyz', s3, X3, f'{sysname} B3LYP/{lv["basis"]} {"optimised" if conv else "NOT converged"}')
                job.put(okey, {'atoms_from_pair': idx, 'charge': charge, 'converged': bool(conv), 'steps': len(hist),
                               'E_Ha': hist[-1]['E_Ha'], 'xyz': f'free/{sysname}.xyz', 'history': hist,
                               'wall_min': round((time.time() - t0) / 60, 2)})
                job.log(f"{okey}: {'converged' if conv else 'NOT converged'} in {len(hist)} steps")
                ok &= bool(conv) or args.smoke
                mol.stdout.close()
            o = job.get(okey)
            geom, charge, frags = job.out / o['xyz'], o['charge'], {}
        for xc in xcs:
            key = f'{sysname}/{xc}'
            if job.done(key):
                continue
            mol = rc.build_mol(geom, charge=charge, spin=0, basis=lv['basis'],
                               basis_overrides=case.get('basis_overrides') if kind == 'pair' else None,
                               verbose=4 if not args.smoke else 3, output=job.out / f'pyscf_{sysname}_{xc}.log')
            chk = job.chk_dir / f'{sysname}_{xc}.chk'
            mf = rc.make_mf(mol, xc=xc, solvent=solvent, conv_tol=lv['conv_tol'], grid_level=lv['grid_level'],
                            lebedev_order=lv['lebedev_order'], chkfile=chk)
            rec = {'system': sysname, 'xc': xc, 'geometry': str(geom), 'charge': charge, 'settings': rc.settings_of(mf)}
            job.put(key, rec, status='running')
            t = time.time()
            rc.run_scf(mf, chk_guess=chk, log=job.log)
            t_scf = time.time() - t
            rec.update({'E_Ha': float(mf.e_tot), 'scf_cycles': mf._cycles, 'scf_wall_s': round(t_scf, 1)})
            if frags:
                fr = hl.frontier(mf, frags)
                rec['frontier'] = fr
            else:
                n = int((mf.mo_occ > 0).sum()); e = mf.mo_energy * rc.HARTREE_EV
                rec['frontier'] = {'HOMO_eV': float(e[n - 1]), 'LUMO_eV': float(e[n]), 'gap_eV': float(e[n] - e[n - 1])}
            td = (mf.TDDFT if rpa else mf.TDA)(**({'equilibrium_solvation': True} if (solvent and args.equilibrium)
                                                       else {}))
            td.nstates = nstates
            td.conv_tol = args.td_conv
            if solvent:
                rec['td_solvent'] = {'eps_response': float(td.with_solvent.eps),
                                     'equilibrium_solvation': bool(td.with_solvent.equilibrium_solvation)}
            t = time.time()
            td.kernel()
            t_td = time.time() - t
            states = analyse_states(td, mf, frags, rpa)
            rec.update({'td_class': type(td).__name__, 'td_wall_s': round(t_td, 1), 'states': states,
                        'all_converged': bool(np.all(np.atleast_1d(td.converged)))})
            bright = [s for s in states if s['f'] > 0.01]
            rec['lowest_bright'] = bright[0] if bright else None
            rec['strongest'] = max(states, key=lambda s: s['f'])
            job.put(key, rec)
            job.r['timings'][key] = round(t_scf + t_td, 1)
            lb = rec['lowest_bright']
            job.log(f"{key}: {len(states)} states in {t_td:.0f} s ({type(td).__name__}); S1 {states[0]['E_eV']:.3f} eV "
                    f"({states[0]['lambda_nm']:.0f} nm, f {states[0]['f']:.4f}" +
                    (f", CT {states[0]['CT_number']:.2f}" if 'CT_number' in states[0] else '') + ")" +
                    (f"; lowest bright S{lb['state']} {lb['E_eV']:.3f} eV ({lb['lambda_nm']:.0f} nm, f {lb['f']:.3f})" if lb else ''))
            ok &= rec['all_converged'] or args.smoke
            mol.stdout.close()
    job.r['summary'] = {k: {'S1_eV': v['states'][0]['E_eV'], 'S1_nm': v['states'][0]['lambda_nm'],
                            'S1_f': v['states'][0]['f'], 'S1_CT': v['states'][0].get('CT_number'),
                            'lowest_bright': v.get('lowest_bright'), 'LUMO_on': v['frontier'].get('LUMO_on'),
                            'all_converged': v.get('all_converged')}
                        for k, v in job.r['items'].items() if v.get('status') == 'done' and 'states' in v}
    job.save()
    return ok


if __name__ == '__main__':
    rc.install_sigterm_handler()
    rc.main_wrapper(parser, run, __file__)
