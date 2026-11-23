"""
Validation of the revision environment against the published BTPC numbers.

Single points at the two stored minima with the level of theory of record (dft/revision/rev_common.py:
B3LYP/6-311++G(d,p), density fitting, grid level 4, IEF-PCM water with Lebedev order 29), compared with the stored
energies.  Pass: |dE| <= 1e-6 Ha against the stored single point (fmo.json E_total_Ha, the conv_tol 1e-10 single
point at opt.xyz; opt_summary.json E_final_Ha agrees with it to < 1e-9 Ha).  Also records the wall time of one SCF
(from the default guess) and of one analytic gradient, the frontier-orbital energies and the Mulliken fragment
composition against fmo.json, and checks that rev_common.make_mf sets exactly what dft/common.py set.

Usage:
  python dft/revision/validate.py --job val --out <dir> [--models gas,pcm] [--no-grad] [--guess minao|chk] [--smoke]
--smoke: STO-3G basis, grid level 2, Lebedev 17, conv_tol 1e-6 -- exercises every code path in about a minute;
         its energies are not comparable and no pass/fail is applied.
"""
import sys, pathlib, argparse, json, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import rev_common as rc          # first: pins the thread pools
import numpy as np

MODELS = {'gas': 'btpc_gas', 'pcm': 'btpc_pcm'}
TOL_HA = 1e-6
TASK_REF = {'gas': -1317.20823094, 'pcm': -1317.25368857}   # the values quoted in the revision plan (8 dp)


def parser():
    p = rc.add_common_args(argparse.ArgumentParser(description=__doc__.split('\n')[1]))
    p.add_argument('--models', default='gas,pcm', help='comma list of gas, pcm')
    p.add_argument('--no-grad', action='store_true', help='skip the gradient timing')
    p.add_argument('--guess', choices=('minao', 'chk'), default='minao',
                   help="initial density: PySCF default (minao) for a true per-SCF timing, or the stored final.chk")
    return p


def stored_refs(case_name):
    d = rc.ROOT / 'molecules' / 'btpc' / 'results' / case_name
    fmo = json.load(open(d / 'fmo.json'))
    opt = json.load(open(d / 'opt_summary.json'))
    return fmo, opt


def run(job, args):
    smoke = args.smoke
    basis = 'sto-3g' if smoke else rc.BASIS
    grid = 2 if smoke else rc.GRID_LEVEL
    leb = 17 if smoke else rc.PCM_LEBEDEV_ORDER
    conv = 1e-6 if smoke else rc.CONV_TOL_SP
    job.r['settings'] = {'xc': rc.XC, 'basis': basis, 'grid_level': grid, 'pcm_lebedev_order': leb,
                         'pcm_method': rc.PCM_METHOD, 'eps_water': rc.SOLVENT_EPS['water'], 'conv_tol': conv,
                         'max_cycle': rc.MAX_CYCLE, 'density_fit': True, 'auxbasis': 'PySCF default',
                         'tolerance_Ha': TOL_HA, 'guess': args.guess}
    job.save()
    all_ok = True
    for model in [m.strip() for m in args.models.split(',') if m.strip()]:
        case_name = MODELS[model]
        if job.done(model):
            it = job.get(model)
            job.log(f"{model}: already done (E = {it['E_Ha']:.10f}, dE = {it.get('dE_Ha')}); skipped")
            all_ok &= it.get('pass') is not False
            continue
        case = rc.load_case(case_name)
        geom = rc.ROOT / rc.STORED[case_name]['xyz']
        job.input_file(f'{model}_geometry', geom)
        job.input_file(f'{model}_case', rc.ROOT / case['case_file'])
        fmo, opt = stored_refs(case_name)
        solvent = case.get('solvent')
        mol = rc.mol_for_case(case, geom, basis=basis, verbose=4 if not smoke else 3,
                              output=job.out / f'pyscf_{model}.log')
        rec = {'case': case_name, 'geometry': str(geom.relative_to(rc.ROOT)), 'solvent': solvent,
               'nao': int(mol.nao), 'natm': int(mol.natm)}

        # 1. same settings as dft/common.py (as launched by run_case.sh)?
        diffs = rc.compare_with_common(mol, solvent)
        rec['settings_vs_common'] = diffs or 'identical'
        job.log(f"{model}: settings vs dft/common.py (grid 4 / Lebedev 29): {diffs or 'identical'}")
        if diffs and not smoke:
            raise RuntimeError(f'rev_common.make_mf differs from dft/common.py: {diffs}')

        # 2. SCF
        chk = job.chk_dir / f'{model}.chk'
        mf = rc.make_mf(mol, solvent=solvent, conv_tol=conv, grid_level=grid, lebedev_order=leb, chkfile=chk)
        rec['settings'] = rc.settings_of(mf)
        guess = None; guess_label = 'minao (PySCF default)'
        if args.guess == 'chk':
            src = rc.ROOT / 'data' / 'r2' / 'molecules' / 'btpc' / 'results' / case_name / 'final.chk'
            guess, guess_label = src, f'stored {src.relative_to(rc.ROOT)}'
        elif chk.exists():
            guess, guess_label = chk, 'own checkpoint from an interrupted run'
        rec['guess'] = guess_label
        job.log(f"{model}: SCF {mol.nao} basis functions, {mol.natm} atoms, solvent {solvent}, guess {guess_label} ...")
        job.put(model, rec, status='running')
        wall = rc.run_scf(mf, chk_guess=guess, log=job.log)
        e = float(mf.e_tot)
        rec.update(rc.runtime_sizes(mf))
        rec.update({'E_Ha': e, 'scf_cycles': mf._cycles, 'scf_wall_s': round(wall, 1),
                    'scf_wall_per_cycle_s': round(wall / max(1, mf._cycles), 2)})
        if hasattr(mf, 'with_solvent') and getattr(mf.with_solvent, 'e', None) is not None:
            rec['E_solvation_Ha'] = float(mf.with_solvent.e)
        job.r['timings'][f'{model}_scf'] = round(wall, 1)
        job.log(f"{model}: E = {e:.10f} Ha in {wall:.0f} s ({mf._cycles} cycles); sizes {rc.runtime_sizes(mf)}")

        # 3. comparison with the stored numbers
        ref = float(fmo['E_total_Ha'])
        rec.update({'E_ref_fmo_Ha': ref, 'E_ref_opt_Ha': float(opt['E_final_Ha']), 'E_ref_plan_Ha': TASK_REF[model],
                    'dE_Ha': e - ref, 'dE_vs_opt_Ha': e - float(opt['E_final_Ha']),
                    'dE_vs_plan_Ha': e - TASK_REF[model]})
        n = int((mf.mo_occ > 0).sum())
        H, L = mf.mo_energy[n - 1] * rc.HARTREE_EV, mf.mo_energy[n] * rc.HARTREE_EV
        rec.update({'HOMO_eV': float(H), 'LUMO_eV': float(L), 'gap_eV': float(L - H),
                    'dHOMO_eV': float(H - fmo['HOMO_eV']), 'dLUMO_eV': float(L - fmo['LUMO_eV'])})
        frags = rc.fragments(case)
        comp = rc.fragment_composition(mf, frags)
        rec['fragment_composition'] = comp
        # stored composition (fmo.py, Mulliken, 0.1 %) keyed by the case's fragment names
        stored = {k: {rc.frag_key(f): v for f, v in d.items()}
                  for k, d in (fmo.get('fragment_composition') or {}).items()}
        rec['fragment_composition_stored'] = stored
        if stored:
            rec['fragment_composition_max_diff_pct'] = max(
                abs(comp[k]['mulliken'][f] - stored[k][f]) for k in stored for f in stored[k])
        if smoke:
            rec['pass'] = None
        else:
            rec['pass'] = bool(abs(e - ref) <= TOL_HA)
            all_ok &= rec['pass']
        job.log(f"{model}: dE = {e - ref:+.3e} Ha vs stored {ref:.10f} (fmo.json); vs plan value "
                f"{TASK_REF[model]:.8f}: {e - TASK_REF[model]:+.2e}; HOMO {H:.4f} (stored {fmo['HOMO_eV']:.4f}) "
                f"LUMO {L:.4f} (stored {fmo['LUMO_eV']:.4f}) eV; LUMO composition {comp['LUMO']['mulliken']} -> "
                f"{'PASS' if rec['pass'] else 'FAIL' if rec['pass'] is False else 'smoke (no criterion)'}")

        # 4. gradient timing
        if not args.no_grad:
            t = time.time()
            g = mf.nuc_grad_method().kernel()
            gw = time.time() - t
            rec.update({'grad_wall_s': round(gw, 1), 'grad_max_Ha_bohr': float(abs(g).max()),
                        'grad_rms_Ha_bohr': float(np.sqrt((g ** 2).mean()))})
            job.r['timings'][f'{model}_grad'] = round(gw, 1)
            job.log(f"{model}: gradient {gw:.0f} s, max |g| {abs(g).max():.2e}, rms {np.sqrt((g**2).mean()):.2e} Ha/Bohr")
        rec['peak_rss_mb'] = round(rc.peak_rss_mb(), 1)
        job.put(model, rec)
        mol.stdout.close()
    job.r['summary'] = {m: {k: job.get(m).get(k) for k in ('E_Ha', 'dE_Ha', 'pass', 'scf_wall_s', 'scf_cycles',
                                                            'grad_wall_s', 'peak_rss_mb')}
                        for m in MODELS if job.get(m)}
    job.save()
    return all_ok


if __name__ == '__main__':
    rc.install_sigterm_handler()
    rc.main_wrapper(parser, run, __file__)
