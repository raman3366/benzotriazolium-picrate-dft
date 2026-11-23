"""
Optimisation of the BTPC ion pair started from the neutral (co-crystal) form (revision item 2).

The transferring proton (case 'proton' block: H4N, donor N4, acceptor O1 for btpc) is moved onto the acceptor,
O1-H 0.98 A along the N4...O1 line, everything else left as in the start structure (default: the stored
minimum of the same model; --from crystal uses the case's X-ray cluster).  The structure is then optimised
without constraints at the level of record (rev_common: B3LYP/6-311++G(d,p), density fitting, grid level 4;
IEF-PCM water with Lebedev order 29 for the pcm model), with Gaussian's default criteria as for every
published geometry.  Reported per model: whether the proton returns to N4 (salt) or stays on O1 (neutral
minimum), the step at which it crossed the N4...O1 midpoint, the final distances, the tight single-point
energy and its difference from the stored ionic minimum, the heavy-atom RMSD to that minimum, charges and
frontier orbitals.

Usage:
  python dft/revision/neutral_start.py --job ID --out DIR [--system btpc] [--models gas,pcm]
         [--from opt|crystal] [--start-gas XYZ] [--start-pcm XYZ] [--r-oh 0.98] [--xc b3lyp] [--disp d3bj]
         [--maxsteps 300] [--smoke]
Cases are <system>_gas and <system>_pcm (molecules/*/cases/).  Outputs in DIR: result.json (items 'gas', 'pcm'),
log.txt, start_<model>.xyz, points/<model>.xyz (final), traj/<model>.xyz (every step), pyscf/<model>.log, chk/.
Restart: a finished model is skipped; an interrupted one resumes from the last geometry of its trajectory.
--smoke: STO-3G, grid 2, Lebedev 17, loose criteria, 4 steps; not for numbers.
"""
import sys, pathlib, argparse
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import rev_common as rc          # first: pins the thread pools
import numpy as np
import proton_scan as ps


def parser():
    p = rc.add_common_args(argparse.ArgumentParser(description=__doc__.split('\n')[1]))
    p.add_argument('--system', default='btpc', help='case prefix: cases <system>_gas / <system>_pcm (default btpc)')
    p.add_argument('--models', default='gas,pcm', help='comma list of gas, pcm')
    p.add_argument('--from', dest='start_from', choices=('opt', 'crystal'), default='opt',
                   help='heavy-atom frame: the stored minimum of the model (default) or the X-ray cluster')
    p.add_argument('--start-gas', default=None, help='explicit start geometry for the gas model (XYZ)')
    p.add_argument('--start-pcm', default=None, help='explicit start geometry for the pcm model (XYZ)')
    p.add_argument('--r-oh', type=float, default=0.98, help='O1-H distance of the neutral start (A)')
    p.add_argument('--xc', default=rc.XC, help='functional (default b3lyp)')
    p.add_argument('--disp', default=None, help='dispersion correction, e.g. d3bj (default none, as published)')
    p.add_argument('--maxsteps', type=int, default=None, help='optimisation steps (default 300; smoke 4)')
    return p


def run(job, args):
    if args.maxsteps is None:
        args.maxsteps = 4 if args.smoke else 300
    job.r['script_sha256']['proton_scan.py'] = rc.sha256_file(pathlib.Path(ps.__file__).resolve())
    models = [m.strip() for m in args.models.split(',') if m.strip()]
    job.r['settings'].update({'system': args.system, 'models': models, 'start_from': args.start_from,
                              'neutral_r_OH_A': args.r_oh})
    job.save()
    if args.smoke:
        job.log('SMOKE settings: STO-3G, grid 2, Lebedev 17, loose criteria -- not for production numbers')
    for model in models:
        if model not in ('gas', 'pcm'):
            raise SystemExit(f'unknown model {model!r}')
        if job.done(model):
            it = job.get(model)
            job.log(f"{model}: done ({it['outcome']}, E = {it['E_Ha']:.8f}); skipped")
            continue
        case = rc.load_case(f'{args.system}_{model}')
        L = ps.level_of_theory(args, case.get('solvent'))
        job.r['settings'][model] = L
        H, N, O = ps.proton_idx(case)
        labs = ps.proton_labels(case)
        explicit = getattr(args, f'start_{model}')
        if explicit:
            src = pathlib.Path(explicit)
        elif args.start_from == 'opt':
            src = ps.stored_minimum(case)
            if src is None:
                raise SystemExit(f"{case['name']} has no stored opt.xyz; use --from crystal or --start-{model}")
        else:
            src = rc.ROOT / case['xyz']
        job.input_file(f'{model}_case', rc.ROOT / case['case_file'])
        job.input_file(f'{model}_start_source', src)
        syms, x_src = rc.read_xyz(src)
        x_start = ps.place_H(x_src, H, O, args.r_oh, toward=N)
        start_file = job.out / f'start_{model}.xyz'
        rc.write_xyz(start_file, syms, x_start,
                     f"neutral start: {labs['H']} on {labs['O']} ({args.r_oh:.2f} A along {labs['O']}->{labs['N']}), "
                     f"frame from {src.name}")
        m_src, m_start = ps.measure(x_src, H, N, O), ps.measure(x_start, H, N, O)
        job.log(f"{model}: {case['name']} ({L['solvent'] or 'gas phase'}); source {src.name}: "
                f"r({labs['N']}-{labs['H']}) {m_src['r_NH_A']:.3f}, r({labs['H']}...{labs['O']}) {m_src['r_HO_A']:.3f}; "
                f"neutral start r({labs['N']}...{labs['H']}) {m_start['r_NH_A']:.3f}, "
                f"r({labs['O']}-{labs['H']}) {m_start['r_HO_A']:.3f} A")
        rec = ps.relax(job, model, case, L, syms, x_start, constraint=None, guess_chk=None,
                       extra={'case': case['name'], 'start_source': str(src.relative_to(rc.ROOT))
                              if str(src.resolve()).startswith(str(rc.ROOT)) else str(src),
                              'start_file': start_file.name, 'source_geometry': m_src, 'labels': labs})
        # outcome
        on_N = rec['r_NH_A'] < rec['r_HO_A']
        rec['proton_on'] = labs['N'] if on_N else labs['O']
        if rec.get('converged'):
            rec['outcome'] = (f"proton returned to {labs['N']} (salt)" if on_N
                              else f"proton stays on {labs['O']} (neutral form is a stationary point here)")
        else:
            rec['outcome'] = (f"proton on {labs['N'] if on_N else labs['O']} at the last step "
                              f"(optimisation NOT converged; smoke)")
        cross = next((h['step'] for h in rec['history'] if h['r_NH_A'] < h['r_HO_A']), None)
        rec['returned_at_step'] = cross
        rec['history_r'] = [(h['step'], h['r_NH_A'], h['r_HO_A']) for h in rec['history']]
        # comparison with the stored ionic minimum of the same model
        smin = ps.stored_minimum(case)
        e_ref = rc.STORED.get(case['name'], {}).get('E_Ha')
        if smin is not None:
            s2, x_min = rc.read_xyz(smin)
            heavy = [i for i, s in enumerate(syms) if s != 'H']
            x_fin = ps.read_frames(job.out / rec['xyz'])[0][1]
            rec['heavy_atom_rmsd_to_stored_min_A'] = round(ps.kabsch_rmsd(x_min, x_fin, heavy)[0], 4)
            rec['all_atom_rmsd_to_stored_min_A'] = round(ps.kabsch_rmsd(x_min, x_fin)[0], 4)
            rec['stored_min_geometry'] = ps.measure(x_min, H, N, O)
        if e_ref is not None and not args.smoke and not args.disp and args.xc == rc.XC:
            rec['E_stored_ionic_min_Ha'] = e_ref
            rec['dE_vs_stored_min_kcal'] = round((rec['E_Ha'] - e_ref) * rc.HARTREE_KCAL, 4)
        if not on_N and rec.get('converged'):
            rec['note'] = ('the proton stayed on the acceptor: confirm with a Hessian (or displacements along the '
                           'proton coordinate) before calling it a minimum')
        job.put(model, rec)
        job.log(f"{model}: {rec['outcome']}; crossed back at step {cross}; final r({labs['N']}-{labs['H']}) "
                f"{rec['r_NH_A']:.3f}, r({labs['H']}...{labs['O']}) {rec['r_HO_A']:.3f}, "
                f"r({labs['N']}...{labs['O']}) {rec['r_NO_A']:.3f} A; E = {rec['E_Ha']:.8f} Ha"
                + (f"; vs stored ionic minimum {rec['dE_vs_stored_min_kcal']:+.4f} kcal/mol"
                   if 'dE_vs_stored_min_kcal' in rec else '')
                + (f"; heavy-atom RMSD to it {rec['heavy_atom_rmsd_to_stored_min_A']:.3f} A"
                   if 'heavy_atom_rmsd_to_stored_min_A' in rec else ''))
    job.r['summary'] = {m: {k: job.get(m).get(k) for k in
                            ('outcome', 'proton_on', 'returned_at_step', 'E_Ha', 'dE_vs_stored_min_kcal', 'r_NH_A',
                             'r_HO_A', 'r_NO_A', 'opt_steps', 'converged_on', 'heavy_atom_rmsd_to_stored_min_A')}
                        for m in models if job.get(m) and job.get(m).get('status') == 'done'}
    job.save()
    return True


if __name__ == '__main__':
    rc.install_sigterm_handler()
    rc.main_wrapper(parser, run, __file__)
