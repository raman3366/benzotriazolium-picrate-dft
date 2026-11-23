"""
Delta-SCF check of where an added electron (N+1) and a removed electron (N-1) go in the BTPC ion pair (plan section 5,
item 3): an orbital-independent test of the LUMO-location claim.

For each (geometry, model) combination -- default all four of {gas minimum, PCM minimum} x {gas phase, IEF-PCM water},
i.e. the two published calculations and the two cross combinations -- the script computes
  N     neutral RKS (from the stored final.chk where one exists),
  N+1   UKS anion doublet at the same geometry (vertical), from one or more initial guesses that put the extra electron
        in the lowest cation-localised virtual, the lowest picrate-localised virtual, or the neutral LUMO,
  N-1   UKS cation doublet, hole in the neutral HOMO (or the highest cation/picrate-localised occupied orbital),
and reports Delta-SCF vertical electron affinity EA = E(N) - E(N+1) and ionisation energy IE = E(N-1) - E(N) (lowest
state of each charge), Koopmans values for comparison, <S^2>, and fragment spin and charge (Mulliken, Loewdin,
Hirshfeld) with the change of fragment charge from N.  With aufbau occupations an initial guess on the "wrong"
fragment may still relax to the other one; that is reported (start fragment vs final spin location).  --mom keeps the
initial occupation (maximum-overlap method) to obtain both localised states and their energy difference.
In IEF-PCM the continuum is fully relaxed for each charge state (equilibrium solvation at fixed geometry).

Level of theory: B3LYP/6-311++G(d,p), density fitting, grid 4, IEF-PCM Lebedev 29 (rev_common.py), conv_tol 1e-10
for N, 1e-9 for the doublets, max_cycle 150 with a level-shifted retry.
Usage:
  python dft/revision/deltascf.py --job dscf --out <dir> [--combos gas/gas,pcm/pcm,gas/pcm,pcm/gas]
         [--plus cation,anion] [--minus homo] [--mom] [--no-hirshfeld] [--case btpc] [--smoke]
  combo = <geometry>/<model>; geometry gas|pcm|label=path.xyz, model gas|pcm.
--smoke: STO-3G, grid 2, Lebedev 17, loose convergence, combos gas/gas,pcm/pcm.
"""
import sys, pathlib, argparse, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import rev_common as rc          # first: pins the thread pools
import numpy as np
import headline as hl
import charges as ch

CONV_DOUBLET = 1e-9


def parser():
    p = rc.add_common_args(argparse.ArgumentParser(description=__doc__.split('\n')[1]))
    p.add_argument('--combos', default=None, help='comma list of geometry/model (default all four)')
    p.add_argument('--plus', default='cation,anion', help='N+1 initial guesses: cation, anion, lumo')
    p.add_argument('--minus', default='homo', help='N-1 initial guesses: homo, cation, anion')
    p.add_argument('--mom', action='store_true', help='also run each guess with the maximum-overlap method')
    p.add_argument('--no-hirshfeld', action='store_true')
    p.add_argument('--case', default='btpc')
    p.add_argument('--xc', default=rc.XC)
    return p


def pick_orbital(fr, kind, plus):
    """Index of the neutral MO that receives the electron (N+1) or loses it (N-1)."""
    if plus:
        if kind == 'lumo':
            return fr['n_occ']
        o = fr['lowest_virtual_on'].get(kind)
    else:
        if kind == 'homo':
            return fr['n_occ'] - 1
        o = fr['highest_occupied_on'].get(kind)
    return None if o is None else o['index']


def uks_guess(C, nocc, k, plus):
    """UKS density from RKS orbitals: alpha gets orbital k as well (N+1), or beta loses orbital k (N-1)."""
    nmo = C.shape[1]
    occ_a = np.zeros(nmo); occ_b = np.zeros(nmo)
    occ_a[:nocc] = 1; occ_b[:nocc] = 1
    if plus:
        occ_a[k] = 1
    else:
        occ_b[k] = 0
    dm = np.array([(C * occ_a) @ C.T, (C * occ_b) @ C.T])
    return dm, (occ_a, occ_b)


def run(job, args):
    lv = hl.level(args.smoke)
    combos = (args.combos or ('gas/gas,pcm/pcm' if args.smoke else 'gas/gas,pcm/pcm,gas/pcm,pcm/gas')).split(',')
    plus_guesses = [g.strip() for g in args.plus.split(',') if g.strip()]
    minus_guesses = [g.strip() for g in args.minus.split(',') if g.strip()]
    job.r['settings'] = {**lv, 'xc': args.xc, 'conv_tol_doublet': 1e-6 if args.smoke else CONV_DOUBLET,
                         'combos': combos, 'plus_guesses': plus_guesses, 'minus_guesses': minus_guesses,
                         'mom': args.mom, 'pcm': f'{rc.PCM_METHOD} water eps {rc.SOLVENT_EPS["water"]}, equilibrium '
                         'solvation for each charge state', 'density_fit': True, 'auxbasis': 'PySCF default',
                         'spin_convention': 'mol.spin = 2S = 1 for both doublets'}
    job.save()
    case = rc.load_case(f'{args.case}_gas')
    frags = rc.fragments(case)
    conv_d = 1e-6 if args.smoke else CONV_DOUBLET
    free_cache = {}

    def analyse(mol, mf, dm):
        out = {k: {'fragments': v['fragments'], **({'spin_fragments': v['spin_fragments']} if 'spin_fragments' in v
                                                   else {})}
               for k, v in hl.populations(mol, dm, frags).items()}
        if not args.no_hirshfeld:
            key = (mol.basis if isinstance(mol.basis, str) else str(mol.basis), args.xc)
            if key not in free_cache:
                free_cache[key] = ch.free_atom_dms(mol, xc=args.xc)
            h = ch.hirshfeld_charges(mol, dm, mf.grids, free_cache[key])
            q = np.array(h['charges'])
            out['hirshfeld'] = {'fragments': {f: float(q[a].sum()) for f, a in frags.items()}}
            if 'spins' in h:
                s = np.array(h['spins'])
                out['hirshfeld']['spin_fragments'] = {f: float(s[a].sum()) for f, a in frags.items()}
        return out

    summary = {}
    for combo in combos:
        gspec, model = combo.strip().split('/')
        glabel, gpath, gmodel = hl.resolve_geom(gspec, args.case)
        solvent = None if model == 'gas' else 'water'
        tag = f'{glabel}/{model}'
        job.input_file(f'geometry_{glabel}', gpath)

        def build(charge, spin):
            mol = rc.mol_for_case(case, gpath, charge=charge, spin=spin, basis=lv['basis'], verbose=4 if not args.smoke else 3,
                                  output=job.out / f'pyscf_{glabel}_{model}.log')
            return mol

        # ---- N (neutral, RKS)
        kN = f'{tag}|N'
        chkN = job.chk_dir / f'{glabel}_{model}_N.chk'
        mol = build(0, 0)
        mfN = rc.make_mf(mol, xc=args.xc, solvent=solvent, conv_tol=lv['conv_tol'], grid_level=lv['grid_level'],
                         lebedev_order=lv['lebedev_order'], chkfile=chkN)
        guess = chkN if chkN.exists() else None
        if guess is None and not args.smoke:
            d = hl.stored_dir(f'{args.case}_{model}') or None
            if d and (d / 'final.chk').exists():
                guess = d / 'final.chk'
        if job.done(kN) and chkN.exists():
            # rebuild the neutral orbitals cheaply from the checkpoint (needed for the guesses)
            wall = rc.run_scf(mfN, chk_guess=chkN, log=job.log)
            recN = job.get(kN)
        else:
            recN = {'geometry': glabel, 'xyz': str(gpath), 'model': model, 'charge': 0, 'spin': 0,
                    'settings': rc.settings_of(mfN)}
            job.put(kN, recN, status='running')
            job.log(f'{kN}: RKS ({mol.nao} AO), guess {guess or "minao"}')
            wall = rc.run_scf(mfN, chk_guess=guess, log=job.log)
            fr = hl.frontier(mfN, frags)
            recN.update({'E_Ha': float(mfN.e_tot), 'scf_cycles': mfN._cycles, 'scf_wall_s': round(wall, 1),
                         'frontier': fr, 'populations': analyse(mol, mfN, mfN.make_rdm1())})
            job.put(kN, recN)
            job.r['timings'][kN] = round(wall, 1)
            job.log(f"{kN}: E = {mfN.e_tot:.10f} Ha ({mfN._cycles} cycles, {wall:.0f} s); LUMO on {fr['LUMO_on']}, "
                    f"HOMO on {fr['HOMO_on']}")
        frN = recN['frontier']
        C = np.asarray(mfN.mo_coeff); nocc = frN['n_occ']

        # ---- N+1 and N-1 (UKS doublets)
        results = {'N+1': [], 'N-1': []}
        for plus, glist in ((True, plus_guesses), (False, minus_guesses)):
            state = 'N+1' if plus else 'N-1'
            used = set()
            for g in glist:
                k = pick_orbital(frN, g, plus)
                if k is None:
                    job.log(f'{tag} {state}: no orbital for guess {g!r}; skipped')
                    continue
                for mode in (['aufbau', 'mom'] if args.mom else ['aufbau']):
                    if (k, mode) in used:
                        job.log(f'{tag} {state}: guess {g!r} is orbital {k}, already done ({mode}); skipped')
                        continue
                    used.add((k, mode))
                    key = f'{tag}|{state}|{g}|{mode}'
                    if job.done(key):
                        results[state].append(job.get(key)); continue
                    charge = -1 if plus else +1
                    molD = build(charge, 1)
                    chk = job.chk_dir / f"{glabel}_{model}_{state.replace('+', 'p').replace('-', 'm')}_{g}_{mode}.chk"
                    mf = rc.make_mf(molD, xc=args.xc, solvent=solvent, conv_tol=conv_d, grid_level=lv['grid_level'],
                                    lebedev_order=lv['lebedev_order'], restricted=False, chkfile=chk)
                    dm0, (oa, ob) = uks_guess(C, nocc, k, plus)
                    if mode == 'mom':
                        from pyscf import scf
                        mf = scf.addons.mom_occ(mf, (C, C), (oa, ob))
                    start_frag = hl.frag_of_share({f: 100 * rc.mo_fragment_share(mol, C[:, k], a) for f, a in frags.items()})
                    rec = {'geometry': glabel, 'model': model, 'state': state, 'charge': charge, 'spin': 1,
                           'guess': g, 'guess_orbital_index': int(k), 'guess_orbital_on': start_frag, 'mode': mode,
                           'settings': rc.settings_of(mf)}
                    job.put(key, rec, status='running')
                    job.log(f'{key}: UKS doublet, electron {"added to" if plus else "removed from"} MO {k} '
                            f'({start_frag}), {mode}')
                    t = time.time()
                    attempts = []
                    try:
                        if chk.exists():
                            rc.run_scf(mf, chk_guess=chk, log=job.log)
                        else:
                            rc.run_scf(mf, dm0=dm0, log=job.log)
                        attempts.append({'level_shift': 0, 'cycles': mf._cycles})
                    except RuntimeError as ex:
                        attempts.append({'level_shift': 0, 'cycles': mf.max_cycle, 'error': str(ex)})
                        job.log(f'{key}: not converged; retry with level shift 0.25 Ha from the last density')
                        dm_last = mf.make_rdm1()
                        mf.level_shift = 0.25
                        rc.run_scf(mf, dm0=dm_last, log=job.log)
                        attempts.append({'level_shift': 0.25, 'cycles': mf._cycles})
                        mf.level_shift = 0
                        rc.run_scf(mf, dm0=mf.make_rdm1(), log=job.log)      # polish without the shift
                        attempts.append({'level_shift': 0, 'cycles': mf._cycles, 'polish': True})
                    wall = time.time() - t
                    s2, mult = mf.spin_square()
                    dm = mf.make_rdm1()
                    pops = analyse(molD, mf, dm)
                    # where did the unpaired electron end up?  (Mulliken spin; Hirshfeld if computed)
                    sp = (pops.get('hirshfeld') or pops['mulliken'])['spin_fragments']
                    end_frag = max(sp, key=lambda f: abs(sp[f]))
                    na = int(round(mf.mo_occ[0].sum())); nb = int(round(mf.mo_occ[1].sum()))
                    # the orbital of the added electron (alpha HOMO) or of the hole (beta LUMO)
                    s_idx, s_spin, s_lab = (na - 1, 0, 'alpha HOMO') if plus else (nb, 1, 'beta LUMO')
                    soc = mf.mo_coeff[s_spin][:, s_idx]
                    somo_sh = {f: round(100 * rc.mo_fragment_share(molD, soc, a), 2) for f, a in frags.items()}
                    dq = {m: {f: pops[m]['fragments'][f] - recN['populations'][m]['fragments'][f] for f in frags}
                          for m in pops if m in recN['populations']}
                    rec.update({'E_Ha': float(mf.e_tot), 'S2': float(s2), 'multiplicity': float(mult),
                                'scf_attempts': attempts, 'scf_wall_s': round(wall, 1), 'populations': pops,
                                'delta_charge_vs_N': dq, 'spin_on': end_frag,
                                'stayed_on_guess_fragment': end_frag == start_frag,
                                'unpaired_orbital': s_lab,
                                'unpaired_orbital_eV': float(mf.mo_energy[s_spin][s_idx] * rc.HARTREE_EV),
                                'unpaired_orbital_shares_mulliken_pct': somo_sh})
                    job.put(key, rec)
                    job.r['timings'][key] = round(wall, 1)
                    job.log(f"{key}: E = {mf.e_tot:.10f} Ha, <S2> {s2:.4f}, spin on {end_frag} "
                            f"({', '.join(f'{f} {v:+.3f}' for f, v in sp.items())}), {wall:.0f} s")
                    molD.stdout.close()
                    results[state].append(job.get(key))
        # ---- Delta-SCF numbers for this combination
        EN = recN['E_Ha']
        s = {'E_N_Ha': EN, 'koopmans_IP_eV': -frN['HOMO_eV'], 'koopmans_EA_eV': -frN['LUMO_eV'],
             'LUMO_on': frN['LUMO_on'], 'HOMO_on': frN['HOMO_on']}
        for state, lst in results.items():
            lst = [r for r in lst if r and r.get('status') == 'done']
            if not lst:
                continue
            best = min(lst, key=lambda r: r['E_Ha'])
            s[f'{state}_states'] = [{'guess': r['guess'], 'mode': r['mode'], 'E_Ha': r['E_Ha'], 'spin_on': r['spin_on'],
                                     'S2': r['S2'], 'E_rel_eV': (r['E_Ha'] - best['E_Ha']) * rc.HARTREE_EV}
                                    for r in lst]
            s[f'{state}_lowest'] = {'guess': best['guess'], 'mode': best['mode'], 'spin_on': best['spin_on']}
            if state == 'N+1':
                s['deltaSCF_EA_eV'] = (EN - best['E_Ha']) * rc.HARTREE_EV
                s['added_electron_on'] = best['spin_on']
            else:
                s['deltaSCF_IE_eV'] = (best['E_Ha'] - EN) * rc.HARTREE_EV
                s['hole_on'] = best['spin_on']
        summary[tag] = s
        job.log(f"{tag}: EA(dSCF) {s.get('deltaSCF_EA_eV', float('nan')):.3f} eV (electron on {s.get('added_electron_on')}),"
                f" IE(dSCF) {s.get('deltaSCF_IE_eV', float('nan')):.3f} eV (hole on {s.get('hole_on')}); Koopmans "
                f"EA {s['koopmans_EA_eV']:.3f}, IP {s['koopmans_IP_eV']:.3f} eV; LUMO on {s['LUMO_on']}")
        mol.stdout.close()
        job.r['summary'] = summary
        job.save()
    return True


if __name__ == '__main__':
    rc.install_sigterm_handler()
    rc.main_wrapper(parser, run, __file__)
