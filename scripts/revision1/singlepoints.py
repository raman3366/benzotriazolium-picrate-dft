"""
Single points with other functionals and DF-MP2 at given geometries (revision items 2 and 8).

Geometries come from proton_scan.py outputs (--scan DIR, repeatable), neutral_start.py outputs (--neutral DIR),
the stored minima (--minima: molecules/<system>/results/<system>_gas|pcm/opt.xyz) or plain XYZ files (--xyz).
Every method uses the basis of record, 6-311++G(d,p) (spherical), with the case's charge, spin and basis
overrides; the DFT settings are those of record (density fitting with PySCF's default auxiliary basis, grid level
4, conv_tol 1e-10, max_cycle 150; IEF-PCM with Lebedev order 29 when --solvent is given).

Methods (--methods, comma list):
  wb97x-d3bj   omega-B97X-D3(BJ): PySCF 'wb97x-d3bj' = libxc HYB_GGA_XC_WB97X_V (466) without VV10, plus D3(BJ)
               from pyscf-dispersion (Najibi & Goerigk's definition); checked on msi 20 Nov 2026
  m062x        M06-2X (libxc 450+236) at grid 4; a grid-5 check at two geometries (the lowest-energy and the most
               transferred one) gives the grid error of the absolute and of the relative energy
  pbe0, cam-b3lyp, b3lyp, b3lyp-d3bj   (b3lyp recomputes the reference; the scan's own B3LYP energies are
               tabulated anyway)
  dfmp2        RI-MP2 on a density-fitted RHF reference (RI-JK with PySCF's default JK-fit basis), frozen core
               (PySCF chemcore: 1s of C, N, O), RI auxiliary basis aug-cc-pVTZ-RI (PySCF's own default for
               6-311++G(d,p) with mp2fit=True; set explicitly, because DFRMP2 would otherwise reuse the JK-fit basis
               of a density-fitted reference).  Also SCS-MP2 (Grimme: 6/5 OS + 1/3 SS).  With --solvent: the
               reference is IEF-PCM HF and the reaction field is frozen at the HF density (PTE scheme).
Default methods: wb97x-d3bj, m062x, pbe0, cam-b3lyp, dfmp2.

Usage:
  python dft/revision/singlepoints.py --job ID --out DIR [--scan DIR ...] [--neutral DIR ...] [--minima]
         [--xyz FILE ...] [--case btpc_gas] [--select lowest|all] [--methods LIST] [--solvent gas|water|<eps>]
         [--grid-check-level 5] [--no-grid-check] [--mp2-aux aug-cc-pvtz-ri] [--max-geoms N] [--smoke]
--select lowest (default): where several scan chains reached the same target, only the lowest-B3LYP geometry.
Outputs in DIR: result.json (items '<geometry>|<method>'), log.txt, table.csv (long), table.md (relative
energies, kcal/mol, per geometry and method), geoms/<geometry>.xyz (the geometries used), pyscf/<method>.log,
chk/.  Restart: finished items are skipped; SCF guesses come from the item's own checkpoint, else the same
method at the previous geometry, else the scan's B3LYP checkpoint.
--smoke: STO-3G, grid 2 (check at 3), Lebedev 17, conv_tol 1e-6, at most 2 geometries; not for numbers.
"""
import sys, pathlib, argparse, json, time, gc, hashlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import rev_common as rc          # first: pins the thread pools
import numpy as np
import proton_scan as ps

DFT = {'wb97x-d3bj': 'wB97X-D3(BJ)', 'm062x': 'M06-2X', 'pbe0': 'PBE0', 'cam-b3lyp': 'CAM-B3LYP',
       'b3lyp': 'B3LYP', 'b3lyp-d3bj': 'B3LYP-D3(BJ)'}
ALL = list(DFT) + ['dfmp2']
DEFAULT_METHODS = 'wb97x-d3bj,m062x,pbe0,cam-b3lyp,dfmp2'
MP2_AUX = 'aug-cc-pvtz-ri'
SCS = (6 / 5, 1 / 3)       # Grimme 2003: OS, SS scaling


def parser():
    p = rc.add_common_args(argparse.ArgumentParser(description=__doc__.split('\n')[1]))
    p.add_argument('--scan', action='append', default=[], help='proton_scan.py output directory (repeatable)')
    p.add_argument('--neutral', action='append', default=[], help='neutral_start.py output directory (repeatable)')
    p.add_argument('--xyz', action='append', default=[], help='XYZ file (repeatable)')
    p.add_argument('--minima', action='store_true', help='add the stored gas and PCM minima of --system')
    p.add_argument('--system', default='btpc', help='system of --minima (default btpc)')
    p.add_argument('--case', default='btpc_gas', help='case giving charge, spin, basis overrides, fragments, proton')
    p.add_argument('--select', choices=('lowest', 'all'), default='lowest')
    p.add_argument('--methods', default=DEFAULT_METHODS, help=f'comma list from {ALL}')
    p.add_argument('--solvent', default=None, help='gas (default), water or a dielectric constant (IEF-PCM)')
    p.add_argument('--grid-check-level', type=int, default=5, help='grid level of the M06-2X check (default 5)')
    p.add_argument('--no-grid-check', action='store_true', help='skip the M06-2X grid check')
    p.add_argument('--mp2-aux', default=MP2_AUX, help=f'RI auxiliary basis of DF-MP2 (default {MP2_AUX})')
    p.add_argument('--max-geoms', type=int, default=None, help='use only the first N geometries (by r(H...O))')
    return p


def _safe(label):
    return label.replace(':', '__').replace('/', '_').replace('|', '_')


def _geom_hash(path):
    s, x = rc.read_xyz(path)
    return hashlib.sha256((' '.join(s) + np.array2string(np.round(x, 5), precision=5)).encode()).hexdigest()[:16]


def collect(args, case):
    """Geometries to compute: list of dicts (label, source, path, r_target_A, chain, E_b3lyp_Ha, chk)."""
    H, N, O = ps.proton_idx(case)
    geoms = []
    for d in args.scan:
        d = pathlib.Path(d).resolve()
        res = json.load(open(d / 'result.json'))
        if res.get('smoke') != bool(args.smoke):
            raise SystemExit(f'{d}: a {"smoke" if res.get("smoke") else "production"} scan; '
                             f'this is a {"smoke" if args.smoke else "production"} run')
        for key, it in res['items'].items():
            if it.get('status') != 'done' or it.get('E_Ha') is None:
                continue
            geoms.append({'label': key, 'source': f"scan {res['job']}", 'path': d / it['xyz'],
                          'r_target_A': it.get('r_target_A'), 'chain': it.get('chain'), 'E_b3lyp_Ha': it['E_Ha'],
                          'chk': (d / it['chk']) if it.get('chk') else None, 'scan_dir': str(d)})
    for d in args.neutral:
        d = pathlib.Path(d).resolve()
        res = json.load(open(d / 'result.json'))
        if bool(res.get('smoke')) != bool(args.smoke):
            raise SystemExit(f'{d}: smoke/production mismatch with this run')
        for key, it in res['items'].items():
            if it.get('status') == 'done' and it.get('xyz'):
                geoms.append({'label': f'neutral_{key}', 'source': f"neutral_start {res['job']}", 'path': d / it['xyz'],
                              'r_target_A': None, 'chain': f'neutral_start_{key}', 'E_b3lyp_Ha': None,
                              'chk': (d / it['chk']) if it.get('chk') else None})
    if args.minima:
        for model in ('gas', 'pcm'):
            c = rc.load_case(f'{args.system}_{model}')
            p = ps.stored_minimum(c)
            if p is not None:
                geoms.append({'label': f'min_{c["name"]}', 'source': 'stored minimum', 'path': p, 'r_target_A': None,
                              'chain': 'stored_min', 'E_b3lyp_Ha': None, 'chk': None})
    for x in args.xyz:
        p = pathlib.Path(x).resolve()
        geoms.append({'label': p.stem, 'source': 'xyz', 'path': p, 'r_target_A': None, 'chain': 'xyz',
                      'E_b3lyp_Ha': None, 'chk': None})
    # identical geometries (the free minimum appears in every chain that starts from it) are computed once
    seen, uniq = {}, []
    for g in geoms:
        h = _geom_hash(g['path'])
        if h in seen:
            continue
        seen[h] = g; uniq.append(g)
    geoms = uniq
    if args.select == 'lowest':
        best = {}
        for g in geoms:
            if g['r_target_A'] is None:
                continue
            t = round(g['r_target_A'], 3)
            if t not in best or (g['E_b3lyp_Ha'] or 0) < (best[t]['E_b3lyp_Ha'] or 0):
                best[t] = g
        keep = {id(g) for g in best.values()}
        geoms = [g for g in geoms if g['r_target_A'] is None or id(g) in keep]
    labels = [g['label'] for g in geoms]
    dup = {l for l in labels if labels.count(l) > 1}
    for g in geoms:
        if g['label'] in dup:
            g['label'] = f"{g['source'].split()[-1]}:{g['label']}"
        s, x = rc.read_xyz(g['path'])
        g.update({k: round(v, 5) for k, v in ps.measure(x, H, N, O).items()})
    geoms.sort(key=lambda g: -g['r_HO_A'])
    if args.max_geoms:
        geoms = geoms[:args.max_geoms]
    return geoms


def run_scf_robust(mf, guesses, log):
    """SCF from the first readable guess; on non-convergence one second-order (Newton) attempt from there."""
    g = next((c for c in guesses if c is not None and pathlib.Path(c).exists()), None)
    try:
        return rc.run_scf(mf, chk_guess=g, log=log), str(g) if g else 'minao', False
    except RuntimeError as e:
        log(f'{e}; second-order (Newton) SCF from the last orbitals')
        t = time.time()
        mo, occ = mf.mo_coeff, mf.mo_occ
        mf2 = mf.newton()
        mf2.max_cycle = 50
        mf2.kernel(mo, occ)
        if not mf2.converged:
            raise RuntimeError('SCF did not converge (DIIS and Newton)')
        mf.__dict__.update({k: getattr(mf2, k) for k in ('e_tot', 'mo_coeff', 'mo_occ', 'mo_energy', 'converged')})
        mf._cycles = getattr(mf, '_cycles', 0) + 50
        return time.time() - t, str(g) if g else 'minao', True


def dft_point(job, g, method, case, L, grid_level, prev_chk, key, logf):
    chk = job.chk_dir / f'{_safe(key)}.chk'
    mol = rc.build_mol(str(job.out / g['local']), case['charge'], case['spin'], basis=L['basis'],
                       basis_overrides=case.get('basis_overrides'), verbose=4, output=logf)
    mf = rc.make_mf(mol, xc=method, solvent=L['eps'], conv_tol=L['conv_tol'], grid_level=grid_level,
                    lebedev_order=L['pcm_lebedev_order'], chkfile=chk, max_cycle=rc.MAX_CYCLE)
    wall, guess, newton = run_scf_robust(mf, [chk, prev_chk, g.get('chk')], job.log)
    rec = {'geometry': g['label'], 'method': method, 'kind': 'dft', 'xc': mf.xc, 'grid_level': grid_level,
           'E_Ha': float(mf.e_tot), 'scf_cycles': mf._cycles, 'scf_wall_s': round(wall, 1), 'guess': guess,
           'newton_fallback': newton, 'scf_settings': rc.settings_of(mf), 'chk': str(chk.relative_to(job.out))}
    summ = getattr(mf, 'scf_summary', {}) or {}
    if 'dispersion' in summ:
        rec['E_disp_Ha'] = float(summ['dispersion'])
    if hasattr(mf, 'with_solvent') and getattr(mf.with_solvent, 'e', None) is not None:
        rec['E_solvation_Ha'] = float(mf.with_solvent.e)
    rec.update(ps.analyse(mf, case)); rec.update(rc.runtime_sizes(mf))
    mol.stdout.close()
    del mf; gc.collect()
    return rec, chk


def mp2_point(job, g, case, L, args, prev_chk, key, logf):
    from pyscf import scf, df
    from pyscf.mp import dfmp2
    from pyscf.data import elements
    chk = job.chk_dir / f'{_safe(key)}.chk'
    mol = rc.build_mol(str(job.out / g['local']), case['charge'], case['spin'], basis=L['basis'],
                       basis_overrides=case.get('basis_overrides'), verbose=4, output=logf)
    mf = scf.RHF(mol).density_fit()
    mf.conv_tol = L['conv_tol']; mf.max_cycle = rc.MAX_CYCLE; mf.chkfile = str(chk)
    if L['eps'] is not None:
        mf = mf.PCM()
        mf.with_solvent.method = rc.PCM_METHOD; mf.with_solvent.eps = L['eps']
        mf.with_solvent.lebedev_order = L['pcm_lebedev_order']
    wall, guess, newton = run_scf_robust(mf, [chk, prev_chk, g.get('chk')], job.log)
    e_hf = float(mf.e_tot)
    nfc = int(elements.chemcore(mol))
    t = time.time()
    pt = dfmp2.DFMP2(mf, frozen=nfc)
    pt.with_df = df.DF(mol, auxbasis=args.mp2_aux)        # RI basis for MP2, not the JK-fit basis of the reference
    pt.max_memory = mol.max_memory
    # Energy only: with_t2=False avoids holding the t2 amplitudes in memory (about 8.6 GB here; the default
    # with_t2=True raised MemoryError on the VM at an 8 GB cap, 22 Nov 2026). The correlation energy is unchanged.
    e_corr = pt.kernel(with_t2=False)[0]
    w_mp2 = time.time() - t
    os_ = float(getattr(pt, 'e_corr_os', getattr(e_corr, 'e_corr_os', np.nan)))
    ss = float(getattr(pt, 'e_corr_ss', getattr(e_corr, 'e_corr_ss', np.nan)))
    auxmol = getattr(pt.with_df, 'auxmol', None)
    occ = mf.mo_occ; n = int((occ > 0).sum())
    rec = {'geometry': g['label'], 'method': 'dfmp2', 'kind': 'mp2', 'E_HF_Ha': e_hf, 'E_corr_Ha': float(e_corr),
           'E_corr_os_Ha': os_, 'E_corr_ss_Ha': ss, 'E_Ha': e_hf + float(e_corr),
           'E_SCS_MP2_Ha': e_hf + SCS[0] * os_ + SCS[1] * ss, 'frozen_core_orbitals': nfc,
           'mp2_auxbasis': args.mp2_aux, 'n_aux_mp2': int(auxmol.nao_nr()) if auxmol is not None else None,
           'hf_auxbasis': 'PySCF default JK-fit', 'hf_scf_cycles': mf._cycles, 'hf_wall_s': round(wall, 1),
           'mp2_wall_s': round(w_mp2, 1), 'guess': guess, 'newton_fallback': newton,
           'HF_HOMO_eV': float(mf.mo_energy[n - 1] * rc.HARTREE_EV), 'HF_LUMO_eV': float(mf.mo_energy[n] * rc.HARTREE_EV),
           'scf_settings': rc.settings_of(mf), 'chk': str(chk.relative_to(job.out))}
    if L['eps'] is not None:
        rec['solvent_scheme'] = 'PTE: MP2 correlation with the IEF-PCM reaction field frozen at the HF density'
        rec['E_solvation_HF_Ha'] = float(mf.with_solvent.e)
    rec.update(rc.runtime_sizes(mf))
    mol.stdout.close()
    del mf, pt; gc.collect()
    return rec, chk


def tables(job, geoms, methods):
    items = job.r['items']
    cols = (['b3lyp_scan'] if any(g.get('E_b3lyp_Ha') is not None for g in geoms) else []) + methods
    if 'dfmp2' in methods:
        cols.append('scs-mp2')

    def energy(g, m):
        if m == 'b3lyp_scan':
            return g.get('E_b3lyp_Ha')
        it = items.get(f"{g['label']}|{'dfmp2' if m == 'scs-mp2' else m}")
        if not it or it.get('status') != 'done':
            return None
        return it.get('E_SCS_MP2_Ha') if m == 'scs-mp2' else it.get('E_Ha')
    ref = next((g for g in geoms if g['label'] == 'min' or g['label'].endswith(':min')), None)
    rows, refE = [], {}
    for m in cols:
        es = [(g, energy(g, m)) for g in geoms if energy(g, m) is not None]
        if not es:
            continue
        refE[m] = energy(ref, m) if ref is not None and energy(ref, m) is not None else min(e for _, e in es)
    for g in geoms:
        row = {'geometry': g['label'], 'r_HO_A': g['r_HO_A'], 'r_NH_A': g['r_NH_A'], 'r_NO_A': g['r_NO_A']}
        for m in cols:
            e = energy(g, m)
            row[m] = None if e is None or m not in refE else round((e - refE[m]) * rc.HARTREE_KCAL, 3)
        rows.append(row)
    with open(job.out / 'table.csv', 'w') as f:
        f.write('geometry,r_HO_A,r_NH_A,r_NO_A,method,E_Ha,dE_kcal\n')
        for g in geoms:
            for m in cols:
                e = energy(g, m)
                if e is not None and m in refE:
                    f.write(f"{g['label']},{g['r_HO_A']},{g['r_NH_A']},{g['r_NO_A']},{m},{e:.10f},"
                            f"{(e - refE[m]) * rc.HARTREE_KCAL:.4f}\n")
    with open(job.out / 'table.md', 'w') as f:
        f.write(f"Relative energies (kcal/mol) per method; zero = {'the free minimum (min)' if ref else 'the lowest geometry'}.\n\n")
        f.write('| geometry | r(H...O) A | r(N-H) A | ' + ' | '.join(cols) + ' |\n')
        f.write('|---|---|---|' + '---|' * len(cols) + '\n')
        for r in rows:
            f.write(f"| {r['geometry']} | {r['r_HO_A']:.3f} | {r['r_NH_A']:.3f} | "
                    + ' | '.join('' if r[m] is None else f'{r[m]:.2f}' for m in cols) + ' |\n')
    return {'reference_geometry': ref['label'] if ref else 'lowest per method', 'E_ref_Ha': refE,
            'columns': cols, 'rows': rows}


def run(job, args):
    case = rc.load_case(args.case)
    smoke = args.smoke
    eps = rc.solvent_eps(None if args.solvent in (None, 'gas', 'none') else args.solvent)
    L = {'basis': 'sto-3g' if smoke else rc.BASIS, 'grid_level': 2 if smoke else rc.GRID_LEVEL,
         'grid_check_level': (3 if smoke else args.grid_check_level),
         'pcm_lebedev_order': 17 if smoke else rc.PCM_LEBEDEV_ORDER, 'solvent': args.solvent if eps else None,
         'eps': eps, 'pcm_method': rc.PCM_METHOD if eps else None, 'conv_tol': 1e-6 if smoke else rc.CONV_TOL_SP,
         'max_cycle': rc.MAX_CYCLE, 'density_fit': True, 'dft_auxbasis': 'PySCF default (aug-cc-pVTZ-JKFIT)',
         'mp2_auxbasis': args.mp2_aux, 'mp2_frozen_core': 'pyscf.data.elements.chemcore', 'scs_mp2': SCS}
    methods = [m.strip().lower() for m in args.methods.split(',') if m.strip()]
    bad = [m for m in methods if m not in ALL]
    if bad:
        raise SystemExit(f'unknown method(s) {bad}; choose from {ALL}')
    if smoke and not args.max_geoms:
        args.max_geoms = 2
    job.r['script_sha256']['proton_scan.py'] = rc.sha256_file(pathlib.Path(ps.__file__).resolve())
    geoms = collect(args, case)
    if not geoms:
        raise SystemExit('no geometries (give --scan, --neutral, --minima or --xyz)')
    (job.out / 'geoms').mkdir(exist_ok=True); (job.out / 'pyscf').mkdir(exist_ok=True)
    for g in geoms:
        s, x = rc.read_xyz(g['path'])
        local = job.out / 'geoms' / f"{_safe(g['label'])}.xyz"
        if not local.exists():
            rc.write_xyz(local, s, x, f"{g['label']} from {g['source']} ({g['path']})")
        g['local'] = str(local.relative_to(job.out))
    job.r['settings'].update(L)
    job.r['settings'].update({'case': case['name'], 'methods': methods, 'select': args.select,
                              'geometries': [{k: (str(v) if isinstance(v, pathlib.Path) else v) for k, v in g.items()}
                                             for g in geoms]})
    job.input_file('case', rc.ROOT / case['case_file'])
    job.save()
    job.log(f"{len(geoms)} geometries x {len(methods)} methods ({', '.join(methods)}); "
            f"{'IEF-PCM eps ' + str(eps) if eps else 'gas phase'}; basis {L['basis']}, grid {L['grid_level']}")
    if smoke:
        job.log('SMOKE settings: STO-3G, grid 2, Lebedev 17, conv_tol 1e-6 -- not for production numbers')

    for method in methods:
        prev_chk = None
        logf = job.out / 'pyscf' / f'{method}.log'
        for g in geoms:
            key = f"{g['label']}|{method}"
            if job.done(key):
                prev_chk = job.out / job.get(key)['chk']
                continue
            job.put(key, {'geometry': g['label'], 'method': method}, status='running')
            t = time.time()
            if method == 'dfmp2':
                rec, prev_chk = mp2_point(job, g, case, L, args, prev_chk, key, logf)
            else:
                rec, prev_chk = dft_point(job, g, method, case, L, L['grid_level'], prev_chk, key, logf)
            rec.update({'r_HO_A': g['r_HO_A'], 'wall_s': round(time.time() - t, 1),
                        'peak_rss_mb': round(rc.peak_rss_mb(), 1)})
            job.put(key, rec)
            job.r['timings'][key] = rec['wall_s']; job.save()
            job.log(f"{key}: E = {rec['E_Ha']:.10f} Ha ({rec['wall_s']:.0f} s"
                    + (f", E_HF {rec['E_HF_Ha']:.8f}, E_corr {rec['E_corr_Ha']:.8f}" if method == 'dfmp2' else
                       f", {rec['scf_cycles']} cycles, LUMO on anion {rec['frontier_composition_mulliken_pct']['LUMO']['anion']:.1f} %")
                    + ')')
        # M06-2X grid check: at the lowest-energy geometry (by the method itself) and the most transferred one
        if method == 'm062x' and not args.no_grid_check:
            done = [g for g in geoms if job.done(f"{g['label']}|m062x")]
            if done:
                lo = min(done, key=lambda g: job.get(f"{g['label']}|m062x")['E_Ha'])
                tr = min(done, key=lambda g: g['r_HO_A'])
                checks = {}
                for g in ([lo] if lo is tr else [lo, tr]):
                    key = f"{g['label']}|m062x@grid{L['grid_check_level']}"
                    if not job.done(key):
                        rec, _ = dft_point(job, g, 'm062x', case, L, L['grid_check_level'],
                                           job.out / job.get(f"{g['label']}|m062x")['chk'], key,
                                           job.out / 'pyscf' / 'm062x_gridcheck.log')
                        rec['dE_vs_grid_default_Ha'] = rec['E_Ha'] - job.get(f"{g['label']}|m062x")['E_Ha']
                        job.put(key, rec)
                    checks[g['label']] = job.get(key)['dE_vs_grid_default_Ha']
                gc_rec = {'grid_default': L['grid_level'], 'grid_check': L['grid_check_level'],
                          'dE_abs_Ha': checks,
                          'dE_abs_kcal': {k: round(v * rc.HARTREE_KCAL, 4) for k, v in checks.items()}}
                if len(checks) == 2:
                    a, b = list(checks.values())
                    gc_rec['relative_energy_error_kcal'] = round((b - a) * rc.HARTREE_KCAL, 4)
                job.r['m062x_grid_check'] = gc_rec; job.save()
                job.log(f"M06-2X grid {L['grid_level']} -> {L['grid_check_level']}: {gc_rec}")
    job.r['summary'] = tables(job, geoms, methods)
    job.save()
    with open(job.out / 'table.md') as f:
        for line in f:
            job.log('  ' + line.rstrip())
    return True


if __name__ == '__main__':
    rc.install_sigterm_handler()
    rc.main_wrapper(parser, run, __file__)
