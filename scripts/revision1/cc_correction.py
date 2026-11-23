"""
Coupled-cluster correction to the DF-MP2 energies of the BTPC ion pair (revision; added 22 Nov 2026).

Why: in the IEF-PCM water model the hybrid functionals keep the salt (N4-H...O1) below the neutral pair (N4...H-O1)
while DF-MP2 and SCS-MP2 reverse the order (sp_min_pcm, sp_neutral_pcm).  CCSD(T) decides between them.

Composite (focal-point) scheme:  E[CCSD(T)] ~ E[DF-MP2 / 6-311++G(d,p)] + delta,
  delta = E_corr[CCSD(T)] - E_corr[MP2], both in a smaller basis (--basis) and in the same space of frozen natural
  orbitals (FNO): the virtual natural orbitals of the full-space MP2 density, kept until --fno-pct of the virtual
  occupation sum, semicanonicalised (PySCF mp.make_fno).  MP2 and CCSD(T) share the orbitals, the frozen core
  (PySCF chemcore: 1s of C, N, O), the density-fitted integrals (RI auxiliary basis from
  df.make_auxbasis(mol, mp2fit=True)) and the reference, so the FNO truncation and the DF error largely cancel in
  delta and cancel further in the difference between two geometries.
Reference: density-fitted RHF (PySCF default JK-fit basis), conv_tol 1e-10; with --solvent water it is IEF-PCM HF
(eps 78.3553, Lebedev 29) and the reaction field stays frozen at the HF density for MP2 and CCSD(T), as for the DF-MP2
numbers of singlepoints.py (PTE scheme): the correlation energies are computed with the solvent-polarised HF orbitals
and Fock matrix.

Per geometry and FNO size (item '<label>|f<frac>'): E_HF, E_corr MP2 in the full space (same RI integrals), the FNO
size and kept occupation, E_corr MP2 in the FNO space (from the CCSD initial amplitudes, RI integrals), E_corr CCSD, E(T),
delta_CCSD = E_CCSD - E_MP2(FNO), delta_T = E(T), delta = delta_CCSD + delta_T, the truncation estimates
dMP2 = MP2(full) - MP2(FNO) and (T) scaled by MP2(full)/MP2(FNO), the T1/D1/D2 diagnostics and |E(T)|/|E_CCSD|, and checks
(Fock off-diagonal norm in the semicanonical basis, CCSD convergence).  Summary: one entry per finished geometry and FNO size; the salt-minus-neutral
difference of delta is formed in the analysis (separate jobs per geometry run in parallel).

Usage:
  python dft/revision/cc_correction.py --job ID --out DIR --xyz FILE [--xyz FILE ...] --labels a[,b]
         [--solvent water] [--basis 6-31g*] [--fno-fracs 0.5,0.65,0.8] [--case btpc_gas]
Production truncation: --fno-fracs gives fixed numbers of virtual NOs (fraction x total virtuals, identical for both
geometries, since they have the same atoms and basis), run from the smallest; the series shows whether the
salt-minus-neutral difference of delta has converged.  (The selftest found that 0.99 of the occupation sum moves the
absolute delta by 1.6-3 kcal/mol in 6-31G* for small molecules, so a single aggressive threshold is not enough.)
  python dft/revision/cc_correction.py --job ID --out DIR --selftest
--selftest: water dimer and formamide-water, 6-31G*, gas and IEF-PCM: (1) FNO keeping all virtual NOs reproduces
plain CCSD(T) and MP2 (no truncation) to 1e-7 Ha; (2) the CCSD initial-amplitude MP2 equals the full-space DF-MP2
with the same RI basis; (3) the delta error of keeping 0.99 of the occupation sum (--fno-pct), with and without
the scaled (T).  Exit 1 on a failed check.
Restart: finished (geometry, FNO size) items are skipped; a restart redoes the SCF (from its checkpoint) and the
full-space MP2 for that geometry.
"""
import sys, pathlib, argparse, time, gc
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import rev_common as rc          # first: pins the thread pools
import numpy as np


def parser():
    p = rc.add_common_args(argparse.ArgumentParser(description=__doc__.split('\n')[1]))
    p.add_argument('--xyz', action='append', default=[], help='geometry (XYZ, Angstrom); repeatable')
    p.add_argument('--labels', default='', help='comma list, one label per --xyz')
    p.add_argument('--case', default='btpc_gas', help='case for charge and spin (default btpc_gas: 0, 0)')
    p.add_argument('--solvent', default=None, help='gas (default) | water | <eps>')
    p.add_argument('--basis', default='6-31g*')
    p.add_argument('--fno-pct', type=float, default=0.99, help='kept fraction of the virtual NO occupation sum (selftest)')
    p.add_argument('--fno-fracs', default='0.5,0.65,0.8',
                   help='production: comma list of kept fractions of the virtual space (fixed NO counts, the same for '
                        'every geometry), run in this order; one item per geometry and fraction')
    p.add_argument('--conv-tol', type=float, default=1e-8, help='CCSD energy convergence (Ha)')
    p.add_argument('--selftest', action='store_true')
    return p


def make_ref(mol, eps, conv_tol=1e-10, chk=None, guess=None, log=print):
    from pyscf import scf
    mf = scf.RHF(mol).density_fit()
    mf.conv_tol = conv_tol; mf.max_cycle = rc.MAX_CYCLE
    if chk:
        mf.chkfile = str(chk)
    if eps is not None:
        mf = mf.PCM()
        mf.with_solvent.method = rc.PCM_METHOD; mf.with_solvent.eps = eps
        mf.with_solvent.lebedev_order = rc.PCM_LEBEDEV_ORDER
    rc.run_scf(mf, chk_guess=guess, log=log)
    return mf


def ri_df(mol):
    """One RI auxiliary basis per geometry (df.make_auxbasis(mp2fit=True)), shared by the full-space MP2 and every
    CC step, so that MP2(full) - MP2(FNO) measures the truncation alone.  Returns the DF object and the basis spec."""
    from pyscf import df
    aux = df.make_auxbasis(mol, mp2fit=True)
    return df.DF(mol, auxbasis=aux), aux


def full_mp2(mf, nfc, rdf):
    """Full-space DF-MP2 with t2 kept (for the natural orbitals), on the shared RI basis."""
    from pyscf import mp
    pt = mp.MP2(mf, frozen=nfc)
    pt.with_df = rdf
    pt.max_memory = mf.mol.max_memory
    pt.kernel(with_t2=True)
    return pt


def no_sets(pt, nfc, nvir_list=(), pct=None):
    """Frozen lists and semicanonical NO coefficients for each requested number of kept virtual NOs (and for a
    kept fraction of the occupation sum if pct is given), built while t2 exists.  Returns (occupations, sets);
    sets[k] = (frozen, coeff, kept occupation fraction, smallest kept occupation)."""
    dm = pt.make_rdm1(with_frozen=False)
    n = np.sort(np.linalg.eigvalsh(dm[pt.nocc:, pt.nocc:]))[::-1]
    nvir = pt.nmo - pt.nocc
    sets = {}
    for k in list(nvir_list) + ([None] if pct is not None else []):
        frozen, coeff = pt.make_fno(pct_occ=pct, nvir_act=k) if k is None else pt.make_fno(nvir_act=k)
        kk = nvir - (len(frozen) - nfc)
        sets['pct' if k is None else kk] = (frozen, coeff, float(n[:kk].sum() / n.sum()), float(n[kk - 1]))
    return n, sets


def cc_in_space(mf, frozen, coeff, rdf, conv_tol=1e-8, log=print):
    """DF-CCSD(T) and MP2 (CCSD initial amplitudes) in the space given by frozen/coeff, on the shared RI basis."""
    from pyscf import cc
    mol = mf.mol
    out = {}
    mycc = cc.CCSD(mf, frozen=frozen, mo_coeff=coeff)
    out['cc_class'] = f'{type(mycc).__module__}.{type(mycc).__name__}'
    if 'dfccsd' not in out['cc_class']:
        raise RuntimeError(f'expected DF-CCSD, got {out["cc_class"]}')
    mycc.with_df = rdf
    mycc.max_memory = mol.max_memory
    mycc.conv_tol = conv_tol; mycc.conv_tol_normt = 1e-6; mycc.max_cycle = 100
    t = time.time()
    eris = mycc.ao2mo()
    nocc_a = eris.nocc
    f = eris.fock
    off = f - np.diag(np.diag(f))
    out['fock_offdiag_max'] = {'oo': float(np.abs(off[:nocc_a, :nocc_a]).max()),
                               'vv': float(np.abs(off[nocc_a:, nocc_a:]).max()),
                               'ov': float(np.abs(off[:nocc_a, nocc_a:]).max())}
    out['naux_cc'] = int(eris.naux) if hasattr(eris, 'naux') else None
    out['n_active_occ'] = int(nocc_a); out['n_active_vir'] = int(f.shape[0] - nocc_a)
    e_mp2 = float(mycc.init_amps(eris)[0])            # amplitudes not kept: kernel rebuilds the same guess
    out['E_corr_mp2_active_Ha'] = e_mp2
    log(f"  DF-CCSD space: {nocc_a} occ x {out['n_active_vir']} vir, naux {out['naux_cc']}; MP2 (same space, RI) "
        f"{e_mp2:.8f} Ha; |F_offdiag| max oo/vv/ov {out['fock_offdiag_max']['oo']:.1e}/"
        f"{out['fock_offdiag_max']['vv']:.1e}/{out['fock_offdiag_max']['ov']:.1e}")
    mycc.kernel(eris=eris)
    if not mycc.converged:
        raise RuntimeError('CCSD did not converge')
    out['E_corr_ccsd_Ha'] = float(mycc.e_corr); out['wall_ccsd_s'] = round(time.time() - t, 1)
    out['diagnostics'] = {'T1': float(mycc.get_t1_diagnostic()), 'D1': float(mycc.get_d1_diagnostic()),
                          'D2': float(mycc.get_d2_diagnostic())}
    log(f"  CCSD {mycc.e_corr:.8f} Ha ({out['wall_ccsd_s']:.0f} s); T1 {out['diagnostics']['T1']:.4f}, "
        f"D1 {out['diagnostics']['D1']:.4f}, D2 {out['diagnostics']['D2']:.4f}")
    t = time.time()
    et = float(mycc.ccsd_t(eris=eris))
    out['E_T_Ha'] = et; out['wall_t_s'] = round(time.time() - t, 1)
    out['diagnostics']['T_over_ccsd'] = abs(et) / abs(out['E_corr_ccsd_Ha'])
    out['delta_ccsd_Ha'] = out['E_corr_ccsd_Ha'] - e_mp2
    out['delta_T_Ha'] = et
    out['delta_Ha'] = out['delta_ccsd_Ha'] + et
    log(f"  (T) {et:.8f} Ha ({out['wall_t_s']:.0f} s); delta = CCSD(T) - MP2 = {out['delta_Ha']:.8f} Ha")
    del eris, mycc; gc.collect()
    return out


def add_truncation_estimates(rec, e_mp2_full):
    """Standard FNO truncation estimates from the shared-RI MP2: dMP2 = MP2(full) - MP2(FNO), and (T) scaled by
    MP2(full)/MP2(FNO).  delta itself needs no correction for the CCSD part to the extent CCSD tracks MP2."""
    rec['E_corr_mp2_full_ri_Ha'] = float(e_mp2_full)
    rec['dMP2_truncation_Ha'] = float(e_mp2_full - rec['E_corr_mp2_active_Ha'])
    rec['T_scaled_Ha'] = rec['E_T_Ha'] * e_mp2_full / rec['E_corr_mp2_active_Ha']
    rec['delta_T_scaled_Ha'] = rec['delta_ccsd_Ha'] + rec['T_scaled_Ha']
    return rec


def selftest(job, args):
    from pyscf.data import elements
    mols = {
        'water_dimer': 'O -1.551 -0.115 0.000; H -1.934 0.763 0.000; H -0.600 0.041 0.000; '
                       'O 1.351 0.111 0.000; H 1.680 -0.374 -0.759; H 1.680 -0.374 0.759',
        'formamide_water': 'C 0.000 0.418 0.000; O 1.204 0.418 0.000; N -0.680 1.590 0.000; H -0.512 -0.562 0.000; '
                           'H -0.190 2.470 0.000; H -1.690 1.560 0.000; O 2.050 2.900 0.000; H 1.680 2.000 0.000; '
                           'H 2.990 2.780 0.000'}
    ok = True
    res = {}
    for name, atom in mols.items():
        for eps in (None, 78.3553):
            tag = f"{name}_{'pcm' if eps else 'gas'}"
            mol = rc.build_mol([(a.split()[0], tuple(float(x) for x in a.split()[1:])) for a in atom.split(';')],
                               0, 0, basis='6-31g*', verbose=0)
            mf = make_ref(mol, eps, log=job.log)
            nfc = int(elements.chemcore(mol))
            rdf, aux = ri_df(mol)
            pt = full_mp2(mf, nfc, rdf)
            nvir = pt.nmo - pt.nocc
            n, sets = no_sets(pt, nfc, [nvir], pct=0.99)
            pt.t2 = None
            job.log(f'{tag}: E_HF {mf.e_tot:.10f}, MP2 full (RI) {pt.e_corr:.10f}')
            a = cc_in_space(mf, nfc, mf.mo_coeff, rdf, conv_tol=1e-10, log=job.log)
            b = cc_in_space(mf, sets[nvir][0], sets[nvir][1], rdf, conv_tol=1e-10, log=job.log)
            c = add_truncation_estimates(cc_in_space(mf, sets['pct'][0], sets['pct'][1], rdf, conv_tol=1e-10,
                                                     log=job.log), pt.e_corr)
            chk = {'plain_vs_fno100_ccsd': abs(a['E_corr_ccsd_Ha'] - b['E_corr_ccsd_Ha']),
                   'plain_vs_fno100_T': abs(a['E_T_Ha'] - b['E_T_Ha']),
                   'plain_vs_fno100_mp2': abs(a['E_corr_mp2_active_Ha'] - b['E_corr_mp2_active_Ha']),
                   'cc_init_mp2_vs_full_dfmp2_same_ri': abs(a['E_corr_mp2_active_Ha'] - pt.e_corr)}
            passed = all(v < 1e-7 for v in chk.values())
            ok &= passed
            res[tag] = {'checks_Ha': chk, 'passed': passed, 'delta_plain_Ha': a['delta_Ha'], 'delta_fno099_Ha': c['delta_Ha'],
                        'delta_fno099_Tscaled_Ha': c['delta_T_scaled_Ha'],
                        'delta_error_fno099_kcal': (c['delta_Ha'] - a['delta_Ha']) * rc.HARTREE_KCAL,
                        'delta_error_fno099_Tscaled_kcal': (c['delta_T_scaled_Ha'] - a['delta_Ha']) * rc.HARTREE_KCAL,
                        'nvir_kept_fno099': c['n_active_vir'], 'nvir_total': nvir, 'diagnostics_plain': a['diagnostics']}
            job.log(f"{tag}: {'PASS' if passed else 'FAIL'} {chk}; delta plain {a['delta_Ha']:.8f}, FNO 0.99 "
                    f"{c['delta_Ha']:.8f} ({res[tag]['delta_error_fno099_kcal']:+.4f} kcal/mol; with scaled (T) "
                    f"{res[tag]['delta_error_fno099_Tscaled_kcal']:+.4f}), {c['n_active_vir']}/{nvir} vir")
    job.r['selftest'] = res; job.save()
    return ok


def run(job, args):
    if args.selftest:
        return selftest(job, args)
    case = rc.load_case(args.case)
    eps = rc.solvent_eps(args.solvent)
    labels = [s.strip() for s in args.labels.split(',') if s.strip()]
    if not args.xyz or len(labels) != len(args.xyz):
        raise SystemExit('give one label per --xyz')
    job.r['settings'].update({'basis': args.basis, 'solvent': args.solvent if eps else None, 'eps': eps,
                              'pcm_method': rc.PCM_METHOD if eps else None,
                              'pcm_lebedev_order': rc.PCM_LEBEDEV_ORDER if eps else None,
                              'reference': 'DF-RHF (PySCF default JK-fit), conv_tol 1e-10' + (', IEF-PCM, reaction field frozen at HF (PTE)' if eps else ''),
                              'correlation': 'DF-CCSD(T) and MP2 in the FNO space; RI basis df.make_auxbasis(mp2fit=True); frozen core chemcore',
                              'cc_conv_tol': args.conv_tol, 'charge': case['charge'],
                              'spin': case['spin']})
    (job.out / 'pyscf').mkdir(exist_ok=True)
    for lab, x in zip(labels, args.xyz):
        job.input_file(f'xyz_{lab}', x)
    job.save()
    job.log(f"{len(labels)} geometries ({', '.join(labels)}); basis {args.basis}; "
            f"{'IEF-PCM eps ' + str(eps) if eps else 'gas phase'}; FNO fractions of the virtual space {args.fno_fracs}")
    from pyscf.data import elements
    fracs = [float(v) for v in args.fno_fracs.split(',') if v.strip()]
    job.r['settings']['fno_fracs'] = fracs
    job.save()
    for lab, x in zip(labels, args.xyz):
        keys = [f'{lab}|f{fr:.3f}' for fr in fracs]
        if all(job.done(k) for k in keys):
            continue
        t0 = time.time()
        mol = rc.build_mol(str(x), case['charge'], case['spin'], basis=args.basis, verbose=4,
                           output=job.out / 'pyscf' / f'{lab}.log')
        chk = job.chk_dir / f'{lab}.chk'
        mf = make_ref(mol, eps, chk=chk, guess=chk if chk.exists() else None, log=job.log)
        nfc = int(elements.chemcore(mol))
        rdf, aux = ri_df(mol)
        job.log(f"{lab}: E_HF = {mf.e_tot:.10f} Ha ({mf._cycles} cycles), nao {mol.nao}; RI basis {aux}; "
                f"full-space MP2 for the NOs")
        t = time.time()
        pt = full_mp2(mf, nfc, rdf)
        nvir_total = pt.nmo - pt.nocc
        ks = {fr: max(1, int(round(fr * nvir_total))) for fr in fracs}
        n, sets = no_sets(pt, nfc, sorted(set(ks.values())))
        e_mp2_full = float(pt.e_corr); nocc_act = int(pt.nocc)
        pt.t2 = None; del pt; gc.collect()                 # t2 released before any CCSD
        job.log(f"{lab}: MP2 full space (RI) {e_mp2_full:.8f} Ha ({time.time() - t:.0f} s), {nocc_act} occ x "
                f"{nvir_total} vir; NO sets {sorted(set(ks.values()))}")
        aux_rec = {'scf_jk': str(getattr(getattr(mf, 'with_df', None), 'auxbasis', None)), 'ri': aux}
        for fr, key in zip(fracs, keys):
            if job.done(key):
                continue
            k = ks[fr]
            frozen, coeff, occ_kept, occ_min = sets[k]
            job.put(key, {'geometry': lab, 'frac': fr, 'nvir_act': k}, status='running')
            t = time.time()
            job.log(f"{key}: FNO-CCSD(T) with {k} of {nvir_total} virtual NOs ({100 * occ_kept:.2f} % of the NO "
                    f"occupation, smallest kept {occ_min:.2e})")
            rec = {'geometry': lab, 'xyz': str(x), 'frac': fr, 'E_HF_Ha': float(mf.e_tot), 'n_frozen_core': nfc,
                   'nao': int(mol.nao), 'nvir_total': int(nvir_total), 'nvir_kept': int(k),
                   'no_occ_kept_frac': occ_kept, 'no_occ_smallest_kept': occ_min, 'no_occ_sum': float(n.sum()),
                   'auxbasis': aux_rec}
            rec.update(cc_in_space(mf, frozen, coeff, rdf, conv_tol=args.conv_tol, log=job.log))
            add_truncation_estimates(rec, e_mp2_full)
            rec.update({'wall_s': round(time.time() - t, 1), 'peak_rss_mb': round(rc.peak_rss_mb(), 1)})
            job.put(key, rec)
            job.r['timings'][key] = rec['wall_s']; job.save()
            job.log(f"{key}: delta_CCSD {rec['delta_ccsd_Ha']:.8f}, delta_T {rec['delta_T_Ha']:.8f}, delta "
                    f"{rec['delta_Ha']:.8f} Ha (with scaled (T) {rec['delta_T_scaled_Ha']:.8f}); MP2 truncation "
                    f"{rec['dMP2_truncation_Ha']:.6f} Ha ({rec['wall_s'] / 3600:.2f} h)")
        del mf, sets; gc.collect()
        job.log(f"{lab}: all fractions done ({(time.time() - t0) / 3600:.2f} h)")
    fields = ('geometry', 'frac', 'nvir_kept', 'nvir_total', 'no_occ_kept_frac', 'delta_ccsd_Ha', 'delta_T_Ha',
              'delta_Ha', 'delta_T_scaled_Ha', 'dMP2_truncation_Ha', 'E_HF_Ha', 'E_corr_mp2_full_ri_Ha',
              'E_corr_mp2_active_Ha', 'E_corr_ccsd_Ha', 'E_T_Ha', 'diagnostics', 'wall_s')
    job.r['summary'] = {k: {f: job.get(k).get(f) for f in fields} for k in job.r['items']
                        if '|f' in k and job.done(k)}
    job.save()
    return True


if __name__ == '__main__':
    rc.install_sigterm_handler()
    rc.main_wrapper(parser, run, __file__)
