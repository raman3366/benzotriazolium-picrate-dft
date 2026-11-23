"""
Hirshfeld and CM5 charges for the BTPC ion pair (plan section 5, item 3: "Hirshfeld or CM5 charges in place of
Mulliken"), compared with the Mulliken and Loewdin values of the paper.

Hirshfeld (Theor. Chim. Acta 44, 129 (1977)) is implemented here in PySCF, so it runs on every host (Multiwfn has no
ARM64 binary):  q_A = Z_A - integral w_A(r) rho(r) dr,  w_A = rho_A^0 / sum_B rho_B^0, with rho_A^0 the spherically
averaged, spin-restricted neutral free-atom density at the same functional and basis (pyscf.scf.atom_ks.get_atm_nrks,
fractional occupations), and the integral taken on the molecule's DFT grid (level 4 in production).  Spin populations
of a UKS density use the same weights.
CM5 (Marenich, Jerome, Cramer, Truhlar, J. Chem. Theory Comput. 8, 527 (2012)) adds the pairwise geometric
correction  q_k^CM5 = q_k^H + sum_k' T_kk' exp[-alpha (r_kk' - R_k - R_k')]  to the Hirshfeld charges (parameters
below, H/C/N/O only, which is all BTPC needs).  When Multiwfn is available (msi), its Hirshfeld (built-in free-atom
densities) and CM5 are computed from a molden file of the same wavefunction and compared atom by atom; Multiwfn's
CM5 - Hirshfeld difference checks the CM5 parameters exactly (it depends on the geometry only).

Sources (--source, comma list):
  stored            data/r2/molecules/<case>/results/<case>_{gas,pcm}/final.chk (the published wavefunctions; Mulliken
                    and Loewdin from them are checked against papers/btpc/data/facts.json)
  label=chk:PATH    any PySCF checkpoint holding mol and scf (e.g. an opt_d3.py output: <out>/chk/final.chk)
  label=xyz:PATH[:gas|pcm]  a new B3LYP single point at that geometry (gas default), then the analysis
Usage:
  python dft/revision/charges.py --job charges --out <dir> [--source stored] [--case btpc] [--multiwfn auto|none|PATH]
--smoke: a new STO-3G single point at the stored gas geometry (grid 2), PySCF Hirshfeld/CM5, and Multiwfn on its molden.

The functions free_atom_dms(), hirshfeld(), cm5() are imported by deltascf.py and opt_d3.py.
"""
import sys, os, pathlib, argparse, json, re, subprocess, shutil, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import rev_common as rc          # first: pins the thread pools
import numpy as np
import headline as hl

FACTS = rc.ROOT / 'papers' / 'btpc' / 'data' / 'facts.json'

# ---------------------------------------------------------------- CM5 parameters (Marenich et al. 2012, Table 2)
CM5_ALPHA = 2.474                                     # 1/Angstrom
# Atomic radii (Angstrom): the mean of the Pyykko-Atsumi (2009) and Cordero et al. (2008) single-bond covalent radii,
# which is what Multiwfn 2026.11.13 uses for CM5: its CM5 - Hirshfeld corrections for the 33-atom ion pair are reproduced
# to 8e-7 e with these radii and the pair parameters below (both radius tables are present in the Multiwfn binary),
# whereas the rounded CRC values (H 0.32, C 0.75, N 0.71, O 0.64) differ from it by up to 2.5e-3 e (20 Nov 2026).
CM5_RADII = {'H': 0.315, 'C': 0.755, 'N': 0.71, 'O': 0.645}
CM5_D = {'H': 0.0056, 'C': -0.0446, 'N': -0.1072, 'O': -0.0802}           # element parameters D_Z
CM5_PAIR = {('H', 'C'): 0.0502, ('H', 'N'): 0.1747, ('H', 'O'): 0.1671,   # pair-specific D_ZZ' (T_ZZ' = -T_Z'Z)
            ('C', 'N'): 0.0556, ('C', 'O'): 0.0234, ('N', 'O'): -0.0346}


def cm5_T(a, b):
    if a == b:
        return 0.0
    if (a, b) in CM5_PAIR:
        return CM5_PAIR[(a, b)]
    if (b, a) in CM5_PAIR:
        return -CM5_PAIR[(b, a)]
    return CM5_D[a] - CM5_D[b]


def cm5(syms, xyz_A, q_hirshfeld):
    """CM5 charges from Hirshfeld charges and the geometry (Angstrom)."""
    missing = sorted(set(syms) - set(CM5_RADII))
    if missing:
        raise ValueError(f'CM5 parameters not tabulated here for {missing}')
    X = np.asarray(xyz_A, float)
    q = np.array(q_hirshfeld, float).copy()
    corr = np.zeros(len(syms))
    for k in range(len(syms)):
        for j in range(len(syms)):
            if j == k:
                continue
            r = np.linalg.norm(X[k] - X[j])
            B = np.exp(-CM5_ALPHA * (r - CM5_RADII[syms[k]] - CM5_RADII[syms[j]]))
            corr[k] += cm5_T(syms[k], syms[j]) * B
    return q + corr, corr


# ---------------------------------------------------------------- Hirshfeld
def free_atom_dms(mol, xc='b3lyp'):
    """Spherically averaged neutral free-atom density matrices {element: dm in that atom's AO block}."""
    from pyscf.scf import atom_ks
    v = mol.verbose
    mol.verbose = 0
    try:
        res = atom_ks.get_atm_nrks(mol, xc=xc)
    finally:
        mol.verbose = v
    out = {}
    for el, (e, mo_e, C, occ) in res.items():
        out[el] = (C * occ) @ C.T
    return out


def hirshfeld(mol, dms, grids, free, max_memory=None):
    """Hirshfeld populations of each density matrix in `dms` (list) on `grids`.
    Returns {'pop': (ndm, natm) array of integral w_A rho_k, 'integral': per-dm total, 'free_integral': per-atom
    integral of rho_A^0 (should equal Z_A), 'promolecule_integral'}."""
    from pyscf import dft
    ni = dft.numint.NumInt()
    nao = mol.nao
    sl = mol.aoslice_by_atom()
    natm = mol.natm
    el = [mol._atom[a][0] for a in range(natm)]
    pops = np.zeros((len(dms), natm)); tot = np.zeros(len(dms)); free_int = np.zeros(natm)
    if grids.coords is None:
        grids.build(with_non0tab=True)
    mm = max_memory or max(500, rc.MAX_MEMORY_MB // 2)
    for ao, mask, weight, coords in ni.block_loop(mol, grids, nao, deriv=0, max_memory=mm):
        rho0 = np.empty((natm, weight.size))
        for a in range(natm):
            p0, p1 = sl[a][2], sl[a][3]
            aoA = np.ascontiguousarray(ao[:, p0:p1])
            rho0[a] = np.einsum('pi,pi->p', aoA @ free[el[a]], aoA)
        pro = rho0.sum(0)
        w = np.divide(rho0, pro, out=np.zeros_like(rho0), where=pro > 1e-30)
        for k, dm in enumerate(dms):
            rho = ni.eval_rho(mol, ao, np.asarray(dm), mask, xctype='LDA', hermi=1)
            rw = rho * weight
            pops[k] += w @ rw
            tot[k] += rw.sum()
        free_int += rho0 @ weight
    return {'pop': pops, 'integral': tot, 'free_integral': free_int}


def hirshfeld_charges(mol, dm, grids, free):
    """Hirshfeld charges (and spins for a UKS (2, nao, nao) density) per atom."""
    dm = np.asarray(dm)
    uks = dm.ndim == 3
    dms = [dm[0] + dm[1], dm[0] - dm[1]] if uks else [dm]
    h = hirshfeld(mol, dms, grids, free)
    Z = mol.atom_charges()
    out = {'charges': (Z - h['pop'][0]).tolist(), 'n_electrons_integrated': float(h['integral'][0]),
           'free_atom_integral_minus_Z_max': float(np.abs(h['free_integral'] - Z).max())}
    if uks:
        out['spins'] = h['pop'][1].tolist()
        out['spin_integrated'] = float(h['integral'][1])
    return out


def full_analysis(mol, dm, grids, free, frags, case):
    """Mulliken, Loewdin, Hirshfeld, CM5 charges per atom and per fragment (spins for UKS)."""
    syms, X = hl.mol_geometry(mol)
    pop = hl.populations(mol, dm, frags)
    h = hirshfeld_charges(mol, dm, grids, free)
    q_cm5, corr = cm5(syms, X, h['charges'])
    res = {'mulliken': pop['mulliken'], 'lowdin': pop['lowdin'],
           'hirshfeld': {'atoms': h['charges'],
                         'fragments': {f: float(np.sum(np.array(h['charges'])[a])) for f, a in frags.items()},
                         'n_electrons_integrated': h['n_electrons_integrated'],
                         'free_atom_integral_minus_Z_max': h['free_atom_integral_minus_Z_max']},
           'cm5': {'atoms': q_cm5.tolist(), 'correction': corr.tolist(),
                   'fragments': {f: float(q_cm5[a].sum()) for f, a in frags.items()}}}
    if 'spins' in h:
        res['hirshfeld']['spin_atoms'] = h['spins']
        res['hirshfeld']['spin_fragments'] = {f: float(np.sum(np.array(h['spins'])[a])) for f, a in frags.items()}
    p = case.get('proton')
    if p:
        res['proton_H'] = {m: float(res[m]['atoms'][p['H']]) for m in ('mulliken', 'lowdin', 'hirshfeld', 'cm5')}
    return res


# ---------------------------------------------------------------- Multiwfn
def find_multiwfn(spec):
    if spec in (None, '', 'none'):
        return None
    if spec != 'auto':
        return spec if os.path.exists(spec) else None
    for c in (os.environ.get('MWFN_BIN'), str(rc.ROOT / 'tools' / 'multiwfn' / 'Multiwfn_noGUI'),
              os.path.expanduser('~/r007/multiwfn/Multiwfn_noGUI'), shutil.which('Multiwfn_noGUI'),
              shutil.which('Multiwfn')):
        if c and os.path.exists(c):
            return c
    return None


def multiwfn_hirshfeld_cm5(binary, molden, workdir, threads):
    """Run Multiwfn (menu 7 -> 16: Hirshfeld with built-in free-atom densities, then CM5) on a molden file."""
    workdir = pathlib.Path(workdir); workdir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, Multiwfnpath=str(pathlib.Path(binary).resolve().parent))
    t = time.time()
    p = subprocess.run([binary, str(molden), '-nt', str(threads)], input='7\n16\n1\nn\n0\nq\n', text=True,
                       capture_output=True, cwd=str(workdir), env=env, timeout=3600)
    out = p.stdout
    (workdir / 'multiwfn_cm5.out').write_text(out)
    rows = re.findall(r'Atom:\s+(\d+)\s*([A-Za-z]+)\s+CM5 charge:\s+(-?\d+\.\d+)\s+Hirshfeld charge:\s+(-?\d+\.\d+)', out)
    if not rows:
        raise RuntimeError(f'Multiwfn output not parsed (rc {p.returncode}); see {workdir}/multiwfn_cm5.out')
    m = re.search(r'Multiwfn.*?[Vv]ersion\s+([\w.\-()]+)', out)
    return {'binary': binary, 'version_line': m.group(0)[:120] if m else None, 'wall_s': round(time.time() - t, 1),
            'cm5': [float(r[2]) for r in rows], 'hirshfeld': [float(r[3]) for r in rows],
            'note': 'raw values before Multiwfn renormalises the sum to the electron count; built-in sphericalized '
                    'free-atom densities'}


def write_molden(mol, mo_coeff, mo_energy, mo_occ, path):
    """Molden file padded to a square coefficient matrix (as dft/fmo.py does, for Multiwfn)."""
    from pyscf.tools import molden
    C = np.asarray(mo_coeff); nao, nmo = C.shape
    e = np.asarray(mo_energy); occ = np.asarray(mo_occ)
    if nmo < nao:
        C2 = np.zeros((nao, nao)); C2[:, :nmo] = C
        e2 = np.zeros(nao); e2[:nmo] = e; e2[nmo:] = e.max() + 1e3
        o2 = np.zeros(nao); o2[:nmo] = occ
        C, e, occ = C2, e2, o2
    molden.from_mo(mol, str(path), C, ene=e, occ=occ)


# ---------------------------------------------------------------- driver
def parser():
    p = rc.add_common_args(argparse.ArgumentParser(description=__doc__.split('\n')[1]))
    p.add_argument('--source', default='stored', help='comma list: stored | label=chk:PATH | label=xyz:PATH[:gas|pcm]')
    p.add_argument('--case', default='btpc')
    p.add_argument('--multiwfn', default='auto', help="auto (msi: ~/r007/multiwfn), none, or the binary's path")
    p.add_argument('--xc', default=rc.XC, help='functional for the free atoms and any new single point')
    return p


def facts_values():
    try:
        d = json.load(open(FACTS))
    except Exception:
        return {}
    out = {}
    for model in ('gas', 'pcm'):
        for k in ('mulliken_charge_picrate', 'lowdin_charge_picrate', 'mulliken_charge_H4N', 'lowdin_charge_H4N'):
            try:                      # facts.json is flat: {"gas.mulliken_charge_picrate": {"value": ...}, ...}
                out[f'{model}.{k}'] = float(d[f'{model}.{k}']['value'])
            except Exception:
                pass
    return out


def expand_sources(spec, case_prefix, smoke):
    srcs = []
    for s in [x.strip() for x in spec.split(',') if x.strip()]:
        if s == 'stored':
            for model in ('gas', 'pcm'):
                d = hl.stored_dir(f'{case_prefix}_{model}')
                srcs.append((model, 'chk', d / 'final.chk' if d else None, model))
        else:
            label, _, rest = s.partition('=')
            kind, _, path = rest.partition(':')
            model = 'gas'
            if kind == 'xyz' and path.rsplit(':', 1)[-1] in ('gas', 'pcm'):
                path, model = path.rsplit(':', 1)
            p = pathlib.Path(path)
            srcs.append((label, kind, p if p.is_absolute() else rc.ROOT / p, model))
    if smoke:
        srcs = [('smoke_gas', 'xyz', rc.ROOT / rc.STORED[f'{case_prefix}_gas']['xyz'] if f'{case_prefix}_gas' in rc.STORED
                 else hl.resolve_geom('gas', case_prefix)[1], 'gas')]
    return srcs


def run(job, args):
    lv = hl.level(args.smoke)
    mwfn = find_multiwfn(args.multiwfn)
    job.r['settings'] = {'hirshfeld': 'PySCF implementation (this script): neutral spherically averaged spin-restricted '
                                      'free atoms, get_atm_nrks(xc), on the molecular DFT grid',
                         'free_atom_xc': args.xc, 'grid_level': lv['grid_level'], 'basis': lv['basis'],
                         'pcm_lebedev_order': lv['lebedev_order'], 'cm5_alpha': CM5_ALPHA, 'cm5_radii': CM5_RADII,
                         'cm5_D': CM5_D, 'cm5_pair': {f'{a}-{b}': v for (a, b), v in CM5_PAIR.items()},
                         'multiwfn': mwfn}
    job.save()
    case = rc.load_case(f'{args.case}_gas')
    frags = rc.fragments(case)
    labels = None
    facts = facts_values() if args.case == 'btpc' else {}
    ok = True
    from pyscf import dft
    from pyscf.lib import chkfile
    for label, kind, path, model in expand_sources(args.source, args.case, args.smoke):
        key = f'{label}'
        if job.done(key):
            job.log(f'{key}: done; skipped')
            continue
        if path is None or not pathlib.Path(path).exists():
            raise FileNotFoundError(f'source {label}: {path} not found')
        job.input_file(f'source_{label}', path)
        rec = {'source': kind, 'path': str(path), 'model': model}
        t0 = time.time()
        if kind == 'chk':
            mol = chkfile.load_mol(str(path))
            mol.max_memory = rc.MAX_MEMORY_MB
            if np.ndim(mol._ecpbas) != 2:          # a Mole restored from a chk keeps an empty 1-D ECP table
                mol._ecpbas = np.zeros((0, 8), dtype=np.int32)
            scf_rec = chkfile.load(str(path), 'scf')
            C, occ, e = scf_rec['mo_coeff'], scf_rec['mo_occ'], scf_rec['mo_energy']
            dm = (np.einsum('xpi,xi,xqi->xpq', C, occ, C) if np.ndim(occ) == 2 else (C * occ) @ C.T)
            rec['E_Ha_in_chk'] = float(scf_rec.get('e_tot', np.nan))
            # the case's stored geometry for comparison, if this is the stored wavefunction
            ref = rc.ROOT / rc.STORED.get(f'{args.case}_{model}', {}).get('xyz', '') if model in ('gas', 'pcm') else None
            if ref and ref.is_file():
                s0, x0 = rc.read_xyz(ref)
                rec['max_dev_from_opt_xyz_A'] = float(np.abs(mol.atom_coords() * rc.BOHR_A - x0).max())
            grids = dft.gen_grid.Grids(mol); grids.level = lv['grid_level']; grids.build(with_non0tab=True)
            mo_for_molden = (C, e, occ) if np.ndim(occ) == 1 else None
        else:   # new single point
            mol = rc.mol_for_case(rc.load_case(f'{args.case}_{model}'), path, basis=lv['basis'], verbose=3,
                                  output=job.out / f'pyscf_{label}.log')
            mf = rc.make_mf(mol, xc=args.xc, solvent=None if model == 'gas' else 'water', conv_tol=lv['conv_tol'],
                            grid_level=lv['grid_level'], lebedev_order=lv['lebedev_order'],
                            chkfile=job.chk_dir / f'{label}.chk')
            rec['settings'] = rc.settings_of(mf)
            wall = rc.run_scf(mf, chk_guess=job.chk_dir / f'{label}.chk', log=job.log)
            rec.update({'E_Ha': float(mf.e_tot), 'scf_wall_s': round(wall, 1), 'scf_cycles': mf._cycles})
            dm = mf.make_rdm1(); grids = mf.grids
            mo_for_molden = (mf.mo_coeff, mf.mo_energy, mf.mo_occ)
        free = free_atom_dms(mol, xc=args.xc)
        res = full_analysis(mol, dm, grids, free, frags, case)
        rec.update({'natm': mol.natm, 'nao': int(mol.nao), 'n_grid_points': int(grids.weights.size),
                    'basis': mol.basis if isinstance(mol.basis, str) else str(mol.basis), **res,
                    'analysis_wall_s': round(time.time() - t0, 1)})
        labels = hl.labels_of(case, mol.natm)
        rec['labels'] = labels
        # check of the wavefunction source against the paper's Mulliken/Loewdin values
        chk_rows = {}
        for m in ('mulliken', 'lowdin'):
            fk = f'{label}.{m}_charge_picrate'
            if fk in facts:
                chk_rows[fk] = {'facts': facts[fk], 'here': res[m]['fragments']['anion'],
                                'diff': res[m]['fragments']['anion'] - facts[fk]}
            fk = f'{label}.{m}_charge_H4N'
            if fk in facts and 'proton_H' in res:
                chk_rows[fk] = {'facts': facts[fk], 'here': res['proton_H'][m], 'diff': res['proton_H'][m] - facts[fk]}
        if chk_rows:
            rec['check_vs_facts'] = chk_rows
            worst = max(abs(v['diff']) for v in chk_rows.values())
            rec['check_vs_facts_pass'] = worst < 1e-4
            ok &= rec['check_vs_facts_pass']
            job.log(f"{key}: Mulliken/Loewdin vs facts.json: max |diff| {worst:.1e} e -> "
                    f"{'PASS' if worst < 1e-4 else 'FAIL'}")
        # Multiwfn comparison on a molden of the same wavefunction
        if mwfn and mo_for_molden is not None:
            molden_path = job.out / f'{label}.molden'
            stored_molden = pathlib.Path(path).with_name('final.molden') if kind == 'chk' else None
            if stored_molden and stored_molden.exists():
                molden_path = stored_molden
            else:
                write_molden(mol, *mo_for_molden, molden_path)
            try:
                mw = multiwfn_hirshfeld_cm5(mwfn, molden_path, job.out / f'mwfn_{label}', rc.THREADS)
                q_h = np.array(res['hirshfeld']['atoms']); q_c = np.array(res['cm5']['atoms'])
                corr_mw = np.array(mw['cm5']) - np.array(mw['hirshfeld'])
                mw.update({'molden': str(molden_path),
                           'max_abs_diff_hirshfeld': float(np.abs(np.array(mw['hirshfeld']) - q_h).max()),
                           'max_abs_diff_cm5': float(np.abs(np.array(mw['cm5']) - q_c).max()),
                           'max_abs_diff_cm5_correction': float(np.abs(corr_mw - np.array(res['cm5']['correction'])).max()),
                           'fragments_hirshfeld': {f: float(np.sum(np.array(mw['hirshfeld'])[a])) for f, a in frags.items()},
                           'fragments_cm5': {f: float(np.sum(np.array(mw['cm5'])[a])) for f, a in frags.items()}})
                rec['multiwfn'] = mw
                job.log(f"{key}: Multiwfn vs PySCF: Hirshfeld max |diff| {mw['max_abs_diff_hirshfeld']:.4f}, "
                        f"CM5 correction max |diff| {mw['max_abs_diff_cm5_correction']:.1e} e")
                if mw['max_abs_diff_cm5_correction'] > 1e-4:
                    job.log(f'{key}: WARNING CM5 correction differs from Multiwfn by more than 1e-4 e')
                    ok = False
            except Exception as ex:
                rec['multiwfn'] = {'error': f'{type(ex).__name__}: {ex}'}
                job.log(f'{key}: Multiwfn comparison failed: {ex}')
        fa = {m: rec[m]['fragments'].get('anion') for m in ('mulliken', 'lowdin', 'hirshfeld', 'cm5')}
        job.log(f"{key}: picrate charge  Mulliken {fa['mulliken']:+.4f}  Loewdin {fa['lowdin']:+.4f}  "
                f"Hirshfeld {fa['hirshfeld']:+.4f}  CM5 {fa['cm5']:+.4f} e  (grid {rec['n_grid_points']} points, "
                f"N integrated {res['hirshfeld']['n_electrons_integrated']:.5f})")
        job.put(key, rec)
    job.r['summary'] = {k: {m: v[m]['fragments'] for m in ('mulliken', 'lowdin', 'hirshfeld', 'cm5')} |
                        ({'proton_H': v['proton_H']} if 'proton_H' in v else {})
                        for k, v in job.r['items'].items() if v.get('status') == 'done'}
    job.save()
    return ok


if __name__ == '__main__':
    rc.install_sigterm_handler()
    rc.main_wrapper(parser, run, __file__)
