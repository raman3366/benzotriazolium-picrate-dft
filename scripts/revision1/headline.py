"""
Headline tests for the BTPC revision (plan section 5, item 3): where does the LUMO live, and is it the geometry or the
dielectric field that moves it?

  (i)  cross single points: the gas-phase minimum in IEF-PCM water and the IEF-PCM minimum in the gas phase;
  (ii) a permittivity series eps = 1, 2, 4, 8, 20, 78.3553 at both minima (eps = 1 is the gas-phase calculation,
       eps = 78.3553 is water), which contains the cross points and the two stored calculations.
For every point: total energy, HOMO-1, HOMO, LUMO, LUMO+1..LUMO+3 energies and their cation/picrate shares (Mulliken
share as dft/fmo.py defines it, plus Loewdin), the gap, dipole, and fragment charges.  For each geometry the eps at
which the LUMO switches from the cation to the picrate is found from the sign change of
    Delta(eps) = E(lowest cation-localised virtual) - E(lowest picrate-localised virtual)      ("localised": share > 50 %)
by linear interpolation in the Pekar-type factor f = (eps - 1)/eps (also given: interpolation in ln eps), and refined
with --refine extra single points at the interpolated eps.

Level of theory (rev_common.py): B3LYP/6-311++G(d,p), density fitting, grid level 4, IEF-PCM Lebedev order 29,
conv_tol 1e-10.  The (gas geometry, eps 1) and (PCM geometry, eps 78.3553) points must reproduce the stored energies
to 1e-6 Ha (B3LYP, btpc only); otherwise the job exits non-zero after finishing all points.

Usage:
  python dft/revision/headline.py --job headline_gas --out <dir> [--geom gas,pcm] [--eps 1,2,4,8,20,78.3553]
         [--xc b3lyp] [--disp d3bj] [--case btpc] [--refine 1] [--guess auto|minao] [--smoke]
  --geom entries: 'gas' / 'pcm' (the case's stored minimum, molecules/<case>/results/<case>_<model>/opt.xyz) or
         'label=path/to/geometry.xyz' (path relative to the repository or absolute; a directory means <dir>/opt.xyz).
  --xc pbe0 / camb3lyp with --eps 1,78.3553 gives the functional check of plan item 8 with the same analysis.
--smoke: STO-3G, grid 2, Lebedev 17, conv_tol 1e-7, eps 1,4,78.3553 -- exercises every code path in minutes.

This module also holds helpers shared by the other builder-2 scripts (frontier-orbital analysis, fmo.py descriptors,
geometry metrics); they are plain functions and run nothing on import.
"""
import sys, pathlib, argparse, json, math, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import rev_common as rc          # first: pins the thread pools
import numpy as np

STORED_TOL_HA = 1e-6
WATER = rc.SOLVENT_EPS['water']
DEBYE = 2.541746473              # a.u. of dipole -> Debye


# ============================================================================ shared helpers (used by other scripts)
def level(smoke, grid=None):
    """Basis, grid, Lebedev order and SCF tolerance: production values, or the reduced --smoke set."""
    if smoke:
        return {'basis': 'sto-3g', 'grid_level': grid or 2, 'lebedev_order': 17, 'conv_tol': 1e-7}
    return {'basis': rc.BASIS, 'grid_level': rc.GRID_LEVEL, 'lebedev_order': rc.PCM_LEBEDEV_ORDER,
            'conv_tol': rc.CONV_TOL_SP}


def resolve_geom(spec, case_prefix='btpc'):
    """'gas' | 'pcm' | 'label=path' | path  ->  (label, absolute path, model of the stored minimum or None)."""
    if spec in ('gas', 'pcm'):
        name = f'{case_prefix}_{spec}'
        if name in rc.STORED:
            return spec, rc.ROOT / rc.STORED[name]['xyz'], spec
        mol_dir = case_prefix
        return spec, rc.ROOT / 'molecules' / mol_dir / 'results' / name / 'opt.xyz', spec
    label, _, path = spec.partition('=') if '=' in spec else (pathlib.Path(spec).stem, '', spec)
    p = pathlib.Path(path)
    if not p.is_absolute():
        p = rc.ROOT / p
    if p.is_dir():
        p = p / 'opt.xyz'
    return label, p, None


def stored_dir(case_name):
    """data/r2/molecules/<id>/results/<case>/ (raw outputs: final.chk, hessian.npy, final.molden) or None."""
    for d in sorted((rc.ROOT / 'data' / 'r2' / 'molecules').glob(f'*/results/{case_name}')):
        return d
    return None


def labels_of(case, natm):
    lab = case.get('cif_labels')
    if lab:
        return [l.split('#')[0] for l in lab]
    return [f'X{i + 1}' for i in range(natm)]


def frag_of_share(shares, thresh=0.5):
    """Fragment holding more than `thresh` of an orbital ({'anion': x, 'cation': y} in percent), else 'mixed'."""
    best = max(shares, key=shares.get)
    return best if shares[best] > 100 * thresh else 'mixed'


def frontier(mf, frags, n_occ=2, n_virt=4, scan=12):
    """Frontier-orbital energies (eV) and fragment shares (percent; Mulliken as dft/fmo.py, and Loewdin) of an RKS
    solution: HOMO-(n_occ-1)..HOMO, LUMO..LUMO+(n_virt-1).  Also the lowest virtual and highest occupied orbital
    localised (> 50 % Mulliken) on each fragment, searched over `scan` orbitals."""
    mol = mf.mol
    S = mol.intor_symmetric('int1e_ovlp'); S_half = rc.lowdin_half(S)
    e = np.asarray(mf.mo_energy) * rc.HARTREE_EV
    C = np.asarray(mf.mo_coeff)
    nocc = int((np.asarray(mf.mo_occ) > 0.5).sum())
    nmo = C.shape[1]

    def shares(k):
        return {m: {f: round(100 * rc.mo_fragment_share(mol, C[:, k], atoms, m, S, S_half), 3)
                    for f, atoms in frags.items()} for m in ('mulliken', 'lowdin')}

    out = {'n_occ': nocc, 'n_mo': nmo, 'orbitals': {}}
    labels = [(f'HOMO-{i}' if i else 'HOMO', nocc - 1 - i) for i in range(n_occ - 1, -1, -1)]
    labels += [(f'LUMO+{i}' if i else 'LUMO', nocc + i) for i in range(n_virt)]
    for lab, k in labels:
        if 0 <= k < nmo:
            sh = shares(k)
            out['orbitals'][lab] = {'index': int(k), 'eV': float(e[k]), **sh,
                                    'on': frag_of_share(sh['mulliken'])}
    out['HOMO_eV'] = float(e[nocc - 1]); out['LUMO_eV'] = float(e[nocc]); out['gap_eV'] = float(e[nocc] - e[nocc - 1])
    out['LUMO_on'] = out['orbitals']['LUMO']['on']; out['HOMO_on'] = out['orbitals']['HOMO']['on']
    lowest_virt, highest_occ = {}, {}
    for k in range(nocc, min(nmo, nocc + scan)):
        sh = shares(k)['mulliken']; f = frag_of_share(sh)
        if f in frags and f not in lowest_virt:
            lowest_virt[f] = {'index': int(k), 'label': 'LUMO' + (f'+{k - nocc}' if k > nocc else ''),
                              'eV': float(e[k]), 'share_pct': sh[f]}
    for k in range(nocc - 1, max(-1, nocc - 1 - scan), -1):
        sh = shares(k)['mulliken']; f = frag_of_share(sh)
        if f in frags and f not in highest_occ:
            highest_occ[f] = {'index': int(k), 'label': 'HOMO' + (f'-{nocc - 1 - k}' if k < nocc - 1 else ''),
                              'eV': float(e[k]), 'share_pct': sh[f]}
    out['lowest_virtual_on'] = lowest_virt
    out['highest_occupied_on'] = highest_occ
    if 'cation' in lowest_virt and 'anion' in lowest_virt:
        out['delta_virt_cation_minus_anion_eV'] = lowest_virt['cation']['eV'] - lowest_virt['anion']['eV']
    return out


def fmo_descriptors(H, L):
    """Koopmans-type descriptors exactly as dft/fmo.py computes them (eV)."""
    IP, EA = -H, -L
    chi = (IP + EA) / 2; mu = -chi; eta = (IP - EA) / 2; S = 1 / (2 * eta); omega = mu ** 2 / (2 * eta)
    return {'HOMO_eV': H, 'LUMO_eV': L, 'gap_eV': L - H, 'IP_eV': IP, 'EA_eV': EA, 'electronegativity_chi_eV': chi,
            'chemical_potential_mu_eV': mu, 'hardness_eta_eV': eta, 'softness_S_eV-1': S,
            'electrophilicity_omega_eV': omega}


def populations(mol, dm, frags):
    """Mulliken and Loewdin charges per atom and per fragment (and spin for a UKS density pair)."""
    S = mol.intor_symmetric('int1e_ovlp'); S_half = rc.lowdin_half(S)
    dm = np.asarray(dm)
    uks = dm.ndim == 3
    tot = dm[0] + dm[1] if uks else dm
    Z = mol.atom_charges()
    sl = mol.aoslice_by_atom()
    res = {}
    for name in ('mulliken', 'lowdin'):
        if name == 'mulliken':
            p = np.einsum('ij,ji->i', tot, S); s = np.einsum('ij,ji->i', dm[0] - dm[1], S) if uks else None
        else:
            p = np.einsum('ij,jk,ki->i', S_half, tot, S_half); s = (np.einsum('ij,jk,ki->i', S_half, dm[0] - dm[1], S_half)
                                                                   if uks else None)
        q_at = np.array([Z[a] - p[sl[a][2]:sl[a][3]].sum() for a in range(mol.natm)])
        r = {'atoms': q_at.tolist(), 'fragments': {f: float(q_at[atoms].sum()) for f, atoms in frags.items()}}
        if uks:
            s_at = np.array([s[sl[a][2]:sl[a][3]].sum() for a in range(mol.natm)])
            r['spin_atoms'] = s_at.tolist()
            r['spin_fragments'] = {f: float(s_at[atoms].sum()) for f, atoms in frags.items()}
        res[name] = r
    return res


def dipole_debye(mf):
    v = mf.dip_moment(unit='Debye', verbose=0)
    return {'total': float(np.linalg.norm(v)), 'vector': [float(x) for x in v]}


def _dihedral(X, a, b, c, d):
    b0 = X[a] - X[b]; b1 = X[c] - X[b]; b2 = X[d] - X[c]
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n; w = b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def _plane_normal(P):
    c = P.mean(0)
    return np.linalg.svd(P - c)[2][-1]


def nearest(X, syms, i, element, among):
    """Index of the nearest atom of `element` to atom i, among the indices `among`."""
    cand = [j for j in among if syms[j] == element and j != i]
    return min(cand, key=lambda j: np.linalg.norm(X[j] - X[i])) if cand else None


def hbond_torsion_atoms(case, syms, X):
    """C(acceptor)-O(acceptor)...N(donor)-C(donor ring): for btpc C1-O1...N4-C12 (indices 4, 0, 18, 21), the
    inter-ion torsion of revision/round1/audit/evidence/dft/dihedrals.py."""
    p = case.get('proton')
    if not p:
        return None
    fr = rc.fragments(case)
    O, N = p['acceptor_O'], p['donor_N']
    cO = nearest(X, syms, O, 'C', fr.get('anion', range(len(syms))))
    cN = nearest(X, syms, N, 'C', fr.get('cation', range(len(syms))))
    return (cO, O, N, cN) if None not in (cO, cN) else None


def geometry_metrics(case, syms, X):
    """Inter-ion geometry: H-bond distances and angle, proton position, the inter-ion torsion and the angle between
    the picrate ring (anion carbons) and the benzotriazolium ring system (cation heavy atoms)."""
    X = np.asarray(X, float)
    lab = labels_of(case, len(syms))
    out = {}
    p = case.get('proton')
    if p:
        H, N, O = p['H'], p['donor_N'], p['acceptor_O']
        dHN = float(np.linalg.norm(X[H] - X[N])); dHO = float(np.linalg.norm(X[H] - X[O]))
        v1, v2 = X[N] - X[H], X[O] - X[H]
        out.update({'d_N...O_A': float(np.linalg.norm(X[N] - X[O])), 'd_H-N_A': dHN, 'd_H...O_A': dHO,
                    'angle_N-H...O_deg': float(np.degrees(np.arccos(np.dot(v1, v2) / np.linalg.norm(v1) /
                                                                     np.linalg.norm(v2)))),
                    'proton_transferred': dHO < dHN,
                    'hbond_labels': f'{lab[N]}-{lab[H]}...{lab[O]}'})
    t = hbond_torsion_atoms(case, syms, X)
    if t:
        out['torsion_atoms'] = [int(i) for i in t]
        out['torsion_labels'] = '-'.join(lab[i] for i in t)
        out['torsion_deg'] = _dihedral(X, *t)
    fr = rc.fragments(case)
    if 'anion' in fr and 'cation' in fr:
        ring_a = [i for i in fr['anion'] if syms[i] == 'C']
        ring_c = [i for i in fr['cation'] if syms[i] != 'H']
        n1, n2 = _plane_normal(X[ring_a]), _plane_normal(X[ring_c])
        out['inter_ring_angle_deg'] = float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(n1, n2)))))))
    return out


def dihedral(X, idx):
    return _dihedral(np.asarray(X, float), *idx)


def mol_geometry(mol):
    syms = [mol.atom_pure_symbol(i) for i in range(mol.natm)]
    return syms, mol.atom_coords() * rc.BOHR_A


def solvent_label(eps):
    return 'gas' if eps is None else f'IEF-PCM eps={eps:g}'


def solvation_energy(mf):
    ws = getattr(mf, 'with_solvent', None)
    return None if ws is None or getattr(ws, 'e', None) is None else float(ws.e)


# ============================================================================ this script
def parser():
    p = rc.add_common_args(argparse.ArgumentParser(description=__doc__.split('\n')[1]))
    p.add_argument('--geom', default='gas,pcm', help="comma list: gas, pcm, label=path.xyz")
    p.add_argument('--eps', default=None, help='comma list of permittivities (1 = gas phase); default 1,2,4,8,20,78.3553')
    p.add_argument('--xc', default=rc.XC, help='functional (default b3lyp; item 8: pbe0, camb3lyp)')
    p.add_argument('--disp', default=None, help="dispersion correction: none (default) or d3bj")
    p.add_argument('--case', default='btpc', help='molecule prefix: btpc (default) or btpcz')
    p.add_argument('--refine', type=int, default=1, help='extra single points at the interpolated switch eps (default 1)')
    p.add_argument('--guess', choices=('auto', 'minao'), default='auto',
                   help='auto: stored final.chk of the same model for the first point, then the previous eps')
    p.add_argument('--test-delta-offset', type=float, default=0.0,
                   help='TEST HOOK, --smoke only: add this many eV to Delta so that the switch/refine path runs')
    return p


def switch_from(points):
    """points: [(eps, delta_eV)] sorted by eps.  Returns the bracket and interpolated eps of the first sign change."""
    pts = [(e, d) for e, d in points if d is not None]
    pts.sort()
    for (e0, d0), (e1, d1) in zip(pts, pts[1:]):
        if d0 == 0:
            return {'eps_switch_f': e0, 'eps_switch_ln': e0, 'bracket': [e0, e0]}
        if d0 * d1 < 0:
            f = lambda e: (e - 1) / e
            x = f(e0) + (0 - d0) * (f(e1) - f(e0)) / (d1 - d0)
            l = math.log(e0) + (0 - d0) * (math.log(e1) - math.log(e0)) / (d1 - d0)
            return {'eps_switch_f': 1 / (1 - x), 'eps_switch_ln': math.exp(l), 'bracket': [e0, e1],
                    'delta_at_bracket_eV': [d0, d1],
                    'direction': 'cation->anion' if d0 < 0 else 'anion->cation'}
    return None


def run(job, args):
    lv = level(args.smoke)
    if args.test_delta_offset and not args.smoke:
        raise SystemExit('--test-delta-offset is a test hook for --smoke runs only')
    disp = None if (args.disp in (None, '', 'none')) else args.disp
    eps_list = ([float(x) for x in args.eps.split(',')] if args.eps else
                ([1.0, 4.0, WATER] if args.smoke else [1.0, 2.0, 4.0, 8.0, 20.0, WATER]))
    eps_list = sorted(set(eps_list))
    case = rc.load_case(f'{args.case}_gas')
    frags = rc.fragments(case)
    job.r['settings'] = {**lv, 'xc': args.xc, 'disp': disp, 'pcm_method': rc.PCM_METHOD, 'eps_list': eps_list,
                         'density_fit': True, 'auxbasis': 'PySCF default', 'max_cycle': rc.MAX_CYCLE,
                         'eps_1_means': 'gas-phase calculation (no PCM)', 'fragment_share': 'Mulliken (as dft/fmo.py) '
                         'and Loewdin; an orbital is localised on a fragment when its Mulliken share > 50 %',
                         'switch_interpolation': 'linear in f=(eps-1)/eps (also ln eps)', 'refine': args.refine,
                         'case': case['name'], 'fragments': frags}
    job.input_file('case', rc.find_case_file(case['name']))
    job.save()
    check_stored = (not args.smoke and args.xc.lower() == 'b3lyp' and disp is None and args.case == 'btpc')
    ok = True

    def point(glabel, gpath, model, eps, guess_chk):
        key = f'{glabel}|eps={eps:g}'
        chk = job.chk_dir / f'{glabel}_eps{eps:g}.chk'
        if job.done(key):
            return job.get(key), chk
        mol = rc.mol_for_case(case, gpath, basis=lv['basis'], verbose=3 if args.smoke else 4,
                              output=job.out / f'pyscf_{glabel}.log')
        solvent = None if eps == 1.0 else eps
        mf = rc.make_mf(mol, xc=args.xc, disp=disp, solvent=solvent, conv_tol=lv['conv_tol'],
                        grid_level=lv['grid_level'], lebedev_order=lv['lebedev_order'], chkfile=chk)
        rec = {'geometry': glabel, 'xyz': str(gpath), 'eps': eps, 'model': solvent_label(solvent),
               'settings': rc.settings_of(mf)}
        job.put(key, rec, status='running')
        g = chk if chk.exists() else guess_chk
        job.log(f'{key}: SCF ({mol.nao} AO, {rec["model"]}), guess {g if g and pathlib.Path(g).exists() else "minao"}')
        wall = rc.run_scf(mf, chk_guess=g, log=job.log)
        fr = frontier(mf, frags)
        if args.test_delta_offset and 'delta_virt_cation_minus_anion_eV' in fr:
            fr['delta_virt_cation_minus_anion_eV'] += args.test_delta_offset
            fr['TEST_delta_offset_eV'] = args.test_delta_offset
        rec.update({'E_Ha': float(mf.e_tot), 'E_solvation_Ha': solvation_energy(mf), 'scf_cycles': mf._cycles,
                    'scf_wall_s': round(wall, 1), **rc.runtime_sizes(mf), 'frontier': fr,
                    'dipole_debye': dipole_debye(mf),
                    'charges': {k: v['fragments'] for k, v in populations(mol, mf.make_rdm1(), frags).items()}})
        if disp:
            rec['E_dispersion_Ha'] = float(np.asarray(mf.scf_summary.get('dispersion', np.nan)))
        job.r['timings'][key] = round(wall, 1)
        orb = fr['orbitals']
        job.log(f"{key}: E = {mf.e_tot:.10f} Ha ({mf._cycles} cycles, {wall:.0f} s); HOMO {fr['HOMO_eV']:.4f} "
                f"({fr['HOMO_on']}), LUMO {fr['LUMO_eV']:.4f} ({fr['LUMO_on']}, cation "
                f"{orb['LUMO']['mulliken'].get('cation', float('nan')):.1f} %), gap {fr['gap_eV']:.4f} eV; "
                f"Delta(cat-an) {fr.get('delta_virt_cation_minus_anion_eV', float('nan')):+.4f} eV")
        job.put(key, rec)
        mol.stdout.close()
        return rec, chk

    for spec in [s.strip() for s in args.geom.split(',') if s.strip()]:
        glabel, gpath, model = resolve_geom(spec, args.case)
        if not gpath.exists():
            raise FileNotFoundError(f'geometry {gpath} for --geom {spec}')
        job.input_file(f'geometry_{glabel}', gpath)
        stored_chk = None
        if args.guess == 'auto' and not args.smoke and model:
            d = stored_dir(f'{args.case}_{model}')
            if d and (d / 'final.chk').exists():
                stored_chk = d / 'final.chk'
        prev = stored_chk
        series = []
        for eps in eps_list:
            rec, chk = point(glabel, gpath, model, eps, prev)
            prev = chk
            series.append(rec)
            ref_key = (f'{args.case}_gas' if (model == 'gas' and eps == 1.0) else
                       f'{args.case}_pcm' if (model == 'pcm' and abs(eps - WATER) < 1e-9) else None)
            if check_stored and ref_key in rc.STORED:
                dE = rec['E_Ha'] - rc.STORED[ref_key]['E_Ha']
                rec['dE_vs_stored_Ha'] = dE
                rec['reproduces_stored'] = abs(dE) <= STORED_TOL_HA
                job.put(f'{glabel}|eps={eps:g}', rec)
                job.log(f"{glabel} eps {eps:g}: dE vs stored {ref_key} = {dE:+.2e} Ha -> "
                        f"{'PASS' if rec['reproduces_stored'] else 'FAIL'}")
                ok &= rec['reproduces_stored']
        # LUMO switch and refinement
        pts = [(r['eps'], r['frontier'].get('delta_virt_cation_minus_anion_eV')) for r in series]
        sw = switch_from(pts)
        for it in range(args.refine if sw else 0):
            e_new = float(f"{sw['eps_switch_f']:.4g}")
            if any(abs(e_new - e) < 1e-6 for e, _ in pts):
                break
            # guess: the converged density of the nearest eps in the series
            near = min(series, key=lambda r: abs(r['eps'] - e_new))
            rec, _ = point(glabel, gpath, model, e_new, job.chk_dir / f"{glabel}_eps{near['eps']:g}.chk")
            rec['refinement'] = it + 1
            job.put(f'{glabel}|eps={e_new:g}', rec)
            pts.append((e_new, rec['frontier'].get('delta_virt_cation_minus_anion_eV')))
            sw = switch_from(pts) or sw
        lumo_on = {f"{e:g}": job.get(f'{glabel}|eps={e:g}')['frontier']['LUMO_on'] for e, _ in sorted(pts)}
        job.r.setdefault('switch', {})[glabel] = sw or {'eps_switch_f': None, 'note': 'no sign change of Delta in the '
                                                         f'eps range {min(eps_list):g}-{max(eps_list):g}',
                                                         'LUMO_on': lumo_on}
        if sw:
            sw['LUMO_on'] = lumo_on
            job.log(f"{glabel}: LUMO switches {sw['direction']} at eps ~ {sw['eps_switch_f']:.3f} (f-interpolation; "
                    f"ln-eps {sw['eps_switch_ln']:.3f}), bracket {sw['bracket']}")
        else:
            job.log(f'{glabel}: no LUMO switch in the eps range; LUMO on {lumo_on}')
        job.save()

    # summary tables
    summ = {}
    for key, it in sorted(job.r['items'].items(), key=lambda kv: (kv[1].get('geometry', ''), kv[1].get('eps', 0))):
        if it.get('status') != 'done':
            continue
        fr = it['frontier']; o = fr['orbitals']
        summ.setdefault(it['geometry'], []).append({
            'eps': it['eps'], 'E_Ha': it['E_Ha'], 'HOMO_eV': fr['HOMO_eV'], 'LUMO_eV': fr['LUMO_eV'],
            'gap_eV': fr['gap_eV'], 'LUMO_on': fr['LUMO_on'],
            'LUMO_cation_pct': o['LUMO']['mulliken'].get('cation'),
            **{f'{k}_eV': o[k]['eV'] for k in ('LUMO+1', 'LUMO+2', 'LUMO+3') if k in o},
            'delta_cat_minus_an_eV': fr.get('delta_virt_cation_minus_anion_eV'),
            'dE_vs_stored_Ha': it.get('dE_vs_stored_Ha')})
    job.r['summary'] = {'series': summ, 'switch': job.r.get('switch'), 'stored_reproduced': ok if check_stored else None}
    job.save()
    return ok


if __name__ == '__main__':
    rc.install_sigterm_handler()
    rc.main_wrapper(parser, run, __file__)
