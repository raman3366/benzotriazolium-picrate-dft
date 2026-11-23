"""Revision facts: the proton position (revision item 2 and the 'Option A' CCSD(T) correction, Nov 2026).

Inputs, read-only (data/revision_runs/queue*/out/<job>/result.json; job purposes in dft/revision/jobs_vm*.json):
  scan_gas_up, scan_gas_down, scan_gas_neutral, scan_gas_oside   B3LYP/6-311++G(d,p) relaxed gas-phase scan of
                                                                 r(H4N...O1) (dft/revision/proton_scan.py)
  sp_scan_{wb97xd3bj,m062x,pbe0,camb3lyp,dfmp2}, sp_oside_all    single points at those geometries (singlepoints.py)
  dvr_gas_all (dvr_gas_b3lyp: first run, superseded)             1-D proton levels on every profile (dvr.py)
  neutral_start_gas, neutral_start_pcm                           B3LYP optimisations started from the neutral form
  sp_min_pcm, sp_neutral_pcm                                     IEF-PCM single points at the two IEF-PCM geometries
  sp_gas_pcmgeoms                                                gas-phase DF-MP2 at the same two geometries
  cc_*                                                           FNO-CCSD(T) - MP2 corrections (cc_correction.py)
Outputs: facts rev.proton.* and, only when called from build_facts.py (not from `python -m rev.common proton`),
papers/btpc/data/derived/rev_proton_{scan_profiles,dvr_levels,dvr_grid,pcm,cc}.csv.

The DVR ground-state densities are not stored in dvr_gas_all/result.json; they are recomputed here from the stored
profiles with the functions of dft/revision/dvr.py (the script, same SHA-256, that produced dvr_gas_all) at the
stored settings, and every recomputed level and O-side weight is checked against the stored ones before use.
"""
from .common import *
import os, csv, types, hashlib, importlib.util
import numpy as np

KC = HA_KCAL
DERIVED = ROOT / 'papers/btpc/data/derived'
SCAN_JOBS = ('scan_gas_up', 'scan_gas_down', 'scan_gas_neutral', 'scan_gas_oside')
SP_SCAN = {'wb97x-d3bj': 'sp_scan_wb97xd3bj', 'm062x': 'sp_scan_m062x', 'pbe0': 'sp_scan_pbe0',
           'cam-b3lyp': 'sp_scan_camb3lyp', 'dfmp2': 'sp_scan_dfmp2'}
SP_OSIDE = 'sp_oside_all'
# (method id in the job files, key slug, name)
METHODS = [('b3lyp', 'b3lyp', 'B3LYP'), ('wb97x-d3bj', 'wb97xd3bj', 'wB97X-D3(BJ)'), ('m062x', 'm062x', 'M06-2X'),
           ('pbe0', 'pbe0', 'PBE0'), ('cam-b3lyp', 'camb3lyp', 'CAM-B3LYP'), ('dfmp2', 'dfmp2', 'DF-MP2'),
           ('scs-mp2', 'scsmp2', 'SCS-MP2')]
SLUG = {m: s for m, s, _ in METHODS}
NAME = {m: n for m, _, n in METHODS}
GEOM_NOTE = ('single point at the B3LYP/6-311++G(d,p) relaxed-scan geometry (6-311++G(d,p), density fitting, '
             'grid level 4; DF-MP2 frozen core, RI basis aug-cc-pVTZ-RI; SCS-MP2 = 6/5 OS + 1/3 SS of the same job)')
CC_RUNS = [  # run id, model, basis, frac, salt (job, item), neutral (job, item), DF-MP2 reference
    ('pcm_631gd_f050', 'pcm', '6-31G*', 0.50, ('cc_salt_631gd', 'salt|f0.500'), ('cc_neutral_631gd', 'neutral|f0.500')),
    ('pcm_631gd_f065', 'pcm', '6-31G*', 0.65, ('cc_salt_631gd', 'salt|f0.650'), ('cc_neutral_631gd_f065', 'neutral|f0.650')),
    ('pcm_631gd_f080', 'pcm', '6-31G*', 0.80, ('cc_salt_631gd_f080', 'salt|f0.800'), ('cc_neutral_631gd_f080', 'neutral|f0.800')),
    ('pcm_631ppgd_f050', 'pcm', '6-31+G(d,p)', 0.50, ('cc_salt_631ppgd', 'salt|f0.500'), ('cc_neutral_631ppgd', 'neutral|f0.500')),
    ('gas_631gd_f050', 'gas', '6-31G*', 0.50, ('cc_gas_631gd', 'salt|f0.500'), ('cc_gas_631gd', 'neutral|f0.500')),
]
DVR_VARIANTS = {('xH', 'H', 'cubic', 'extrapolate'): 'main', ('xH', 'H', 'pchip', 'extrapolate'): 'pchip',
                ('xH', 'H', 'cubic', 'box'): 'box_walls', ('xH', 'reduced', 'cubic', 'extrapolate'): 'reduced_mass',
                ('half_delta', 'H', 'cubic', 'extrapolate'): 'half_delta', ('rHO', 'H', 'cubic', 'extrapolate'): 'rHO',
                ('mwpath', 'unit (mass-weighted path)', 'cubic', 'extrapolate'): 'mwpath'}

_R = {}


def R(job):
    if job not in _R:
        _R[job] = result(job)
    return _R[job][0]


def S(job, path=''):
    if job not in _R:
        R(job)
    return _R[job][1] + (f':{path}' if path else '')


def _py(v):
    if isinstance(v, dict):
        return {k: _py(x) for k, x in v.items()}
    if isinstance(v, (list, tuple, np.ndarray)):
        return [_py(x) for x in v]
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return float(f'{float(v):.12g}')
    return v


def read_xyz(p):
    L = open(p).read().splitlines(); n = int(L[0].split()[0])
    return [l.split()[0] for l in L[2:2 + n]], np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])


def sha256(p):
    return sha256_as_run(p)   # dataset: checksum as the run wrote the file (rev/common.py)


def _published(facts):
    if 'gas.e_total' in facts:                         # inside build_facts.py
        return lambda k: facts[k]['value']
    pub = json.load(open(ROOT / 'papers/btpc/data/facts.json'))
    return lambda k: pub[k]['value']


def _done(job):
    d = R(job)
    assert d['status'] == 'done' and not d.get('smoke') and not d.get('error'), (job, d['status'], d.get('error'))
    return d


# ------------------------------------------------------------------ geometry helpers (CIF labels of the case file)
def _labels():
    cj = json.load(open(ROOT / 'molecules/btpc/cases/btpc_pcm.json'))
    lab = [l.split('#')[0] for l in cj['cif_labels']]
    return {l: i for i, l in enumerate(lab)}


def _plane_normal(X, idx):
    P = X[idx] - X[idx].mean(0)
    return np.linalg.svd(P)[2][2]


def interplanar(X, L):
    pic = [L[a] for a in ('C1', 'C2', 'C3', 'C4', 'C5', 'C6')]
    cat = [L[a] for a in ('N4', 'N6', 'C7', 'C12', 'N5', 'C11', 'C10', 'C9', 'C8')]
    c = abs(float(_plane_normal(X, pic) @ _plane_normal(X, cat)))
    return float(np.degrees(np.arccos(min(1.0, c))))


def dist(X, L, a, b):
    return float(np.linalg.norm(X[L[a]] - X[L[b]]))


def angle(X, L, a, b, c):
    v1, v2 = X[L[a]] - X[L[b]], X[L[c]] - X[L[b]]
    return float(np.degrees(np.arccos(np.clip(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2), -1, 1))))


# ------------------------------------------------------------------ the gas-phase profile
def scan_points():
    """All finished scan points and the profile: the lowest-energy point per target (as dvr.py and
    singlepoints.py --select lowest), sorted from the N side (large r(H4N...O1)) to the O side."""
    allp = []
    for job in SCAN_JOBS:
        d = _done(job)
        st = d['settings']
        assert st['xc'] == 'b3lyp' and st['disp'] is None and st['basis'] == '6-311++g(d,p)' and st['grid_level'] == 4
        assert st['solvent'] is None and st['labels'] == {'H': 'H4N', 'N': 'N4', 'O': 'O1'}
        for key, it in d['items'].items():
            if it.get('status') != 'done' or it.get('E_Ha') is None:
                continue
            assert it['converged'] is True, (job, key)
            allp.append({'label': key, 'job': job, 'chain': it.get('chain'),
                         't': 'min' if it.get('r_target_A') is None else round(it['r_target_A'], 3),
                         'E': it['E_Ha'], 'r_HO': it['r_HO_A'], 'r_NH': it['r_NH_A'], 'r_NO': it['r_NO_A'],
                         'NHO': it['NHO_deg'], 'steps': it.get('opt_steps', 0),
                         'dev': (it.get('constraint') or {}).get('deviation_A'),
                         'lumo_anion': ((it.get('frontier_composition_mulliken_pct') or {}).get('LUMO') or {}).get('anion')})
    best, by_t = {}, {}
    for p in allp:
        by_t.setdefault(p['t'], []).append(p)
        if p['t'] not in best or p['E'] < best[p['t']]['E']:
            best[p['t']] = p
    prof = sorted(best.values(), key=lambda p: -p['r_HO'])
    return allp, by_t, prof


def sp_energies(prof):
    """E (Ha) per method at the profile geometries: {method: [E per point]}; checks that every single-point job used
    the same geometry (label and r(H4N...O1)) as the profile."""
    E = {'b3lyp': [p['E'] for p in prof]}
    srcs = {'b3lyp': [S(p['job'], f"items.{p['label']}.E_Ha") for p in prof]}
    for m, job in SP_SCAN.items():
        d = _done(job); o = _done(SP_OSIDE)
        geoms = {g['label']: g for g in d['settings']['geometries']}
        geoms_o = {g['label']: g for g in o['settings']['geometries']}
        for mm, fld in ([(m, 'E_Ha')] + ([('scs-mp2', 'E_SCS_MP2_Ha'), ('hf', 'E_HF_Ha')] if m == 'dfmp2' else [])):
            vals, ss = [], []
            for p in prof:
                if p['job'] == 'scan_gas_oside':
                    it, g, j = o['items'][f"{p['label']}|{m}"], geoms_o[p['label']], SP_OSIDE
                else:
                    it, g, j = d['items'][f"{p['label']}|{m}"], geoms[p['label']], job
                assert it['status'] == 'done' and not it.get('newton_fallback'), (j, p['label'], m)
                assert abs(g['r_HO_A'] - p['r_HO']) < 1e-4 and abs(g['r_NH_A'] - p['r_NH']) < 1e-4, (j, p['label'])
                vals.append(it[fld]); ss.append(S(j, f"items.{p['label']}|{m}.{fld}"))
            E[mm] = vals; srcs[mm] = ss
    return E, srcs


# ------------------------------------------------------------------ DVR (recomputed densities)
def _dvr_module():
    env = dict(os.environ)
    spec = importlib.util.spec_from_file_location('rev_dvr_script', paths.resolve('dft/revision/dvr.py'))
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    finally:                                            # dvr.py pins thread variables at import; undo that
        os.environ.clear(); os.environ.update(env)
    return mod


def dvr_recompute(stored):
    """Main-variant DVR (the stored settings) for every method: grid, potential, |psi_n|^2; checked against the
    stored levels (2-decimal cm-1) and O-side weights (5 decimals)."""
    dv = _dvr_module()
    st = stored['settings']
    assert (st['coord'], st['mass'], st['interp'], st['walls']) == ('xH', 'H', 'cubic', 'extrapolate'), st
    assert sha256(paths.resolve('dft/revision/dvr.py')) == stored['script_sha256']['dvr.py'], 'dvr.py changed since dvr_gas_all'
    args = types.SimpleNamespace(scan=[str(job_dir(j)) for j in SCAN_JOBS],
                                 sp=[str(job_dir(j)) for j in list(SP_SCAN.values()) + [SP_OSIDE]],
                                 smoke=False, methods='all', chains='all', npts=st['npts'], pad=st['pad_A'],
                                 nlevels=st['nlevels'])
    profiles, proton = dv.load_profiles(args, lambda *a: None)
    assert proton == (st['proton_atoms_0based']['H'], st['proton_atoms_0based']['donor_N'],
                      st['proton_atoms_0based']['acceptor_O'])
    out = {}
    for m, _, _ in METHODS:
        prof = dv.make_profile(profiles[m]['pts'], proton)
        x = prof['coords']['xH']; V = prof['E'] - prof['E'].min()
        x_b = x / dv.BOHR_A
        g, v, f = dv.build_potential(x_b, V, 'cubic', 'extrapolate', args.pad / dv.BOHR_A, args.npts)
        E, P = dv.solve(g, v, dv.MASS['H'] * dv.AMU_ME, args.nlevels)
        xs, vs, *_ = dv.extrema(f, x_b[0], x_b[-1])
        vmin = float(vs.min())
        o_side = np.interp(g, x_b, prof['delta_A']) > 0
        outside = (g < x_b[0]) | (g > x_b[-1])
        lv = stored['items'][m]['levels']
        for n in range(args.nlevels):
            assert abs((E[n] - vmin) * dv.HARTREE_CM - lv[n]['E_cm']) < 0.006, (m, n)
            assert abs(P[o_side, n].sum() - lv[n]['P_O_side']) < 6e-6, (m, n)
            assert abs(P[outside, n].sum() - lv[n]['P_outside_scan']) < 6e-6, (m, n)
        xp = [p['x_A'] for p in stored['items'][m]['points']]
        assert np.allclose(x, xp, atol=6e-6), m
        out[m] = {'labels': prof['labels'], 'x_pts_A': x, 'r_pts_A': prof['r_HO_A'], 'g_A': g * dv.BOHR_A,
                  'V_kcal': (v - vmin) * KC, 'E_kcal': (E - vmin) * KC, 'E_cm': (E - vmin) * dv.HARTREE_CM,
                  'P': P, 'P_O': P[o_side].sum(0), 'P_out': P[outside].sum(0), 'undershoot_kcal': -vmin * KC}
    return out


def r_on_grid(gx, xp, rp):
    """r(H4N...O1) on the DVR grid: linear in xH between scan points, linear extrapolation in the pads."""
    xp, rp = np.asarray(xp), np.asarray(rp)
    r = np.interp(gx, xp, rp)
    lo, hi = gx < xp[0], gx > xp[-1]
    r[lo] = rp[0] + (gx[lo] - xp[0]) * (rp[1] - rp[0]) / (xp[1] - xp[0])
    r[hi] = rp[-1] + (gx[hi] - xp[-1]) * (rp[-1] - rp[-2]) / (xp[-1] - xp[-2])
    return r


# ------------------------------------------------------------------ CSV writers
def _write_csv(name, header, rows):
    p = DERIVED / name
    tmp = p.with_suffix('.csv.tmp')
    with open(tmp, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)
    os.replace(tmp, p)


def _g(v, nd=10):
    return '' if v is None else (f'{v:.{nd}g}' if isinstance(v, float) else v)


# ======================================================================================== build
def build(add, facts):
    pub = _published(facts)
    writing = 'standalone' not in getattr(add, '__qualname__', '')
    L = _labels()

    def A(key, value, unit, case, source, note):
        add(f'rev.proton.{key}', _py(value), unit, case, source, note)

    # ---------------------------------------------------------------- 1. gas-phase relaxed scan (B3LYP)
    allp, by_t, prof = scan_points()
    n = len(prof)
    lab = [p['label'] for p in prof]
    rHO = np.array([p['r_HO'] for p in prof]); rNH = np.array([p['r_NH'] for p in prof])
    rNO = np.array([p['r_NO'] for p in prof]); delta = rNH - rHO
    scan_src = ', '.join(S(j) for j in SCAN_JOBS)
    lvl = R('scan_gas_up')['settings']
    A('scan.level', f"B3LYP/{lvl['basis'].replace('g(d,p)', 'G(d,p)')}, no dispersion, density fitting, grid level "
      f"{lvl['grid_level']}, Gas phase; every coordinate relaxed except r(H4N...O1), held by a geomeTRIC "
      f"distance constraint; {lvl['opt_criteria_name']}", '', 'rev_gas', scan_src + ':settings',
      'dft/revision/proton_scan.py; four chains: up (1.46-1.70 A) and down (1.40-1.08 A) from the published gas-phase '
      'minimum, neutral (0.98-1.28 A) from a start with H4N on O1, and the O-side extension (0.90, 0.94 A) added on '
      '22 Nov after the first DVR run found the O-side well unbracketed')
    A('scan.coordinate', 'r(H4N...O1)', '', 'rev_gas', scan_src + ':settings.scan_coordinate',
      'the scanned (constrained) coordinate: the H4N...O1 distance; r(N4-H4N), r(N4...O1) and the angle relax')
    A('scan.n_optimisations', sum(1 for p in allp if p['chain'] != 'min'), '', 'rev_gas', scan_src + ':items',
      'constrained optimisations in the four scan jobs (every one converged on the Gaussian default criteria); '
      'plus the single point at the published minimum (chain min)')
    A('scan.n_points', n, '', 'rev_gas', scan_src + ':items',
      'points of the profile: one per target distance (the lower-energy chain where two chains reached the same '
      'target) plus the published gas-phase minimum')
    A('scan.labels', lab, '', 'rev_gas', scan_src + ':items', 'profile points, N side (large r(H4N...O1)) first; '
      'the order of every rev.proton.scan.* and rev.proton.gas.*.dE_kcal list')
    A('scan.jobs', [p['job'] for p in prof], '', 'rev_gas', scan_src, 'job of each profile point')
    for k, arr, u, what in (('r_HO_A', rHO, 'Å', 'r(H4N...O1)'), ('r_NH_A', rNH, 'Å', 'r(N4-H4N)'),
                            ('r_NO_A', rNO, 'Å', 'r(N4...O1)'), ('NHO_deg', [p['NHO'] for p in prof], 'deg', 'angle N4-H4N...O1'),
                            ('delta_A', delta, 'Å', 'r(N4-H4N) - r(H4N...O1); > 0: proton nearer O1')):
        A(f'scan.{k}', list(arr), u, 'rev_gas', scan_src + f':items.<label>.{k if k != "delta_A" else "r_NH_A - r_HO_A"}',
          f'{what} at each profile point (order of rev.proton.scan.labels); B3LYP optimised geometries')
    A('scan.r_HO_range_A', [float(rHO.min()), float(rHO.max())], 'Å', 'rev_gas', scan_src, 'range of the scan')
    dev = max(abs(p['dev']) for p in allp if p['dev'] is not None)
    A('scan.max_constraint_dev_A', dev, 'Å', 'rev_gas', scan_src + ':items.<label>.constraint.deviation_A',
      'largest |r(H4N...O1) - target| over the 20 constrained optimisations')
    A('scan.max_opt_steps', max(p['steps'] for p in allp), '', 'rev_gas', scan_src + ':items.<label>.opt_steps',
      'largest number of optimisation steps of one scan point (limit 100)')
    hyst = {t: (max(q['E'] for q in v) - min(q['E'] for q in v)) * KC for t, v in by_t.items()
            if t != 'min' and len(v) > 1}
    A('scan.hysteresis_max_kcal', max(hyst.values()), 'kcal mol-1', 'rev_gas', scan_src + ':items',
      f"largest energy difference between the down chain (from the salt) and the neutral chain (from the neutral "
      f"start) at the {len(hyst)} shared targets {sorted(hyst)} A: no hysteresis, the two branches meet")
    emin = [q['E'] for q in by_t['min']]
    A('scan.min_vs_published_Ha', abs(emin[0] - pub('gas.e_total')), 'Ha', 'rev_gas',
      S('scan_gas_up', 'items.min.E_Ha') + ' vs gas.e_total',
      'reproduction check: tight single point at the published gas-phase minimum vs the published energy')
    k = next(i for i in range(n - 1) if delta[i] <= 0 < delta[i + 1])
    w = -delta[k] / (delta[k + 1] - delta[k])
    A('scan.transfer_r_HO_A', rHO[k] + w * (rHO[k + 1] - rHO[k]), 'Å', 'rev_gas', 'rev.proton.scan.r_HO_A, delta_A',
      f'r(H4N...O1) where r(N4-H4N) = r(H4N...O1) on the relaxed path (linear interpolation between {lab[k]} and '
      f'{lab[k + 1]})')
    A('scan.transfer_r_NO_A', rNO[k] + w * (rNO[k + 1] - rNO[k]), 'Å', 'rev_gas', 'rev.proton.scan.r_NO_A',
      'r(N4...O1) at that point (same interpolation)')
    A('scan.lumo_on_picrate_pct', [p['lumo_anion'] for p in prof], '%', 'rev_gas',
      scan_src + ':items.<label>.frontier_composition_mulliken_pct.LUMO.anion',
      'B3LYP LUMO, Mulliken share on the picrate atoms (H4N counted with the cation fragment), per profile point: '
      'on the cation down to r(H4N...O1) = 1.28 A, on the picrate (picric acid) from 1.23 A')

    # ---------------------------------------------------------------- 2. per-method profiles
    E, Esrc = sp_energies(prof)
    i_min = lab.index('min')
    i_low = {m: int(np.argmin(E[m])) for m in E if m != 'hf'}
    assert len(set(i_low.values())) == 1, i_low
    il = i_low['b3lyp']
    A('gas.lowest_point_r_HO_A', rHO[il], 'Å', 'rev_gas', 'rev.proton.scan.r_HO_A',
      f'the lowest-energy profile point ({lab[il]}, H4N on O1) is the same at every level, B3LYP, wB97X-D3(BJ), '
      'M06-2X, PBE0, CAM-B3LYP, DF-MP2 and SCS-MP2 (HF alone has its lowest point at 0.94 A)')
    A('gas.ionic_point_r_HO_A', rHO[i_min], 'Å', 'rev_gas', 'rev.proton.scan.r_HO_A',
      'the published B3LYP gas-phase ionic minimum (profile point min), the reference of every dE_kcal list')
    ionic_well = {}
    for m in [mm for mm, _, _ in METHODS] + ['hf']:
        s = SLUG.get(m, m)
        dE = (np.array(E[m]) - E[m][i_min]) * KC
        nm = NAME.get(m, 'HF')
        esrc = (scan_src + ':items.<label>.E_Ha' if m == 'b3lyp' else
                f"{S(SP_SCAN.get(m, 'sp_scan_dfmp2'))}, {S(SP_OSIDE)}: items.<label>|{m if m in SP_SCAN else 'dfmp2'}"
                f".{ {'scs-mp2': 'E_SCS_MP2_Ha', 'hf': 'E_HF_Ha'}.get(m, 'E_Ha') }")
        A(f'gas.{s}.dE_kcal', list(dE), 'kcal mol-1', 'rev_gas', esrc,
          f'{nm}: E(point) - E(published B3LYP ionic minimum geometry), per profile point (order of '
          f'rev.proton.scan.labels); ' + ('B3LYP tight single points of the scan job itself' if m == 'b3lyp' else GEOM_NOTE)
          + ('; HF = the RHF reference of the DF-MP2 job, not a level of record' if m == 'hf' else ''))
        if m == 'hf':
            continue
        de_ion = -dE[il]
        A(f'gas.{s}.dE_ionic_minus_neutral_kcal', de_ion, 'kcal mol-1', 'rev_gas', f'rev.proton.gas.{s}.dE_kcal',
          f'{nm}: E(salt, published B3LYP ionic minimum geometry, r(H4N...O1) {rHO[i_min]:.3f} A) - E(neutral, lowest '
          f'profile point, r(H4N...O1) {rHO[il]:.2f} A); positive = the neutral form lies lower. Electronic energies '
          'along the B3LYP path, no zero-point or thermal terms')
        # interior minima of the data on the N side (proton on N4: delta < 0)
        Vk = dE
        dmin = [i for i in range(1, n - 1) if Vk[i] < Vk[i - 1] and Vk[i] < Vk[i + 1] and delta[i] < 0]
        ionic_well[m] = bool(dmin)
        A(f'gas.{s}.ionic_well', 'yes' if dmin else 'no', '', 'rev_gas', f'rev.proton.gas.{s}.dE_kcal',
          f'{nm}: does the profile have a local minimum with the proton on N4 (a point lower than both neighbours, '
          f'r(N4-H4N) < r(H4N...O1))? Judged on the computed points, not on the interpolating spline (the spline of '
          f'dvr.py shows spurious 0.0003 kcal/mol dips for M06-2X and PBE0)')
        if dmin:
            i0 = dmin[0]
            j = i0 + int(np.argmax(Vk[i0:il + 1]))
            A(f'gas.{s}.ionic_min_r_HO_A', rHO[i0], 'Å', 'rev_gas', 'rev.proton.scan.r_HO_A',
              f'{nm}: r(H4N...O1) of the lowest N-side point ({lab[i0]})')
            A(f'gas.{s}.barrier_r_HO_A', rHO[j], 'Å', 'rev_gas', 'rev.proton.scan.r_HO_A',
              f'{nm}: r(H4N...O1) of the highest point between that minimum and the O side ({lab[j]})')
            A(f'gas.{s}.well_depth_kcal', Vk[j] - Vk[i0], 'kcal mol-1', 'rev_gas', f'rev.proton.gas.{s}.dE_kcal',
              f'{nm}: depth of the ionic well = E({lab[j]}) - E({lab[i0]}), barrier to proton transfer from N4 to O1 '
              f'on the computed points')
    # M06-2X grid check
    gc1, gc2 = R('sp_scan_m062x')['m062x_grid_check'], R(SP_OSIDE)['m062x_grid_check']
    A('gas.m062x.grid5_abs_max_kcal', max(abs(v) for v in list(gc1['dE_abs_Ha'].values()) + list(gc2['dE_abs_Ha'].values())) * KC,
      'kcal mol-1', 'rev_gas', S('sp_scan_m062x', 'm062x_grid_check') + '; ' + S(SP_OSIDE, 'm062x_grid_check'),
      'M06-2X: largest change of the absolute energy from grid level 4 to 5 at the three checked points '
      f"({', '.join(list(gc1['dE_abs_Ha']) + list(gc2['dE_abs_Ha']))})")
    A('gas.m062x.grid5_rel_kcal', abs(gc2['relative_energy_error_kcal']), 'kcal mol-1', 'rev_gas',
      S(SP_OSIDE, 'm062x_grid_check.relative_energy_error_kcal'),
      'M06-2X: change of the relative energy of the two O-side points from grid level 4 to 5')

    # ---------------------------------------------------------------- 3. DVR
    dvr = _done('dvr_gas_all')
    dsrc = S('dvr_gas_all')
    st = dvr['settings']
    A('dvr.settings', f"sinc-DVR (Colbert-Miller) on the {n}-point profile of each level; coordinate xH = arc "
      f"length of the path of H4N after superposing consecutive scan geometries on all other atoms, zero at the "
      f"lowest point; mass of H ({st['masses_u']['H']} u); cubic not-a-knot spline, continued {st['pad_A']:.2f} A "
      f"beyond each end with the end slope and curvature, then a hard wall; {st['npts']} grid points; "
      f"{st['nlevels']} levels", '', 'rev_gas', dsrc + ':settings', 'dft/revision/dvr.py')
    A('dvr.npts', st['npts'], '', 'rev_gas', dsrc + ':settings.npts', 'DVR grid points (checked against 2N-1)')
    A('dvr.pad_A', st['pad_A'], 'Å', 'rev_gas', dsrc + ':settings.pad_A', 'extension of the potential beyond the scan')
    it0 = dvr['items']['b3lyp']
    A('dvr.dx_A', it0['dx_A'], 'Å', 'rev_gas', dsrc + ':items.b3lyp.dx_A', 'DVR grid spacing')
    A('dvr.x_range_A', it0['scan_range_A'], 'Å', 'rev_gas', dsrc + ':items.b3lyp.scan_range_A',
      'xH of the first (r(H4N...O1) 1.70 A) and last (0.90 A) scan point')
    A('scan.x_A', [p['x_A'] for p in it0['points']], 'Å', 'rev_gas', dsrc + ':items.b3lyp.points[].x_A',
      'DVR coordinate xH of each profile point (order of rev.proton.scan.labels); the same for every level')
    A('dvr.m_eff_path_u', it0['m_eff_path_u'], '', 'rev_gas', dsrc + ':items.b3lyp.m_eff_path_u',
      'effective mass (u) of the whole relaxed path = (mass-weighted path length / proton path length)^2: the heavy '
      'atoms move a lot along the path, so the 1-D model with the mass of H is the light limit and the mwpath '
      'variant (unit mass on the mass-weighted path) the heavy limit')
    sel = dvr['selftest']
    A('dvr.selftest', f"{sum(v['pass'] for v in sel.values())} of {len(sel)} passed ({', '.join(sel)})", '',
      'rev_gas', dsrc + ':selftest', 'analytic checks of the DVR code run in the same job: harmonic and Morse levels, '
      'a spline-sampled harmonic potential, symmetric and asymmetric double wells')
    assert all(v['pass'] for v in sel.values())
    first = _done('dvr_gas_b3lyp')['items']['b3lyp']
    A('dvr.first_run_p_outside_pct', 100 * first['levels'][0]['P_outside_scan'], '%', 'rev_gas',
      S('dvr_gas_b3lyp', 'items.b3lyp.levels[0].P_outside_scan'),
      'superseded first DVR run (14 points, 0.98-1.70 A): its lowest point was the O end, so 24 % of the ground state '
      'lay outside the scan; this led to scan_gas_oside (0.90, 0.94 A). Not to be quoted as a result')
    rec = dvr_recompute(dvr)
    for m, s, nm in METHODS:
        it = dvr['items'][m]
        p = f'items.{m}'
        assert it['n_points'] == n and [q['label'] for q in it['points']] == lab, m
        assert np.allclose([q['dE_kcal'] for q in it['points']], (np.array(E[m]) - min(E[m])) * KC, atol=6e-5), m
        lv = it['levels']
        A(f'dvr.{s}.zpe_cm', it['ZPE_cm'], 'cm-1', 'rev_gas', dsrc + f':{p}.ZPE_cm',
          f'{nm}: ground proton level v = 0 above the minimum of the interpolated potential')
        A(f'dvr.{s}.zpe_kcal', it['ZPE_kcal'], 'kcal mol-1', 'rev_gas', dsrc + f':{p}.ZPE_kcal', f'{nm}: the same in kcal/mol')
        A(f'dvr.{s}.fundamental_cm', it['fundamental_cm'], 'cm-1', 'rev_gas', dsrc + f':{p}.fundamental_cm',
          f'{nm}: first excitation E(v=1) - E(v=0) along the proton path')
        A(f'dvr.{s}.v1_kcal', lv[1]['E_kcal'], 'kcal mol-1', 'rev_gas', dsrc + f':{p}.levels[1].E_kcal',
          f'{nm}: level v = 1 above the potential minimum')
        rc_ = ' (value from the recomputation in rev/proton.py, unrounded; equal to the stored 5-decimal value within 6e-6)'
        A(f'dvr.{s}.p_O_v0_pct', 100 * rec[m]['P_O'][0], '%', 'rev_gas', dsrc + f':{p}.levels[0].P_O_side',
          f'{nm}: probability of the ground state on the O side (r(N4-H4N) > r(H4N...O1), the neutral form)' + rc_)
        A(f'dvr.{s}.p_O_v1_pct', 100 * rec[m]['P_O'][1], '%', 'rev_gas', dsrc + f':{p}.levels[1].P_O_side',
          f'{nm}: the same for v = 1' + rc_)
        A(f'dvr.{s}.p_outside_v0_pct', 100 * rec[m]['P_out'][0], '%', 'rev_gas', dsrc + f':{p}.levels[0].P_outside_scan',
          f'{nm}: ground-state weight outside the scanned range (in the extrapolated pads); dvr.py flags > 1 %' + rc_)
        A(f'dvr.{s}.grid_conv_cm', it['grid_convergence_cm'], 'cm-1', 'rev_gas', dsrc + f':{p}.grid_convergence_cm',
          f"{nm}: largest level change from {st['npts']} to {2 * st['npts'] - 1} grid points (dvr.py warns above "
          '1 cm-1; the change is below 0.4 % of the ZPE)')
        A(f'dvr.{s}.spline_undershoot_kcal', it['spline_undershoot_kcal'], 'kcal mol-1', 'rev_gas',
          dsrc + f':{p}.spline_undershoot_kcal', f'{nm}: depth of the spline minimum below the lowest computed point '
          '(the DVR zero); dvr.py flags > 0.2')
        A(f'dvr.{s}.transfer_above_v0_kcal', it['V_transfer_point_kcal'] - it['ZPE_kcal'], 'kcal mol-1', 'rev_gas',
          dsrc + f':{p}.V_transfer_point_kcal - ZPE_kcal',
          f'{nm}: potential at the point r(N4-H4N) = r(H4N...O1) minus the ground level; positive = the ground level '
          'lies below the half-way point of the transfer')
        omin = [q for q in it['minima'] if q['side'] == 'O']
        assert len(omin) == 1
        A(f'dvr.{s}.neutral_min_r_HO_A', omin[0]['r_HO_A'], 'Å', 'rev_gas', dsrc + f':{p}.minima[side O].r_HO_A',
          f'{nm}: r(H4N...O1) of the minimum of the interpolated potential (O side), interpolated from xH')
        well = ionic_well[m]
        assert well == any(q['r_HO_A'] > rHO[il] + 0.1 for q in it['data_minima']), m
        if well:                              # a real ionic minimum (in the computed points, not only on the spline)
            sm, bar = it['second_minimum'], it['barrier']
            assert sm['side'] == 'N' and bar is not None, m
            A(f'dvr.{s}.ionic_min_above_v0_kcal', sm['V_kcal'] - it['ZPE_kcal'], 'kcal mol-1', 'rev_gas',
              dsrc + f':{p}.second_minimum.V_kcal - ZPE_kcal',
              f"{nm}: bottom of the ionic well (spline, r(H4N...O1) {sm['r_HO_A']:.3f} A) above the ground level; "
              'positive = the ionic minimum lies above the zero-point level, so it holds no bound proton state of its own')
            A(f'dvr.{s}.barrier_above_v0_kcal', bar['V_kcal'] - it['ZPE_kcal'], 'kcal mol-1', 'rev_gas',
              dsrc + f':{p}.barrier.V_kcal - ZPE_kcal', f"{nm}: barrier top (spline, r(H4N...O1) {bar['r_HO_A']:.3f} A) "
              'above the ground level')
            A(f'dvr.{s}.v1_minus_barrier_kcal', lv[1]['E_kcal'] - bar['V_kcal'], 'kcal mol-1', 'rev_gas',
              dsrc + f':{p}.levels[1].E_kcal - barrier.V_kcal', f'{nm}: level v = 1 relative to the barrier top')
            A(f'gas.{s}.well_depth_spline_kcal', bar['V_kcal'] - sm['V_kcal'], 'kcal mol-1', 'rev_gas',
              dsrc + f':{p}.barrier.V_kcal - second_minimum.V_kcal',
              f'{nm}: ionic-well depth on the interpolating spline (compare rev.proton.gas.{s}.well_depth_kcal)')
            A(f'gas.{s}.ionic_min_r_HO_spline_A', sm['r_HO_A'], 'Å', 'rev_gas', dsrc + f':{p}.second_minimum.r_HO_A',
              f'{nm}: position of the ionic minimum on the spline')
            A(f'gas.{s}.barrier_r_HO_spline_A', bar['r_HO_A'], 'Å', 'rev_gas', dsrc + f':{p}.barrier.r_HO_A',
              f'{nm}: position of the barrier top on the spline')
        A(f'dvr.{s}.verdict', (f"double well: ionic minimum {it['second_minimum']['V_kcal'] - it['ZPE_kcal']:.2f} "
                               f"kcal/mol above the ground level" if well else
                               'single well: no ionic minimum; the ground state sits in the neutral (O1-H) well')
          + f"; ground state {100 * rec[m]['P_O'][0]:.1f} % on the O side", '', 'rev_gas', dsrc + f':{p}',
          f'{nm}: summary; ionic minimum judged on the computed points (dvr.py verdict for M06-2X and PBE0 reports a '
          'spline dip of 0.0003 kcal/mol as a second minimum: a spline artefact, not used)')
        # sensitivity
        sz, sp_ = {'main': it['ZPE_cm']}, {'main': 100 * rec[m]['P_O'][0]}
        for v in it['sensitivity']:
            name = DVR_VARIANTS[(v['coord'], v['mass'], v['interp'], v['walls'])]
            assert not v.get('skipped'), (m, name)
            sz[name] = v['ZPE_cm']; sp_[name] = 100 * v['P_O_side_v0']
        A(f'dvr.{s}.sens_zpe_cm', sz, 'cm-1', 'rev_gas', dsrc + f':{p}.sensitivity[].ZPE_cm',
          f'{nm}: ZPE per DVR variant: main (xH, mass of H, cubic, extrapolated walls), pchip interpolation, box walls '
          'at the scan ends, reduced mass H against N+O, coordinate [r(N-H) - r(H...O)]/2, coordinate -r(H...O), '
          'mass-weighted path with unit mass (heavy-atom limit)')
        A(f'dvr.{s}.sens_p_O_v0_pct', sp_, '%', 'rev_gas', dsrc + f':{p}.sensitivity[].P_O_side_v0',
          f'{nm}: ground-state O-side probability per DVR variant (names as sens_zpe_cm)')
        lo = min(sp_, key=sp_.get)
        A(f'dvr.{s}.sens_p_O_v0_min_pct', sp_[lo], '%', 'rev_gas', f'rev.proton.dvr.{s}.sens_p_O_v0_pct',
          f'{nm}: lowest O-side probability over the variants (variant {lo}); the -r(H...O) coordinate with the mass '
          'of H compresses the path to half its length and so overstates the zero-point spread')
        ex = {k2: v2 for k2, v2 in sp_.items() if k2 != 'rHO'}
        A(f'dvr.{s}.sens_p_O_v0_min_excl_rHO_pct', min(ex.values()), '%', 'rev_gas', f'rev.proton.dvr.{s}.sens_p_O_v0_pct',
          f'{nm}: lowest O-side probability over the variants other than rHO (variant {min(ex, key=ex.get)})')
        # recomputed levels must equal the stored ones (checked in dvr_recompute); undershoot as stored
        assert abs(rec[m]['undershoot_kcal'] - it['spline_undershoot_kcal']) < 6e-5, m

    # ---------------------------------------------------------------- 4. B3LYP optimisations from the neutral form
    for model, case in (('gas', 'rev_gas'), ('pcm', 'rev_pcm')):
        job = f'neutral_start_{model}'
        d = _done(job); it = d['items'][model]
        assert it['converged'] and it['proton_on'] == 'O1'
        ref = 'gas.e_total' if model == 'gas' else 'pcm.e_total'
        assert abs(it['E_stored_ionic_min_Ha'] - pub(ref)) < 1e-8, (model, it['E_stored_ionic_min_Ha'], pub(ref))
        sy, X = read_xyz(job_dir(job) / it['xyz'])
        _, Xs = read_xyz(ROOT / f'molecules/btpc/results/btpc_{model}/opt.xyz')
        assert abs(interplanar(Xs, L) - pub(f'{model}.interplanar_picrate_cation')) < 1e-6
        assert abs(dist(X, L, 'H4N', 'O1') - it['r_HO_A']) < 1e-6
        js = S(job, f'items.{model}')
        mname = 'Gas phase' if model == 'gas' else 'IEF-PCM (water)'
        A(f'neutral_start_{model}.outcome', f"proton stays on O1: neutral minimum (O1-H4N...N4), converged in "
          f"{it['opt_steps']} steps", '', case, js + '.outcome',
          f"{mname}, B3LYP/6-311++G(d,p) (as published), no constraint; start = the published {model} minimum with H4N "
          f"moved onto O1 (O1-H {d['settings']['neutral_r_OH_A']:.2f} A along O1->N4); {it['converged_on']}. No Hessian "
          'was computed at this structure')
        A(f'neutral_start_{model}.dE_vs_salt_kcal', it['dE_vs_stored_min_kcal'], 'kcal mol-1', case,
          js + '.dE_vs_stored_min_kcal', f'{mname}, B3LYP: E(neutral minimum) - E(published ionic minimum, {ref}); '
          'negative = the neutral form is lower')
        for k2, fld, u, wh in (('r_OH_A', 'r_HO_A', 'Å', 'O1-H4N (covalent)'),
                               ('r_N_H_A', 'r_NH_A', 'Å', 'H4N...N4 (hydrogen bond)'),
                               ('r_NO_A', 'r_NO_A', 'Å', 'N4...O1'),
                               ('angle_deg', 'NHO_deg', 'deg', 'angle O1-H4N...N4')):
            A(f'neutral_start_{model}.{k2}', it[fld], u, case, js + f'.{fld}', f'{mname}, B3LYP neutral minimum: {wh}')
        A(f'neutral_start_{model}.opt_steps', it['opt_steps'], '', case, js + '.opt_steps', f'{mname}: optimisation steps')
        A(f'neutral_start_{model}.heavy_rmsd_A', it['heavy_atom_rmsd_to_stored_min_A'], 'Å', case,
          js + '.heavy_atom_rmsd_to_stored_min_A', f'{mname}: heavy-atom RMSD (Kabsch) to the published ionic minimum')
        A(f'neutral_start_{model}.interplanar_deg', interplanar(X, L), 'deg', case, rel(job_dir(job) / it['xyz']),
          f'{mname}, neutral minimum: angle between the picrate C6 ring and the 9 cation ring atoms (definition of '
          f'{model}.interplanar_picrate_cation, {pub(model + ".interplanar_picrate_cation"):.1f} deg at the salt)')
        A(f'neutral_start_{model}.r_C1O1_A', dist(X, L, 'C1', 'O1'), 'Å', case, rel(job_dir(job) / it['xyz']),
          f'{mname}, neutral minimum: C1-O1 (phenol C-O; {pub(model + ".r_O1_C1"):.3f} A in the salt)')
        A(f'neutral_start_{model}.lumo_on_picrate_pct', it['frontier_composition_mulliken_pct']['LUMO']['anion'], '%',
          case, js + '.frontier_composition_mulliken_pct.LUMO.anion',
          f'{mname}, neutral minimum: Mulliken share of the LUMO on the picrate atoms (H4N counted with the cation)')
    ng = R('neutral_start_gas')['items']['gas']
    p98 = next(p for p in prof if p['t'] == 0.98)
    A('neutral_start_gas.dE_vs_scan_kcal', (ng['E_Ha'] - p98['E']) * KC, 'kcal mol-1', 'rev_gas',
      S('neutral_start_gas', 'items.gas.E_Ha') + ' - ' + S(p98['job'], f"items.{p98['label']}.E_Ha"),
      'free neutral minimum minus the constrained scan point r(H4N...O1) = 0.98 A: the scan point lies just off the '
      'minimum (O1-H 0.998 A when free)')

    # ---------------------------------------------------------------- 5. IEF-PCM: E(neutral) - E(salt) by method
    smin, sneu, sgas = _done('sp_min_pcm'), _done('sp_neutral_pcm'), _done('sp_gas_pcmgeoms')
    for dd_, sv in ((smin, 'water'), (sneu, 'water'), (sgas, None)):
        assert dd_['settings']['basis'] == '6-311++g(d,p)' and dd_['settings']['solvent'] == sv
    for dd_, g_lab, ref in ((smin, 'opt', ROOT / 'molecules/btpc/results/btpc_pcm/opt.xyz'),
                            (sneu, 'neutral_pcm', job_dir('neutral_start_pcm') / 'points/pcm.xyz'),
                            (sgas, 'opt', ROOT / 'molecules/btpc/results/btpc_pcm/opt.xyz'),
                            (sgas, 'neutral_pcm', job_dir('neutral_start_pcm') / 'points/pcm.xyz')):
        g = {q['label']: q for q in dd_['settings']['geometries']}[g_lab]
        jd = job_dir(dd_['job'])
        assert np.abs(read_xyz(jd / g['local'])[1] - read_xyz(ref)[1]).max() < 1e-7, (dd_['job'], g_lab)
    for run in CC_RUNS:
        for jb, _ in run[4:6]:
            inp = R(jb)['inputs']
            for kk, v in inp.items():
                loc = ROOT / v['path'] if not v['path'].startswith(('/', '<run-dir>')) else job_dir('sp_neutral_pcm') / 'geoms/neutral_pcm.xyz'
                assert sha256(loc) == v['sha256'], (jb, kk)
    nsp = R('neutral_start_pcm')['items']['pcm']
    pcm = {'b3lyp': (nsp['E_Ha'] - nsp['E_stored_ionic_min_Ha'],
                     S('neutral_start_pcm', 'items.pcm.E_Ha') + ' - pcm.e_total', nsp['E_stored_ionic_min_Ha'], nsp['E_Ha'])}
    for m in ('wb97x-d3bj', 'm062x', 'pbe0', 'cam-b3lyp', 'dfmp2'):
        a, b = smin['items'][f'opt|{m}'], sneu['items'][f'neutral_pcm|{m}']
        assert a['status'] == b['status'] == 'done'
        pcm[m] = (b['E_Ha'] - a['E_Ha'], S('sp_neutral_pcm', f'items.neutral_pcm|{m}.E_Ha') + ' - ' +
                  S('sp_min_pcm', f'items.opt|{m}.E_Ha'), a['E_Ha'], b['E_Ha'])
    a, b = smin['items']['opt|dfmp2'], sneu['items']['neutral_pcm|dfmp2']
    for m, fld in (('scs-mp2', 'E_SCS_MP2_Ha'), ('hf', 'E_HF_Ha')):
        pcm[m] = (b[fld] - a[fld], S('sp_neutral_pcm', f'items.neutral_pcm|dfmp2.{fld}') + ' - ' +
                  S('sp_min_pcm', f'items.opt|dfmp2.{fld}'), a[fld], b[fld])
    pgeo = ('two B3LYP/6-311++G(d,p) IEF-PCM (water) geometries: the published salt minimum (molecules/btpc/results/'
            'btpc_pcm/opt.xyz) and the neutral minimum of neutral_start_pcm; IEF-PCM eps 78.3553, Lebedev 29, '
            '6-311++G(d,p)')
    for m, (dE, src_, _, _) in pcm.items():
        s = SLUG.get(m, m)
        extra = {'dfmp2': '; DF-MP2 on an IEF-PCM HF reference with the reaction field frozen at the HF density (PTE), '
                          'frozen core, RI basis aug-cc-pVTZ-RI',
                 'scs-mp2': '; SCS-MP2 from the same DF-MP2 job (6/5 OS + 1/3 SS)',
                 'hf': '; HF = the IEF-PCM RHF reference of the DF-MP2 job',
                 'b3lyp': '; B3LYP: the optimised energies themselves'}.get(m, '; single points')
        A(f'pcm.{s}.dE_neutral_minus_salt_kcal', dE * KC, 'kcal mol-1', 'rev_pcm', src_,
          f"IEF-PCM (water), {NAME.get(m, 'HF')}: E(neutral pair) - E(salt) at the {pgeo}; positive = the salt is lower"
          + extra)
    A('pcm.dfmp2.dE_corr_part_kcal', (pcm['dfmp2'][0] - pcm['hf'][0]) * KC, 'kcal mol-1', 'rev_pcm',
      'rev.proton.pcm.dfmp2 - rev.proton.pcm.hf', 'IEF-PCM (water): MP2 correlation part of E(neutral) - E(salt) '
      '(DF-MP2 minus its HF reference)')
    A('pcm.n_dft_salt_lower', sum(1 for m in ('b3lyp', 'wb97x-d3bj', 'm062x', 'pbe0', 'cam-b3lyp') if pcm[m][0] > 0),
      '', 'rev_pcm', 'rev.proton.pcm.*.dE_neutral_minus_salt_kcal', 'number of the five functionals (B3LYP, '
      'wB97X-D3(BJ), M06-2X, PBE0, CAM-B3LYP) that put the salt below the neutral pair in IEF-PCM (water)')
    gm = {}
    for m, fld in (('dfmp2', 'E_Ha'), ('scs-mp2', 'E_SCS_MP2_Ha'), ('hf', 'E_HF_Ha')):
        a, b = sgas['items']['opt|dfmp2'], sgas['items']['neutral_pcm|dfmp2']
        gm[m] = (b[fld] - a[fld]) * KC
        A(f'gas_pcmgeom.{SLUG.get(m, m)}.dE_neutral_minus_salt_kcal', gm[m], 'kcal mol-1', 'rev_gas',
          S('sp_gas_pcmgeoms', f'items.neutral_pcm|dfmp2.{fld}') + ' - ' + S('sp_gas_pcmgeoms', f'items.opt|dfmp2.{fld}'),
          f"Gas phase {NAME.get(m, 'HF')}/6-311++G(d,p) at the two IEF-PCM geometries (sp_gas_pcmgeoms): separates "
          'the solvent reaction field from the intrinsic (gas-phase) preference at fixed geometry')
    A('pcm.dfmp2.solvent_shift_kcal', pcm['dfmp2'][0] * KC - gm['dfmp2'], 'kcal mol-1', 'rev_pcm',
      'rev.proton.pcm.dfmp2.dE_neutral_minus_salt_kcal - rev.proton.gas_pcmgeom.dfmp2.dE_neutral_minus_salt_kcal',
      'DF-MP2: how far the IEF-PCM reaction field moves E(neutral) - E(salt) toward the salt at fixed geometries')

    # ---------------------------------------------------------------- 6. CCSD(T) corrections
    dfmp2_pcm = pcm['dfmp2'][0] * KC
    cc_rows, run_dd = [], {}
    diag = {'T1': [], 'D1': [], 'D2': []}
    for run, model, basis, frac, js_, jn_ in CC_RUNS:
        case = f'rev_{model}'
        recs = {}
        for sp_name, (jb, key) in (('salt', js_), ('neutral', jn_)):
            d = _done(jb); it = d['items'][key]
            st_ = d['settings']
            assert it['status'] == 'done' and it['frac'] == frac and it['geometry'] == sp_name
            assert st_['basis'] == {'6-31G*': '6-31g*', '6-31+G(d,p)': '6-31+g**'}[basis]
            assert (st_['solvent'] == 'water') == (model == 'pcm')
            assert 'dfccsd' in it['cc_class'] and max(it['fock_offdiag_max'].values()) < 1e-6
            recs[sp_name] = it
            src_ = S(jb, f'items.{key}')
            pre = f'cc.{run}.{sp_name}'
            what = (f"{'IEF-PCM (water), PTE reference' if model == 'pcm' else 'Gas phase'}, {basis}, FNO "
                    f'{100 * frac:.0f} % of the virtual space, {sp_name}')
            A(f'{pre}.nvir_kept', it['nvir_kept'], '', case, src_ + '.nvir_kept', f'{what}: virtual natural orbitals kept')
            A(f'{pre}.nvir_total', it['nvir_total'], '', case, src_ + '.nvir_total',
              f'{what}: virtual orbitals in the basis' + ('; 400 for the neutral, 401 for the salt: PySCF removed one '
                                                         'near-linearly-dependent AO combination at the neutral '
                                                         'geometry (diffuse functions); both keep 200 NOs'
                                                         if basis == '6-31+G(d,p)' else ''))
            A(f'{pre}.no_occ_pct', 100 * it['no_occ_kept_frac'], '%', case, src_ + '.no_occ_kept_frac',
              f'{what}: share of the virtual NO occupation sum (MP2 density) kept')
            for fld, k2, wh in (('delta_ccsd_Ha', 'delta_ccsd_kcal', 'E_corr(CCSD) - E_corr(MP2), same FNO space'),
                                ('delta_T_Ha', 'delta_T_kcal', '(T) energy'),
                                ('delta_Ha', 'delta_kcal', 'delta = E_corr(CCSD(T)) - E_corr(MP2) in the FNO space'),
                                ('delta_T_scaled_Ha', 'delta_Tscaled_kcal', 'delta with (T) scaled by MP2(full)/MP2(FNO)')):
                A(f'{pre}.{k2}', it[fld] * KC, 'kcal mol-1', case, src_ + f'.{fld}', f'{what}: {wh}')
            for dk in ('T1', 'D1', 'D2'):
                A(f'{pre}.{dk}', it['diagnostics'][dk], '', case, src_ + f'.diagnostics.{dk}', f'{what}: CCSD {dk} diagnostic')
                diag[dk].append(it['diagnostics'][dk])
            cc_rows.append([run, model, basis, frac, sp_name, jb, key, it['nvir_kept'], it['nvir_total'],
                            it['no_occ_kept_frac'], it['E_HF_Ha'], it['E_corr_mp2_full_ri_Ha'], it['E_corr_mp2_active_Ha'],
                            it['E_corr_ccsd_Ha'], it['E_T_Ha'], it['delta_ccsd_Ha'], it['delta_T_Ha'], it['delta_Ha'],
                            it['delta_T_scaled_Ha'], it['diagnostics']['T1'], it['diagnostics']['D1'],
                            it['diagnostics']['D2'], it['wall_s']])
        s_, n_ = recs['salt'], recs['neutral']
        assert s_['nvir_kept'] == n_['nvir_kept']
        what = (f"{'IEF-PCM (water)' if model == 'pcm' else 'Gas phase'}, {basis}, FNO {100 * frac:.0f} %")
        dd = {}
        for fld, k2, wh in (('delta_ccsd_Ha', 'dd_ccsd_kcal', 'CCSD part'), ('delta_T_Ha', 'dd_T_kcal', '(T) part'),
                            ('delta_Ha', 'dd_kcal', 'CCSD(T) - MP2 correction'),
                            ('delta_T_scaled_Ha', 'dd_Tscaled_kcal', 'correction with scaled (T)')):
            dd[k2] = (n_[fld] - s_[fld]) * KC
            A(f'cc.{run}.{k2}', dd[k2], 'kcal mol-1', case,
              f'rev.proton.cc.{run}.neutral.{k2.replace("dd_", "delta_")} - salt',
              f'{what}: {wh} of E(neutral) - E(salt), delta(neutral) - delta(salt); positive = toward the salt')
        run_dd[run] = dd
        ref_mp2 = dfmp2_pcm if model == 'pcm' else gm['dfmp2']
        refk = 'rev.proton.pcm.dfmp2' if model == 'pcm' else 'rev.proton.gas_pcmgeom.dfmp2'
        for k2, ddk in (('composite_kcal', 'dd_kcal'), ('composite_Tscaled_kcal', 'dd_Tscaled_kcal')):
            A(f'cc.{run}.{k2}', ref_mp2 + dd[ddk], 'kcal mol-1', case, f'{refk}.dE_neutral_minus_salt_kcal + rev.proton.cc.{run}.{ddk}',
              f'{what}: composite E(neutral) - E(salt) = DF-MP2/6-311++G(d,p) + CCSD(T) correction'
              + (' (scaled (T))' if 'Tscaled' in k2 else '') + '; negative = the neutral pair is lower')
    # FNO series (6-31G*, IEF-PCM): linear extrapolation in x = 1 - kept occupation fraction
    ser = ['pcm_631gd_f050', 'pcm_631gd_f065', 'pcm_631gd_f080']
    xs = []
    for run in ser:
        r_ = next(c for c in CC_RUNS if c[0] == run)
        xs.append(1 - 0.5 * (R(r_[4][0])['items'][r_[4][1]]['no_occ_kept_frac'] + R(r_[5][0])['items'][r_[5][1]]['no_occ_kept_frac']))
    xs = np.array(xs)
    fr_ = [R(c[4][0])['items'][c[4][1]]['frac'] for r in ser for c in CC_RUNS if c[0] == r]
    A('cc.fno_pct', [100 * v for v in fr_], '%', 'rev_pcm',
      ', '.join(S(j) + ':items.salt|f*.frac' for j in ('cc_salt_631gd', 'cc_salt_631gd_f080')),
      'FNO sizes of the 6-31G* IEF-PCM series: kept share of the virtual orbitals (fixed NO counts, equal for salt and '
      'neutral)')
    A('cc.extrap_x', list(xs), '', 'rev_pcm', 'rev.proton.cc.pcm_631gd_f0{50,65,80}.{salt,neutral}.no_occ_pct',
      'x = 1 - kept fraction of the virtual NO occupation (mean of salt and neutral, which differ by < 2e-4), for '
      'FNO 50, 65, 80 % (6-31G*, IEF-PCM)')
    ext = {}
    for k2 in ('dd_kcal', 'dd_Tscaled_kcal'):
        y = np.array([run_dd[r][k2] for r in ser])
        slope2 = (y[1] - y[2]) / (xs[1] - xs[2])
        e2 = y[2] - slope2 * xs[2]
        c3 = np.polyfit(xs, y, 1)
        e3 = c3[1]
        res3 = float(np.abs(np.polyval(c3, xs) - y).max())
        tag = '' if k2 == 'dd_kcal' else '_Tscaled'
        A(f'cc.extrap_2pt{tag}_kcal', e2, 'kcal mol-1', 'rev_pcm', f'rev.proton.cc.pcm_631gd_f065/f080.{k2}, rev.proton.cc.extrap_x',
          f'CCSD(T) correction{" (scaled (T))" if tag else ""} extrapolated to the full virtual space (x = 0) by the '
          f'line through the 65 and 80 % points (slope {slope2:.2f} kcal/mol per unit x); 6-31G*, IEF-PCM (water)')
        A(f'cc.extrap_3pt{tag}_kcal', e3, 'kcal mol-1', 'rev_pcm', f'rev.proton.cc.pcm_631gd_f0{{50,65,80}}.{k2}, rev.proton.cc.extrap_x',
          f'the same with a least-squares line through the 50, 65 and 80 % points (slope {c3[0]:.2f}; largest residual '
          f'{res3:.2f} kcal/mol: the series is not monotonic, 50 % lies off the line)')
        ext[k2] = (e2, e3)
        if k2 == 'dd_kcal':
            for nm_, e_ in (('2pt', e2), ('3pt', e3)):
                A(f'cc.extrap_{nm_}_composite_kcal', dfmp2_pcm + e_, 'kcal mol-1', 'rev_pcm',
                  f'rev.proton.pcm.dfmp2.dE_neutral_minus_salt_kcal + rev.proton.cc.extrap_{nm_}_kcal',
                  f'IEF-PCM (water): composite E(neutral) - E(salt) with the {nm_} extrapolated correction')
    lo_, hi_ = min(run_dd['pcm_631gd_f080']['dd_kcal'], *ext['dd_kcal']), max(run_dd['pcm_631gd_f080']['dd_kcal'], *ext['dd_kcal'])
    A('cc.dd_range_kcal', [lo_, hi_], 'kcal mol-1', 'rev_pcm', 'rev.proton.cc.pcm_631gd_f080.dd_kcal, extrap_2pt_kcal, extrap_3pt_kcal',
      'range of the CCSD(T) correction to E(neutral) - E(salt) in IEF-PCM (water): the 80 % value and the two '
      'extrapolations to the full virtual space (the 80 % value lies inside)')
    best = dfmp2_pcm + run_dd['pcm_631gd_f080']['dd_kcal']
    A('cc.best_kcal', best, 'kcal mol-1', 'rev_pcm', 'rev.proton.cc.pcm_631gd_f080.composite_kcal',
      'IEF-PCM (water) CCSD(T) estimate of E(neutral) - E(salt) at the B3LYP IEF-PCM geometries: DF-MP2/6-311++G(d,p) '
      '+ FNO-CCSD(T)/6-31G* correction at 80 % of the virtual space (222 of 277 NOs); negative = the neutral pair is lower')
    A('cc.best_range_kcal', [dfmp2_pcm + lo_, dfmp2_pcm + hi_], 'kcal mol-1', 'rev_pcm',
      'rev.proton.pcm.dfmp2.dE_neutral_minus_salt_kcal + rev.proton.cc.dd_range_kcal',
      'IEF-PCM (water): composite E(neutral) - E(salt) over the range of the correction')
    A('cc.flip_needed_kcal', -dfmp2_pcm, 'kcal mol-1', 'rev_pcm', '-rev.proton.pcm.dfmp2.dE_neutral_minus_salt_kcal',
      'correction that would be needed to bring the salt level with the neutral pair at the DF-MP2 baseline')
    pr = [r[0] for r in CC_RUNS if r[1] == 'pcm']
    mx = max(pr, key=lambda r: run_dd[r]['dd_kcal'])
    A('cc.dd_max_kcal', run_dd[mx]['dd_kcal'], 'kcal mol-1', 'rev_pcm', f'rev.proton.cc.{mx}.dd_kcal',
      f'largest CCSD(T) correction among the IEF-PCM runs ({mx})')
    A('cc.basis_check_kcal', run_dd['pcm_631ppgd_f050']['dd_kcal'] - run_dd['pcm_631gd_f050']['dd_kcal'], 'kcal mol-1',
      'rev_pcm', 'rev.proton.cc.pcm_631ppgd_f050.dd_kcal - rev.proton.cc.pcm_631gd_f050.dd_kcal',
      'basis check at FNO 50 %: 6-31+G(d,p) minus 6-31G* (the 6-31+G(d,p) run keeps 97.3 % of the NO occupation, '
      'the 6-31G* run 93.0 %, so the two differ in truncation too)')
    A('cc.gas_minus_pcm_kcal', run_dd['gas_631gd_f050']['dd_kcal'] - run_dd['pcm_631gd_f050']['dd_kcal'], 'kcal mol-1',
      'rev_pcm', 'rev.proton.cc.gas_631gd_f050.dd_kcal - rev.proton.cc.pcm_631gd_f050.dd_kcal',
      'gas-phase check (6-31G*, FNO 50 %, same two IEF-PCM geometries): the correction without the continuum minus '
      'with it')
    A('cc.ccsd_change_65_80_kcal', run_dd['pcm_631gd_f080']['dd_ccsd_kcal'] - run_dd['pcm_631gd_f065']['dd_ccsd_kcal'],
      'kcal mol-1', 'rev_pcm', 'rev.proton.cc.pcm_631gd_f080.dd_ccsd_kcal - pcm_631gd_f065.dd_ccsd_kcal',
      'change of the CCSD part from 65 to 80 % (converged)')
    A('cc.T_change_65_80_kcal', run_dd['pcm_631gd_f080']['dd_T_kcal'] - run_dd['pcm_631gd_f065']['dd_T_kcal'],
      'kcal mol-1', 'rev_pcm', 'rev.proton.cc.pcm_631gd_f080.dd_T_kcal - pcm_631gd_f065.dd_T_kcal',
      'change of the (T) part from 65 to 80 % (not converged)')
    for dk, v in diag.items():
        A(f'cc.{dk}_range', [min(v), max(v)], '', 'rev_pcm', 'rev.proton.cc.*.{salt,neutral}.' + dk,
          f'range of the CCSD {dk} diagnostic over all ten CC calculations (salt and neutral alike)')
    A('cc.settings', 'focal-point correction delta = E_corr[CCSD(T)] - E_corr[MP2] in a reduced basis and a frozen-'
      'natural-orbital (FNO) virtual space (MP2 natural orbitals, semicanonical), DF-CCSD(T) and MP2 on the same '
      'RI integrals (cc-pVDZ-RI for 6-31G*, heavy-aug-cc-pVDZ-RI for 6-31+G(d,p)), frozen core; reference DF-RHF, '
      'in IEF-PCM (water) with the reaction field frozen at the HF density (PTE), as the DF-MP2 numbers; added to '
      'DF-MP2/6-311++G(d,p)', '', 'rev_pcm', S('cc_salt_631gd', 'settings'), 'dft/revision/cc_correction.py')

    # ---------------------------------------------------------------- 7. derived CSV files (build_facts.py only)
    if not writing:
        print('rev.proton: standalone run, derived CSV files not written')
        return
    DERIVED.mkdir(parents=True, exist_ok=True)
    hdr = ['label', 'job', 'chain', 'r_HO_A', 'r_NH_A', 'r_NO_A', 'NHO_deg', 'delta_A', 'x_A', 'lumo_on_picrate_pct']
    meths = [m for m, _, _ in METHODS] + ['hf']
    hdr += [f'dE_kcal_{SLUG.get(m, m)}' for m in meths] + [f'E_Ha_{SLUG.get(m, m)}' for m in meths]
    rows = []
    for i, p in enumerate(prof):
        rows.append([p['label'], p['job'], p['chain'], _g(p['r_HO']), _g(p['r_NH']), _g(p['r_NO']), _g(p['NHO']),
                     _g(float(delta[i])), _g(it0['points'][i]['x_A']), _g(p['lumo_anion'])]
                    + [_g(float((E[m][i] - E[m][i_min]) * KC)) for m in meths] + [f'{E[m][i]:.10f}' for m in meths])
    _write_csv('rev_proton_scan_profiles.csv', hdr, rows)
    rows = []
    for m, s, _ in METHODS:
        r_ = rec[m]
        for k2 in range(len(r_['E_cm'])):
            rows.append([s, k2, f"{r_['E_cm'][k2]:.4f}", f"{r_['E_kcal'][k2]:.6f}", f"{r_['P_O'][k2]:.6f}",
                         f"{r_['P_out'][k2]:.6f}"])
    _write_csv('rev_proton_dvr_levels.csv', ['method', 'n', 'E_cm', 'E_kcal', 'P_O_side', 'P_outside_scan'], rows)
    g0 = rec['b3lyp']['g_A']
    for m in rec:
        assert np.allclose(rec[m]['g_A'], g0) and np.allclose(rec[m]['x_pts_A'], rec['b3lyp']['x_pts_A'])
    rg = r_on_grid(g0, rec['b3lyp']['x_pts_A'], rec['b3lyp']['r_pts_A'])
    hdr = ['x_A', 'r_HO_A', 'in_scan']
    for m, s, _ in METHODS:
        hdr += [f'V_kcal_{s}', f'psi2_v0_{s}']
    xp = rec['b3lyp']['x_pts_A']
    rows = []
    for i, x in enumerate(g0):
        row = [f'{x:.6f}', f'{rg[i]:.6f}', int(xp[0] <= x <= xp[-1])]
        for m, s, _ in METHODS:
            row += [f"{min(rec[m]['V_kcal'][i], 999.0):.5f}", f"{rec[m]['P'][i, 0]:.5e}"]
        rows.append(row)
    _write_csv('rev_proton_dvr_grid.csv', hdr, rows)
    rows = [[SLUG.get(m, m), f'{v[2]:.10f}', f'{v[3]:.10f}', f'{v[0] * KC:.6f}', v[1]] for m, v in pcm.items()]
    _write_csv('rev_proton_pcm.csv', ['method', 'E_salt_Ha', 'E_neutral_Ha', 'dE_neutral_minus_salt_kcal', 'source'], rows)
    _write_csv('rev_proton_cc.csv', ['run', 'model', 'basis', 'fno_frac', 'pair', 'job', 'item', 'nvir_kept', 'nvir_total',
                                     'no_occ_kept_frac', 'E_HF_Ha', 'E_corr_mp2_full_Ha', 'E_corr_mp2_fno_Ha',
                                     'E_corr_ccsd_Ha', 'E_T_Ha', 'delta_ccsd_Ha', 'delta_T_Ha', 'delta_Ha',
                                     'delta_T_scaled_Ha', 'T1', 'D1', 'D2', 'wall_s'],
               [[_g(v, 12) if isinstance(v, float) else v for v in r] for r in cc_rows])
    print('rev.proton: wrote derived/rev_proton_{scan_profiles,dvr_levels,dvr_grid,pcm,cc}.csv')
