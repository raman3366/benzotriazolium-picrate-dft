"""
One-dimensional proton vibrational levels on a proton-scan profile (revision item 2): sinc-DVR.

Pure numpy/scipy (matplotlib only for the figures, with a plain-SVG fallback); no PySCF, so it runs on the Pi,
piserver or any machine with the scan outputs.  It answers the plan's question: does any second minimum along
the proton coordinate lie below the proton's zero-point level?

Profiles: B3LYP from proton_scan.py outputs (--scan DIR, repeatable; where chains reached the same target the
lowest energy is used, and the free minimum is included), and every method of singlepoints.py outputs (--sp DIR,
repeatable: wb97x-d3bj, m062x, pbe0, cam-b3lyp, dfmp2 and scs-mp2) at those geometries.

Proton coordinate (--coord) and mass (--mass); every combination below is also run as a sensitivity check:
  xH          arc length of the proton's path, after superposing consecutive geometries on all other atoms
              (default; with the mass of H, 1.00782503 u)
  half_delta  (r(N-H) - r(H...O)) / 2, the proton's displacement for a fixed N...O line (mass of H)
  rHO         -r(H...O), the scanned distance itself (mass of H)
  mwpath      mass-weighted arc length over all atoms (u^1/2 A) with unit mass: the effective mass of the relaxed
              path, i.e. the heavy atoms follow the proton adiabatically (the heavy-mass limit)
  --mass reduced: mu = m_H (m_N + m_O) / (m_H + m_N + m_O), H against a rigid N...O pair
Profile: cubic spline (not-a-knot) through the points (pchip as a check).  Outside the scanned range the potential
continues with its end slope (made to rise outward) and end curvature (>= 0), with a hard wall --pad beyond the
end; the 'box' variant puts hard walls at the scanned ends.  The weight of each level outside the scanned range
is reported; above 1 % the scan range limits the result (flagged).

DVR: Colbert & Miller, J. Chem. Phys. 96, 1982 (1992), uniform grid, T_ii = pi^2/3 * h, T_ij = 2(-1)^(i-j)/(i-j)^2 * h,
h = hbar^2 / (2 m dx^2); atomic units; converged against a grid of half the spacing (reported).
Reported per profile: levels (cm-1 and kcal/mol above the potential minimum), ZPE and fundamental, minima and
barriers of the profile, the transfer point (r(N-H) = r(H...O)), each level's weight on the O side and outside the
scan, classical turning points of the ground level, and a plain verdict.

Usage:
  python dft/revision/dvr.py --job ID --out DIR --scan DIR [--scan DIR ...] [--sp DIR ...] [--methods all|LIST]
         [--chains all|up,down,neutral] [--coord xH] [--mass H] [--interp cubic] [--walls extrapolate]
         [--npts 801] [--pad 0.30] [--nlevels 6] [--smoke]
--smoke (or --selftest): analytic checks first (harmonic and Morse levels, a spline-sampled harmonic, a symmetric
double well), then the given profiles with the same settings.  Outputs: result.json, log.txt, levels.csv,
dvr_<method>.png|svg, profiles.png|svg.
"""
import os, sys
for _k in ('THREADS', 'OMP_NUM_THREADS'):
    _v = os.environ.get(_k, '').strip()
    if _v.isdigit() and int(_v) > 0:
        THREADS = int(_v); break
else:
    THREADS = 1
for _v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = str(THREADS)
import json, time, argparse, pathlib, hashlib, socket, platform, datetime, traceback, signal   # noqa: E402
import numpy as np                                                                               # noqa: E402
import scipy                                                                                     # noqa: E402
from scipy.interpolate import CubicSpline, PchipInterpolator                                     # noqa: E402
from scipy.linalg import eigh                                                                    # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
IST = datetime.timezone(datetime.timedelta(hours=5, minutes=30), 'IST')
BOHR_A = 0.529177210903            # CODATA 2018
AMU_ME = 1822.888486209            # CODATA 2018, u in electron masses
HARTREE_CM = 219474.6313632        # CODATA 2018
HARTREE_KCAL = 627.509474
MASS = {'H': 1.00782503207, 'C': 12.0, 'N': 14.0030740048, 'O': 15.99491461956, 'Cl': 34.96885268,
        'S': 31.97207100, 'Zn': 63.9291422}
V_WALL = 1.0                        # Hartree: hard wall
COORDS = ('xH', 'half_delta', 'rHO', 'mwpath')
METHOD_ORDER = ['b3lyp', 'wb97x-d3bj', 'm062x', 'pbe0', 'cam-b3lyp', 'dfmp2', 'scs-mp2', 'b3lyp-d3bj']
METHOD_LABEL = {'b3lyp': 'B3LYP', 'wb97x-d3bj': 'ωB97X-D3(BJ)', 'm062x': 'M06-2X', 'pbe0': 'PBE0',
                'cam-b3lyp': 'CAM-B3LYP', 'dfmp2': 'DF-MP2', 'scs-mp2': 'SCS-MP2', 'b3lyp-d3bj': 'B3LYP-D3(BJ)'}
# categorical palette (dataviz reference instance, light mode, fixed order; text in ink, never series colour)
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'


def now_ist():
    return datetime.datetime.now(IST).isoformat(timespec='seconds')


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def save_json(path, obj):
    path = pathlib.Path(path); tmp = path.with_name(path.name + '.tmp')
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=2, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o))
    os.replace(tmp, path)


# ---------------------------------------------------------------- minimal job record (rev_common.Job without PySCF)
class Job:
    def __init__(self, args, script):
        self.args = args
        self.out = pathlib.Path(args.out).resolve(); self.out.mkdir(parents=True, exist_ok=True)
        self.path = self.out / 'result.json'
        self._logf = open(self.out / 'log.txt', 'a', buffering=1)
        prev = {}
        if self.path.exists():
            try:
                prev = json.load(open(self.path))
            except Exception:
                prev = {}
        if prev and bool(prev.get('smoke')) != bool(args.smoke):
            raise SystemExit(f'{self.path} holds a {"smoke" if prev.get("smoke") else "production"} run; use another --out')
        here = pathlib.Path(script).resolve()
        sw = {'python': sys.version.split()[0], 'numpy': np.__version__, 'scipy': scipy.__version__,
              'executable': sys.executable}
        try:
            import matplotlib
            sw['matplotlib'] = matplotlib.__version__
        except Exception:
            sw['matplotlib'] = None
        self.r = {'job': args.job, 'script': str(here.relative_to(ROOT)) if str(here).startswith(str(ROOT)) else str(here),
                  'argv': sys.argv, 'smoke': bool(args.smoke), 'status': 'running',
                  'started': prev.get('started') or now_ist(), 'restarts': prev.get('restarts', -1) + 1,
                  'last_start': now_ist(), 'finished': None,
                  'host': {'hostname': socket.gethostname(), 'machine': platform.machine(), 'cpu_count': os.cpu_count()},
                  'threads': THREADS, 'software': sw, 'script_sha256': {here.name: sha256_file(here)},
                  'settings': {}, 'inputs': prev.get('inputs', {}), 'items': prev.get('items', {}),
                  'timings': prev.get('timings', {}), 'error': None}
        self.save()
        self.log(f"job {args.job} ({self.r['script']}) on {self.r['host']['hostname']}: threads {THREADS}, "
                 f"smoke {bool(args.smoke)}, restart #{self.r['restarts']}; numpy {np.__version__}, scipy {scipy.__version__}")

    def log(self, msg):
        line = f"[{datetime.datetime.now(IST).strftime('%Y-%m-%d %I:%M:%S %p IST')}] {msg}"
        print(line, flush=True); self._logf.write(line + '\n')

    def save(self):
        save_json(self.path, self.r)

    def get(self, key):
        return self.r['items'].get(key)

    def put(self, key, rec, status='done'):
        rec = dict(rec); rec['status'] = status; rec['at'] = now_ist()
        self.r['items'][key] = rec; self.save()

    def input_file(self, label, path):
        p = pathlib.Path(path)
        self.r['inputs'][label] = {'path': str(p), 'sha256': sha256_file(p)}; self.save()

    def finish(self, ok, error=None):
        self.r['status'] = 'done' if ok else 'failed'; self.r['finished'] = now_ist()
        if error:
            self.r['error'] = error
        self.save()
        self.log(f"job {self.r['job']} {'finished' if ok else 'FAILED'}")


# ---------------------------------------------------------------- geometry
def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0])
    syms = [l.split()[0] for l in lines[2:2 + n]]
    return syms, np.array([[float(v) for v in l.split()[1:4]] for l in lines[2:2 + n]])


def superpose(A, B, idx, w=None):
    """B superposed onto A (Kabsch) using atoms idx, optionally weighted."""
    idx = np.asarray(idx); w = np.ones(len(A)) if w is None else np.asarray(w, float)
    wa = w[idx] / w[idx].sum()
    ca, cb = (A[idx] * wa[:, None]).sum(0), (B[idx] * wa[:, None]).sum(0)
    U, S, Vt = np.linalg.svd(((B[idx] - cb) * wa[:, None]).T @ (A[idx] - ca))
    d = np.sign(np.linalg.det(U @ Vt))
    return (B - cb) @ (U @ np.diag([1, 1, d]) @ Vt) + ca


def path_coordinates(pts, H, N, O):
    """pts sorted from the N side (large r(H...O)) to the O side.  Returns dict of coordinate arrays (zero at the
    lowest-energy point), r_HO, r_NH, delta and the overall effective mass of the relaxed path."""
    X = [p['xyz'] for p in pts]; syms = pts[0]['syms']
    m = np.array([MASS.get(s, MASS.get(s[:1], 12.0)) for s in syms])
    rHO = np.array([np.linalg.norm(x[H] - x[O]) for x in X]); rNH = np.array([np.linalg.norm(x[H] - x[N]) for x in X])
    i0 = int(np.argmin([p['E'] for p in pts]))
    others = [i for i in range(len(syms)) if i != H]
    allat = list(range(len(syms)))
    sH, sM = [0.0], [0.0]
    for k in range(len(X) - 1):
        B = superpose(X[k], X[k + 1], others)
        sH.append(sH[-1] + float(np.linalg.norm(B[H] - X[k][H])))
        Bm = superpose(X[k], X[k + 1], allat, m)
        sM.append(sM[-1] + float(np.sqrt((m * ((Bm - X[k]) ** 2).sum(1)).sum())))
    sH, sM = np.array(sH), np.array(sM)
    delta = rNH - rHO
    c = {'xH': sH - sH[i0], 'half_delta': 0.5 * (delta - delta[i0]), 'rHO': -(rHO - rHO[i0]), 'mwpath': sM - sM[i0]}
    meff = float((sM[-1] - sM[0]) ** 2 / (sH[-1] - sH[0]) ** 2) if sH[-1] > sH[0] else None
    seg = np.diff(sM) ** 2 / np.maximum(np.diff(sH) ** 2, 1e-30)
    return c, {'r_HO_A': rHO, 'r_NH_A': rNH, 'delta_A': delta, 'i_min': i0,
               'm_eff_path_u': meff, 'm_eff_segments_u': seg}


# ---------------------------------------------------------------- the DVR
def kinetic(n, dx, mass_me):
    i = np.arange(n)
    d = i[:, None] - i[None, :]
    with np.errstate(divide='ignore'):
        T = np.where(d == 0, np.pi ** 2 / 3.0, 2.0 * (-1.0) ** np.abs(d) / np.where(d == 0, 1, d) ** 2)
    return T / (2.0 * mass_me * dx ** 2)


def build_potential(x_b, V, interp, walls, pad_b, npts):
    """x_b (bohr, increasing), V (Ha).  Returns grid, potential on it and the interpolant."""
    f = CubicSpline(x_b, V, bc_type='not-a-knot') if interp == 'cubic' else PchipInterpolator(x_b, V)
    xL, xR = x_b[0], x_b[-1]
    g = np.linspace(xL - pad_b, xR + pad_b, npts)
    v = np.empty_like(g)
    inside = (g >= xL) & (g <= xR)
    v[inside] = f(g[inside])
    if walls == 'box':
        v[~inside] = V_WALL
    else:
        d1, d2 = f.derivative(1), f.derivative(2)
        sL, kL = min(float(d1(xL)), 0.0), max(float(d2(xL)), 0.0)
        sR, kR = max(float(d1(xR)), 0.0), max(float(d2(xR)), 0.0)
        left, right = g < xL, g > xR
        dl, dr = xL - g[left], g[right] - xR
        v[left] = f(xL) - sL * dl + 0.5 * kL * dl ** 2
        v[right] = f(xR) + sR * dr + 0.5 * kR * dr ** 2
        v[np.isclose(g, g[0]) | np.isclose(g, g[-1])] = V_WALL
        f.unbound_ends = [e for e, bad in (('N end', d1(xL) > 0), ('O end', d1(xR) < 0)) if bad]
    return g, np.minimum(v, V_WALL), f


def solve(g, v, mass_me, nlevels):
    dx = g[1] - g[0]
    E, C = eigh(kinetic(len(g), dx, mass_me) + np.diag(v), subset_by_index=[0, nlevels - 1])
    return E, C ** 2          # columns: probabilities per grid point (sum 1)


def extrema(f, xL, xR, n=4001):
    xs = np.linspace(xL, xR, n); vs = f(xs)
    mins = [i for i in range(1, n - 1) if vs[i] < vs[i - 1] and vs[i] <= vs[i + 1]]
    maxs = [i for i in range(1, n - 1) if vs[i] > vs[i - 1] and vs[i] >= vs[i + 1]]
    edge = []
    if vs[0] < vs[1]:
        edge.append(0)
    if vs[-1] < vs[-2]:
        edge.append(n - 1)
    return xs, vs, mins, maxs, edge


def run_variant(prof, coord, mass_kind, interp, walls, args, nlevels):
    """One DVR calculation.  prof: dict from make_profile.  Returns a record."""
    x = prof['coords'][coord]
    if np.any(np.diff(x) <= 0):
        return {'coord': coord, 'mass': mass_kind, 'interp': interp, 'walls': walls,
                'skipped': f'coordinate {coord} is not monotonic along the path'}
    if coord == 'mwpath':
        mass_u = 1.0; mass_kind = 'unit (mass-weighted path)'
    elif mass_kind == 'reduced':
        mH, mN, mO = MASS['H'], MASS['N'], MASS['O']
        mass_u = mH * (mN + mO) / (mH + mN + mO)
    else:
        mass_u = MASS['H']
    V = prof['E'] - prof['E'].min()
    x_b = x / BOHR_A
    pad_b = args.pad / BOHR_A * (np.sqrt(MASS['H']) if coord == 'mwpath' else 1.0)
    g, v, f = build_potential(x_b, V, interp, walls, pad_b, args.npts)
    E, P = solve(g, v, mass_u * AMU_ME, nlevels)
    g2, v2, _ = build_potential(x_b, V, interp, walls, pad_b, 2 * args.npts - 1)
    E2, _ = solve(g2, v2, mass_u * AMU_ME, nlevels)
    conv_cm = float(np.abs(E2 - E).max() * HARTREE_CM)
    # profile analysis inside the scanned range
    xs, vs, mins, maxs, edge = extrema(f, x_b[0], x_b[-1])
    vmin = float(vs.min()) if len(vs) else 0.0
    dl = np.interp(g, x_b, prof['delta_A'])                 # delta = r(N-H) - r(H...O); > 0: proton nearer O
    o_side = dl > 0
    outside = (g < x_b[0]) | (g > x_b[-1])
    to_A = BOHR_A / (np.sqrt(MASS['H']) if coord == 'mwpath' else 1.0)   # back to A (mwpath: A at the mass of H)
    # transfer point: delta = 0 (linear in the data)
    xt = None
    sgn = np.sign(prof['delta_A'])
    for k in range(len(x_b) - 1):
        if sgn[k] <= 0 < sgn[k + 1]:
            d0, d1 = prof['delta_A'][k], prof['delta_A'][k + 1]
            xt = float(x_b[k] + (x_b[k + 1] - x_b[k]) * (-d0) / (d1 - d0)); break
    levels = []
    for n in range(nlevels):
        e = float(E[n] - vmin)
        levels.append({'n': n, 'E_cm': round(e * HARTREE_CM, 2), 'E_kcal': round(e * HARTREE_KCAL, 4),
                       'P_O_side': round(float(P[o_side, n].sum()), 5),
                       'P_outside_scan': round(float(P[outside, n].sum()), 5)})
    minima = [{'x_A': round(float(xs[i] * to_A), 4), 'r_HO_A': round(float(np.interp(xs[i], x_b, prof['r_HO_A'])), 4),
               'V_kcal': round(float((vs[i] - vmin) * HARTREE_KCAL), 4),
               'side': 'O' if np.interp(xs[i], x_b, prof['delta_A']) > 0 else 'N'} for i in mins]
    maxima = [{'x_A': round(float(xs[i] * to_A), 4), 'r_HO_A': round(float(np.interp(xs[i], x_b, prof['r_HO_A'])), 4),
               'V_kcal': round(float((vs[i] - vmin) * HARTREE_KCAL), 4)} for i in maxs]
    E0 = float(E[0] - vmin)
    # classical turning points of the ground level around the global minimum
    ig = int(np.argmin(v)); tl = tr = None
    for i in range(ig, -1, -1):
        if v[i] - vmin > E0:
            tl = float(g[i] * to_A); break
    for i in range(ig, len(g)):
        if v[i] - vmin > E0:
            tr = float(g[i] * to_A); break
    gl = min(minima, key=lambda m: m['V_kcal']) if minima else None
    second = [m for m in minima if m is not gl]
    barrier = None
    if gl and second:
        m2 = min(second, key=lambda m: m['V_kcal'])
        lo, hi = sorted([gl['x_A'], m2['x_A']])
        between = [m for m in maxima if lo < m['x_A'] < hi]
        barrier = max(between, key=lambda m: m['V_kcal']) if between else None
    rec = {'coord': coord, 'mass': mass_kind, 'mass_u': mass_u, 'interp': interp, 'walls': walls,
           'npts': args.npts, 'dx_A': float((g[1] - g[0]) * to_A), 'grid_A': [float(g[0] * to_A), float(g[-1] * to_A)],
           'scan_range_A': [float(x_b[0] * to_A), float(x_b[-1] * to_A)], 'grid_convergence_cm': conv_cm,
           'ZPE_cm': round(E0 * HARTREE_CM, 2), 'ZPE_kcal': round(E0 * HARTREE_KCAL, 4),
           'fundamental_cm': round(float(E[1] - E[0]) * HARTREE_CM, 2) if nlevels > 1 else None,
           'levels': levels, 'minima': minima, 'maxima': maxima,
           'boundary_minimum': [('N end' if i == 0 else 'O end') for i in edge],
           'transfer_point_x_A': None if xt is None else round(xt * to_A, 4),
           'V_transfer_point_kcal': None if xt is None else round(float((f(xt) - vmin) * HARTREE_KCAL), 4),
           'ground_turning_points_A': [tl, tr], 'barrier': barrier,
           'second_minimum': (min(second, key=lambda m: m['V_kcal']) if second else None)}
    rec['second_minimum_below_ZPE'] = (None if not second else bool(rec['second_minimum']['V_kcal'] < rec['ZPE_kcal']))
    rec['range_limited'] = bool(levels[0]['P_outside_scan'] > 0.01)
    # checks on the data themselves, so that a spline artefact cannot pass for a minimum
    Vk = V * HARTREE_KCAL
    dmin = [k for k in range(1, len(Vk) - 1) if Vk[k] < Vk[k - 1] and Vk[k] < Vk[k + 1]]
    rec['data_minima'] = [{'index': k, 'r_HO_A': round(float(prof['r_HO_A'][k]), 4),
                           'V_kcal': round(float(Vk[k] - Vk.min()), 4)} for k in dmin]
    rec['spline_undershoot_kcal'] = round(float((V.min() - vmin) * HARTREE_KCAL), 4)
    w = []
    k0 = int(np.argmin(Vk))
    if k0 in (0, len(Vk) - 1):
        w.append(f"the lowest point is at the {'N' if k0 == 0 else 'O'} end of the scan: the scan does not bracket "
                 f"the minimum")
    if rec['spline_undershoot_kcal'] > 0.2:
        w.append(f"the interpolant dips {rec['spline_undershoot_kcal']:.2f} kcal/mol below the lowest point "
                 f"(overshoot; compare the pchip variant)")
    if len(mins) != len(dmin):
        w.append(f"the interpolant has {len(mins)} interior minima, the points themselves {len(dmin)}")
    for e in getattr(f, 'unbound_ends', []):
        w.append(f'the profile does not rise toward the {e}; the extension there is flat up to the hard wall')
    if rec['range_limited']:
        w.append(f"{100 * levels[0]['P_outside_scan']:.1f} % of the ground state lies outside the scanned range")
    if conv_cm > 1.0:
        w.append(f'DVR grid not converged ({conv_cm:.2f} cm-1 between N and 2N-1 points)')
    rec['warnings'] = w
    rec['_plot'] = {'g': g * to_A, 'v': (v - vmin) * HARTREE_KCAL, 'E': (E - vmin) * HARTREE_KCAL, 'P': P,
                    'x': x_b * to_A, 'V': (V - vmin) * HARTREE_KCAL, 'V_data': (V - V.min()) * HARTREE_KCAL,
                    'xt': None if xt is None else xt * to_A}
    return rec


def verdict(rec, method):
    if rec.get('skipped'):
        return rec['skipped']
    z, p0 = rec['ZPE_kcal'], rec['levels'][0]['P_O_side']
    s = f"{METHOD_LABEL.get(method, method)}: ZPE {rec['ZPE_cm']:.0f} cm-1 ({z:.2f} kcal/mol); "
    if rec['second_minimum'] is None:
        s += 'single well: no second minimum along the proton coordinate'
        if rec['V_transfer_point_kcal'] is not None:
            above = rec['V_transfer_point_kcal'] - z
            s += (f"; the transfer point (r(N-H) = r(H...O)) lies {rec['V_transfer_point_kcal']:.2f} kcal/mol above the "
                  f"minimum, {abs(above):.2f} kcal/mol {'above' if above > 0 else 'BELOW'} the ground level")
    else:
        m2 = rec['second_minimum']
        s += (f"double well: second minimum on the {m2['side']} side at r(H...O) {m2['r_HO_A']:.3f} A, "
              f"{m2['V_kcal']:.2f} kcal/mol above the global minimum")
        if rec['barrier']:
            s += f", barrier {rec['barrier']['V_kcal']:.2f} kcal/mol"
        s += f"; it lies {'BELOW' if rec['second_minimum_below_ZPE'] else 'above'} the zero-point level"
    s += f"; ground-state weight on the O side {100 * p0:.1f} %"
    if rec.get('warnings'):
        s += ' | WARNINGS: ' + '; '.join(rec['warnings'])
    return s


# ---------------------------------------------------------------- profiles from the result files
def load_profiles(args, log):
    """{method: {'pts': [...], 'source': ...}} and the proton atoms."""
    profiles, proton = {}, None
    want = None if args.methods in ('all', '') else {m.strip().lower() for m in args.methods.split(',')}
    chains = None if args.chains in ('all', '') else {c.strip() for c in args.chains.split(',')} | {'min'}
    scan_pts = []
    for d in args.scan:
        d = pathlib.Path(d).resolve(); res = json.load(open(d / 'result.json'))
        if bool(res.get('smoke')) != bool(args.smoke):
            raise SystemExit(f'{d}: smoke/production mismatch with this run')
        st = res.get('settings', {})
        pa = st.get('proton_atoms_0based')
        if pa:
            p = (pa['H'], pa['donor_N'], pa['acceptor_O'])
            if proton and p != proton:
                raise SystemExit(f'{d}: different proton atoms {p} vs {proton}')
            proton = p
        for key, it in res['items'].items():
            if it.get('status') != 'done' or it.get('E_Ha') is None:
                continue
            if chains and it.get('chain') not in chains:
                continue
            syms, xyz = read_xyz(d / it['xyz'])
            scan_pts.append({'label': key, 'chain': it.get('chain'), 'r_target': it.get('r_target_A'),
                             'E': it['E_Ha'], 'syms': syms, 'xyz': xyz, 'src': str(d), 'xc': st.get('xc', 'b3lyp'),
                             'disp': st.get('disp')})
    if scan_pts:
        best = {}
        for p in scan_pts:
            t = 'min' if p['r_target'] is None else round(p['r_target'], 3)
            if t not in best or p['E'] < best[t]['E']:
                best[t] = p
        xc = scan_pts[0]['xc'] + (f"-{scan_pts[0]['disp']}" if scan_pts[0]['disp'] else '')
        profiles[xc.lower()] = {'pts': list(best.values()), 'source': 'proton_scan'}
    for d in args.sp:
        d = pathlib.Path(d).resolve(); res = json.load(open(d / 'result.json'))
        if bool(res.get('smoke')) != bool(args.smoke):
            raise SystemExit(f'{d}: smoke/production mismatch with this run')
        geoms = {g['label']: g for g in res.get('settings', {}).get('geometries', [])}
        for key, it in res['items'].items():
            if it.get('status') != 'done' or '@grid' in key or '|' not in key:
                continue
            glab, meth = key.split('|', 1)
            g = geoms.get(glab)
            if g is None or not str(g.get('source', '')).startswith('scan'):
                continue                         # only scan geometries make a profile
            syms, xyz = read_xyz(d / g['local'])
            for m, e in ((meth, it.get('E_Ha')), ('scs-mp2', it.get('E_SCS_MP2_Ha')) if meth == 'dfmp2' else (None, None)):
                if m is None or e is None:
                    continue
                profiles.setdefault(m, {'pts': [], 'source': f'singlepoints {res["job"]}'})['pts'].append(
                    {'label': glab, 'chain': g.get('chain'), 'r_target': g.get('r_target_A'), 'E': e,
                     'syms': syms, 'xyz': xyz, 'src': str(d)})
    if want:
        profiles = {m: p for m, p in profiles.items() if m in want}
    if proton is None:
        case = json.load(open(ROOT / 'molecules' / 'btpc' / 'cases' / 'btpc_gas.json'))
        proton = (case['proton']['H'], case['proton']['donor_N'], case['proton']['acceptor_O'])
        log(f'proton atoms not in the scan settings; using btpc_gas {proton}')
    return profiles, proton


def make_profile(pts, proton):
    H, N, O = proton
    pts = sorted(pts, key=lambda p: -np.linalg.norm(p['xyz'][H] - p['xyz'][O]))
    c, info = path_coordinates(pts, H, N, O)
    return {'labels': [p['label'] for p in pts], 'E': np.array([p['E'] for p in pts]), 'coords': c,
            'r_HO_A': info['r_HO_A'], 'r_NH_A': info['r_NH_A'], 'delta_A': info['delta_A'],
            'm_eff_path_u': info['m_eff_path_u'], 'm_eff_segments_u': info['m_eff_segments_u'], 'i_min': info['i_min']}


# ---------------------------------------------------------------- figures
def _colour(method):
    """Fixed colour per method (follows the entity, not its rank among the plotted series)."""
    return SERIES[METHOD_ORDER.index(method) % len(SERIES)] if method in METHOD_ORDER else INK2


def _mpl():
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        return plt
    except Exception:
        return None


def _style(ax):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.yaxis.label.set_color(INK); ax.xaxis.label.set_color(INK)
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)


def coord_label(coord):
    return {'xH': 'proton displacement along the path, x$_H$ (Å)',
            'half_delta': '[r(N–H) − r(H···O)] / 2 (Å)',
            'rHO': '−Δr(H···O) (Å)',
            'mwpath': 'mass-weighted path / √m$_H$ (Å)'}[coord]


def plot_levels(path, rec, method, labs, nshow):
    d = rec['_plot']; ymax = max(float(d['V'].max()), float(d['E'][min(nshow, len(d['E'])) - 1]) * 1.15) + 1.0
    plt = _mpl()
    if plt is None:
        return svg_levels(path.with_suffix('.svg'), rec, method, nshow, ymax)
    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=200)
    _style(ax)
    ax.plot(d['g'], np.where(d['v'] < ymax, d['v'], np.nan), color=SERIES[0], lw=2, zorder=3)
    ax.plot(d['x'], d['V'], 'o', ms=5, mfc=SERIES[0], mec='white', mew=1.2, zorder=4)
    scale = 0.18 * ymax / max(float(d['P'][:, :nshow].max()), 1e-12)
    for n in range(min(nshow, len(d['E']))):
        e = float(d['E'][n]); inside = d['v'] <= e
        if inside.any():
            ax.hlines(e, d['g'][inside].min(), d['g'][inside].max(), color=INK2, lw=1, zorder=2)
        ax.fill_between(d['g'], e, e + scale * d['P'][:, n], color=SERIES[0], alpha=0.18, lw=0, zorder=1)
        ax.text(d['g'][-1], e, f' v={n}', va='center', ha='left', fontsize=8, color=INK2)
    if d['xt'] is not None:
        ax.axvline(d['xt'], color=INK2, lw=1, ls=(0, (4, 3)))
        ax.text(d['xt'], ymax * 0.97, f" r({labs['N']}–H) = r(H···{labs['O']})", fontsize=8,
                color=INK2, va='top')
    ax.set_xlim(d['g'][0], d['g'][-1]); ax.set_ylim(-0.5, ymax)
    ax.set_xlabel(coord_label(rec['coord'])); ax.set_ylabel('energy above the minimum (kcal/mol)')
    ax.set_title(f"{METHOD_LABEL.get(method, method)}: proton levels on the relaxed scan (sinc-DVR)",
                 fontsize=10, color=INK, loc='left')
    ax.text(0.01, -0.2, f"{labs['N']} side ←   → {labs['O']} side;  ZPE {rec['ZPE_cm']:.0f} cm⁻¹, "
            f"P(O side, v=0) {100 * rec['levels'][0]['P_O_side']:.1f} %", transform=ax.transAxes, fontsize=8,
            color=INK2)
    fig.tight_layout(); fig.savefig(path, facecolor='white'); plt.close(fig)
    return path


def plot_profiles(path, results, primary_coord):
    plt = _mpl()
    series = [(m, results[m]['_plot']) for m in METHOD_ORDER + sorted(set(results) - set(METHOD_ORDER))
              if m in results and '_plot' in results[m]]
    if not series:
        return None
    if plt is None:
        return svg_profiles(path.with_suffix('.svg'), series)
    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=200)
    _style(ax)
    for m, d in series:
        ax.plot(d['x'], d['V_data'], '-o', color=_colour(m), lw=2, ms=4, mec='white', mew=0.8,
                label=METHOD_LABEL.get(m, m))
    if len(series) <= 4:                       # direct labels at the right end, nudged apart
        ends = sorted(((float(d['V_data'][-1]), m, float(d['x'][-1])) for m, d in series))
        span = max(float(d['V_data'].max()) for _, d in series) or 1.0
        ys = []
        for y, m, x in ends:
            y = max(y, ys[-1] + 0.045 * span) if ys else y
            ys.append(y)
            ax.text(x, y, ' ' + METHOD_LABEL.get(m, m), fontsize=8, color=INK, va='center')
    ax.set_xlabel(coord_label(primary_coord)); ax.set_ylabel('energy above the lowest point (kcal/mol)')
    ax.set_title('Proton-transfer profiles (B3LYP relaxed scan geometries)', fontsize=10, color=INK, loc='left')
    leg = ax.legend(fontsize=8, frameon=False, loc='upper left')
    for t in leg.get_texts():
        t.set_color(INK)
    fig.tight_layout(); fig.savefig(path, facecolor='white'); plt.close(fig)
    return path


def _svg(path, w, h, body):
    path.write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" font-family="sans-serif" '
                    f'font-size="11"><rect width="100%" height="100%" fill="white"/>{body}</svg>\n')
    return path


def _axes(x0, x1, y0, y1, W=640, Hh=420, m=(60, 20, 20, 50)):
    L, R, T, B = m
    sx = lambda x: L + (x - x0) / (x1 - x0) * (W - L - R)
    sy = lambda y: Hh - B - (y - y0) / (y1 - y0) * (Hh - T - B)
    body = [f'<line x1="{L}" y1="{Hh - B}" x2="{W - R}" y2="{Hh - B}" stroke="{INK2}"/>',
            f'<line x1="{L}" y1="{T}" x2="{L}" y2="{Hh - B}" stroke="{INK2}"/>']
    for t in np.linspace(x0, x1, 6):
        body.append(f'<text x="{sx(t):.1f}" y="{Hh - B + 15}" text-anchor="middle" fill="{INK2}">{t:.2f}</text>')
    for t in np.linspace(y0, y1, 6):
        body.append(f'<text x="{L - 5}" y="{sy(t) + 4:.1f}" text-anchor="end" fill="{INK2}">{t:.1f}</text>')
    return sx, sy, body


def svg_levels(path, rec, method, nshow, ymax):
    d = rec['_plot']; sx, sy, body = _axes(d['g'][0], d['g'][-1], 0, ymax)
    m = d['v'] < ymax
    body.append('<polyline fill="none" stroke="%s" stroke-width="2" points="%s"/>' % (
        SERIES[0], ' '.join(f'{sx(a):.1f},{sy(b):.1f}' for a, b in zip(d['g'][m], d['v'][m]))))
    for a, b in zip(d['x'], d['V']):
        body.append(f'<circle cx="{sx(a):.1f}" cy="{sy(b):.1f}" r="4" fill="{SERIES[0]}" stroke="white"/>')
    for n in range(min(nshow, len(d['E']))):
        e = float(d['E'][n]); inside = d['v'] <= e
        if inside.any():
            body.append(f'<line x1="{sx(d["g"][inside].min()):.1f}" x2="{sx(d["g"][inside].max()):.1f}" '
                        f'y1="{sy(e):.1f}" y2="{sy(e):.1f}" stroke="{INK2}"/>')
    body.append(f'<text x="60" y="14" fill="{INK}">{METHOD_LABEL.get(method, method)}: proton levels (sinc-DVR), '
                f'kcal/mol vs {rec["coord"]} (A)</text>')
    return _svg(path, 640, 420, ''.join(body))


def svg_profiles(path, series):
    x0 = min(float(d['x'].min()) for _, d in series); x1 = max(float(d['x'].max()) for _, d in series)
    y1 = max(float(d['V_data'].max()) for _, d in series) + 1
    sx, sy, body = _axes(x0, x1, 0, y1)
    for m, d in series:
        body.append('<polyline fill="none" stroke="%s" stroke-width="2" points="%s"/>' % (
            _colour(m), ' '.join(f'{sx(a):.1f},{sy(b):.1f}' for a, b in zip(d['x'], d['V_data']))))
        body.append(f'<text x="{sx(d["x"][-1]) + 4:.1f}" y="{sy(d["V_data"][-1]):.1f}" fill="{INK}">'
                    f'{METHOD_LABEL.get(m, m)}</text>')
    return _svg(path, 640, 420, ''.join(body))


# ---------------------------------------------------------------- analytic self-tests (--smoke / --selftest)
def selftest(job, args):
    res, ok = {}, True
    mH = MASS['H'] * AMU_ME
    # 1. harmonic oscillator, 3000 cm-1, mass of H
    w = 3000 / HARTREE_CM; k = mH * w ** 2
    g = np.linspace(-1.2, 1.2, args.npts) / BOHR_A
    E, P = solve(g, np.minimum(0.5 * k * g ** 2, V_WALL), mH, 5)
    exact = w * (np.arange(5) + 0.5)
    err = float(np.abs(E - exact).max() / exact[0])
    res['harmonic'] = {'E_cm': (E * HARTREE_CM).round(4).tolist(), 'exact_cm': (exact * HARTREE_CM).round(4).tolist(),
                       'max_rel_err': err, 'pass': err < 1e-8}
    # 2. Morse oscillator (omega 3000 cm-1, D 0.18 Ha), exact levels
    D = 0.18; a = w * np.sqrt(mH / (2 * D))
    gm = np.linspace(-0.8, 3.0, args.npts) / BOHR_A
    Em, _ = solve(gm, np.minimum(D * (1 - np.exp(-a * gm)) ** 2, V_WALL), mH, 5)
    n = np.arange(5) + 0.5
    exm = w * n - (w * n) ** 2 / (4 * D)
    errm = float(np.abs(Em - exm).max() / exm[0])
    res['morse'] = {'E_cm': (Em * HARTREE_CM).round(4).tolist(), 'exact_cm': (exm * HARTREE_CM).round(4).tolist(),
                    'max_rel_err': errm, 'pass': errm < 1e-7}
    # 3. the profile path: a harmonic potential sampled at 14 irregular points, splined and extrapolated
    xs = np.array([-0.30, -0.22, -0.15, -0.10, -0.06, -0.03, 0.0, 0.03, 0.07, 0.12, 0.18, 0.25, 0.33, 0.42])
    Vs = 0.5 * k * (xs / BOHR_A) ** 2
    gg, vv, _ = build_potential(xs / BOHR_A, Vs, 'cubic', 'extrapolate', 0.3 / BOHR_A, args.npts)
    Es, Ps = solve(gg, vv, mH, 3)
    errs = float(np.abs(Es - exact[:3]).max() / exact[0])
    res['sampled_harmonic'] = {'E_cm': (Es * HARTREE_CM).round(3).tolist(), 'max_rel_err': errs, 'pass': errs < 1e-3,
                               'P_outside_scan_v0': float(Ps[(gg < xs[0] / BOHR_A) | (gg > xs[-1] / BOHR_A), 0].sum())}
    # 4. symmetric double well: ground state shared 50/50, tunnelling splitting > 0
    b = 0.35 / BOHR_A; Vb = 6 / HARTREE_KCAL
    gd = np.linspace(-1.0, 1.0, args.npts) / BOHR_A
    Ed, Pd = solve(gd, np.minimum(Vb * ((gd / b) ** 2 - 1) ** 2, V_WALL), mH, 4)
    wl = float(Pd[gd < -1e-12, 0].sum() + 0.5 * Pd[np.abs(gd) <= 1e-12, 0].sum())     # centre point split
    res['double_well'] = {'barrier_kcal': 6.0, 'E_cm': (Ed * HARTREE_CM).round(3).tolist(),
                          'splitting_cm': float((Ed[1] - Ed[0]) * HARTREE_CM), 'P_left_v0': wl,
                          'pass': abs(wl - 0.5) < 1e-6 and Ed[1] > Ed[0]}
    # 5. the whole analysis path on an asymmetric double well sampled like a scan (14 points): two minima, the
    #    barrier, the O-side weight and the verdict, against a direct DVR on the analytic potential
    b5, Vb5, c5 = 0.30, 8 / HARTREE_KCAL, 3 / HARTREE_KCAL          # wells at -/+0.30 A, O well ~3.5 kcal/mol higher
    fV = lambda xa: Vb5 * ((xa / b5) ** 2 - 1) ** 2 + c5 * (xa / b5)
    xs5 = np.array([-0.62, -0.52, -0.44, -0.37, -0.31, -0.25, -0.18, -0.10, 0.0, 0.10, 0.20, 0.28, 0.34, 0.42, 0.52])
    prof5 = {'E': fV(xs5), 'coords': {'xH': xs5}, 'delta_A': xs5, 'r_HO_A': 1.44 - xs5}
    r5 = run_variant(prof5, 'xH', 'H', 'cubic', 'extrapolate', args, 3)
    gx = np.linspace(-1.2, 1.2, args.npts)
    Ex, Px = solve(gx / BOHR_A, np.minimum(fV(gx) - fV(gx).min(), V_WALL), mH, 3)
    e0x = float(Ex[0] * HARTREE_CM)
    res['sampled_double_well'] = {
        'minima': r5['minima'], 'barrier': r5['barrier'], 'ZPE_cm': r5['ZPE_cm'], 'ZPE_direct_cm': round(e0x, 2),
        'P_O_v0': r5['levels'][0]['P_O_side'], 'P_O_v0_direct': round(float(Px[gx > 0, 0].sum()), 5),
        'second_minimum_below_ZPE': r5['second_minimum_below_ZPE'], 'verdict': verdict(r5, 'test'),
        'pass': (len(r5['minima']) == 2 and r5['barrier'] is not None and abs(r5['ZPE_cm'] - e0x) / e0x < 0.01
                 and abs(r5['levels'][0]['P_O_side'] - float(Px[gx > 0, 0].sum())) < 0.01)}
    for kk, v in res.items():
        ok &= bool(v['pass'])
        job.log(f"selftest {kk}: {'PASS' if v['pass'] else 'FAIL'} "
                + ', '.join(f'{a}={b}' for a, b in v.items() if a not in ('pass',)))
    job.r['selftest'] = res; job.save()
    return ok


# ---------------------------------------------------------------- main
def parser():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    p.add_argument('--job', required=True); p.add_argument('--out', required=True)
    p.add_argument('--smoke', action='store_true', help='analytic self-tests, then the given (smoke) profiles')
    p.add_argument('--selftest', action='store_true', help='run the analytic self-tests also in production mode')
    p.add_argument('--scan', action='append', default=[], help='proton_scan.py output directory (repeatable)')
    p.add_argument('--sp', action='append', default=[], help='singlepoints.py output directory (repeatable)')
    p.add_argument('--methods', default='all', help="profiles to use: all, or a comma list (b3lyp, wb97x-d3bj, ...)")
    p.add_argument('--chains', default='all', help="scan chains to use: all, or e.g. up,down")
    p.add_argument('--coord', choices=COORDS, default='xH')
    p.add_argument('--mass', choices=('H', 'reduced'), default='H')
    p.add_argument('--interp', choices=('cubic', 'pchip'), default='cubic')
    p.add_argument('--walls', choices=('extrapolate', 'box'), default='extrapolate')
    p.add_argument('--npts', type=int, default=801, help='DVR grid points (checked against 2N-1)')
    p.add_argument('--pad', type=float, default=0.30, help='extension beyond the scanned range (A)')
    p.add_argument('--nlevels', type=int, default=6)
    p.add_argument('--keep-scratch', action='store_true', help='(accepted for the common interface; unused)')
    return p


SENSITIVITY = [('xH', 'H', 'pchip', 'extrapolate'), ('xH', 'H', 'cubic', 'box'), ('xH', 'reduced', 'cubic', 'extrapolate'),
               ('half_delta', 'H', 'cubic', 'extrapolate'), ('rHO', 'H', 'cubic', 'extrapolate'),
               ('mwpath', 'unit', 'cubic', 'extrapolate')]


def run(job, args):
    job.r['settings'] = {'coord': args.coord, 'mass': args.mass, 'interp': args.interp, 'walls': args.walls,
                         'npts': args.npts, 'pad_A': args.pad, 'nlevels': args.nlevels, 'chains': args.chains,
                         'methods': args.methods, 'sensitivity': SENSITIVITY,
                         'masses_u': MASS, 'constants': {'bohr_A': BOHR_A, 'u_me': AMU_ME, 'hartree_cm': HARTREE_CM,
                                                         'hartree_kcal': HARTREE_KCAL}}
    job.save()
    ok = True
    if args.smoke or args.selftest:
        ok = selftest(job, args)
    if not args.scan and not args.sp:
        if not args.smoke and not args.selftest:
            raise SystemExit('give --scan and/or --sp (or --selftest)')
        return ok
    for d in args.scan + args.sp:
        job.input_file(f'result:{d}', pathlib.Path(d) / 'result.json')
    profiles, proton = load_profiles(args, job.log)
    job.r['settings']['proton_atoms_0based'] = {'H': proton[0], 'donor_N': proton[1], 'acceptor_O': proton[2]}
    labs = {'N': 'N4', 'O': 'O1'}
    for d in args.scan:
        st = json.load(open(pathlib.Path(d) / 'result.json')).get('settings', {})
        if st.get('labels'):
            labs = st['labels']; break
    rows, results = [], {}
    for method in METHOD_ORDER + sorted(set(profiles) - set(METHOD_ORDER)):
        if method not in profiles:
            continue
        pts = profiles[method]['pts']
        if len(pts) < 4:
            job.log(f'{method}: only {len(pts)} points; at least 4 are needed; skipped')
            continue
        prof = make_profile(pts, proton)
        t = time.time()
        main = run_variant(prof, args.coord, args.mass, args.interp, args.walls, args, args.nlevels)
        if main.get('skipped'):
            raise RuntimeError(f"{method}: {main['skipped']}")
        sens = []
        for c, mk, ip, wl in SENSITIVITY:
            if (c, mk, ip, wl) == (args.coord, args.mass, args.interp, args.walls):
                continue
            s = run_variant(prof, c, mk, ip, wl, args, min(args.nlevels, 3))
            s.pop('_plot', None)
            sens.append({k: s.get(k) for k in ('coord', 'mass', 'mass_u', 'interp', 'walls', 'skipped', 'ZPE_cm',
                                               'fundamental_cm', 'second_minimum', 'second_minimum_below_ZPE',
                                               'range_limited', 'grid_convergence_cm', 'warnings')}
                        | ({'P_O_side_v0': s['levels'][0]['P_O_side'], 'P_outside_v0': s['levels'][0]['P_outside_scan']}
                           if not s.get('skipped') else {}))
        plot = main.pop('_plot')
        results[method] = dict(main, _plot=plot)
        figs = []
        f = plot_levels(job.out / f"dvr_{method}.png", dict(main, _plot=plot), method, labs, min(args.nlevels, 4))
        if f:
            figs.append(pathlib.Path(f).name)
        rec = dict(main)
        rec.update({'method': method, 'source': profiles[method]['source'], 'n_points': len(pts),
                    'points': [{'label': l, 'r_HO_A': round(float(r), 4), 'r_NH_A': round(float(n), 4),
                                'x_A': round(float(x), 5), 'dE_kcal': round(float((e - prof['E'].min()) * HARTREE_KCAL), 4)}
                               for l, r, n, x, e in zip(prof['labels'], prof['r_HO_A'], prof['r_NH_A'],
                                                        prof['coords'][args.coord], prof['E'])],
                    'm_eff_path_u': prof['m_eff_path_u'],
                    'm_eff_segments_u': [round(float(v), 3) for v in prof['m_eff_segments_u']],
                    'sensitivity': sens, 'figures': figs, 'wall_s': round(time.time() - t, 2)})
        rec['verdict'] = verdict(main, method)
        job.put(method, rec)
        job.log(rec['verdict'])
        job.log(f"  levels (cm-1): {[l['E_cm'] for l in main['levels']]}; P(O side): "
                f"{[l['P_O_side'] for l in main['levels']]}; grid check {main['grid_convergence_cm']:.1e} cm-1; "
                f"m_eff of the path {prof['m_eff_path_u']:.2f} u")
        for s in sens:
            job.log(f"  sensitivity {s['coord']}/{s['mass']}/{s['interp']}/{s['walls']}: "
                    + (s['skipped'] if s.get('skipped') else
                       f"ZPE {s['ZPE_cm']:.0f} cm-1, fundamental {s['fundamental_cm']:.0f} cm-1, "
                       f"P(O, v=0) {100 * s['P_O_side_v0']:.1f} %, second minimum below ZPE: {s['second_minimum_below_ZPE']}"))
        for l in main['levels']:
            rows.append((method, l['n'], l['E_cm'], l['E_kcal'], l['P_O_side'], l['P_outside_scan']))
    with open(job.out / 'levels.csv', 'w') as f:
        f.write('method,n,E_cm,E_kcal,P_O_side,P_outside_scan\n')
        for r in rows:
            f.write(','.join(str(v) for v in r) + '\n')
    f = plot_profiles(job.out / 'profiles.png', results, args.coord)
    job.r['summary'] = {m: {k: job.get(m).get(k) for k in ('verdict', 'ZPE_cm', 'fundamental_cm',
                                                            'second_minimum_below_ZPE', 'range_limited')}
                        for m in results}
    job.r['figures'] = sorted(p.name for p in job.out.glob('*.png')) + sorted(p.name for p in job.out.glob('*.svg'))
    job.save()
    return ok


def main():
    args = parser().parse_args()
    job = Job(args, __file__)

    def _h(signum, frame):
        raise KeyboardInterrupt(f'signal {signum}')
    signal.signal(signal.SIGTERM, _h)
    try:
        ok = run(job, args)
        job.finish(ok is not False)
        sys.exit(0 if ok is not False else 1)
    except SystemExit:
        raise
    except BaseException as e:
        tb = traceback.format_exc()
        job.log(f'ERROR {type(e).__name__}: {e}\n{tb}')
        job.finish(False, {'type': type(e).__name__, 'message': str(e), 'traceback': tb})
        sys.exit(1)


if __name__ == '__main__':
    main()
