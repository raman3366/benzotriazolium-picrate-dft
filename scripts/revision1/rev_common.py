"""
Common setup for the BTPC revision calculations (dft/revision/*.py).

Import this module FIRST in every revision script, before numpy or pyscf: it pins every thread pool from the
environment before any BLAS library is loaded (project notes, not deposited: "Pin every thread pool").

Interface contract shared by all revision scripts
    python dft/revision/<name>.py --job <job_id> --out <dir> [options] [--smoke]
    env THREADS          threads for PySCF, OpenMP and every BLAS pool (default: OMP_NUM_THREADS, else all CPUs)
    env PYSCF_TMPDIR     PySCF scratch (default <out>/scratch, removed after a successful run); never a tmpfs
    env MAX_MEMORY_MB    PySCF max_memory (default min(6000, 70 % of RAM))
    writes <out>/result.json (inputs, settings incl. grid and Lebedev order, software versions, energies,
    timings) and <out>/log.txt; restartable (finished items are kept in result.json and skipped, SCF
    checkpoints in <out>/chk/ are reused); exits non-zero on failure.

Level of theory of record (dft/common.py make_mf as launched by run_case.sh, PySCF 2.14.0):
    RKS, xc 'b3lyp' (libxc HYB_GGA_XC_B3LYP, VWN-RPA = Gaussian's B3LYP), basis 6-311++G(d,p) (spherical) with
    the case's per-element overrides, density fitting via mf.density_fit() with PySCF's default auxiliary basis
    (aug-cc-pVTZ-JKFIT for 6-311++G(d,p), verified 20 Nov 2026), integration grid level 4, max_cycle 150,
    IEF-PCM water eps 78.3553 with Lebedev order 29 and PySCF's other PCM defaults (SWIG surface,
    modified-Bondi radii x 1.2, no probe radius).  dft/common.py defaults to grid 3 / Lebedev 17; run_case.sh
    exported 4 / 29 for every published calculation, so this module sets 4 / 29 explicitly and records them.
"""
import os, sys, json, time, socket, hashlib, platform, pathlib, datetime, traceback, contextlib, argparse, shutil

# ---------------------------------------------------------------- threads (before numpy / pyscf are imported)
def _threads_from_env():
    for k in ('THREADS', 'OMP_NUM_THREADS', 'DFT_THREADS'):
        v = os.environ.get(k, '').strip()
        if v.isdigit() and int(v) > 0:
            return int(v)
    return os.cpu_count() or 1


def _check_shadowing():
    """A file in dft/revision/ named like a standard-library module (queue.py, json.py, ...) shadows it for every
    script run from here.  Refuse to run rather than fail obscurely later (geomeTRIC needs the real `queue`)."""
    here = pathlib.Path(__file__).resolve().parent
    names = getattr(sys, 'stdlib_module_names', ())
    bad = sorted(p.name for p in here.glob('*.py') if p.stem in names)
    if bad:
        raise SystemExit(f'[rev_common] {bad} in {here} shadow standard-library modules; rename them')


_check_shadowing()
THREADS = _threads_from_env()
if any(m in sys.modules for m in ('numpy', 'pyscf', 'scipy')):
    print('[rev_common] WARNING: numpy/scipy/pyscf were imported before rev_common; '
          'BLAS thread pools may not be pinned', file=sys.stderr)
for _v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
           'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = str(THREADS)

import numpy as np                      # noqa: E402
import pyscf                            # noqa: E402
from pyscf import gto, dft, lib         # noqa: E402

lib.num_threads(THREADS)

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent      # repository root
HARTREE_EV = 27.211386245988
HARTREE_KCAL = 627.509474
BOHR_A = 0.529177210903
IST = datetime.timezone(datetime.timedelta(hours=5, minutes=30), 'IST')

# ---------------------------------------------------------------- level of theory of record
XC = 'b3lyp'
BASIS = '6-311++g(d,p)'
GRID_LEVEL = 4                 # published setting (run_case.sh); dft/common.py default is 3
PCM_LEBEDEV_ORDER = 29         # published setting (run_case.sh); dft/common.py default is 17
PCM_METHOD = 'IEF-PCM'
SOLVENT_EPS = {'water': 78.3553, 'dmso': 46.826, 'methanol': 32.613, 'ethanol': 24.852}
MAX_CYCLE = 150
CONV_TOL_SP = 1e-10            # single points (fmo.py, frequencies.py, _fd_worker.py)
CONV_TOL_OPT_GAS = 1e-9        # optimize.py, gas
CONV_TOL_OPT_PCM = 1e-8        # optimize.py, solvent
STORED = {                     # papers/btpc/data/facts.json gas.e_total / pcm.e_total (opt_summary.json E_final_Ha)
    'btpc_gas': {'E_Ha': -1317.20823093777, 'xyz': 'molecules/btpc/results/btpc_gas/opt.xyz'},
    'btpc_pcm': {'E_Ha': -1317.25368856, 'xyz': 'molecules/btpc/results/btpc_pcm/opt.xyz'},
}


def _ram_mb():
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                if line.startswith('MemTotal:'):
                    return int(line.split()[1]) // 1024
    except OSError:
        pass
    return 8000


MAX_MEMORY_MB = int(os.environ.get('MAX_MEMORY_MB') or min(6000, int(0.7 * _ram_mb())))


def now_ist():
    return datetime.datetime.now(IST).isoformat(timespec='seconds')


# ---------------------------------------------------------------- scratch
def setup_scratch(out_dir):
    """PySCF scratch: PYSCF_TMPDIR if given, else <out>/scratch.  Returns (path, created_here)."""
    given = os.environ.get('PYSCF_TMPDIR')
    path = pathlib.Path(given) if given else pathlib.Path(out_dir) / 'scratch'
    path.mkdir(parents=True, exist_ok=True)
    lib.param.TMPDIR = str(path)
    import tempfile
    tempfile.tempdir = str(path)
    os.environ['TMPDIR'] = str(path)
    return path, not given


# ---------------------------------------------------------------- cases and geometries
def find_case_file(name):
    """molecules/<id>/cases/<name>.json (same search as dft/common.py)."""
    p = pathlib.Path(name)
    if p.suffix == '.json' and p.exists():
        return p.resolve()
    for c in sorted((ROOT / 'molecules').glob(f'*/cases/{name}.json')):
        return c
    raise FileNotFoundError(f'no case file for {name!r} under molecules/*/cases/')


def load_case(name):
    path = find_case_file(name)
    case = json.load(open(path))
    case['name'] = path.stem
    case['case_file'] = str(path.relative_to(ROOT)) if str(path).startswith(str(ROOT)) else str(path)
    return case


def read_xyz(path):
    lines = open(path).read().strip().splitlines()
    n = int(lines[0])
    syms, xyz = [], []
    for line in lines[2:2 + n]:
        p = line.split()
        syms.append(p[0]); xyz.append([float(v) for v in p[1:4]])
    return syms, np.array(xyz)


def write_xyz(path, syms, xyz, comment=''):
    tmp = pathlib.Path(str(path) + '.tmp')
    with open(tmp, 'w') as f:
        f.write(f'{len(syms)}\n{comment}\n')
        for s, (x, y, z) in zip(syms, xyz):
            f.write(f'{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n')
    os.replace(tmp, path)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fragments(case):
    """{'cation': [...], 'anion': [...]} from the case's fragment names (btpc: 'benzotriazolium cation',
    'picrate anion').  Any other fragments are returned under their own names."""
    return {frag_key(name): [int(a) for a in atoms] for name, atoms in (case.get('fragments') or {}).items()}


def frag_key(name):
    """'benzotriazolium cation' -> 'cation', 'picrate anion' -> 'anion'; other names unchanged."""
    low = name.lower()
    return 'cation' if 'cation' in low else 'anion' if ('anion' in low or 'picrate' in low) else name


# ---------------------------------------------------------------- level of theory
def basis_spec(basis=BASIS, overrides=None):
    ov = dict(overrides or {})
    return basis if not ov else dict({'default': basis}, **ov)


def build_mol(geom, charge=0, spin=0, basis=BASIS, basis_overrides=None, verbose=3, max_memory=None, output=None):
    """geom: path to an XYZ file (read by PySCF, exactly as dft/common.py does) or a list of
    (symbol, (x, y, z)) in Angstrom.  spin = 2S (number of unpaired electrons).  output: PySCF log file
    (appended; default stdout)."""
    atom = str(geom) if isinstance(geom, (str, pathlib.Path)) else geom
    mol = gto.M(atom=atom, basis=basis_spec(basis, basis_overrides), charge=charge, spin=spin,
                verbose=0, max_memory=max_memory or MAX_MEMORY_MB, symmetry=False)
    mol.verbose = verbose
    if output:
        mol.stdout = open(output, 'a', buffering=1)
    if verbose >= 4:
        mol.dump_input()
    return mol


def mol_for_case(case, geom=None, charge=None, spin=None, basis=None, verbose=3, output=None):
    """Mole for a case at its stored geometry (or `geom`), with the case's basis overrides."""
    g = geom if geom is not None else ROOT / case['xyz']
    return build_mol(g, case['charge'] if charge is None else charge, case['spin'] if spin is None else spin,
                     basis=basis or BASIS, basis_overrides=case.get('basis_overrides'), verbose=verbose,
                     output=output)


def solvent_eps(solvent):
    """None -> gas phase; a name ('water') -> its eps; a number -> that eps (permittivity series)."""
    if solvent is None or (isinstance(solvent, str) and solvent.lower() in ('', 'none', 'gas')):
        return None
    if isinstance(solvent, (int, float)):
        return float(solvent)
    try:
        return float(solvent)
    except ValueError:
        return SOLVENT_EPS[solvent.lower()]


def make_mf(mol, xc=XC, disp=None, solvent=None, conv_tol=CONV_TOL_SP, grid_level=GRID_LEVEL,
            lebedev_order=PCM_LEBEDEV_ORDER, restricted=None, density_fit=True, chkfile=None,
            max_cycle=MAX_CYCLE):
    """RKS/UKS exactly as dft/common.py make_mf builds the published calculations.

    xc: libxc name ('b3lyp', 'cam-b3lyp', 'pbe0', 'm062x', 'wb97x-d3bj', ...).  disp: None or 'd3bj' (appended to
    xc as PySCF's '<xc>-d3bj', handled by pyscf-dispersion).  solvent: None, 'water' or a dielectric constant.
    restricted: None -> RKS for closed shells, UKS otherwise; False forces UKS.
    """
    if restricted is None:
        restricted = mol.spin == 0
    cls = dft.RKS if restricted else dft.UKS
    mf = cls(mol)
    if density_fit:
        mf = mf.density_fit()              # PySCF's default auxiliary basis, as in dft/common.py
    mf.xc = xc if not disp else f'{xc}-{disp}'
    mf.grids.level = int(grid_level)
    mf.conv_tol = conv_tol
    mf.max_cycle = max_cycle
    if chkfile:
        mf.chkfile = str(chkfile)
    eps = solvent_eps(solvent)
    if eps is not None:
        mf = mf.PCM()
        mf.with_solvent.method = PCM_METHOD
        mf.with_solvent.eps = eps
        mf.with_solvent.lebedev_order = int(lebedev_order)
    return mf


def settings_of(mf):
    """Machine-readable record of the level of theory actually set on an SCF object."""
    mol = mf.mol
    s = {'method': type(mf).__name__, 'xc': getattr(mf, 'xc', None), 'disp': getattr(mf, 'disp', None),
         'basis': mol.basis if isinstance(mol.basis, str) else dict(mol.basis), 'cart': bool(mol.cart),
         'nao': int(mol.nao), 'nelectron': int(mol.nelectron), 'charge': int(mol.charge), 'spin': int(mol.spin),
         'grid_level': int(mf.grids.level) if hasattr(mf, 'grids') else None,
         'conv_tol': mf.conv_tol, 'max_cycle': mf.max_cycle, 'max_memory_mb': mol.max_memory,
         'density_fit': hasattr(mf, 'with_df'), 'solvent': None}
    if hasattr(mf, 'with_df'):
        s['auxbasis'] = mf.with_df.auxbasis       # None = PySCF default (resolved name recorded after build)
        try:
            from pyscf.df import addons
            s['auxbasis_default'] = addons.make_auxbasis(mol)
        except Exception as e:                   # pragma: no cover
            s['auxbasis_default'] = f'? {e}'
    ws = getattr(mf, 'with_solvent', None)
    if ws is not None:
        s['solvent'] = {'model': ws.method, 'eps': ws.eps, 'lebedev_order': ws.lebedev_order,
                        'vdw_scale': ws.vdw_scale, 'r_probe': ws.r_probe,
                        'radii_table': None if ws.radii_table is None else 'custom',
                        'surface_discretization_method': ws.surface_discretization_method,
                        'equilibrium_solvation': ws.equilibrium_solvation}
    return s


def runtime_sizes(mf):
    """Sizes known after a kernel: grid points, auxiliary functions, PCM surface points, retained MOs."""
    r = {}
    try:
        r['n_grid_points'] = int(mf.grids.weights.size)
    except Exception:
        pass
    try:
        r['n_aux'] = int(mf.with_df.auxmol.nao)
    except Exception:
        pass
    try:
        r['n_pcm_points'] = int(mf.with_solvent.surface['grid_coords'].shape[0])
    except Exception:
        pass
    try:
        C = mf.mo_coeff[0] if isinstance(mf.mo_coeff, (tuple, list)) or np.ndim(mf.mo_coeff) == 3 else mf.mo_coeff
        r['n_mo'] = int(C.shape[1])
    except Exception:
        pass
    return r


def run_scf(mf, dm0=None, chk_guess=None, log=None):
    """kernel() with an optional initial density from a checkpoint; raises if the SCF does not converge.
    Returns wall seconds; mf._cycles holds the number of SCF cycles."""
    if dm0 is None and chk_guess and os.path.exists(chk_guess):
        try:
            dm0 = mf.from_chk(str(chk_guess))
            if log:
                log(f'initial density from {chk_guess}')
        except Exception as e:
            if log:
                log(f'could not read {chk_guess} ({e}); default guess')
            dm0 = None
    cycles = []
    old_cb = mf.callback
    mf.callback = lambda env: (cycles.append(env.get('cycle')), old_cb(env) if old_cb else None)
    t = time.time()
    mf.kernel(dm0=dm0)
    wall = time.time() - t
    mf.callback = old_cb
    mf._cycles = len(cycles)
    if not mf.converged:
        raise RuntimeError(f'SCF did not converge in {mf.max_cycle} cycles (E = {mf.e_tot})')
    return wall


def compare_with_common(mol, solvent=None):
    """Build dft/common.py's make_mf (with GRID_LEVEL=4 / PCM_LEBEDEV_ORDER=29, as run_case.sh exported them) and
    this module's make_mf for the same Mole, and list every setting that differs.  Empty list = identical."""
    import importlib.util
    os.environ['GRID_LEVEL'] = str(GRID_LEVEL)
    os.environ['PCM_LEBEDEV_ORDER'] = str(PCM_LEBEDEV_ORDER)
    spec = importlib.util.spec_from_file_location('_common_of_record', ROOT / 'dft' / 'common.py')
    common = importlib.util.module_from_spec(spec)
    saved_threads = lib.num_threads()
    spec.loader.exec_module(common)
    lib.num_threads(saved_threads)
    a = common.make_mf(mol, conv_tol=CONV_TOL_SP, solvent=solvent)
    b = make_mf(mol, solvent=solvent, conv_tol=CONV_TOL_SP)
    diffs = []
    ca, cb = [c.__name__ for c in type(a).__mro__], [c.__name__ for c in type(b).__mro__]
    if ca != cb:
        diffs.append(('class', ca, cb))
    for attr in ('xc', 'conv_tol', 'max_cycle', 'conv_tol_grad', 'init_guess', 'level_shift', 'damp', 'diis_space',
                 'small_rho_cutoff', 'nlc'):
        va, vb = getattr(a, attr, None), getattr(b, attr, None)
        if va != vb:
            diffs.append((attr, va, vb))
    for attr in ('level', 'atom_grid', 'prune', 'radi_method', 'becke_scheme', 'atomic_radii', 'radii_adjust'):
        va, vb = getattr(a.grids, attr, None), getattr(b.grids, attr, None)
        if (va is not vb) and va != vb:
            diffs.append((f'grids.{attr}', str(va), str(vb)))
    for attr in ('auxbasis',):
        va, vb = getattr(a.with_df, attr, None), getattr(b.with_df, attr, None)
        if va != vb:
            diffs.append((f'with_df.{attr}', va, vb))
    if solvent:
        for attr in ('method', 'eps', 'lebedev_order', 'vdw_scale', 'r_probe', 'radii_table',
                     'surface_discretization_method', 'max_cycle', 'conv_tol', 'equilibrium_solvation', 'frozen'):
            va, vb = getattr(a.with_solvent, attr, None), getattr(b.with_solvent, attr, None)
            if va != vb:
                diffs.append((f'with_solvent.{attr}', va, vb))
    return diffs


# ---------------------------------------------------------------- fragment analysis
def _coeff_and_overlap(mol, S=None):
    return S if S is not None else mol.intor_symmetric('int1e_ovlp')


def ao_indices(mol, atoms):
    sl = mol.aoslice_by_atom()
    return np.concatenate([np.arange(sl[a][2], sl[a][3]) for a in atoms]) if atoms else np.array([], int)


def mo_fragment_share(mol, c, atoms, method='mulliken', S=None, S_half=None):
    """Fraction (0..1) of one MO (coefficient vector c) on a set of atoms.
    'mulliken': sum_{mu in frag} c_mu (S c)_mu  (as dft/fmo.py);  'lowdin': sum_{mu in frag} ((S^1/2 c)_mu)^2."""
    S = _coeff_and_overlap(mol, S)
    idx = ao_indices(mol, atoms)
    if method == 'mulliken':
        return float((c[idx] * (S @ c)[idx]).sum())
    if method == 'lowdin':
        if S_half is None:
            S_half = lowdin_half(S)
        v = S_half @ c
        return float((v[idx] ** 2).sum())
    raise ValueError(method)


def lowdin_half(S):
    w, U = np.linalg.eigh(S)
    return (U * np.sqrt(w)) @ U.T


def fragment_composition(mf, frags, orbitals=None, methods=('mulliken', 'lowdin')):
    """{orbital label: {method: {fragment: percent}}} for the frontier orbitals of an RKS/UKS solution.
    orbitals: dict label -> (spin index or None, MO index); default HOMO-1..LUMO+1 (alpha for UKS)."""
    mol = mf.mol
    S = mol.intor_symmetric('int1e_ovlp'); S_half = lowdin_half(S)
    C, occ = mf.mo_coeff, mf.mo_occ
    unrestricted = np.ndim(occ) == 2
    if orbitals is None:
        if unrestricted:
            orbitals = {}
            for s, tag in ((0, 'a'), (1, 'b')):
                n = int((occ[s] > 0.5).sum())
                orbitals.update({f'HOMO({tag})': (s, n - 1), f'LUMO({tag})': (s, n)})
        else:
            n = int((occ > 0.5).sum())
            orbitals = {'HOMO-1': (None, n - 2), 'HOMO': (None, n - 1), 'LUMO': (None, n), 'LUMO+1': (None, n + 1)}
    out = {}
    for label, (s, k) in orbitals.items():
        c = C[s][:, k] if s is not None else C[:, k]
        out[label] = {m: {f: round(100 * mo_fragment_share(mol, c, atoms, m, S, S_half), 2)
                          for f, atoms in frags.items()} for m in methods}
    return out


def fragment_populations(mf, frags):
    """Mulliken and Loewdin fragment charges (and spin populations for UKS)."""
    mol = mf.mol
    S = mol.intor_symmetric('int1e_ovlp'); S_half = lowdin_half(S)
    dm = mf.make_rdm1()
    unrestricted = np.ndim(dm) == 3
    dm_tot = dm[0] + dm[1] if unrestricted else dm
    Z = mol.atom_charges()
    res = {}
    for name, P in (('mulliken', None), ('lowdin', None)):
        if name == 'mulliken':
            pop_ao = np.einsum('ij,ji->i', dm_tot, S)
            spin_ao = np.einsum('ij,ji->i', dm[0] - dm[1], S) if unrestricted else None
        else:
            pop_ao = np.diag(S_half @ dm_tot @ S_half)
            spin_ao = np.diag(S_half @ (dm[0] - dm[1]) @ S_half) if unrestricted else None
        r = {}
        for f, atoms in frags.items():
            idx = ao_indices(mol, atoms)
            r[f] = {'charge': float(Z[atoms].sum() - pop_ao[idx].sum())}
            if unrestricted:
                r[f]['spin'] = float(spin_ao[idx].sum())
        res[name] = r
    return res


# ---------------------------------------------------------------- versions
def software_versions():
    v = {'python': sys.version.split()[0], 'pyscf': pyscf.__version__, 'numpy': np.__version__,
         'executable': sys.executable, 'prefix': sys.prefix}
    for mod in ('scipy', 'h5py', 'geometric'):
        try:
            v[mod] = __import__(mod).__version__
        except Exception:
            v[mod] = None
    try:
        from pyscf.dft import libxc
        v['libxc'] = libxc.libxc_version()
    except Exception:
        v['libxc'] = None
    try:
        import importlib.metadata as md
        v['pyscf-dispersion'] = md.version('pyscf-dispersion')
    except Exception:
        v['pyscf-dispersion'] = None
    # conda build strings of the packages that matter (identifies e.g. the CPU-level build of PySCF)
    meta = pathlib.Path(sys.prefix) / 'conda-meta'
    if meta.is_dir():
        builds = {}
        for name in ('pyscf', 'pyscf-dispersion', 'libxc-c', 'libcint', 'libopenblas', 'libblas', 'numpy',
                     'scipy', 'geometric', 'python', '_x86_64-microarch-level', 'simple-dftd3', 'dftd4'):
            for p in meta.glob(f'{name}-[0-9]*.json'):
                d = json.load(open(p))
                if d.get('name') == name:
                    builds[name] = f"{d['version']} {d['build']}"
        v['conda_builds'] = builds
    return v


def host_info():
    cpu = None
    try:
        for line in open('/proc/cpuinfo'):
            if line.startswith('model name'):
                cpu = line.split(':', 1)[1].strip(); break
    except OSError:
        pass
    if cpu is None:                       # ARM: no 'model name' line; use the board model (e.g. Raspberry Pi 5)
        try:
            cpu = open('/proc/device-tree/model').read().strip('\x00 \n')
        except OSError:
            cpu = platform.processor() or None
    return {'hostname': socket.gethostname(), 'machine': platform.machine(), 'cpu': cpu,
            'cpu_count': os.cpu_count(), 'ram_mb': _ram_mb()}


def peak_rss_mb():
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


# ---------------------------------------------------------------- job runner: args, log, result.json
def add_common_args(p):
    p.add_argument('--job', required=True, help='job id (recorded in result.json)')
    p.add_argument('--out', required=True, help='output directory (result.json, log.txt, chk/)')
    p.add_argument('--smoke', action='store_true',
                   help='reduced settings for a quick end-to-end test (never for production numbers)')
    p.add_argument('--keep-scratch', action='store_true', help='do not delete <out>/scratch after success')
    return p


def _json_default(o):
    if hasattr(o, 'tolist'):
        return o.tolist()
    if isinstance(o, (set, tuple)):
        return list(o)
    if isinstance(o, pathlib.Path):
        return str(o)
    return str(o)


def save_json(path, obj):
    path = pathlib.Path(path)
    tmp = path.with_name(path.name + '.tmp')
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=2, default=_json_default)
    os.replace(tmp, path)


class Job:
    """One calculation job: owns <out>/result.json and <out>/log.txt.

    result.json layout:
      job, script, argv, smoke, status (running|done|failed), started/finished (IST), host, threads,
      max_memory_mb, scratch, software, script_sha256, settings, inputs, items {key: {...}}, timings {label: s},
      error (on failure).
    `items` holds the finished units of work (one SCF, one scan point, ...); scripts call job.done(key) to skip
    them on restart and job.put(key, record) to store them as they finish."""

    def __init__(self, args, script, settings=None):
        self.args = args
        self.out = pathlib.Path(args.out).resolve()
        self.out.mkdir(parents=True, exist_ok=True)
        self.chk_dir = self.out / 'chk'
        self.chk_dir.mkdir(exist_ok=True)
        self.path = self.out / 'result.json'
        self._logf = open(self.out / 'log.txt', 'a', buffering=1)
        prev = {}
        if self.path.exists():
            try:
                prev = json.load(open(self.path))
            except Exception:
                prev = {}
        if prev and bool(prev.get('smoke')) != bool(args.smoke):
            raise SystemExit(f'{self.path} holds a {"smoke" if prev.get("smoke") else "production"} run; '
                             f'use another --out for a {"smoke" if args.smoke else "production"} run')
        self.scratch, self._own_scratch = setup_scratch(self.out)
        here = pathlib.Path(script).resolve()
        self.r = {
            'job': args.job, 'script': str(here.relative_to(ROOT)) if str(here).startswith(str(ROOT)) else str(here),
            'argv': sys.argv, 'smoke': bool(args.smoke), 'status': 'running',
            'started': prev.get('started') or now_ist(), 'restarts': prev.get('restarts', -1) + 1,
            'last_start': now_ist(), 'finished': None,
            'host': host_info(), 'threads': THREADS, 'max_memory_mb': MAX_MEMORY_MB, 'scratch': str(self.scratch),
            'software': software_versions(),
            'script_sha256': {here.name: sha256_file(here),
                              'rev_common.py': sha256_file(pathlib.Path(__file__).resolve())},
            'settings': settings or {}, 'inputs': prev.get('inputs', {}), 'items': prev.get('items', {}),
            'timings': prev.get('timings', {}), 'error': None,
        }
        self.save()
        self.log(f"job {args.job} ({self.r['script']}) on {self.r['host']['hostname']}: threads {THREADS}, "
                 f"max_memory {MAX_MEMORY_MB} MB, scratch {self.scratch}, smoke {bool(args.smoke)}, "
                 f"restart #{self.r['restarts']}; PySCF {pyscf.__version__}")

    # -- logging
    def log(self, msg):
        line = f"[{datetime.datetime.now(IST).strftime('%Y-%m-%d %I:%M:%S %p IST')}] {msg}"
        print(line, flush=True)
        self._logf.write(line + '\n')

    # -- result bookkeeping
    def save(self):
        save_json(self.path, self.r)

    def done(self, key):
        it = self.r['items'].get(key)
        return bool(it) and it.get('status') == 'done'

    def get(self, key):
        return self.r['items'].get(key)

    def put(self, key, record, status='done'):
        record = dict(record); record['status'] = status; record['at'] = now_ist()
        self.r['items'][key] = record
        self.save()

    def input_file(self, label, path):
        p = pathlib.Path(path)
        self.r['inputs'][label] = {'path': str(p.relative_to(ROOT)) if str(p.resolve()).startswith(str(ROOT)) else str(p),
                                   'sha256': sha256_file(p)}
        self.save()

    @contextlib.contextmanager
    def timer(self, label):
        t = time.time()
        try:
            yield
        finally:
            self.r['timings'][label] = round(time.time() - t, 2)
            self.save()

    # -- end of run
    def finish(self, ok=True, error=None):
        self.r['status'] = 'done' if ok else 'failed'
        self.r['finished'] = now_ist()
        self.r['peak_rss_mb'] = round(peak_rss_mb(), 1)
        if error:
            self.r['error'] = error
        self.save()
        self.log(f"job {self.r['job']} {'finished' if ok else 'FAILED'} (peak RSS {self.r['peak_rss_mb']:.0f} MB)")
        if ok and self._own_scratch and not getattr(self.args, 'keep_scratch', False):
            shutil.rmtree(self.scratch, ignore_errors=True)


def main_wrapper(build_parser, run, script):
    """Parse args, create the Job, call run(job, args); record failures and exit non-zero on any error."""
    args = build_parser().parse_args()
    job = Job(args, script)
    try:
        ok = run(job, args)
        job.finish(ok is not False)
        sys.exit(0 if ok is not False else 1)
    except SystemExit:
        raise
    except BaseException as e:            # KeyboardInterrupt and SIGTERM-as-exception included
        tb = traceback.format_exc()
        job.log(f'ERROR {type(e).__name__}: {e}\n{tb}')
        job.finish(False, {'type': type(e).__name__, 'message': str(e), 'traceback': tb})
        sys.exit(1)


def install_sigterm_handler():
    """Turn SIGTERM (queue stop, kill) into an exception so result.json records the failure."""
    import signal

    def _h(signum, frame):
        raise KeyboardInterrupt(f'signal {signum}')
    signal.signal(signal.SIGTERM, _h)
