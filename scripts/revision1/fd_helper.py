"""
Extra workers for a running finite-difference Hessian of c0_case.py (--freq fd), added 22 Nov 2026.

c0_case.py computes the FD Hessian in chunks of Cartesian columns (cols_XXX_YYY.npz in <out>/fd/), one chunk after
another in a single process, and skips every chunk whose .npz already exists.  This helper computes chunks from the END
of the list in a separate process, with the same command (dft/_fd_worker.py <case> c0 c1 <part> <step>) and the same
environment of record (threads, memory, GRID_LEVEL, PCM_LEBEDEV_ORDER, NO_AUTO_PCM), copied from the running main worker
when there is one and otherwise set to the values below.  Each result is written to a private <name>.<tag>.part.npz and
renamed atomically to <name>.npz, as c0_case.py does, so c0_case.py treats it as done and in the end only assembles.

A helper takes chunk j only if: no .npz exists for it; no other helper has claimed it (<name>.claim, created with
O_EXCL); and j >= m + --gap, where m is the lowest chunk that is neither finished nor claimed (the one c0_case.py is
computing, since it works upward from 0).  With the gap the main process cannot reach a chunk a helper is still
computing (it has at least gap - 1 whole chunks to do first).  It stops when no chunk qualifies.

Usage (one queue job per extra worker; THREADS from the queue):
  python dft/revision/fd_helper.py --fd-dir <out>/fd --case btpcz_pcm --tag fdh1 [--chunk 3] [--gap 3] [--step 0.005]
"""
import os, sys, time, argparse, subprocess, pathlib, datetime
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'dft'))
KEYS = ('MAX_MEMORY_MB', 'DFT_MAX_MEMORY_MB', 'GRID_LEVEL', 'PCM_LEBEDEV_ORDER', 'NO_AUTO_PCM')


def now():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=5, minutes=30))).strftime('%Y-%m-%d %I:%M:%S %p IST')


def log(msg):
    print(f'[{now()}] {msg}', flush=True)


def main_worker_env(case, fd):
    """Environment of the main c0_case.py worker (the _fd_worker.py process writing <name>.part.npz in fd)."""
    for p in pathlib.Path('/proc').iterdir():
        if not p.name.isdigit():
            continue
        try:
            argv = (p / 'cmdline').read_bytes().split(b'\0')
            argv = [a.decode() for a in argv if a]
            if len(argv) >= 6 and argv[1].endswith('dft/_fd_worker.py') and argv[2] == case \
                    and argv[5].startswith(str(fd) + '/') and pathlib.Path(argv[5]).name.count('.') == 2:   # cols_a_b.part.npz
                env = dict(l.decode().split('=', 1) for l in (p / 'environ').read_bytes().split(b'\0') if b'=' in l)
                return {k: env[k] for k in KEYS if k in env}, int(p.name)
        except (OSError, PermissionError, ValueError):
            continue
    return None, None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--fd-dir', required=True)
    ap.add_argument('--case', required=True)
    ap.add_argument('--tag', required=True, help='unique name of this helper (the queue job id)')
    ap.add_argument('--chunk', type=int, default=3, help='columns per chunk; must equal c0_case.py --fd-chunk')
    ap.add_argument('--gap', type=int, default=3)
    ap.add_argument('--step', type=float, default=0.005, help='Bohr; must equal c0_case.py --fd-step')
    ap.add_argument('--grid', default='4'); ap.add_argument('--lebedev', default='29')
    a = ap.parse_args()
    fd = pathlib.Path(a.fd_dir).resolve()
    from common import load_case
    case = load_case(a.case)
    natm = int(open(case['dir'] / 'opt.xyz').readline().split()[0])
    ncol = 3 * natm
    edges = list(range(0, ncol, a.chunk)) + [ncol]
    chunks = list(zip(edges[:-1], edges[1:]))
    t = os.environ.get('THREADS', '4')
    mem = os.environ.get('MAX_MEMORY_MB', '8000')
    rec = {'MAX_MEMORY_MB': mem, 'DFT_MAX_MEMORY_MB': mem, 'GRID_LEVEL': a.grid, 'PCM_LEBEDEV_ORDER': a.lebedev,
           'NO_AUTO_PCM': '1'}
    menv, mpid = main_worker_env(a.case, fd)
    if menv:
        if menv != rec:
            log(f'environment of the main worker (pid {mpid}) differs from the defaults: {menv} vs {rec}; using the main worker\'s')
        rec = menv
    else:
        log('no main worker running; using the default environment of record')
    log(f'{a.tag}: case {a.case}, {natm} atoms, {len(chunks)} chunks of {a.chunk} columns, step {a.step} Bohr, '
        f'{t} threads, env {rec}')

    def done(i):
        c0, c1 = chunks[i]
        f = fd / f'cols_{c0:03d}_{c1:03d}.npz'
        if not f.exists():
            return False
        try:
            return sorted(int(x) for x in np.load(f).files) == list(range(c0, c1))
        except Exception:
            return False

    def claimed(i):
        c0, c1 = chunks[i]
        return (fd / f'cols_{c0:03d}_{c1:03d}.claim').exists()

    n_done = 0
    while True:
        missing = [i for i in range(len(chunks)) if not done(i)]
        free = [i for i in missing if not claimed(i)]
        if not free:
            break
        m = min(free)
        cand = [i for i in free if i >= m + a.gap]
        if not cand:
            log(f'{a.tag}: no chunk at least {a.gap} ahead of chunk {m} (main); stopping')
            break
        j = max(cand)
        c0, c1 = chunks[j]
        name = f'cols_{c0:03d}_{c1:03d}'
        try:
            os.close(os.open(fd / f'{name}.claim', os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        (fd / f'{name}.claim').write_text(f'{a.tag} {now()}\n')
        part = fd / f'{name}.{a.tag}.part.npz'
        part.unlink(missing_ok=True)
        scr = pathlib.Path(os.environ.get('PYSCF_TMPDIR', '/tmp')) / name
        scr.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ)
        env.update({'THREADS': t, 'DFT_THREADS': t, 'OMP_NUM_THREADS': t, 'MKL_NUM_THREADS': t,
                    'OPENBLAS_NUM_THREADS': t, 'NUMEXPR_NUM_THREADS': t, 'VECLIB_MAXIMUM_THREADS': t,
                    'MWFN_THREADS': t, 'PYSCF_TMPDIR': str(scr), 'TMPDIR': str(scr), **rec})
        log(f'{a.tag}: chunk {name} (index {j}; main at {m}) started')
        t0 = time.time()
        r = subprocess.run([sys.executable, str(ROOT / 'dft' / '_fd_worker.py'), a.case, str(c0), str(c1), str(part),
                            repr(a.step)], cwd=str(ROOT), env=env)
        if r.returncode != 0 or not part.exists():
            (fd / f'{name}.claim').unlink(missing_ok=True)
            log(f'{a.tag}: chunk {name} FAILED (exit {r.returncode})')
            sys.exit(1)
        os.replace(part, fd / f'{name}.npz')
        n_done += 1
        log(f'{a.tag}: chunk {name} done ({2 * (c1 - c0)} gradients, {(time.time() - t0) / 60:.1f} min)')
    log(f'{a.tag}: finished, {n_done} chunk(s) computed')


if __name__ == '__main__':
    main()
