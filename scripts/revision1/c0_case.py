"""
Restartable driver for revision item C0 (the P2_1/n ion pair, cases btpcz_gas and btpcz_pcm): the published pipeline
stages optimize -> proton_check -> frequencies -> fmo, made resumable so that an interruption (reboot, power cut,
`jobqueue.py stop --kill`, a runner restart) costs minutes, not the whole job.

Why not `run_case.sh <case> -- --stage <s>` in a loop (the first C0 command): pipeline.run_stage skips nothing;
dft/optimize.py without --restart re-runs the whole optimisation and truncates opt_traj.xyz; PySCF 2.14's
geometric_solver.optimize() returns the last geometry at maxsteps with no flag, so optimize.py would call a stalled
optimisation "converged"; dft/frequencies.py with one worker writes its .npz only after all 198 gradients, into a fresh
mkdtemp, so an interruption at gradient 190 loses everything.  dft/ itself is left untouched.

Stages (each skipped when its output in molecules/<mol>/results/<case>/ exists and is not older than opt.xyz):
  preopt        B3LYP/6-31G(d) pre-optimisation from the crystal pair with GAU_LOOSE criteria and SCF 1e-8, in the
                case's own model (gas or IEF-PCM), as dft/optimize.py stage 0 -> preopt.xyz.  Driven through
                dft/revision/opt_d3.py (--disp none), which restarts from its own trajectory.  Like optimize.py, the
                pre-optimised geometry is used even if the loose criteria were not met within the cap (it is only a
                start); that is recorded.  Skipped when preopt.xyz exists.
  optimize      B3LYP/6-311++G(d,p) optimisation from preopt.xyz with geomeTRIC's defaults (= Gaussian's) and SCF 1e-9
                gas / 1e-8 PCM, as dft/optimize.py, through opt_d3.py: explicit converged flag, restart from its own
                trajectory, stall stop (--negligible-stop, see opt_d3.py), cumulative step and wall-time caps.  Only an
                accepted optimisation is handed on: opt.xyz, opt_summary.json (optimize.py's keys plus converged /
                converged_on), geometry_vs_start.csv, opt_traj.xyz and opt.chk in the case folder.
  proton_check  dft/pipeline.py --stage proton_check through run_case.sh (published code, unchanged).
  frequencies   --freq fd: the finite-difference Hessian of dft/frequencies.py, computed by calling dft/_fd_worker.py
                (unchanged) on chunks of --fd-chunk Cartesian columns one after another into <out>/fd/ (a chunk is
                kept once complete and skipped on restart; a manifest ties the chunks to the sha256 of opt.xyz and the
                step), after the same reference SCF (conv_tol 1e-10 from opt.chk, writes final.chk).  hessian.npy and
                freq.json are then assembled exactly as frequencies.py does (column layout, symmetrisation,
                thermo.harmonic_analysis, dft/vibrations.py projection, thermo.thermo at 298.15 K / 1 atm).  No
                machine-wide lock is taken (FREQ_NOLOCK is irrelevant here).
                --freq softmode: dft/revision/softmode_check.py at opt.xyz (B3LYP, no dispersion) along the --nmodes
                softest modes of --ref-hessian (for btpcz_pcm: the btpcz_gas FD Hessian, same atom order): 2N+1
                gradients instead of 6N.  Writes softmode_check.json; the pipeline's frequencies stage is recorded
                as skipped (no freq.json, so RESULTS.md has no frequency table for that case).
                --freq none: nothing (status untouched).
  fmo           dft/pipeline.py --stage fmo through run_case.sh (published code, unchanged).
The Multiwfn/plot/render/report stages stay with `run_case.sh <case> -- --from mwfn` (jobs c0_*_analysis).

Settings of record: every subprocess of published code gets the environment run_case.sh exports (GRID_LEVEL=4,
PCM_LEBEDEV_ORDER=29, DFT_THREADS = OMP/MKL/OpenBLAS threads = THREADS, DFT_MAX_MEMORY_MB = MAX_MEMORY_MB,
PYSCF_TMPDIR below the job's scratch); they are recorded in result.json.  run_case.sh needs tools/env (job c0_prep).

Usage:
  python dft/revision/c0_case.py --job ID --out DIR --case btpcz_gas|btpcz_pcm
         [--stages preopt,optimize,proton_check,frequencies,fmo] [--freq fd|softmode|none]
         [--fd-chunk 3] [--fd-step 0.005] [--ref-hessian PATH] [--nmodes 10] [--soft-step 0.02]
         [--preopt-basis 6-31g(d)] [--preopt-maxsteps 300] [--opt-maxsteps 150] [--opt-max-hours H]
         [--negligible-stop 10] [--smoke]
--smoke: reduced settings for the subprocesses (opt_d3 / softmode_check --smoke, GRID_LEVEL 2, Lebedev 17); use it
         only with a throwaway case (e.g. a small molecule with STO-3G basis_overrides), never with btpcz: results go
         into the case folder.
"""
import sys, os, json, time, shutil, pathlib, argparse, subprocess, datetime, importlib.util
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# ---------------------------------------------------------------- internal mode: the reference SCF of frequencies.py
if len(sys.argv) >= 4 and sys.argv[1] == '--internal-refscf':
    # Run in a child process so that its density-fitting tensors are freed before the FD workers start.
    _ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
    sys.path.insert(0, str(_ROOT / 'dft'))
    import time as _t
    from common import load_case as _lc, mol_for_case as _mfc, make_mf as _mm, run_scf as _rs, save_json as _sj
    _case = _lc(sys.argv[2]); _d = _case['dir']
    _mol = _mfc(_case, _d / 'opt.xyz')
    _mf = _mm(_mol, conv_tol=1e-10, chkfile=_d / 'final.chk', solvent=_case.get('solvent'))
    _t0 = _t.time()
    _rs(_mf, init_chk=_d / 'opt.chk')
    _sj(sys.argv[3], {'E_Ha': float(_mf.e_tot), 'wall_s': round(_t.time() - _t0, 1), 'nao': int(_mol.nao),
                      'natm': int(_mol.natm), 'grid_level': int(_mf.grids.level),
                      'lebedev_order': (int(_mf.with_solvent.lebedev_order) if hasattr(_mf, 'with_solvent') else None),
                      'max_memory_mb': _mol.max_memory, 'chk': str(_d / 'final.chk')})
    print(f"[c0 refscf] {_case['name']}: SCF E = {_mf.e_tot:.10f} Ha ({_t.time() - _t0:.0f} s); {_mol.nao} basis "
          f"functions, {_mol.natm} atoms", flush=True)
    sys.exit(0)

import rev_common as rc          # first: pins the thread pools
import numpy as np

STAGES = ('preopt', 'optimize', 'proton_check', 'frequencies', 'fmo')


# ---------------------------------------------------------------- helpers
def published_common():
    """dft/common.py (the code of record) as a module, for bonds_from_geometry / hbond_geometry / mol_for_case."""
    os.environ.setdefault('DFT_THREADS', str(rc.THREADS))
    if str(rc.ROOT / 'dft') not in sys.path:
        sys.path.append(str(rc.ROOT / 'dft'))
    spec = importlib.util.spec_from_file_location('common', rc.ROOT / 'dft' / 'common.py')
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault('common', mod)
    spec.loader.exec_module(mod)
    return mod


def case_dir(case):
    p = rc.find_case_file(case['name'])
    d = p.parent.parent / 'results' / case['name']
    d.mkdir(parents=True, exist_ok=True)
    return d


def stage_note(d, text):
    """Same line format as dft/common.stage_note (local time; IST on the hosts)."""
    with open(d / 'STAGES.md', 'a') as f:
        f.write(f"- {datetime.datetime.now():%Y-%m-%d %H:%M}  {text}\n")


def set_stage(d, case_name, stage, rec):
    """Record a stage in status.json the way dft/pipeline.py run_stage does."""
    p = d / 'status.json'
    st = json.load(open(p)) if p.exists() else {'case': case_name, 'stages': {}}
    st.setdefault('stages', {})[stage] = rec
    rc.save_json(p, st)


def stage_status(d, stage):
    p = d / 'status.json'
    if not p.exists():
        return None
    return (json.load(open(p)).get('stages', {}).get(stage) or {}).get('status')


def fresh(d, names):
    """All files exist and none is older than opt.xyz."""
    ref = d / 'opt.xyz'
    t0 = ref.stat().st_mtime if ref.exists() else 0
    return all((d / n).exists() and (d / n).stat().st_mtime >= t0 - 1e-3 for n in names)


def now_local():
    return datetime.datetime.now().isoformat(timespec='seconds')


class Runner:
    def __init__(self, job, args, case):
        self.job, self.args, self.case = job, args, case
        self.d = case_dir(case)
        self.prefix, self.model = case['name'].rsplit('_', 1)
        self.grid = 2 if args.smoke else rc.GRID_LEVEL
        self.leb = 17 if args.smoke else rc.PCM_LEBEDEV_ORDER
        (job.out / 'logs').mkdir(exist_ok=True)

    # -- environment of record for published code (what run_case.sh exports) and for revision subprocesses
    def env(self, tag):
        t = str(rc.THREADS)
        scr = self.job.scratch / tag
        scr.mkdir(parents=True, exist_ok=True)
        e = dict(os.environ)
        e.update({'THREADS': t, 'DFT_THREADS': t, 'OMP_NUM_THREADS': t, 'MKL_NUM_THREADS': t,
                  'OPENBLAS_NUM_THREADS': t, 'NUMEXPR_NUM_THREADS': t, 'VECLIB_MAXIMUM_THREADS': t,
                  'MAX_MEMORY_MB': str(rc.MAX_MEMORY_MB), 'DFT_MAX_MEMORY_MB': str(rc.MAX_MEMORY_MB),
                  'GRID_LEVEL': str(self.grid), 'PCM_LEBEDEV_ORDER': str(self.leb),
                  'PYSCF_TMPDIR': str(scr), 'TMPDIR': str(scr), 'NO_AUTO_PCM': '1', 'MWFN_THREADS': t})
        return e, scr

    def run(self, cmd, tag, log_name=None):
        """Run a subprocess with the environment of record; output appended to <out>/logs/<tag>.log."""
        env, scr = self.env(tag)
        lf = self.job.out / 'logs' / f'{log_name or tag}.log'
        with open(lf, 'a') as f:
            f.write(f"\n===== {now_local()}  $ {' '.join(map(str, cmd))}\n"); f.flush()
            t = time.time()
            rc_ = subprocess.run([str(c) for c in cmd], cwd=str(rc.ROOT), env=env, stdout=f,
                                 stderr=subprocess.STDOUT).returncode
        wall = time.time() - t
        if rc_ == 0:
            shutil.rmtree(scr, ignore_errors=True)
        return rc_, wall, lf

    def opt_d3(self, tag, extra):
        cmd = [sys.executable, rc.ROOT / 'dft' / 'revision' / 'opt_d3.py', '--job', f'{self.job.r["job"]}_{tag}',
               '--out', self.job.out / tag, '--case', self.prefix, '--model', self.model, '--disp', 'none',
               '--seed-hessian', 'none', '--no-fmo'] + extra + (['--smoke'] if self.args.smoke else [])
        rc_, wall, lf = self.run(cmd, tag)
        res = self.job.out / tag / 'result.json'
        r = json.load(open(res)) if res.exists() else {}
        return rc_, wall, r, lf

    # ---------------------------------------------------------------- stages
    def preopt(self):
        d, a = self.d, self.args
        if (d / 'preopt.xyz').exists():
            self.job.log('preopt: preopt.xyz exists; skipped (as dft/optimize.py reuses it)')
            return True
        extra = ['--start', 'xray', '--criteria', 'loose', '--scf-tol', '1e-8', '--maxsteps', str(a.preopt_maxsteps),
                 '--negligible-stop', str(a.negligible_stop)]
        if not a.smoke:
            extra += ['--basis', a.preopt_basis]
        if a.preopt_max_hours:
            extra += ['--max-hours', str(a.preopt_max_hours)]
        self.job.log(f"preopt: B3LYP/{'sto-3g (smoke)' if a.smoke else a.preopt_basis}, GAU_LOOSE, from the crystal pair")
        rc_, wall, r, lf = self.opt_d3('preopt', extra)
        o = (r.get('items') or {}).get('opt') or {}
        usable = (self.job.out / 'preopt' / 'opt.xyz').exists() and (
            rc_ == 0 or (o.get('status') == 'stopped' and o.get('stop_reason') in ('maxsteps', 'max_hours')))
        if not usable:
            raise RuntimeError(f'pre-optimisation failed (exit {rc_}); see {lf} and {self.job.out / "preopt"}/log.txt')
        syms, X = rc.read_xyz(self.job.out / 'preopt' / 'opt.xyz')
        how = 'converged (GAU_LOOSE)' if o.get('converged') else (o.get('converged_on') or
                                                                   f"NOT converged ({o.get('stop_reason')}); used as start")
        rc.write_xyz(d / 'preopt.xyz', syms, X,
                     f"pre-optimised at B3LYP/{'sto-3g' if a.smoke else a.preopt_basis} ({how}), "
                     f"E = {o.get('E_final_Ha', float('nan')):.10f} Ha; dft/revision/c0_case.py via opt_d3.py")
        msg = (f"pre-optimisation at B3LYP/{'sto-3g' if a.smoke else a.preopt_basis}: {how} after {o.get('steps')} "
               f"steps, {o.get('wall_min_total', wall / 60):.1f} min (dft/revision/c0_case.py)")
        stage_note(d, msg)
        self.job.put('preopt', {'steps': o.get('steps'), 'converged': o.get('converged'), 'how': how,
                                'E_Ha': o.get('E_final_Ha'), 'wall_s': round(wall, 1)})
        self.job.log('preopt: ' + msg)
        return True

    def optimize(self):
        d, a, case = self.d, self.args, self.case
        if (d / 'opt.xyz').exists() and (d / 'opt_summary.json').exists():
            self.job.log('optimize: opt.xyz and opt_summary.json exist; skipped')
            return True
        start = d / 'preopt.xyz' if (d / 'preopt.xyz').exists() else rc.ROOT / case['xyz']
        extra = ['--start', str(start), '--maxsteps', str(a.opt_maxsteps), '--negligible-stop', str(a.negligible_stop)]
        if a.opt_max_hours:
            extra += ['--max-hours', str(a.opt_max_hours)]
        t0 = now_local()
        set_stage(d, case['name'], 'optimize', {'status': 'running', 'start': t0})
        stage_note(d, f"optimisation started from {start.name} (dft/revision/c0_case.py via opt_d3.py); "
                      f"solvent={case.get('solvent') or 'none (gas phase)'}")
        self.job.log(f'optimize: from {start.name}, cap {a.opt_maxsteps} evaluations'
                     f"{f', {a.opt_max_hours} h' if a.opt_max_hours else ''}, stall stop {a.negligible_stop}")
        rc_, wall, r, lf = self.opt_d3('opt', extra)
        o = (r.get('items') or {}).get('opt') or {}
        if rc_ != 0 or not (o.get('accepted') or a.smoke):
            set_stage(d, case['name'], 'optimize', {'status': 'failed', 'start': t0, 'end': now_local(),
                                                    'wall_min': round(wall / 60, 2),
                                                    'note': f"exit {rc_}: {o.get('stop_reason') or 'see ' + str(lf)}"})
            stage_note(d, f"optimize FAILED / stopped ({o.get('stop_reason')}, {o.get('steps')} evaluations); "
                          f"rerun the job to continue from {self.job.out / 'opt' / 'opt_traj.xyz'}")
            raise RuntimeError(f"optimisation not accepted (exit {rc_}, stop reason {o.get('stop_reason')}); "
                               f"rerun to continue; see {lf}")
        self.hand_over(o, wall, t0)
        return True

    def hand_over(self, o, wall, t0):
        """opt_d3 output -> the files dft/optimize.py would have written in the case folder."""
        d, case, a = self.d, self.case, self.args
        cm = published_common()
        od = self.job.out / 'opt'
        syms, xyz = rc.read_xyz(od / 'opt.xyz')
        E = float(o['E_final_Ha'])
        how = ('optimised (all five geomeTRIC/Gaussian criteria)' if o.get('converged')
               else f"accepted on {o.get('converged_on')}")
        for src, dst in ((od / 'opt_traj.xyz', d / 'opt_traj.xyz'), (od / 'chk' / 'opt.chk', d / 'opt.chk')):
            if src.exists():
                shutil.copy2(src, str(dst) + '.tmp'); os.replace(str(dst) + '.tmp', dst)
        s0, x0 = cm.read_xyz(rc.ROOT / case['xyz'])
        lab = case.get('cif_labels') or cm.atom_labels(syms)
        rows = []
        for i, j, d0 in cm.bonds_from_geometry(s0, x0):
            rows.append((f'{lab[i]}-{lab[j]}', d0, float(np.linalg.norm(xyz[i] - xyz[j]))))
        with open(d / 'geometry_vs_start.csv', 'w') as f:
            f.write('bond,start_A,dft_A,diff_A\n')
            for n, x, y in rows:
                f.write(f'{n},{x:.4f},{y:.4f},{y - x:+.4f}\n')
        diffs = np.array([y - x for _, x, y in rows])
        pre = self.job.get('preopt') or {}
        hist = [{'step': h['step'], 'E_Ha': h['E_Ha'], 'gmax': h.get('gmax_cart', h['gmax']),
                 'grms': h.get('grms_cart', h['grms'])} for h in o.get('history', [])]
        summary = {'case': case['name'], 'E_final_Ha': E, 'steps': o.get('steps'),
                   'wall_min': o.get('wall_min_total', round(wall / 60, 2)),
                   'preopt_basis': ('sto-3g' if a.smoke else a.preopt_basis) if (d / 'preopt.xyz').exists() else None,
                   'preopt_steps': pre.get('steps'), 'n_basis': o.get('nao'),
                   'bond_rmsd_A': float(np.sqrt((diffs ** 2).mean())), 'bond_maxdev_A': float(abs(diffs).max()),
                   'history': hist,
                   # revision keys (dft/optimize.py has none of these; finalize_opt.py writes converged_on)
                   'converged': bool(o.get('converged')), 'converged_on': o.get('converged_on'),
                   'optimiser': 'geomeTRIC 1.1.1 via dft/revision/opt_d3.py (driven by dft/revision/c0_case.py)',
                   'history_gradient_norms': 'Cartesian components (as dft/optimize.py); opt_d3 result.json also '
                                             'has geomeTRIC per-atom norms',
                   'grid_level': self.grid, 'pcm_lebedev_order': self.leb if case.get('solvent') else None,
                   'opt_result': str((od / 'result.json').relative_to(rc.ROOT)) if str(od).startswith(str(rc.ROOT))
                   else str(od / 'result.json')}
        for name, (i, j) in case.get('distances', {}).items():
            summary.setdefault('distances_A', {})[name] = {'start': float(np.linalg.norm(x0[i] - x0[j])),
                                                         'dft': float(np.linalg.norm(xyz[i] - xyz[j]))}
        for name, (D, H, A) in case.get('hbonds', {}).items():
            summary.setdefault('hbonds', {})[name] = {'start': cm.hbond_geometry(x0, D, H, A),
                                                      'dft': cm.hbond_geometry(xyz, D, H, A)}
        cm.save_json(d / 'opt_summary.json', summary)
        # opt.xyz last: its presence (with opt_summary.json) marks the stage done
        rc.write_xyz(d / 'opt.xyz', syms, xyz, f"{case['title']} -- B3LYP/{cm.basis_description(case)} {how}, "
                                                 f"E = {E:.8f} Ha")
        note = (f"optimisation {'converged' if o.get('converged') else 'accepted on negligible forces'} in "
                f"{o.get('steps')} steps, {summary['wall_min']:.1f} min, E = {E:.8f} Ha"
                + ('' if o.get('converged') else f" ({o.get('converged_on')}); minimum to be confirmed by frequencies"))
        stage_note(d, note)
        set_stage(d, case['name'], 'optimize', {'status': 'done', 'start': t0, 'end': now_local(),
                                                'wall_min': round(wall / 60, 2), 'note': note})
        self.job.put('optimize', {'E_Ha': E, 'steps': o.get('steps'), 'converged': bool(o.get('converged')),
                                  'converged_on': o.get('converged_on'),
                                  'opt_xyz_sha256': rc.sha256_file(d / 'opt.xyz')})
        self.job.log('optimize: ' + note)

    def pipeline_stage(self, stage, outputs):
        d, case = self.d, self.case
        if fresh(d, outputs) and stage_status(d, stage) == 'done':
            self.job.log(f'{stage}: {", ".join(outputs)} up to date; skipped')
            return True
        if not (rc.ROOT / 'tools' / 'env' / 'bin' / 'python').exists():
            raise RuntimeError('tools/env/bin/python missing: run the c0_prep job (links tools/env to the environment)')
        self.job.log(f'{stage}: ./run_case.sh {case["name"]} {rc.THREADS} -- --stage {stage}')
        rc_, wall, lf = self.run(['bash', rc.ROOT / 'run_case.sh', case['name'], str(rc.THREADS), '--', '--stage', stage],
                                 stage)
        if rc_ != 0 or not fresh(d, outputs):
            raise RuntimeError(f'{stage} failed (exit {rc_}); see {lf} and {d}/logs/{stage}.log')
        self.job.put(stage, {'wall_s': round(wall, 1), 'outputs': outputs,
                             'opt_xyz_sha256': rc.sha256_file(d / 'opt.xyz')})
        self.job.log(f'{stage}: done in {wall / 60:.1f} min')
        return True

    # -- frequencies
    def frequencies(self):
        a = self.args
        if a.freq == 'none':
            self.job.log('frequencies: --freq none; skipped')
            return True
        if a.freq == 'softmode':
            return self.softmode()
        return self.fd_hessian()

    def softmode(self):
        d, a, case = self.d, self.args, self.case
        if fresh(d, ['softmode_check.json']):
            self.job.log('frequencies (soft-mode check): softmode_check.json up to date; skipped')
            return True
        if not a.ref_hessian:
            raise SystemExit('--freq softmode needs --ref-hessian (for btpcz_pcm: the btpcz_gas hessian.npy)')
        ref = pathlib.Path(a.ref_hessian)
        ref = ref if ref.is_absolute() else rc.ROOT / ref
        if not ref.exists():
            raise RuntimeError(f'reference Hessian {ref} not found (has the gas-phase C0 job finished and been copied?)')
        t0 = now_local()
        cmd = [sys.executable, rc.ROOT / 'dft' / 'revision' / 'softmode_check.py', '--job', f'{self.job.r["job"]}_softmode',
               '--out', self.job.out / 'softmode', '--geom', d / 'opt.xyz', '--model', self.model, '--case', self.prefix,
               '--disp', 'none', '--ref-hessian', ref, '--step', str(a.soft_step)]
        if a.nmodes:
            cmd += ['--nmodes', str(a.nmodes)]
        if a.smoke:
            cmd += ['--smoke']
        self.job.log(f'frequencies: soft-mode check at opt.xyz along the softest modes of {ref}')
        rc_, wall, lf = self.run(cmd, 'softmode')
        res = self.job.out / 'softmode' / 'result.json'
        if rc_ != 0 or not res.exists():
            raise RuntimeError(f'soft-mode check failed (exit {rc_}); see {lf}')
        r = json.load(open(res))
        s = r.get('summary') or {}
        out = {'case': case['name'], 'method': 'soft-mode check (dft/revision/softmode_check.py) in place of the '
                                               'full finite-difference Hessian',
               'reference_hessian': str(ref), 'reference': r.get('reference'), 'settings': r.get('settings'),
               'summary': s, 'residual_gradient': r.get('residual_gradient'),
               'n_gradients': sum(1 for k in (r.get('items') or {}) if k == 'x0' or k.startswith('mode')),
               'result_json': str(res), 'opt_xyz_sha256': rc.sha256_file(d / 'opt.xyz')}
        rc.save_json(d / 'softmode_check.json', out)
        note = (f"frequencies replaced by a soft-mode check ({len(s.get('projected_freqs_cm-1', []))} modes of "
                f"{ref.name}): projected {np.round(s.get('projected_freqs_cm-1', []), 1).tolist()} cm-1; "
                f"{s.get('verdict')}")
        stage_note(d, note)
        set_stage(d, case['name'], 'frequencies', {'status': 'skipped', 'start': t0, 'end': now_local(),
                                                   'wall_min': round(wall / 60, 2), 'note': note + ' (softmode_check.json)'})
        self.job.put('frequencies', {'method': 'softmode', 'wall_s': round(wall, 1), 'verdict': s.get('verdict')})
        self.job.log('frequencies: ' + note)
        return True

    def fd_hessian(self):
        d, a, case = self.d, self.args, self.case
        if fresh(d, ['freq.json', 'hessian.npy']) and stage_status(d, 'frequencies') == 'done':
            self.job.log('frequencies: freq.json and hessian.npy up to date; skipped')
            return True
        fd = self.job.out / 'fd'
        sha = rc.sha256_file(d / 'opt.xyz')
        syms, _ = rc.read_xyz(d / 'opt.xyz')
        natm = len(syms); ncol = 3 * natm
        manifest = {'case': case['name'], 'opt_xyz_sha256': sha, 'step_bohr': a.fd_step, 'ncol': ncol,
                    'grid_level': self.grid, 'pcm_lebedev_order': self.leb}
        mpath = fd / 'manifest.json'
        if mpath.exists() and json.load(open(mpath)) != manifest:
            stale = fd.with_name(f'fd.stale-{datetime.datetime.now():%Y%m%d-%H%M%S}')
            fd.rename(stale)
            self.job.log(f'frequencies: chunks in {fd} were for another opt.xyz / step / grid; moved to {stale.name}')
        fd.mkdir(exist_ok=True)
        rc.save_json(mpath, manifest)
        t_stage = now_local()
        if stage_status(d, 'frequencies') != 'running':
            set_stage(d, case['name'], 'frequencies', {'status': 'running', 'start': t_stage})
            stage_note(d, f'Hessian by finite differences in sequential chunks of {a.fd_chunk} columns '
                          f'({2 * ncol} gradients, step {a.fd_step} Bohr, {rc.THREADS} threads; dft/revision/c0_case.py)')
        # reference SCF (frequencies.py: conv_tol 1e-10 from opt.chk, chkfile final.chk)
        ref = fd / 'ref.json'
        if not (ref.exists() and (d / 'final.chk').exists()
                and (d / 'final.chk').stat().st_mtime >= (d / 'opt.xyz').stat().st_mtime):
            self.job.log('frequencies: reference SCF (conv_tol 1e-10 from opt.chk -> final.chk)')
            rc_, wall, lf = self.run([sys.executable, pathlib.Path(__file__).resolve(), '--internal-refscf',
                                      case['name'], ref], 'fd_refscf')
            if rc_ != 0:
                raise RuntimeError(f'reference SCF failed (exit {rc_}); see {lf}')
        refr = json.load(open(ref))
        self.job.log(f"frequencies: reference E = {refr['E_Ha']:.10f} Ha ({refr['wall_s']:.0f} s, grid "
                     f"{refr['grid_level']}, Lebedev {refr['lebedev_order']})")
        # chunks
        k = max(1, a.fd_chunk)
        edges = list(range(0, ncol, k)) + [ncol]
        prev = self.job.get('fd') or {}
        times = dict(prev.get('chunk_wall_s', {}))
        for c0, c1 in zip(edges[:-1], edges[1:]):
            name = f'cols_{c0:03d}_{c1:03d}'
            fin = fd / f'{name}.npz'
            if fin.exists():
                try:
                    z = np.load(fin)
                    if sorted(int(x) for x in z.files) == list(range(c0, c1)):
                        continue
                except Exception:
                    pass
                fin.unlink()
            part = fd / f'{name}.part.npz'
            part.unlink(missing_ok=True)
            t = time.time()
            rc_, wall, lf = self.run([sys.executable, rc.ROOT / 'dft' / '_fd_worker.py', case['name'], c0, c1, part,
                                      repr(a.fd_step)], 'fd', log_name='fd_worker')
            if rc_ != 0 or not part.exists():
                raise RuntimeError(f'FD chunk {name} failed (exit {rc_}); see {lf}')
            os.replace(part, fin)
            times[name] = round(time.time() - t, 1)
            done = len([p for p in fd.glob('cols_*_*.npz') if not p.name.endswith('.part.npz')])
            self.job.put('fd', {'chunk_wall_s': times, 'chunks_done': done, 'chunks': len(edges) - 1}, status='running')
            self.job.log(f'frequencies: chunk {name} done ({2 * (c1 - c0)} gradients, {times[name] / 60:.1f} min); '
                         f'{done}/{len(edges) - 1} chunks')
        self.assemble(fd, refr, edges, times, t_stage)
        return True

    def assemble(self, fd, refr, edges, times, t_stage):
        """hessian.npy + freq.json exactly as dft/frequencies.py builds them."""
        d, a, case = self.d, self.args, self.case
        cm = published_common()
        from pyscf.hessian import thermo
        import vibrations                                   # dft/vibrations.py (published)
        mol = cm.mol_for_case(case, d / 'opt.xyz', verbose=0)
        ncol = 3 * mol.natm
        hess = np.zeros((mol.natm, mol.natm, 3, 3))
        got = 0
        for c0, c1 in zip(edges[:-1], edges[1:]):
            z = np.load(fd / f'cols_{c0:03d}_{c1:03d}.npz')
            for key in z.files:
                col = int(key); at, kk = divmod(col, 3)
                hess[at, :, kk, :] = z[key]; got += 1
        if got != ncol:
            raise RuntimeError(f'expected {ncol} Hessian columns, collected {got}')
        hess = 0.5 * (hess + hess.transpose(1, 0, 3, 2))
        res = thermo.harmonic_analysis(mol, hess, imaginary_freq=True)
        pyscf_freqs = np.sort(np.array(res['freq_wavenumber']))
        freqs, n_imag_proj, tr_res = vibrations.frequencies(mol, hess)
        if len(freqs) != len(pyscf_freqs) or not np.allclose(np.sort(freqs), pyscf_freqs, atol=1.0):
            stage_note(d, f'frequency projection disagreed with PySCF; using projected set '
                          f'({n_imag_proj} imaginary vs {int((pyscf_freqs < 0).sum())} reported)')

        class _Model:                                       # thermo.thermo reads only .mol and .e_tot
            pass
        model = _Model(); model.mol = mol; model.e_tot = refr['E_Ha']
        th = thermo.thermo(model, res['freq_au'], 298.15, 101325)
        hwall = sum(times.values())
        method = (f'finite-difference (gradients, step {a.fd_step} Bohr, dft/_fd_worker.py in {len(edges) - 1} '
                  f'sequential chunks of {a.fd_chunk} columns x {rc.THREADS} threads; dft/revision/c0_case.py)')
        out = {'case': case['name'], 'E_Ha': float(refr['E_Ha']), 'n_imaginary': n_imag_proj, 'hessian_method': method,
               'n_basis': int(refr['nao']), 'n_atoms': int(mol.natm),
               'frequencies_cm-1': [float(x) for x in np.sort(freqs)], 'hessian_wall_min': hwall / 60,
               'frequencies_pyscf_cm-1': [float(x) for x in pyscf_freqs],
               'n_imaginary_pyscf': int((pyscf_freqs < 0).sum()), 'tr_residual_cm-1': tr_res,
               'ZPE_Ha': float(th['ZPE'][0]), 'H298_Ha': float(th['H_tot'][0]), 'G298_Ha': float(th['G_tot'][0]),
               # revision keys
               'grid_level': refr['grid_level'], 'pcm_lebedev_order': refr['lebedev_order'],
               'fd_chunks': len(edges) - 1, 'fd_dir': str(fd)}
        np.save(d / 'hessian.npy', hess)
        cm.save_json(d / 'freq.json', out)
        note = (f"frequencies ({method}): {n_imag_proj} imaginary, lowest {np.sort(freqs)[0]:.1f} cm-1, "
                f"{hwall / 60:.0f} min")
        stage_note(d, note)
        set_stage(d, case['name'], 'frequencies', {'status': 'done', 'start': t_stage, 'end': now_local(),
                                                   'wall_min': round(hwall / 60, 2), 'note': ''})
        self.job.put('fd', {'chunk_wall_s': times, 'chunks_done': len(edges) - 1, 'chunks': len(edges) - 1})
        self.job.put('frequencies', {'method': 'fd', 'n_imaginary': n_imag_proj,
                                     'lowest_cm-1': [float(x) for x in np.sort(freqs)[:6]],
                                     'tr_residual_cm-1': tr_res, 'hessian_wall_min': round(hwall / 60, 1),
                                     'opt_xyz_sha256': rc.sha256_file(d / 'opt.xyz')})
        self.job.log('frequencies: ' + note + f'; lowest {np.round(np.sort(freqs)[:6], 1)} cm-1; '
                     f'TR residual {tr_res:.2f} cm-1')


def parser():
    p = rc.add_common_args(argparse.ArgumentParser(description=__doc__.split('\n')[1]))
    p.add_argument('--case', required=True, help='btpcz_gas or btpcz_pcm (any <prefix>_<gas|pcm> case)')
    p.add_argument('--stages', default=','.join(STAGES))
    p.add_argument('--freq', choices=('fd', 'softmode', 'none'), default='fd')
    p.add_argument('--fd-chunk', type=int, default=3, help='Cartesian columns per _fd_worker.py call (default 3)')
    p.add_argument('--fd-step', type=float, default=0.005, help='Bohr (frequencies.py default 0.005)')
    p.add_argument('--ref-hessian', default=None, help='--freq softmode: reference hessian.npy (same atom order)')
    p.add_argument('--nmodes', type=int, default=None, help='--freq softmode: modes (default 10; smoke 2)')
    p.add_argument('--soft-step', type=float, default=0.02, help='--freq softmode: largest atomic displacement, Bohr')
    p.add_argument('--preopt-basis', default='6-31g(d)', help="dft/optimize.py's PREOPT_BASIS default")
    p.add_argument('--preopt-maxsteps', type=int, default=300, help='as dft/optimize.py')
    p.add_argument('--preopt-max-hours', type=float, default=None)
    p.add_argument('--opt-maxsteps', type=int, default=150, help='total evaluations of the main optimisation')
    p.add_argument('--opt-max-hours', type=float, default=None, help='total wall hours of the main optimisation')
    p.add_argument('--negligible-stop', type=int, default=10, help='stall window for opt_d3.py (0 = off)')
    return p


def run(job, args):
    case = rc.load_case(args.case)
    if '_' not in case['name'] or case['name'].rsplit('_', 1)[1] not in ('gas', 'pcm'):
        raise SystemExit(f"case name must be <prefix>_gas or <prefix>_pcm, got {case['name']}")
    if bool(case.get('solvent')) != (case['name'].endswith('_pcm')):
        raise SystemExit(f"{case['name']}: the solvent key and the _gas/_pcm suffix disagree")
    if args.smoke and case['name'].split('_')[0] in ('btpc', 'btpcz'):
        raise SystemExit(f"--smoke would write reduced-settings results into the case folder of {case['name']}; "
                         f"use a throwaway case (see the docstring)")
    stages = [s.strip() for s in args.stages.split(',') if s.strip()]
    bad = [s for s in stages if s not in STAGES]
    if bad:
        raise SystemExit(f'unknown stage(s) {bad}; choose from {STAGES}')
    R = Runner(job, args, case)
    job.input_file('case', rc.find_case_file(case['name']))
    job.r['settings'] = {'case': case['name'], 'case_dir': str(R.d), 'stages': stages, 'freq': args.freq,
                         'fd_chunk_columns': args.fd_chunk, 'fd_step_bohr': args.fd_step,
                         'ref_hessian': args.ref_hessian, 'nmodes': args.nmodes, 'soft_step_bohr': args.soft_step,
                         'preopt_basis': 'sto-3g' if args.smoke else args.preopt_basis,
                         'preopt_maxsteps': args.preopt_maxsteps, 'opt_maxsteps': args.opt_maxsteps,
                         'opt_max_hours': args.opt_max_hours, 'negligible_stop_window': args.negligible_stop,
                         'basis': 'case basis (6-311++G(d,p) + overrides)', 'xc': rc.XC,
                         'grid_level': R.grid, 'pcm_lebedev_order': R.leb if case.get('solvent') else None,
                         'solvent': case.get('solvent'), 'density_fit': True,
                         'env_for_published_code': {k: v for k, v in R.env('_probe')[0].items()
                                                    if k in ('GRID_LEVEL', 'PCM_LEBEDEV_ORDER', 'DFT_THREADS',
                                                             'DFT_MAX_MEMORY_MB', 'OMP_NUM_THREADS', 'NO_AUTO_PCM')}}
    shutil.rmtree(job.scratch / '_probe', ignore_errors=True)
    job.save()
    for s in STAGES:
        if s not in stages:
            continue
        if s in ('proton_check', 'frequencies', 'fmo') and not (R.d / 'opt.xyz').exists():
            raise RuntimeError(f'{s}: no opt.xyz in {R.d}; run the optimize stage first')
        if s == 'preopt':
            R.preopt()
        elif s == 'optimize':
            R.optimize()
        elif s == 'proton_check':
            R.pipeline_stage('proton_check', ['proton_check.json'])
        elif s == 'frequencies':
            R.frequencies()
        elif s == 'fmo':
            R.pipeline_stage('fmo', ['fmo.json', 'final.molden', 'esp.cube'])
    job.r['summary'] = {s: (job.get(s) or {}) for s in ('preopt', 'optimize', 'proton_check', 'frequencies', 'fmo')}
    for s in ('optimize', 'proton_check', 'frequencies', 'fmo'):
        job.r['summary'][s] = {k: v for k, v in job.r['summary'][s].items() if k != 'history'}
    job.save()
    return True


if __name__ == '__main__':
    rc.install_sigterm_handler()
    rc.main_wrapper(parser, run, __file__)
