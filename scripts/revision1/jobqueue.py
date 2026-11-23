"""
A small local job queue for the BTPC revision calculations (one machine, no daemon, no root).

Named jobqueue.py, not queue.py: a queue.py in dft/revision/ shadows Python's standard-library `queue` module for every
script run from this directory (the script's directory is sys.path[0]).  That was seen on 20 Nov 2026: geomeTRIC's
NetworkX import failed ("Cannot import optional NetworkX module, topology tools won't work"), which would break the TRIC
coordinates of every optimisation.

  python dft/revision/jobqueue.py run    --jobs jobs.yaml --root <queue dir> [--slots N] [--threads T]
                                       [--scratch-root DIR] [--poll S]
  python dft/revision/jobqueue.py status --root <queue dir> [--json]
  python dft/revision/jobqueue.py retry  --root <queue dir> ID [ID ...]   (failed/blocked -> pending)
  python dft/revision/jobqueue.py stop   --root <queue dir> [--kill]      (no new starts; --kill also SIGTERMs jobs)

Jobs file (YAML or JSON): either a list of jobs or {"defaults": {...}, "jobs": [...]}.  Each job:
  id           unique name (also the scratch directory name)
  command      shell command (bash -c); placeholders {id} {scratch} {root} {repo} {python} are substituted
  threads      threads the job uses (exported as THREADS, OMP/MKL/OpenBLAS/NUMEXPR_NUM_THREADS); default 1
  depends_on   list of job ids that must finish successfully first; default []
  est_hours    estimated wall hours (status/ETA only); default 0
  cwd          working directory; default the repository root
  env          extra environment variables {NAME: value}
  retries      automatic re-runs after a non-zero exit; default 0 (an interrupted job is always re-queued)

Scheduling: a job starts when its dependencies are done, a slot is free (--slots, default 1) and its threads fit in
the thread budget (--threads, default all CPUs).  Each job gets PYSCF_TMPDIR=<scratch-root>/<id> (default
<root>/scratch/<id>), removed after success and kept after a failure, and its output in <root>/logs/<id>.log.
A job that depends on a failed job is marked blocked.

Crash safety: jobs run in their own session, so they survive the runner; each writes its exit code to
<root>/rc/<id>.rc when it ends.  The command runs in a subshell, so an `exit` inside it (e.g. `cmd || exit $?`) cannot
skip that file (21 Nov 2026: before, a C0 loop's `|| exit $?` left no rc file, and a restarted runner then re-queued a
failed job as if interrupted).  A job whose wrapper shell is killed by a signal (`stop --kill`, a reboot) is treated
as interrupted and re-queued, not counted as a failed attempt.  A restarted runner adopts jobs that are still alive, records the exit code of
those that finished meanwhile, and re-queues those that died without one (the revision scripts are restartable).
State lives in <root>/status.json (atomic writes); <root>/runner.lock stops two runners on one queue.
Times are shown in IST.
"""
import os, sys, json, time, shlex, signal, fcntl, argparse, pathlib, datetime, subprocess, shutil

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
IST = datetime.timezone(datetime.timedelta(hours=5, minutes=30), 'IST')
STATES = ('pending', 'running', 'done', 'failed', 'blocked')


def now():
    return datetime.datetime.now(IST)


def iso(t=None):
    return (t or now()).isoformat(timespec='seconds')


def parse_iso(s):
    return datetime.datetime.fromisoformat(s) if s else None


def load_jobs(path):
    text = open(path).read()
    if str(path).endswith(('.yaml', '.yml')):
        import yaml
        data = yaml.safe_load(text)
    else:
        data = json.loads(text)
    defaults = {}
    if isinstance(data, dict):
        defaults = data.get('defaults', {}) or {}
        data = data.get('jobs', [])
    jobs, seen = [], set()
    for j in data:
        j = dict(defaults, **j)
        if 'id' not in j or 'command' not in j:
            raise SystemExit(f'job without id/command: {j}')
        if j['id'] in seen:
            raise SystemExit(f'duplicate job id {j["id"]}')
        seen.add(j['id'])
        j.setdefault('threads', 1); j.setdefault('depends_on', []); j.setdefault('est_hours', 0.0)
        j.setdefault('retries', 0); j.setdefault('env', {})
        if isinstance(j['depends_on'], str):
            j['depends_on'] = [j['depends_on']]
        jobs.append(j)
    ids = {j['id'] for j in jobs}
    for j in jobs:
        missing = [d for d in j['depends_on'] if d not in ids]
        if missing:
            raise SystemExit(f'job {j["id"]} depends on unknown job(s) {missing}')
    # cycle check
    order, state = [], {}
    deps = {j['id']: j['depends_on'] for j in jobs}

    def visit(n, stack):
        if state.get(n) == 1:
            raise SystemExit(f'dependency cycle: {" -> ".join(stack + [n])}')
        if state.get(n) == 2:
            return
        state[n] = 1
        for d in deps[n]:
            visit(d, stack + [n])
        state[n] = 2; order.append(n)
    for n in deps:
        visit(n, [])
    return jobs


class Queue:
    def __init__(self, root, scratch_root=None):
        self.root = pathlib.Path(root).resolve()
        for sub in ('logs', 'rc'):
            (self.root / sub).mkdir(parents=True, exist_ok=True)
        self.scratch_root = pathlib.Path(scratch_root).resolve() if scratch_root else self.root / 'scratch'
        self.status_path = self.root / 'status.json'
        self.st = json.load(open(self.status_path)) if self.status_path.exists() else {'jobs': {}, 'order': []}

    # ---------------------------------------------------------------- persistence
    def save(self):
        self.st['updated'] = iso()
        tmp = self.status_path.with_name('status.json.tmp')
        with open(tmp, 'w') as f:
            json.dump(self.st, f, indent=2)
        os.replace(tmp, self.status_path)

    def log(self, msg):
        line = f"[{now().strftime('%Y-%m-%d %I:%M:%S %p IST')}] {msg}"
        print(line, flush=True)
        with open(self.root / 'runner.log', 'a') as f:
            f.write(line + '\n')

    def merge_jobs(self, jobs):
        """Add new jobs from the jobs file and refresh the definition of jobs that have not started."""
        for j in jobs:
            cur = self.st['jobs'].get(j['id'])
            if cur is None:
                self.st['jobs'][j['id']] = {'def': j, 'state': 'pending', 'attempts': 0}
                self.st['order'].append(j['id'])
            elif cur['state'] in ('pending', 'blocked'):
                cur['def'] = j
        self.save()

    # ---------------------------------------------------------------- process helpers
    @staticmethod
    def proc_start_ticks(pid):
        try:
            with open(f'/proc/{pid}/stat') as f:
                return int(f.read().rsplit(')', 1)[1].split()[19])
        except (OSError, IndexError, ValueError):
            return None

    def alive(self, rec):
        pid = rec.get('pid')
        if not pid:
            return False
        ticks = self.proc_start_ticks(pid)
        if ticks is None or ticks != rec.get('pid_start_ticks'):
            return False
        try:
            with open(f'/proc/{pid}/stat') as f:
                if f.read().rsplit(')', 1)[1].split()[0] == 'Z':
                    return False
        except OSError:
            return False
        return True

    def rc_of(self, jid):
        p = self.root / 'rc' / f'{jid}.rc'
        if p.exists():
            try:
                return int(p.read_text().strip())
            except ValueError:
                return 255
        return None

    # ---------------------------------------------------------------- job lifecycle
    def start(self, jid):
        rec = self.st['jobs'][jid]; j = rec['def']
        scratch = self.scratch_root / jid
        scratch.mkdir(parents=True, exist_ok=True)
        rcfile = self.root / 'rc' / f'{jid}.rc'
        rcfile.unlink(missing_ok=True)
        cmd = str(j['command'])
        for k, v in (('{id}', jid), ('{scratch}', str(scratch)), ('{root}', str(self.root)),
                     ('{repo}', str(ROOT)), ('{python}', sys.executable)):
            cmd = cmd.replace(k, v)
        t = str(int(j['threads']))
        env = dict(os.environ)
        env.update({k: str(v) for k, v in (j.get('env') or {}).items()})
        env.update({'THREADS': t, 'OMP_NUM_THREADS': t, 'MKL_NUM_THREADS': t, 'OPENBLAS_NUM_THREADS': t,
                    'NUMEXPR_NUM_THREADS': t, 'VECLIB_MAXIMUM_THREADS': t,
                    'PYSCF_TMPDIR': str(scratch), 'TMPDIR': str(scratch),
                    'QUEUE_JOB_ID': jid, 'QUEUE_ROOT': str(self.root)})
        wrapper = (f'(\n{cmd}\n)\nrc=$?\necho $rc > {shlex.quote(str(rcfile))}.tmp && '
                   f'mv {shlex.quote(str(rcfile))}.tmp {shlex.quote(str(rcfile))}\nexit $rc\n')
        logf = open(self.root / 'logs' / f'{jid}.log', 'a')
        rec['attempts'] = rec.get('attempts', 0) + 1
        logf.write(f'\n===== queue: start {jid} attempt {rec["attempts"]} at {iso()} =====\n$ {cmd}\n'); logf.flush()
        p = subprocess.Popen(['bash', '-c', wrapper], cwd=j.get('cwd') or str(ROOT), env=env, stdout=logf,
                             stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True)
        logf.close()
        rec.update({'state': 'running', 'pid': p.pid, 'pgid': p.pid, 'pid_start_ticks': self.proc_start_ticks(p.pid),
                    'start': iso(), 'end': None, 'rc': None, 'command': cmd, 'scratch': str(scratch),
                    'host': os.uname().nodename})
        self.save()
        self._children = getattr(self, '_children', {})
        self._children[jid] = p
        self.log(f'start {jid} (pid {p.pid}, {t} threads, est {j["est_hours"]} h)')

    def finish(self, jid, rc, note=''):
        rec = self.st['jobs'][jid]; j = rec['def']
        rec['end'] = iso(); rec['rc'] = rc
        st = parse_iso(rec.get('start'))
        rec['wall_h'] = round((now() - st).total_seconds() / 3600, 3) if st else None
        if rc == 0:
            rec['state'] = 'done'
            if rec.get('note', '').startswith(('exit ', 'interrupted')):
                rec['note'] = f"done after {rec['attempts']} attempt(s); earlier: {rec['note']}"
            shutil.rmtree(rec.get('scratch') or '', ignore_errors=True)
        elif rc is None:                                  # died without an exit code: interrupted
            rec['state'] = 'pending'; rec['note'] = f'interrupted at {iso()}; re-queued'
        elif rec['attempts'] <= int(j.get('retries', 0)):
            rec['state'] = 'pending'; rec['note'] = f'exit {rc}; retry {rec["attempts"]}/{j.get("retries")}'
        else:
            rec['state'] = 'failed'
        if note:
            rec['note'] = note
        rec.pop('pid', None)
        self.save()
        self.log(f'end {jid}: rc {rc} -> {rec["state"]} ({rec.get("wall_h")} h){" " + note if note else ""}')

    def reap(self):
        """Collect finished jobs (children of this runner, or adopted ones via their rc file)."""
        ch = getattr(self, '_children', {})
        for jid, rec in self.st['jobs'].items():
            if rec['state'] != 'running':
                continue
            p = ch.get(jid)
            if p is not None:
                if p.poll() is None:
                    continue
                rc = self.rc_of(jid)
                if rc is None and p.returncode is not None and p.returncode >= 0:
                    rc = p.returncode
                # rc None: the wrapper shell itself was killed by a signal -> interrupted, re-queued
                self.finish(jid, rc)
                ch.pop(jid, None)
            elif not self.alive(rec):
                rc = self.rc_of(jid)
                self.finish(jid, rc, '' if rc is not None else 'process gone without an exit code')

    def adopt(self):
        for jid, rec in self.st['jobs'].items():
            if rec['state'] == 'running':
                if self.alive(rec):
                    self.log(f'adopting {jid} (pid {rec["pid"]} still running)')
                else:
                    rc = self.rc_of(jid)
                    self.finish(jid, rc, 'finished while no runner was active' if rc is not None
                                else 'runner restart: process gone without an exit code')

    def block_dependents(self):
        changed = True
        while changed:
            changed = False
            for jid, rec in self.st['jobs'].items():
                if rec['state'] != 'pending':
                    continue
                bad = [d for d in rec['def']['depends_on'] if self.st['jobs'][d]['state'] in ('failed', 'blocked')]
                if bad:
                    rec['state'] = 'blocked'; rec['note'] = f'dependency {bad} failed'
                    changed = True
                    self.log(f'blocked {jid}: dependency {bad} failed')
        self.save()

    def runnable(self):
        return [jid for jid in self.st['order'] if self.st['jobs'][jid]['state'] == 'pending'
                and all(self.st['jobs'][d]['state'] == 'done' for d in self.st['jobs'][jid]['def']['depends_on'])]

    def running(self):
        return [jid for jid, r in self.st['jobs'].items() if r['state'] == 'running']

    def run(self, slots, threads, poll):
        self.st['runner'] = {'pid': os.getpid(), 'host': os.uname().nodename, 'started': iso(), 'slots': slots,
                             'threads': threads, 'scratch_root': str(self.scratch_root)}
        self.save()
        self.adopt()
        stop = self.root / 'STOP'
        while True:
            self.reap()
            self.block_dependents()
            run_ids = self.running()
            used = sum(int(self.st['jobs'][j]['def']['threads']) for j in run_ids)
            if not stop.exists():
                for jid in self.runnable():
                    need = int(self.st['jobs'][jid]['def']['threads'])
                    if len(run_ids) >= slots:
                        break
                    if used + need > threads and run_ids:
                        continue                      # wait for threads (a lone over-wide job still runs)
                    self.start(jid); run_ids.append(jid); used += need
            run_ids = self.running()
            if not run_ids and (stop.exists() or not self.runnable()):
                break
            time.sleep(poll)
        counts = {s: sum(1 for r in self.st['jobs'].values() if r['state'] == s) for s in STATES}
        self.st['runner']['ended'] = iso(); self.save()
        self.log(f'runner exits: {counts}{" (STOP file present)" if stop.exists() else ""}')
        return counts


# ---------------------------------------------------------------- status
def status(root, as_json=False):
    q = Queue(root)
    jobs = q.st['jobs']
    if as_json:
        print(json.dumps(q.st, indent=2)); return
    t = now()
    rows = []
    for jid in q.st['order']:
        r = jobs[jid]; d = r['def']
        if r['state'] == 'running':
            el = (t - parse_iso(r['start'])).total_seconds() / 3600
            alive = q.alive(r)
            wall = f'{el:.2f}' + ('' if alive else ' (dead?)')
        else:
            wall = f"{r['wall_h']:.2f}" if r.get('wall_h') is not None else ''
        rows.append((jid, r['state'], str(d['threads']), f"{float(d['est_hours']):.2f}", wall,
                     ','.join(d['depends_on']) or '-', '' if r.get('rc') is None else str(r['rc']),
                     (parse_iso(r['start']).strftime('%d %b %I:%M %p') if r.get('start') else ''),
                     r.get('note', '')))
    hdr = ('job', 'state', 'thr', 'est_h', 'wall_h', 'depends_on', 'rc', 'started (IST)', 'note')
    w = [max(len(h), *(len(x[i]) for x in rows)) if rows else len(h) for i, h in enumerate(hdr)]
    print('  '.join(h.ljust(w[i]) for i, h in enumerate(hdr)))
    for x in rows:
        print('  '.join(c.ljust(w[i]) for i, c in enumerate(x)))
    counts = {s: sum(1 for r in jobs.values() if r['state'] == s) for s in STATES}
    left = sum(float(r['def']['est_hours']) for r in jobs.values() if r['state'] in ('pending', 'running'))
    total = sum(float(r['def']['est_hours']) for r in jobs.values()) or 1.0
    done_h = sum(float(r['def']['est_hours']) for r in jobs.values() if r['state'] == 'done')
    rn = q.st.get('runner', {})
    slots = rn.get('slots') or 1
    runner_alive = False
    if rn.get('pid') and not rn.get('ended'):
        try:
            os.kill(rn['pid'], 0); runner_alive = rn.get('host') == os.uname().nodename
        except OSError:
            pass
    print(f"\n{counts}; progress {100 * done_h / total:.0f} % of estimated hours; remaining est {left:.1f} h "
          f"(~{left / slots:.1f} h on {slots} slot(s)); runner {'alive (pid %s)' % rn.get('pid') if runner_alive else 'not running'}"
          f"; updated {q.st.get('updated', '?')}")


def main():
    ap = argparse.ArgumentParser(description='local job queue for the BTPC revision')
    sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run'); r.add_argument('--jobs', required=True); r.add_argument('--root', required=True)
    r.add_argument('--slots', type=int, default=1); r.add_argument('--threads', type=int, default=os.cpu_count())
    r.add_argument('--scratch-root'); r.add_argument('--poll', type=float, default=5.0)
    s = sub.add_parser('status'); s.add_argument('--root', required=True); s.add_argument('--json', action='store_true')
    t = sub.add_parser('retry'); t.add_argument('--root', required=True); t.add_argument('ids', nargs='+')
    k = sub.add_parser('stop'); k.add_argument('--root', required=True); k.add_argument('--kill', action='store_true')
    a = ap.parse_args()

    if a.cmd == 'status':
        return status(a.root, a.json)
    if a.cmd == 'stop':
        q = Queue(a.root)
        (q.root / 'STOP').write_text(iso() + '\n')
        if a.kill:
            for jid, rec in q.st['jobs'].items():
                if rec['state'] == 'running' and q.alive(rec):
                    os.killpg(rec['pgid'], signal.SIGTERM)
                    print(f'SIGTERM sent to {jid} (process group {rec["pgid"]})')
        print(f'STOP file written; the runner starts no new jobs (remove {q.root / "STOP"} to resume)')
        return
    lockf = open(pathlib.Path(a.root).resolve().joinpath('runner.lock') if pathlib.Path(a.root).exists()
                 else (pathlib.Path(a.root).mkdir(parents=True, exist_ok=True) or pathlib.Path(a.root).resolve() / 'runner.lock'), 'w')
    try:
        fcntl.flock(lockf, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit(f'another runner holds {a.root}/runner.lock')
    if a.cmd == 'retry':
        q = Queue(a.root)
        for jid in a.ids:
            rec = q.st['jobs'][jid]
            if rec['state'] in ('failed', 'blocked'):
                rec['state'] = 'pending'; rec['attempts'] = 0; rec['note'] = f'retry requested {iso()}'
                print(f'{jid}: pending')
            else:
                print(f'{jid}: state {rec["state"]}, unchanged')
        for jid, rec in q.st['jobs'].items():          # re-evaluated by the next runner (re-blocked if still needed)
            if rec['state'] == 'blocked':
                rec['state'] = 'pending'; rec['note'] = ''
        q.save()
        return
    q = Queue(a.root, a.scratch_root)
    (q.root / 'STOP').unlink(missing_ok=True)
    q.merge_jobs(load_jobs(a.jobs))
    counts = q.run(a.slots, a.threads, a.poll)
    sys.exit(0 if counts['failed'] == 0 and counts['blocked'] == 0 and counts['pending'] == 0 else 1)


if __name__ == '__main__':
    main()
