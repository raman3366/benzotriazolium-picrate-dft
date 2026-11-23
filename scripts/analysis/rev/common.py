"""Shared helpers for the revision fact modules (read-only access to data/revision_runs)."""
import json, pathlib, sys, importlib, os, hashlib
# dataset: the deposited files are read through the workspace-layout view of scripts/paths.py, in which
# data/revision_runs/<queue>/out/<job>/ is calculations/revision1/<job>/ (docs/revision.md)
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))   # scripts/ (paths.py)
import paths
ROOT = paths.legacy_root(quiet=True)
RUNS = ROOT / 'data/revision_runs'
QUEUES = ('queue_vm2', 'queue')          # later run first, so a job id present in both resolves to the newer one
HA_EV = 27.211386245988; HA_KCAL = 627.509474; HA_KJ = 2625.499639; BOHR = 0.529177210903; AU2CM = 219474.6313702

def job_dir(job):
    for q in QUEUES:
        p = RUNS / q / 'out' / job
        if (p / 'result.json').exists():
            return p
    raise FileNotFoundError(f'no result.json for job {job} under {RUNS}/{{{",".join(QUEUES)}}}/out')

def result(job):
    """(dict, repo-relative source path) of one revision job's result.json."""
    p = job_dir(job) / 'result.json'
    return json.load(open(p)), str(p.relative_to(ROOT))

def rel(p):
    return str(pathlib.Path(os.path.abspath(p)).relative_to(ROOT))   # dataset: the view is made of links, so no resolve()

_SCRUBBED = None
def sha256_as_run(p):
    """SHA-256 of a deposited file as the run wrote it (dataset). Text files in which paths, host names or the user
    name were replaced for deposit (docs/provenance.md, rules D1-D4) are listed in docs/scrubbed_files.tsv with
    both checksums; for those the checksum before the replacement is returned, for every other file its own."""
    global _SCRUBBED
    if _SCRUBBED is None:
        _SCRUBBED = {}
        for line in open(paths.REPO / 'docs' / 'scrubbed_files.tsv', encoding='utf-8'):
            f = line.rstrip('\n').split('\t')
            if len(f) >= 3 and f[0] != 'path':
                _SCRUBBED[f[0]] = (f[1], f[2])
    h = hashlib.sha256(open(p, 'rb').read()).hexdigest()
    before, after = _SCRUBBED.get(os.path.relpath(os.path.realpath(p), paths.REPO), (None, None))
    return before if h == after else h

def standalone(modname):
    """Run one module into a fresh dict and print its facts; nothing is written."""
    facts = {}
    def add(key, value, unit='', case=None, source='', note=''):
        assert key.startswith('rev.'), key
        assert key not in facts, key
        facts[key] = {'value': value, 'unit': unit, 'case': case, 'source': source, 'note': note}
    m = importlib.import_module('rev.' + modname)
    m.build(add, facts)
    for k, v in facts.items():
        val = v['value']
        s = json.dumps(val, ensure_ascii=False)
        print(f"{k} = {s[:160]}{'…' if len(s) > 160 else ''} [{v['unit']}] ({v['case']})")
    print(f'{len(facts)} facts from rev.{modname}')

if __name__ == '__main__':
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
    standalone(sys.argv[1])
