"""The one place that maps the original project layout onto this dataset.

The analysis, table and figure scripts were written in the project workspace, whose layout differs from this
repository. For example, the workspace file `molecules/btpc/results/btpc_pcm/opt.xyz` is
`calculations/btpc_pcm/geometry/opt.xyz` here. The scripts quote the workspace paths in their outputs
(the "source" of every entry in analysis/facts.json, the inputs listed in tables/INDEX.md), so those paths are
kept, and this module translates them:

* `resolve(old)` returns the dataset path of a workspace path (a string such as
  'data/r2/molecules/btpc/results/btpc_gas/logs/optimize.log').
* `legacy_root()` builds `.layout/` (git-ignored): a tree of relative symbolic links that shows the dataset
  files under their workspace paths. The analysis, table and figure scripts run against that tree, so they read
  and write the dataset files without any other change.
* The DFT pipeline (scripts/dft) runs from the repository root and writes a recomputation to `results/<case>/`
  (git-ignored), never into `calculations/`.

Two things are not part of the dataset and are found here:
* the CIF of CCDC 835950: `structure/BTPC.cif`, or the file named by the environment variable BTPC_CIF;
* third-party programs installed by scripts/setup_linux.sh: `tools/` (git-ignored).

docs/provenance.md lists the same mapping as a table.
"""
from __future__ import annotations

import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / 'scripts'
CASES_DIR = REPO / 'cases'
CALC = REPO / 'calculations'
STRUCTURE = REPO / 'structure'
ANALYSIS = REPO / 'analysis'
TABLES = REPO / 'tables'
FIGURES = REPO / 'figures'
TOOLS = REPO / 'tools'            # third-party software (setup_linux.sh); not in the dataset
RESULTS = REPO / 'results'        # where a recomputation with scripts/run_case.sh writes; not in the dataset
LEGACY = REPO / '.layout'         # symbolic-link view in the workspace layout (legacy_root)

MOLECULE = 'btpc'
CASES = ('btpc_gas', 'btpc_pcm')
CLUSTER = STRUCTURE / 'btpc_cluster.xyz'   # rebuilt from the CIF by scripts/dft/build_cluster.py
# The five literature records whose claims the scripts parse (pKa values, earlier H-bond distances): an extract,
# verbatim, of the workspace's literature list papers/btpc/lit/candidates.jsonl.
LITERATURE = ANALYSIS / 'literature_claims.jsonl'
# Driver logs of the two runs (workspace data/r2/logs/) and the case each belongs to.
DRIVER_LOGS = {'btpc_gas.out': 'btpc_gas', 'btpc_mep_retry.out': 'btpc_gas',
               'btpc_pcm.out': 'btpc_pcm', 'btpc_pcm_freq_step0.01.out': 'btpc_pcm'}
# Files of the workspace's git-tracked results folder molecules/btpc/results/<case>/ (a subset of the raw
# results archive data/r2/molecules/btpc/results/<case>/, byte-identical where they overlap).
TRACKED = ('RESULTS.md', 'STAGES.md', 'fmo.json', 'fmo_table.md', 'freq.json', 'freq_step0.005.json',
           'geometry_vs_start.csv', 'mep_extrema.csv', 'mep_range.json', 'mep_surface.json', 'nci_spikes.csv',
           'opt.xyz', 'opt_summary.json', 'proton_check.json', 'status.json')
# File of a results folder -> sub-folder of calculations/<case>/ (patterns are matched in this order).
_SUBDIR = (
    ('geometry', ('preopt.xyz', 'opt.xyz', 'opt_traj.xyz', 'opt_summary.json', 'geometry_vs_start.csv',
                  'proton_check.json')),
    ('wavefunction', ('final.molden', 'final.chk')),
    ('frequencies', ('hessian.npy', 'freq.json', 'freq_step0.005.json', 'softmode_check.json')),
    ('orbitals', ('fmo.json', 'fmo_table.md')),
    ('grids', ('*.cube', '*_plane.txt', 'rdg_scatter.txt')),
    ('multiwfn', ('mwfn_*.in', 'surfanalysis.txt', 'vtx.txt', 'mep_extrema.csv', 'mep_range.json',
                  'mep_surface.json', 'nci_spikes.csv', 'rdg_scatter.npz', '*_ring.npz')),
    ('run', ('RESULTS.md', 'STAGES.md', 'status.json')),
)

# Revision calculations (docs/revision.md). Each job folder calculations/revision1/<job>/ came from the queue folder
# data/revision_runs/<queue>/out/<job>/ of the project: these two jobs ran on the second VM (queue_vm2), all others on
# the first (queue). The second crystal form (structure/btpcz/README.md) has the cases below; its results were in
# data/revision_runs/repo_results/btpcz/<case>/.
REVISION = CALC / 'revision1'
REVISION_VM2_JOBS = ('cc_neutral_631gd_f080', 'cc_salt_631gd_f080')
BTPCZ_CASES = ('btpcz_gas', 'btpcz_pcm')


def cif_path() -> Path:
    """The CIF of CCDC 835950 (not distributed here; see structure/README.md)."""
    return Path(os.environ.get('BTPC_CIF') or STRUCTURE / 'BTPC.cif').expanduser().absolute()


def calc_file(case: str, name: str) -> Path:
    """Dataset path of `name`, a path inside the workspace results folder of `case` (e.g. 'logs/fmo.log')."""
    if name.startswith('logs/'):
        return CALC / case / 'run' / name
    if name == 'mwfn_home/settings.ini':
        return CALC / case / 'multiwfn' / 'settings.ini'
    if '/' not in name:
        for sub, pats in _SUBDIR:
            if any(Path(name).match(p) for p in pats):
                return CALC / case / sub / name
    raise KeyError(f'no dataset location for {case}/{name}')


def _revision_runs(r: str) -> Path:
    """Dataset path of a path below data/revision_runs/ of the project (the revision calculations)."""
    q, _, rest = r.partition('/')
    if q in ('queue', 'queue_vm2'):
        vm = 'vm2' if q == 'queue_vm2' else 'vm1'
        if rest.startswith('out/'):
            job, _, name = rest[len('out/'):].partition('/')
            return REVISION / job / name if name else REVISION / job
        if rest.startswith('logs/') and rest.endswith('.log'):
            return REVISION / rest[len('logs/'):-len('.log')] / 'queue.log'
        if rest in ('status.json', 'runner.log', 'runner.out'):
            return REVISION / 'runner' / vm / rest
    if q in ('vm_extras', 'vm2_extras'):
        return REPO / 'environment' / 'revision1' / ('vm2' if q == 'vm2_extras' else 'vm1') / rest
    if q == 'repo_results' and rest.startswith('btpcz/'):
        case, _, name = rest[len('btpcz/'):].partition('/')
        return calc_file(case, name)
    raise KeyError(f'no dataset location for workspace path data/revision_runs/{r!r}')


def resolve(old: str | os.PathLike) -> Path:
    """Dataset path of a workspace path (relative to the workspace root)."""
    o = Path(old).as_posix().lstrip('./')
    if o.startswith('dft/revision/'):
        return SCRIPTS / 'revision1' / o[len('dft/revision/'):]
    if o.startswith('data/revision_runs/'):
        return _revision_runs(o[len('data/revision_runs/'):])
    if o.startswith('molecules/btpcz/'):
        r = o[len('molecules/btpcz/'):]
        if r.startswith('cases/'):
            return CASES_DIR / r[len('cases/'):]
        if r == 'input/btpcz_cluster.xyz':
            return STRUCTURE / 'btpcz' / 'btpcz_cluster.xyz'
        if r == 'input/2230422.cif':
            return STRUCTURE / '2230422.cif'              # not deposited (structure/btpcz/README.md)
        if r.startswith('results/'):
            case, _, name = r[len('results/'):].partition('/')
            return calc_file(case, name)
    if o.startswith('data/r2/'):
        o = o[len('data/r2/'):]
        if o.startswith('logs/') and o[5:] in DRIVER_LOGS:
            return CALC / DRIVER_LOGS[o[5:]] / 'run' / o[5:]
    m = f'molecules/{MOLECULE}/'
    if o.startswith(m + 'cases/'):
        return CASES_DIR / o[len(m + 'cases/'):]
    if o == m + 'input/btpc_cluster.xyz':
        return CLUSTER
    if o == m + 'input/BTPC.cif':
        return cif_path()
    if o.startswith(m + 'results/'):
        case, _, name = o[len(m + 'results/'):].partition('/')
        return calc_file(case, name)
    if o == 'papers/btpc/lit/candidates.jsonl':
        return LITERATURE
    for prefix, new in (('papers/btpc/data/derived/', ANALYSIS / 'derived'), ('papers/btpc/tables/', TABLES),
                        ('papers/btpc/figures/scripts/', SCRIPTS / 'figures'), ('papers/btpc/figures/', FIGURES),
                        ('dft/', SCRIPTS / 'dft'), ('papers/_shared/figkit/', SCRIPTS / 'figures' / 'figkit'),
                        ('tools/', TOOLS)):
        if o.startswith(prefix):
            return new / o[len(prefix):]
    if o.startswith('papers/btpc/data/'):
        name = o[len('papers/btpc/data/'):]
        return (SCRIPTS / 'analysis' / name) if name.endswith('.py') else (ANALYSIS / name)
    if o in ('run_case.sh', 'setup_linux.sh'):
        return SCRIPTS / o
    raise KeyError(f'no dataset location for workspace path {old!r}')


def _old_new_pairs():
    """Every (workspace path, dataset path) the legacy view contains."""
    m = f'molecules/{MOLECULE}'
    pairs = []
    for base in (m, 'data/r2/' + m):
        for c in CASES:
            pairs.append((f'{base}/cases/{c}.json', CASES_DIR / f'{c}.json'))
        pairs.append((f'{base}/input/btpc_cluster.xyz', CLUSTER))
        pairs.append((f'{base}/input/BTPC.cif', cif_path()))
    for c in CASES:
        root = CALC / c
        for p in sorted(root.rglob('*')):
            if not p.is_file():
                continue
            rel = p.relative_to(root)
            if rel.parts[0] == 'run' and len(rel.parts) == 2 and rel.name in DRIVER_LOGS:
                pairs.append((f'data/r2/logs/{rel.name}', p))
                continue
            if rel.parts[0] == 'run' and rel.parts[1] == 'logs':
                name = 'logs/' + rel.name
            elif rel.parts[0] == 'multiwfn' and rel.name == 'settings.ini':
                name = 'mwfn_home/settings.ini'
            elif rel.parts[0] == 'grids' and rel.suffix == '.xz':
                name = rel.stem                       # the unpacked grid (scripts/unpack_grids.sh)
                p = p.with_suffix('')
            elif rel.parts[0] == 'grids':
                continue                              # the unpacked copy itself, linked via its .xz above
            else:
                name = rel.name
            pairs.append((f'data/r2/{m}/results/{c}/{name}', p))
            if name in TRACKED:
                pairs.append((f'{m}/results/{c}/{name}', p))
    # revision calculations: the job folders under their queue folders, the VM environment records, the second
    # crystal form's results, cases and starting geometry (docs/revision.md)
    if REVISION.is_dir():
        for d in sorted(REVISION.iterdir()):
            if d.is_dir() and d.name != 'runner':
                q = 'queue_vm2' if d.name in REVISION_VM2_JOBS else 'queue'
                pairs.append((f'data/revision_runs/{q}/out/{d.name}', d))
    pairs += [('data/revision_runs/vm_extras', REPO / 'environment' / 'revision1' / 'vm1'),
              ('data/revision_runs/vm2_extras', REPO / 'environment' / 'revision1' / 'vm2')]
    for c in BTPCZ_CASES:
        pairs.append((f'molecules/btpcz/cases/{c}.json', CASES_DIR / f'{c}.json'))
        root = CALC / c
        for p in sorted(root.rglob('*')) if root.is_dir() else []:
            if not p.is_file():
                continue
            rel = p.relative_to(root)
            if rel.parts[0] == 'run' and rel.parts[1] == 'logs':
                name = 'logs/' + rel.name
            elif rel.parts[0] == 'multiwfn' and rel.name == 'settings.ini':
                name = 'mwfn_home/settings.ini'
            elif rel.parts[0] == 'grids' and rel.suffix == '.xz':
                name = rel.stem
                p = p.with_suffix('')
            elif rel.parts[0] == 'grids':
                continue
            else:
                name = rel.name
            pairs.append((f'data/revision_runs/repo_results/btpcz/{c}/{name}', p))
    pairs.append(('molecules/btpcz/input/btpcz_cluster.xyz', STRUCTURE / 'btpcz' / 'btpcz_cluster.xyz'))
    pairs.append(('papers/btpc/data/experimental.json', ANALYSIS / 'experimental.json'))
    pairs.append(('papers/btpc/data/rev', SCRIPTS / 'analysis' / 'rev'))
    for name in ('facts.json', 'facts.md', 'derived'):
        pairs.append((f'papers/btpc/data/{name}', ANALYSIS / name))
    for name in ('build_facts.py', 'build_facts_md.py', 'wfn_topology.py'):
        pairs.append((f'papers/btpc/data/{name}', SCRIPTS / 'analysis' / name))
    pairs.append(('papers/btpc/lit/candidates.jsonl', LITERATURE))
    pairs += [('papers/btpc/tables', TABLES), ('papers/btpc/figures', FIGURES),
              ('papers/_shared/figkit', SCRIPTS / 'figures' / 'figkit'),
              ('papers/_shared/tools', SCRIPTS / 'tables' / 'pandoc'),
              ('dft', SCRIPTS / 'dft'), ('tools', TOOLS),
              ('run_case.sh', SCRIPTS / 'run_case.sh'), ('setup_linux.sh', SCRIPTS / 'setup_linux.sh')]
    return pairs


def legacy_root(quiet: bool = False) -> Path:
    """Build (or refresh) the symbolic-link view `.layout/` and return its root."""
    FIGURES.mkdir(exist_ok=True)
    (ANALYSIS / 'derived').mkdir(parents=True, exist_ok=True)
    missing = []
    for old, new in _old_new_pairs():
        link = LEGACY / old
        link.parent.mkdir(parents=True, exist_ok=True)
        target = os.path.relpath(new, link.parent)
        if link.is_symlink():
            if os.readlink(link) == target:
                pass
            else:
                link.unlink()
        if not link.is_symlink():
            if link.exists():
                raise FileExistsError(f'{link} exists and is not a link')
            os.symlink(target, link)
        if not new.exists() and new.suffix in ('.cube', '.txt') and 'grids' in new.parts:
            missing.append(new.relative_to(REPO))
    if missing and not quiet:
        print(f'[paths] {len(missing)} grids are not unpacked yet (run scripts/unpack_grids.sh), e.g. {missing[0]}')
    return LEGACY
