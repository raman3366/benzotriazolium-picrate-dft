"""
Import-level check that dft/pipeline.py would accept a molecule's case files (no SCF, no optimisation).

  python dft/revision/tests/check_cases_pipeline.py btpcz [btpc]

For each <id>_gas / <id>_pcm case: loads it with dft/common.py load_case (as every pipeline stage does), imports
dft/pipeline.py, checks every case key the dft/*.py stages read, builds the Mole at the input geometry for the main
basis and the 6-31G(d) pre-optimisation basis, builds (does not run) common.make_mf, compares it with
rev_common.make_mf, checks the watched hydrogen bond on the input geometry, and checks that pipeline.make_pcm_case
on the gas case reproduces the stored _pcm case.  Side effect, as with any pipeline call: load_case creates the empty
directory molecules/<id>/results/<case>/ (removed again here if this script created it).
"""
import sys, os, json, shutil, pathlib, importlib, tempfile
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import rev_common as rc                     # pins threads
ROOT = rc.ROOT
os.environ.setdefault('GRID_LEVEL', '4'); os.environ.setdefault('PCM_LEBEDEV_ORDER', '29')
sys.path.insert(0, str(ROOT / 'dft'))
import common                               # noqa: E402  (dft/common.py, the pipeline's own module)
import pipeline                             # noqa: E402  (import only: defines the stages)
import numpy as np                          # noqa: E402

REQUIRED = ['title', 'xyz', 'charge', 'spin', 'cif_labels', 'fragments', 'hbonds', 'distances', 'planes', 'proton',
            'chemical_name', 'formula', 'source_doi', 'notes']


def check(name):
    out = {'case': name}
    res_dir = ROOT / 'molecules' / name.rsplit('_', 1)[0] / 'results' / name
    existed, parent_existed = res_dir.exists(), res_dir.parent.exists()
    case = common.load_case(name)
    try:
        missing = [k for k in REQUIRED if k not in case]
        assert not missing, f'missing keys {missing}'
        syms, xyz = common.read_xyz(ROOT / case['xyz'])
        n = len(syms)
        assert len(case['cif_labels']) == n
        assert sorted(i for v in case['fragments'].values() for i in v) == list(range(n))
        for k, idx in list(case['hbonds'].items()) + list(case['distances'].items()) + list(case['planes'].items()):
            assert all(0 <= i < n for i in idx), k
        p = case['proton']
        hb = common.hbond_geometry(xyz, p['donor_N'], p['H'], p['acceptor_O'])
        out['watched_hbond'] = {'label': f"{case['cif_labels'][p['donor_N']]}-{case['cif_labels'][p['H']]}..."
                                         f"{case['cif_labels'][p['acceptor_O']]}", **{k: round(v, 3) for k, v in hb.items()}}
        mol = common.mol_for_case(case, ROOT / case['xyz'], verbose=0)
        pre = common.build_mol(ROOT / case['xyz'], case['charge'], case['spin'], basis='6-31g(d)',
                               basis_overrides=case.get('basis_overrides'), verbose=0)
        mf = common.make_mf(mol, solvent=case.get('solvent'))
        out.update({'natoms': n, 'nao': int(mol.nao), 'nao_preopt_631gd': int(pre.nao), 'nelectron': int(mol.nelectron),
                    'solvent': case.get('solvent'), 'mf': type(mf).__name__,
                    'basis': common.basis_description(case),
                    'settings_vs_rev_common': rc.compare_with_common(mol, case.get('solvent')) or 'identical'})
        assert mol.nelectron % 2 == 0 and case['spin'] == 0
        if not case.get('solvent'):
            # make_pcm_case writes <id>_pcm.json next to the gas case: do it in a temporary copy and compare
            with tempfile.TemporaryDirectory() as td:
                tmp_cases = pathlib.Path(td) / 'cases'; tmp_cases.mkdir()
                shutil.copy(case['case_file'], tmp_cases / pathlib.Path(case['case_file']).name)
                orig = pipeline.find_case_file
                pipeline.find_case_file = lambda nm: tmp_cases / f'{nm}.json'
                try:
                    pcm_name = pipeline.make_pcm_case(case)
                finally:
                    pipeline.find_case_file = orig
                made = json.load(open(tmp_cases / f'{pcm_name}.json'))
                stored = json.load(open(common.find_case_file(pcm_name)))
                out['make_pcm_case_matches_stored_pcm'] = made == stored
                if made != stored:
                    out['make_pcm_case_diff'] = sorted(k for k in set(made) | set(stored) if made.get(k) != stored.get(k))
        out['ok'] = True
    finally:
        if not existed and res_dir.exists() and not any(res_dir.iterdir()):
            res_dir.rmdir()
        if not parent_existed and res_dir.parent.exists() and not any(res_dir.parent.iterdir()):
            res_dir.parent.rmdir()
    return out


if __name__ == '__main__':
    ok = True
    for mid in sys.argv[1:] or ['btpcz']:
        for suffix in ('gas', 'pcm'):
            try:
                r = check(f'{mid}_{suffix}')
            except Exception as e:
                r = {'case': f'{mid}_{suffix}', 'ok': False, 'error': f'{type(e).__name__}: {e}'}
            ok &= r.get('ok', False)
            print(json.dumps(r, default=str))
    sys.exit(0 if ok else 1)
