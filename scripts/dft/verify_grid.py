"""Which integration grid / PCM cavity did each finished case actually use?

report.py hard-codes "integration grid level 3" into every RESULTS.md and nothing records the
setting in machine-readable form, so the delivered methods paragraph cannot be trusted.  Recompute
the SCF at the stored optimised geometry under each candidate setting and see which one reproduces
the stored total energy.  The match is the setting that was really used.
"""
import sys, os, json, pathlib, itertools
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import paths

case_name = sys.argv[1]
os.environ.setdefault('DFT_THREADS', '8')

import numpy as np

results = []
for grid in (3, 4):
    lebs = (17, 29)
    for leb in lebs:
        os.environ['GRID_LEVEL'] = str(grid)
        os.environ['PCM_LEBEDEV_ORDER'] = str(leb)
        for m in [k for k in list(sys.modules) if k in ('common',)]:
            del sys.modules[m]
        import common
        case = common.load_case(case_name)
        if not case.get('solvent') and leb != lebs[0]:
            continue                      # cavity is irrelevant without a solvent
        mol = common.mol_for_case(case, paths.calc_file(case_name, 'opt.xyz'))   # the archived geometry
        mf = common.make_mf(mol, conv_tol=1e-9, solvent=case.get('solvent'))
        e = mf.kernel()
        results.append((grid, leb if case.get('solvent') else None, float(e)))
        print(f'  grid {grid}  lebedev {leb if case.get("solvent") else "-":>3}  E = {e:.8f}', flush=True)

stored = json.load(open(paths.calc_file(case_name, 'opt_summary.json')))['E_final_Ha']
print(f'\n  stored E_final_Ha = {stored:.8f}')
best = min(results, key=lambda r: abs(r[2] - stored))
print(f'  MATCH -> grid {best[0]}, lebedev {best[1]}, |dE| = {abs(best[2]-stored)*1e6:.2f} microHa')
