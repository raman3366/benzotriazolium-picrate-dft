"""
Smoke test of the dft/revision/rev_common.py paths that validate.py does not exercise (minutes, reduced settings).

  THREADS=2 python dft/revision/tests/smoke_rev_common.py --job smoke_rev_common --out <dir> [--smoke]

STO-3G ion pair (btpc_gas geometry), grid level 2, Lebedev 17, conv_tol 1e-6:
  1. ΔSCF-style UKS cation (N-1, charge +1, doublet) and anion (N+1, charge -1, doublet): fragment Mulliken/Loewdin
     charges and spin populations, alpha/beta frontier-orbital shares;
  2. dispersion: B3LYP-D3(BJ) and wB97X-D3(BJ) energies and the dispersion part of the analytic gradient;
  3. a numeric permittivity (eps = 4) through the solvent argument;
  4. the restart path: a second call with the same --out skips every finished item.
Always runs with reduced settings (the --smoke flag is accepted for the common interface and recorded).
"""
import sys, pathlib, argparse, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import rev_common as rc
import numpy as np

KW = dict(conv_tol=1e-6, grid_level=2, lebedev_order=17)


def parser():
    return rc.add_common_args(argparse.ArgumentParser(description='rev_common smoke test'))


def run(job, args):
    case = rc.load_case('btpc_gas')
    geom = rc.ROOT / 'molecules/btpc/results/btpc_gas/opt.xyz'
    job.input_file('geometry', geom)
    frags = rc.fragments(case)
    job.r['settings'] = {'basis': 'sto-3g', **KW, 'fragments': {k: len(v) for k, v in frags.items()}}
    skipped = 0
    for key, charge, spin, xc, disp, solvent in [
            ('uks_cation_Nminus1', 1, 1, 'b3lyp', None, None),
            ('uks_anion_Nplus1', -1, 1, 'b3lyp', None, 'water'),
            ('b3lyp_d3bj', 0, 0, 'b3lyp', 'd3bj', None),
            ('wb97x_d3bj', 0, 0, 'wb97x', 'd3bj', None),
            ('b3lyp_eps4', 0, 0, 'b3lyp', None, 4.0)]:
        if job.done(key):
            skipped += 1; job.log(f'{key}: done earlier, skipped'); continue
        mol = rc.mol_for_case(case, geom, charge=charge, spin=spin, basis='sto-3g', verbose=0)
        mf = rc.make_mf(mol, xc=xc, disp=disp, solvent=solvent, chkfile=job.chk_dir / f'{key}.chk', **KW)
        wall = rc.run_scf(mf, chk_guess=job.chk_dir / f'{key}.chk')
        rec = {'charge': charge, 'spin': spin, 'settings': rc.settings_of(mf), 'E_Ha': float(mf.e_tot),
               'scf_wall_s': round(wall, 1), 'cycles': mf._cycles, **rc.runtime_sizes(mf)}
        if spin:
            rec['populations'] = rc.fragment_populations(mf, frags)
            rec['frontier'] = rc.fragment_composition(mf, frags)
            ss, mult = mf.spin_square()
            rec['S2'] = float(ss)
        if disp:
            rec['E_disp_Ha'] = float(mf.scf_summary.get('dispersion', np.nan))
            t = time.time(); g = mf.nuc_grad_method().kernel()
            rec.update({'grad_wall_s': round(time.time() - t, 1), 'grad_max': float(abs(g).max())})
            mf0 = rc.make_mf(mol, xc=xc, solvent=solvent, **KW); rc.run_scf(mf0, dm0=mf.make_rdm1())
            g0 = mf0.nuc_grad_method().kernel()
            rec['grad_disp_max'] = float(abs(g - g0).max())       # D3(BJ) part of the gradient, must be > 0
            assert rec['grad_disp_max'] > 1e-6, 'dispersion gradient missing'
            assert rec['E_disp_Ha'] < 0, 'dispersion energy missing'
        if solvent is not None:
            rec['eps'] = mf.with_solvent.eps
            rec['E_solvation_Ha'] = float(mf.with_solvent.e)
        job.put(key, rec)
        extra = ''
        if spin:
            mp = rec['populations']['mulliken']
            extra = (f"; Mulliken charge/spin cation {mp['cation']['charge']:+.3f}/{mp['cation']['spin']:+.3f}, "
                     f"anion {mp['anion']['charge']:+.3f}/{mp['anion']['spin']:+.3f}; <S^2> {rec['S2']:.3f}")
        if disp:
            extra = f"; E_disp {rec['E_disp_Ha']:.6f} Ha, max |g_disp| {rec['grad_disp_max']:.2e}"
        job.log(f"{key}: {type(mf).__name__} xc={mf.xc} E = {mf.e_tot:.8f} ({wall:.0f} s){extra}")
    job.r['skipped_on_this_run'] = skipped
    job.save()
    return True


if __name__ == '__main__':
    rc.install_sigterm_handler()
    rc.main_wrapper(parser, run, __file__)
