"""
Soft-mode Hessian check at an optimised geometry (plan section 5, item 4): is it a true minimum along the floppiest
directions, without a full Hessian?

The N softest vibrational normal modes q_k of a reference Hessian (default: the stored finite-difference Hessian of
the same model, data/r2/molecules/btpc/results/btpc_<model>/hessian.npy) are taken at the target geometry (overall
translations and rotations projected out, mass-weighted with average isotope masses).  For each mode the analytic
gradient is computed at x0 +/- h d_k (d_k = M^-1/2 q_k scaled so that the largest atomic displacement is h, default
0.02 Bohr): 2N gradients, plus one at x0 for the residual gradient.  Central differences give the mass-weighted
Hessian columns H q_k, from which
  * the projected N x N Hessian K = Q^T H Q (symmetrised; its antisymmetric part measures the numerical noise) and its
    frequencies -- the Rayleigh-Ritz values of the true Hessian in the soft subspace;
  * the residual ||(1 - Q Q^T) H q_k|| (coupling of each mode to the rest of the spectrum);
  * the reference Hessian updated with the new columns (symmetric rank-2N update), whose lowest frequencies estimate
    the full spectrum including that coupling.
A negative eigenvalue of K proves a negative curvature (an imaginary frequency); all-positive values with small
residuals show that the soft subspace, which holds the modes that were imaginary or near zero before, is now stable.

Gradients: the level of theory of the optimisation (--xc, --disp, --model), conv_tol 1e-10, grid 4, Lebedev 29,
PySCF's default gradients without grid-weight response (as dft/frequencies.py); --grid-response adds it.
Usage:
  python dft/revision/softmode_check.py --job soft_pcm --out <dir> --geom <xyz or opt_d3 output dir> --model pcm
         [--disp d3bj|none] [--ref-hessian auto|PATH.npy] [--nmodes 10] [--step 0.02] [--case btpc] [--smoke]
--smoke: STO-3G, grid 2, Lebedev 17, conv_tol 1e-7, 2 modes (5 gradients).
"""
import sys, pathlib, argparse, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import rev_common as rc          # first: pins the thread pools
import numpy as np
import headline as hl
from opt_d3 import tr_basis, complement, hessian_matrix

AMU_ME = 1822.888486209          # atomic mass unit in electron masses
AU_CM = 219474.6313632           # hartree -> cm-1


def freqs_cm(evals):
    ev = np.asarray(evals, float)
    return np.sign(ev) * np.sqrt(np.abs(ev)) * AU_CM


def vib_modes(H_cart, X_bohr, masses_amu):
    """Mass-weighted vibrational modes of a Cartesian Hessian at X (TR projected): (eigenvalues a.u., modes 3N x m,
    TR basis)."""
    m = np.repeat(np.asarray(masses_amu, float) * AMU_ME, 3)
    Hm = H_cart / np.sqrt(np.outer(m, m))
    T = tr_basis(np.asarray(X_bohr), masses_amu)
    N = complement(T)
    w, V = np.linalg.eigh(N.T @ Hm @ N)
    return w, N @ V, T, m


def parser():
    p = rc.add_common_args(argparse.ArgumentParser(description=__doc__.split('\n')[1]))
    p.add_argument('--geom', required=True, help='xyz file, an opt_d3.py output directory, or stored-gas / stored-pcm')
    p.add_argument('--model', choices=('gas', 'pcm'), required=True)
    p.add_argument('--case', default='btpc')
    p.add_argument('--xc', default=rc.XC)
    p.add_argument('--disp', default='d3bj', help='d3bj (default, for the D3(BJ) minima) or none')
    p.add_argument('--ref-hessian', default='auto')
    p.add_argument('--nmodes', type=int, default=None)
    p.add_argument('--step', type=float, default=0.02, help='largest atomic displacement, Bohr (default 0.02)')
    p.add_argument('--grid-response', action='store_true')
    return p


def run(job, args):
    lv = hl.level(args.smoke)
    disp = None if args.disp in (None, '', 'none') else args.disp
    solvent = 'water' if args.model == 'pcm' else None
    nmodes = args.nmodes or (2 if args.smoke else 10)
    case = rc.load_case(f'{args.case}_{args.model}')
    if args.geom.startswith('stored-'):
        gpath = hl.resolve_geom(args.geom.split('-', 1)[1], args.case)[1]
    else:
        gpath = pathlib.Path(args.geom)
        gpath = gpath if gpath.is_absolute() else rc.ROOT / gpath
        if gpath.is_dir():
            gpath = gpath / 'opt.xyz'
    job.input_file('geometry', gpath)
    hp = (hl.stored_dir(f'{args.case}_{args.model}') / 'hessian.npy' if args.ref_hessian == 'auto'
          else pathlib.Path(args.ref_hessian))
    job.input_file('ref_hessian', hp)
    job.r['settings'] = {**lv, 'xc': args.xc, 'disp': disp, 'model': args.model, 'nmodes': nmodes,
                         'step_max_atom_bohr': args.step, 'grid_response': args.grid_response,
                         'masses': 'average isotope masses (PySCF atom_mass_list(isotope_avg=True))',
                         'differences': 'central, analytic gradients', 'pcm_method': rc.PCM_METHOD if solvent else None}
    job.save()

    mol = rc.mol_for_case(case, gpath, basis=lv['basis'], verbose=3, output=job.out / 'pyscf.log')
    X0 = mol.atom_coords()                       # Bohr
    masses = mol.atom_mass_list(isotope_avg=True)
    H_ref = hessian_matrix(hp)
    if H_ref.shape != (3 * mol.natm,) * 2:
        raise ValueError(f'reference Hessian shape {H_ref.shape} does not match {mol.natm} atoms')
    w_ref, modes, T, mvec = vib_modes(H_ref, X0, masses)
    Q = modes[:, :nmodes]                        # mass-weighted unit vectors, softest first
    ref_freqs = freqs_cm(w_ref[:max(nmodes, 12)])
    # reference frequencies at the reference geometry (stored opt.xyz of the same model), for the record
    rec_ref = {'path': str(hp), 'lowest_freqs_at_target_geometry_cm-1': ref_freqs.tolist()}
    try:
        Xr = rc.read_xyz(hl.resolve_geom(args.model, args.case)[1])[1] / rc.BOHR_A
        rec_ref['lowest_freqs_at_reference_geometry_cm-1'] = freqs_cm(vib_modes(H_ref, Xr, masses)[0][:12]).tolist()
        rec_ref['rms_geometry_difference_A'] = float(np.sqrt(((Xr - X0) ** 2).sum(1).mean()) * rc.BOHR_A)
    except Exception as ex:
        rec_ref['note'] = f'reference geometry not compared: {ex}'
    job.r['reference'] = rec_ref
    job.log(f'reference Hessian {hp}: softest {nmodes} modes at the target geometry: {np.round(ref_freqs[:nmodes], 1)} cm-1')

    chk0 = job.chk_dir / 'x0.chk'

    def gradient(key, X):
        if job.done(key):
            return np.array(job.get(key)['gradient'])
        m2 = mol.copy(); m2.set_geom_(X, unit='Bohr')
        mf = rc.make_mf(m2, xc=args.xc, disp=disp, solvent=solvent, conv_tol=lv['conv_tol'],
                        grid_level=lv['grid_level'], lebedev_order=lv['lebedev_order'],
                        chkfile=chk0 if key == 'x0' else None)
        t = time.time()
        rc.run_scf(mf, chk_guess=chk0 if chk0.exists() else None)
        gm = mf.nuc_grad_method()
        if args.grid_response:
            gm.grid_response = True
        g = gm.kernel()
        wall = time.time() - t
        job.put(key, {'E_Ha': float(mf.e_tot), 'gradient': g.tolist(), 'scf_cycles': mf._cycles,
                      'wall_s': round(wall, 1), 'net_force_Ha_bohr': float(np.linalg.norm(g.sum(0)))})
        job.r['timings'][key] = round(wall, 1)
        job.log(f'{key}: E = {mf.e_tot:.10f} Ha, |sum g| {np.linalg.norm(g.sum(0)):.1e}, {wall:.0f} s '
                f'({mf._cycles} SCF cycles)')
        return g

    g0 = gradient('x0', X0)
    res0 = hl.geometry_metrics(case, *hl.mol_geometry(mol))
    ga = np.linalg.norm(g0, axis=1)
    job.r['residual_gradient'] = {'rms_atom_norm': float(np.sqrt((ga ** 2).mean())), 'max_atom_norm': float(ga.max()),
                                  'projected_on_soft_modes': (Q.T @ (g0.ravel() / np.sqrt(mvec))).tolist()}
    cols = np.zeros((3 * mol.natm, nmodes)); scales = []
    for k in range(nmodes):
        d = (Q[:, k] / np.sqrt(mvec)).reshape(-1, 3)          # Cartesian direction of mode k (Bohr per unit q)
        s = args.step / np.linalg.norm(d, axis=1).max()
        scales.append(float(s))
        gp = gradient(f'mode{k + 1}+', X0 + s * d)
        gm = gradient(f'mode{k + 1}-', X0 - s * d)
        cols[:, k] = ((gp - gm).ravel() / (2 * s)) / np.sqrt(mvec)  # = H_mw q_k
    # projected Hessian and diagnostics
    P_vib = np.eye(3 * mol.natm) - T @ T.T
    cols = P_vib @ cols
    K = Q.T @ cols
    asym = float(np.linalg.norm(K - K.T) / max(np.linalg.norm(K), 1e-30))
    Ks = 0.5 * (K + K.T)
    wK, UK = np.linalg.eigh(Ks)
    resid = np.linalg.norm(cols - Q @ K, axis=0)
    # reference Hessian (mass-weighted, TR projected) updated with the new columns
    Hm_ref = P_vib @ (H_ref / np.sqrt(np.outer(mvec, mvec))) @ P_vib
    D = cols - Hm_ref @ Q
    Hm_new = Hm_ref + D @ Q.T + Q @ D.T - Q @ (Q.T @ D) @ Q.T
    Hm_new = 0.5 * (Hm_new + Hm_new.T)
    N = complement(T)
    w_new = np.linalg.eigvalsh(N.T @ Hm_new @ N)
    out = {'projected_freqs_cm-1': freqs_cm(wK).tolist(),
           'projected_eigvec_on_ref_modes': UK.tolist(),
           'diagonal_freqs_cm-1': freqs_cm(np.diag(Ks)).tolist(),
           'reference_freqs_cm-1': freqs_cm(w_ref[:nmodes]).tolist(),
           'residual_norm_outside_subspace': resid.tolist(),
           'residual_relative': (resid / np.maximum(np.linalg.norm(cols, axis=0), 1e-30)).tolist(),
           'asymmetry_rel': asym,
           'updated_full_lowest_freqs_cm-1': freqs_cm(w_new[:max(12, nmodes)]).tolist(),
           'displacement_scales': scales, 'n_imaginary_projected': int((wK < 0).sum()),
           'n_imaginary_updated_full': int((w_new < 0).sum()), 'geometry': res0}
    out['verdict'] = ('no imaginary frequency in the soft subspace' if out['n_imaginary_projected'] == 0
                      else f"{out['n_imaginary_projected']} negative curvature(s) in the soft subspace: NOT a minimum")
    job.r['summary'] = out
    job.log(f"projected frequencies (cm-1): {np.round(out['projected_freqs_cm-1'], 1)}; reference "
            f"{np.round(out['reference_freqs_cm-1'], 1)}; updated full lowest {np.round(out['updated_full_lowest_freqs_cm-1'][:6], 1)};"
            f" asymmetry {asym:.2e}; {out['verdict']}")
    job.save()
    mol.stdout.close()
    return True          # the verdict is a result, not a failure; it is in result.json summary.verdict


if __name__ == '__main__':
    rc.install_sigterm_handler()
    rc.main_wrapper(parser, run, __file__)
