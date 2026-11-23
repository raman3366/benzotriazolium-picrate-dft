#!/usr/bin/env python3
"""rev_struct_tables.py: the revision tables on structure (round 1, COMPTC-D-26-02710).

Run (from the repository root, after build_facts.py):
    python scripts/tables/rev_struct_tables.py

Writes, in the format of make_tables.py (one pandoc pipe table per file, caption line
'Table: <caption> {#tbl:<stem>}'):
    tables/tbl_rev_minima.md     main: hydrogen bond, inter-ion torsion, ring planes, lowest wavenumber and relative
                                 energy of the crystal pairs and the B3LYP / B3LYP-D3(BJ) minima of both crystal forms
    tables/tblS_rev_torsions.md  SI: torsion angles of the P2_1/c pair (Reviewer 1, point 5), with Klyne-Prelog ranges
    tables/tblS_rev_torscan.md   SI: relaxed gas-phase scans of C1-O1...N4-C12, B3LYP and B3LYP-D3(BJ)
    tables/tblS_rev_softmode.md  SI: soft-mode Hessian checks
Every number comes from papers/btpc/data/facts.json (keys rev.struct.* and the published gas.*, pcm.*, xray.*);
the script stops on a missing key or a wrong unit (make_tables.py helpers) and lints the output with the house rules.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_tables as MT                                  # noqa: E402  (helpers: val, fnum, fmt, wn, render, lint)
from make_tables import fmt, wn, MINUS, NDASH            # noqa: E402
HERE = MT.HERE                             # dataset: tables/ (make_tables.py, scripts/paths.py)

S = 'rev.struct.'
CS = '*P*2~1~/*c*'
CN = '*P*2~1~/*n*'
GAS, PCM = MT.MODEL['gas'], MT.MODEL['pcm']


def f(key, unit=None):
    return MT.fnum(key, unit)


def s(key):
    return MT.fstr(key, '')


def lst(key, unit):
    v = MT.val(key, unit)
    if not isinstance(v, list):
        MT.fail(f"key '{key}' is not a list")
    return v


def kp(t):
    """Klyne-Prelog range (IUPAC 1996): sp |t| < 30, sc 30-90, ac 90-150, ap >= 150; sign for sc and ac."""
    a = abs(t)
    if a < 30:
        return 'sp'
    if a >= 150:
        return 'ap'
    return (('+' if t > 0 else MINUS) + ('sc' if a < 90 else 'ac'))


def start(stem):
    MT.CTX.stem, MT.CTX.keys, MT.CTX.files, MT.CTX.checks = stem, [], [str(MT.FACTS_PATH.relative_to(MT.REPO))], []


def write(stem, where, title, spec):
    t = MT.Table(stem=stem, where=where, title=title, keys=list(MT.CTX.keys), files=list(MT.CTX.files),
                 checks=list(MT.CTX.checks), **spec)
    text = MT.render(t)
    MT.lint(text, stem)
    (HERE / f'{stem}.md').write_text(text, encoding='utf-8')
    print(f'{stem:<22} {len(t.rows):3d} rows  {len(t.keys):3d} keys')


# =================================================================================================
def t_minima():
    start('tbl_rev_minima')
    rows = []
    hb4 = 'hb_N4_H4N_O1'
    hz = 'hb_' + s(S + 'btpcz.defining_hbond').replace('-', '_').replace('...', '_')
    MT.check(hz == 'hb_N6_H6_O1', 'P2_1/n defining hydrogen bond is N6-H6...O1')

    def geo(pfx, hb, tau_key, ring_key):
        return [fmt(f(f'{pfx}{hb}_dh', 'Å'), 3), fmt(f(f'{pfx}{hb}_ha', 'Å'), 3), fmt(f(f'{pfx}{hb}_da', 'Å'), 3),
                fmt(f(f'{pfx}{hb}_angle', 'deg'), 1), fmt(f(tau_key, 'deg'), 1), fmt(f(ring_key, 'deg'), 1)]

    # --- P2_1/c pair
    rows.append([f'{CS} pair (CCDC 835950), N4–H4N···O1'] + [''] * 8)
    rows.append(['X-ray'] + geo('xray.', hb4, S + 'tors_xray_C1_O1_N4_C12', 'xray.interplanar_picrate_cation') + ['–', '–'])
    MT.check(f('gas.n_imag', '') == 0, 'published gas-phase minimum has no imaginary mode')
    rows.append([f'{GAS}, B3LYP'] + geo('gas.', hb4, S + 'tors_gas_C1_O1_N4_C12', 'gas.interplanar_picrate_cation')
                + [wn(f('gas.freq_lowest', 'cm-1')), 'ref.'])
    MT.check(s(S + 'softmode_d3_gas.assessment') == 'informative', 'gas D3(BJ) soft-mode check is informative')
    MT.check(f(S + 'softmode_d3_gas.n_imag_projected', '') == 0, 'gas D3(BJ) soft subspace has no negative curvature')
    rows.append([f'{GAS}, B3LYP-D3(BJ)'] + geo(S + 'd3_gas.', hb4, S + 'd3_gas.tors_C1_O1_N4_C12', S + 'd3_gas.interplanar_picrate_cation')
                + [wn(f(S + 'softmode_d3_gas.lowest_projected', 'cm-1')) + '^a^', '–'])
    MT.check(MT.from_value('pcm.opt_converged', r'^no') is not None, 'published IEF-PCM optimisation not converged')
    rows.append([f'{PCM}, B3LYP^b^'] + geo('pcm.', hb4, S + 'tors_pcm_C1_O1_N4_C12', 'pcm.interplanar_picrate_cation')
                + [wn(f('pcm.freq_lowest', 'cm-1')), 'ref.'])
    MT.check(s(S + 'softmode_d3_pcm.assessment') == 'not informative', 'PCM D3(BJ) soft-mode check is not informative')
    rows.append([f'{PCM}, B3LYP-D3(BJ)'] + geo(S + 'd3_pcm.', hb4, S + 'd3_pcm.tors_C1_O1_N4_C12', S + 'd3_pcm.interplanar_picrate_cation')
                + ['n.d.^c^', '–'])
    # --- P2_1/n pair
    rows.append([f'{CN} pair (COD 2230422), N6–H6···O1'] + [''] * 8)
    rows.append(['X-ray'] + geo(S + 'btpcz_xray.', hz, S + 'btpcz_xray.tors_C1_O1_N6_C12', S + 'btpcz_xray.interplanar_picrate_cation')
                + ['–', '–'])
    MT.check(s(S + 'btpcz_gas.proton_on') == 'O1', 'P2_1/n gas-phase structure has the proton on O1')
    rows.append([f'{GAS}, B3LYP^d^'] + geo(S + 'btpcz_gas.', hz, S + 'btpcz_gas.tors_C1_O1_N6_C12', S + 'btpcz_gas.interplanar_picrate_cation')
                + [wn(f(S + 'btpcz_gas.freq_lowest', 'cm-1')), fmt(f(S + 'btpcz_gas.e_minus_btpc', 'kcal mol-1'), 2)])
    MT.check(f(S + 'btpcz_pcm.n_imag', '') == 0 and s(S + 'btpcz_pcm.proton_on') == 'N6', 'P2_1/n IEF-PCM: ion pair, no imaginary mode')
    rows.append([f'{PCM}, B3LYP'] + geo(S + 'btpcz_pcm.', hz, S + 'btpcz_pcm.tors_C1_O1_N6_C12', S + 'btpcz_pcm.interplanar_picrate_cation')
                + [wn(f(S + 'btpcz_pcm.freq_lowest', 'cm-1')), fmt(f(S + 'btpcz_pcm.e_minus_btpc', 'kcal mol-1'), 2)])
    nm = int(f(S + 'softmode_d3_gas.n_modes', ''))
    caption = (
        'Hydrogen bond and relative orientation of the ions in the crystal pairs of the two crystal forms of '
        'benzotriazolium picrate and in the pairs optimised from them with B3LYP and B3LYP-D3(BJ)/6-311++G(d,p) in '
        f'the gas phase and in IEF-PCM (water). N is the cation nitrogen that donates to the phenolate O1: N4 in the '
        f'{CS} form and N6 in the {CN} form, the nitrogen bonded to C12 in both. *τ* is the pseudo-torsion '
        'C1–O1···N–C12 about the non-bonded O1···N axis, which measures the turn of one ion against the other about the '
        'hydrogen bond; a mirror-image structure has the opposite sign. Ring planes: angle between the least-squares '
        'planes of the six picrate ring carbons and of the nine ring atoms of the cation. Lowest ν̃: lowest harmonic '
        'wavenumber of the full finite-difference Hessian, unscaled (i, imaginary). Δ*E*: energy of the '
        f'{CN} structure relative to the {CS} structure (ref.) at the same level and in the same model; '
        'both pairs contain the same atoms.')
    notes = [
        f'^a^ Lowest eigenvalue of the Hessian projected on the {MT.words(nm)}-dimensional subspace of the softest modes '
        f'of the published gas-phase Hessian (soft-mode check in the Supplementary Information); no full Hessian '
        f'was computed.',
        f'^b^ Published structure, accepted on the gradient criteria only; its displacement and energy criteria were '
        f'not met.',
        f'^c^ Not determined. The soft-mode check used the Hessian of the published IEF-PCM structure, which lies '
        f'{fmt(f(S + "softmode_d3_pcm.ref_geometry_rmsd", "Å"), 2)} Å (best-fit RMSD) from this one, so the modes it '
        f'probed are not the soft modes of the new structure (projected wavenumbers '
        f'{wn(f(S + "softmode_d3_pcm.lowest_projected", "cm-1"), 0)}{NDASH}{wn(f(S + "softmode_d3_pcm.highest_projected", "cm-1"), 0)} cm⁻¹).',
        f'^d^ The proton moved to O1 during the optimisation, so this is the neutral pair: *d*(N–H) and '
        f'*d*(H···O1) are N6···H6 and O1–H6. Its last optimisation steps alternated between two SCF solutions '
        f'(mean |Δ*E*| per step {fmt(f(S + "btpcz_gas.opt_tail_mean_abs_dE", "kcal mol-1"), 2)} kcal mol⁻¹) and its '
        f'finite-difference Hessian is {fmt(f(S + "btpcz_gas.fd_hessian_asymmetry_pct", "%"), 1)} % asymmetric before '
        f'symmetrisation, against {fmt(f(S + "btpcz_pcm.fd_hessian_asymmetry_pct", "%"), 1)} % for the IEF-PCM '
        f'structure, so the imaginary mode is not a reliable verdict on this structure.']
    return dict(header=['Structure', '*d*(N–H) (Å)', '*d*(H···O1) (Å)', '*d*(N···O1) (Å)', '∠N–H···O1 (°)', '*τ* (°)',
                        'Ring planes (°)', 'Lowest ν̃ (cm⁻¹)', 'Δ*E* (kcal mol⁻¹)'],
                align='lrrrrrrrr', rows=rows, caption=caption, notes=notes)


# =================================================================================================
TORS_BLOCKS = [
    ('Inter-ion orientation (°)', [('C1', 'O1', 'N4', 'C12'), ('C6', 'C1', 'O1', 'H4N')]),
    ('Nitro groups, O–N–C–C (°)', [('O2', 'N1', 'C2', 'C1'), ('O3', 'N1', 'C2', 'C3'), ('O4', 'N2', 'C4', 'C3'),
                                   ('O5', 'N2', 'C4', 'C5'), ('O7', 'N3', 'C6', 'C1'), ('O6', 'N3', 'C6', 'C5')]),
    ('Picrate ring (°)', [('O1', 'C1', 'C2', 'N1'), ('O1', 'C1', 'C6', 'N3'), ('C6', 'C1', 'C2', 'C3')]),
    ('Cation (°)', [('C12', 'N4', 'N5', 'N6'), ('N4', 'N5', 'N6', 'C7'), ('H4N', 'N4', 'N5', 'N6'), ('H6N', 'N6', 'N5', 'N4')]),
]
TORS_MODELS = ['xray', 'gas', 'pcm', 'd3gas', 'd3pcm']
NONB = {('O1', 'N4'), ('N4', 'O1'), ('O1', 'H4N'), ('H4N', 'O1')}
ORTHO_PARA = {'N1': ' (*ortho*)', 'N2': ' (*para*)', 'N3': ' (*ortho*)'}


def tname(t):
    out = t[0]
    for a, b in zip(t, t[1:]):
        out += ('···' if (a, b) in NONB else NDASH) + b
    if t[1] in ORTHO_PARA:
        out += ORTHO_PARA[t[1]]
    return out


def t_torsions():
    start('tblS_rev_torsions')
    rows = []
    n_tors = 0
    for title, ts in TORS_BLOCKS:
        rows.append([f'*{title}*'] + [''] * 5)
        for t in ts:
            cells = []
            for mdl in TORS_MODELS:
                v = f(S + f"tors_{mdl}_{'_'.join(t)}", 'deg')
                cells.append(f'{fmt(v, 1)} ({kp(v)})')
            rows.append([tname(t)] + cells)
            n_tors += 1
    rows.append(['*Ring–ring relation*'] + [''] * 5)
    ring = {'xray': 'xray.interplanar_picrate_cation', 'gas': 'gas.interplanar_picrate_cation',
            'pcm': 'pcm.interplanar_picrate_cation', 'd3gas': S + 'd3_gas.interplanar_picrate_cation',
            'd3pcm': S + 'd3_pcm.interplanar_picrate_cation'}
    cent = {k: v.replace('interplanar_picrate_cation', 'centroid_distance_picrate_cation') for k, v in ring.items()}
    rows.append(['Angle between the ring planes (°)'] + [fmt(f(ring[m], 'deg'), 1) for m in TORS_MODELS])
    rows.append(['Ring centroid distance (Å)'] + [fmt(f(cent[m], 'Å'), 3) for m in TORS_MODELS])
    rows.append(['O1···N4 line to the picrate ring plane (°)'] + [fmt(f(S + f'line_O1_N4_vs_picrate_{m}', 'deg'), 1) for m in TORS_MODELS])
    rows.append(['O1···N4 line to the cation plane (°)'] + [fmt(f(S + f'line_O1_N4_vs_cation_{m}', 'deg'), 1) for m in TORS_MODELS])
    caption = (
        f'Torsion angles of the benzotriazolium picrate ion pair in the crystal ({CS}, 298 K), in the published '
        'B3LYP/6-311++G(d,p) structures and in the B3LYP-D3(BJ)/6-311++G(d,p) minima, gas phase and IEF-PCM (water), '
        f'with the relation between the two rings ({n_tors} torsions). The Klyne–Prelog range of each torsion is given in '
        'parentheses: sp, synperiplanar (|*τ*| < 30°); sc, synclinal (30–90°); ac, anticlinal (90–150°); ap, '
        'antiperiplanar (150–180°); the sign is given for sc and ac. C1–O1···N4–C12 and C6–C1–O1···H4N are '
        'pseudo-torsions through the non-bonded O1···N4 and O1···H4N contacts; the first measures the turn of the two '
        'ions about the N4–H4N···O1 bond. X-ray torsions were computed from the deposited coordinates of the modelled '
        'pair; the CIF contains no torsion list, so they carry no standard uncertainties, and those involving H4N or '
        'H6N rest on the short X-ray N–H positions. The published IEF-PCM structure was accepted on the gradient '
        'criteria only; the B3LYP-D3(BJ) IEF-PCM structure met all five criteria, but its Hessian was not '
        'checked (@tbl:tblS_rev_softmode). Ring planes are the least-squares planes of the six picrate ring carbons and of the nine '
        'ring atoms of the cation. Atom labels follow the CIF.')
    return dict(header=['Torsion or angle', 'X-ray', f'{GAS}, B3LYP', f'{PCM}, B3LYP', f'{GAS}, B3LYP-D3(BJ)',
                        f'{PCM}, B3LYP-D3(BJ)'],
                align='lrrrrr', rows=rows, caption=caption)


# =================================================================================================
def t_torscan():
    start('tblS_rev_torscan')
    B, D = S + 'scan_b3lyp.', S + 'scan_d3.'
    grid = lst(B + 'grid', 'deg')
    if lst(D + 'grid', 'deg') != grid:
        MT.fail('B3LYP and D3(BJ) grids differ')
    cols = [(B + 'dE', 'kcal mol-1', 3), (D + 'dE', 'kcal mol-1', 3), (B + 'd_N4_O1', 'Å', 3), (D + 'd_N4_O1', 'Å', 3),
            (B + 'd_N4_H4N', 'Å', 3), (D + 'd_N4_H4N', 'Å', 3), (B + 'interplanar', 'deg', 1), (D + 'interplanar', 'deg', 1)]
    data = [lst(k, u) for k, u, _ in cols]
    for k in (B, D):
        MT.check(f(k + 'n_converged', '') == f(k + 'n_points', '') == len(grid), f'{k}: every point converged')
        MT.check(f(k + 'n_proton_on_N4', '') == len(grid), f'{k}: proton on N4 at every point')
        tau = lst(k + 'tau', 'deg')
        MT.check(max(abs(a - b) for a, b in zip(tau, grid)) < 0.01, f'{k}: constrained torsions reached to 0.01°')
    rows = []
    for i, g in enumerate(grid):
        rows.append([fmt(g, 0)] + [fmt(d[i], nd) for d, (_, _, nd) in zip(data, cols)])
    rt = f(S + 'rt_298_kcal', 'kcal mol-1')
    caption = (
        'Relaxed scans of the pseudo-torsion C1–O1···N4–C12 of the ion pair in the gas phase with B3LYP and '
        'B3LYP-D3(BJ)/6-311++G(d,p): at each value every other coordinate was optimised to the Gaussian default '
        'criteria. Δ*E* is relative to the unconstrained minimum of the same functional (B3LYP, the published '
        f'gas-phase minimum at *τ* = {fmt(f(B + "tau_ref", "deg"), 1)}°; B3LYP-D3(BJ), the new minimum at '
        f'{fmt(f(D + "tau_ref", "deg"), 1)}°); *RT* at 298.15 K is {fmt(rt, 3)} kcal mol⁻¹. The X-ray torsion is '
        f'{fmt(f(S + "tors_xray_C1_O1_N4_C12", "deg"), 1)}°. All {MT.words(len(grid))} points of both scans converged and H4N '
        'stays on N4 at every point. Ring planes: angle between the least-squares planes of the picrate ring and of '
        'the cation.')
    return dict(header=['*τ* (°)', 'Δ*E*, B3LYP (kcal mol⁻¹)', 'Δ*E*, D3(BJ) (kcal mol⁻¹)', 'N4···O1, B3LYP (Å)',
                        'N4···O1, D3(BJ) (Å)', 'N4–H4N, B3LYP (Å)', 'N4–H4N, D3(BJ) (Å)', 'Ring planes, B3LYP (°)',
                        'Ring planes, D3(BJ) (°)'],
                align='rrrrrrrrr', rows=rows, caption=caption)


# =================================================================================================
def t_softmode():
    start('tblS_rev_softmode')
    checks = [('d3_gas', f'{CS}, {GAS}, B3LYP-D3(BJ)', f'the published B3LYP gas-phase Hessian'),
              ('d3_pcm', f'{CS}, {PCM}, B3LYP-D3(BJ)', f'the published B3LYP IEF-PCM Hessian'),
              ('btpcz_pcm', f'{CN}, {PCM}, B3LYP', f'the {CN} gas-phase Hessian')]
    rows = []
    summary = []
    nm0 = None
    for tag, title, refname in checks:
        P = S + f'softmode_{tag}.'
        nm = int(f(P + 'n_modes', ''))
        nm0 = nm0 or nm
        ref, dg, pr = lst(P + 'ref_freqs', 'cm-1'), lst(P + 'diagonal_freqs', 'cm-1'), lst(P + 'projected_freqs', 'cm-1')
        if not len(ref) == len(dg) == len(pr) == nm:
            MT.fail(f'{tag}: {len(ref)}/{len(dg)}/{len(pr)} values for {nm} modes')
        rows.append([title, '', '', ''])
        for k in range(nm):
            rows.append([str(k + 1), wn(ref[k]), wn(dg[k]), wn(pr[k])])
        rr = lst(P + 'residual_relative_range', '')
        summary.append(
            f'{title}: reference {refname}, structures '
            f'{fmt(f(P + "ref_geometry_rmsd", "Å"), 3)} Å apart; updated Hessian lowest '
            f'{wn(f(P + "updated_lowest", "cm-1"))} cm⁻¹ ({int(f(P + "n_imag_updated", ""))} imaginary); '
            f'asymmetry {fmt(f(P + "asymmetry_pct", "%"), 2)} %; share of *Hq*~*k*~ outside the subspace '
            f'{fmt(rr[0], 2)}–{fmt(rr[1], 2)}; {s(P + "assessment")}')
    step = f(S + 'softmode_d3_gas.step_bohr', 'bohr')
    caption = (
        f'Soft-mode Hessian checks. For each structure the {MT.words(nm0)} softest vibrational modes *q*~*k*~ of a '
        'reference Hessian were taken at the structure, and analytic gradients were computed at displacements of '
        f'±{fmt(step, 2)} Bohr (largest atomic displacement) along each; central differences give the columns *Hq*~*k*~. '
        'Reference ν̃: wavenumbers of the reference Hessian at the structure; along mode: curvature along each '
        'reference mode alone; projected: eigenvalues of the Hessian projected on the soft subspace (Rayleigh–Ritz '
        'values; a negative value would prove an imaginary frequency). Unscaled; i, imaginary. The check speaks for a '
        'structure only when the reference Hessian belongs to a nearby structure, so that its soft modes are still the '
        'soft modes there (rule used: best-fit RMSD ≤ 0.1 Å and highest projected value ≤ 1.5 × highest reference '
        'value). Summary per check: ' + '. '.join(summary) + '. The last two checks probe directions that are '
        'stiff at the structure checked and say nothing about its softest modes; the updated Hessians, which mix a '
        'Hessian of one structure with gradients of another, have spurious imaginary values. The '
        f'{CN} IEF-PCM minimum was confirmed by its full finite-difference Hessian instead (lowest '
        f'{wn(f(S + "btpcz_pcm.freq_lowest", "cm-1"))} cm⁻¹, none imaginary).')
    return dict(header=['Mode', 'Reference ν̃ (cm⁻¹)', 'Along mode (cm⁻¹)', 'Projected (cm⁻¹)'],
                align='lrrr', rows=rows, caption=caption)


def main():
    write('tbl_rev_minima', 'main', 'Hydrogen bond, inter-ion torsion and lowest mode of all minima, both crystal forms', t_minima())
    write('tblS_rev_torsions', 'SI', 'Torsion angles of the P2_1/c pair with Klyne-Prelog ranges (R1.5)', t_torsions())
    write('tblS_rev_torscan', 'SI', 'Relaxed C1-O1...N4-C12 scans, B3LYP and B3LYP-D3(BJ)', t_torscan())
    write('tblS_rev_softmode', 'SI', 'Soft-mode Hessian checks', t_softmode())


if __name__ == '__main__':
    try:
        main()
    except MT.TableError as e:
        print(f'ERROR {e}', file=sys.stderr)
        sys.exit(1)
