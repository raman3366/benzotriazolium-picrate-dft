#!/usr/bin/env python3
"""rev_proton_tables.py: the proton-position tables of the revision (COMPTC-D-26-02710, round 1).

Run (from the repository root, after papers/btpc/data/build_facts.py):
    python scripts/tables/rev_proton_tables.py

Writes papers/btpc/tables/<stem>.md for
    tbl_rev_proton   (main)  proton position by level: gas-phase scan + DVR columns, IEF-PCM E(neutral) - E(salt)
    tblS_rev_cc      (SI)    FNO-CCSD(T) correction: every basis / FNO size, salt and neutral, extrapolations
    tblS_rev_scan    (SI)    the 16 points of the relaxed gas-phase scan with the energies of every level
Every number comes from facts.json (rev.proton.* keys, papers/btpc/data/rev/proton.py); the script stops on a missing
key or a wrong unit (the helpers of make_tables.py), runs the cross-checks listed in each function and applies the
house-style lint of make_tables.py. make_tables.py itself, its INDEX.md and the other tables are not touched.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_tables as mt                     # noqa: E402  (fail-loud facts access, render, lint)
from make_tables import Table, CTX, fnum, fstr, val, fmt, check, close, fail, MINUS, NDASH   # noqa: E402
HERE = mt.HERE                             # dataset: tables/ (make_tables.py, scripts/paths.py)

P = 'rev.proton.'
KC = 'kcal mol-1'
LEVELS = [('b3lyp', 'B3LYP'), ('wb97xd3bj', 'ωB97X-D3(BJ)'), ('m062x', 'M06-2X'), ('pbe0', 'PBE0'),
          ('camb3lyp', 'CAM-B3LYP'), ('dfmp2', 'DF-MP2'), ('scsmp2', 'SCS-MP2')]


def lst(key, unit=None, n=None):
    v = val(P + key, unit)
    if not isinstance(v, list) or (n is not None and len(v) != n):
        fail(f"key '{P + key}' is not a list of {n} values: {v!r}")
    return [float(x) for x in v]


def f(key, unit=None):
    return fnum(P + key, unit)


def s(key):
    return fstr(P + key)


def sg(x, nd=2):
    return fmt(x, nd, sign=True)


def rng(a, b, nd=2):
    return f'{fmt(a, nd)} to {fmt(b, nd)}'


# =================================================================================================
def t_rev_proton():
    n = int(f('scan.n_points', ''))
    rHO = lst('scan.r_HO_A', 'Å', n)
    i_min = rHO.index(min(rHO, key=lambda r: abs(r - f('gas.ionic_point_r_HO_A', 'Å'))))
    i_low = rHO.index(min(rHO, key=lambda r: abs(r - f('gas.lowest_point_r_HO_A', 'Å'))))
    rows, zpe_cm = [], []
    for k, name in LEVELS:
        dE = lst(f'gas.{k}.dE_kcal', KC, n)
        close(dE[i_min], 0.0, 1e-12, f'{k}: profile zero at the ionic minimum')
        check(min(range(n), key=lambda i: dE[i]) == i_low, f'{k}: lowest profile point = rev.proton.gas.lowest_point')
        de = f(f'gas.{k}.dE_ionic_minus_neutral_kcal', KC)
        close(de, -dE[i_low], 1e-9, f'{k}: dE(ionic - neutral) = -dE_kcal[lowest]')
        well = s(f'gas.{k}.ionic_well')
        if well == 'yes':
            depth = fmt(f(f'gas.{k}.well_depth_kcal', KC), 2)
            check(f(f'dvr.{k}.ionic_min_above_v0_kcal', KC) > 0, f'{k}: ionic minimum above the DVR ground level')
        elif well == 'no':
            depth = 'none'
        else:
            fail(f'{k}: ionic_well = {well!r}')
        zpe_cm.append(f(f'dvr.{k}.zpe_cm', 'cm-1'))
        rows.append([name, fmt(de, 2), depth, fmt(f(f'dvr.{k}.zpe_kcal', KC), 2), fmt(f(f'dvr.{k}.p_O_v0_pct', '%'), 1),
                     sg(f(f'pcm.{k}.dE_neutral_minus_salt_kcal', KC))])
    best = f('cc.best_kcal', KC)
    lo, hi = lst('cc.best_range_kcal', KC, 2)
    check(lo <= best <= hi, 'CCSD(T) estimate inside its range')
    close(best, f('pcm.dfmp2.dE_neutral_minus_salt_kcal', KC) + f('cc.pcm_631gd_f080.dd_kcal', KC), 1e-9,
          'best = DF-MP2 + 80 % correction')
    close(lo, f('cc.extrap_2pt_composite_kcal', KC), 1e-9, 'range low end = two-point extrapolation')
    close(hi, f('cc.extrap_3pt_composite_kcal', KC), 1e-9, 'range high end = three-point extrapolation')
    rows.append(['DF-MP2 + Δ[CCSD(T)]^a^', NDASH, NDASH, NDASH, NDASH, f'{sg(best)} ({rng(lo, hi)})'])
    r_lo, r_hi = lst('scan.r_HO_range_A', 'Å', 2)
    nk, nt = f('cc.pcm_631gd_f080.salt.nvir_kept', ''), f('cc.pcm_631gd_f080.salt.nvir_total', '')
    check(nk == f('cc.pcm_631gd_f080.neutral.nvir_kept', '') and nt == f('cc.pcm_631gd_f080.neutral.nvir_total', ''),
          'same FNO space for salt and neutral')
    check(f('pcm.n_dft_salt_lower', '') == 5, 'all five functionals keep the salt lower in IEF-PCM')
    neu = f('neutral_start_gas.dE_vs_salt_kcal', KC)
    check(neu < 0, 'free gas-phase neutral minimum below the salt')
    notes = [
        f'Gas phase: relaxed B3LYP scan of *r*(H4N···O1) from {fmt(r_hi, 2)} to {fmt(r_lo, 2)} Å ({n} points); the other '
        'levels are single points at the scan geometries. Δ*E*(ionic − neutral) is the energy of the published ionic '
        f'minimum (*r*(H4N···O1) = {fmt(rHO[i_min], 3)} Å) relative to the lowest point of the profile '
        f'({fmt(rHO[i_low], 2)} Å, proton on O1), which is the same point at every level; positive values mean that the '
        f'neutral form is lower. The fully optimised B3LYP neutral minimum lies {fmt(abs(neu), 2)} kcal mol⁻¹ below the '
        'ionic one. The ionic well is the rise from the lowest point with the proton on N4 to the highest point on the '
        'way to O1; "none" means that the energy falls all the way from the ionic geometry to the neutral one. '
        'ZPE is the ground level of the proton (one-dimensional sinc-DVR along its path on the scan, mass of H) above '
        f'the potential minimum ({fmt(min(zpe_cm), 0)}–{fmt(max(zpe_cm), 0)} cm⁻¹), and *P*(O side) is its probability '
        'with *r*(N4–H4N) > *r*(H4N···O1).',
        'IEF-PCM (water): energy of the neutral pair (proton on O1) minus that of the salt, both at their B3LYP '
        'IEF-PCM minima; positive values mean that the salt is lower. DF-MP2 and SCS-MP2 use an IEF-PCM HF reference '
        'with the reaction field frozen at the HF density.',
        '^a^ DF-MP2 plus the CCSD(T) − MP2 difference from frozen-natural-orbital CCSD(T)/6-31G* with '
        f'{fmt(nk, 0)} of {fmt(nt, 0)} virtual natural orbitals; in parentheses the range from the two extrapolations to '
        'the full virtual space (Supplementary Information).',
    ]
    return dict(header=['Level', 'Gas phase: Δ*E*(ionic − neutral) (kcal mol⁻¹)',
                        'Gas phase: ionic-well depth (kcal mol⁻¹)', 'Gas phase: proton ZPE (kcal mol⁻¹)',
                        'Gas phase: *P*(O side) (%)', 'IEF-PCM (water): *E*(neutral) − *E*(salt) (kcal mol⁻¹)'],
                align='lrrrrr', rows=rows, notes=notes,
                caption='Proton position in the benzotriazolium picrate pair by level of theory, '
                        'all with the 6-311++G(d,p) basis set at B3LYP/6-311++G(d,p) geometries.')


# =================================================================================================
CC_ROWS = [('pcm_631gd_f050', 'IEF-PCM (water), 6-31G*', 0), ('pcm_631gd_f065', 'IEF-PCM (water), 6-31G*', 1),
           ('pcm_631gd_f080', 'IEF-PCM (water), 6-31G*', 2), ('pcm_631ppgd_f050', 'IEF-PCM (water), 6-31+G(d,p)', 0),
           ('gas_631gd_f050', 'Gas phase, 6-31G*', 0)]


def t_rev_cc():
    fno = lst('cc.fno_pct', '%', 3)
    rows = []
    for run, label, ifno in CC_ROWS:
        pre = f'cc.{run}'
        nv = {}
        for k, sp in enumerate(('salt', 'neutral')):
            q = f'{pre}.{sp}'
            d = [f(f'{q}.delta_ccsd_kcal', KC), f(f'{q}.delta_T_kcal', KC), f(f'{q}.delta_kcal', KC),
                 f(f'{q}.delta_Tscaled_kcal', KC)]
            close(d[0] + d[1], d[2], 1e-9, f'{q}: delta = delta_CCSD + delta_(T)')
            nv[sp] = (f(f'{q}.nvir_kept', ''), f(f'{q}.nvir_total', ''))
            rows.append([f'{label}, {fmt(fno[ifno], 0)} %' if k == 0 else '', sp.capitalize(),
                         f'{fmt(nv[sp][0], 0)} of {fmt(nv[sp][1], 0)}', fmt(f(f'{q}.no_occ_pct', '%'), 2),
                         sg(d[0]), sg(d[1]), sg(d[2]), sg(d[3]), fmt(f(f'{q}.T1', ''), 3), fmt(f(f'{q}.D1', ''), 3), ''])
        check(nv['salt'][0] == nv['neutral'][0], f'{run}: same number of kept NOs')
        check(abs(nv['salt'][0] / nv['salt'][1] * 100 - fno[ifno]) < 0.5, f'{run}: kept NOs = FNO share of the virtuals')
        dd = [f(f'{pre}.dd_ccsd_kcal', KC), f(f'{pre}.dd_T_kcal', KC), f(f'{pre}.dd_kcal', KC), f(f'{pre}.dd_Tscaled_kcal', KC)]
        for j, nm in enumerate(('delta_ccsd_kcal', 'delta_T_kcal', 'delta_kcal', 'delta_Tscaled_kcal')):
            close(dd[j], f(f'{pre}.neutral.{nm}', KC) - f(f'{pre}.salt.{nm}', KC), 1e-9, f'{run}: Δδ = neutral - salt ({nm})')
        ref = f('pcm.dfmp2.dE_neutral_minus_salt_kcal' if run.startswith('pcm') else
                'gas_pcmgeom.dfmp2.dE_neutral_minus_salt_kcal', KC)
        comp = f(f'{pre}.composite_kcal', KC)
        close(comp, ref + dd[2], 1e-9, f'{run}: composite = DF-MP2 + Δδ')
        rows.append(['', 'Neutral − salt', '', '', sg(dd[0]), sg(dd[1]), sg(dd[2]), sg(dd[3]), '', '', sg(comp)])
    for nm, which in (('2pt', f'{fmt(fno[1], 0)} and {fmt(fno[2], 0)} %'), ('3pt', f'{fmt(fno[0], 0)}, {fmt(fno[1], 0)} and {fmt(fno[2], 0)} %')):
        e = f(f'cc.extrap_{nm}_kcal', KC)
        c = f(f'cc.extrap_{nm}_composite_kcal', KC)
        close(c, f('pcm.dfmp2.dE_neutral_minus_salt_kcal', KC) + e, 1e-9, f'extrapolated composite ({nm})')
        rows.append([f'IEF-PCM (water), 6-31G*, extrapolated from {which}', 'Neutral − salt', '', '', '', '', sg(e),
                     sg(f(f'cc.extrap_{nm}_Tscaled_kcal', KC)), '', '', sg(c)])
    t1, d1, d2 = lst('cc.T1_range', '', 2), lst('cc.D1_range', '', 2), lst('cc.D2_range', '', 2)
    nt_s, nt_n = f('cc.pcm_631ppgd_f050.salt.nvir_total', ''), f('cc.pcm_631ppgd_f050.neutral.nvir_total', '')
    check(nt_s - nt_n == 1, '6-31+G(d,p): one more virtual orbital for the salt (linear dependence at the neutral)')
    notes = [
        'δ = *E*~corr~[CCSD(T)] − *E*~corr~[MP2], both in the same space of frozen natural orbitals (FNO: virtual '
        'natural orbitals of the MP2 density, kept in order of occupation) with the same density-fitting integrals and '
        'frozen 1s cores; δ~CCSD~ and δ~(T)~ are its CCSD and (T) parts, and the scaled (T) is multiplied by the ratio '
        'of the full-space to the FNO-space MP2 correlation energy. NO occupation: share of the virtual natural-orbital '
        'occupation that is kept. Reference: density-fitted RHF; in IEF-PCM (water) the reaction field is frozen at the '
        'HF density for MP2 and CCSD(T), as for the DF-MP2 values. Geometries: the B3LYP/6-311++G(d,p) IEF-PCM '
        'minima of the salt and of the neutral pair, also for the gas-phase rows.',
        'Neutral − salt rows give Δδ = δ(neutral) − δ(salt); positive values move the balance toward the salt. Last '
        'column: DF-MP2/6-311++G(d,p) *E*(neutral) − *E*(salt) '
        f"({sg(f('pcm.dfmp2.dE_neutral_minus_salt_kcal', KC))} kcal mol⁻¹ in IEF-PCM (water), "
        f"{sg(f('gas_pcmgeom.dfmp2.dE_neutral_minus_salt_kcal', KC))} kcal mol⁻¹ in the gas phase at the same geometries) "
        'plus Δδ. Extrapolated rows: Δδ continued linearly in 1 − (kept NO occupation fraction) to the full virtual '
        'space, through the two smaller truncations or by least squares through all three.',
        f'*T*~1~ ranges from {fmt(t1[0], 3)} to {fmt(t1[1], 3)}, *D*~1~ from {fmt(d1[0], 3)} to {fmt(d1[1], 3)} and '
        f'*D*~2~ from {fmt(d2[0], 3)} to {fmt(d2[1], 3)}, alike for salt and neutral pair. In 6-31+G(d,p) one '
        f'near-linear-dependent combination of basis functions was removed at the neutral geometry ({fmt(nt_n, 0)} '
        f'instead of {fmt(nt_s, 0)} virtual orbitals); both pairs keep the same number of natural orbitals.',
    ]
    return dict(header=['Model, basis, FNO size', 'Pair', 'Virtual NOs kept', 'NO occupation (%)',
                        'δ~CCSD~ (kcal mol⁻¹)', 'δ~(T)~ (kcal mol⁻¹)', 'δ (kcal mol⁻¹)', 'δ, scaled (T) (kcal mol⁻¹)',
                        '*T*~1~', '*D*~1~', '*E*(neutral) − *E*(salt) (kcal mol⁻¹)'],
                align='llrrrrrrrrr', rows=rows, notes=notes,
                caption='Coupled-cluster correction to the DF-MP2 energy of the neutral pair relative to the salt '
                        'of benzotriazolium picrate.')


# =================================================================================================
def t_rev_scan():
    n = int(f('scan.n_points', ''))
    labels = val(P + 'scan.labels', '')
    rHO, rNH, rNO = lst('scan.r_HO_A', 'Å', n), lst('scan.r_NH_A', 'Å', n), lst('scan.r_NO_A', 'Å', n)
    delta = lst('scan.delta_A', 'Å', n)
    lumo = lst('scan.lumo_on_picrate_pct', '%', n)
    dE = {k: lst(f'gas.{k}.dE_kcal', KC, n) for k, _ in LEVELS}
    start = {'up': 'salt', 'down': 'salt', 'min': 'published minimum', 'neutral': 'neutral'}
    rows = []
    for i in range(n):
        close(delta[i], rNH[i] - rHO[i], 1e-9, f'point {labels[i]}: delta = r(N-H) - r(H...O)')
        if i:
            check(rHO[i] < rHO[i - 1], 'points ordered from the N side to the O side')
        rows.append([fmt(rHO[i], 3), fmt(rNH[i], 3), fmt(rNO[i], 3), start[labels[i].split('_')[0]]]
                    + [fmt(dE[k][i], 2) for k, _ in LEVELS] + [fmt(lumo[i], 1)])
    t_r = f('scan.transfer_r_HO_A', 'Å')
    hyst = f('scan.hysteresis_max_kcal', KC)
    dev = f('scan.max_constraint_dev_A', 'Å')
    check(dev < 1e-6, 'constraints held to better than 1e-6 Å')
    g5 = f('gas.m062x.grid5_rel_kcal', KC)
    notes = [
        'Energies in kcal mol⁻¹ relative to the published ionic minimum (the row marked published minimum) at each '
        'level; negative values lie below it. Every point is a B3LYP optimisation with *r*(H4N···O1) held fixed and '
        'all other coordinates free, started either from the salt or from a structure with the proton on O1 '
        '(started from); where both starts reached the same distance, the lower one is listed, and the two agreed to '
        f'within {fmt(hyst, 3)} kcal mol⁻¹. The other levels are single points at these geometries with the '
        '6-311++G(d,p) basis set (DF-MP2 and SCS-MP2 from the same calculation). '
        f'*r*(N4–H4N) = *r*(H4N···O1) at *r*(H4N···O1) = {fmt(t_r, 3)} Å (linear interpolation). Last column: share of '
        'the B3LYP LUMO on the picrate atoms (Mulliken, H4N counted with the cation). A finer M06-2X integration grid '
        f'(level 5 instead of 4) changes the energy difference of the two shortest points by {fmt(g5, 3)} kcal mol⁻¹ '
        f"and no checked absolute energy by more than {fmt(f('gas.m062x.grid5_abs_max_kcal', KC), 3)} kcal mol⁻¹.",
    ]
    return dict(header=['*r*(H4N···O1) (Å)', '*r*(N4–H4N) (Å)', '*r*(N4···O1) (Å)', 'Started from']
                + [name for _, name in LEVELS] + ['LUMO on picrate (%)'],
                align='rrrl' + 'r' * len(LEVELS) + 'r', rows=rows, notes=notes,
                caption='Relaxed scan of the proton coordinate of the benzotriazolium picrate pair in the gas phase, '
                        'B3LYP/6-311++G(d,p) geometries, with the energies of every level.')


TABLES = [('tbl_rev_proton', 'main', 'Proton position by level: gas-phase scan + DVR, IEF-PCM salt vs neutral', t_rev_proton),
          ('tblS_rev_cc', 'SI', 'FNO-CCSD(T) correction to E(neutral) - E(salt)', t_rev_cc),
          ('tblS_rev_scan', 'SI', 'Relaxed gas-phase proton scan, energies of every level', t_rev_scan)]


def main():
    for stem, where, title, fn in TABLES:
        CTX.stem, CTX.keys, CTX.files, CTX.checks = stem, [], [str(mt.FACTS_PATH.relative_to(mt.REPO))], []
        spec = fn()
        t = Table(stem=stem, where=where, title=title, keys=list(CTX.keys), files=list(CTX.files),
                  checks=list(CTX.checks), **spec)
        text = mt.render(t)
        mt.lint(text, stem)
        (HERE / f'{stem}.md').write_text(text, encoding='utf-8')
        print(f'{stem:<16} {where:<4} {len(t.rows):3d} rows  {len(t.keys):3d} keys  {len(t.checks):3d} checks')


if __name__ == '__main__':
    try:
        main()
    except mt.TableError as e:
        print(f'ERROR {e}', file=sys.stderr)
        sys.exit(1)
