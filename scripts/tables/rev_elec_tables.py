#!/usr/bin/env python3
"""rev_elec_tables.py: tables of the revision's electronic-structure results (BTPC, COMPTC-D-26-02710, round 1).

Run from the repository root, after papers/btpc/data/build_facts.py has written the rev.elec.* keys:
    python scripts/tables/rev_elec_tables.py

Writes (one pandoc pipe table per file, caption line 'Table: <caption> {#tbl:<stem>}'):
    tbl_rev_eps.md           main  frontier orbitals versus the permittivity of the continuum, both minima, + Delta-SCF
    tbl_rev_tddft.md         main  TDA excited states of the pair (Gas phase, IEF-PCM) and of the free ions
    tblS_rev_charges.md      SI    Mulliken, Loewdin, Hirshfeld and CM5 fragment charges at six points
    tblS_rev_functionals.md  SI    B3LYP, PBE0 and CAM-B3LYP frontier orbitals at eps = 1 and 78.3553
Every number comes from papers/btpc/data/facts.json (rev.elec.* keys from papers/btpc/data/rev/electronic.py); the
script stops on a missing key or a wrong unit.  It reuses make_tables.py's fail-loud facts access, number formatting,
pipe-table writer and style lint; it does not touch make_tables.py, its tables or INDEX.md.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_tables as mt                                                    # noqa: E402
from make_tables import (Table, TableError, CTX, MODEL, MINUS, NDASH, fnum, fstr, fint, fmt, render, lint,  # noqa: E402
                         check, close, method, FACTS, FACTS_PATH, REPO)
HERE = mt.HERE                             # dataset: tables/ (make_tables.py, scripts/paths.py)

EN = NDASH                                     # empty cell
XC = {'b3lyp': 'B3LYP', 'pbe0': 'PBE0', 'camb3lyp': 'CAM-B3LYP'}


def etag(eps):
    return 'e' + f'{eps:g}'.replace('.', 'p')


def eps_str(e):
    return f'{e:g}'


def series_eps(series):
    """eps values of a permittivity series, read from the keys present in facts.json."""
    pre = f'rev.elec.eps.{series}.e'
    tags = sorted({k[len(pre) - 1:].split('.')[0] for k in FACTS if k.startswith(pre) and k.endswith('.lumo_ev')})
    eps = sorted(float(t[1:].replace('p', '.')) for t in tags)
    check(len(eps) == fint(f'rev.elec.eps.{series}.n_points'), f'{series}: n_points keys present')
    return eps


# =================================================================================================
def t_eps():
    wat = fnum('rev.elec.eps.water')
    rows = []
    block = {'gasgeom': f'*Gas-phase minimum*', 'pcmgeom': f'*IEF-PCM minimum*'}
    dscf = {('gasgeom', 1.0): 'gasgeom_vac', ('gasgeom', wat): 'gasgeom_water',
            ('pcmgeom', 1.0): 'pcmgeom_vac', ('pcmgeom', wat): 'pcmgeom_water'}
    for s in ('gasgeom', 'pcmgeom'):
        rows.append([block[s]] + [''] * 7)
        for e in series_eps(s):
            k = f'rev.elec.eps.{s}.{etag(e)}'
            h, l, d = fnum(k + '.homo_ev', 'eV'), fnum(k + '.lumo_ev', 'eV'), fnum(k + '.delta_cat_pic_ev', 'eV')
            sh = fnum(k + '.lumo_picrate_pct', '%')
            on = fstr(k + '.lumo_on')
            check((d < 0) == (on != 'picrate'), f'{s} eps {e:g}: sign of Delta matches the LUMO location')
            row = [eps_str(e), fmt(h, 3), fmt(l, 3), fmt(sh, 1), fmt(d, 3, sign=True)]
            if (s, e) in dscf:
                q = f'rev.elec.dscf.{dscf[(s, e)]}'
                row += [fmt(fnum(q + '.ea_ev', 'eV'), 3), fmt(fnum(q + '.ie_ev', 'eV'), 3),
                        fmt(fnum(q + '.added_electron_picrate_pct', '%'), 1)]
                close(fnum(q + '.koopmans_ea_ev', 'eV'), -l, 1e-4, f'{s} eps {e:g}: Delta-SCF N state = series point')
            else:
                row += [EN, EN, EN]
            rows.append(row)
    sw = {s: fnum(f'rev.elec.eps.{s}.switch_eps') for s in ('gasgeom', 'pcmgeom', 'btpcz_pcmgeom')}
    lo = {s: fnum(f'rev.elec.eps.{s}.switch_eps_lo') for s in ('gasgeom', 'pcmgeom')}
    hi = {s: fnum(f'rev.elec.eps.{s}.switch_eps_hi') for s in ('gasgeom', 'pcmgeom')}
    ref = {s: fnum(f'rev.elec.eps.{s}.refine_eps') for s in ('gasgeom', 'pcmgeom')}
    for s in ('gasgeom', 'pcmgeom'):
        check(lo[s] < sw[s] < hi[s] and ref[s] == lo[s], f'{s}: switch eps inside its bracket, refinement point = lower end')
        check(fstr(f'rev.elec.eps.{s}.switch') == 'cation->picrate', f'{s}: LUMO moves from cation to picrate')
        check(fstr(f'rev.elec.eps.{s}.homo_on_all') == 'picrate', f'{s}: HOMO on the picrate at every eps')
    hmin = min(fnum(f'rev.elec.eps.{s}.homo_picrate_pct_min', '%') for s in ('gasgeom', 'pcmgeom'))
    hmax = max(fnum(f'rev.elec.eps.{s}.homo_picrate_pct_max', '%') for s in ('gasgeom', 'pcmgeom'))
    hole = min(fnum(f'rev.elec.dscf.{c}.hole_picrate_pct', '%') for c in dscf.values())
    split = {s: fnum(f'rev.elec.eps.{s}.near_switch_split_ev', 'eV') for s in ('gasgeom', 'pcmgeom')}
    lsh = {s: fnum(f'rev.elec.eps.{s}.{etag(ref[s])}.lumo_picrate_pct', '%') for s in ('gasgeom', 'pcmgeom')}
    check(abs(fnum('rev.elec.dscf.gasgeom_vac.added_electron_picrate_pct', '%') - 50) < 25,
          'in vacuum at the gas-phase minimum the added electron is shared between the ions')
    caption = (f'Frontier orbitals of the ion pair at fixed geometry as the permittivity *ε* of the continuum is varied, '
               f'{method()}, at the gas-phase and at the IEF-PCM (water) minimum; *ε* = 1 is the calculation without '
               f'continuum and *ε* = {eps_str(wat)} is IEF-PCM (water).')
    notes = [
        f'LUMO on picrate: Mulliken share of the LUMO on the picrate atoms; the HOMO is on the picrate at every point '
        f'({fmt(hmin, 1)}{EN}{fmt(hmax, 1)} %). Δ = *E*(lowest cation-localised virtual) − *E*(lowest '
        f'picrate-localised virtual), so a negative Δ means a cation LUMO. Δ changes sign at *ε* = '
        f'{fmt(sw["gasgeom"], 2)} at the gas-phase minimum and at *ε* = {fmt(sw["pcmgeom"], 2)} at the IEF-PCM minimum '
        f'(linear interpolation in (*ε* − 1)/*ε*); the points at *ε* = {eps_str(ref["gasgeom"])} and '
        f'{eps_str(ref["pcmgeom"])} were added next to these crossings, where the two acceptor levels lie '
        f'{fmt(split["gasgeom"], 3)} and {fmt(split["pcmgeom"], 3)} eV apart and mix (LUMO {fmt(lsh["gasgeom"], 1)} '
        f'and {fmt(lsh["pcmgeom"], 1)} % on the picrate). At the IEF-PCM minimum of the *P*2~1~/*n* form the '
        f'crossing is at *ε* = {fmt(sw["btpcz_pcmgeom"], 2)}.',
        'EA~ΔSCF~ = *E*(N) − *E*(N+1) and IE~ΔSCF~ = *E*(N−1) − *E*(N): vertical values from separate unrestricted '
        'calculations of the anion and cation doublets at the same geometry, with the continuum relaxed for each '
        'charge state. Added electron on picrate: Hirshfeld share of the unpaired spin of the N+1 doublet on the '
        f'picrate atoms; the hole of the N−1 doublet is on the picrate in all four cases (at least {fmt(hole, 1)} %).',
    ]
    header = ['*ε*', '*E*~HOMO~ (eV)', '*E*~LUMO~ (eV)', 'LUMO on picrate (%)', 'Δ (eV)', 'EA~ΔSCF~ (eV)',
              'IE~ΔSCF~ (eV)', 'Added electron on picrate (%)']
    return dict(header=header, align='lrrrrrrr', rows=rows, caption=caption, notes=notes)


# =================================================================================================
def t_tddft():
    systems = [('pair_gas', f'Ion pair, {MODEL["gas"]}'), ('pair_pcm', f'Ion pair, {MODEL["pcm"]}'),
               ('picrate_pcm', f'Free picrate anion, {MODEL["pcm"]}'), ('bta_pcm', f'1*H*-Benzotriazole, {MODEL["pcm"]}')]
    roles = [('s1', None), ('bright', 'first bright'), ('strongest', 'strongest'), ('ctstate', 'lowest CT')]
    rows = []
    for sysname, label in systems:
        for xc in ('b3lyp', 'camb3lyp'):
            base = f'rev.elec.td.{sysname}.{xc}'
            pair = sysname.startswith('pair')
            listed = {}
            for r, name in roles:
                if f'{base}.{r}_state' not in FACTS and r == 'ctstate' and pair:
                    mt.fail(f'{base}: no picrate → cation CT state among the computed singlets')
                if f'{base}.{r}_state' not in FACTS:
                    continue
                n = fint(f'{base}.{r}_state')
                listed.setdefault(n, (r, []))
                if name:
                    listed[n][1].append(name)
            first = True
            for n in sorted(listed):
                r, names = listed[n]
                e, lam, f = fnum(f'{base}.{r}_ev', 'eV'), fnum(f'{base}.{r}_nm', 'nm'), fnum(f'{base}.{r}_f')
                st = f'S~{n}~' + (f' ({", ".join(names)})' if names else '')
                ct = fmt(fnum(f'{base}.{r}_ct_number'), 2) if pair else EN
                rows.append([label if first else '', XC[xc] if first else '', st, fmt(e, 3), f'{lam:.0f}', fmt(f, 3),
                             ct, fstr(f'{base}.{r}_character')])
                first = False
    eo = fnum('rev.elec.td.eps_optical')
    ns = fint('rev.elec.td.nstates')
    fb, ctt = fnum('rev.elec.td.bright_threshold_f'), fnum('rev.elec.td.ct_threshold')
    check(fstr('rev.elec.td.pair_gas.b3lyp.s1_character') == 'picrate → cation', 'B3LYP gas S1 is the CT state')
    caption = (f'Singlet excited states from TD-DFT in the Tamm–Dancoff approximation ({ns} states per calculation, '
               f'6-311++G(d,p)): the ion pair at its gas-phase minimum in the gas phase and at its IEF-PCM minimum in '
               f'IEF-PCM (water), and the free picrate anion and neutral 1*H*-benzotriazole optimised in IEF-PCM (water).')
    notes = [
        f'Listed per calculation: S~1~, the first state with *f* > {fb:g} (first bright), the strongest of the {ns} '
        'states and the lowest picrate → cation charge-transfer (CT) state, where these differ. CT number: share of '
        f'the transition density (Löwdin partition) that connects the two ions; a state is classed as CT above {ctt:g}, '
        'and its direction follows from where the hole and the electron sit. In IEF-PCM (water) the continuum '
        f'responds to the excitation with the optical permittivity {fmt(eo, 2)} (non-equilibrium) while the ground state '
        'is equilibrated with the static one. CAM-B3LYP: [@yanai2004hybrid].',
    ]
    header = ['System and model', 'Functional', 'State', '*E* (eV)', '*λ* (nm)', '*f*', 'CT number', 'Character']
    return dict(header=header, align='lllrrrrl', rows=rows, caption=caption, notes=notes)


# =================================================================================================
def t_charges():
    cols = [('gas', MODEL['gas']), ('pcm', MODEL['pcm']), ('gasgeom_water', f'Gas-phase minimum in {MODEL["pcm"]}'),
            ('pcmgeom_vac', 'IEF-PCM minimum in the gas phase'), ('btpcz_gas', '*P*2~1~/*n* form, gas phase'),
            ('btpcz_pcm', f'*P*2~1~/*n* form, {MODEL["pcm"]}')]
    schemes = [('mulliken', 'Mulliken'), ('lowdin', 'Löwdin'), ('hirshfeld', 'Hirshfeld'), ('cm5', 'CM5')]
    rows = []

    def cell(key):
        return fmt(fnum(key, 'e'), 3) if key in FACTS else EN
    for m, name in schemes:
        rows.append([f'Picrate, {name}'] + [cell(f'rev.elec.chg.{c}.picrate_{m}_e') for c, _ in cols])
    for m, name in (('hirshfeld', 'Hirshfeld'), ('cm5', 'CM5')):
        rows.append([f'Picrate, {name} (Multiwfn)'] + [cell(f'rev.elec.chg.{c}.picrate_{m}_multiwfn_e') for c, _ in cols])
    for m, name in schemes:
        rows.append([f'Inter-ion H, {name}'] + [cell(f'rev.elec.chg.{c}.proton_{m}_e') for c, _ in cols])
    for c in ('gas', 'pcm'):
        check(fnum(f'rev.elec.chg.{c}.check_vs_published_max_e', 'e') < 1e-6,
              f'{c}: Mulliken/Löwdin recomputed from the checkpoint = published values')
    ph = {m: fnum(f'rev.elec.chg.btpcz_gas.picric_acid_{m}_e', 'e') for m, _ in schemes}
    fe = fnum('rev.elec.chg.hirshfeld_field_effect_gasgeom_e', 'e')
    ge = fnum('rev.elec.chg.hirshfeld_geom_effect_vac_e', 'e')
    close(fe, fnum('rev.elec.chg.gasgeom_water.picrate_hirshfeld_e', 'e') - fnum('rev.elec.chg.gas.picrate_hirshfeld_e', 'e'),
          1e-12, 'field effect = cross point - gas')
    caption = (f'Net charge of the picrate atoms and charge of the inter-ion hydrogen (e), {method()}, at the two '
               'minima of the *P*2~1~/*c* pair (Gas phase, IEF-PCM (water)), at the two cross single points (each '
               'minimum in the other model) and at the two minima of the *P*2~1~/*n* form.')
    notes = [
        'Hirshfeld charges [@hirshfeld1977bonded] were integrated on the molecular DFT grid with spherically averaged '
        'B3LYP free-atom densities; CM5 [@marenich2012charge] adds the CM5 pair correction to them. The rows marked '
        'Multiwfn repeat both with that program’s built-in free-atom densities as a check. The cross single points have '
        'no CM5 values. The inter-ion hydrogen is H4N in the *P*2~1~/*c* pair and H6 in the *P*2~1~/*n* form. At the '
        'gas-phase minimum of the *P*2~1~/*n* form that hydrogen sits on O1, so the picrate atoms there belong to '
        f'picric acid, whose net charge including H6 is {fmt(ph["hirshfeld"], 3)} e (Hirshfeld) and '
        f'{fmt(ph["cm5"], 3)} e (CM5). Moving from the gas-phase to the IEF-PCM geometry in vacuum changes the '
        f'Hirshfeld picrate charge by {fmt(ge, 3)} e; placing the gas-phase geometry in the continuum changes it by '
        f'{fmt(fe, 3)} e.',
    ]
    header = ['Charge (e)'] + [n for _, n in cols]
    return dict(header=header, align='l' + 'r' * len(cols), rows=rows, caption=caption, notes=notes)


# =================================================================================================
def t_functionals():
    wat = fnum('rev.elec.eps.water')
    rows = []
    for xc in ('b3lyp', 'pbe0', 'camb3lyp'):
        for g, gname in (('gas', 'Gas-phase minimum'), ('pcm', 'IEF-PCM minimum')):
            for k, e in enumerate((1.0, wat)):
                key = (f'rev.elec.eps.{g}geom.{etag(e)}' if xc == 'b3lyp' else f'rev.elec.fx.{xc}.{g}geom.{etag(e)}')
                h, l, gap = fnum(key + '.homo_ev', 'eV'), fnum(key + '.lumo_ev', 'eV'), fnum(key + '.gap_ev', 'eV')
                close(gap, l - h, 1e-9, f'{xc} {g} eps {e:g}: gap = LUMO - HOMO')
                sw = fmt(fnum(f'rev.elec.fx.{xc}.{g}geom.switch_eps_2pt'), 2) if k == 0 else ''
                rows.append([XC[xc] if (g == 'gas' and k == 0) else '', gname if k == 0 else '', eps_str(e),
                             fmt(h, 3), fmt(l, 3), fmt(gap, 3), fmt(fnum(key + '.lumo_picrate_pct', '%'), 1),
                             fmt(fnum(key + '.delta_cat_pic_ev', 'eV'), 3, sign=True), sw])
    check(fstr('rev.elec.fx.lumo_location_agrees') == 'yes', 'the three functionals agree on the LUMO location')
    spread = fnum('rev.elec.fx.delta_spread_max_ev', 'eV')
    full = {g: fnum(f'rev.elec.eps.{g}geom.switch_eps') for g in ('gas', 'pcm')}
    caption = ('Frontier orbitals of the ion pair with three functionals at the B3LYP minima (single points, '
               f'6-311++G(d,p)), without continuum (*ε* = 1) and in IEF-PCM (water) (*ε* = {eps_str(wat)}).')
    notes = [
        'PBE0 [@adamo1999reliable] and CAM-B3LYP [@yanai2004hybrid] shift all orbital energies but leave Δ = '
        '*E*(lowest cation-localised virtual) − *E*(lowest picrate-localised virtual) within '
        f'{fmt(spread, 2)} eV of one another, and all three put the LUMO on the cation without continuum and on '
        'the picrate in water at both minima. Switch *ε*: where Δ changes sign by linear interpolation in '
        '(*ε* − 1)/*ε* between the two points of each row pair; the full B3LYP series places the same crossings at '
        f'*ε* = {fmt(full["gas"], 2)} and {fmt(full["pcm"], 2)}, so this two-point estimate runs low and serves only to '
        'compare the functionals.',
    ]
    header = ['Functional', 'Geometry', '*ε*', '*E*~HOMO~ (eV)', '*E*~LUMO~ (eV)', 'Δ*E*~gap~ (eV)',
              'LUMO on picrate (%)', 'Δ (eV)', 'Switch *ε* (two-point)']
    return dict(header=header, align='lllrrrrrr', rows=rows, caption=caption, notes=notes)


TABLES = [('tbl_rev_eps', 'main', 'Frontier orbitals versus the permittivity of the continuum, with ΔSCF', t_eps),
          ('tbl_rev_tddft', 'main', 'TD-DFT (TDA) excited states of the pair and of the free ions', t_tddft),
          ('tblS_rev_charges', 'SI', 'Mulliken, Löwdin, Hirshfeld and CM5 fragment charges', t_charges),
          ('tblS_rev_functionals', 'SI', 'B3LYP, PBE0 and CAM-B3LYP frontier orbitals', t_functionals)]


def main():
    for stem, where, title, fn in TABLES:
        CTX.stem, CTX.keys, CTX.files, CTX.checks = stem, [], [str(FACTS_PATH.relative_to(REPO))], []
        spec = fn()
        t = Table(stem=stem, where=where, title=title, keys=list(CTX.keys), files=list(CTX.files),
                  checks=list(CTX.checks), **spec)
        text = render(t)
        lint(text, stem)
        (HERE / f'{stem}.md').write_text(text, encoding='utf-8')
        print(f'{where:>4}  {stem:<22} {len(t.rows):3d} rows  {len(t.keys):3d} keys  {len(t.checks):3d} checks')


if __name__ == '__main__':
    try:
        main()
    except TableError as e:
        print(f'ERROR {e}', file=sys.stderr)
        sys.exit(1)
