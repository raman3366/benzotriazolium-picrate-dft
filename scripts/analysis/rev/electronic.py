"""Revision facts: electronic structure (COMPTC-D-26-02710, revision round 1; plan section 5, items 3, 7 and 8).

Where does the LUMO of the ion pair live, and is it the geometry or the dielectric field that decides?  Read-only
inputs are the result.json files of these jobs (data/revision_runs/queue*/out/<job>/), all B3LYP/6-311++G(d,p) with
density fitting, grid level 4, IEF-PCM Lebedev order 29 unless stated:

  headline_gasgeom, headline_pcmgeom  permittivity series eps = 1, 2, 4, 8, 20, 78.3553 (eps = 1: no PCM) at the
                                      published P2_1/c gas-phase and IEF-PCM minima, plus one refinement point at
                                      the interpolated switch (1.64 and 3.153: NOT crystal permittivities)
  headline_btpcz                      the same at the P2_1/n minima; its gas-phase minimum is the neutral pair
                                      (proton on O1), its IEF-PCM minimum the salt (refinement point 2.963)
  deltascf_gasgeom, deltascf_pcmgeom  Delta-SCF N+1 / N-1 (UKS doublets) at both geometries, in vacuum and in water
  charges_stored, charges_btpcz       Hirshfeld and CM5 (with Mulliken, Loewdin) at the four minima
  td_pair_b3lyp, td_pair_camb3lyp     TDA, 20 singlets, pair at the gas minimum (gas) and IEF-PCM minimum (PCM)
  td_free                             free picrate and neutral 1H-benzotriazole optimised in IEF-PCM, TDA B3LYP/CAM
  sp8_pbe0, sp8_camb3lyp              PBE0 and CAM-B3LYP at both P2_1/c minima, eps = 1 and 78.3553

Key layout (rev.elec.*):
  eps.<series>.<etag>.<q>   series gasgeom | pcmgeom | btpcz_gasgeom | btpcz_pcmgeom; etag e1, e1p64, e2, ...,
                            e78p3553 (the eps printed with '.' -> 'p'); q homo_ev, lumo_ev, gap_ev, lumo_on,
                            lumo_picrate_pct, delta_cat_pic_ev, mulliken_picrate_e, dipole_debye
  eps.<series>.<q>          switch_eps (+ _ln, _lo, _hi, _linfit), fit_*, shifts of the localised levels, ...
  eps.field_* / eps.geom_*  geometry-versus-field decomposition of Delta
  dscf.<combo>.<q>          combo gasgeom_vac | gasgeom_water | pcmgeom_water | pcmgeom_vac
  chg.<tag>.<q>             tag gas | pcm | btpcz_gas | btpcz_pcm | gasgeom_water | pcmgeom_vac
  td.<system>.<xc>.<q>      system pair_gas | pair_pcm | picrate_pcm | bta_pcm; xc b3lyp | camb3lyp; q
                            <sel>_state/_ev/_nm/_f/_character/_ct_number with sel s1 | bright |
                            strongest | ctstate (lowest picrate -> cation CT state)
  fx.<xc>.<geom>.<etag>.<q> PBE0 / CAM-B3LYP points; fx.<xc>.<geom>.switch_eps_2pt also for b3lyp
"Picrate" and "cation" are the atom sets of the case files (picrate atoms 0-17, cation atoms 18-32).  At the neutral
P2_1/n gas minimum they are picric acid without its O-H hydrogen and benzotriazole carrying that hydrogen.

Derived CSVs (written only in the full build, i.e. when facts already holds the published keys):
  papers/btpc/data/derived/rev_elec_eps_series.csv, rev_elec_deltascf.csv, rev_elec_charges.csv,
  rev_elec_tddft.csv, rev_elec_functionals.csv
"""
from .common import *
import csv, math, os, re

import numpy as np

FACTS_JSON = ROOT / 'papers/btpc/data/facts.json'
DERIVED = ROOT / 'papers/btpc/data/derived'
CT_MIN = 0.5                    # a TD state is 'charge transfer' when the inter-fragment share of Omega exceeds this
LOC_MIN = 50.0                  # headline.py: an orbital is localised on a fragment when its Mulliken share > 50 %
FRAG = {'anion': 'picrate', 'cation': 'cation', 'mixed': 'mixed'}
XC_NAME = {'b3lyp': 'B3LYP', 'pbe0': 'PBE0', 'camb3lyp': 'CAM-B3LYP'}

SERIES = [  # name, job, geometry label inside the job, case, plain description
    ('gasgeom', 'headline_gasgeom', 'gas', 'rev_gas', 'P2_1/c pair at its gas-phase minimum'),
    ('pcmgeom', 'headline_pcmgeom', 'pcm', 'rev_pcm', 'P2_1/c pair at its IEF-PCM minimum'),
    ('btpcz_gasgeom', 'headline_btpcz', 'gas', 'btpcz_gas',
     'P2_1/n pair at its gas-phase minimum (neutral: proton on O1)'),
    ('btpcz_pcmgeom', 'headline_btpcz', 'pcm', 'btpcz_pcm', 'P2_1/n pair at its IEF-PCM minimum (salt)'),
]
DSCF = [  # name, job, combo (geometry/model), case, description
    ('gasgeom_vac', 'deltascf_gasgeom', 'gas/gas', 'rev_gas', 'gas-phase minimum in vacuum'),
    ('gasgeom_water', 'deltascf_gasgeom', 'gas/pcm', 'rev_gas', 'gas-phase minimum in IEF-PCM water'),
    ('pcmgeom_water', 'deltascf_pcmgeom', 'pcm/pcm', 'rev_pcm', 'IEF-PCM minimum in IEF-PCM water'),
    ('pcmgeom_vac', 'deltascf_pcmgeom', 'pcm/gas', 'rev_pcm', 'IEF-PCM minimum in vacuum'),
]
CHG = [  # tag, job, item, case, description
    ('gas', 'charges_stored', 'gas', 'rev_gas', 'P2_1/c gas-phase minimum (published final.chk)'),
    ('pcm', 'charges_stored', 'pcm', 'rev_pcm', 'P2_1/c IEF-PCM minimum (published final.chk)'),
    ('btpcz_gas', 'charges_btpcz', 'gas', 'btpcz_gas', 'P2_1/n gas-phase minimum (neutral pair)'),
    ('btpcz_pcm', 'charges_btpcz', 'pcm', 'btpcz_pcm', 'P2_1/n IEF-PCM minimum'),
]
SCHEMES = ('mulliken', 'lowdin', 'hirshfeld', 'cm5')
SCHEME_NAME = {'mulliken': 'Mulliken', 'lowdin': 'Loewdin', 'hirshfeld': 'Hirshfeld', 'cm5': 'CM5'}
TD = [  # system, job, item key, xc, case
    ('pair_gas', 'td_pair_b3lyp', 'pair-gas/b3lyp', 'b3lyp', 'rev_gas'),
    ('pair_gas', 'td_pair_camb3lyp', 'pair-gas/camb3lyp', 'camb3lyp', 'rev_gas'),
    ('pair_pcm', 'td_pair_b3lyp', 'pair-pcm/b3lyp', 'b3lyp', 'rev_pcm'),
    ('pair_pcm', 'td_pair_camb3lyp', 'pair-pcm/camb3lyp', 'camb3lyp', 'rev_pcm'),
    ('picrate_pcm', 'td_free', 'picrate-pcm/b3lyp', 'b3lyp', 'free_picrate_pcm'),
    ('picrate_pcm', 'td_free', 'picrate-pcm/camb3lyp', 'camb3lyp', 'free_picrate_pcm'),
    ('bta_pcm', 'td_free', 'bta-pcm/b3lyp', 'b3lyp', 'free_bta_pcm'),
    ('bta_pcm', 'td_free', 'bta-pcm/camb3lyp', 'camb3lyp', 'free_bta_pcm'),
]
TD_DESC = {'pair_gas': 'ion pair at the gas-phase minimum, Gas phase',
           'pair_pcm': 'ion pair at the IEF-PCM minimum, IEF-PCM (water), non-equilibrium LR-PCM',
           'picrate_pcm': 'free picrate anion optimised (B3LYP) in IEF-PCM (water)',
           'bta_pcm': 'free neutral 1H-benzotriazole optimised (B3LYP) in IEF-PCM (water)'}
FX = [('pbe0', 'sp8_pbe0'), ('camb3lyp', 'sp8_camb3lyp')]


# ============================================================================ helpers
def etag(eps):
    return 'e' + f'{eps:g}'.replace('.', 'p')


def fvar(eps):
    """Pekar-type factor f = (eps - 1)/eps (0 in vacuum, -> 1 for a conductor)."""
    return (eps - 1.0) / eps


def switch_f(points):
    """First sign change of Delta(eps) in a sorted [(eps, delta)] list, interpolated linearly in f (and in ln eps),
    exactly as dft/revision/headline.py switch_from does; None when Delta keeps its sign."""
    pts = sorted(p for p in points if p[1] is not None)
    for (e0, d0), (e1, d1) in zip(pts, pts[1:]):
        if d0 * d1 < 0:
            x = fvar(e0) + (0 - d0) * (fvar(e1) - fvar(e0)) / (d1 - d0)
            l = math.log(e0) + (0 - d0) * (math.log(e1) - math.log(e0)) / (d1 - d0)
            return {'f': 1 / (1 - x), 'ln': math.exp(l), 'lo': e0, 'hi': e1, 'd_lo': d0, 'd_hi': d1}
    return None


_PUB = None


def pub(facts, key):
    """A published fact: from the facts being built if present, else from papers/btpc/data/facts.json on disk."""
    global _PUB
    if key in facts:
        return facts[key]['value']
    if _PUB is None:
        _PUB = json.load(open(FACTS_JSON))
    return _PUB[key]['value']


def check_settings(d, job, xc='b3lyp'):
    s = d['settings']
    assert d['status'] == 'done' and not d.get('smoke'), (job, d['status'])
    assert s['basis'] == '6-311++g(d,p)' and s['grid_level'] == 4 and s['lebedev_order'] == 29, (job, s)
    xs = s['xc'] if isinstance(s['xc'], list) else [s['xc']]
    assert xc in xs, (job, s['xc'])
    return s


def write_csv(name, header, rows):
    p = DERIVED / name
    tmp = p.with_suffix('.csv.tmp')
    with open(tmp, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(['' if v is None else (f'{v:.15g}' if isinstance(v, float) else v) for v in r])
    os.replace(tmp, p)


def loc(sh):
    """'picrate' | 'cation' | 'mixed' from a {'anion': %, 'cation': %} share dict (> 50 % rule of headline.py)."""
    best = max(sh, key=sh.get)
    return FRAG[best] if sh[best] > LOC_MIN else 'mixed'


# ============================================================================ build
def build(add, facts):
    full_build = 'gas.homo_ev' in facts            # standalone runs (python -m rev.common) write no files
    A = lambda key, value, unit='', case=None, source='', note='': add(
        'rev.elec.' + key, (float(value) if isinstance(value, (np.floating, float)) else
                            int(value) if isinstance(value, np.integer) else value), unit, case, source, note)
    csv_rows = {}

    # ------------------------------------------------------------------ 1. permittivity series (headline jobs)
    series = {}
    eps_rows = []
    lv_src = None
    for name, job, g, case, desc in SERIES:
        d, src = result(job)
        s = check_settings(d, job)
        assert s['xc'] == 'b3lyp' and s['disp'] is None and s['conv_tol'] == 1e-10, (job, s)
        lv_src = lv_src or (s, src)
        items = sorted((it for it in d['items'].values() if it['geometry'] == g and it['status'] == 'done'),
                       key=lambda it: it['eps'])
        key_of = {it['eps']: k for k, it in d['items'].items() if it['geometry'] == g}
        pts = []
        for it in items:
            eps = it['eps']; fr = it['frontier']; o = fr['orbitals']; k = key_of[eps]
            assert it['settings']['xc'] == 'b3lyp'
            assert (it['settings']['solvent'] is None) == (eps == 1.0), (job, k)
            if eps != 1.0:
                assert abs(it['settings']['solvent']['eps'] - eps) < 1e-12 and it['settings']['solvent']['model'] == 'IEF-PCM'
            lv = fr['lowest_virtual_on']; ho = fr['highest_occupied_on']
            p = {'eps': eps, 'f': fvar(eps), 'key': k, 'refinement': it.get('refinement'), 'E_Ha': it['E_Ha'],
                 'homo': fr['HOMO_eV'], 'lumo': fr['LUMO_eV'], 'gap': fr['gap_eV'],
                 'lumo_pic': o['LUMO']['mulliken']['anion'], 'lumo_pic_lowdin': o['LUMO']['lowdin']['anion'],
                 'homo_pic': o['HOMO']['mulliken']['anion'], 'lumo_on': FRAG[fr['LUMO_on']],
                 'homo_on': FRAG[fr['HOMO_on']], 'lumo1': o['LUMO+1']['eV'],
                 'lumo1_pic': o['LUMO+1']['mulliken']['anion'],
                 'cat_virt': lv['cation']['eV'], 'cat_virt_label': lv['cation']['label'],
                 'pic_virt': lv['anion']['eV'], 'pic_virt_label': lv['anion']['label'],
                 'cat_occ': ho.get('cation', {}).get('eV'), 'cat_occ_label': ho.get('cation', {}).get('label'),
                 'delta': fr['delta_virt_cation_minus_anion_eV'],
                 'mull': it['charges']['mulliken']['anion'], 'lowd': it['charges']['lowdin']['anion'],
                 'dip': it['dipole_debye']['total'], 'dE_stored': it.get('dE_vs_stored_Ha')}
            # the definitions behind the stored labels, re-checked here
            assert abs(p['delta'] - (p['cat_virt'] - p['pic_virt'])) < 1e-9
            assert abs(p['gap'] - (p['lumo'] - p['homo'])) < 1e-9
            assert p['lumo_on'] == loc(o['LUMO']['mulliken']) and p['homo_on'] == loc(o['HOMO']['mulliken'])
            pts.append(p)
        series[name] = pts
        base = f'eps.{name}'
        for p in pts:
            t = f"{base}.{etag(p['eps'])}"
            srcp = f"{src}:items['{p['key']}']"
            what = (f"{desc}, B3LYP, " + ('no continuum (eps = 1)' if p['eps'] == 1 else
                                          f"IEF-PCM eps = {p['eps']:g}") +
                    (' (refinement point placed by headline.py at the switch eps interpolated from the series)'
                     if p['refinement'] else ''))
            A(f'{t}.homo_ev', p['homo'], 'eV', case, srcp + '.frontier.HOMO_eV', f'HOMO energy; {what}')
            A(f'{t}.lumo_ev', p['lumo'], 'eV', case, srcp + '.frontier.LUMO_eV', f'LUMO energy; {what}')
            A(f'{t}.gap_ev', p['gap'], 'eV', case, srcp + '.frontier.gap_eV', f'HOMO-LUMO gap; {what}')
            A(f'{t}.lumo_on', p['lumo_on'], '', case, srcp + '.frontier.LUMO_on',
              f"fragment holding > 50 % (Mulliken) of the LUMO; {what}; LUMO+1 is {p['lumo1']:.4f} eV with "
              f"{p['lumo1_pic']:.1f} % on the picrate")
            A(f'{t}.lumo_picrate_pct', p['lumo_pic'], '%', case, srcp + ".frontier.orbitals.LUMO.mulliken.anion",
              f"Mulliken share of the LUMO on the picrate atoms (Loewdin {p['lumo_pic_lowdin']:.2f} %); {what}")
            A(f'{t}.delta_cat_pic_ev', p['delta'], 'eV', case,
              srcp + '.frontier.delta_virt_cation_minus_anion_eV',
              f"Delta = E(lowest cation-localised virtual, {p['cat_virt_label']}) - E(lowest picrate-localised "
              f"virtual, {p['pic_virt_label']}); negative: the LUMO is on the cation; {what}")
            A(f'{t}.mulliken_picrate_e', p['mull'], 'e', case, srcp + '.charges.mulliken.anion',
              f"Mulliken charge summed over the picrate atoms (Loewdin {p['lowd']:.4f} e); {what}")
            A(f'{t}.dipole_debye', p['dip'], 'D', case, srcp + '.dipole_debye.total', f'dipole moment; {what}')
            eps_rows.append([name, case, g, p['eps'], p['f'], p['refinement'] or '', p['E_Ha'], p['homo'],
                             p['homo_on'], p['homo_pic'], p['lumo'], p['lumo_on'], p['lumo_pic'], p['lumo_pic_lowdin'],
                             p['lumo1'], p['lumo1_pic'], p['gap'], p['cat_virt_label'], p['cat_virt'],
                             p['pic_virt_label'], p['pic_virt'], p['delta'], p['cat_occ_label'], p['cat_occ'],
                             p['mull'], p['lowd'], p['dip'], f"{src}:items['{p['key']}']"])

        A(f'{base}.n_points', len(pts), '', case, src + ':items', f'points in the series ({desc}): eps = ' +
          ', '.join(f"{p['eps']:g}" for p in pts))
        refs = [p['eps'] for p in pts if p['refinement']]
        homo_on = sorted({p['homo_on'] for p in pts}); lumo_on = [p['lumo_on'] for p in pts]
        A(f'{base}.homo_on_all', homo_on[0] if len(homo_on) == 1 else 'varies', '', case, src + ':items[*].frontier.HOMO_on',
          f'fragment holding the HOMO at every eps of the series ({desc}); Mulliken HOMO share on the picrate '
          f"{min(p['homo_pic'] for p in pts):.1f}-{max(p['homo_pic'] for p in pts):.1f} %")
        for q, fn in (('min', min), ('max', max)):
            A(f'{base}.homo_picrate_pct_{q}', fn(p['homo_pic'] for p in pts), '%', case,
              src + ':items[*].frontier.orbitals.HOMO.mulliken.anion',
              f'{q}imum over the series of the Mulliken share of the HOMO on the picrate atoms; {desc}')
        # switch eps, recomputed here from the points with the formula of headline.py and checked against the job
        sw = switch_f([(p['eps'], p['delta']) for p in pts])
        jsw = d['switch'][g]
        if sw is None:
            assert jsw.get('eps_switch_f') is None
            A(f'{base}.switch', 'none', '', case, src + f":switch['{g}']",
              f"Delta keeps its sign over eps {pts[0]['eps']:g}-{pts[-1]['eps']:g}: LUMO on the "
              f"{lumo_on[0]} at every eps ({desc})")
            assert len(set(lumo_on)) == 1
        else:
            assert abs(sw['f'] - jsw['eps_switch_f']) < 1e-9 and abs(sw['ln'] - jsw['eps_switch_ln']) < 1e-9
            assert [sw['lo'], sw['hi']] == jsw['bracket']
            A(f'{base}.switch', jsw['direction'].replace('anion', 'picrate'), '', case, src + f":switch['{g}'].direction",
              f'direction in which the LUMO moves as eps grows ({desc})')
            A(f'{base}.switch_eps', sw['f'], '', case, src + f":switch['{g}'].eps_switch_f",
              f"eps at which Delta changes sign, linear interpolation in f = (eps-1)/eps between eps {sw['lo']:g} "
              f"(Delta {sw['d_lo']:+.4f} eV) and {sw['hi']:g} (Delta {sw['d_hi']:+.4f} eV); recomputed here from the "
              f'points; {desc}')
            A(f'{base}.switch_eps_ln', sw['ln'], '', case, src + f":switch['{g}'].eps_switch_ln",
              f'the same crossing interpolated linearly in ln eps (sensitivity of the interpolation); {desc}')
            A(f'{base}.switch_eps_lo', sw['lo'], '', case, src + f":switch['{g}'].bracket[0]",
              f'largest computed eps with the LUMO still on the cation side (Delta < 0); {desc}')
            A(f'{base}.switch_eps_hi', sw['hi'], '', case, src + f":switch['{g}'].bracket[1]",
              f'smallest computed eps with the LUMO on the picrate (Delta > 0); {desc}')
        if refs:
            r = next(p for p in pts if p['refinement'])
            A(f'{base}.refine_eps', r['eps'], '', case, src + f":items['{r['key']}'].refinement",
              'eps of the extra single point that headline.py (--refine 1) placed at the switch eps interpolated '
              'from the six-point series (eps 1, 2, 4, 8, 20, 78.3553); a numerical refinement point, not a '
              f'measured or crystal permittivity; {desc}')
            A(f'{base}.near_switch_split_ev', r['lumo1'] - r['lumo'], 'eV', case,
              src + f":items['{r['key']}'].frontier.orbitals",
              f"E(LUMO+1) - E(LUMO) at eps {r['eps']:g}, next to the crossing: the two acceptor levels mix there "
              f"(LUMO {r['lumo_pic']:.1f} % and LUMO+1 {r['lumo1_pic']:.1f} % on the picrate), i.e. the switch is "
              f'an avoided crossing of a cation-localised and a picrate-localised level; {desc}')
        # Delta is close to linear in f: Delta = a + b f (least squares over all points)
        fx = np.array([p['f'] for p in pts]); dy = np.array([p['delta'] for p in pts])
        b, a = np.polyfit(fx, dy, 1)
        r2 = 1 - ((dy - (a + b * fx)) ** 2).sum() / ((dy - dy.mean()) ** 2).sum()
        A(f'{base}.fit_slope_ev', b, 'eV', case, src + ':items[*].frontier.delta_virt_cation_minus_anion_eV',
          f'slope b of the least-squares line Delta = a + b f over the {len(pts)} points, f = (eps-1)/eps (the '
          f'Born/Pekar factor): how far the field moves the cation level against the picrate level; {desc}')
        A(f'{base}.fit_intercept_ev', a, 'eV', case, src + ':items[*].frontier.delta_virt_cation_minus_anion_eV',
          f'intercept a (vacuum value of the fitted line) of Delta = a + b f; {desc}')
        A(f'{base}.fit_r2', r2, '', case, src + ':items[*].frontier.delta_virt_cation_minus_anion_eV',
          f'coefficient of determination of Delta = a + b f over {len(pts)} points; {desc}')
        if a * (a + b * fvar(pts[-1]['eps'])) < 0:
            A(f'{base}.switch_eps_linfit', 1 / (1 + a / b), '', case,
              src + ':items[*].frontier.delta_virt_cation_minus_anion_eV',
              f'eps where the fitted line crosses zero: f* = -a/b, eps = 1/(1 - f*); {desc}')
        # how far each localised level moves between vacuum and water (fixed geometry)
        p1, pw = pts[0], pts[-1]
        assert p1['eps'] == 1.0
        for q, k, lab in (('cation_virtual_shift_ev', 'cat_virt', 'lowest cation-localised virtual'),
                          ('picrate_virtual_shift_ev', 'pic_virt', 'lowest picrate-localised virtual'),
                          ('homo_shift_ev', 'homo', 'HOMO')):
            A(f'{base}.{q}', pw[k] - p1[k], 'eV', case, src + f":items['{pw['key']}'] - items['{p1['key']}']",
              f"E({lab}) at eps {pw['eps']:g} minus at eps 1 ({p1[k]:.4f} -> {pw[k]:.4f} eV); {desc}")
        if name in ('gasgeom', 'pcmgeom'):
            ref = pts[0] if name == 'gasgeom' else pts[-1]
            m = 'gas' if name == 'gasgeom' else 'pcm'
            A(f'{base}.dE_vs_published_ha', ref['dE_stored'], 'Ha', case, src + f":items['{ref['key']}'].dE_vs_stored_Ha",
              f"E({'eps 1' if m == 'gas' else 'eps 78.3553'}) minus the published {m}.e_total: the series reproduces "
              'the published calculation (tolerance 1e-6 Ha)')
            assert abs(ref['dE_stored']) < 1e-6
            assert abs(ref['homo'] - pub(facts, f'{m}.homo_ev')) < 1e-3 and abs(ref['lumo'] - pub(facts, f'{m}.lumo_ev')) < 1e-3
            assert abs(ref['lumo_pic'] - pub(facts, f'{m}.lumo_picrate_pct')) < 0.051
            assert abs(ref['mull'] - pub(facts, f'{m}.mulliken_charge_picrate')) < 1e-6
            assert abs(ref['dip'] - pub(facts, f'{m}.dipole_debye')) < 1e-3
        if name == 'btpcz_gasgeom':
            pc = json.load(open(ROOT / 'data/revision_runs/repo_results/btpcz/btpcz_gas/proton_check.json'))
            assert pc['proton_transferred'] and d['inputs']['geometry_gas']['path'].endswith('btpcz_gas/opt.xyz')
            A(f'{base}.proton_on', 'O1', '', case, 'data/revision_runs/repo_results/btpcz/btpcz_gas/proton_check.json',
              f"the P2_1/n gas-phase minimum is the neutral pair: H6 on the picrate O1 (O1-H {pc['d_H_O_A']:.3f} A, "
              f"H...N6 {pc['d_H_N_A']:.3f} A, bonds broken {pc['bonds_broken']}, formed {pc['bonds_formed']}); its "
              "'picrate' atom set is picric acid without H6 and its 'cation' atom set is benzotriazole plus H6")
    wat = series['gasgeom'][-1]['eps']
    for nm in series:
        assert series[nm][0]['eps'] == 1.0 and series[nm][-1]['eps'] == wat
    gs, ps = {p['eps']: p for p in series['gasgeom']}, {p['eps']: p for p in series['pcmgeom']}
    dsrc = f"{result('headline_gasgeom')[1]} and {result('headline_pcmgeom')[1]}"
    A('eps.field_shift_gasgeom_ev', gs[wat]['delta'] - gs[1.0]['delta'], 'eV', 'rev_gas', dsrc,
      f'Delta(gas-phase geometry, eps {wat:g}) - Delta(gas-phase geometry, eps 1): the effect of the continuum alone')
    A('eps.field_shift_pcmgeom_ev', ps[wat]['delta'] - ps[1.0]['delta'], 'eV', 'rev_pcm', dsrc,
      f'Delta(IEF-PCM geometry, eps {wat:g}) - Delta(IEF-PCM geometry, eps 1): the effect of the continuum alone')
    A('eps.geom_shift_vac_ev', ps[1.0]['delta'] - gs[1.0]['delta'], 'eV', None, dsrc,
      'Delta(IEF-PCM geometry, eps 1) - Delta(gas-phase geometry, eps 1): the effect of the geometry alone, in vacuum; '
      'negative: the IEF-PCM geometry holds the LUMO more firmly on the cation')
    A('eps.geom_shift_water_ev', ps[wat]['delta'] - gs[wat]['delta'], 'eV', None, dsrc,
      f'Delta(IEF-PCM geometry, eps {wat:g}) - Delta(gas-phase geometry, eps {wat:g}): the effect of the geometry alone, '
      'in the continuum')
    A('eps.btpcz_minus_btpc_pcmgeom_switch_eps',
      facts['rev.elec.eps.btpcz_pcmgeom.switch_eps']['value'] - facts['rev.elec.eps.pcmgeom.switch_eps']['value'],
      '', None, dsrc + ' and ' + result('headline_btpcz')[1],
      'switch eps at the P2_1/n IEF-PCM minimum minus that at the P2_1/c IEF-PCM minimum')
    s, src = lv_src
    A('eps.protocol', f"{s['xc'].upper()}/{s['basis'].replace('g(d,p)', 'G(d,p)')}, density fitting "
      f"({s['auxbasis']} auxiliary basis), grid level {s['grid_level']}, {s['pcm_method']} Lebedev order "
      f"{s['lebedev_order']}, SCF conv_tol {s['conv_tol']:g} Ha; eps = 1 is a calculation without continuum; "
      'fixed geometries; each SCF started from the converged density of the previous eps', '', None,
      src + ':settings', 'settings of the permittivity series (headline.py), identical in the four series')
    A('eps.water', wat, '', None, src + ':settings.eps_list[-1]', 'eps of water in PySCF (the published IEF-PCM model)')
    A('eps.eps_list', [float(e) for e in s['eps_list']], '', None, src + ':settings.eps_list',
      'permittivities of the series before refinement (eps = 1: no continuum); one refinement point is added per '
      'geometry at the interpolated switch')
    csv_rows['rev_elec_eps_series.csv'] = (
        ['series', 'case', 'geometry', 'eps', 'f_(eps-1)/eps', 'refinement', 'E_Ha', 'HOMO_eV', 'HOMO_on',
         'HOMO_picrate_pct', 'LUMO_eV', 'LUMO_on', 'LUMO_picrate_pct', 'LUMO_picrate_pct_lowdin', 'LUMO+1_eV',
         'LUMO+1_picrate_pct', 'gap_eV', 'lowest_cation_virtual', 'lowest_cation_virtual_eV',
         'lowest_picrate_virtual', 'lowest_picrate_virtual_eV', 'delta_cat_minus_pic_eV', 'highest_cation_occupied',
         'highest_cation_occupied_eV', 'mulliken_picrate_e', 'lowdin_picrate_e', 'dipole_D', 'source'], eps_rows)

    # ------------------------------------------------------------------ 2. Delta-SCF
    dscf_rows = []
    s2max, spread = 0.0, 0.0
    hirsh_n = {}
    for name, job, combo, case, desc in DSCF:
        d, src = result(job)
        check_settings(d, job)
        sm = d['summary'][combo]
        it = d['items']
        n = it[f'{combo}|N']
        plus = [(k, v) for k, v in it.items() if k.startswith(f'{combo}|N+1|') and v['status'] == 'done']
        minus = [(k, v) for k, v in it.items() if k.startswith(f'{combo}|N-1|') and v['status'] == 'done']
        kp, vp = min(plus, key=lambda kv: kv[1]['E_Ha'])
        km, vm = min(minus, key=lambda kv: kv[1]['E_Ha'])
        ea = (n['E_Ha'] - vp['E_Ha']) * HA_EV
        ie = (vm['E_Ha'] - n['E_Ha']) * HA_EV
        assert abs(ea - sm['deltaSCF_EA_eV']) < 1e-6 and abs(ie - sm['deltaSCF_IE_eV']) < 1e-6
        kea, kip = -n['frontier']['LUMO_eV'], -n['frontier']['HOMO_eV']
        assert abs(kea - sm['koopmans_EA_eV']) < 1e-9
        spread = max(spread, (max(v['E_Ha'] for _, v in plus) - min(v['E_Ha'] for _, v in plus)) * HA_EV)
        for _, v in plus + minus:
            s2max = max(s2max, v['S2'])
            assert v['settings']['spin'] == 1
        model_water = combo.endswith('/pcm')
        assert (n['settings']['solvent'] is not None) == model_water
        geo_series = 'gasgeom' if combo.startswith('gas/') else 'pcmgeom'
        sp = {p['eps']: p for p in series[geo_series]}[wat if model_water else 1.0]
        assert abs(sp['E_Ha'] - n['E_Ha']) < 1e-6 and abs(sp['lumo'] + kea) < 1e-4, (name, sp['E_Ha'], n['E_Ha'])
        hp, hm = vp['populations']['hirshfeld'], vm['populations']['hirshfeld']
        hn = n['populations']['hirshfeld']['fragments']['anion']
        hirsh_n[name] = (hn, n['populations'])
        base = f'dscf.{name}'
        what = f'{desc}; B3LYP, vertical (fixed geometry)' + (', continuum fully relaxed for each charge state'
                                                              if model_water else '')
        A(f'{base}.ea_ev', ea, 'eV', case, f"{src}:items['{combo}|N'].E_Ha - items['{kp}'].E_Ha",
          f'Delta-SCF vertical electron affinity E(N) - E(N+1), lowest N+1 doublet; {what}')
        A(f'{base}.ie_ev', ie, 'eV', case, f"{src}:items['{km}'].E_Ha - items['{combo}|N'].E_Ha",
          f'Delta-SCF vertical ionisation energy E(N-1) - E(N); {what}')
        A(f'{base}.koopmans_ea_ev', kea, 'eV', case, f"{src}:items['{combo}|N'].frontier.LUMO_eV", f'-E(LUMO) of N; {what}')
        A(f'{base}.koopmans_ip_ev', kip, 'eV', case, f"{src}:items['{combo}|N'].frontier.HOMO_eV", f'-E(HOMO) of N; {what}')
        A(f'{base}.ea_minus_koopmans_ev', ea - kea, 'eV', case, src, f'Delta-SCF EA minus -E(LUMO); {what}')
        A(f'{base}.ie_minus_koopmans_ev', ie - kip, 'eV', case, src, f'Delta-SCF IE minus -E(HOMO); {what}')
        A(f'{base}.added_electron_picrate_pct', 100 * hp['spin_fragments']['anion'], '%', case,
          f"{src}:items['{kp}'].populations.hirshfeld.spin_fragments.anion",
          f"share of the unpaired spin of the N+1 doublet on the picrate atoms, Hirshfeld (Mulliken "
          f"{100 * vp['populations']['mulliken']['spin_fragments']['anion']:.1f} %, Loewdin "
          f"{100 * vp['populations']['lowdin']['spin_fragments']['anion']:.1f} %); <S^2> {vp['S2']:.4f}; the "
          f"guesses with the extra electron in the cation-localised and in the picrate-localised virtual both "
          f"converged to this state; the neutral LUMO is on the {FRAG[n['frontier']['LUMO_on']]}; {what}")
        A(f'{base}.added_charge_picrate_hirshfeld_e', hp['fragments']['anion'] - hn, 'e', case,
          f"{src}:items['{kp}'].populations.hirshfeld.fragments.anion - items['{combo}|N']...",
          f'change of the Hirshfeld picrate charge from N to N+1 (where the added electron goes); {what}')
        A(f'{base}.hole_picrate_pct', 100 * hm['spin_fragments']['anion'], '%', case,
          f"{src}:items['{km}'].populations.hirshfeld.spin_fragments.anion",
          f"share of the unpaired spin of the N-1 doublet (hole) on the picrate atoms, Hirshfeld (Mulliken "
          f"{100 * vm['populations']['mulliken']['spin_fragments']['anion']:.1f} %); <S^2> {vm['S2']:.4f}; {what}")
        for st, kk, v in (('N', f'{combo}|N', n), ('N+1', kp, vp), ('N-1', km, vm)):
            P = v['populations']
            dscf_rows.append([name, case, combo, st, kk, v['E_Ha'], v.get('S2'),
                              *[P[m]['fragments']['anion'] for m in ('mulliken', 'lowdin', 'hirshfeld')],
                              *[(P[m].get('spin_fragments') or {}).get('anion') for m in ('mulliken', 'lowdin', 'hirshfeld')],
                              ea if st == 'N+1' else (ie if st == 'N-1' else None), src])
    A('dscf.s2_max', s2max, '', None, 'deltascf_gasgeom, deltascf_pcmgeom result.json items[*].S2',
      'largest <S^2> of the eight doublets (0.75 for a pure doublet): spin contamination is small')
    A('dscf.guess_spread_max_ev', spread, 'eV', None, 'deltascf_gasgeom, deltascf_pcmgeom result.json items[*N+1*].E_Ha',
      'largest energy difference between the two N+1 solutions started with the electron on the cation and on the '
      'picrate (same combination): both guesses reach the same state in every case (aufbau occupations, no MOM)')
    csv_rows['rev_elec_deltascf.csv'] = (
        ['combo', 'case', 'geometry/model', 'state', 'item', 'E_Ha', 'S2', 'mulliken_picrate_e', 'lowdin_picrate_e',
         'hirshfeld_picrate_e', 'mulliken_spin_picrate', 'lowdin_spin_picrate', 'hirshfeld_spin_picrate',
         'deltaSCF_EA_or_IE_eV', 'source'], dscf_rows)

    # ------------------------------------------------------------------ 3. charges
    chg_rows = []
    chg = {}
    for tag, job, item, case, desc in CHG:
        d, src = result(job)
        assert d['status'] == 'done'
        it = d['items'][item]
        assert it['basis'] == '6-311++g(d,p)' and it['status'] == 'done'
        base = f'chg.{tag}'
        lab = it['labels'][int(json.load(open(ROOT / 'molecules' / ('btpcz' if 'btpcz' in tag else 'btpc') / 'cases' /
                                               f"{'btpcz' if 'btpcz' in tag else 'btpc'}_{item}.json"))['proton']['H'])]
        chg[tag] = {}
        for m in SCHEMES:
            q = it[m]['fragments']['anion']
            chg[tag][m] = q
            assert abs(q + it[m]['fragments']['cation']) < 2e-5, (tag, m)
            extra = (f"; the electron count integrates to {it['hirshfeld']['n_electrons_integrated']:.5f} on the "
                     f"molecular grid ({it['n_grid_points']} points)" if m == 'hirshfeld' else
                     f"; CM5 = Hirshfeld + the CM5 pair correction (alpha {d['settings']['cm5_alpha']})" if m == 'cm5'
                     else '')
            A(f'{base}.picrate_{m}_e', q, 'e', case, f"{src}:items['{item}'].{m}.fragments.anion",
              f'{SCHEME_NAME[m]} charge summed over the picrate atoms; {desc}; B3LYP/6-311++G(d,p){extra}')
            A(f'{base}.proton_{m}_e', it['proton_H'][m], 'e', case, f"{src}:items['{item}'].proton_H.{m}",
              f'{SCHEME_NAME[m]} charge of the inter-ion hydrogen {lab}; {desc}')
            chg_rows.append([tag, case, m, q, it[m]['fragments']['cation'], it['proton_H'][m], src])
        mw = it['multiwfn']
        A(f'{base}.picrate_hirshfeld_multiwfn_e', mw['fragments_hirshfeld']['anion'], 'e', case,
          f"{src}:items['{item}'].multiwfn.fragments_hirshfeld.anion",
          'Hirshfeld picrate charge from Multiwfn (its built-in sphericalised free-atom densities, raw values) as a '
          f"cross-check of the PySCF implementation; largest per-atom difference {mw['max_abs_diff_hirshfeld']:.4f} e")
        A(f'{base}.picrate_cm5_multiwfn_e', mw['fragments_cm5']['anion'], 'e', case,
          f"{src}:items['{item}'].multiwfn.fragments_cm5.anion",
          f"CM5 picrate charge from Multiwfn; the CM5 corrections agree to {mw['max_abs_diff_cm5_correction']:.1e} e")
        chg_rows.append([tag, case, 'hirshfeld_multiwfn', mw['fragments_hirshfeld']['anion'],
                         mw['fragments_hirshfeld']['cation'], None, src])
        chg_rows.append([tag, case, 'cm5_multiwfn', mw['fragments_cm5']['anion'], mw['fragments_cm5']['cation'], None, src])
        if tag == 'btpcz_gas':
            for m in SCHEMES:
                A(f'{base}.picric_acid_{m}_e', chg[tag][m] + it['proton_H'][m], 'e', case,
                  f"{src}:items['{item}'].{m}.fragments.anion + proton_H.{m}",
                  f'{SCHEME_NAME[m]} net charge of picric acid (picrate atoms plus H6, which sits on O1 at this '
                  'neutral minimum): the charge moved from benzotriazole to picric acid in the neutral pair')
        if 'check_vs_facts' in it:
            assert it['check_vs_facts_pass']
            A(f'{base}.check_vs_published_max_e', max(abs(v['diff']) for v in it['check_vs_facts'].values()), 'e',
              case, f"{src}:items['{item}'].check_vs_facts",
              'largest difference between the Mulliken/Loewdin charges recomputed from the stored checkpoint and the '
              'published facts (' + ', '.join(it['check_vs_facts']) + ')')
            for k, v in it['check_vs_facts'].items():
                assert abs(pub(facts, k) - v['here']) < 1e-6, k
    for m in SCHEMES:
        A(f'chg.pcm_minus_gas_{m}_e', chg['pcm'][m] - chg['gas'][m], 'e', None, 'rev.elec.chg.pcm/gas.picrate_*',
          f'{SCHEME_NAME[m]} picrate charge, IEF-PCM minimum minus gas-phase minimum (P2_1/c)')
    # cross points from the N states of the Delta-SCF jobs (no CM5 there)
    for name, tag in (('gasgeom_water', 'gasgeom_water'), ('pcmgeom_vac', 'pcmgeom_vac'),
                      ('gasgeom_vac', None), ('pcmgeom_water', None)):
        hn, P = hirsh_n[name]
        job = 'deltascf_gasgeom' if name.startswith('gas') else 'deltascf_pcmgeom'
        combo = dict((n_, c_) for n_, _, c_, _, _ in DSCF)[name]
        src = result(job)[1]
        if tag is None:            # diagonal points: the same density as charges_stored, recomputed independently
            t2 = 'gas' if name == 'gasgeom_vac' else 'pcm'
            for m in ('mulliken', 'lowdin', 'hirshfeld'):
                assert abs(P[m]['fragments']['anion'] - chg[t2][m]) < 1e-4, (name, m)
            continue
        case = 'rev_gas' if name.startswith('gas') else 'rev_pcm'
        desc = dict((n_, d_) for n_, _, _, _, d_ in DSCF)[name]
        for m in ('mulliken', 'lowdin', 'hirshfeld'):
            A(f'chg.{tag}.picrate_{m}_e', P[m]['fragments']['anion'], 'e', case,
              f"{src}:items['{combo}|N'].populations.{m}.fragments.anion",
              f'{SCHEME_NAME[m]} picrate charge of the closed-shell pair, {desc} (cross single point; no CM5 here)')
            chg_rows.append([tag, case, m, P[m]['fragments']['anion'], P[m]['fragments']['cation'], None, src])
    A('chg.hirshfeld_field_effect_gasgeom_e', facts['rev.elec.chg.gasgeom_water.picrate_hirshfeld_e']['value']
      - chg['gas']['hirshfeld'], 'e', 'rev_gas', 'rev.elec.chg.gasgeom_water / chg.gas',
      'Hirshfeld picrate charge: gas-phase geometry in IEF-PCM water minus in vacuum (the field alone)')
    A('chg.hirshfeld_geom_effect_vac_e', facts['rev.elec.chg.pcmgeom_vac.picrate_hirshfeld_e']['value']
      - chg['gas']['hirshfeld'], 'e', None, 'rev.elec.chg.pcmgeom_vac / chg.gas',
      'Hirshfeld picrate charge: IEF-PCM geometry minus gas-phase geometry, both in vacuum (the geometry alone)')
    csv_rows['rev_elec_charges.csv'] = (['point', 'case', 'scheme', 'picrate_e', 'cation_e', 'proton_H_e', 'source'],
                                        chg_rows)

    # ------------------------------------------------------------------ 4. TD-DFT
    td_rows = []
    s1 = {}
    nst = set()
    for sysname, job, key, xc, case in TD:
        d, src = result(job)
        st_ = check_settings(d, job, xc)
        assert st_['response'] == 'TDA'
        it = d['items'][key]
        assert it['status'] == 'done' and it['all_converged'] and len(it['states']) == st_['nstates']
        nst.add(st_['nstates'])
        pair = sysname.startswith('pair')
        states = it['states']
        if 'pcm' in sysname:
            assert it['td_solvent']['equilibrium_solvation'] is False
            eps_opt = it['td_solvent']['eps_response']

        def char(s):
            if not pair:
                return 'within picrate' if sysname.startswith('picrate') else 'within benzotriazole'
            om = s['omega']
            if s['CT_number'] > CT_MIN:
                return 'picrate → cation' if om['anion->cation'] >= om['cation->anion'] else 'cation → picrate'
            if om['anion->anion'] > CT_MIN:
                return 'within picrate'
            if om['cation->cation'] > CT_MIN:
                return 'within cation'
            return 'mixed'
        base = f'td.{sysname}.{xc}'
        what = f'{XC_NAME[xc]} TDA, {TD_DESC[sysname]}'
        lb = it['lowest_bright']
        assert states[lb['state'] - 1]['state'] == lb['state']
        strongest = max(states, key=lambda s: s['f'])
        assert strongest['state'] == it['strongest']['state']
        sel = [('s1', states[0], 'lowest singlet S1'),
               ('bright', states[lb['state'] - 1], 'lowest bright singlet (first with f > 0.01, tddft.py)'),
               ('strongest', strongest, f"strongest of the {len(states)} singlets computed")]
        if pair:
            cts = [s for s in states if s['omega']['anion->cation'] > CT_MIN]
            if cts:
                sel.append(('ctstate', cts[0], f'lowest picrate → cation charge-transfer singlet (Omega share '
                                          f'anion->cation > {CT_MIN})'))
        for tagq, s, lab in sel:
            sp_ = f"{src}:items['{key}'].states[{s['state'] - 1}]"
            tr = '; '.join(f"{t['from']}->{t['to']} {t['weight_pct']} %" for t in s['transitions'][:3])
            ctn = (f", CT number {s['CT_number']:.3f} (hole {100 * s['hole']['anion']:.0f} % / electron "
                   f"{100 * s['electron']['anion']:.0f} % on the picrate)" if pair else '')
            A(f'{base}.{tagq}_state', s['state'], '', case, sp_ + '.state', f'{lab}: state number; {what}')
            A(f'{base}.{tagq}_ev', s['E_eV'], 'eV', case, sp_ + '.E_eV', f'{lab}: excitation energy; {what}; {tr}{ctn}')
            A(f'{base}.{tagq}_nm', s['lambda_nm'], 'nm', case, sp_ + '.lambda_nm', f'{lab}: wavelength; {what}')
            A(f'{base}.{tagq}_f', s['f'], '', case, sp_ + '.f', f'{lab}: oscillator strength (length gauge); {what}')
            A(f'{base}.{tagq}_character', char(s), '', case, sp_ + '.omega',
              f'{lab}: character from the Omega matrix (Loewdin partition of the transition density; CT when the '
              f'inter-fragment share exceeds {CT_MIN}); {what}; {tr}')
            if pair:
                A(f'{base}.{tagq}_ct_number', s['CT_number'], '', case, sp_ + '.CT_number',
                  f'{lab}: charge-transfer number (inter-fragment share of Omega); {what}')
        A(f'{base}.max_ev', states[-1]['E_eV'], 'eV', case, f"{src}:items['{key}'].states[-1].E_eV",
          f'highest of the {len(states)} singlets computed (the spectrum is covered up to here); {what}')
        fr = it['frontier']
        if not pair:
            for q in ('HOMO', 'LUMO'):
                A(f'{base}.{q.lower()}_ev', fr[f'{q}_eV'], 'eV', case, f"{src}:items['{key}'].frontier.{q}_eV",
                  f'{q} energy of the ground state; {what}')
        s1[(sysname, xc)] = states[0]['E_eV']
        for s in states:
            td_rows.append([sysname, xc, case, s['state'], s['E_eV'], s['lambda_nm'], s['f'],
                            s.get('CT_number'), (s.get('omega') or {}).get('anion->cation'),
                            (s.get('omega') or {}).get('cation->anion'), (s.get('hole') or {}).get('anion'),
                            (s.get('electron') or {}).get('anion'), char(s),
                            '; '.join(f"{t['from']}->{t['to']} {t['weight_pct']}" for t in s['transitions']), src])
    A('td.eps_optical', eps_opt, '', 'rev_pcm', result('td_pair_b3lyp')[1] + ":items['pair-pcm/b3lyp'].td_solvent.eps_response",
      'optical permittivity of the non-equilibrium linear-response IEF-PCM (PySCF 2.14 default, that of water); the '
      'ground state is equilibrated with eps 78.3553')
    assert len(nst) == 1
    tdsrc = paths.resolve('dft/revision/tddft.py').read_text()
    m = re.search(r"bright = \[s for s in states if s\['f'\] > ([0-9.]+)\]", tdsrc)
    A('td.bright_threshold_f', float(m.group(1)), '', None, 'dft/revision/tddft.py (lowest_bright)',
      'oscillator strength above which tddft.py calls a state bright (its lowest_bright entry)')
    A('td.ct_threshold', CT_MIN, '', None, 'papers/btpc/data/rev/electronic.py CT_MIN',
      'definition used here: a state is charge transfer when the inter-fragment share of its Omega matrix (the CT '
      'number) exceeds this value; the lowest CT state is the first with Omega(picrate->cation) above it')
    A('td.nstates', nst.pop(), '', None, 'td_* result.json settings.nstates', 'singlets per TDA calculation (Tamm-Dancoff approximation)')
    for xc in ('b3lyp', 'camb3lyp'):
        A(f'td.pair_pcm_minus_free_picrate_s1_{xc}_ev', s1[('pair_pcm', xc)] - s1[('picrate_pcm', xc)], 'eV', 'rev_pcm',
          f'rev.elec.td.pair_pcm.{xc}.s1_ev - td.picrate_pcm.{xc}.s1_ev',
          f'{XC_NAME[xc]}: S1 of the pair in IEF-PCM minus S1 of the free picrate in IEF-PCM (both within the picrate): '
          'the shift of the picrate band caused by the partner cation')
    for sysname in ('pair_gas', 'pair_pcm', 'picrate_pcm', 'bta_pcm'):
        A(f'td.{sysname}.cam_minus_b3lyp_s1_ev', s1[(sysname, 'camb3lyp')] - s1[(sysname, 'b3lyp')], 'eV',
          dict((s_, c_) for s_, _, _, _, c_ in TD)[sysname], f'rev.elec.td.{sysname}.camb3lyp/b3lyp.s1_ev',
          f'S1 with CAM-B3LYP minus S1 with B3LYP; {TD_DESC[sysname]}')
    csv_rows['rev_elec_tddft.csv'] = (
        ['system', 'xc', 'case', 'state', 'E_eV', 'lambda_nm', 'f', 'CT_number', 'omega_picrate_to_cation',
         'omega_cation_to_picrate', 'hole_on_picrate', 'electron_on_picrate', 'character', 'transitions_pct', 'source'],
        td_rows)

    # ------------------------------------------------------------------ 5. other functionals (sp8 jobs)
    fx_rows = []
    deltas = {}
    for g, sname in (('gas', 'gasgeom'), ('pcm', 'pcmgeom')):
        for p in series[sname]:
            if p['eps'] in (1.0, wat):
                deltas.setdefault((g, p['eps']), {})['b3lyp'] = p['delta']
                fx_rows.append(['b3lyp', g, p['eps'], p['homo'], p['lumo'], p['gap'], p['lumo_on'], p['lumo_pic'],
                                p['delta'], result('headline_' + sname)[1]])
        e1 = {p['eps']: p for p in series[sname]}
        sw2 = switch_f([(1.0, e1[1.0]['delta']), (wat, e1[wat]['delta'])])
        A(f'fx.b3lyp.{g}geom.switch_eps_2pt', sw2['f'], '', f'rev_{g}', result('headline_' + sname)[1],
          'B3LYP switch eps interpolated linearly in f from the eps = 1 and 78.3553 points only, for a like-for-like '
          'comparison with PBE0 and CAM-B3LYP (the full series gives rev.elec.eps.' + sname + '.switch_eps)')
    for xc, job in FX:
        d, src = result(job)
        s = check_settings(d, job, xc)
        for g in ('gas', 'pcm'):
            case = f'rev_{g}'
            pts = sorted((it for it in d['items'].values() if it['geometry'] == g), key=lambda it: it['eps'])
            assert [it['eps'] for it in pts] == [1.0, wat]
            for it in pts:
                fr = it['frontier']; o = fr['orbitals']; k = f"{g}|eps={it['eps']:g}"
                assert it['settings']['xc'] == xc and d['items'][k] is it
                t = f"fx.{xc}.{g}geom.{etag(it['eps'])}"
                srcp = f"{src}:items['{k}']"
                what = (f"{XC_NAME[xc]}, P2_1/c {'gas-phase' if g == 'gas' else 'IEF-PCM'} minimum, " +
                        ('no continuum' if it['eps'] == 1 else f"IEF-PCM eps {it['eps']:g}"))
                dl = fr['delta_virt_cation_minus_anion_eV']
                deltas.setdefault((g, it['eps']), {})[xc] = dl
                A(f'{t}.homo_ev', fr['HOMO_eV'], 'eV', case, srcp + '.frontier.HOMO_eV', f'HOMO energy; {what}')
                A(f'{t}.lumo_ev', fr['LUMO_eV'], 'eV', case, srcp + '.frontier.LUMO_eV', f'LUMO energy; {what}')
                A(f'{t}.gap_ev', fr['gap_eV'], 'eV', case, srcp + '.frontier.gap_eV', f'HOMO-LUMO gap; {what}')
                A(f'{t}.lumo_on', FRAG[fr['LUMO_on']], '', case, srcp + '.frontier.LUMO_on', f'LUMO location; {what}')
                A(f'{t}.lumo_picrate_pct', o['LUMO']['mulliken']['anion'], '%', case,
                  srcp + '.frontier.orbitals.LUMO.mulliken.anion', f'Mulliken share of the LUMO on the picrate; {what}')
                A(f'{t}.delta_cat_pic_ev', dl, 'eV', case, srcp + '.frontier.delta_virt_cation_minus_anion_eV',
                  f'Delta = E(lowest cation-localised virtual) - E(lowest picrate-localised virtual); {what}')
                assert FRAG[fr['HOMO_on']] == 'picrate'
                fx_rows.append([xc, g, it['eps'], fr['HOMO_eV'], fr['LUMO_eV'], fr['gap_eV'], FRAG[fr['LUMO_on']],
                                o['LUMO']['mulliken']['anion'], dl, src])
            sw = switch_f([(it['eps'], it['frontier']['delta_virt_cation_minus_anion_eV']) for it in pts])
            assert abs(sw['f'] - d['switch'][g]['eps_switch_f']) < 1e-9
            A(f'fx.{xc}.{g}geom.switch_eps_2pt', sw['f'], '', case, src + f":switch['{g}'].eps_switch_f",
              f'{XC_NAME[xc]} switch eps interpolated linearly in f between eps 1 and {wat:g} (two points only; '
              f"the same estimate with B3LYP is rev.elec.fx.b3lyp.{g}geom.switch_eps_2pt)")
            # the CAM-B3LYP ground state of the TD job is the same calculation: check
            if xc == 'camb3lyp':
                tdk = 'pair-gas/camb3lyp' if g == 'gas' else 'pair-pcm/camb3lyp'
                tdf = result('td_pair_camb3lyp')[0]['items'][tdk]['frontier']
                ref = pts[0] if g == 'gas' else pts[-1]
                assert abs(tdf['LUMO_eV'] - ref['frontier']['LUMO_eV']) < 2e-3 and tdf['LUMO_on'] == ref['frontier']['LUMO_on']
    locs = {(r[0], r[1], r[2]): r[6] for r in fx_rows}
    agree = all(locs[(x, g, 1.0)] == 'cation' and locs[(x, g, wat)] == 'picrate'
                for x in ('b3lyp', 'pbe0', 'camb3lyp') for g in ('gas', 'pcm'))
    A('fx.lumo_location_agrees', 'yes' if agree else 'no', '', None, 'rev.elec.fx.*.lumo_on and rev.elec.eps.*.lumo_on',
      'B3LYP, PBE0 and CAM-B3LYP all put the LUMO on the cation at eps 1 and on the picrate at eps 78.3553, at both '
      'P2_1/c minima' if agree else 'the three functionals disagree on the LUMO location at some point')
    spread_d = max(max(v.values()) - min(v.values()) for v in deltas.values())
    assert all(len(v) == 3 for v in deltas.values())
    A('fx.delta_spread_max_ev', spread_d, 'eV', None, 'rev.elec.fx.*.delta_cat_pic_ev and rev.elec.eps.*.delta_cat_pic_ev',
      'largest spread (max - min) of Delta over B3LYP, PBE0 and CAM-B3LYP at any of the four (geometry, eps) points')
    csv_rows['rev_elec_functionals.csv'] = (
        ['xc', 'geometry', 'eps', 'HOMO_eV', 'LUMO_eV', 'gap_eV', 'LUMO_on', 'LUMO_picrate_pct',
         'delta_cat_minus_pic_eV', 'source'], fx_rows)

    if full_build:
        DERIVED.mkdir(exist_ok=True)
        for fn, (hdr, rows) in csv_rows.items():
            write_csv(fn, hdr, rows)
