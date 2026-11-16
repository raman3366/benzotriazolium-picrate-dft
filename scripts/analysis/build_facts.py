"""
Build papers/btpc/data/facts.json -- the facts dossier for the benzotriazolium picrate (BTPC) paper.

Every value is read from a project file or derived from one here with an explicit formula; nothing is typed in
from memory except the provenance notes.  Re-run with (from the repository root; docs/REPRODUCE.md):
    python scripts/analysis/build_facts.py
Paths in this file and in facts.json are workspace paths; scripts/paths.py maps them to this dataset.
Inputs are read-only (molecules/btpc/..., data/r2/molecules/btpc/...).  Outputs:
    papers/btpc/data/facts.json
    papers/btpc/data/derived/mo_fragment_composition_<case>.csv   (Mulliken fragment shares of every MO)
    papers/btpc/data/derived/crystal_hbonds.csv                     (contacts computed from the CIF by symmetry)
No quantum-chemistry calculation is run: the only "new" numbers are geometry, linear algebra on the stored
Hessians / wavefunctions, and statistics on stored grids.
"""
import sys, os, re, json, math, csv, itertools, pathlib
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))   # scripts/ (paths.py)
import paths
ROOT = paths.legacy_root()   # the dataset files under their workspace paths (scripts/paths.py)
os.chdir(ROOT)
sys.path.insert(0, str(ROOT / 'dft'))
import build_cluster as bc          # CIF parser + symmetry helpers used to build the cluster (pure numpy)

OUT = ROOT / 'papers/btpc/data'
(OUT / 'derived').mkdir(parents=True, exist_ok=True)
T = 'molecules/btpc/results/'          # tracked results
RAW = 'data/r2/molecules/btpc/results/'  # raw results archive
CASES = ('btpc_pcm', 'btpc_gas')
PFX = {'btpc_pcm': 'pcm', 'btpc_gas': 'gas'}
HA_EV = 27.211386245988; HA_KCAL = 627.509474; HA_KJ = 2625.499639; BOHR = 0.529177210903
AU2CM = 219474.6313702; AMU2ME = 1822.888486209
MASS = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999}   # masses PySCF used for freq.json (fact-check: these reproduce freq.json to <1e-10 cm-1; 1.00794/12.0107/... left a 0.10 cm-1 residual)

facts = {}
def add(key, value, unit='', case=None, source='', note=''):
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    if isinstance(value, float):
        value = float(f'{value:.15g}')
    assert key not in facts, key
    facts[key] = {'value': value, 'unit': unit, 'case': case, 'source': source, 'note': note}

def jload(p): return json.load(open(p))
def read_xyz(p):
    L = open(p).read().split('\n'); n = int(L[0].split()[0])
    return [l.split()[0] for l in L[2:2 + n]], np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])

case = jload('molecules/btpc/cases/btpc_pcm.json')
case_gas = jload('molecules/btpc/cases/btpc_gas.json')
LABF = case['cif_labels']; LAB = [l.split('#')[0] for l in LABF]; L = {l: i for i, l in enumerate(LAB)}
PIC = sorted(case['fragments']['picrate anion']); CAT = sorted(case['fragments']['benzotriazolium cation'])
SYM, XR = read_xyz(case['xyz'])                 # X-ray cluster (start geometry)
GEOM = {'xray': XR}
for c in CASES:
    GEOM[PFX[c]] = read_xyz(T + c + '/opt.xyz')[1]

def dist(Y, a, b): return float(np.linalg.norm(Y[L[a]] - Y[L[b]]))
def ang(Y, a, b, c):
    v1 = Y[L[a]] - Y[L[b]]; v2 = Y[L[c]] - Y[L[b]]
    return math.degrees(math.acos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))
def plane(Y, idx):
    c = Y[idx].mean(0); u, s, vt = np.linalg.svd(Y[idx] - c); return c, vt[2]

# =====================================================================================  1. compound / crystal
CIF = 'molecules/btpc/input/BTPC.cif'
cif = open(CIF).read()
def cv(tag):
    m = re.search(rf'^{re.escape(tag)}\s+(.+)$', cif, re.M); return m.group(1).strip().strip("'")
def cnum(tag):
    s = cv(tag); return float(s.split('(')[0])
def csu(tag):
    return cv(tag)

add('compound.project_id', 'btpc', '', None, 'molecules/btpc/ (directory and case-file names btpc_gas, btpc_pcm)', 'Project identifier; the CIF itself calls the compound BTZP.')
add('compound.cif_common_name', cv('_chemical_name_common'), '', None, CIF + ':_chemical_name_common', '_chemical_name_systematic is "?" in the CIF.')
add('compound.chemical_name_casefile', case['chemical_name'], '', None, 'molecules/btpc/cases/btpc_pcm.json:chemical_name',
    'Name assigned by dft/curate_cases.py. The CIF places H on N4 and N6 (the two outer triazole N) and none on N5, i.e. the '
    '1H,3H-benzotriazolium cation (systematically 1H-1,2,3-benzotriazol-3-ium); check the name against the source paper before use.')
add('compound.formula_sum', cv('_chemical_formula_sum'), '', None, CIF + ':_chemical_formula_sum')
add('compound.formula_moiety', cv('_chemical_formula_moiety'), '', None, CIF + ':_chemical_formula_moiety', 'picrate anion C6H2N3O7(-), benzotriazolium cation C6H6N3(+)')
add('compound.ccdc', 835950, '', None, CIF + ':_database_code_depnum_ccdc_archive')
add('compound.source_doi', '10.1016/j.ijleo.2012.06.018', '', None, CIF + ':_citation_doi', 'Original crystal-structure paper (Optik). CIF gives _citation_year 2013.')
add('compound.source_year', 2013, '', None, CIF + ':_citation_year')
add('xray.audit_block_doi', cv('_audit_block_doi'), '', None, CIF + ':_audit_block_doi', 'DOI of the CCDC data block')
add('xray.deposition_date', '2011-07-21', '', None, CIF + ':_audit_update_record')
add('xray.mr', cnum('_chemical_formula_weight'), 'g mol-1', None, CIF + ':_chemical_formula_weight',
    'Multiwfn reports 348.22841 Da from its own atomic masses (data/r2/.../logs/mwfn_mep.out); quote the CIF value.')
add('xray.crystal_system', cv('_symmetry_cell_setting').lower(), '', None, CIF + ':_symmetry_cell_setting')
add('xray.space_group', 'P2_1/c', '', None, CIF + ':_symmetry_space_group_name_H-M', 'CIF string "P2(1)/c"; Hall -P 2ybc; space-group No. 14 (standard numbering, not in CIF).')
add('xray.space_group_hall', cv('_symmetry_space_group_name_Hall'), '', None, CIF + ':_symmetry_space_group_name_Hall')
add('xray.symmetry_operations', 'x,y,z; -x,y+1/2,-z+1/2; -x,-y,-z; x,-y-1/2,z-1/2', '', None, CIF + ':_symmetry_equiv_pos_as_xyz',
    'Operation numbers 1-4 in this order are the first part of the n_klm symmetry codes used below.')
for k, t, u in [('a', '_cell_length_a', 'Å'), ('b', '_cell_length_b', 'Å'), ('c', '_cell_length_c', 'Å'),
                ('alpha', '_cell_angle_alpha', 'deg'), ('beta', '_cell_angle_beta', 'deg'), ('gamma', '_cell_angle_gamma', 'deg'),
                ('volume', '_cell_volume', 'Å^3')]:
    add(f'xray.cell_{k}', cnum(t), u, None, f'{CIF}:{t}', f'with s.u.: {csu(t)}')
add('xray.z', int(cnum('_cell_formula_units_Z')), '', None, CIF + ':_cell_formula_units_Z')
add('xray.z_prime', 1, '', None, CIF + ':_atom_site loop', 'derived: one cation + one anion (33 atoms) in the asymmetric unit, Z = 4 with 4 general positions')
add('xray.temperature', cnum('_diffrn_ambient_temperature'), 'K', None, CIF + ':_diffrn_ambient_temperature', 'CIF: 298(2) K (also cell measurement temperature)')
add('xray.radiation', 'Mo Kα', '', None, CIF + ':_diffrn_radiation_type', 'graphite monochromator, fine-focus sealed tube')
add('xray.wavelength', cnum('_diffrn_radiation_wavelength'), 'Å', None, CIF + ':_diffrn_radiation_wavelength')
add('xray.density_calc', cnum('_exptl_crystal_density_diffrn'), 'g cm-3', None, CIF + ':_exptl_crystal_density_diffrn')
add('xray.mu', cnum('_exptl_absorpt_coefficient_mu'), 'mm-1', None, CIF + ':_exptl_absorpt_coefficient_mu')
add('xray.f000', int(cnum('_exptl_crystal_F_000')), '', None, CIF + ':_exptl_crystal_F_000')
add('xray.crystal_description', 'yellow needle, 0.30 x 0.22 x 0.18 mm', '', None, CIF + ':_exptl_crystal_description/colour/size_*')
add('xray.theta_range', '2.23-26.00', 'deg', None, CIF + ':_diffrn_reflns_theta_min/max', 'completeness to theta_max 0.999')
add('xray.reflections_measured', int(cnum('_diffrn_reflns_number')), '', None, CIF + ':_diffrn_reflns_number')
add('xray.reflections_unique', int(cnum('_reflns_number_total')), '', None, CIF + ':_reflns_number_total')
add('xray.reflections_observed', int(cnum('_reflns_number_gt')), '', None, CIF + ':_reflns_number_gt', 'I > 2sigma(I)')
add('xray.r_int', cnum('_diffrn_reflns_av_R_equivalents'), '', None, CIF + ':_diffrn_reflns_av_R_equivalents')
add('xray.n_parameters', int(cnum('_refine_ls_number_parameters')), '', None, CIF + ':_refine_ls_number_parameters')
add('xray.n_restraints', int(cnum('_refine_ls_number_restraints')), '', None, CIF + ':_refine_ls_number_restraints')
add('xray.r1_gt', cnum('_refine_ls_R_factor_gt'), '', None, CIF + ':_refine_ls_R_factor_gt', 'R1 for I > 2sigma(I)')
add('xray.r1_all', cnum('_refine_ls_R_factor_all'), '', None, CIF + ':_refine_ls_R_factor_all')
add('xray.wr2_gt', cnum('_refine_ls_wR_factor_gt'), '', None, CIF + ':_refine_ls_wR_factor_gt')
add('xray.wr2_all', cnum('_refine_ls_wR_factor_ref'), '', None, CIF + ':_refine_ls_wR_factor_ref', 'wR2 on all data (ref)')
add('xray.gof', cnum('_refine_ls_goodness_of_fit_ref'), '', None, CIF + ':_refine_ls_goodness_of_fit_ref')
add('xray.diff_density_max', cnum('_refine_diff_density_max'), 'e Å-3', None, CIF + ':_refine_diff_density_max')
add('xray.diff_density_min', cnum('_refine_diff_density_min'), 'e Å-3', None, CIF + ':_refine_diff_density_min')
add('xray.software', 'SHELXS-97 (solution), SHELXL-97 (refinement); Bruker SMART/SAINT (collection/reduction)', '', None,
    CIF + ':_computing_*', 'The absorption-correction record names CrysAlisPro / SCALE3 ABSPACK (Oxford Diffraction) while _computing_* names Bruker SMART/SAINT: an internal inconsistency of the deposited CIF.')
add('xray.absorption_correction', 'multi-scan, Tmin 0.9595, Tmax 0.9754', '', None, CIF + ':_exptl_absorpt_correction_*')
add('xray.hydrogen_treatment', 'mixed; every H atom has refined coordinates (with s.u.) and a refined Uiso', '', None,
    CIF + ':_refine_ls_hydrogen_treatment, _atom_site loop', 'N-H 0.91(4) (N4-H4N) and 0.83(4) (N6-H6N) Å; C-H 0.81(3)-0.96(4) Å: X-ray H positions are short, as usual.')
add('xray.hbond_table_deposited', 'none', '', None, CIF, 'The CIF has no _geom_hbond_* loop. All hydrogen-bond geometries below are computed here from the deposited coordinates and symmetry (derived).')

# ---- crystal hydrogen bonds by symmetry (derived)
cell, atoms, ops = bc.read_cif(CIF)
M = bc.frac_to_cart_matrix(cell)
els0 = [a['el'] for a in atoms]; labs0 = [a['label'] for a in atoms]; F0 = np.array([a['f'] for a in atoms]); X0 = F0 @ M.T
bs0 = bc.bonds(els0, X0); hpar = {}
for i, j in bs0:
    if els0[i] == 'H': hpar[i] = j
    if els0[j] == 'H': hpar[j] = i
def opstr(op, sh):
    R, t = bc.parse_op(op); out = []
    for r in range(3):
        terms = ''.join(('+' if R[r, k] > 0 else '-') + 'xyz'[k] for k in range(3) if R[r, k])
        tt = t[r] + sh[r]
        from fractions import Fraction
        f = Fraction(tt).limit_denominator(4)
        s = terms.lstrip('+') + ('' if f == 0 else (f'+{f}' if f > 0 else f'{f}'))
        out.append(s)
    return ', '.join(out)
xhb = []
for oi, op in enumerate(ops):
    R, t = bc.parse_op(op)
    for sh in itertools.product(range(-2, 3), repeat=3):
        Xf = (F0 @ R.T + t + np.array(sh)) @ M.T
        code = f'{oi + 1}_{"".join(str(s + 5) for s in sh)}'
        for h, d in hpar.items():
            for k in range(len(atoms)):
                if els0[k] not in ('N', 'O') or (code == '1_555' and k in (h, d)): continue
                dha = float(np.linalg.norm(X0[h] - Xf[k]))
                if dha < 2.9:
                    a_ = bc.angle(X0[d], X0[h], Xf[k])
                    if a_ > 100:
                        xhb.append(dict(D=labs0[d], H=labs0[h], A=labs0[k], code=code, op=opstr(op, sh),
                                        DH=float(np.linalg.norm(X0[d] - X0[h])), HA=dha, DA=float(np.linalg.norm(X0[d] - Xf[k])), angle=a_))
xhb.sort(key=lambda r: (r['D'][0] != 'N', r['HA']))
with open(OUT / 'derived/crystal_hbonds.csv', 'w') as f:
    w = csv.writer(f); w.writerow(['D', 'H', 'A', 'symmetry_code_of_A', 'symmetry_operation_of_A', 'D-H_A', 'H...A_A', 'D...A_A', 'D-H...A_deg', 'kind'])
    for r in xhb:
        kind = 'intramolecular' if r['code'] == '1_555' else 'intermolecular'
        w.writerow([r['D'], r['H'], r['A'], r['code'], r['op'], f"{r['DH']:.3f}", f"{r['HA']:.3f}", f"{r['DA']:.3f}", f"{r['angle']:.1f}", kind])
k = 0
for r in xhb:
    k += 1
    kind = 'intramolecular (picrate), S(5)-type C-H...O(nitro) contact with a near-100 deg angle; not a hydrogen bond by the usual angle criterion' \
        if r['code'] == '1_555' else 'intermolecular'
    in_cluster = (r['D'], r['A'], r['code']) in {('N4', 'O1', '4_465'), ('N4', 'O7', '4_465'), ('C11', 'O2', '4_465')}
    base = f"xray.hb{k:02d}_{r['D']}_{r['H']}_{r['A']}"
    add(base, f"{r['D']}-{r['H']}...{r['A']}", '', None, CIF + ' (atom_site + symmetry)',
        f"derived; acceptor at symmetry {r['code']} = ({r['op']}); donor at x,y,z; {kind}; "
        f"{'present in the DFT cluster' if in_cluster else 'not in the DFT cluster'}")
    add(base + '_dh', r['DH'], 'Å', None, CIF, 'derived from deposited coordinates')
    add(base + '_ha', r['HA'], 'Å', None, CIF, 'derived')
    add(base + '_da', r['DA'], 'Å', None, CIF, 'derived')
    add(base + '_angle', r['angle'], 'deg', None, CIF, 'derived')
add('xray.hb_search_criteria', 'H...A < 2.9 Å and D-H...A > 100 deg, acceptors N/O, all images within +-2 cells', '', None,
    'papers/btpc/data/build_facts.py', 'The cluster builder (dft/build_cluster.py) uses stricter criteria: H...A < 2.65 Å, angle > 130 deg, heavy-atom clash < 2.40 Å.')
# asymmetric-unit closest contact and stacking
frs = bc.fragments(els0, X0)
dmin = min((float(np.linalg.norm(X0[i] - X0[j])), labs0[i], labs0[j]) for i in frs[0] for j in frs[1])
add('xray.asu_closest_interion_contact', f'{dmin[1]}...{dmin[2]}', '', None, CIF,
    'derived; the two moieties of the asymmetric unit are not hydrogen bonded (N5 carries no H), which is why the cluster was rebuilt from a symmetry image')
add('xray.asu_closest_interion_distance', dmin[0], 'Å', None, CIF, 'derived; matches the 2.99 Å quoted in the case-file notes')
# triazole stacking with inversion image 3_565
tri0 = [labs0.index(x) for x in ('N4', 'N5', 'N6', 'C7', 'C12')]
R3, t3 = bc.parse_op(ops[2]); Xi = (F0 @ R3.T + t3 + np.array([0, 1, 0])) @ M.T
c1, n1 = plane(X0, tri0); c2, n2 = plane(Xi, tri0)
add('xray.stacking_triazole_cg_cg', float(np.linalg.norm(c1 - c2)), 'Å', None, CIF,
    'derived; centroid-centroid distance between the triazole ring and its inversion image (3_565 = -x, 1-y, -z); rings parallel')
add('xray.stacking_triazole_perp', abs(float(np.dot(c2 - c1, n1))), 'Å', None, CIF, 'derived; perpendicular separation of the two parallel triazole planes')
add('xray.stacking_note', 'no picrate...cation or picrate...picrate ring centroid pair closer than 4.3 Å', '', None, CIF, 'derived (search over all images)')

# ---- CIF bond lengths/angles with s.u. for the chemically important bonds
geom_bonds = {}
for m in re.finditer(r'^(\w+)\s+(\w+)\s+(\d+\.\d+\(\d+\))\s+\.\s+\?', cif, re.M):
    geom_bonds[(m.group(1), m.group(2))] = m.group(3)
for (a, b) in [('O1', 'C1'), ('N4', 'N5'), ('N6', 'N5'), ('N4', 'C12'), ('N6', 'C7'), ('N4', 'H4N'), ('N6', 'H6N'),
               ('C6', 'C1'), ('C1', 'C2'), ('N1', 'O2'), ('N1', 'O3'), ('N2', 'O4'), ('N2', 'O5'), ('N3', 'O6'), ('N3', 'O7'),
               ('C2', 'N1'), ('C4', 'N2'), ('N3', 'C6')]:
    s = geom_bonds[(a, b)]
    add(f'xray.bond_{a}_{b}', float(s.split('(')[0]), 'Å', None, CIF + f':_geom_bond {a} {b}', f'CIF value with s.u. {s}')
geom_angles = {}
for m in re.finditer(r'^(\w+)\s+(\w+)\s+(\w+)\s+(\d+\.?\d*\(\d+\))\s+\.\s+\.\s+\?', cif, re.M):
    geom_angles[(m.group(1), m.group(2), m.group(3))] = m.group(4)
for t_ in [('C2', 'C1', 'C6'), ('N4', 'N5', 'N6'), ('N5', 'N4', 'C12'), ('N5', 'N6', 'C7')]:
    s = geom_angles[t_]
    add(f'xray.angle_{"_".join(t_)}', float(s.split('(')[0]), 'deg', None, CIF + f':_geom_angle {" ".join(t_)}', f'CIF value with s.u. {s}')

# =====================================================================================  2. structural model
add('model.n_atoms', len(SYM), '', None, case['xyz'])
add('model.n_electrons', 178, '', None, RAW + 'btpc_pcm/logs/mwfn_mep.out', 'Multiwfn: 178 electrons (89 alpha, 89 beta); fmo.json n_electrons 178, n_occ 89')
add('model.charge', case['charge'], 'e', None, 'molecules/btpc/cases/btpc_pcm.json:charge', 'neutral ion pair: picrate(-1) + benzotriazolium(+1)')
add('model.multiplicity', 1, '', None, 'molecules/btpc/cases/btpc_pcm.json:spin', 'spin = 0 (2S), closed-shell restricted KS')
add('model.fragments', 'picrate anion: 18 atoms (atom indices 0-17, C6H2N3O7); benzotriazolium cation: 15 atoms (indices 18-32, C6H6N3)', '',
    None, 'molecules/btpc/cases/btpc_pcm.json:fragments')
add('model.atom_order', 'CIF atom order within each ion, picrate first (indices 0-17), then cation (18-32); the CIF itself interleaves the two ions (O1, N3, O7, N4, C6, ...). Cation labels carry the internal suffix #4_666 (never show it in the paper)', '', None, 'molecules/btpc/cases/btpc_pcm.json:cif_labels')
add('model.symmetry_picrate', '1_555 (x, y, z)', '', None, 'molecules/btpc/cases/btpc_pcm.json:notes')
add('model.symmetry_cation', '4_666 (x+1, -y+1/2, z+1/2)', '', None, 'molecules/btpc/cases/btpc_pcm.json:cif_labels / dft/build_cluster.py',
    'derived: operation 4 (x, -y-1/2, z-1/2) plus lattice shift (1,1,1). Equivalent, seen from a cation at x,y,z: the picrate at 4_465 = (x-1, -y+1/2, z-1/2).')
add('model.defining_hbond', 'N4-H4N...O1 (cation N-H to picrate phenolate O)', '', None, 'molecules/btpc/cases/btpc_pcm.json:hbonds',
    'Chosen by dft/build_cluster.py (exhaustive search over symmetry images, H...A < 2.65 Å, angle > 130 deg); the bifurcated N4-H4N...O7 (123 deg) is also inside the cluster but below the builder angle cutoff.')
add('model.missing_hbond', 'N6-H6N...O6 (to a second, inversion-related picrate at 1-x,1-y,1-z) is not represented', '', None,
    'molecules/btpc/cases/btpc_pcm.json:notes', 'In the cluster N6-H6N is a free N-H; see xray.hb for its crystal geometry. "Inversion-related" refers to the asymmetric-unit frame (donor cation at x,y,z): that picrate is the inversion image 3_666 of the asymmetric-unit picrate. Relative to the picrate actually in the cluster it is the 2_1-screw image 2_746 (2-x, y-1/2, 3/2-z), H6N...O6 2.103 Å (fact-check, derived).')
add('model.n_basis', 606, '', None, RAW + 'btpc_pcm/logs/frequencies.log', 'also fmo.json n_basis, Multiwfn "Basis functions: 606", molden 606 functions')
add('model.basis_per_element', 'C, N, O: 22 functions each (5s4p1d, spherical d); H: 7 functions (4s1p)', '', None, RAW + 'btpc_pcm/final.molden:[GTO]',
    'derived by counting shells in the molden; 25 heavy atoms x 22 + 8 H x 7 = 606')
add('model.basis_overrides', 'none', '', None, 'molecules/btpc/cases/btpc_pcm.json', 'no basis_overrides key (unlike PZC, which uses def2-TZVP on Zn)')
add('pcm.n_lindep_removed', 3, '', 'btpc_pcm', RAW + 'btpc_pcm/logs/frequencies.log',
    '"3 small eigenvectors of overlap matrix removed" at the final PCM geometry (603 MOs; molden padded to 606 with zero vectors). 4 were removed at some earlier PCM frames (optimize.log: 41 x 4, 30 x 3).')
add('gas.n_lindep_removed', 4, '', 'btpc_gas', RAW + 'btpc_gas/logs/frequencies.log', '602 MOs kept (fmo.json mo_energies_eV has 602 entries)')
add('model.preopt_basis_functions', 366, '', None, 'basis-set count: 6-31G(d) with spherical d functions, 25 heavy atoms x 14 + 8 H x 2', 'B3LYP/6-31G(d) pre-optimisation size; derived check: 25 x 14 + 8 x 2 = 366 with spherical d')

# =====================================================================================  3. protocol
add('protocol.functional', 'B3LYP (PySCF xc="b3lyp" = libxc HYB_GGA_XC_B3LYP, ID 402, VWN-RPA local correlation, same definition as Gaussian)', '', None,
    'dft/common.py:XC; PySCF 2.14.0 pyscf/dft/libxc.py (XC_CODES "B3LYP": 402, "B3LYPG": 402 "used by Gaussian")',
    'No B3LYP_WITH_VWN5 override exists anywhere in the repository.')
add('protocol.basis', '6-311++G(d,p), spherical (5d) functions, all atoms', '', None, 'dft/common.py:BASIS; final.molden [5d]')
add('protocol.density_fitting', 'yes (RI-JK via mf.density_fit())', '', None, 'dft/common.py:make_mf')
add('protocol.aux_basis', 'PySCF automatic default (not logged)', '', None, 'dft/common.py:make_mf; PySCF 2.14.0 pyscf/df/addons.py:make_auxbasis',
    'INFERRED, NOT VERIFIED FROM A RUN LOG (SCF verbosity 3 hides the choice). Per the PySCF 2.14.0 source, Pople 6-311++G names are reduced to "6311++g", which DEFAULT_AUXBASIS maps to aug-cc-pVTZ-JKFIT for hybrid functionals, unless BSE metadata supplies a JK-fit set first. Write "the default auxiliary basis chosen by PySCF" unless confirmed.')
add('protocol.grid_level', 4, '', None, 'verify_grid.log (dft/verify_grid.py run, not deposited); run_case.sh (GRID_LEVEL default 4)',
    'Known issue: every RESULTS.md says "grid level 3" - wrong. verify_grid.py matched grid 4 for btpc_gas (stored E reproduced to 0.00 microHa; grid 3 differs by 1.5 microHa). btpc_pcm was NOT in the verify run; grid 4 for it rests only on the run_case.sh defaults. The "grid 4" in the btpc_pcm opt.xyz comment is a literal string written by dft/finalize_opt.py, not a recorded setting, so it is not independent evidence.')
add('protocol.grid_points', 'H: 60 radial x 434 angular (Lebedev order 35); C, N, O: 90 radial x 590 angular (order 41); Treutler-Ahlrichs radial, Becke partition, NWChem pruning',
    '', None, 'PySCF 2.14.0 pyscf/dft/gen_grid.py (RAD_GRIDS/ANG_ORDER row for level 4)', 'from the PySCF source, not from a run log')
add('protocol.solvation_model', 'IEF-PCM', '', 'btpc_pcm', 'dft/common.py:make_mf (with_solvent.method = "IEF-PCM")')
add('protocol.solvent', 'water', '', 'btpc_pcm', 'molecules/btpc/cases/btpc_pcm.json:solvent')
add('protocol.solvent_eps', 78.3553, '', 'btpc_pcm', 'dft/common.py:SOLVENT_EPS["water"]')
add('protocol.pcm_lebedev_order', 29, '', 'btpc_pcm', 'run_case.sh (PCM_LEBEDEV_ORDER default 29)',
    '302 points per atomic sphere. dft/common.py (PCM_LEBEDEV_ORDER default) still says order 17 (the 3ATB Mac setting) - superseded. verify_grid.log confirms Lebedev 29 for 3atb_pair_pcm and pzc_pcm; btpc_pcm itself was not verified.')
add('protocol.pcm_radii', 'modified Bondi radii (H 1.10, C 1.70, N 1.55, O 1.52 Å) scaled by 1.2 -> H 1.32, C 2.04, N 1.86, O 1.824 Å', '', 'btpc_pcm',
    'PySCF 2.14.0 pyscf/solvent/pcm.py (vdw_scale 1.2, modified_Bondi) and pyscf/data/radii.py (VDW)', 'PySCF defaults (not overridden in dft/common.py); from source, not logged')
add('protocol.pcm_discretisation', 'SWIG (switching-Gaussian) surface, no probe radius, electrostatic term only (no cavitation/dispersion/repulsion)', '', 'btpc_pcm',
    'PySCF 2.14.0 pyscf/solvent/pcm.py (surface_discretization_method "SWIG", r_probe 0.0)')
add('protocol.scf_conv_tol_opt_gas', 1e-9, 'Ha', 'btpc_gas', 'dft/optimize.py (conv_tol 1e-9 gas / 1e-8 solvent)')
add('protocol.scf_conv_tol_opt_pcm', 1e-8, 'Ha', 'btpc_pcm', 'dft/optimize.py')
add('protocol.scf_conv_tol_final', 1e-10, 'Ha', None, 'dft/frequencies.py, dft/fmo.py, dft/_fd_worker.py', 'single points for frequencies, FMO/analysis wavefunction and every displaced gradient')
add('protocol.scf_max_cycle', 150, '', None, 'dft/common.py:make_mf')
add('protocol.optimizer', 'geomeTRIC 1.1.1 via pyscf.geomopt.geometric_solver, TRIC coordinates: 99 (= 3N) delocalised internals built from 170 primitives (34 bonds, 35 angles, 17 out-of-plane, 72 dihedrals, translation + rotation of each ion), trust radius 0.1 (max 0.3) Å, maxiter 300',
    '', None, RAW + 'btpc_pcm/logs/optimize.log (banner and "Optimization Info")')
add('protocol.opt_criteria', '|dE| < 1.0e-6 Ha, RMS grad < 3.0e-4, max grad < 4.5e-4 Ha/Bohr, RMS disp < 1.2e-3, max disp < 1.8e-3 Å (all five required)', '', None,
    RAW + 'btpc_pcm/logs/optimize.log ("Convergence Criteria")', 'geomeTRIC defaults, equal to Gaussian default thresholds. btpc_gas met all five; btpc_pcm did not (see pcm.opt_*).')
add('protocol.preopt', 'B3LYP/6-31G(d) with loose criteria (dE 1e-6, grms 1.7e-3, gmax 2.5e-3, drms 6.7e-3, dmax 1.0e-2) - pipeline default', '', None,
    'dft/optimize.py (PREOPT_BASIS, PRE_CONV)', 'For BTPC the final grid-4 runs skipped this stage and started from preopt.xyz files carried over from earlier grid-3 runs on the OVH VM (see gas.opt_start, pcm.opt_start); the 6-31G(d) stage happened in those earlier runs (not in the archive).')
add('protocol.hessian_method', 'central finite differences of analytic gradients (6N = 198 gradients), each displaced SCF started from the converged density', '', None,
    'dft/frequencies.py; dft/_fd_worker.py', 'H = (g(+h) - g(-h)) / 2h, symmetrised 0.5(H + H^T)')
add('protocol.hessian_projection', 'translations and rotations projected out in mass-weighted Cartesians (QR-orthonormalised projector), own routine dft/vibrations.py', '',
    None, 'dft/vibrations.py; freq.json tr_residual_cm-1', 'residual 2.7e-5 (PCM), 2.8e-5 (gas) cm-1. PySCF harmonic_analysis listed the PCM imaginary mode as 0+3.2j and counted it as real (n_imaginary_pyscf = 0).')
add('protocol.frequency_scaling', 'none (harmonic, unscaled)', '', None, 'dft/frequencies.py', 'No IR/Raman intensities were computed.')
add('protocol.thermo_T', 298.15, 'K', None, 'dft/frequencies.py (thermo.thermo(mf, freq, 298.15, 101325))')
add('protocol.thermo_p', 101325, 'Pa', None, 'dft/frequencies.py', 'PySCF rigid-rotor/harmonic-oscillator thermochemistry, no low-frequency correction')
add('protocol.fmo_descriptors', 'IP = -E_HOMO; EA = -E_LUMO; chi = (IP+EA)/2; mu = -chi; eta = (IP-EA)/2; S = 1/(2 eta); omega = mu^2/(2 eta)', '', None,
    'dft/fmo.py', 'Koopmans-type. Softness convention S = 1/(2 eta) (some papers use S = 1/eta, which doubles S). Values from the PCM wavefunction include the continuum reaction field.')
add('protocol.fragment_composition_method', 'Mulliken partition of each MO over the two fragments: sum over fragment AOs of c_mu (S c)_mu', '', None, 'dft/fmo.py:frag_share')
add('protocol.cubes', 'HOMO/LUMO/density 80^3, ESP 64^3, margin 4.5 Bohr (2.38 Å) beyond the outermost nuclei (PySCF cubegen); render isovalues: orbitals +-0.02 a.u., MEP on 0.001 a.u. density, NCI RDG = 0.5 coloured -0.035..+0.02 a.u.',
    '', None, 'dft/fmo.py; dft/render.py; cube headers', 'PySCF cubegen takes margin in Bohr: every cube origin equals the smallest nuclear coordinate minus 4.5 Bohr on each axis (e.g. PCM x: -8.1537 - 4.5 = -12.6537 Bohr). The "4 A"/"4.5 Å" wording in dft/fmo.py is wrong. The 0.001 a.u. surface extends at most ~1.8 Å beyond the nuclei (Multiwfn vertex ranges), so it is not clipped.')
add('protocol.mep_surface', 'Multiwfn main function 12: electron-density isosurface 0.001 a.u., grid spacing 0.25 Bohr, marching tetrahedra; ESP evaluated exactly at surface vertices', '', None,
    RAW + 'btpc_pcm/logs/mwfn_mep.out (menu lines "iso: 0.00100", "Spacing of grid points ... 0.250000")', 'Colour range of the rendered MEP = surface min..max (mep_range.json). Reported MEP numbers do not depend on the 64^3 ESP cube (render only).')
add('protocol.mep_esp_code_pcm', 'iESPcode 2 (LIBRETA fast)', '', 'btpc_pcm', RAW + 'btpc_pcm/mwfn_home/settings.ini; logs/mwfn_mep.out', 'Multiwfn exited with code 59 after writing complete output; accepted by dft/mwfn.py')
add('protocol.mep_esp_code_gas', 'iESPcode 1 (Multiwfn original ESP code)', '', 'btpc_gas', RAW + 'btpc_gas/mwfn_home/settings.ini; data/r2/logs/btpc_mep_retry.out',
    'LIBRETA fast (2) and slow (3) both segfaulted (exit 174) on the gas wavefunction; code 1 ran 863 s, exit 59, output complete; statistics then parsed separately (STAGES 2026-09-17 00:31).')
add('protocol.mep_note', 'The ESP is that of the (PCM-polarised, for btpc_pcm) molecular charge density alone; the PCM surface charges are not part of the mapped potential', '', None,
    'dft/mwfn.py (Multiwfn reads only final.molden)', 'derived from how the analysis was run')
add('protocol.nci_settings', 'Multiwfn 20 -> 1 (NCI), "High-quality grid" covering the whole system with 1.5 Bohr extension; RDG_maxrho = 0.1 a.u. (points with rho >= 0.1 set to RDG 100); promolecular approximation not used (full wavefunction)',
    '', None, RAW + 'btpc_pcm/logs/mwfn_nci.out; mwfn_home/settings.ini', 'Default RDG_maxrho is 0.05; raised to 0.1 (setup_linux.sh) so the short N-H...O bond survives (its critical-point density is 0.087 a.u. in the gas phase, gas.bcp_H4N_O1_rho).')
add('protocol.nci_grid_pcm', '126 x 148 x 95 = 1,771,560 points, spacing 0.158553 Bohr (0.0839 Å)', '', 'btpc_pcm', RAW + 'btpc_pcm/logs/mwfn_nci.out; rdg.cube header')
add('protocol.nci_grid_gas', '125 x 147 x 96 = 1,764,000 points, spacing 0.153832 Bohr (0.0814 Å)', '', 'btpc_gas', RAW + 'btpc_gas/logs/mwfn_nci.out; rdg.cube header')
add('protocol.nci_spike_method', 'dft/mwfn.py:nci_post: sign(l2)rho binned from -0.200 to +0.050 a.u. in 0.001 bins; a bin is a candidate if it holds >20 points and its minimum RDG < 0.3; contiguous candidate bins are merged and represented by the bin with the lowest RDG; type label: < -0.015 "attractive (H-bond)", < +0.010 "van der Waals", otherwise "steric repulsion"',
    '', None, 'dft/mwfn.py', 'The label is purely a threshold on sign(l2)rho; assignments to atoms are in *.nci_region_* facts (derived here).')
add('protocol.elf_lol_settings', 'Multiwfn 4 -> 9 (ELF) / 10 (LOL), plane through three atoms, 300 x 300 points, aug2D = 6.0 Bohr, Becke definition (ELFLOL_type 0, ELF_addminimal 1)',
    '', None, RAW + 'btpc_pcm/mwfn_elf_*.in; mwfn_home/settings.ini; logs/mwfn.log')
add('protocol.elf_plane_picrate', 'picrate ring plane through C6, C4, C2 (1-based Multiwfn atoms 4, 9, 6)', '', None, 'dft/curate_cases.py; data/r2/.../mwfn_elf_picrate_ring.in')
add('protocol.elf_plane_cation', 'benzotriazolium plane through N4, N5, C7 (1-based atoms 19, 23, 21)', '', None, 'dft/curate_cases.py; data/r2/.../mwfn_elf_benzotriazolium_ring.in',
    'triazole and benzo rings are coplanar (9-atom rms deviation 0.0009 Å PCM, 0.0013 Å gas), so one map covers the whole cation')
add('protocol.host', 'IBM Cloud mx3d-48x480 VM (48 vCPU / 24 physical Intel Sapphire Rapids cores, 472 GB RAM, local NVMe scratch)', '', None,
    'IBM Cloud instance profile mx3d-48x480 (public profile specification); run_case.sh comment on the scratch volume; run paths <run-dir> and <scratch-dir> in the logs; install log "threads 48" (not deposited)')
add('software.pyscf', '2.14.0', '', None, RAW + 'btpc_pcm/final.molden header "made by pyscf v[2.14.0]"; density.cube header; install log (not deposited)')
add('software.geometric', '1.1.1', '', None, RAW + 'btpc_pcm/logs/optimize.log ("geomeTRIC started. Version: 1.1.1"); install log (not deposited)')
add('software.multiwfn', '2026.9.13 (official Linux noGUI binary)', '', None, RAW + 'btpc_pcm/logs/mwfn_mep.out banner; install log (not deposited)', 'Multiwfn asks to cite J. Comput. Chem. 33, 580 (2012) and J. Chem. Phys. 161, 082503 (2024); MEP module J. Mol. Graph. Model. 38, 314 (2012); LIBRETA ESP PCCP 23, 20323 (2021); MPI Carbon 171, 514 (2021); NCI original paper J. Am. Chem. Soc. 132, 6498 (2010) plus the reviews Comprehensive Computational Chemistry vol. 2, 240 (2024, DOI 10.1016/B978-0-12-821978-2.00076-3) and Angew. Chem. Int. Ed. e202504895 (2025) (as printed in the logs; the NCI entries are in logs/mwfn_nci.out)')
add('software.paraview', '6.1.1 (pvbatch, headless)', '', None, 'install log (not deposited); environment/software.md (ParaView-6.1.1-MPI-Linux-Python3.12-x86_64)',
    'The install log shows a libgomp.so.1 load error for pvbatch at install time; the renders later succeeded (logs/render.log "ok"), so it was fixed outside the log.')
add('software.python', '3.11.16', '', None, 'environment/compute-env.yml (conda package list of the compute environment)')
add('software.numpy_scipy_matplotlib', 'NumPy 2.4.6, SciPy 1.17.1, matplotlib 3.11.1', '', None, 'environment/compute-env.yml (conda package list of the compute environment)')
add('software.libxc', 'version not recorded', '', None, 'environment/compute-env.yml; install log (searched, absent; not deposited)', 'PySCF 2.14.0 conda-forge build; libxc version absent from the archive')

# =====================================================================================  4. per case
def rms(v): return float(np.sqrt(np.mean(np.square(v))))
nb = {i: set() for i in range(len(SYM))}
COV = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66}
BONDS = [(i, j) for i in range(len(SYM)) for j in range(i + 1, len(SYM)) if np.linalg.norm(XR[i] - XR[j]) < 1.25 * (COV[SYM[i]] + COV[SYM[j]])]
for i, j in BONDS: nb[i].add(j); nb[j].add(i)
ANGLES = [(i, j, k) for j in range(len(SYM)) for i in nb[j] for k in nb[j] if i < k]
heavy = lambda t: all(SYM[a] != 'H' for a in t)
def dij(Y, i, j): return float(np.linalg.norm(Y[i] - Y[j]))
def aijk(Y, i, j, k):
    v1 = Y[i] - Y[j]; v2 = Y[k] - Y[j]; return math.degrees(math.acos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))
def kabsch(A, B):
    a = A - A.mean(0); b = B - B.mean(0); U, S_, Vt = np.linalg.svd(a.T @ b); dd = np.sign(np.linalg.det(Vt.T @ U.T))
    Rm = Vt.T @ np.diag([1, 1, dd]) @ U.T; return float(np.sqrt((((Rm @ a.T).T - b) ** 2).sum(1).mean()))
pring = [L[a] for a in ('C1', 'C2', 'C3', 'C4', 'C5', 'C6')]
cat9 = [L[a] for a in ('N4', 'N5', 'N6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12')]
NITRO = {'N1': ('N1', 'O2', 'O3'), 'N2': ('N2', 'O4', 'O5'), 'N3': ('N3', 'O6', 'O7')}
KEYB = [('O1', 'C1'), ('C1', 'C2'), ('C1', 'C6'), ('C2', 'C3'), ('C3', 'C4'), ('C4', 'C5'), ('C5', 'C6'), ('C2', 'N1'), ('C4', 'N2'), ('C6', 'N3'),
        ('N1', 'O2'), ('N1', 'O3'), ('N2', 'O4'), ('N2', 'O5'), ('N3', 'O6'), ('N3', 'O7'), ('N4', 'N5'), ('N5', 'N6'), ('N4', 'C12'), ('N6', 'C7'),
        ('C7', 'C12'), ('C7', 'C8'), ('C8', 'C9'), ('C9', 'C10'), ('C10', 'C11'), ('C11', 'C12'), ('N4', 'H4N'), ('N6', 'H6N')]
KEYA = [('C2', 'C1', 'C6'), ('N4', 'N5', 'N6'), ('N5', 'N4', 'C12'), ('N5', 'N6', 'C7'), ('C1', 'O1', 'H4N'), ('C1', 'O1', 'N4')]
INTRA_CONTACTS = [('N4', 'H4N', 'O1'), ('N4', 'H4N', 'O7'), ('C11', 'H11', 'O2'), ('C11', 'H11', 'O1'), ('C5', 'H5', 'O6'), ('C5', 'H5', 'O5'), ('C3', 'H3', 'O3'), ('C3', 'H3', 'O4')]
OO = [('O1', 'O7'), ('O1', 'O2'), ('N5', 'O1')]

def geom_facts(tag, Y, cs, src):
    """structural facts for one geometry (tag = xray / pcm / gas)"""
    note_x = 'derived from the X-ray cluster molecules/btpc/input/btpc_cluster.xyz (CIF coordinates, symmetry-applied)' if tag == 'xray' else f'derived from {src}'
    for a, b in KEYB:
        add(f'{tag}.r_{a}_{b}', dij(Y, L[a], L[b]), 'Å', cs, src, note_x)
    for t_ in KEYA:
        add(f'{tag}.a_{"_".join(t_)}', aijk(Y, *[L[x] for x in t_]), 'deg', cs, src, note_x)
    for D, H, A in INTRA_CONTACTS:
        base = f'{tag}.hb_{D}_{H}_{A}'
        kind = 'inter-ion' if ((L[D] in PIC) != (L[A] in PIC)) else 'intramolecular (picrate C-H...O(nitro), angle ~90-106 deg: geometric contact, not a hydrogen bond)'
        add(base + '_dh', dij(Y, L[D], L[H]), 'Å', cs, src, f'{note_x}; {kind}')
        add(base + '_ha', dij(Y, L[H], L[A]), 'Å', cs, src, note_x)
        add(base + '_da', dij(Y, L[D], L[A]), 'Å', cs, src, note_x)
        add(base + '_angle', aijk(Y, L[D], L[H], L[A]), 'deg', cs, src, note_x)
    for a, b in OO:
        add(f'{tag}.r_{a}_{b}', dij(Y, L[a], L[b]), 'Å', cs, src, note_x + '; non-bonded contact')
    c1, n1 = plane(Y, pring); c2, n2 = plane(Y, cat9)
    add(f'{tag}.interplanar_picrate_cation', math.degrees(math.acos(min(1, abs(float(np.dot(n1, n2)))))), 'deg', cs, src,
        note_x + '; angle between the least-squares planes of the picrate C6 ring and the 9 ring atoms of the cation')
    add(f'{tag}.centroid_distance_picrate_cation', float(np.linalg.norm(c1 - c2)), 'Å', cs, src, note_x + '; picrate ring centroid to cation 9-atom centroid')
    for nm, (N, Oa, Ob) in NITRO.items():
        cN, nN = plane(Y, [L[N], L[Oa], L[Ob]])
        pos = {'N1': 'ortho, on C2', 'N2': 'para, on C4', 'N3': 'ortho, on C6 (next to O1)'}[nm]
        add(f'{tag}.nitro_twist_{nm}', math.degrees(math.acos(min(1, abs(float(np.dot(nN, n1)))))), 'deg', cs, src,
            note_x + f'; angle between the O-{nm}-O plane and the picrate ring plane ({pos})')
    return c1, n1, c2, n2

gx = geom_facts('xray', XR, None, 'molecules/btpc/input/btpc_cluster.xyz')

freq = {c: jload(T + c + '/freq.json') for c in CASES}
fmo = {c: jload(T + c + '/fmo.json') for c in CASES}
opts = {c: jload(T + c + '/opt_summary.json') for c in CASES}
msurf = {c: jload(T + c + '/mep_surface.json') for c in CASES}

# ---- Hessian re-analysis (normal modes)
def modes_for(c):
    s, xyz = read_xyz(T + c + '/opt.xyz'); x = xyz / BOHR
    H = np.load(RAW + c + '/hessian.npy'); n = len(s); m = np.array([MASS[e] for e in s]) * AMU2ME
    Hm = H.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n); w = np.repeat(m ** -0.5, 3); Hm = Hm * w[:, None] * w[None, :]; Hm = 0.5 * (Hm + Hm.T)
    com = (m[:, None] * x).sum(0) / m.sum(); r = x - com
    TR = np.zeros((6, 3 * n))
    for i in range(n):
        sq = np.sqrt(m[i]); TR[0, 3 * i] = TR[1, 3 * i + 1] = TR[2, 3 * i + 2] = sq
        for k_, ax in enumerate(np.eye(3)): TR[3 + k_, 3 * i:3 * i + 3] = sq * np.cross(ax, r[i])
    q, _ = np.linalg.qr(TR.T); P = np.eye(3 * n) - q @ q.T
    ev, V = np.linalg.eigh(P @ Hm @ P); f = np.sign(ev) * np.sqrt(abs(ev)) * AU2CM
    order = np.argsort(abs(f)); tr = set(order[:6].tolist()); keep = [i for i in np.argsort(f) if i not in tr]
    return s, x, m, f[keep], V[:, keep], q
def rigid_basis(x, m, groups):
    B = []
    for g in groups:
        gm = m[g]; com = (gm[:, None] * x[g]).sum(0) / gm.sum()
        for k_ in range(3):
            v = np.zeros((len(x), 3)); v[g, k_] = np.sqrt(gm); B.append(v.ravel())
        for ax in np.eye(3):
            v = np.zeros((len(x), 3))
            for i in g: v[i] = np.sqrt(m[i]) * np.cross(ax, x[i] - com)
            B.append(v.ravel())
    Q, _ = np.linalg.qr(np.array(B).T); return Q
def nh_modes(x, m, f, V, H, D):
    c_, n_ = plane(x, cat9); h, d = L[H], L[D]; b = x[h] - x[d]; b /= np.linalg.norm(b); ip = np.cross(n_, b); ip /= np.linalg.norm(ip)
    out = {'str': None, 'ip': None, 'oop': []}
    best = {'str': 0, 'ip': 0}
    for k_ in range(len(f)):
        v = V[:, k_].reshape(-1, 3); wH = float((v[h] ** 2).sum())
        if wH < 0.15: continue
        vh = v[h] / np.sqrt(wH); fs, fi, fo = [float(np.dot(vh, e) ** 2) for e in (b, ip, n_)]
        if fs > 0.5 and wH * fs > best['str']: best['str'] = wH * fs; out['str'] = (f[k_], wH, k_)
        if fi > 0.5 and wH * fi > best['ip']: best['ip'] = wH * fi; out['ip'] = (f[k_], wH, k_)
        if fo > 0.5: out['oop'].append((f[k_], wH, k_))
    return out

# ---- orbital composition from the molden (validated overlap)
from math import comb
def parse_molden(p):
    sec = None; atoms_ = []; shells = []; mos = []; cur = None
    for line in open(p):
        s = line.strip()
        if s.startswith('['): sec = s.split(']')[0][1:].lower(); continue
        q = s.split()
        if sec == 'atoms' and len(q) == 6: atoms_.append(np.array([float(v) for v in q[3:6]]))
        elif sec == 'gto':
            if len(q) == 2 and q[1] == '0': cur = int(q[0]) - 1
            elif len(q) == 3 and q[0] in 'spdf': shells.append({'atom': cur, 'l': 'spdf'.index(q[0]), 'exp': [], 'c': []})
            elif len(q) == 2 and shells: shells[-1]['exp'].append(float(q[0].replace('D', 'E'))); shells[-1]['c'].append(float(q[1].replace('D', 'E')))
        elif sec == 'mo':
            if s.startswith('Ene='): mos.append({'E': float(s.split('=')[1]), 'c': []})
            elif s.startswith('Occup='): mos[-1]['occ'] = float(s.split('=')[1])
            elif q and q[0].isdigit(): mos[-1]['c'].append(float(q[1]))
    return atoms_, shells, mos
def dfact(n): return 1 if n <= 0 else n * dfact(n - 2)
def s1d(i, j, PA, PB, p):
    tot = 0.0
    for k1 in range(i + 1):
        for k2 in range(j + 1):
            kk = k1 + k2
            if kk % 2: continue
            tot += comb(i, k1) * comb(j, k2) * PA ** (i - k1) * PB ** (j - k2) * dfact(kk - 1) / (2 * p) ** (kk / 2)
    return tot * math.sqrt(math.pi / p)
def pnorm(a, lx, ly, lz):
    return (2 * a / math.pi) ** 0.75 * (4 * a) ** ((lx + ly + lz) / 2) / math.sqrt(dfact(2 * lx - 1) * dfact(2 * ly - 1) * dfact(2 * lz - 1))
CART = {0: [(0, 0, 0)], 1: [(1, 0, 0), (0, 1, 0), (0, 0, 1)], 2: [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)]}
SPH_D = np.array([[-0.5, -0.5, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1], [math.sqrt(3) / 2, -math.sqrt(3) / 2, 0, 0, 0, 0], [0, 0, 0, math.sqrt(3), 0, 0]])
def overlap(atoms_, shells):
    nc = [len(CART[sh['l']]) for sh in shells]; off = np.cumsum([0] + nc); N = off[-1]; Sc = np.zeros((N, N))
    for I, a in enumerate(shells):
        A = atoms_[a['atom']]
        for J in range(I, len(shells)):
            b = shells[J]; B = atoms_[b['atom']]; AB2 = float(((A - B) ** 2).sum()); blk = np.zeros((nc[I], nc[J]))
            for ea, ca in zip(a['exp'], a['c']):
                for eb, cb in zip(b['exp'], b['c']):
                    p = ea + eb; P = (ea * A + eb * B) / p; K = math.exp(-ea * eb / p * AB2)
                    if K < 1e-14: continue
                    for u, (ax, ay, az) in enumerate(CART[a['l']]):
                        for v, (bx, by, bz) in enumerate(CART[b['l']]):
                            blk[u, v] += ca * cb * pnorm(ea, ax, ay, az) * pnorm(eb, bx, by, bz) * K * s1d(ax, bx, P[0] - A[0], P[0] - B[0], p) * \
                                s1d(ay, by, P[1] - A[1], P[1] - B[1], p) * s1d(az, bz, P[2] - A[2], P[2] - B[2], p)
            Sc[off[I]:off[I + 1], off[J]:off[J + 1]] = blk; Sc[off[J]:off[J + 1], off[I]:off[I + 1]] = blk.T
    Tm = []; aof = []
    for I, sh in enumerate(shells):
        if sh['l'] < 2:
            for k_ in range(nc[I]):
                row = np.zeros(N); row[off[I] + k_] = 1; Tm.append(row); aof.append(sh['atom'])
        else:
            for cb_ in SPH_D:
                row = np.zeros(N); row[off[I]:off[I + 1]] = cb_; Tm.append(row); aof.append(sh['atom'])
    Tm = np.array(Tm); S = Tm @ Sc @ Tm.T; dd = np.sqrt(np.diag(S)); return S / dd[:, None] / dd[None, :], np.array(aof)

def read_cube(p):
    with open(p) as fh:
        fh.readline(); fh.readline()
        nat, *org = fh.readline().split(); org = np.array([float(v) for v in org[:3]])
        ax = [fh.readline().split() for _ in range(3)]; n = [int(a[0]) for a in ax]; sp = np.array([[float(v) for v in a[1:4]] for a in ax])
        for _ in range(abs(int(nat))): fh.readline()
        data = np.fromstring(fh.read(), sep=' ')
    return org, sp, n, data.reshape(n)

def orbname(k, nocc):
    return ('HOMO' if k == nocc - 1 else f'HOMO-{nocc - 1 - k}') if k < nocc else ('LUMO' if k == nocc else f'LUMO+{k - nocc}')

MODES = {}
for c in CASES:
    p = PFX[c]; src_opt = T + c + '/opt.xyz'; Y = GEOM[p]
    fq, fm, op, ms = freq[c], fmo[c], opts[c], msurf[c]
    # ---- energies, optimisation
    add(f'{p}.e_total', op['E_final_Ha'], 'Ha', c, T + c + '/opt_summary.json:E_final_Ha',
        f"freq.json E_Ha {fq['E_Ha']:.10f}, fmo.json E_total_Ha {fm['E_total_Ha']:.10f} (agree to < 1e-9 Ha)")
    add(f'{p}.e_total_ev', op['E_final_Ha'] * HA_EV, 'eV', c, T + c + '/opt_summary.json:E_final_Ha', 'derived: E x 27.211386245988')
    if c == 'btpc_gas':
        add('gas.opt_steps', op['steps'], '', c, T + c + '/opt_summary.json:steps', 'number of energy+gradient evaluations (history entries 0-2); geomeTRIC "Converged!" after step 2')
        add('gas.opt_wall', op['wall_min'], 'min', c, T + c + '/opt_summary.json:wall_min', '4 threads on the 48-vCPU compute VM')
        add('gas.opt_converged', 'yes, all five geomeTRIC/Gaussian criteria', '', c, RAW + c + '/logs/optimize.log',
            'final step: dE -6.7e-7 Ha, grad rms/max 3.50e-5/6.83e-5, disp rms/max 5.61e-4/1.18e-3 Å')
        add('gas.opt_start', 'preopt.xyz "seeded from the converged grid-3 geometry; re-optimised here at grid level 4"', '', c, RAW + c + '/preopt.xyz (comment line)',
            'start already near-stationary at grid 4 (step 0 max grad 1.1e-4 Ha/Bohr), hence 3 evaluations; opt_summary preopt_basis null / preopt_steps 0 for this run')
        add('gas.opt_date', '2026-09-16 12:05-12:19 UTC', '', c, T + c + '/STAGES.md')
    else:
        add('pcm.opt_steps', op['steps'], '', c, T + c + '/opt_summary.json:steps',
            '70 optimisation steps = 71 energy+gradient evaluations (frames 0-70 of opt_traj.xyz). NB: the gas "3 steps" counts evaluations, this counts steps.')
        add('pcm.opt_wall', 'about 13.4 h (2026-09-16 14:37 -> 2026-09-17 ~04:01 UTC), 4 threads', '', c, T + c + '/STAGES.md; opt_summary.json wall_min = NaN',
            'derived from STAGES timestamps and the optimize.log/opt_traj.xyz modification time; ~11.4 min per step')
        add('pcm.opt_converged', 'no - accepted on forces only', '', c, T + c + '/opt_summary.json:converged_on; ' + RAW + c + '/logs/optimize.log',
            'dft/finalize_opt.py took the last frame when forces were inside Gaussian thresholds (rms 1.40e-4, max 3.85e-4 Ha/Bohr). The energy and displacement criteria were NOT met: last step dE -8.6e-6 Ha (> 1e-6), disp rms/max 1.99e-2/3.76e-2 Å (> 1.2e-3/1.8e-3), trust radius growing, energy still falling. The run was stopped during step 71.')
        olog = re.sub(r'\x1b\[[0-9;]*m', '', open(RAW + c + '/logs/optimize.log', errors='replace').read())
        ls = re.findall(r'Step\s+(\d+) : Displace = ([\d.e+-]+)/([\d.e+-]+) \(rms/max\) Trust = ([\d.e+-]+).*?Grad = ([\d.e+-]+)/([\d.e+-]+) \(rms/max\) E \(change\) = (-?[\d.]+) \(([-+\d.e]+)\)', olog)[-1]
        src_log = RAW + c + f'/logs/optimize.log (Step {ls[0]} line)'
        add('pcm.opt_last_step_dE', float(ls[7]), 'Ha', c, src_log, 'E (change) printed by geomeTRIC for the last completed step; criterion 1e-6')
        add('pcm.opt_last_step_disp_rms', float(ls[1]), 'Å', c, src_log, 'criterion 1.2e-3')
        add('pcm.opt_last_step_disp_max', float(ls[2]), 'Å', c, src_log, 'criterion 1.8e-3')
        add('pcm.opt_last_trust_radius', float(ls[3]), 'Å', c, src_log, 'trust radius had been growing for 10 steps (8.7e-4 at step 60)')
        add('pcm.opt_final_grad_rms', float(ls[4]), 'Ha/Bohr', c, src_log, 'geomeTRIC definition; criterion 3.0e-4')
        add('pcm.opt_final_grad_max', float(ls[5]), 'Ha/Bohr', c, src_log, 'geomeTRIC definition; criterion 4.5e-4')
        add('pcm.opt_start', 'preopt.xyz "harvested from the grid-3 IEF-PCM optimisation trajectory on OVH (frame 16)"', '', c, RAW + c + '/preopt.xyz (comment line)',
            'grid-3 run itself started from the X-ray cluster; not in the archive')
        tr = open(RAW + c + '/opt_traj.xyz').read().split('\n')
        e0 = float(tr[1].split('=')[1].split()[0])
        add('pcm.opt_energy_drop_total', (op['E_final_Ha'] - e0) * HA_KCAL, 'kcal mol-1', c, RAW + c + '/opt_traj.xyz (step 0 vs step 70)', f'derived: (E70 - E0) x 627.509474; E0 = {e0:.8f} Ha')
    # ---- proton position / H-bond
    pc = jload(T + c + '/proton_check.json')
    add(f'{p}.proton_transferred', 'no' if not pc['proton_transferred'] else 'yes', '', c, T + c + '/proton_check.json:proton_transferred',
        'criterion: transferred if H4N...O1 < N4-H4N; no covalent bond broken or formed (connectivity_changed false)')
    geom_facts(p, Y, c, src_opt)
    add(f'{p}.nh_minus_ho_asym', dij(Y, L['H4N'], L['O1']) - dij(Y, L['N4'], L['H4N']), 'Å', c, src_opt, 'derived: d(H4N...O1) - d(N4-H4N); 0 would be a centred proton')
    # ---- RMSD vs X-ray
    db = np.array([dij(Y, i, j) - dij(XR, i, j) for i, j in BONDS]); dbh = np.array([dij(Y, i, j) - dij(XR, i, j) for i, j in BONDS if heavy((i, j))])
    da = np.array([aijk(Y, *t_) - aijk(XR, *t_) for t_ in ANGLES]); dah = np.array([aijk(Y, *t_) - aijk(XR, *t_) for t_ in ANGLES if heavy(t_)])
    add(f'{p}.bond_rmsd_all', rms(db), 'Å', c, T + c + '/opt_summary.json:bond_rmsd_A', f'34 bonds incl. 8 C-H/N-H; recomputed here {rms(db):.4f} = file {op["bond_rmsd_A"]:.4f}; dominated by short X-ray X-H')
    add(f'{p}.bond_max_dev_all', float(abs(db).max()), 'Å', c, T + c + '/opt_summary.json:bond_maxdev_A', 'C5-H5 (X-ray 0.808 Å)')
    add(f'{p}.bond_rmsd_heavy', rms(dbh), 'Å', c, src_opt, 'derived: 26 bonds between non-H atoms, DFT - X-ray (cluster xyz)')
    add(f'{p}.bond_mae_heavy', float(abs(dbh).mean()), 'Å', c, src_opt, 'derived: mean |DFT - X-ray| over 26 heavy-atom bonds')
    order = np.argsort(-abs(dbh)); hb_list = [b for b in BONDS if heavy(b)]
    add(f'{p}.bond_largest_dev_heavy', '; '.join(f'{LAB[hb_list[k_][0]]}-{LAB[hb_list[k_][1]]} {dbh[k_]:+.3f}' for k_ in order[:5]), 'Å', c, src_opt, 'derived: five largest |DFT - X-ray| heavy-atom bonds')
    add(f'{p}.bond_mean_signed_dev_heavy', float(dbh.mean()), 'Å', c, src_opt, 'derived: mean(DFT - X-ray); positive = DFT longer')
    add(f'{p}.angle_rmsd_all', rms(da), 'deg', c, src_opt, 'derived: 52 bond angles from X-ray connectivity')
    add(f'{p}.angle_rmsd_heavy', rms(dah), 'deg', c, src_opt, 'derived: 36 angles with no H atom')
    ang_h = [t_ for t_ in ANGLES if heavy(t_)]; oa = np.argsort(-abs(dah))
    add(f'{p}.angle_largest_dev_heavy', '; '.join(f'{LAB[ang_h[k_][0]]}-{LAB[ang_h[k_][1]]}-{LAB[ang_h[k_][2]]} {dah[k_]:+.2f}' for k_ in oa[:5]), 'deg', c, src_opt, 'derived')
    hv = [i for i in range(len(SYM)) if SYM[i] != 'H']
    add(f'{p}.kabsch_rmsd_heavy', kabsch(Y[hv], XR[hv]), 'Å', c, src_opt, 'derived: best-fit superposition RMSD of all 25 heavy atoms vs X-ray cluster (measures the relative ion orientation)')
    add(f'{p}.kabsch_rmsd_picrate_heavy', kabsch(Y[[i for i in PIC if SYM[i] != "H"]], XR[[i for i in PIC if SYM[i] != "H"]]), 'Å', c, src_opt, 'derived: picrate heavy atoms only (nitro rotations)')
    add(f'{p}.kabsch_rmsd_cation_heavy', kabsch(Y[[i for i in CAT if SYM[i] != "H"]], XR[[i for i in CAT if SYM[i] != "H"]]), 'Å', c, src_opt, 'derived: cation heavy atoms only')
    # ---- frequencies
    fr = fq['frequencies_cm-1']
    add(f'{p}.n_imag', fq['n_imaginary'], '', c, T + c + '/freq.json:n_imaginary', 'from the project projection (dft/vibrations.py)')
    add(f'{p}.n_imag_pyscf', fq['n_imaginary_pyscf'], '', c, T + c + '/freq.json:n_imaginary_pyscf', 'PySCF harmonic_analysis count (misses 3.2i in PCM because it stores it as a complex number with zero real part)')
    add(f'{p}.hessian_method', fq['hessian_method'], '', c, T + c + '/freq.json:hessian_method')
    add(f'{p}.hessian_wall', fq['hessian_wall_min'], 'min', c, T + c + '/freq.json:hessian_wall_min')
    add(f'{p}.n_modes', len(fr), '', c, T + c + '/freq.json', '3N - 6 = 93')
    add(f'{p}.freq_lowest10', ', '.join(f'{v:.2f}' for v in fr[:10]), 'cm-1', c, T + c + '/freq.json:frequencies_cm-1[0:10]', 'negative = imaginary')
    add(f'{p}.freq_lowest', fr[0], 'cm-1', c, T + c + '/freq.json:frequencies_cm-1[0]')
    add(f'{p}.freq_highest', fr[-1], 'cm-1', c, T + c + '/freq.json:frequencies_cm-1[-1]', 'N6-H6N stretch (free N-H)')
    add(f'{p}.n_modes_below_100', sum(1 for v in fr if v < 100), '', c, T + c + '/freq.json', 'derived; counts every listed value < 100 cm-1, including an imaginary mode (PCM: 10 = 1 imaginary + 9 real). The thermochemistry uses only the real modes (ZPE = 0.5 x sum over 92 real modes for PCM).')
    add(f'{p}.zpe', fq['ZPE_Ha'], 'Ha', c, T + c + '/freq.json:ZPE_Ha', f'= {fq["ZPE_Ha"] * HA_KCAL:.2f} kcal mol-1 = {fq["ZPE_Ha"] * HA_EV:.3f} eV; equals 0.5 x sum of the projected real frequencies to 2e-11 Ha')
    add(f'{p}.zpe_kcal', fq['ZPE_Ha'] * HA_KCAL, 'kcal mol-1', c, T + c + '/freq.json:ZPE_Ha', 'derived: x 627.509474')
    add(f'{p}.h298', fq['H298_Ha'], 'Ha', c, T + c + '/freq.json:H298_Ha', '298.15 K, 1 atm, RRHO, unscaled')
    add(f'{p}.g298', fq['G298_Ha'], 'Ha', c, T + c + '/freq.json:G298_Ha', '298.15 K, 1 atm, RRHO; entropy of the many <100 cm-1 modes treated harmonically')
    add(f'{p}.thermal_corr_h', fq['H298_Ha'] - fq['E_Ha'], 'Ha', c, T + c + '/freq.json', 'derived: H298 - E')
    add(f'{p}.thermal_corr_g', fq['G298_Ha'] - fq['E_Ha'], 'Ha', c, T + c + '/freq.json', 'derived: G298 - E')
    add(f'{p}.entropy', (fq['H298_Ha'] - fq['G298_Ha']) * HA_KCAL * 1000 / 298.15, 'cal mol-1 K-1', c, T + c + '/freq.json', 'derived: (H - G)/T')
    s_, x_, m_, f_, V_, qTR = modes_for(c); MODES[c] = (s_, x_, m_, f_, V_, qTR)
    add(f'{p}.freq_recomputed_check', float(abs(f_ - np.array(fr)).max()), 'cm-1', c, RAW + c + '/hessian.npy',
        'derived: max |recomputed - freq.json| over 93 modes from the stored Hessian with the masses PySCF used (H 1.008, C 12.011, N 14.007, O 15.999); numerical noise only')
    for H, D in (('H4N', 'N4'), ('H6N', 'N6')):
        o = nh_modes(x_, m_, f_, V_, H, D)
        add(f'{p}.nu_{D}_{H}_stretch', fr[o['str'][2]], 'cm-1', c, T + c + f"/freq.json:frequencies_cm-1[{o['str'][2]}]",
            f'value from freq.json; assignment derived from the eigenvectors of {RAW}{c}/hessian.npy: H weight in the mass-weighted eigenvector {o["str"][1]:.2f}, >50% of the H motion along the {D}-{H} bond; unscaled harmonic')
        add(f'{p}.nu_{D}_{H}_inplane_bend', fr[o['ip'][2]], 'cm-1', c, T + c + f"/freq.json:frequencies_cm-1[{o['ip'][2]}]",
            f'assignment derived from hessian.npy: mode with the largest in-plane H motion perpendicular to the bond (H weight {o["ip"][1]:.2f}); coupled with ring modes')
        add(f'{p}.nu_{D}_{H}_oop_wag', ', '.join(f'{fr[k2]:.1f}' for _, _, k2 in o['oop']), 'cm-1', c, T + c + '/freq.json',
            'assignment derived from hessian.npy: modes with H weight > 0.15 and >50% out-of-plane H motion')
    Q = rigid_basis(x_, m_, [PIC, CAT])
    for k_ in range(3):
        add(f'{p}.mode{k_ + 1}_rigid_fraction', float(np.linalg.norm(Q.T @ V_[:, k_]) ** 2), '', c, RAW + c + '/hessian.npy',
            f'derived: fraction of normal mode {k_ + 1} ({f_[k_]:.2f} cm-1) lying in the 12-dim space of rigid translations/rotations of the two ions')
    # ---- FMO
    nocc = fm['n_occ']; E = fm['mo_energies_eV']
    add(f'{p}.homo_ev', fm['HOMO_eV'], 'eV', c, T + c + '/fmo.json:HOMO_eV')
    add(f'{p}.lumo_ev', fm['LUMO_eV'], 'eV', c, T + c + '/fmo.json:LUMO_eV')
    add(f'{p}.gap_ev', fm['gap_eV'], 'eV', c, T + c + '/fmo.json:gap_eV', 'E_LUMO - E_HOMO')
    add(f'{p}.homo_m1_ev', fm['HOMO-1_eV'], 'eV', c, T + c + '/fmo.json:HOMO-1_eV')
    add(f'{p}.lumo_p1_ev', fm['LUMO+1_eV'], 'eV', c, T + c + '/fmo.json:LUMO+1_eV')
    add(f'{p}.ip_ev', fm['IP_eV'], 'eV', c, T + c + '/fmo.json:IP_eV', 'IP = -E_HOMO (Koopmans)')
    add(f'{p}.ea_ev', fm['EA_eV'], 'eV', c, T + c + '/fmo.json:EA_eV', 'EA = -E_LUMO')
    add(f'{p}.chi_ev', fm['electronegativity_chi_eV'], 'eV', c, T + c + '/fmo.json:electronegativity_chi_eV', 'chi = (IP + EA)/2')
    add(f'{p}.mu_ev', fm['chemical_potential_mu_eV'], 'eV', c, T + c + '/fmo.json:chemical_potential_mu_eV', 'mu = -(IP + EA)/2')
    add(f'{p}.eta_ev', fm['hardness_eta_eV'], 'eV', c, T + c + '/fmo.json:hardness_eta_eV', 'eta = (IP - EA)/2')
    add(f'{p}.softness_ev_inv', fm['softness_S_eV-1'], 'eV-1', c, T + c + '/fmo.json:softness_S_eV-1', 'S = 1/(2 eta)  [convention; 1/eta would give twice this]')
    add(f'{p}.omega_ev', fm['electrophilicity_omega_eV'], 'eV', c, T + c + '/fmo.json:electrophilicity_omega_eV', 'omega = mu^2/(2 eta)')
    add(f'{p}.dipole_debye', fm['dipole_Debye'], 'D', c, T + c + '/fmo.json:dipole_Debye', 'total dipole of the neutral ion pair (origin-independent); PCM value from the polarised density')
    add(f'{p}.dipole_vector_debye', ', '.join(f'{v:.3f}' for v in fm['dipole_vector_Debye']), 'D', c, T + c + '/fmo.json:dipole_vector_Debye', 'x, y, z in the opt.xyz frame')
    add(f'{p}.n_occ', nocc, '', c, T + c + '/fmo.json:n_occ')
    for nm in ('HOMO-1', 'HOMO', 'LUMO', 'LUMO+1'):
        fc = fm['fragment_composition'][nm]
        add(f'{p}.{nm.lower().replace("-", "_m").replace("+", "_p")}_picrate_pct', fc['picrate anion'], '%', c, T + c + f'/fmo.json:fragment_composition.{nm}',
            f"Mulliken; cation {fc['benzotriazolium cation']} %")
    atoms_, shells, mos = parse_molden(RAW + c + '/final.molden')
    S, aof = overlap(atoms_, shells)
    nreal = sum(1 for mo in mos if np.abs(mo['c']).sum() > 0)
    Cm = np.array([mo['c'] for mo in mos]).T[:, :nreal]
    orth = float(abs(Cm.T @ S @ Cm - np.eye(nreal)).max())
    add(f'{p}.molden_orthonormality_check', orth, '', c, RAW + c + '/final.molden', 'derived: max |C^T S C - I| with the overlap matrix rebuilt from the molden basis; validates the composition analysis below')
    isp = np.array([a in PIC for a in aof]); SC = S @ Cm
    share = np.array([float((Cm[isp, k_] * SC[isp, k_]).sum()) for k_ in range(nreal)])
    occ = np.array([mo.get('occ', 0) for mo in mos])[:nreal]
    with open(OUT / f'derived/mo_fragment_composition_{c}.csv', 'w') as fcsv:
        w = csv.writer(fcsv); w.writerow(['index0', 'label', 'energy_eV', 'occupation', 'picrate_pct', 'cation_pct'])
        for k_ in range(nreal):
            w.writerow([k_, orbname(k_, nocc), f'{mos[k_]["E"] * HA_EV:.5f}', f'{occ[k_]:.1f}', f'{100 * share[k_]:.2f}', f'{100 * (1 - share[k_]):.2f}'])
    for nm in ('HOMO-1', 'HOMO', 'LUMO', 'LUMO+1'):
        k_ = {'HOMO-1': nocc - 2, 'HOMO': nocc - 1, 'LUMO': nocc, 'LUMO+1': nocc + 1}[nm]
        assert abs(100 * share[k_] - fm['fragment_composition'][nm]['picrate anion']) < 0.06, (c, nm, share[k_])
    rows = []
    for k_ in range(nocc - 8, nocc + 8):
        rows.append(f'{orbname(k_, nocc)} {E[k_]:.3f} eV ({100 * share[k_]:.1f}% picrate)')
    add(f'{p}.orbitals_near_gap', '; '.join(rows), '', c, T + c + '/fmo.json:mo_energies_eV + ' + RAW + c + '/final.molden',
        'energies from fmo.json; fragment shares derived here (Mulliken with rebuilt overlap; reproduces the fmo.json HOMO-1..LUMO+1 values to 0.05%)')
    cat_occ = [k_ for k_ in range(nocc) if share[k_] < 0.5][-1]; cat_virt = [k_ for k_ in range(nocc, nreal) if share[k_] < 0.5][0]
    pic_occ = [k_ for k_ in range(nocc) if share[k_] > 0.5][-1]; pic_virt = [k_ for k_ in range(nocc, nreal) if share[k_] > 0.5][0]
    for tag_, k_ in (('cation_highest_occ', cat_occ), ('cation_lowest_virt', cat_virt), ('picrate_highest_occ', pic_occ), ('picrate_lowest_virt', pic_virt)):
        add(f'{p}.{tag_}_ev', E[k_], 'eV', c, RAW + c + '/final.molden; ' + T + c + '/fmo.json:mo_energies_eV',
            f'derived: {orbname(k_, nocc)}, {100 * share[k_]:.1f}% on picrate (Mulliken)')
    add(f'{p}.ct_gap_picrate_to_cation_ev', E[cat_virt] - E[pic_occ], 'eV', c, RAW + c + '/final.molden; fmo.json',
        'derived: E(lowest cation-localised virtual) - E(highest picrate-localised occupied); inter-ion charge-transfer orbital gap')
    add(f'{p}.local_gap_picrate_ev', E[pic_virt] - E[pic_occ], 'eV', c, RAW + c + '/final.molden; fmo.json', 'derived: picrate-localised LUMO - HOMO')
    add(f'{p}.local_gap_cation_ev', E[cat_virt] - E[cat_occ], 'eV', c, RAW + c + '/final.molden; fmo.json', 'derived: cation-localised lowest virtual - highest occupied')
    Z = {'O': 8, 'N': 7, 'C': 6, 'H': 1}
    Dm = (Cm * occ) @ Cm.T; pop = np.einsum('ij,ji->i', Dm, S)
    qpic = sum(Z[SYM[a]] - pop[aof == a].sum() for a in PIC)
    add(f'{p}.mulliken_charge_picrate', float(qpic), 'e', c, RAW + c + '/final.molden', 'derived: Mulliken fragment charge (cation = -this). Diffuse basis: indicative only; atomic Mulliken charges not reported.')
    # ---- MEP
    add(f'{p}.mep_min_au', ms['global_min_au'], 'a.u.', c, T + c + '/mep_surface.json:global_min_au')
    add(f'{p}.mep_max_au', ms['global_max_au'], 'a.u.', c, T + c + '/mep_surface.json:global_max_au')
    ext = list(csv.DictReader(open(T + c + '/mep_extrema.csv')))
    gmin = [e for e in ext if e['global'] == 'yes' and e['type'] == 'minimum'][0]; gmax = [e for e in ext if e['global'] == 'yes' and e['type'] == 'maximum'][0]
    # unrounded kcal mol-1 values from Multiwfn's own surfanalysis.txt (revision round 2, figures_r2 #2 / QA r4-m1):
    # mep_extrema.csv stores them rounded to 3 decimals, and rounding those again to 2 put six table cells 0.01 off
    sa_ = {'minimum': {}, 'maximum': {}}; kind_ = None
    for ln_ in open(RAW + c + '/surfanalysis.txt'):
        if 'surface minima' in ln_: kind_ = 'minimum'; continue
        if 'surface maxima' in ln_: kind_ = 'maximum'; continue
        t_ = ln_.replace('*', ' ').split()
        if kind_ and len(t_) == 7 and t_[0].isdigit():
            sa_[kind_][int(t_[0])] = float(t_[3])
    add(f'{p}.mep_min_kcal', sa_['minimum'][int(gmin['index'])], 'kcal mol-1', c, RAW + c + '/surfanalysis.txt (global minimum; mep_extrema.csv)', f"nearest atoms {gmin['nearest_atoms']}; unrounded (the CSV gives {gmin['ESP_kcal_mol']})")
    add(f'{p}.mep_max_kcal', sa_['maximum'][int(gmax['index'])], 'kcal mol-1', c, RAW + c + '/surfanalysis.txt (global maximum; mep_extrema.csv)', f"nearest atoms {gmax['nearest_atoms'].replace('#4_666', '')}; unrounded (the CSV gives {gmax['ESP_kcal_mol']})")
    add(f'{p}.mep_min_site', gmin['nearest_atoms'].split(' (')[0], '', c, T + c + '/mep_extrema.csv', 'nitro O of the para (N2) group, not the phenolate O1')
    add(f'{p}.mep_max_site', gmax['nearest_atoms'].split(' (')[0].replace('#4_666', ''), '', c, T + c + '/mep_extrema.csv', 'the free N6-H6N of the cation (its crystal H-bond partner is absent from the cluster)')
    add(f'{p}.mep_n_minima', ms['n_minima'], '', c, T + c + '/mep_surface.json:n_minima')
    add(f'{p}.mep_n_maxima', ms['n_maxima'], '', c, T + c + '/mep_surface.json:n_maxima')
    for e in sorted(ext, key=lambda e: float(e['ESP_kcal_mol'])):
        near = e['nearest_atoms'].replace('#4_666', '')
        a1 = near.split(' (')[0]; frag = 'picrate' if L[a1] in PIC else 'cation'
        kc_ = sa_[e['type']][int(e['index'])]
        assert abs(kc_ - float(e['ESP_kcal_mol'])) < 6e-4, (c, e, kc_)
        add(f"{p}.mep_{e['type'][:3]}{int(e['index']):02d}", kc_, 'kcal mol-1', c, RAW + c + f"/surfanalysis.txt:{e['type']} {e['index']} (mep_extrema.csv)",
            f"{e['ESP_au']} a.u.; at ({e['x_A']}, {e['y_A']}, {e['z_A']}) Å; nearest {near}; {frag}{'; GLOBAL' if e['global'] == 'yes' else ''}; "
            f"kcal mol-1 value unrounded from Multiwfn surfanalysis.txt (mep_extrema.csv rounds it to {float(e['ESP_kcal_mol']):.3f})")
    for key, unit, lab_ in [('area_A2', 'Å^2', 'surface_area'), ('positive_area_A2', 'Å^2', 'positive_area'), ('negative_area_A2', 'Å^2', 'negative_area'),
                            ('volume_A3', 'Å^3', 'volume'), ('overall_avg_kcal', 'kcal mol-1', 'avg_all'), ('positive_avg_kcal', 'kcal mol-1', 'avg_pos'),
                            ('negative_avg_kcal', 'kcal mol-1', 'avg_neg'), ('sigma2_tot_kcal2', '(kcal mol-1)^2', 'sigma2_tot'), ('nu', '', 'nu'),
                            ('Pi_kcal', 'kcal mol-1', 'Pi'), ('MPI_eV', 'eV', 'mpi'), ('polar_area_pct', '%', 'polar_area_pct')]:
        add(f'{p}.mep_{lab_}', ms[key], unit, c, T + c + f'/mep_surface.json:{key}', 'Multiwfn surface statistics on the 0.001 a.u. isosurface')
    mo_txt = re.sub(r'\x1b\[[0-9;]*m', '', open(RAW + c + '/logs/mwfn_mep.out', errors='replace').read())
    def grab(pat): return float(re.search(pat, mo_txt).group(1))
    add(f'{p}.mep_sigma2_pos', grab(r'Positive variance:.*?\(\s*([\d.]+) \(kcal'), '(kcal mol-1)^2', c, RAW + c + '/logs/mwfn_mep.out:Positive variance')
    add(f'{p}.mep_sigma2_neg', grab(r'Negative variance:.*?\(\s*([\d.]+) \(kcal'), '(kcal mol-1)^2', c, RAW + c + '/logs/mwfn_mep.out:Negative variance')
    add(f'{p}.mep_sigma2_nu', grab(r'Product of sigma\^2_tot and nu:.*?\(\s*([\d.]+) \(kcal'), '(kcal mol-1)^2', c, RAW + c + '/logs/mwfn_mep.out:Product of sigma^2_tot and nu')
    add(f'{p}.mep_mpi_kcal', grab(r'Molecular polarity index \(MPI\):.*?\(\s*([\d.]+) kcal'), 'kcal mol-1', c, RAW + c + '/logs/mwfn_mep.out:MPI')
    add(f'{p}.mep_nonpolar_area', grab(r'Nonpolar surface area.*?:\s*([\d.]+) Angstrom'), 'Å^2', c, RAW + c + '/logs/mwfn_mep.out', '|ESP| <= 10 kcal mol-1')
    add(f'{p}.mep_polar_area', grab(r'Polar surface area.*?:\s*([\d.]+) Angstrom'), 'Å^2', c, RAW + c + '/logs/mwfn_mep.out', '|ESP| > 10 kcal mol-1')
    add(f'{p}.mep_skew_all', grab(r'Overall skewness:\s*([-\d.]+)'), '', c, RAW + c + '/logs/mwfn_mep.out')
    add(f'{p}.mep_skew_pos', grab(r'Positive skewness:\s*([-\d.]+)'), '', c, RAW + c + '/logs/mwfn_mep.out')
    add(f'{p}.mep_skew_neg', grab(r'Negative skewness:\s*([-\d.]+)'), '', c, RAW + c + '/logs/mwfn_mep.out')
    add(f'{p}.mep_sphericity', grab(r'Sphericity:\s*([\d.]+)'), '', c, RAW + c + '/logs/mwfn_mep.out')
    add(f'{p}.mep_density_mv', grab(r'\(M/V\):\s*([\d.]+)'), 'g cm-3', c, RAW + c + '/logs/mwfn_mep.out', 'Multiwfn mass/volume estimate inside the 0.001 a.u. surface (cf. X-ray 1.645)')
    add(f'{p}.mep_n_vertices', int(re.findall(r'After elimination, V=\s*(\d+)', mo_txt)[-1]), '', c, RAW + c + '/logs/mwfn_mep.out', 'surface vertices after merging; = rows of vtx.txt')
    add(f'{p}.mep_wall', grab(r'Calculation of mapped function took up wall clock time\s+(\d+)'), 's', c, RAW + c + '/logs/mwfn_mep.out')
    # ---- NCI
    sp = list(csv.DictReader(open(T + c + '/nci_spikes.csv')))
    D_ = np.loadtxt(RAW + c + '/rdg_scatter.txt')
    for k_, row in enumerate(sp):
        cen = float(row['sign_l2_rho_au']); sel = (D_[:, 3] >= cen - 0.0005) & (D_[:, 3] < cen + 0.0005); rmin = D_[sel, 4].min()
        pts = D_[sel & (D_[:, 4] < rmin + 0.05)]; tags = {}
        for pt in pts:
            dd = np.linalg.norm(Y - pt[:3], axis=1); a, b = np.argsort(dd)[:2]
            key = f'{LAB[a]}/{LAB[b]}' + (' (inter-ion)' if (a in PIC) != (b in PIC) else '')
            tags[key] = tags.get(key, 0) + 1
        where = ', '.join(sorted(tags, key=lambda t_: -tags[t_]))
        add(f'{p}.nci_spike{k_ + 1}', cen, 'a.u.', c, T + c + f'/nci_spikes.csv:row {k_ + 1}',
            f"sign(l2)rho bin centre; min RDG {row['min_RDG']}; {row['n_points']} points in bin; pipeline label '{row['type']}'; lowest-RDG points lie between {where} (derived from rdg_scatter.txt)")
    del D_
    org, spc, n_, s2 = read_cube(RAW + c + '/sl2r.cube'); _, _, _, rg = read_cube(RAW + c + '/rdg.cube')
    from scipy import ndimage
    lab_, nl = ndimage.label(rg < 0.5); idx = np.indices(n_).reshape(3, -1).T; regs = []
    for k_ in range(1, nl + 1):
        sel = np.where(lab_.ravel() == k_)[0]
        if len(sel) < 20: continue
        sv = s2.ravel()[sel]; rv = rg.ravel()[sel]; pts = (org + idx[sel] @ spc) * BOHR; j = np.argmin(rv)
        dd = np.linalg.norm(Y - pts[j], axis=1); a, b = np.argsort(dd)[:2]
        regs.append((float(sv[j]), float(rv[j]), float(sv.min()), float(sv.max()), len(sel), LAB[a], LAB[b], (a in PIC) != (b in PIC)))
    regs.sort()
    for k_, rg_ in enumerate(regs):
        add(f'{p}.nci_region{k_ + 1:02d}', rg_[0], 'a.u.', c, RAW + c + '/sl2r.cube + rdg.cube',
            f'derived: connected RDG<0.5 region ({rg_[4]} grid points); value = sign(l2)rho at its lowest-RDG point (RDG {rg_[1]:.3f}); range {rg_[2]:+.4f}..{rg_[3]:+.4f}; '
            f'lowest-RDG point nearest {rg_[5]} and {rg_[6]}' + (' (inter-ion)' if rg_[7] else ' (within one ion)'))
    # ---- ELF / LOL planes
    for fn in ('elf', 'lol'):
        for pl in ('picrate_ring', 'benzotriazolium_ring'):
            z = np.load(RAW + c + f'/{fn}_{pl}.npz'); F = z['F']
            add(f'{p}.{fn}_{pl}_range', f'{F.min():.4f}..{F.max():.4f}', '', c, RAW + c + f'/{fn}_{pl}.npz:F', f'300 x 300 map; also logs/mwfn.log')
    from scipy.interpolate import RegularGridInterpolator
    for fn in ('elf', 'lol'):
        z = np.load(RAW + c + f'/{fn}_benzotriazolium_ring.npz'); X_, Y_, F = z['X'], z['Y'], z['F']
        if np.allclose(X_[:, 0], X_[0, 0]): ip_ = RegularGridInterpolator((Y_[:, 0], X_[0, :]), F); fv = lambda uv: ip_(uv[:, ::-1])
        else: ip_ = RegularGridInterpolator((X_[:, 0], Y_[0, :]), F); fv = lambda uv: ip_(uv)
        axy = z['atoms_xy']
        tri = [L[a] for a in ('N4', 'N5', 'N6', 'C7', 'C12')]; ben = [L[a] for a in ('C7', 'C8', 'C9', 'C10', 'C11', 'C12')]
        add(f'{p}.{fn}_triazole_centroid', float(fv(axy[tri].mean(0)[None])[0]), '', c, RAW + c + f'/{fn}_benzotriazolium_ring.npz', 'derived: interpolated at the ring centroid')
        add(f'{p}.{fn}_benzo_centroid', float(fv(axy[ben].mean(0)[None])[0]), '', c, RAW + c + f'/{fn}_benzotriazolium_ring.npz', 'derived')
        z = np.load(RAW + c + f'/{fn}_picrate_ring.npz'); X_, Y_, F = z['X'], z['Y'], z['F']
        if np.allclose(X_[:, 0], X_[0, 0]): ip_ = RegularGridInterpolator((Y_[:, 0], X_[0, :]), F); fv = lambda uv: ip_(uv[:, ::-1])
        else: ip_ = RegularGridInterpolator((X_[:, 0], Y_[0, :]), F); fv = lambda uv: ip_(uv)
        add(f'{p}.{fn}_picrate_centroid', float(fv(z['atoms_xy'][pring].mean(0)[None])[0]), '', c, RAW + c + f'/{fn}_picrate_ring.npz', 'derived')
    z = np.load(RAW + c + '/elf_benzotriazolium_ring.npz'); cen_, nrm_ = z['center'], z['normal']
    add(f'{p}.cation_plane_offset_O1', float(np.dot(Y[L['O1']] - cen_, nrm_)), 'Å', c, RAW + c + '/elf_benzotriazolium_ring.npz + opt.xyz',
        'derived: signed distance of O1 from the ELF/LOL cation plane (H4N lies within 0.01 Å of it), i.e. how well the cation-plane map shows the H4N...O1 region')
    add(f'{p}.cube_grid', '80 x 80 x 80 (HOMO, LUMO, density), 64 x 64 x 64 (ESP), Bohr units, same origin', '', c, RAW + c + '/*.cube headers')

# ---- extra PCM saddle-point diagnostics
c = 'btpc_pcm'; s_, x_, m_, f_, V_, qTR = MODES[c]
f5 = jload(T + c + '/freq_step0.005.json')
add('saddle.freq_imag_step0p01', freq[c]['frequencies_cm-1'][0], 'cm-1', c, T + c + '/freq.json:frequencies_cm-1[0]', 'authoritative Hessian: step 0.01 Bohr, 10 workers x 2 threads, 2026-09-18 16:05-21:51 UTC')
add('saddle.freq_imag_step0p005', f5['frequencies_cm-1'][0], 'cm-1', c, T + c + '/freq_step0.005.json:frequencies_cm-1[0]', 'first successful Hessian: step 0.005 Bohr, 12 workers x 2 threads, 2026-09-17 13:56-19:42 UTC')
dfr = np.array(freq[c]['frequencies_cm-1']) - np.array(f5['frequencies_cm-1'])
add('saddle.step_change_max', float(abs(dfr).max()), 'cm-1', c, T + c + '/freq.json vs freq_step0.005.json', f'derived: max |nu(0.01) - nu(0.005)| over 93 modes (mode at {f5["frequencies_cm-1"][int(np.argmax(abs(dfr)))]:.1f} cm-1); mean |diff| {abs(dfr).mean():.3f}')
add('saddle.step_change_zpe', (freq[c]['ZPE_Ha'] - f5['ZPE_Ha']) * HA_KCAL, 'kcal mol-1', c, 'freq.json vs freq_step0.005.json', 'derived')
add('saddle.step_change_g', (freq[c]['G298_Ha'] - f5['G298_Ha']) * HA_KCAL, 'kcal mol-1', c, 'freq.json vs freq_step0.005.json', 'derived')
add('saddle.first_attempt', 'failed: 2026-09-17 12:21-12:48 UTC, 19 workers x 2 threads, exit code 1', '', c, T + c + '/STAGES.md', 'no Hessian from this attempt')
add('saddle.imag_mode_rigid_fraction', float(np.linalg.norm(rigid_basis(x_, m_, [PIC, CAT]).T @ V_[:, 0]) ** 2), '', c, RAW + c + '/hessian.npy',
    'derived: 99.9% of the -3.2 cm-1 eigenvector is rigid-body motion of the two ions relative to each other (inter-ion libration), not an internal distortion')
w0 = (V_[:, 0].reshape(-1, 3) ** 2).sum(1)
add('saddle.imag_mode_picrate_weight', float(w0[PIC].sum()), '', c, RAW + c + '/hessian.npy', 'derived: mass-weighted participation on picrate (cation 0.66); largest atomic weights C10, N5, N6, C9 (cation)')
def read_traj(pth):
    Ls = open(pth).read().split('\n'); frames = []; i = 0
    while i < len(Ls) and Ls[i].strip():
        n = int(Ls[i]); frames.append((Ls[i + 1], np.array([[float(v) for v in l.split()[1:4]] for l in Ls[i + 2:i + 2 + n]]))); i += n + 2
    return frames
frames = read_traj(RAW + c + '/opt_traj.xyz'); xf = x_ * BOHR; mw = m_ / AMU2ME
def align(A, B, w):
    ca = (w[:, None] * A).sum(0) / w.sum(); cb = (w[:, None] * B).sum(0) / w.sum(); a = (A - ca) * np.sqrt(w)[:, None]; b = (B - cb) * np.sqrt(w)[:, None]
    U, S_, Vt = np.linalg.svd(a.T @ b); dd = np.sign(np.linalg.det(Vt.T @ U.T)); Rm = Vt.T @ np.diag([1, 1, dd]) @ U.T; return (Rm @ (A - ca).T).T + cb
for k0 in (0, 40, 60):
    A = align(frames[k0][1], xf, mw); dx = xf - A; q = (dx * np.sqrt(mw)[:, None]).ravel(); q = q - qTR @ (qTR.T @ q); q /= np.linalg.norm(q)
    ov = (V_.T @ q) ** 2; E0 = float(frames[k0][0].split('=')[1].split()[0]); E1 = float(frames[-1][0].split('=')[1].split()[0])
    add(f'saddle.drift_{k0}_to_70_overlap_imag', float(ov[0]), '', c, RAW + c + '/opt_traj.xyz + hessian.npy',
        f'derived: squared overlap of the mass-weighted displacement (frame {k0} -> 70, aligned, rigid rotation removed) with the -3.2 cm-1 mode; '
        f'overlap with the 9.0 cm-1 mode {ov[1]:.3f}; max atom displacement {np.linalg.norm(dx, axis=1).max():.2f} Å; energy change {(E1 - E0) * HA_KCAL:.3f} kcal mol-1')
pring_i = pring
for k0 in (0, 25, 30, 40, 45, 50, 60, 70):
    Yk = frames[k0][1]; c1, n1 = plane(Yk, pring); c2, n2 = plane(Yk, cat9)
    ia = math.degrees(math.acos(min(1, abs(float(np.dot(n1, n2)))))); Ek = float(frames[k0][0].split('=')[1].split()[0])
    add(f'saddle.traj_step{k0:02d}', f"E {Ek:.8f} Ha; N4-H4N {dij(Yk, L['N4'], L['H4N']):.3f}; H4N...O1 {dij(Yk, L['H4N'], L['O1']):.3f}; N4...O1 {dij(Yk, L['N4'], L['O1']):.3f} Å; "
        f"angle {aijk(Yk, L['N4'], L['H4N'], L['O1']):.1f} deg; interplanar {ia:.1f} deg; H11...O2 {dij(Yk, L['H11'], L['O2']):.3f} Å", '', c, RAW + c + '/opt_traj.xyz',
        'derived: geometry along the PCM optimisation (shows the H-bond is stable while the ions rotate)')
add('saddle.interpretation', 'The PCM stationary-point search ended on an almost flat inter-ion libration coordinate: over the last 30 steps the ions rotated by ~38 deg for 0.30 kcal/mol, the residual imaginary frequency (3.2i cm-1) belongs to that same rigid-body motion and survived a doubled FD step. The intra-ion geometry and the N4-H4N...O1 bond are converged; the relative orientation of the two rings is not determined to better than a few tens of degrees.',
    '', c, 'derived from the saddle.* facts', 'interpretation for the writers; every number it rests on is in saddle.*')

# =====================================================================================  5. gas vs PCM differences
def diff(key, unit='', scale=1.0, note=''):
    g, pv = facts[f'gas.{key}']['value'], facts[f'pcm.{key}']['value']
    add(f'diff.{key}', (pv - g) * scale, unit or facts[f'pcm.{key}']['unit'], None, f'pcm.{key}, gas.{key}',
        f'derived: PCM - gas' + (f' (x {scale})' if scale != 1 else '') + (f'; {note}' if note else ''))
diff('e_total', 'kcal mol-1', HA_KCAL, 'total-energy difference of the two optimised structures; contains the PCM electrostatic solvation term and the geometry change; NOT a solvation free energy (no non-electrostatic terms, different geometries)')
add('diff.e_total_ha', facts['pcm.e_total']['value'] - facts['gas.e_total']['value'], 'Ha', None, 'pcm.e_total, gas.e_total', 'derived: PCM - gas')
diff('g298', 'kcal mol-1', HA_KCAL, 'same caveats as diff.e_total')
diff('zpe', 'kcal mol-1', HA_KCAL)
for k_ in ['r_N4_H4N', 'hb_N4_H4N_O1_ha', 'hb_N4_H4N_O1_da', 'hb_N4_H4N_O1_angle', 'hb_N4_H4N_O7_ha', 'hb_C11_H11_O2_ha', 'r_O1_C1', 'r_N4_N5', 'r_N5_N6', 'r_N6_H6N',
           'r_C1_C2', 'r_C1_C6', 'a_C2_C1_C6', 'interplanar_picrate_cation', 'centroid_distance_picrate_cation', 'nitro_twist_N1', 'nitro_twist_N2', 'nitro_twist_N3',
           'nh_minus_ho_asym', 'bond_rmsd_all', 'bond_rmsd_heavy', 'angle_rmsd_heavy', 'kabsch_rmsd_heavy',
           'homo_ev', 'lumo_ev', 'gap_ev', 'homo_m1_ev', 'lumo_p1_ev', 'ip_ev', 'ea_ev', 'chi_ev', 'mu_ev', 'eta_ev', 'softness_ev_inv', 'omega_ev', 'dipole_debye',
           'cation_highest_occ_ev', 'cation_lowest_virt_ev', 'picrate_highest_occ_ev', 'picrate_lowest_virt_ev', 'ct_gap_picrate_to_cation_ev', 'local_gap_picrate_ev',
           'mulliken_charge_picrate',
           'mep_min_kcal', 'mep_max_kcal', 'mep_surface_area', 'mep_positive_area', 'mep_negative_area', 'mep_volume', 'mep_avg_all', 'mep_avg_pos', 'mep_avg_neg',
           'mep_sigma2_tot', 'mep_sigma2_pos', 'mep_sigma2_neg', 'mep_nu', 'mep_sigma2_nu', 'mep_Pi', 'mep_mpi', 'mep_polar_area_pct',
           'nu_N4_H4N_stretch', 'nu_N6_H6N_stretch', 'freq_lowest', 'freq_highest', 'h298', 'entropy']:
    diff(k_)
add('diff.nci_hbond_spike', facts['pcm.nci_spike1']['value'] - facts['gas.nci_spike1']['value'], 'a.u.', None, 'pcm.nci_spike1, gas.nci_spike1',
    'derived: PCM - gas of the most negative sign(l2)rho spike (both at H4N...O1); |value| ~ rho at the H...O contact')
add('diff.nh_redshift_pcm', facts['pcm.nu_N6_H6N_stretch']['value'] - facts['pcm.nu_N4_H4N_stretch']['value'], 'cm-1', 'btpc_pcm', 'pcm.nu_N6_H6N_stretch, pcm.nu_N4_H4N_stretch',
    'derived: nu(free N6-H6N) - nu(H-bonded N4-H4N) within the PCM model; an internal reference for the H-bond red shift')
add('diff.nh_redshift_gas', facts['gas.nu_N6_H6N_stretch']['value'] - facts['gas.nu_N4_H4N_stretch']['value'], 'cm-1', 'btpc_gas', 'gas.nu_N6_H6N_stretch, gas.nu_N4_H4N_stretch', 'derived: as above, gas phase')
add('diff.xray_vs_pcm_n_o', facts['pcm.hb_N4_H4N_O1_da']['value'] - facts['xray.hb_N4_H4N_O1_da']['value'], 'Å', None, 'pcm/xray hb_N4_H4N_O1_da', 'derived: DFT - X-ray N4...O1')
add('diff.xray_vs_gas_n_o', facts['gas.hb_N4_H4N_O1_da']['value'] - facts['xray.hb_N4_H4N_O1_da']['value'], 'Å', None, 'gas/xray hb_N4_H4N_O1_da', 'derived: DFT - X-ray N4...O1')

# =====================================================================================  6. known inconsistencies (machine-readable)
add('inconsistency.grid_level_text', 'RESULTS.md (both cases) says "integration grid level 3"; the runs used level 4', '', None, 'molecules/btpc/results/*/RESULTS.md vs verify_grid.log (not deposited)', 'authoritative: level 4')
add('inconsistency.pcm_hessian_step_text', 'btpc_pcm RESULTS.md methods say step 0.01 Bohr (true for freq.json) while the dft/frequencies.py comment says 0.005 is what "every completed case used"', '', 'btpc_pcm',
    'molecules/btpc/results/btpc_pcm/RESULTS.md; dft/frequencies.py comment', 'authoritative: freq.json hessian_method (0.01 Bohr); freq_step0.005.json is the 0.005 check')
add('inconsistency.pcm_convergence_text', 'handover README and opt pipeline text claim full convergence with Gaussian criteria; btpc_pcm met only the force criteria', '', 'btpc_pcm',
    'hand-over package README of Sep 2026 (not deposited) vs optimize.log', 'authoritative: optimize.log')
add('inconsistency.lebedev_text', 'dft/common.py defaults to Lebedev order 17; runs used 29', '', 'btpc_pcm', 'run_case.sh; verify_grid.log (not deposited)', 'authoritative: 29 (not directly verified for btpc_pcm)')
add('inconsistency.pcm_results_md_minimum', 'btpc_pcm RESULTS.md: "1 imaginary -> NOT a minimum" and "(nan min)" for the optimisation time; SUMMARY/charts say lowest -3.2', '', 'btpc_pcm',
    'molecules/btpc/results/btpc_pcm/RESULTS.md', 'consistent with freq.json; wall time recovered from STAGES.md')
add('inconsistency.fmo_label_threshold', 'STAGES.md prints PCM HOMO/LUMO composition as 99.9/0.1 for both while RESULTS.md table shows HOMO-1 0/100 - both consistent with fmo.json (HOMO-1 is 99.8% cation)', '', 'btpc_pcm', 'STAGES.md, RESULTS.md, fmo.json', 'no conflict; noted to avoid confusion')
add('inconsistency.nci_scatter_units', 'dft/mwfn.py comment says rdg_scatter.txt x,y,z are in Bohr; Multiwfn output states Angstrom and the values match Å', '', None,
    'dft/mwfn.py:nci_post vs logs/mwfn_nci.out ("Column 1/2/3: X/Y/Z in Angstrom")', 'authoritative: Angstrom')
add('inconsistency.chart_title', 'results/charts/btpc/SUMMARY.md title truncated to "1H-benzotriazol-1-ium 2"', '', None, 'results/charts/btpc/SUMMARY.md (pipeline chart summary, not deposited)', 'cosmetic')
add('inconsistency.gas_freq_policy', 'The project notes (not deposited) say frequencies only on the primary (PCM) model; btpc_gas nevertheless has a full FD Hessian (0 imaginary)', '', 'btpc_gas',
    'molecules/btpc/results/btpc_gas/freq.json; STAGES.md', 'authoritative: freq.json exists and is valid')

# =====================================================================================  7. writer's additions (manuscript, 23 Sep 2026)
# Values the manuscript and captions print that had no key of their own (OUTLINE.md section 9).  Simple
# arithmetic on existing keys is stored with its "expr" (numbers_audit.py re-evaluates it); everything else is
# recomputed here from the stored result files, with the formula in the note.
def addx(key, expr, unit='', case=None, note=''):
    names = sorted(set(re.findall(r'[A-Za-z_][\w]*\.[\w.]+', expr)), key=lambda n: (-len(n), n))   # longest first; ties by name (deterministic)
    e = expr
    for i, n in enumerate(names):
        e = re.sub(r'(?<![\w.])' + re.escape(n) + r'(?![\w.])', f'__v{i}', e)
    val = eval(e, {'__builtins__': {}, 'abs': abs}, {f'__v{i}': facts[n]['value'] for i, n in enumerate(names)})
    add(key, val, unit, case, ', '.join(names), 'derived: ' + expr + (f'; {note}' if note else ''))
    facts[key]['expr'] = expr

addx('xray.nh_minus_ho_asym', 'xray.hb_N4_H4N_O1_ha - xray.hb_N4_H4N_O1_dh', 'Å', None,
     'X-ray analogue of *.nh_minus_ho_asym; the X-ray N-H is short, so the value is an upper bound')
for p_ in ('xray', 'pcm', 'gas'):
    addx(f'{p_}.nn_bond_asym', f'{p_}.r_N5_N6 - {p_}.r_N4_N5', 'Å', None if p_ == 'xray' else f'btpc_{p_}',
         'asymmetry of the two triazole N-N bonds of the cation (N5-N6 minus N4-N5)')
addx('gas.lumo_cation_pct', '100 - gas.lumo_picrate_pct', '%', 'btpc_gas', 'Mulliken share of the gas-phase LUMO on the cation')
addx('pcm.lumo_cation_pct', '100 - pcm.lumo_picrate_pct', '%', 'btpc_pcm', 'Mulliken share of the IEF-PCM LUMO on the cation')
addx('pcm.homo_m1_cation_pct', '100 - pcm.homo_m1_picrate_pct', '%', 'btpc_pcm', 'IEF-PCM HOMO-1 is the cation pi orbital')
addx('gas.homo_cation_pct', '100 - gas.homo_picrate_pct', '%', 'btpc_gas', '')
addx('saddle.imag_mode_rigid_pct', '100 * saddle.imag_mode_rigid_fraction', '%', 'btpc_pcm', 'rigid-body share of the 3.2i cm-1 mode')
addx('gas.mode1_rigid_pct', '100 * gas.mode1_rigid_fraction', '%', 'btpc_gas', 'rigid-body share of the lowest gas-phase mode (12.2 cm-1)')
addx('saddle.imag_mode_cation_pct', '100 * (1 - saddle.imag_mode_picrate_weight)', '%', 'btpc_pcm', 'mass-weighted amplitude of the 3.2i cm-1 mode on the cation')
addx('saddle.drift_0_to_70_pct', '100 * saddle.drift_0_to_70_overlap_imag', '%', 'btpc_pcm', 'share of the step 0 -> 70 displacement along the 3.2i cm-1 mode')
addx('saddle.drift_40_to_70_pct', '100 * saddle.drift_40_to_70_overlap_imag', '%', 'btpc_pcm', 'share of the step 40 -> 70 displacement along the 3.2i cm-1 mode')
# drift of the last ten steps along the 9.0 cm-1 libration (same recipe as saddle.drift_*_overlap_imag)
A60 = align(frames[60][1], xf, mw); dx60 = xf - A60; q60 = (dx60 * np.sqrt(mw)[:, None]).ravel(); q60 = q60 - qTR @ (qTR.T @ q60); q60 /= np.linalg.norm(q60)
add('saddle.drift_60_to_70_pct_9cm', 100 * float(((V_.T @ q60) ** 2)[1]), '%', c, RAW + c + '/opt_traj.xyz + hessian.npy',
    'derived: 100 x squared overlap of the mass-weighted displacement frame 60 -> 70 (aligned, rigid rotation removed) with the 9.0 cm-1 mode (second lowest)')
# largest excursion of the hydrogen-bond distances from their step-0 values over the 71 frames
exc = {}
for nm_, (a_, b_) in {'N4-H4N': ('N4', 'H4N'), 'H4N...O1': ('H4N', 'O1'), 'N4...O1': ('N4', 'O1')}.items():
    d_ = np.array([dij(Yk, L[a_], L[b_]) for _, Yk in frames]); exc[nm_] = (float(abs(d_ - d_[0]).max()), float(d_.max() - d_.min()))
kmax = max(exc, key=lambda k_: exc[k_][0])
add('saddle.hbond_max_excursion', exc[kmax][0], 'Å', c, RAW + c + '/opt_traj.xyz',
    f'derived: largest |d(step k) - d(step 0)| of N4-H4N, H4N...O1, N4...O1 over the 71 frames ({kmax}); ' +
    '; '.join(f'{k_} excursion {v_[0]:.3f}, full range {v_[1]:.3f}' for k_, v_ in exc.items()))
# ring angle along the trajectory, same plane definition as saddle.traj_stepNN
ia_ = []
for _, Yk in frames:
    c1, n1 = plane(Yk, pring); c2, n2 = plane(Yk, cat9); ia_.append(math.degrees(math.acos(min(1, abs(float(np.dot(n1, n2)))))))
add('saddle.ring_rotation_0_to_70', ia_[-1] - ia_[0], 'deg', c, RAW + c + '/opt_traj.xyz',
    f'derived: inter-ring angle at step 70 minus step 0 ({ia_[0]:.1f} -> {ia_[-1]:.1f} deg)')
add('saddle.ring_rotation_40_to_70', ia_[-1] - ia_[40], 'deg', c, RAW + c + '/opt_traj.xyz',
    f'derived: inter-ring angle at step 70 minus step 40 ({ia_[40]:.1f} -> {ia_[-1]:.1f} deg)')
E_ = [float(fr[0].split('=')[1].split()[0]) for fr in frames]
add('saddle.energy_drop_40_to_70', (E_[-1] - E_[40]) * HA_KCAL, 'kcal mol-1', c, RAW + c + '/opt_traj.xyz', 'derived: E(step 70) - E(step 40)')
# step test: mean change and the mode with the largest change
kbig = int(np.argmax(abs(dfr)))
add('saddle.step_change_mean', float(abs(dfr).mean()), 'cm-1', c, T + c + '/freq.json vs freq_step0.005.json', 'derived: mean |nu(0.01) - nu(0.005)| over the 93 modes')
add('saddle.step_change_max_mode', freq[c]['frequencies_cm-1'][kbig], 'cm-1', c, T + c + '/freq.json:frequencies_cm-1',
    f'derived: wavenumber in the authoritative 0.01 Bohr Hessian of the mode with the largest change (mode {kbig + 1}); '
    f'the same mode lies at {f5["frequencies_cm-1"][kbig]:.1f} cm-1 in the 0.005 Bohr Hessian; signed change {dfr[kbig]:.3f} cm-1')
add('saddle.step_change_max_signed', float(dfr[kbig]), 'cm-1', c, T + c + '/freq.json vs freq_step0.005.json', 'derived: nu(0.01) - nu(0.005) of that mode')
# C-H stretching wavenumbers from the Hessian eigenvectors (all modes between 3000 and 3400 cm-1 are C-H stretches)
CH = {'C3': 'H3', 'C5': 'H5', 'C8': 'H8', 'C9': 'H9', 'C10': 'H10', 'C11': 'H11'}
chv = {}
for c_ in CASES:
    p = PFX[c_]; s_c, x_c, m_c, f_c, V_c, _ = MODES[c_]
    sel = [k_ for k_ in range(len(f_c)) if 3000 < f_c[k_] < 3400]
    assert len(sel) == 6, (c_, len(sel))
    own = {}
    for k_ in sel:
        v = V_c[:, k_].reshape(-1, 3); w = {h: float((v[L[h]] ** 2).sum()) for h in CH.values()}
        own[k_] = max(w, key=w.get)
    pic_modes = {h: [k_ for k_ in sel if own[k_] == h] for h in ('H3', 'H5')}
    assert all(len(v_) == 1 for v_ in pic_modes.values()), pic_modes
    benzo = sorted(float(f_c[k_]) for k_ in sel if own[k_] not in ('H3', 'H5'))
    assert len(benzo) == 4
    chv[p] = {'H3': float(f_c[pic_modes['H3'][0]]), 'H5': float(f_c[pic_modes['H5'][0]]), 'benzo': benzo}
    add(f'{p}.nu_C3_H3_stretch', chv[p]['H3'], 'cm-1', c_, RAW + c_ + '/hessian.npy', 'derived: the mode in 3000-3400 cm-1 with the largest H3 weight (picrate C-H)')
    add(f'{p}.nu_C5_H5_stretch', chv[p]['H5'], 'cm-1', c_, RAW + c_ + '/hessian.npy', 'derived: the mode in 3000-3400 cm-1 with the largest H5 weight (picrate C-H)')
    add(f'{p}.nu_benzo_CH_stretches', [round(v_, 1) for v_ in benzo], 'cm-1', c_, RAW + c_ + '/hessian.npy',
        'derived: the four C-H stretches of the benzo ring (H8-H11), coupled, ascending')
    add(f'{p}.benzo_ch_spread', benzo[-1] - benzo[0], 'cm-1', c_, RAW + c_ + '/hessian.npy', 'derived: highest minus lowest of the four benzo C-H stretches')
shifts = [chv['pcm']['H3'] - chv['gas']['H3'], chv['pcm']['H5'] - chv['gas']['H5']] + [a_ - b_ for a_, b_ in zip(chv['pcm']['benzo'], chv['gas']['benzo'])]
add('diff.ch_stretch_max_shift', float(max(abs(v_) for v_ in shifts)), 'cm-1', None, 'pcm/gas nu_C3_H3_stretch, nu_C5_H5_stretch, nu_benzo_CH_stretches',
    'derived: largest |PCM - gas| over the six C-H stretches (benzo modes paired in ascending order); shifts ' + ', '.join(f'{v_:+.1f}' for v_ in shifts))
# O7 offsets from the ELF/LOL planes (figure captions)
for c_ in CASES:
    p = PFX[c_]; Y = GEOM[p]
    z = np.load(RAW + c_ + '/elf_benzotriazolium_ring.npz')
    add(f'{p}.cation_plane_offset_O7', float(abs(np.dot(Y[L['O7']] - z['center'], z['normal']))), 'Å', c_, RAW + c_ + '/elf_benzotriazolium_ring.npz + opt.xyz',
        'derived: distance of the nitro O7 from the ELF/LOL cation plane (Fig. 12 marks its projection in the gas phase)')
    z = np.load(RAW + c_ + '/elf_picrate_ring.npz')
    add(f'{p}.picrate_plane_offset_O7', float(abs(np.dot(Y[L['O7']] - z['center'], z['normal']))), 'Å', c_, RAW + c_ + '/elf_picrate_ring.npz + opt.xyz',
        'derived: distance of the nitro O7 from the ELF/LOL picrate-ring plane (SI ELF/LOL figure)')
# optimisation thresholds as numbers (parsed from protocol.opt_criteria, which quotes optimize.log)
_oc = facts['protocol.opt_criteria']['value']
for k_, pat, u_ in (('energy', r'\|dE\| < ([0-9.eE+-]+)', 'Ha'), ('grad_rms', r'RMS grad < ([0-9.eE+-]+)', 'Ha/Bohr'),
                    ('grad_max', r'max grad < ([0-9.eE+-]+)', 'Ha/Bohr'), ('disp_rms', r'RMS disp < ([0-9.eE+-]+)', 'Å'),
                    ('disp_max', r'max disp < ([0-9.eE+-]+)', 'Å')):
    add(f'protocol.opt_threshold_{k_}', float(re.search(pat, _oc).group(1)), u_, None, facts['protocol.opt_criteria']['source'],
        'parsed from protocol.opt_criteria (geomeTRIC/Gaussian default convergence threshold)')
add('protocol.pcm_lebedev_points', 302, '', 'btpc_pcm', 'run_case.sh (PCM_LEBEDEV_ORDER 29); PySCF 2.14.0 pyscf/dft/gen_grid.py LEBEDEV_ORDER table',
    'number of Lebedev points of order 29 per atomic sphere of the PCM cavity (see protocol.pcm_lebedev_order)')
# lower end of the gas-phase N4-H4N...O1 NCI region (from the range recorded with gas.nci_region01)
_rng = re.search(r'range ([+-][0-9.]+)\.\.([+-][0-9.]+)', facts['gas.nci_region01']['note'])
add('gas.nci_region01_min', float(_rng.group(1)), 'a.u.', 'btpc_gas', facts['gas.nci_region01']['source'],
    'derived: most negative sign(l2)rho inside the gas-phase N4-H4N...O1 RDG < 0.5 region (range recorded with gas.nci_region01); close to the RDG_maxrho = 0.1 a.u. screen')
# plotting constants printed in captions (source: the figure scripts)
add('plot.nci_colour_min', -0.035, 'a.u.', None, 'papers/_shared/figkit rdg.NCI_RANGE (fig_nci.py VRANGE)', 'lower end of the linear sign(l2)rho colour scale (blue)')
add('plot.nci_colour_max', 0.02, 'a.u.', None, 'papers/_shared/figkit rdg.NCI_RANGE (fig_nci.py VRANGE)', 'upper end of the linear sign(l2)rho colour scale (red)')
addx('plot.nci_colour_mid', '(plot.nci_colour_min + plot.nci_colour_max) / 2', 'a.u.', None, 'the green middle of the linear NCI scale')
add('plot.mep_colour_min', -55, 'kcal mol-1', None, 'papers/btpc/figures/scripts/fig_mep.py (vr, outer multiples of 5 of all surface values)', 'lower end of the MEP colour scale, both models')
add('plot.mep_colour_max', 110, 'kcal mol-1', None, 'papers/btpc/figures/scripts/fig_mep.py (vr, outer multiples of 5 of all surface values)', 'upper end of the MEP colour scale, both models')
add('plot.mep_hist_bin', 2.5, 'kcal mol-1', None, 'papers/btpc/figures/scripts/fig_mep.py (bins)', 'bin width of the ESP surface-area histogram')
add('plot.mep_label_threshold', 30, 'kcal mol-1', None, 'papers/btpc/figures/scripts/fig_mep.py LABEL_ABOVE', 'extrema with |V| above this are labelled in the extrema strip')
add('plot.pdos_fwhm', 0.3, 'eV', None, 'papers/btpc/figures/scripts/fig_levels_pdos.py FWHM', 'Gaussian broadening (full width at half maximum) of the DOS')
add('plot.orbital_box_margin', 7.0, 'bohr', None, 'papers/btpc/figures/scripts/figS_orbitals_extra.py MARGIN_BOHR', 'margin of the orbital evaluation box beyond the outermost nuclei (3.70 Å)')


# =====================================================================================  8. revision round 1 (24 Sep 2026)
# Quantities the round-1 audits asked for (papers/btpc/audit/*_r1.md; response in audit/response_r1.md).  Local
# post-processing of stored files only: the CIF, opt_traj.xyz, hessian.npy, freq.json and final.molden.
R1 = []
def add8(key, value, unit='', case=None, source='', note=''):
    add(key, value, unit, case, source, note); facts[key]['added'] = 'r1'; R1.append(key)

# ---- 8a. CIF bond lengths with s.u. for the three Table 2 bonds that had none (figures_r1 #6)
for (a, b) in [('C2', 'C3'), ('C6', 'C5'), ('C7', 'C12')]:
    s8 = geom_bonds[(a, b)]
    add8(f'xray.bond_{a}_{b}', float(s8.split('(')[0]), 'Å', None, CIF + f':_geom_bond {a} {b}', f'CIF value with s.u. {s8}')
    facts[f'xray.bond_{a}_{b}']['su'] = s8
# ---- 8a2. versions of the packages that made the paper's figures and tables (citations_r1 #18: ParaView made only
#            the old hand-over renders; the paper figures come from the paper environment on the Raspberry Pi)
import importlib.metadata as _md
_pk = [('Matplotlib', 'matplotlib'), ('NumPy', 'numpy'), ('SciPy', 'scipy'), ('PyVista', 'pyvista'), ('VTK', 'vtk'),
       ('scikit-image', 'scikit-image'), ('gemmi', 'gemmi'), ('RDKit', 'rdkit')]
add8('software.figures_env', '; '.join(f'{n} {_md.version(d)}' for n, d in _pk), '', None,
     'tools/paper-env (installed package metadata, read with importlib.metadata when facts.json is built)',
     'environment of the figure and table scripts (papers/btpc/figures/scripts, tables/make_tables.py); the calculations used '
     'the versions in software.numpy_scipy_matplotlib on the compute host')
# ---- 8b. equivalent isotropic displacement parameters of the oxygen atoms (referee m12)
ueq = {}
for m8 in re.finditer(r'^(O\d)\s+O\s+\S+\s+\S+\s+\S+\s+(\d+\.\d+)\((\d+)\)\s+Uani', cif, re.M):
    ueq[m8.group(1)] = (float(m8.group(2)), f'{m8.group(2)}({m8.group(3)})')
assert set(ueq) == {f'O{i}' for i in range(1, 8)}, ueq
for o_ in ('O2', 'O3'):
    add8(f'xray.ueq_{o_}', ueq[o_][0], 'Å^2', None, CIF + ':_atom_site_U_iso_or_equiv', f'U(eq) {ueq[o_][1]} Å²; the two largest of all 33 atoms (nitro group N1 on C2)')
oth = sorted(v[0] for k8, v in ueq.items() if k8 not in ('O2', 'O3'))
add8('xray.ueq_O_other', [oth[0], oth[-1]], 'Å^2', None, CIF + ':_atom_site_U_iso_or_equiv', 'range of U(eq) over the other five O atoms (O1, O4-O7)')

# ---- 8c. trajectory quantities (figures_r1 #2, referee M6)
d8 = {nm_: np.array([dij(Yk, L[a_], L[b_]) for _, Yk in frames]) for nm_, (a_, b_) in
      {'H4N...O7': ('H4N', 'O7'), 'N4...O7': ('N4', 'O7')}.items()}
add8('saddle.h4n_o7_max_excursion', float(abs(d8['H4N...O7'] - d8['H4N...O7'][0]).max()), 'Å', 'btpc_pcm', RAW + 'btpc_pcm/opt_traj.xyz',
     f"derived: largest |d(step k) - d(step 0)| of H4N...O7 over the 71 frames (step 0 {d8['H4N...O7'][0]:.3f}, step 70 {d8['H4N...O7'][-1]:.3f}, "
     f"maximum {d8['H4N...O7'].max():.3f} Å); N4...O7 {float(abs(d8['N4...O7'] - d8['N4...O7'][0]).max()):.3f} Å")
ang8 = np.array([aijk(Yk, L['N4'], L['H4N'], L['O1']) for _, Yk in frames])
add8('saddle.hb_angle_step0', float(ang8[0]), 'deg', 'btpc_pcm', RAW + 'btpc_pcm/opt_traj.xyz', 'derived: angle N4-H4N...O1 at step 0 (inter-ring angle 5.8 deg)')
addx('saddle.hb_angle_opening', 'pcm.hb_N4_H4N_O1_angle - saddle.hb_angle_step0', 'deg', 'btpc_pcm',
     'opening of the N4-H4N...O1 angle during the IEF-PCM optimisation, while the rings turned from 5.8 to 58.2 deg')
facts['saddle.hb_angle_opening']['added'] = 'r1'; R1.append('saddle.hb_angle_opening')
nd8 = {k8: np.array([dij(Yk, L['N4'], L['H4N']) if k8 == 'NH' else dij(Yk, L['N4'], L['O1']) for _, Yk in frames]) for k8 in ('NH', 'NO')}
add8('saddle.nh_no_max_excursion', float(max(abs(nd8['NH'] - nd8['NH'][0]).max(), abs(nd8['NO'] - nd8['NO'][0]).max())), 'Å', 'btpc_pcm',
     RAW + 'btpc_pcm/opt_traj.xyz', 'derived: largest |d(step k) - d(step 0)| of N4-H4N and N4...O1 together over the 71 frames')

# ---- 8d. character of the modes described in the text and in Fig. 6 (numbers_r1 #9)
for c8 in CASES:
    p8 = PFX[c8]; s8_, x8, m8_, f8, V8, _ = MODES[c8]
    n4, h4 = L['N4'], L['H4N']
    k_n4 = int(np.argmin(abs(f8 - facts[f'{p8}.nu_N4_H4N_stretch']['value'])))
    disp = (V8[:, k_n4].reshape(-1, 3) / np.sqrt(m8_)[:, None])
    b8 = x8[h4] - x8[n4]; b8 /= np.linalg.norm(b8)
    wH = float((V8[:, k_n4].reshape(-1, 3)[h4] ** 2).sum())
    add8(f'{p8}.nu_N4_H4N_stretch_h_weight', wH, '', c8, RAW + c8 + '/hessian.npy',
         'derived: share of the mass-weighted eigenvector of the N4-H4N stretch carried by H4N')
    off = math.degrees(math.acos(abs(float(disp[h4] @ b8)) / np.linalg.norm(disp[h4])))
    add8(f'{p8}.nu_N4_H4N_stretch_offaxis_deg', off, 'deg', c8, RAW + c8 + '/hessian.npy',
         'derived: angle between the Cartesian displacement of H4N in the N4-H4N stretch and the N4-H4N bond (fig_vibrations.py check_offaxis_caption)')
    # rigid-body decomposition of the lowest mode, ion by ion (same least-squares fit as figures/scripts/_vib.py)
    q8 = V8[:, 0].reshape(-1, 3); mm8 = m8_; out8 = {}
    for nm8, idx8 in (('picrate', PIC), ('cation', CAT)):
        mi, xi = mm8[idx8], x8[idx8]; com8 = (mi[:, None] * xi).sum(0) / mi.sum(); r8 = xi - com8
        A8 = np.zeros((3 * len(idx8), 6))
        for a8 in range(len(idx8)):
            sq = np.sqrt(mi[a8]); rx = r8[a8]
            A8[3 * a8:3 * a8 + 3, :3] = sq * np.eye(3)
            A8[3 * a8:3 * a8 + 3, 3:] = sq * np.array([[0, rx[2], -rx[1]], [-rx[2], 0, rx[0]], [rx[1], -rx[0], 0]])
        sol8, *_ = np.linalg.lstsq(A8, q8[idx8].ravel(), rcond=None)
        ft, fr8 = A8[:, :3] @ sol8[:3], A8[:, 3:] @ sol8[3:]
        out8[nm8] = dict(t=float(ft @ ft), r=float(fr8 @ fr8), w=sol8[3:])
    lab8 = 'saddle.imag_mode' if c8 == 'btpc_pcm' else 'gas.mode1'
    for nm8 in ('picrate', 'cation'):
        tot8 = out8[nm8]['t'] + out8[nm8]['r']
        add8(f'{lab8}_{nm8}_rotation_share', out8[nm8]['r'] / tot8, '', c8, RAW + c8 + '/hessian.npy',
             f'derived: rotational part of the rigid-body motion of the {nm8} in the lowest mode (rotation / (rotation + translation) of its least-squares rigid fit); '
             f'translation {out8[nm8]["t"]:.3f}, rotation {out8[nm8]["r"]:.3f} of the whole mode')
    wp8, wc8 = out8['picrate']['w'], out8['cation']['w']
    ang_axes = math.degrees(math.acos(np.clip(wp8 @ wc8 / np.linalg.norm(wp8) / np.linalg.norm(wc8), -1, 1)))
    add8(f'{lab8}_axes_angle_deg', ang_axes, 'deg', c8, RAW + c8 + '/hessian.npy',
         'derived: angle between the rotation axes (angular velocities) of the two ions in the lowest mode; > 90 deg = the ions turn in opposite senses')
    link8 = x8[cat9].mean(0) - x8[pring].mean(0)
    a_cl = math.degrees(math.acos(np.clip(abs(wc8 @ link8) / np.linalg.norm(wc8) / np.linalg.norm(link8), -1, 1)))
    add8(f'{lab8}_cation_axis_to_centroid_line_deg', a_cl, 'deg', c8, RAW + c8 + '/hessian.npy + opt.xyz',
         'derived: angle between the rotation axis of the cation in the lowest mode and the line joining the picrate and cation ring centroids')

for k8, e8 in (('gas.nu_N4_H4N_stretch_h_pct', '100 * gas.nu_N4_H4N_stretch_h_weight'),
               ('pcm.nu_N4_H4N_stretch_h_pct', '100 * pcm.nu_N4_H4N_stretch_h_weight'),
               ('saddle.imag_mode_picrate_rotation_pct', '100 * saddle.imag_mode_picrate_rotation_share'),
               ('saddle.imag_mode_cation_rotation_pct', '100 * saddle.imag_mode_cation_rotation_share'),
               ('gas.mode3_rigid_pct', '100 * gas.mode3_rigid_fraction'), ('pcm.mode3_rigid_pct', '100 * pcm.mode3_rigid_fraction')):
    addx(k8, e8, '%', 'btpc_gas' if k8.startswith('gas') else 'btpc_pcm', 'the same share in per cent')
    facts[k8]['added'] = 'r1'; R1.append(k8)

# ---- 8e. quasi-RRHO entropy (Grimme's free-rotor interpolation for the soft modes; referee m7)
KB, HPL, CLIGHT, NA, RGAS = 1.380649e-23, 6.62607015e-34, 2.99792458e10, 6.02214076e23, 8.314462618
TK, PPA, BAV, NU0 = 298.15, 101325.0, 1.0e-44, 100.0
def s_harm(nu):
    x = HPL * CLIGHT * nu / (KB * TK); return RGAS * (x / math.expm1(x) - math.log(-math.expm1(-x)))
def s_rot_free(nu):
    mu = HPL / (8 * math.pi ** 2 * CLIGHT * nu); mup = mu * BAV / (mu + BAV)
    return RGAS * (0.5 + math.log(math.sqrt(8 * math.pi ** 3 * mup * KB * TK / HPL ** 2)))
for c8 in CASES:
    p8 = PFX[c8]; s8_, x8, m8_, f8, V8, _ = MODES[c8]
    mkg = np.array([MASS[e] for e in s8_]) * 1.66053906660e-27; xm = x8 * BOHR * 1e-10
    com8 = (mkg[:, None] * xm).sum(0) / mkg.sum(); r8 = xm - com8
    I8 = np.linalg.eigvalsh(sum(mi * (np.dot(ri, ri) * np.eye(3) - np.outer(ri, ri)) for mi, ri in zip(mkg, r8)))
    M8 = mkg.sum()
    s_tr = RGAS * (math.log((2 * math.pi * M8 * KB * TK / HPL ** 2) ** 1.5 * KB * TK / PPA) + 2.5)
    s_ro = RGAS * (math.log(math.sqrt(math.pi * np.prod(8 * math.pi ** 2 * I8 * KB * TK / HPL ** 2))) + 1.5)
    nus = [float(v) for v in freq[c8]['frequencies_cm-1'] if v > 0]
    s_v = sum(s_harm(v) for v in nus)
    s_q = sum((1 / (1 + (NU0 / v) ** 4)) * s_harm(v) + (1 - 1 / (1 + (NU0 / v) ** 4)) * s_rot_free(v) for v in nus)
    S_rrho = (s_tr + s_ro + s_v) / 4.184; S_q = (s_tr + s_ro + s_q) / 4.184
    assert abs(S_rrho - facts[f'{p8}.entropy']['value']) < 0.05, (c8, S_rrho, facts[f'{p8}.entropy']['value'])
    add8(f'{p8}.entropy_rrho_check', S_rrho, 'cal mol-1 K-1', c8, T + c8 + '/freq.json + opt.xyz',
         f'derived: translational + rotational + harmonic vibrational entropy recomputed here (sigma = 1, {len(nus)} real modes); reproduces {p8}.entropy (from H298 - G298)')
    add8(f'{p8}.entropy_qrrho', S_q, 'cal mol-1 K-1', c8, T + c8 + '/freq.json + opt.xyz',
         'derived: quasi-RRHO entropy, each vibrational term interpolated between the harmonic and a free-rotor entropy with the weight '
         '1/(1 + (100 cm-1 / nu)^4) and B_av = 1e-44 kg m^2 (Grimme 2012); translation and rotation as in the RRHO value')
    addx(f'{p8}.entropy_qrrho_lowering', f'{p8}.entropy - {p8}.entropy_qrrho', 'cal mol-1 K-1', c8, 'RRHO minus quasi-RRHO entropy')
    facts[f'{p8}.entropy_qrrho_lowering']['added'] = 'r1'; R1.append(f'{p8}.entropy_qrrho_lowering')

# ---- 8f. Loewdin populations beside the Mulliken ones (referee M4)
for c8 in CASES:
    p8 = PFX[c8]
    atoms8, shells8, mos8 = parse_molden(RAW + c8 + '/final.molden')
    S8, aof8 = overlap(atoms8, shells8)
    nreal8 = sum(1 for mo in mos8 if np.abs(mo['c']).sum() > 0)
    C8 = np.array([mo['c'] for mo in mos8]).T[:, :nreal8]; occ8 = np.array([mo.get('occ', 0) for mo in mos8])[:nreal8]
    D8 = (C8 * occ8) @ C8.T
    ev8, U8 = np.linalg.eigh(S8); Sh8 = (U8 * np.sqrt(ev8)) @ U8.T
    popM = np.einsum('ij,ji->i', D8, S8); popL = np.diag(Sh8 @ D8 @ Sh8)
    assert abs(popM.sum() - 178) < 1e-6 and abs(popL.sum() - 178) < 1e-6, (popM.sum(), popL.sum())
    Zc = {'O': 8, 'N': 7, 'C': 6, 'H': 1}
    qat = lambda pop, a: Zc[SYM[a]] - pop[aof8 == a].sum()
    qpicL = sum(qat(popL, a) for a in PIC)
    add8(f'{p8}.lowdin_charge_picrate', float(qpicL), 'e', c8, RAW + c8 + '/final.molden',
         f'derived: Loewdin fragment charge, diag(S^1/2 D S^1/2), overlap rebuilt from the molden basis (smallest overlap eigenvalue {ev8.min():.1e}); cation = -this')
    add8(f'{p8}.mulliken_charge_H4N', float(qat(popM, L['H4N'])), 'e', c8, RAW + c8 + '/final.molden', 'derived: Mulliken atomic charge of H4N (counted with the cation)')
    add8(f'{p8}.lowdin_charge_H4N', float(qat(popL, L['H4N'])), 'e', c8, RAW + c8 + '/final.molden', 'derived: Loewdin atomic charge of H4N (counted with the cation)')
for k8 in ('lowdin_charge_picrate', 'mulliken_charge_H4N', 'lowdin_charge_H4N'):
    addx(f'diff.{k8}', f'pcm.{k8} - gas.{k8}', 'e', None, 'IEF-PCM minus gas phase')
    facts[f'diff.{k8}']['added'] = 'r1'; R1.append(f'diff.{k8}')

# ---- 8g. density topology of the N4-H4N...O1 contact (referee M5): bond critical point and ELF, from final.molden
sys.path.insert(0, str(OUT))
import wfn_topology as WT
for c8 in CASES:
    p8 = PFX[c8]; w8 = WT.Wfn(RAW + c8 + '/final.molden'); X8 = w8.mf.coords
    bc = w8.bcp(X8[L['H4N']], X8[L['O1']])
    src8 = RAW + c8 + '/final.molden (papers/btpc/data/wfn_topology.py)'
    reg = facts[f'{p8}.nci_region01']['value']
    assert abs(bc['sign_l2_rho'] - reg) < 0.002 and bc['lambda'][1] < 0 < bc['lambda'][2], (c8, bc, reg)
    add8(f'{p8}.bcp_H4N_O1_rho', bc['rho'], 'a.u.', c8, src8,
         f'derived: electron density at the (3,-1) critical point between H4N and O1 (Newton search from the density minimum on the H4N-O1 line; |grad rho| {bc["grad"]:.0e}); '
         f'lambda = {bc["lambda"][0]:.4f}, {bc["lambda"][1]:.4f}, {bc["lambda"][2]:.4f} a.u.; sign(lambda2)rho agrees with the lowest-RDG point of NCI region 01 ({reg:+.4f}) to 0.002')
    add8(f'{p8}.bcp_H4N_O1_lap', bc['lap'], 'a.u.', c8, src8, 'derived: Laplacian of the density at that critical point (e bohr^-5)')
    add8(f'{p8}.bcp_H4N_O1_G', bc['G'], 'a.u.', c8, src8, 'derived: positive-definite kinetic energy density G = 1/2 sum n_i |grad phi_i|^2 (hartree bohr^-3)')
    add8(f'{p8}.bcp_H4N_O1_V', bc['V'], 'a.u.', c8, src8, 'derived: potential energy density from the local virial theorem, V = lap/4 - 2G')
    add8(f'{p8}.bcp_H4N_O1_H', bc['H'], 'a.u.', c8, src8, 'derived: total energy density H = G + V; H < 0 with lap > 0 marks a partly covalent (shared-shell) contribution')
    add8(f'{p8}.bcp_H4N_O1_VG', bc['VG'], '', c8, src8, 'derived: |V|/G; 1 < |V|/G < 2 is the intermediate (transit) region between closed-shell and shared interactions')
    add8(f'{p8}.bcp_H4N_O1_elf', bc['elf'], '', c8, src8, 'derived: ELF (Becke, with the 1e-5 addition of Multiwfn ELF_addminimal) at the critical point')
    add8(f'{p8}.bcp_H4N_O1_dist_H', float(np.linalg.norm(bc['x'] - X8[L['H4N']]) * BOHR), 'Å', c8, src8, 'derived: distance of the critical point from H4N')
    # ELF profile along the H4N -> O1 line: its minimum between the H4N and O1 valence regions
    t8 = np.linspace(0.05, 0.45, 161); seg8 = X8[L['H4N']][None] + t8[:, None] * (X8[L['O1']] - X8[L['H4N']])[None]
    e8 = w8.elf(seg8); k8 = int(np.argmin(e8))
    add8(f'{p8}.elf_min_H4N_O1', float(e8[k8]), '', c8, src8,
         f'derived: minimum of ELF on the straight H4N -> O1 line between the H4N and O1 valence regions, {t8[k8] * np.linalg.norm(X8[L["O1"]] - X8[L["H4N"]]) * BOHR:.3f} Å from H4N')
    # validation of the evaluator against the stored PySCF cube and the Multiwfn ELF plane (written to the notes)
    org8, sp8, n8, cube8 = read_cube(RAW + c8 + '/density.cube')
    rng8 = np.random.default_rng(7); idx8 = np.argwhere(cube8 > 1e-3); idx8 = idx8[rng8.choice(len(idx8), 400, replace=False)]
    rdev = float(np.abs(w8.local(org8 + idx8 @ sp8)[0] / cube8[tuple(idx8.T)] - 1).max())
    pl8 = np.loadtxt(RAW + c8 + '/elf_benzotriazolium_ring_plane.txt'); sel8 = pl8[rng8.choice(len(pl8), 400, replace=False)]; sel8 = sel8[sel8[:, 5] > 0.02]
    edev = float(np.abs(w8.elf(sel8[:, :3] / BOHR) - sel8[:, 5]).max())
    assert rdev < 1e-3 and edev < 1e-3, (c8, rdev, edev)
    add8(f'{p8}.wfn_rho_rel_dev_max', rdev, '', c8, src8 + ' vs ' + RAW + c8 + '/density.cube',
         'derived: largest relative deviation of the evaluated density from the PySCF cube over 400 random points with rho > 1e-3 a.u.')
    add8(f'{p8}.wfn_elf_abs_dev_max', edev, '', c8, src8 + ' vs ' + RAW + c8 + '/elf_benzotriazolium_ring_plane.txt',
         'derived: largest absolute deviation of the evaluated ELF from the Multiwfn plane map (random plane points with ELF > 0.02)')
    add8(f'{p8}.bcp_nci_sign_l2_rho_diff', abs(bc['sign_l2_rho'] - reg), 'a.u.', c8, src8 + '; ' + f'{p8}.nci_region01',
         'derived: |sign(lambda2)rho at the critical point - value at the lowest-RDG point of NCI region 01|')
    add8(f'{p8}.wfn_topology_validation', f'density vs density.cube: max relative deviation {rdev:.1e} (400 points with rho > 1e-3 a.u.); '
         f'ELF vs Multiwfn elf_benzotriazolium_ring_plane.txt: max absolute deviation {edev:.1e} ({len(sel8)} points)', '', c8, src8,
         'validation of the orbital evaluator used for the bcp_* and elf_min_* keys')
print('r1 additions:', len(R1))

# =====================================================================================  9. revision round 2 (24 Sep 2026)
# Quantities the round-2 audits asked for (papers/btpc/audit/*_r2.md; response in audit/response_r2.md).  Local
# post-processing of stored files only: hessian.npy + opt.xyz, the CIF, dft/_fd_worker.py and the figure scripts.
R2 = []
def add9(key, value, unit='', case=None, source='', note=''):
    add(key, value, unit, case, source, note); facts[key]['added'] = 'r2'; R2.append(key)
def addx9(key, expr, unit='', case=None, note=''):
    addx(key, expr, unit, case, note); facts[key]['added'] = 'r2'; R2.append(key)

# ---- 9a. numerical resolution of the finite-difference Hessians (referee_r2 M1): curvature of the stored,
#          unprojected mass-weighted Hessian along the overall rigid rotations and translations of the pair.
#          An exact Hessian at a stationary point gives zero for both.
for c9 in CASES:
    p9 = PFX[c9]; s9, xyz9 = read_xyz(T + c9 + '/opt.xyz'); x9 = xyz9 / BOHR
    H9 = np.load(RAW + c9 + '/hessian.npy'); n9 = len(s9); m9 = np.array([MASS[e] for e in s9]) * AMU2ME
    Hm9 = H9.transpose(0, 2, 1, 3).reshape(3 * n9, 3 * n9); w9 = np.repeat(m9 ** -0.5, 3); Hm9 = Hm9 * w9[:, None] * w9[None, :]
    Hm9 = 0.5 * (Hm9 + Hm9.T)
    r9 = x9 - (m9[:, None] * x9).sum(0) / m9.sum()
    Tb, Rb = np.zeros((3, 3 * n9)), np.zeros((3, 3 * n9))
    for i9 in range(n9):
        sq9 = np.sqrt(m9[i9]); Tb[:, 3 * i9:3 * i9 + 3] = sq9 * np.eye(3)
        for k9, ax9 in enumerate(np.eye(3)): Rb[k9, 3 * i9:3 * i9 + 3] = sq9 * np.cross(ax9, r9[i9])
    for nm9, B9 in (('rot', Rb), ('trans', Tb)):
        Q9, _ = np.linalg.qr(B9.T); ev9 = np.linalg.eigvalsh(Q9.T @ Hm9 @ Q9)
        cm9 = [round(float(np.sign(e9) * np.sqrt(abs(e9)) * AU2CM), 2) for e9 in ev9]
        add9(f'{p9}.hessian_{nm9}_curvature_cm', cm9, 'cm-1', c9, RAW + c9 + '/hessian.npy + ' + T + c9 + '/opt.xyz',
             f'derived: eigenvalues of Q^T M^-1/2 H M^-1/2 Q for the unprojected Hessian, Q = the three orthonormalised rigid '
             f'{"rotations about the centre of mass" if nm9 == "rot" else "translations"} of the whole pair, as wavenumbers (negative = imaginary); '
             'zero for an exact Hessian at a stationary point, so the spread measures the numerical resolution of the Hessian '
             '(integration grid not rotated with the molecule, gradients without grid-weight response, residual gradient)')
    add9(f'{p9}.hessian_rot_curvature_max_abs', float(max(abs(v9) for v9 in facts[f'{p9}.hessian_rot_curvature_cm']['value'])), 'cm-1', c9,
         f'{p9}.hessian_rot_curvature_cm', 'derived: largest |value| of the three rotational curvatures')
add9('protocol.fd_gradient_grid_response', 'no: nuc_grad_method().kernel() with PySCF defaults (grid_response = False), displaced SCF conv_tol 1e-10',
     '', None, 'dft/_fd_worker.py', 'the analytic gradients behind both finite-difference Hessians omit the response of the integration-grid weights')

# ---- 9b. energy density per electron at the H4N...O1 critical point, and the gas/IEF-PCM ratio of H (referee_r2 M2, numbers_r2 #3)
for p9 in ('gas', 'pcm'):
    addx9(f'{p9}.bcp_H4N_O1_H_over_rho', f'{p9}.bcp_H4N_O1_H / {p9}.bcp_H4N_O1_rho', 'a.u.', f'btpc_{p9}',
          'H/rho at the H4N...O1 critical point, hartree per electron (Espinosa et al. 2002 use it to grade covalency within the transit region)')
addx9('diff.bcp_H_ratio_gas_pcm', 'gas.bcp_H4N_O1_H / pcm.bcp_H4N_O1_H', '', None, 'how many times more negative H is in the gas phase')

# ---- 9c. displacement parameters of the triazole nitrogens (referee_r2 m5)
for m9_ in re.finditer(r'^(N[456])\s+N\s+\S+\s+\S+\s+\S+\s+(\d+\.\d+)\((\d+)\)\s+Uani', cif, re.M):
    add9(f'xray.ueq_{m9_.group(1)}', float(m9_.group(2)), 'Å^2', None, CIF + ':_atom_site_U_iso_or_equiv',
         f'U(eq) {m9_.group(2)}({m9_.group(3)}) Å²; triazole nitrogen of the cation (N5 is the unprotonated central atom)')
assert all(f'xray.ueq_{a9}' in facts for a9 in ('N4', 'N5', 'N6'))

# ---- 9d. plotting constants printed in captions (numbers_r2 #14); read from the scripts, not typed
def _const9(path, pat):
    m9c = re.search(pat, open(paths.resolve(path)).read(), re.M); assert m9c, (path, pat); return float(m9c.group(1))
FS9 = 'papers/btpc/figures/scripts/'
# plot.crystal_stack_tilt: not needed since the Fig. 2 caption no longer prints the tilt (it paired 38 deg rotations in
# the text with the plotting constant)
add9('plot.levels_degen_ev', _const9(FS9 + 'fig_levels_pdos.py', r'^DEGEN\s*=\s*([0-9.]+)'), 'eV', None, FS9 + 'fig_levels_pdos.py DEGEN',
     'levels closer than this are drawn side by side')
add9('plot.bond_band', _const9(FS9 + 'fig_xray_comparison.py', r'parity_plot\(.*band=([0-9.]+)'), 'Å', None,
     FS9 + 'fig_xray_comparison.py (parity_plot band); figS_bond_angle_deviations.py (specs)', 'half-width of the grey band in the bond-length plots')
add9('plot.angle_band', _const9(FS9 + 'figS_bond_angle_deviations.py', r"\('b',.*?,\s*([0-9.]+),\s*'°'"), 'deg', None,
     FS9 + 'figS_bond_angle_deviations.py (specs)', 'half-width of the grey band in the angle plot')
# ---- 9e. the Delta pKa range printed in the abstract, Fig. 1 and the Introduction (numbers_r2 #5, figures_r2 #9);
#          parsed from the abstract-checked literature claims (papers/btpc/lit/candidates.jsonl), never typed in
_cl9 = {}
for _l9 in open('papers/btpc/lit/candidates.jsonl'):
    if _l9.strip():
        _d9 = json.loads(_l9); _cl9[_d9['key']] = ' '.join(_d9.get('claims', []))
_pk_acid = float(re.search(r'determined spectrophotometrically as (\d\.\d+)', _cl9['davis1963spectrophotometric']).group(1))
_acid_lo, _acid_hi = map(float, re.search(r'spanning roughly (\d\.\d+)-(\d\.\d+)', _cl9['davis1963spectrophotometric']).groups())
_pk_new = float(re.search(r'pKa = (\d\.\d+)', _cl9['wang2000photochemistry']).group(1))
_pk_old = float(re.search(r'Basicity pKa of benzotriazole (\d\.\d+)', _cl9['alkorta2012theoretical']).group(1))
_src9 = 'papers/btpc/lit/candidates.jsonl (claims of wang2000photochemistry, davis1963spectrophotometric, alkorta2012theoretical)'
add9('lit.pk_picric_acid', _pk_acid, '', None, _src9, 'thermodynamic pK of picric acid in water at 25 C (davis1963spectrophotometric)')
add9('lit.pk_picric_acid_spread', [_acid_lo, _acid_hi], '', None, _src9, 'spread of earlier picric acid pK determinations (davis1963spectrophotometric)')
add9('lit.pka_benzotriazolium_flash', _pk_new, '', None, _src9, 'pKa of the N-protonated benzotriazolium cation, flash photolysis, I = 1.0 M (wang2000photochemistry)')
add9('lit.pka_benzotriazolium_older', _pk_old, '', None, _src9, 'basicity pKa of benzotriazole from the older review quoted by Alkorta et al. (alkorta2012theoretical)')
add9('lit.dpka_picric_acid_flash', _pk_new - _pk_acid, '', None, _src9,
     f'derived: pKa(BTAH+) {_pk_new} (flash photolysis, I = 1.0 M) - pK(picric acid) {_pk_acid} (25 C)')
add9('lit.dpka_picric_acid_flash_range', [round(_pk_new - _acid_hi, 2), round(_pk_new - _acid_lo, 2)], '', None, _src9,
     f'derived: {_pk_new} minus the spread {_acid_lo}-{_acid_hi} of published picric acid values')
add9('lit.dpka_picric_acid_older', _pk_old - _pk_acid, '', None, _src9,
     f'derived: basicity pKa {_pk_old} of the older review quoted by Alkorta et al. - pK(picric acid) {_pk_acid}')
add9('lit.dpka_picric_acid_older_range', [round(_pk_old - _acid_hi, 2), round(_pk_old - _acid_lo, 2)], '', None, _src9,
     f'derived: {_pk_old} minus the spread {_acid_lo}-{_acid_hi}')
add9('lit.dpka_picric_acid_range', [round(_pk_new - _acid_hi, 2), round(_pk_old - _acid_lo, 2)], '', None, _src9,
     'derived: lowest (flash-photolysis value, highest picric acid value) to highest (older value, lowest picric acid '
     'value) Delta pKa over the published values; printed as "about -0.4 to 1.4"')
# ---- 9f. finite-difference steps as numbers (the text quotes them; before they matched only by accident)
for c9 in CASES:
    p9 = PFX[c9]
    add9(f'{p9}.hessian_step_bohr', float(re.search(r'step ([\d.]+) Bohr', facts[f'{p9}.hessian_method']['value']).group(1)),
         'bohr', c9, f'{p9}.hessian_method', 'finite-difference displacement step of the reported Hessian')
add9('saddle.hessian_step_first_bohr', float(re.search(r'step ([\d.]+) Bohr', facts['saddle.freq_imag_step0p005']['note']).group(1)),
     'bohr', 'btpc_pcm', 'saddle.freq_imag_step0p005', 'step of the first IEF-PCM Hessian (the check)')
print('r2 additions:', len(R2))

json.dump(facts, open(OUT / 'facts.json', 'w'), indent=1, ensure_ascii=False)
print(f'wrote {len(facts)} facts to {OUT / "facts.json"}')
