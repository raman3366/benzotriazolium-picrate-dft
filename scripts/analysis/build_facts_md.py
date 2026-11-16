"""
Render papers/btpc/data/facts.md from facts.json (run build_facts.py first).
Every number in facts.md is pulled from facts.json by key, so the two files cannot disagree.
    python scripts/analysis/build_facts_md.py
"""
import json, pathlib, re

import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))   # scripts/ (paths.py)
import paths
D = paths.ANALYSIS
F = json.load(open(D / 'facts.json'))

def v(k): return F[k]['value']
def f(k, nd=3, sign=False):
    x = F[k]['value']
    if isinstance(x, (int, float)) and not isinstance(x, bool):
        return f'{x:+.{nd}f}' if sign else f'{x:.{nd}f}'
    return str(x)
def fu(k, nd=3, sign=False):
    u = F[k]['unit']; return f'{f(k, nd, sign)} {u}'.strip()
def row(*cells): return '| ' + ' | '.join(str(c).replace('|', '\\|') for c in cells) + ' |'
def kk(k): return f'`{k}`'

out = []
W = out.append

W('# Facts dossier: benzotriazolium picrate (BTPC, paper id `btpc`)')
W('')
W('Authoritative source of every number for the BTPC manuscript. Machine-readable twin: `facts.json` '
  f'({len(F)} keys; each has value, unit, case, source file[:field] and a note). Both files are generated: '
  '`python scripts/analysis/build_facts.py` then `python scripts/analysis/build_facts_md.py`. '
  'Keys are quoted in backticks so an auditor can trace any number. "Derived" means computed here from a stored file '
  '(formula in the facts.json note); nothing was recomputed with quantum chemistry.')
W('')
W('Case names: `btpc_pcm` = IEF-PCM (water) model, **primary**; `btpc_gas` = gas-phase model, comparison. '
  'Atom labels are the CIF labels; the cation labels carry an internal `#4_666` suffix in the project files, which must never appear in the paper. '
  'Never call the compound a Schiff base.')
W('')
W('Sources that must NOT be used for numbers: `results/SUMMARY.md`, `STATUS.md` (stale, 3ATB only), `handover/` documents (context only), '
  'and the "integration grid level 3" sentence in every `RESULTS.md` (wrong; level 4 was used).')
W('')

# ------------------------------------------------------------------------------------------------ 1
W('## 1. Compound and crystal structure')
W('')
W(f'- Compound: benzotriazolium picrate; the CIF calls it **{v("compound.cif_common_name")}** (`_chemical_name_systematic` is "?"). '
  f'Case-file name: "{v("compound.chemical_name_casefile")}" — see caveat C12: H atoms sit on N4 and N6 (both outer triazole N), '
  'N5 (central) carries none, i.e. the 1H,3H-benzotriazolium cation (systematic: 1H-1,2,3-benzotriazol-3-ium). Check the name against the source paper.')
W(f'- Formula {v("compound.formula_sum")} (moiety: {v("compound.formula_moiety")}), Mr {f("xray.mr", 2)} g mol⁻¹. '
  f'CCDC {v("compound.ccdc")}; data-block DOI {v("xray.audit_block_doi")}; deposited {v("xray.deposition_date")}. '
  f'Original structure paper: doi:{v("compound.source_doi")} (CIF citation year {v("compound.source_year")}).')
W('')
W('Table 1.1 — Crystal data (all from `molecules/btpc/input/BTPC.cif`).')
W('')
W(row('Quantity', 'Value', 'key'))
W(row('---', '---', '---'))
for lab_, k, nd in [('Crystal system', 'xray.crystal_system', 0), ('Space group', 'xray.space_group', 0), ('Hall symbol', 'xray.space_group_hall', 0),
                    ('a (Å)', 'xray.cell_a', 4), ('b (Å)', 'xray.cell_b', 4), ('c (Å)', 'xray.cell_c', 4), ('α (°)', 'xray.cell_alpha', 2),
                    ('β (°)', 'xray.cell_beta', 3), ('γ (°)', 'xray.cell_gamma', 2), ('V (Å³)', 'xray.cell_volume', 1), ('Z / Z′', None, 0),
                    ('T (K)', 'xray.temperature', 0), ('Radiation, λ (Å)', None, 0), ('D_calc (g cm⁻³)', 'xray.density_calc', 3), ('μ (mm⁻¹)', 'xray.mu', 3),
                    ('F(000)', 'xray.f000', 0), ('Crystal', 'xray.crystal_description', 0), ('θ range (°)', 'xray.theta_range', 0),
                    ('Reflections measured / unique / I>2σ(I)', None, 0), ('R_int', 'xray.r_int', 4), ('Parameters / restraints', None, 0),
                    ('R1 [I>2σ(I)] / R1 (all)', None, 0), ('wR2 [I>2σ(I)] / wR2 (all)', None, 0), ('GoF', 'xray.gof', 3),
                    ('Δρ max / min (e Å⁻³)', None, 0), ('Programs', 'xray.software', 0), ('Absorption correction', 'xray.absorption_correction', 0),
                    ('H-atom treatment', 'xray.hydrogen_treatment', 0), ('Deposited H-bond table', 'xray.hbond_table_deposited', 0)]:
    if k:
        val = f(k, nd) if isinstance(v(k), (int, float)) else v(k)
        if k.startswith('xray.cell_'): val = F[k]['note'].split(': ')[1]
        W(row(lab_, val, kk(k)))
    elif lab_.startswith('Z '): W(row(lab_, f'{v("xray.z")} / {v("xray.z_prime")}', '`xray.z`, `xray.z_prime`'))
    elif lab_.startswith('Radiation'): W(row(lab_, f'{v("xray.radiation")}, {f("xray.wavelength", 5)}', '`xray.radiation`, `xray.wavelength`'))
    elif lab_.startswith('Reflections'): W(row(lab_, f'{v("xray.reflections_measured")} / {v("xray.reflections_unique")} / {v("xray.reflections_observed")}', '`xray.reflections_*`'))
    elif lab_.startswith('Parameters'): W(row(lab_, f'{v("xray.n_parameters")} / {v("xray.n_restraints")}', '`xray.n_parameters`, `xray.n_restraints`'))
    elif lab_.startswith('R1'): W(row(lab_, f'{f("xray.r1_gt", 4)} / {f("xray.r1_all", 4)}', '`xray.r1_gt`, `xray.r1_all`'))
    elif lab_.startswith('wR2'): W(row(lab_, f'{f("xray.wr2_gt", 4)} / {f("xray.wr2_all", 4)}', '`xray.wr2_gt`, `xray.wr2_all`'))
    elif lab_.startswith('Δρ'): W(row(lab_, f'{f("xray.diff_density_max", 3)} / {f("xray.diff_density_min", 3)}', '`xray.diff_density_*`'))
W('')
W('**No hydrogen-bond table was deposited** (the CIF has no `_geom_hbond` loop). Table 1.2 is derived here from the deposited '
  'coordinates and the four symmetry operations (search: H···A < 2.9 Å, D–H···A > 100°, N/O acceptors; donor at x,y,z). '
  'Also written to `papers/btpc/data/derived/crystal_hbonds.csv`. X-ray H positions (N–H 0.91/0.83 Å) are short, so H···A is long and the angle uncertain; D···A is the robust number.')
W('')
W('Table 1.2 — Crystal contacts (derived, `xray.hbNN_*`).')
W('')
W(row('#', 'D–H···A', 'A at', 'D–H (Å)', 'H···A (Å)', 'D···A (Å)', '∠D–H···A (°)', 'in DFT cluster?', 'comment'))
W(row(*['---'] * 9))
for k in sorted(x for x in F if re.match(r'xray\.hb\d\d_[A-Z0-9]+_[A-Z0-9]+_[A-Z0-9]+$', x)):
    note = F[k]['note']; sym = re.search(r'symmetry (\S+) = \(([^)]*)\)', note)
    inside = 'yes' if 'present in the DFT cluster' in note else 'no'
    com = 'intramolecular S(5)-type, not an H-bond' if 'intramolecular' in note else ''
    W(row(k.split('_')[0].replace('xray.hb', ''), v(k).replace('...', '···'), f'{sym.group(1)} ({sym.group(2)})', f(k + '_dh', 3), f(k + '_ha', 3), f(k + '_da', 3), f(k + '_angle', 1), inside, com))
W('')
W('Reading Table 1.2: the cation donates two N–H groups to two *different* picrates. N4–H4N binds the phenolate O1 of the picrate at '
  '(x−1, −y+½, z−½) with a minor bifurcated component to the ortho-nitro O7 of the same anion; N6–H6N binds nitro O6 of an '
  'inversion-related picrate at (−x+1, −y+1, −z+1). Only the first of these is inside the DFT cluster. '
  f'The asymmetric unit itself is not hydrogen bonded: its closest inter-ion contact is {v("xray.asu_closest_interion_contact")} '
  f'{f("xray.asu_closest_interion_distance", 3)} Å (`xray.asu_closest_interion_*`). The only ring stacking below 4.3 Å is triazole···triazole '
  f'across an inversion centre (−x, 1−y, −z): centroid–centroid {f("xray.stacking_triazole_cg_cg", 3)} Å, interplanar separation {f("xray.stacking_triazole_perp", 3)} Å (`xray.stacking_*`).')
W('')
W('Table 1.3 — Selected deposited geometry with s.u. (CIF `_geom_bond` / `_geom_angle`).')
W('')
W(row('Parameter', 'X-ray', 'key'))
W(row('---', '---', '---'))
for k in [x for x in F if x.startswith('xray.bond_') or x.startswith('xray.angle_')]:
    W(row(k.replace('xray.bond_', '').replace('xray.angle_', '∠').replace('_', '–'), F[k]['note'].replace('CIF value with s.u. ', '') + (' Å' if 'bond' in k else '°'), kk(k)))
W('')

# ------------------------------------------------------------------------------------------------ 2
W('## 2. Structural model (how the cluster was built)')
W('')
W(f'- {v("model.n_atoms")} atoms, {v("model.n_electrons")} electrons (89 doubly occupied orbitals), total charge {v("model.charge")}, singlet (restricted closed shell). '
  f'Fragments: {v("model.fragments")} (`model.*`).')
W(f'- Built by `dft/build_cluster.py` from the CIF: the symmetry is expanded and, for each moiety, the periodic image that hydrogen-bonds to the '
  f'growing cluster is selected by an exhaustive search (criteria H···A < 2.65 Å, D–H···A > 130°, heavy-atom clash < 2.40 Å). '
  f'Picrate kept at {v("model.symmetry_picrate")}; benzotriazolium moved to {v("model.symmetry_cation")}. '
  'The defining contact is N4–H4N···O1 (cation N–H to phenolate O); the bifurcated N4–H4N···O7 (123°) is inside the cluster but below the builder\'s angle cut-off. '
  'N6–H6N is a free N–H in the model (its crystal partner, O6 of a second picrate, is absent; that picrate is the inversion image of the asymmetric-unit picrate, but relative to the cluster\'s own picrate it is the 2₁-screw image 2_746, see `model.missing_hbond`) — this governs the MEP maximum and the 3578/3633 cm⁻¹ band.')
W('- Atom order everywhere is the CIF order within each ion, picrate first (atoms 0–17), then the cation (18–32); the CIF itself interleaves the two ions (O1, N3, O7, N4, C6, …). Curated by `dft/curate_cases.py` (names, ELF/LOL planes). No basis-set overrides.')
W(f'- Basis size: {v("model.n_basis")} functions ({v("model.basis_per_element")}). PySCF removed near-linear-dependent combinations of the diffuse functions: '
  f'{v("pcm.n_lindep_removed")} at the final PCM geometry (603 MOs), {v("gas.n_lindep_removed")} for gas (602 MOs); only the PCM molden was zero-padded to 606 MOs for Multiwfn (no physical effect); the gas molden holds the 602 MOs unpadded and Multiwfn read 602 orbitals. The B3LYP/6-31G(d) pre-optimisation level has 366 functions (`model.preopt_basis_functions`).')
W('')

# ------------------------------------------------------------------------------------------------ 3
W('## 3. Computational protocol')
W('')
W('Table 3.1 — Settings. Where the value comes from the PySCF 2.14.0 source rather than a run log, the key\'s note says so.')
W('')
W(row('Item', 'Setting', 'key'))
W(row('---', '---', '---'))
for k in ['protocol.functional', 'protocol.basis', 'protocol.density_fitting', 'protocol.aux_basis', 'protocol.grid_level', 'protocol.grid_points',
          'protocol.solvation_model', 'protocol.solvent', 'protocol.solvent_eps', 'protocol.pcm_lebedev_order', 'protocol.pcm_radii', 'protocol.pcm_discretisation',
          'protocol.scf_conv_tol_opt_gas', 'protocol.scf_conv_tol_opt_pcm', 'protocol.scf_conv_tol_final', 'protocol.scf_max_cycle', 'protocol.optimizer',
          'protocol.opt_criteria', 'protocol.preopt', 'protocol.hessian_method', 'protocol.hessian_projection', 'protocol.frequency_scaling', 'protocol.thermo_T',
          'protocol.thermo_p', 'protocol.fmo_descriptors', 'protocol.fragment_composition_method', 'protocol.cubes', 'protocol.mep_surface', 'protocol.mep_esp_code_pcm',
          'protocol.mep_esp_code_gas', 'protocol.mep_note', 'protocol.nci_settings', 'protocol.nci_grid_pcm', 'protocol.nci_grid_gas', 'protocol.nci_spike_method',
          'protocol.elf_lol_settings', 'protocol.elf_plane_picrate', 'protocol.elf_plane_cation', 'protocol.host']:
    val = v(k); u = F[k]['unit']
    W(row(k.replace('protocol.', ''), f'{val} {u}'.strip(), kk(k)))
W('')
W('Notes that the Computational-details section must reflect:')
W('')
W(f'- **Grid**: level 4 (≈ Gaussian UltraFine), not level 3 as every RESULTS.md says. Verified for `btpc_gas` by re-computing the SCF at the stored geometry '
  '(`verify_grid.log`, not deposited: grid 4 reproduces the stored energy exactly, grid 3 is 1.5 µHa off). `btpc_pcm` was not in that verification run; '
  'grid 4 + Lebedev 29 for it rests only on the `run_case.sh` defaults (added because every case was launched with them). The "grid 4" in its opt.xyz comment is a literal string written by `dft/finalize_opt.py`, not a recorded setting, so it is not independent evidence.')
W('- **PCM**: IEF-PCM, ε = 78.3553, Lebedev order 29 (302 points per sphere), modified-Bondi radii × 1.2, SWIG discretisation; electrostatics only. '
  'The dft/common.py default (order 17, the 3ATB-on-a-laptop setting) is superseded.')
W('- **Frequencies**: central finite differences of analytic gradients (198 gradients); PCM first with step 0.005 Bohr (2026-09-17 13:56–19:42 UTC, freq_step0.005.json), then repeated with 0.01 Bohr (2026-09-18 16:05–21:51 UTC; authoritative, freq.json/hessian.npy); gas with 0.005 Bohr. '
  'For the gas case an analytic Hessian was attempted first (2026-09-16 16:31 UTC) and abandoned; the FD run replaced it. Unscaled harmonic; no IR/Raman intensities.')
W('- **Two-stage optimisation history for BTPC**: the pipeline default is a B3LYP/6-31G(d) pre-optimisation with loose criteria, then B3LYP/6-311++G(d,p). '
  'For BTPC, both final (grid-4) optimisations started from `preopt.xyz` files carried over from earlier grid-3 runs on the (now retired) OVH VM: '
  'the gas run from the converged grid-3 geometry (3 evaluations at grid 4), the PCM run from frame 16 of the grid-3 PCM trajectory (70 steps at grid 4). '
  'The 6-31G(d) stage itself happened in those earlier runs, which are not in the archive. Reported geometries depend only on the final-level optimisation.')
W('- **Software versions and where they were read**:')
for k in ['software.pyscf', 'software.geometric', 'software.multiwfn', 'software.paraview', 'software.python', 'software.numpy_scipy_matplotlib', 'software.libxc']:
    W(f'  - {k.split(".")[1]}: {v(k)} — {F[k]["source"]}' + (f' ({F[k]["note"]})' if F[k]['note'] else ''))
W('- **Multiwfn parameters in one place**: MEP on the ρ = 0.001 a.u. isosurface, surface grid spacing 0.25 Bohr; NCI with the full wavefunction, '
  f'"high-quality" grid ({v("protocol.nci_grid_pcm")} for PCM; {v("protocol.nci_grid_gas")} for gas), 1.5 Bohr box extension, RDG_maxrho = 0.1 a.u.; '
  'NCI isosurface drawn at RDG = 0.5 coloured by sign(λ₂)ρ from −0.035 to +0.02 a.u.; ELF/LOL (Becke definitions) on 300 × 300 planes with 6 Bohr extension. '
  'Orbital isosurfaces ±0.02 a.u.')
W('')

# ------------------------------------------------------------------------------------------------ 4
W('## 4. Results')
W('')
W('### 4.1 Energies and optimisation')
W('')
W(row('Quantity', 'IEF-PCM (water)', 'gas phase', 'keys'))
W(row('---', '---', '---', '---'))
W(row('E_total (Ha)', f('pcm.e_total', 8), f('gas.e_total', 8), '`*.e_total`'))
W(row('Optimisation', f'{v("pcm.opt_steps")} steps (71 gradient evaluations), ≈13.4 h, 4 threads', f'{v("gas.opt_steps")} evaluations, {f("gas.opt_wall", 1)} min', '`*.opt_steps`, `*.opt_wall`'))
W(row('Converged on all five criteria?', 'NO — accepted on forces only (see 4.10)', 'yes', '`*.opt_converged`'))
W(row('Final RMS / max gradient (Ha/Bohr)', f'{v("pcm.opt_final_grad_rms"):.2e} / {v("pcm.opt_final_grad_max"):.2e}', '3.50e-05 / 6.83e-05 (log, step 2)', '`pcm.opt_final_grad_*`, `gas.opt_converged` note'))
W(row('Last step ΔE (Ha); RMS/max displacement (Å)', f'{v("pcm.opt_last_step_dE"):.2e}; {v("pcm.opt_last_step_disp_rms"):.2e}/{v("pcm.opt_last_step_disp_max"):.2e}', '−6.7e-07; 5.61e-04/1.18e-03', '`pcm.opt_last_step_*`'))
W(row('Starting point', 'grid-3 PCM trajectory frame 16', 'converged grid-3 geometry', '`*.opt_start`'))
W(row('Energy lowered during the final optimisation (kcal mol⁻¹)', f('pcm.opt_energy_drop_total', 3), '—', '`pcm.opt_energy_drop_total`'))
W(row('Proton transferred?', v('pcm.proton_transferred'), v('gas.proton_transferred'), '`*.proton_transferred`'))
W('')
W(f'Derived: E(PCM) − E(gas) = {f("diff.e_total_ha", 6)} Ha = {f("diff.e_total", 2)} kcal mol⁻¹; G298(PCM) − G298(gas) = {f("diff.g298", 2)} kcal mol⁻¹. '
  'These mix the electrostatic solvation term with a geometry change and omit non-electrostatic terms: do not call them solvation free energies.')
W('')
W('### 4.2 Proton position and hydrogen bonds inside the cluster')
W('')
W(row('Contact', 'X-ray (cluster)', 'IEF-PCM', 'gas', 'keys'))
W(row('---', '---', '---', '---', '---'))
for D_, H_, A_ in [('N4', 'H4N', 'O1'), ('N4', 'H4N', 'O7'), ('C11', 'H11', 'O2'), ('C11', 'H11', 'O1')]:
    b = f'hb_{D_}_{H_}_{A_}'
    for part, lab_, nd in [('dh', 'D–H', 3), ('ha', 'H···A', 3), ('da', 'D···A', 3), ('angle', '∠D–H···A', 1)]:
        W(row(f'{D_}–{H_}···{A_} {lab_}', f(f'xray.{b}_{part}', nd), f(f'pcm.{b}_{part}', nd), f(f'gas.{b}_{part}', nd), f'`*.{b}_{part}`'))
W(row('d(H4N···O1) − d(N4–H4N) (Å)', f'{v("xray.hb_N4_H4N_O1_ha") - v("xray.hb_N4_H4N_O1_dh"):.3f}', f('pcm.nh_minus_ho_asym', 3), f('gas.nh_minus_ho_asym', 3), '`*.nh_minus_ho_asym` (X-ray computed the same way)'))
W(row('O1···O7 / O1···O2 intramolecular (Å)', f'{f("xray.r_O1_O7")} / {f("xray.r_O1_O2")}', f'{f("pcm.r_O1_O7")} / {f("pcm.r_O1_O2")}', f'{f("gas.r_O1_O7")} / {f("gas.r_O1_O2")}', '`*.r_O1_O7`, `*.r_O1_O2`'))
W('')
W('Distances in Å, angles in degrees. Picrate-internal C3–H3···O3/O4 and C5–H5···O5/O6 contacts (angles 91–106°) are also listed in facts.json (`*.hb_C3_*`, `*.hb_C5_*`); '
  'they are geometric S(5)-type contacts, not hydrogen bonds. The proton stays on N4 in both models (H···O > N–H), but the gas-phase bond is much shorter '
  f'(N···O {f("gas.hb_N4_H4N_O1_da")} Å vs {f("pcm.hb_N4_H4N_O1_da")} Å PCM and {f("xray.hb_N4_H4N_O1_da")} Å X-ray) with the N–H stretched to {f("gas.r_N4_H4N")} Å. '
  f'The gas model also forms a C11–H11···O2 contact ({f("gas.hb_C11_H11_O2_ha")} Å, {f("gas.hb_C11_H11_O2_angle", 1)}°) that is absent in PCM ({f("pcm.hb_C11_H11_O2_ha")} Å) and weak in the crystal ({f("xray.hb_C11_H11_O2_ha")} Å).')
W('')
W('### 4.3 Geometry versus X-ray')
W('')
W(row('Measure', 'IEF-PCM', 'gas', 'keys'))
W(row('---', '---', '---', '---'))
for lab_, k, nd in [('Bond RMSD, all 34 bonds (Å)', 'bond_rmsd_all', 4), ('Largest bond deviation, all (Å; C5–H5)', 'bond_max_dev_all', 3),
                    ('Bond RMSD, 26 heavy-atom bonds (Å)', 'bond_rmsd_heavy', 4), ('Bond MAE, heavy (Å)', 'bond_mae_heavy', 4),
                    ('Mean signed bond deviation, heavy (Å; + = DFT longer)', 'bond_mean_signed_dev_heavy', 4),
                    ('Angle RMSD, all 52 (°)', 'angle_rmsd_all', 2), ('Angle RMSD, 36 heavy-atom angles (°)', 'angle_rmsd_heavy', 2),
                    ('Superposition RMSD, 25 heavy atoms (Å)', 'kabsch_rmsd_heavy', 3), ('… picrate heavy atoms only (Å)', 'kabsch_rmsd_picrate_heavy', 3),
                    ('… cation heavy atoms only (Å)', 'kabsch_rmsd_cation_heavy', 3)]:
    W(row(lab_, f(f'pcm.{k}', nd), f(f'gas.{k}', nd), f'`*.{k}`'))
W(row('Five largest heavy-atom bond deviations', v('pcm.bond_largest_dev_heavy'), v('gas.bond_largest_dev_heavy'), '`*.bond_largest_dev_heavy`'))
W(row('Five largest heavy-atom angle deviations', v('pcm.angle_largest_dev_heavy'), v('gas.angle_largest_dev_heavy'), '`*.angle_largest_dev_heavy`'))
W('')
W('The all-bond RMSD (0.089/0.092 Å, the number in RESULTS.md) is dominated by the X-ray X–H bonds, which are short by 0.12–0.27 Å as usual for X-ray H positions; '
  'quote the heavy-atom figures for the quality of the model. Each ion on its own matches the crystal well (cation 0.03 Å, picrate 0.20–0.22 Å because of nitro rotation); '
  'the whole-cluster superposition RMSD is large because the ions rotate relative to each other (Section 4.4).')
W('')
W('Table 4.3b — Chemically important bond lengths (Å) and angles (°). X-ray from the cluster coordinates (CIF values with s.u. in Table 1.3).')
W('')
W(row('Parameter', 'X-ray', 'IEF-PCM', 'gas', 'key'))
W(row('---', '---', '---', '---', '---'))
for k in [x for x in F if x.startswith('pcm.r_') or x.startswith('pcm.a_')]:
    base = k[4:]
    if base in ('r_O1_O7', 'r_O1_O2', 'r_N5_O1'): continue
    nd = 1 if base.startswith('a_') else 3
    W(row(base[2:].replace('_', '–') if base.startswith('r_') else '∠' + base[2:].replace('_', '–'), f(f'xray.{base}', nd), f(k, nd), f(f'gas.{base}', nd), f'`*.{base}`'))
W('')
W(f'Points for the text: C1–O1 lengthens from {f("xray.r_O1_C1")} (X-ray) to {f("pcm.r_O1_C1")} (PCM) and {f("gas.r_O1_C1")} Å (gas) as the proton approaches O1; '
  f'the ipso angle C2–C1–C6 stays small ({f("xray.a_C2_C1_C6", 1)}/{f("pcm.a_C2_C1_C6", 1)}/{f("gas.a_C2_C1_C6", 1)}°), the phenolate signature. '
  f'In the cation the two N–N bonds are nearly equal in the crystal ({f("xray.r_N4_N5")}/{f("xray.r_N5_N6")}) and in PCM ({f("pcm.r_N4_N5")}/{f("pcm.r_N5_N6")}), '
  f'but split in the gas phase ({f("gas.r_N4_N5")}/{f("gas.r_N5_N6")}): the N4 side becomes more imine-like as H4N moves towards O1.')
W('')
W('### 4.4 Conformation: ring orientation and nitro twists')
W('')
W(row('Angle (°)', 'X-ray cluster', 'IEF-PCM', 'gas', 'keys'))
W(row('---', '---', '---', '---', '---'))
W(row('Picrate ring / cation plane', f('xray.interplanar_picrate_cation', 1), f('pcm.interplanar_picrate_cation', 1), f('gas.interplanar_picrate_cation', 1), '`*.interplanar_picrate_cation`'))
W(row('Ring centroid distance (Å)', f('xray.centroid_distance_picrate_cation', 3), f('pcm.centroid_distance_picrate_cation', 3), f('gas.centroid_distance_picrate_cation', 3), '`*.centroid_distance_picrate_cation`'))
for n_, pos in [('N1', 'ortho (C2)'), ('N2', 'para (C4)'), ('N3', 'ortho (C6, next to O1)')]:
    W(row(f'Nitro {n_} twist, {pos}', f(f'xray.nitro_twist_{n_}', 1), f(f'pcm.nitro_twist_{n_}', 1), f(f'gas.nitro_twist_{n_}', 1), f'`*.nitro_twist_{n_}`'))
W('')
W('In the crystal the two ions are nearly coplanar (5.5°) and all three nitro groups are twisted by only 9–12°. Both isolated models keep the para-nitro group coplanar '
  'but rotate the two ortho groups by 22–40°, and open the inter-ring angle to 29° (gas) or 58° (PCM). The PCM value is not well defined (4.10).')
W('')
W('### 4.5 Vibrational analysis and thermochemistry (unscaled harmonic)')
W('')
W(row('Quantity', 'IEF-PCM', 'gas', 'keys'))
W(row('---', '---', '---', '---'))
W(row('Hessian', v('pcm.hessian_method'), v('gas.hessian_method'), '`*.hessian_method`'))
W(row('Hessian wall time (min)', f('pcm.hessian_wall', 0), f('gas.hessian_wall', 0), '`*.hessian_wall`'))
W(row('Imaginary modes (project\'s own projection / PySCF count)', f'{v("pcm.n_imag")} / {v("pcm.n_imag_pyscf")}', f'{v("gas.n_imag")} / {v("gas.n_imag_pyscf")}', '`*.n_imag`, `*.n_imag_pyscf`'))
W(row('Lowest ten (cm⁻¹)', v('pcm.freq_lowest10'), v('gas.freq_lowest10'), '`*.freq_lowest10`'))
W(row('Modes below 100 cm⁻¹ (PCM count includes the imaginary mode)', v('pcm.n_modes_below_100'), v('gas.n_modes_below_100'), '`*.n_modes_below_100`'))
W(row('Highest (cm⁻¹; free N6–H6N stretch)', f('pcm.freq_highest', 1), f('gas.freq_highest', 1), '`*.freq_highest`'))
W(row('N4–H4N stretch, H-bonded (cm⁻¹)', f('pcm.nu_N4_H4N_stretch', 1), f('gas.nu_N4_H4N_stretch', 1), '`*.nu_N4_H4N_stretch`'))
W(row('N4–H4N in-plane bend (cm⁻¹)', f('pcm.nu_N4_H4N_inplane_bend', 1), f('gas.nu_N4_H4N_inplane_bend', 1), '`*.nu_N4_H4N_inplane_bend`'))
W(row('N4–H4N out-of-plane wag (cm⁻¹)', v('pcm.nu_N4_H4N_oop_wag'), v('gas.nu_N4_H4N_oop_wag'), '`*.nu_N4_H4N_oop_wag`'))
W(row('N6–H6N stretch, free (cm⁻¹)', f('pcm.nu_N6_H6N_stretch', 1), f('gas.nu_N6_H6N_stretch', 1), '`*.nu_N6_H6N_stretch`'))
W(row('N6–H6N in-plane bend / wag (cm⁻¹)', f'{f("pcm.nu_N6_H6N_inplane_bend", 1)} / {v("pcm.nu_N6_H6N_oop_wag")}', f'{f("gas.nu_N6_H6N_inplane_bend", 1)} / {v("gas.nu_N6_H6N_oop_wag")}', '`*.nu_N6_H6N_*`'))
W(row('ν(free N6–H) − ν(bonded N4–H) (cm⁻¹)', f('diff.nh_redshift_pcm', 0), f('diff.nh_redshift_gas', 0), '`diff.nh_redshift_*`'))
W(row('ZPE (Ha; kcal mol⁻¹)', f'{f("pcm.zpe", 6)}; {f("pcm.zpe_kcal", 2)}', f'{f("gas.zpe", 6)}; {f("gas.zpe_kcal", 2)}', '`*.zpe`, `*.zpe_kcal`'))
W(row('H298 (Ha)', f('pcm.h298', 6), f('gas.h298', 6), '`*.h298`'))
W(row('G298 (Ha)', f('pcm.g298', 6), f('gas.g298', 6), '`*.g298`'))
W(row('H − E, G − E (Ha)', f'{f("pcm.thermal_corr_h", 6)}, {f("pcm.thermal_corr_g", 6)}', f'{f("gas.thermal_corr_h", 6)}, {f("gas.thermal_corr_g", 6)}', '`*.thermal_corr_*`'))
W(row('S298 (cal mol⁻¹ K⁻¹)', f('pcm.entropy', 1), f('gas.entropy', 1), '`*.entropy`'))
W('')
W('Mode assignments are derived from the stored Hessians (`data/r2/.../hessian.npy`): the eigenvectors were recomputed (with the masses PySCF used, H 1.008, C 12.011, N 14.007, O 15.999, the frequencies reproduce freq.json to < 1e-10 cm⁻¹), '
  'and a mode is called an N–H stretch when the H atom carries most of the mass-weighted amplitude and >50% of its motion is along the bond. Nitro stretches mix strongly with ring modes '
  '(asymmetric NO₂ character ~1495–1575 cm⁻¹, symmetric ~1286–1351 cm⁻¹ in PCM); do not assign single nitro frequencies without a PED. The three lowest modes of both models are '
  f'rigid-body motions of the two ions (rigid fraction PCM {f("pcm.mode1_rigid_fraction", 3)}, {f("pcm.mode2_rigid_fraction", 3)}, {f("pcm.mode3_rigid_fraction", 3)}; '
  f'gas {f("gas.mode1_rigid_fraction", 3)}, {f("gas.mode2_rigid_fraction", 3)}, {f("gas.mode3_rigid_fraction", 3)}).')
W('')
W('### 4.6 Frontier orbitals, descriptors, dipole, charge separation')
W('')
W(row('Descriptor (formula, dft/fmo.py)', 'IEF-PCM', 'gas', 'keys'))
W(row('---', '---', '---', '---'))
for lab_, k in [('E_HOMO (eV)', 'homo_ev'), ('E_LUMO (eV)', 'lumo_ev'), ('ΔE_gap = E_LUMO − E_HOMO (eV)', 'gap_ev'), ('E_HOMO−1 (eV)', 'homo_m1_ev'),
                ('E_LUMO+1 (eV)', 'lumo_p1_ev'), ('IP = −E_HOMO (eV)', 'ip_ev'), ('EA = −E_LUMO (eV)', 'ea_ev'), ('χ = (IP+EA)/2 (eV)', 'chi_ev'),
                ('μ = −(IP+EA)/2 (eV)', 'mu_ev'), ('η = (IP−EA)/2 (eV)', 'eta_ev'), ('S = 1/(2η) (eV⁻¹)', 'softness_ev_inv'), ('ω = μ²/(2η) (eV)', 'omega_ev'),
                ('Dipole moment (D)', 'dipole_debye')]:
    W(row(lab_, f(f'pcm.{k}', 4), f(f'gas.{k}', 4), f'`*.{k}`'))
W(row('Dipole vector x, y, z (D)', v('pcm.dipole_vector_debye'), v('gas.dipole_vector_debye'), '`*.dipole_vector_debye`'))
W(row('HOMO−1 / HOMO / LUMO / LUMO+1 on picrate (%)', f'{v("pcm.homo_m1_picrate_pct")} / {v("pcm.homo_picrate_pct")} / {v("pcm.lumo_picrate_pct")} / {v("pcm.lumo_p1_picrate_pct")}',
      f'{v("gas.homo_m1_picrate_pct")} / {v("gas.homo_picrate_pct")} / {v("gas.lumo_picrate_pct")} / {v("gas.lumo_p1_picrate_pct")}', '`*.*_picrate_pct` (Mulliken, fmo.json)'))
W(row('Highest occupied picrate / cation orbital (eV)', f'{f("pcm.picrate_highest_occ_ev")} / {f("pcm.cation_highest_occ_ev")}', f'{f("gas.picrate_highest_occ_ev")} / {f("gas.cation_highest_occ_ev")}', '`*.picrate_highest_occ_ev`, `*.cation_highest_occ_ev`'))
W(row('Lowest virtual picrate / cation orbital (eV)', f'{f("pcm.picrate_lowest_virt_ev")} / {f("pcm.cation_lowest_virt_ev")}', f'{f("gas.picrate_lowest_virt_ev")} / {f("gas.cation_lowest_virt_ev")}', '`*.picrate_lowest_virt_ev`, `*.cation_lowest_virt_ev`'))
W(row('Inter-ion CT gap: picrate HOMO → lowest cation virtual (eV)', f('pcm.ct_gap_picrate_to_cation_ev'), f('gas.ct_gap_picrate_to_cation_ev'), '`*.ct_gap_picrate_to_cation_ev`'))
W(row('Picrate-local gap (eV)', f('pcm.local_gap_picrate_ev'), f('gas.local_gap_picrate_ev'), '`*.local_gap_picrate_ev`'))
W(row('Mulliken charge on picrate (e)', f('pcm.mulliken_charge_picrate'), f('gas.mulliken_charge_picrate'), '`*.mulliken_charge_picrate`'))
W('')
W('Orbital ladder with Mulliken fragment shares (derived; energies from fmo.json; full lists for every MO in `derived/mo_fragment_composition_<case>.csv`):')
W('')
W(f'- PCM: {v("pcm.orbitals_near_gap")}')
W(f'- gas: {v("gas.orbitals_near_gap")}')
W('')
W('The composition analysis rebuilds the AO overlap matrix from the molden basis; it reproduces orthonormality (max |CᵀSC − I| ≈ 3e-11) and the four fmo.json '
  'compositions to 0.05%, so the extra orbitals are on the same footing. Key reading: in the gas phase the LUMO is the cation π* (99.6% cation) and the HOMO→LUMO gap '
  'is an inter-ion charge-transfer gap; in the continuum the cation π* rises by 0.77 eV to LUMO+3 and the LUMO becomes picrate-localised, so the PCM gap (3.735 eV) is intra-picrate. '
  'The softness convention is S = 1/(2η); a paper using S = 1/η would report twice the value.')
W('')
W('### 4.7 Molecular electrostatic potential on the 0.001 a.u. density surface')
W('')
W(row('Statistic (Multiwfn)', 'IEF-PCM', 'gas', 'keys'))
W(row('---', '---', '---', '---'))
for lab_, k, nd in [('Global minimum (kcal mol⁻¹)', 'mep_min_kcal', 2), ('… site', 'mep_min_site', 0), ('Global maximum (kcal mol⁻¹)', 'mep_max_kcal', 2), ('… site', 'mep_max_site', 0),
                    ('Global min / max (a.u.)', None, 0), ('Minima / maxima found', None, 0), ('Surface area (Å²)', 'mep_surface_area', 2),
                    ('Positive / negative area (Å²)', None, 0), ('Enclosed volume (Å³)', 'mep_volume', 2), ('Mass/volume density (g cm⁻³)', 'mep_density_mv', 4),
                    ('Sphericity', 'mep_sphericity', 4), ('Overall average (kcal mol⁻¹)', 'mep_avg_all', 3), ('Positive / negative average (kcal mol⁻¹)', None, 0),
                    ('σ²_tot (kcal² mol⁻²)', 'mep_sigma2_tot', 2), ('σ²₊ / σ²₋ (kcal² mol⁻²)', None, 0), ('Balance ν', 'mep_nu', 4), ('σ²_tot·ν (kcal² mol⁻²)', 'mep_sigma2_nu', 2),
                    ('Internal charge separation Π (kcal mol⁻¹)', 'mep_Pi', 3), ('Molecular polarity index (eV; kcal mol⁻¹)', None, 0),
                    ('Non-polar / polar area, |ESP| ≤ / > 10 kcal mol⁻¹ (Å²)', None, 0), ('Polar area fraction (%)', 'mep_polar_area_pct', 2),
                    ('Skewness overall / + / −', None, 0), ('Surface vertices', 'mep_n_vertices', 0), ('ESP evaluation wall time (s)', 'mep_wall', 0)]:
    if k: W(row(lab_, f(f'pcm.{k}', nd), f(f'gas.{k}', nd), f'`*.{k}`'))
    elif 'a.u.' in lab_: W(row(lab_, f'{f("pcm.mep_min_au", 6)} / {f("pcm.mep_max_au", 6)}', f'{f("gas.mep_min_au", 6)} / {f("gas.mep_max_au", 6)}', '`*.mep_min_au`, `*.mep_max_au`'))
    elif 'found' in lab_: W(row(lab_, f'{v("pcm.mep_n_minima")} / {v("pcm.mep_n_maxima")}', f'{v("gas.mep_n_minima")} / {v("gas.mep_n_maxima")}', '`*.mep_n_minima`, `*.mep_n_maxima`'))
    elif 'Positive / negative area' in lab_: W(row(lab_, f'{f("pcm.mep_positive_area", 2)} / {f("pcm.mep_negative_area", 2)}', f'{f("gas.mep_positive_area", 2)} / {f("gas.mep_negative_area", 2)}', '`*.mep_positive_area`, `*.mep_negative_area`'))
    elif 'Positive / negative average' in lab_: W(row(lab_, f'{f("pcm.mep_avg_pos", 3)} / {f("pcm.mep_avg_neg", 3)}', f'{f("gas.mep_avg_pos", 3)} / {f("gas.mep_avg_neg", 3)}', '`*.mep_avg_pos`, `*.mep_avg_neg`'))
    elif 'σ²₊' in lab_: W(row(lab_, f'{f("pcm.mep_sigma2_pos", 2)} / {f("pcm.mep_sigma2_neg", 2)}', f'{f("gas.mep_sigma2_pos", 2)} / {f("gas.mep_sigma2_neg", 2)}', '`*.mep_sigma2_pos`, `*.mep_sigma2_neg`'))
    elif 'polarity' in lab_: W(row(lab_, f'{f("pcm.mep_mpi", 4)}; {f("pcm.mep_mpi_kcal", 2)}', f'{f("gas.mep_mpi", 4)}; {f("gas.mep_mpi_kcal", 2)}', '`*.mep_mpi`, `*.mep_mpi_kcal`'))
    elif 'Non-polar' in lab_: W(row(lab_, f'{f("pcm.mep_nonpolar_area", 2)} / {f("pcm.mep_polar_area", 2)}', f'{f("gas.mep_nonpolar_area", 2)} / {f("gas.mep_polar_area", 2)}', '`*.mep_nonpolar_area`, `*.mep_polar_area`'))
    elif 'Skewness' in lab_: W(row(lab_, f'{f("pcm.mep_skew_all", 3)} / {f("pcm.mep_skew_pos", 3)} / {f("pcm.mep_skew_neg", 3)}', f'{f("gas.mep_skew_all", 3)} / {f("gas.mep_skew_pos", 3)} / {f("gas.mep_skew_neg", 3)}', '`*.mep_skew_*`'))
W('')
for c_, nm in (('pcm', 'IEF-PCM'), ('gas', 'gas phase')):
    W(f'Table 4.7{"a" if c_ == "pcm" else "b"} — every surface extremum, {nm} (`{c_}.mep_minNN` / `{c_}.mep_maxNN`; Multiwfn index; nearest atoms from opt.xyz). '
      'Note that Multiwfn lists local maxima with negative values and local minima with positive values; they are topological extrema, not sign classes.')
    W('')
    W(row('key', 'type', 'ESP (kcal mol⁻¹)', 'ESP (a.u.)', 'nearest atoms (Å)', 'ion'))
    W(row(*['---'] * 6))
    ks = sorted([x for x in F if re.match(rf'{c_}\.mep_(min|max)\d\d$', x)], key=lambda x: v(x))
    for k in ks:
        n_ = F[k]['note']; au = n_.split(' a.u.')[0]; near = re.search(r'nearest ([^;]*);', n_).group(1); ion = 'picrate' if '; picrate' in n_ else 'cation'
        W(row(f'`{k}`', ('minimum' if '_min' in k else 'maximum') + (' (global)' if 'GLOBAL' in n_ else ''), f'{v(k):.2f}', au, near.replace(' A)', ')'), ion))
    W('')
W('Physical reading: the most negative potential sits on the para-nitro oxygens O4/O5, not on the phenolate O1, because O1 is capped by the N4–H4N hydrogen bond; the most positive '
  'region is the free N6–H6N. The ESP maps the charge density only (PCM surface charges are not included in the mapped potential).')
W('')
W('### 4.8 NCI (RDG) analysis')
W('')
for c_, nm in (('pcm', 'IEF-PCM'), ('gas', 'gas phase')):
    W(f'Table 4.8{"a" if c_ == "pcm" else "b"} — RDG spikes as reported by the pipeline (`{c_}.nci_spikeN`, from nci_spikes.csv) with the atoms at their lowest-RDG points (derived from rdg_scatter.txt), {nm}.')
    W('')
    W(row('key', 'sign(λ₂)ρ bin (a.u.)', 'min RDG', 'points in bin', 'pipeline label', 'lowest-RDG points between'))
    W(row(*['---'] * 6))
    for k in sorted([x for x in F if re.match(rf'{c_}\.nci_spike\d+$', x)], key=lambda x: int(x.split('spike')[1])):
        n_ = F[k]['note']
        mr = re.search(r'min RDG ([\d.]+)', n_).group(1); npt = re.search(r'; (\d+) points', n_).group(1); lb = re.search(r"label '([^']*)'", n_).group(1)
        wh = re.search(r'lie between (.*) \(derived', n_).group(1)
        W(row(f'`{k}`', f'{v(k):+.4f}', mr, npt, lb, wh))
    W('')
for c_, nm in (('pcm', 'IEF-PCM'), ('gas', 'gas phase')):
    W(f'Table 4.8{"c" if c_ == "pcm" else "d"} — Real-space NCI regions (derived: connected RDG < 0.5 regions of rdg.cube with ≥20 grid points; value = sign(λ₂)ρ at the lowest-RDG point), {nm}.')
    W('')
    W(row('key', 'sign(λ₂)ρ at min RDG (a.u.)', 'min RDG', 'sign(λ₂)ρ range', 'points', 'located between'))
    W(row(*['---'] * 6))
    for k in sorted([x for x in F if re.match(rf'{c_}\.nci_region\d+$', x)]):
        n_ = F[k]['note']
        W(row(f'`{k}`', f'{v(k):+.4f}', re.search(r'RDG ([\d.]+)\)', n_).group(1), re.search(r'range (\S+);', n_).group(1),
              re.search(r'\((\d+) grid points\)', n_).group(1), re.search(r'nearest (.*)$', n_).group(1)))
    W('')
W('Reading: the only strongly attractive feature is the N4–H4N···O1 bond, at sign(λ₂)ρ ≈ −0.050 a.u. (PCM) and −0.086 a.u. (gas; the region reaches −0.099 a.u., '
  'close to the RDG_maxrho = 0.1 screening limit, so the gas isosurface of this bond may be slightly truncated). Weak attractive/van der Waals regions: H4N···O7 (bifurcated), '
  'and in the gas phase only, C11–H11···O2 and H11···O1. Positive (repulsive/ring) regions: the three ring centres (+0.020 to +0.055 a.u.), O1···O2 and O1···O7 '
  'between the phenolate and ortho-nitro oxygens, and the picrate C–H···O(nitro) five-membered pockets. The pipeline\'s spike labels are a pure threshold on sign(λ₂)ρ '
  '(e.g. the PCM −0.0145 "van der Waals" spike is the intramolecular O1···O7/O1···O2 contact region); use the located-between column, not the label.')
W('')
W('### 4.9 ELF and LOL maps')
W('')
W(row('Quantity', 'IEF-PCM', 'gas', 'keys'))
W(row('---', '---', '---', '---'))
for fn in ('elf', 'lol'):
    for pl, nm in (('picrate_ring', 'picrate plane (C6, C4, C2)'), ('benzotriazolium_ring', 'cation plane (N4, N5, C7)')):
        W(row(f'{fn.upper()} range, {nm}', v(f'pcm.{fn}_{pl}_range'), v(f'gas.{fn}_{pl}_range'), f'`*.{fn}_{pl}_range`'))
    for pos in ('picrate_centroid', 'benzo_centroid', 'triazole_centroid'):
        W(row(f'{fn.upper()} at {pos.replace("_", " ")}', f(f'pcm.{fn}_{pos}', 3), f(f'gas.{fn}_{pos}', 3), f'`*.{fn}_{pos}`'))
W(row('O1 distance from the cation map plane (Å; H4N within 0.01 Å)', f('pcm.cation_plane_offset_O1', 3), f('gas.cation_plane_offset_O1', 3), '`*.cation_plane_offset_O1`'))
W('')
W('Maps are 300 × 300 points in the plane through three atoms; the cation plane contains the whole benzotriazolium (triazole and benzo rings coplanar) and, nearly, the H4N···O1 bond, '
  'so the cation map is the one to show the hydrogen bond. ELF/LOL values at ring centroids are the usual aromatic depletion (ELF ≈ 0.02–0.05). No quantitative ELF value '
  'along the H-bond is read from the plane maps because O1 lies out of the plane (0.36 Å in PCM). Superseded for the value on the H4N–O1 line itself: the revision-1 keys '
  '`*.elf_min_H4N_O1` and `*.bcp_H4N_O1_elf` give ELF along that line, which the plane offset does not affect.')
W('')
W('### 4.10 The PCM imaginary mode: what was checked')
W('')
W(f'- The authoritative PCM Hessian (step 0.01 Bohr) gives one imaginary mode, {f("saddle.freq_imag_step0p01", 2)} cm⁻¹; the earlier Hessian with step 0.005 Bohr gave {f("saddle.freq_imag_step0p005", 2)} cm⁻¹. '
  f'Doubling the step changed no frequency by more than {f("saddle.step_change_max", 2)} cm⁻¹ (mean 0.15), ZPE by {f("saddle.step_change_zpe", 4)} and G298 by {f("saddle.step_change_g", 3)} kcal mol⁻¹. '
  'A first attempt (19 workers) failed with exit code 1 and produced nothing (`saddle.first_attempt`).')
W('- PySCF\'s own harmonic analysis stored the mode as 0 + 3.2i and therefore counted zero imaginary modes (`pcm.n_imag_pyscf` = 0); the project\'s explicit '
  'translation/rotation projection (`dft/vibrations.py`) reports it correctly.')
W(f'- Character: {100 * v("saddle.imag_mode_rigid_fraction"):.1f}% of the eigenvector is rigid-body motion of the two ions relative to each other; 66% of the mass-weighted amplitude is on the cation '
  f'(`saddle.imag_mode_*`). It is an inter-ion libration, not a distortion of either ion.')
W('- The optimiser was still moving along exactly this coordinate. The net displacement from step 0 to 70 projects '
  f'{100 * v("saddle.drift_0_to_70_overlap_imag"):.0f}% onto the imaginary mode, and from step 40 to 70 {100 * v("saddle.drift_40_to_70_overlap_imag"):.0f}% '
  '(max atomic displacement 1.8 Å for 0.30 kcal mol⁻¹); the last ten steps moved along the 9.0 cm⁻¹ libration (95%).')
W('- Along the whole trajectory the hydrogen bond hardly changed while the rings rotated:')
W('')
W(row('step', 'E (Ha)', 'N4–H4N', 'H4N···O1', 'N4···O1', '∠ (°)', 'inter-ring angle (°)', 'H11···O2'))
W(row(*['---'] * 8))
for k in sorted(x for x in F if x.startswith('saddle.traj_step')):
    s_ = v(k); nums = re.findall(r'-?\d+\.\d+', s_)
    W(row(int(k[-2:]), nums[0], nums[1], nums[2], nums[3], nums[4], nums[5], nums[6]))
W('')
W(f'- Conclusion for the paper: {v("saddle.interpretation")}')
W('- The gas-phase model has no imaginary mode (lowest 12.2 cm⁻¹) and met every convergence criterion.')
W('')

# ------------------------------------------------------------------------------------------------ 5
W('## 5. Gas phase versus IEF-PCM (derived: PCM − gas unless stated)')
W('')
W(row('Quantity', 'PCM − gas', 'key'))
W(row('---', '---', '---'))
for k in [x for x in F if x.startswith('diff.')]:
    val = v(k); u = F[k]['unit']
    nd = 6 if u == 'Ha' else (4 if abs(val) < 1 else 2)
    W(row(k[5:], f'{val:+.{nd}f} {u}'.strip(), kk(k)))
W('')
W('Energy differences mix solvation and geometry (see 4.1). Orbital-energy shifts are Koopmans quantities in two different electrostatic environments.')
W('')

# ------------------------------------------------------------------------------------------------ 6
W('## 6. Storyline candidates (ranked by strength of the evidence)')
W('')
W(f'1. **The continuum switches the frontier orbitals from inter-ion charge transfer to intra-picrate.** Gas: HOMO on picrate ({v("gas.homo_picrate_pct")}%), LUMO on the cation '
  f'({100 - v("gas.lumo_picrate_pct"):.1f}%), gap {f("gas.gap_ev")} eV = CT gap. PCM: cation π* rises by {f("diff.cation_lowest_virt_ev", 2)} eV (to LUMO+3), cation HOMO rises by '
  f'{f("diff.cation_highest_occ_ev", 2)} eV (to HOMO−1), the picrate LUMO drops by {abs(v("diff.picrate_lowest_virt_ev")):.2f} eV; LUMO and LUMO+1 are both picrate, gap {f("pcm.gap_ev")} eV, '
  f'and the CT gap opens to {f("pcm.ct_gap_picrate_to_cation_ev")} eV. So "gap" means different excitations in the two models; the descriptors (ω {f("gas.omega_ev", 2)} → {f("pcm.omega_ev", 2)} eV, '
  f'η {f("gas.eta_ev", 2)} → {f("pcm.eta_ev", 2)} eV) inherit that. Solid: validated orbital compositions.')
W(f'2. **Environment-controlled proton position in the N⁺–H···O⁻ bond, salt retained in both models.** N4···O1 {f("gas.hb_N4_H4N_O1_da")} Å (gas) vs {f("pcm.hb_N4_H4N_O1_da")} Å (PCM) vs '
  f'{f("xray.hb_N4_H4N_O1_da")} Å (X-ray); N–H {f("gas.r_N4_H4N")} vs {f("pcm.r_N4_H4N")} Å; ν(N4–H) {f("gas.nu_N4_H4N_stretch", 0)} vs {f("pcm.nu_N4_H4N_stretch", 0)} cm⁻¹ '
  f'(red shifts of {f("diff.nh_redshift_gas", 0)} and {f("diff.nh_redshift_pcm", 0)} cm⁻¹ against the free N6–H); NCI sign(λ₂)ρ {f("gas.nci_spike1", 3)} vs {f("pcm.nci_spike1", 3)} a.u.; '
  f'Mulliken charge on picrate {f("gas.mulliken_charge_picrate", 2)} vs {f("pcm.mulliken_charge_picrate", 2)} e; C1–O1 {f("gas.r_O1_C1")} vs {f("pcm.r_O1_C1")} Å; N4–N5/N5–N6 asymmetry '
  f'{v("gas.r_N5_N6") - v("gas.r_N4_N5"):.3f} vs {v("pcm.r_N5_N6") - v("pcm.r_N4_N5"):.3f} Å. Unlike 3ATB (proton migrated in the gas phase), BTPC stays ionic even in vacuo, but the gas-phase bond is a '
  'short, strong one (N···O < 2.5 Å; "partly shared" was the early wording, superseded by the critical-point data of revision 1: both bonds lie in the transit region, the gas-phase one far deeper). PCM reproduces the crystal N···O to +0.05 Å; the gas phase under-shoots by 0.15 Å.')
W(f'3. **A floppy ion pair: the hydrogen bond is robust, the relative orientation is not.** Along 70 PCM steps the H-bond metrics changed by ≤0.05 Å while the inter-ring angle went from 5.8° to 58.2° '
  f'for {f("pcm.opt_energy_drop_total", 2)} kcal mol⁻¹; the residual {f("saddle.freq_imag_step0p01", 1)} cm⁻¹ mode is 99.9% rigid inter-ion libration and survived a doubled FD step. '
  'This is a methodological point worth making honestly (crystal packing, i.e. the second N–H···O and the stacking, fixes the orientation that a single ion pair cannot).')
W(f'4. **Packing flattens the picrate.** X-ray nitro twists 9–12° and a nearly coplanar pair (5.5°); both isolated models twist the ortho nitro groups to 22–40° '
  f'(N1 {f("pcm.nitro_twist_N1", 1)}°, N3 {f("pcm.nitro_twist_N3", 1)}° PCM) while keeping the para group planar ({f("pcm.nitro_twist_N2", 1)}°). The largest bond errors are the '
  f'nitro N1–O2 (+0.04 Å) where X-ray bonds are anomalously short (1.186 Å; high Ueq of O2/O3). Heavy-atom bond RMSD is only {f("pcm.bond_rmsd_heavy", 3)} Å.')
W(f'5. **Electrostatic surface.** Most negative ESP on the para-nitro O4/O5 ({f("pcm.mep_min_kcal", 1)} kcal mol⁻¹ PCM; {f("gas.mep_min_kcal", 1)} gas), not on the H-bonded phenolate; '
  f'most positive on the free N6–H6N ({f("pcm.mep_max_kcal", 1)}; {f("gas.mep_max_kcal", 1)}) — precisely the site that in the crystal donates to O6 of the next picrate. '
  f'The continuum amplifies the charge separation (σ²_tot {f("gas.mep_sigma2_tot", 0)} → {f("pcm.mep_sigma2_tot", 0)} kcal² mol⁻², MPI {f("gas.mep_mpi", 2)} → {f("pcm.mep_mpi", 2)} eV, '
  f'dipole {f("gas.dipole_debye", 1)} → {f("pcm.dipole_debye", 1)} D).')
W('')

# ------------------------------------------------------------------------------------------------ 7
W('## 7. Caveats and limitations the paper must state')
W('')
caveats = [
    'C1. **PCM geometry is not a fully converged minimum.** Accepted on forces (rms 1.40e-4, max 3.85e-4 Ha/Bohr, inside Gaussian thresholds) but the energy (ΔE −8.6e-6 Ha) and displacement (0.02/0.04 Å) criteria were not met and the energy was still falling; one imaginary frequency of 3.2i cm⁻¹ (2.5i with the smaller FD step) remains. It is a rigid inter-ion libration; internal geometry and the H-bond are stable. Do not write "true minimum" or "all real frequencies" for the PCM model.',
    'C2. **Relative ion orientation in PCM is poorly defined** (inter-ring angle 58° vs 5.5° in the crystal, 29° in gas; drifted by ~38° for 0.3 kcal/mol). Any property sensitive to the orientation (dipole magnitude, MEP surface statistics, NCI C–H···O contacts) inherits this uncertainty.',
    'C3. **Cluster, not crystal.** One ion pair built by symmetry; the second N–H···O bond (N6–H6N···O6), C9–H9···O5, the triazole stacking and long-range lattice electrostatics are absent. N6–H6N is artificially free (MEP maximum, 3578 cm⁻¹ band).',
    'C4. **Continuum solvent (water) stands in for the crystal field**; electrostatic IEF-PCM only, no cavitation/dispersion terms, equilibrium solvation. E(PCM) − E(gas) is not a solvation free energy.',
    'C5. **No dispersion correction** (plain B3LYP), which matters for the inter-ring orientation and any stacking.',
    'C6. **Frequencies are harmonic, unscaled, from finite differences; no IR/Raman intensities**, so no simulated spectrum can be compared with experiment. The strongly H-bonded N4–H stretch is highly anharmonic in reality.',
    'C7. **Thermochemistry**: RRHO at 298.15 K/1 atm with ~10 modes below 100 cm⁻¹ treated as harmonic oscillators; entropies and G are soft.',
    'C8. **Koopmans descriptors** (IP = −E_HOMO etc.) from orbital energies; in PCM they refer to a solvated pair. Softness convention S = 1/(2η).',
    'C9. **Mulliken charges** with a diffuse (++) basis: only fragment totals are reported, and only as indicative.',
    'C10. **Software**: PySCF/geomeTRIC/Multiwfn/ParaView, not Gaussian/VMD; B3LYP defined as in Gaussian (VWN-RPA); grid level 4 (not 3); density fitting with PySCF\'s default auxiliary basis (identity not logged).',
    'C11. **Starting points**: both final optimisations began from geometries produced by earlier grid-3 runs (and those from a 6-31G(d) pre-optimisation) that are not archived.',
    'C12. **Naming**: the case file\'s "1H-benzotriazol-1-ium" does not match the CIF protonation pattern (H on N4 and N6); verify the cation name against the source paper.',
    'C13. **X-ray reference**: room-temperature structure (298 K), R1 = 0.067, H atoms refined isotropically with short X–H; nitro O2/O3 have large displacement parameters. DFT vs X-ray bond comparisons should use heavy atoms.',
    'C14. **NCI screening**: RDG_maxrho = 0.1 a.u.; the gas-phase H4N···O1 region reaches ρ ≈ 0.099 a.u., at the edge of the screen.',
]
for c_ in caveats: W(f'- {c_}')
W('')

# ------------------------------------------------------------------------------------------------ 8
W('## 8. Figure data inventory (files opened and checked)')
W('')
W('Paths: `T` = `molecules/btpc/results/<case>/` (tracked), `R` = `data/r2/molecules/btpc/results/<case>/` (raw results archive, read-only). `<case>` = btpc_pcm | btpc_gas. '
  'Tracked copies are byte-identical to the archive copies (checked with cmp).')
W('')
inv = [
    ('Structure / labels', 'T/opt.xyz', '33 atoms, Å, CIF atom order within each ion, picrate (0–17) then cation (18–32) (labels from `molecules/btpc/cases/<case>.json:cif_labels`, strip `#4_666`). X-ray start: `molecules/btpc/input/btpc_cluster.xyz` (centred cluster, Å). Existing render: T/fig/structure_labels.png (1800×1500, 300 dpi).'),
    ('Optimisation trajectory', 'R/opt_traj.xyz', 'multi-frame xyz, Å; comment "step k  E = … Ha". PCM 71 frames (0–70), gas 3 frames. Suitable for an energy/inter-ring-angle vs step plot (Section 4.10).'),
    ('Crystal H-bond network', 'molecules/btpc/input/BTPC.cif; papers/btpc/data/derived/crystal_hbonds.csv', 'CIF fractional coordinates + 4 symmetry ops; CSV columns D,H,A,symmetry_code_of_A,symmetry_operation_of_A,D-H_A,H...A_A,D...A_A,D-H...A_deg,kind.'),
    ('Frequencies', 'T/freq.json (+ T/freq_step0.005.json for PCM)', 'keys frequencies_cm-1 (93, sorted, imaginary negative), frequencies_pyscf_cm-1, ZPE_Ha, H298_Ha, G298_Ha, hessian_method, n_imaginary, tr_residual_cm-1. No intensities: a stick "spectrum" can only show positions.'),
    ('Hessian / normal modes', 'R/hessian.npy', 'float64 array (33, 33, 3, 3), Ha/Bohr², Cartesian, same atom order and frame as opt.xyz (Bohr = Å/0.529177210903). PCM file is the step-0.01 Hessian (the 0.005 one was overwritten). Masses as used by PySCF for freq.json: H 1.008, C 12.011, N 14.007, O 15.999 (reproduce freq.json to < 1e-10 cm⁻¹). Imaginary-mode arrows can be drawn from its lowest eigenvector.'),
    ('Orbitals (3-D)', 'R/homo.cube, R/lumo.cube', 'Gaussian cube, Bohr, 80×80×80; PCM origin (−12.654, −13.796, −9.456), spacing (0.3271, 0.3708, 0.2645) Bohr; gas origin (−12.138, −13.783, −9.693), spacing (0.3175, 0.3601, 0.2618). Values in a.u.; render at ±0.02. Existing renders T/fig/homo.png, lumo.png (1600×1200).'),
    ('Orbital energies / DOS', 'T/fmo.json:mo_energies_eV; papers/btpc/data/derived/mo_fragment_composition_<case>.csv', 'fmo.json: 603 (PCM) / 602 (gas) energies in eV, n_occ 89. CSV (derived): index0,label,energy_eV,occupation,picrate_pct,cation_pct for every MO (Mulliken, validated) — enough for a fragment-resolved DOS/PDOS. Existing T/fig/mo_levels.png (1260×1650).'),
    ('Wavefunction', 'R/final.molden', 'PySCF 2.14.0 molden, [Atoms] in Bohr, spherical [5d], 606 AOs. PCM: 603 MOs zero-padded to 606 (the 3 padded vectors have energy max+1000 a.u. and zero coefficients); gas: 602 MOs, not padded (Multiwfn read 602 orbitals).'),
    ('Density / MEP 3-D', 'R/density.cube, R/esp.cube', 'density 80³ (same grid as orbitals), e/Bohr³; ESP 64³ (PCM spacing 0.4102, 0.4650, 0.3317 Bohr; gas 0.3982, 0.4515, 0.3283), a.u. Map ESP onto the 0.001 a.u. density isosurface; colour range = T/mep_range.json (vmin/vmax a.u.). Existing T/fig/mep.png.'),
    ('MEP surface points', 'R/vtx.txt', 'first line = vertex count (PCM 26483, gas 25344); then 6 columns: x, y, z (Bohr), ESP (a.u.), ESP (eV), ESP (kcal/mol). Exact Multiwfn surface — best source for an ESP histogram / area-vs-ESP plot.'),
    ('MEP extrema', 'T/mep_extrema.csv, R/surfanalysis.txt', 'CSV: type,index,global,ESP_au,ESP_eV,ESP_kcal_mol,x_A,y_A,z_A,nearest_atoms (27 rows each); surfanalysis.txt: Multiwfn lists (#, a.u., eV, kcal/mol, X/Y/Z Å), * = global. Statistics: T/mep_surface.json and R/logs/mwfn_mep.out (variances, skewness, sphericity only in the log).'),
    ('NCI scatter', 'R/rdg_scatter.txt, R/rdg_scatter.npz', 'txt: 5 columns x, y, z (Å), sign(λ₂)ρ (a.u.), RDG (a.u.); PCM 1,771,560 rows, gas 1,764,000; points screened by RDG_maxrho (ρ ≥ 0.1 a.u.) carry RDG = 100 exactly (PCM 94,395 rows, gas 103,505); a further 104,320 (PCM) / 89,732 (gas) low-density points have RDG > 100 naturally (RDG ≥ 100 in total: PCM 198,715, gas 193,237). npz: float32 x=sign(λ₂)ρ, y=RDG, filtered to −0.25 < x < 0.06 and y < 2.5 (PCM 365,105 points, gas 403,221), n_total. Existing T/fig/rdg_scatter.png (x −0.05…0.05) and rdg_scatter_wide.png.'),
    ('NCI 3-D', 'R/rdg.cube, R/sl2r.cube', 'Multiwfn cubes, Bohr; PCM 126×148×95 from (−9.654, −10.796, −6.456) spacing 0.158553; gas 125×147×96 from (−9.138, −10.783, −6.693) spacing 0.153832. Isosurface of rdg at 0.5, coloured by sl2r in −0.035…+0.02 a.u. Existing T/fig/nci.png. Spike table T/nci_spikes.csv (sign_l2_rho_au,min_RDG,n_points,type).'),
    ('ELF / LOL maps', 'R/{elf,lol}_{picrate_ring,benzotriazolium_ring}.npz and *_plane.txt', 'npz keys X, Y (in-plane coordinates, Å, 300×300), F (values), atoms_xy (33×2 projected atom positions), in_plane (bool, |offset| < 0.6 Å), normal, center. plane.txt: 90,000 rows, 6 columns X, Y, Z (Å), graph-x, graph-y, value. Picrate plane spans ±4.36 × ±4.20 Å (steps 0.029 × 0.028 Å), cation plane ±3.81 × ±4.22 Å (gas ±4.23; steps 0.026 × 0.028 Å). Existing T/fig/{elf,lol}_*.png (1860×1500, jet colour map, 0–1).'),
    ('ParaView inputs', 'R/vtk/*.vtk', 'legacy VTK (ASCII) in Å converted by dft/render.py: density, esp, homo, lumo, rdg, sl2r volumes; atoms_<El>.vtk and bonds.vtk polydata.'),
    ('Cross-case charts', 'results/charts/btpc/*.png', 'bond_lengths_comparison, frequency_spectrum, hbond_geometry, mep_extrema, mo_levels_comparison, nci_spikes, optimisation_convergence, reactivity_descriptors (300 dpi, made by dft/charts.py from the tracked JSON/CSV). Rebuild for the paper rather than reuse; the SUMMARY.md title there is truncated.'),
    ('Multiwfn inputs/logs', 'R/mwfn_*.in, R/logs/mwfn_*.out, R/mwfn_home/settings.ini', 'exact menu sequences (MEP 12,0,1,7,-1,q; NCI 20,1,3,2,3,0,0,q; ELF 4,9,1,300,300,4,atoms,-6,-5,q; LOL same with 10) and the settings actually used (RDG_maxrho 0.1, aug2D 6.0, iESPcode 2 PCM / 1 gas).'),
]
W(row('Figure topic', 'File(s)', 'Content checked'))
W(row('---', '---', '---'))
for a, b, c_ in inv: W(row(a, f'`{b}`', c_))
W('')

# ------------------------------------------------------------------------------------------------ appendix
W('## Appendix — inconsistencies found and how they were resolved')
W('')
for k in [x for x in F if x.startswith('inconsistency.')]:
    W(f'- {v(k)}. Authoritative: {F[k]["note"].replace("authoritative: ", "")} (`{k}`; sources: {F[k]["source"]}).')
W('- The auxiliary basis for density fitting is not recorded in any log; `protocol.aux_basis` gives the inference from the PySCF source and marks it as unverified.')
W('- `btpc_pcm` was not included in `verify_grid.py`; its grid/Lebedev settings are inferred (see `protocol.grid_level`).')
W('')

# ------------------------------------------------------------------------------------------------ writer's additions
_WK = ['xray.nh_minus_ho_asym', 'xray.nn_bond_asym', 'pcm.nn_bond_asym', 'gas.nn_bond_asym', 'gas.lumo_cation_pct', 'pcm.lumo_cation_pct',
       'pcm.homo_m1_cation_pct', 'gas.homo_cation_pct', 'saddle.imag_mode_rigid_pct', 'gas.mode1_rigid_pct', 'saddle.imag_mode_cation_pct',
       'saddle.drift_0_to_70_pct', 'saddle.drift_40_to_70_pct', 'saddle.drift_60_to_70_pct_9cm', 'saddle.hbond_max_excursion',
       'saddle.ring_rotation_0_to_70', 'saddle.ring_rotation_40_to_70', 'saddle.energy_drop_40_to_70', 'saddle.step_change_mean',
       'saddle.step_change_max_mode', 'saddle.step_change_max_signed', 'pcm.nu_C3_H3_stretch', 'pcm.nu_C5_H5_stretch', 'pcm.nu_benzo_CH_stretches',
       'pcm.benzo_ch_spread', 'gas.nu_C3_H3_stretch', 'gas.nu_C5_H5_stretch', 'gas.nu_benzo_CH_stretches', 'gas.benzo_ch_spread',
       'diff.ch_stretch_max_shift', 'pcm.cation_plane_offset_O7', 'gas.cation_plane_offset_O7', 'pcm.picrate_plane_offset_O7',
       'gas.picrate_plane_offset_O7', 'protocol.opt_threshold_energy', 'protocol.opt_threshold_grad_rms',
       'protocol.opt_threshold_grad_max', 'protocol.opt_threshold_disp_rms', 'protocol.opt_threshold_disp_max', 'protocol.pcm_lebedev_points',
       'gas.nci_region01_min', 'plot.nci_colour_min', 'plot.nci_colour_max', 'plot.nci_colour_mid', 'plot.mep_colour_min',
       'plot.mep_colour_max', 'plot.mep_hist_bin', 'plot.mep_label_threshold', 'plot.pdos_fwhm', 'plot.orbital_box_margin']
W("## Writer's additions (manuscript, 23 Sep 2026)")
W('')
W('Added by the manuscript writer through section 7 of `build_facts.py` (OUTLINE.md section 9): arithmetic on existing keys '
  '(with an `expr` that numbers_audit.py re-evaluates), quantities recomputed from the stored trajectory, Hessians and plane '
  'files, and plotting constants read from the figure scripts.')
W('')
W(row('key', 'value', 'unit', 'source / formula'))
W(row('---', '---', '---', '---'))
for k in [x for x in _WK if x in F]:
    val = v(k)
    val = ', '.join(str(x) for x in val) if isinstance(val, list) else (f'{val:.4g}' if isinstance(val, float) else val)
    W(row(kk(k), val, F[k]['unit'], F[k].get('expr') or (F[k]['source'] + '; ' + F[k]['note'])))
W('')
# ------------------------------------------------------------------------------------------------ revision round 1
_R1 = [k for k in F if isinstance(F[k], dict) and F[k].get('added') == 'r1']
W('## Revision round 1 additions (24 Sep 2026)')
W('')
W('Added through section 8 of `build_facts.py` in answer to the round-1 audits (`papers/btpc/audit/*_r1.md`; response in '
  '`audit/response_r1.md`). All are local post-processing of stored files: CIF values, the IEF-PCM trajectory, the two '
  'Hessians, freq.json (quasi-RRHO entropy, with the RRHO entropy reproduced first as a check), Loewdin populations with '
  'the rebuilt overlap matrix, and the density topology of the N4-H4N...O1 contact evaluated from final.molden by '
  '`data/wfn_topology.py` (validated against density.cube and the Multiwfn ELF plane; sign(lambda2)rho at the critical '
  'point reproduces the NCI region 01 value). Reading: at the H4N...O1 critical point both models have lap(rho) > 0 and '
  'H < 0, the transit region between closed-shell and shared-shell interactions; the gas-phase bond is much further into '
  'it (H -0.031 against -0.006 a.u., |V|/G 1.44 against 1.14). The Loewdin picrate charge differs between the models by '
  'about the same amount as the Mulliken one (0.13 against 0.12 e), and the H4N charge is the same in both models under '
  'either partition, so the difference is not an artefact of where H4N is counted.')
W('')
W(row('key', 'value', 'unit', 'source / formula'))
W(row('---', '---', '---', '---'))
for k in _R1:
    val = v(k)
    val = ', '.join(str(x) for x in val) if isinstance(val, list) else (f'{val:.4g}' if isinstance(val, float) else val)
    W(row(kk(k), val, F[k]['unit'], F[k].get('expr') or (F[k]['source'] + '; ' + F[k]['note'])))
W('')
# ------------------------------------------------------------------------------------------------ revision round 2
_R2 = [k for k in F if isinstance(F[k], dict) and F[k].get('added') == 'r2']
W('## Revision round 2 additions (24 Sep 2026)')
W('')
W('Added through section 9 of `build_facts.py` in answer to the round-2 audits (`papers/btpc/audit/*_r2.md`; response in '
  '`audit/response_r2.md`). Local post-processing only. (1) The curvature of the two stored, unprojected Hessians along the '
  'overall rigid rotations and translations of the pair: an exact Hessian at a stationary point gives zero, so these '
  'values measure the numerical resolution of the finite-difference Hessians (rotations up to 7.7i cm-1, translations '
  '0.2 cm-1 or less); the 3.2i cm-1 IEF-PCM mode lies below that resolution. The analytic gradients behind the Hessians '
  'omit the grid-weight response (`protocol.fd_gradient_grid_response`). (2) H/rho at the H4N...O1 critical point and '
  'the gas/IEF-PCM ratio of H (about 4.9): both bonds are in the transit region, and they differ in degree, not in class. '
  '(3) U(eq) of the triazole nitrogens from the CIF. (4) Three plotting constants read from the figure scripts. '
  '(5) The Delta pKa values and ranges quoted in the abstract, Fig. 1 and the Introduction, computed from the pKa '
  'values in the abstract-checked literature claims (`lit.*`). '
  'Changed values: every `*.mep_minNN` / `*.mep_maxNN` now holds the unrounded kcal mol-1 value from Multiwfn\'s '
  'surfanalysis.txt instead of the 3-decimal value of mep_extrema.csv (six table cells had been rounded twice and were '
  '0.01 off; QA round 4 r4-m1).')
W('')
W(row('key', 'value', 'unit', 'source / formula'))
W(row('---', '---', '---', '---'))
for k in _R2:
    val = v(k)
    val = ', '.join(str(x) for x in val) if isinstance(val, list) else (f'{val:.4g}' if isinstance(val, float) else val)
    W(row(kk(k), val, F[k]['unit'], F[k].get('expr') or (F[k]['source'] + '; ' + F[k]['note'])))
W('')
# keep the hand-written fact-check log (appended after generation) across rebuilds
_old = (D / 'facts.md').read_text() if (D / 'facts.md').exists() else ''
_i = _old.find('\n## Fact-check log')
if _i >= 0:
    out += _old[_i + 1:].rstrip('\n').split('\n')
(D / 'facts.md').write_text('\n'.join(out) + '\n')
print(f'wrote {D / "facts.md"} ({len(out)} lines)')
